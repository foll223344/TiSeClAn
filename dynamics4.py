import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from sklearn.cluster import DBSCAN
from tslearn.metrics import dtw
from gseapy import enrichr
import mygene
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
import xlsxwriter
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import warnings
import os
import umap.umap_ as umap 
from sklearn.manifold import TSNE
import seaborn as sns
from scipy.signal import savgol_filter
from sklearn.utils import resample
from sklearn.neighbors import NearestNeighbors
from Bio.Blast import NCBIWWW, NCBIXML
from gseapy import get_library_name
import networkx as nx
from PIL import Image

# Suppress warnings
warnings.filterwarnings('ignore')

# Add command line arguments
parser = argparse.ArgumentParser(description='Gene Expression Time Series Analysis')
parser.add_argument('--counts_file', type=str, required=True, help='Path to the counts file (sorted.norm.CPMs.tsv)')
parser.add_argument('--metadata_file', type=str, required=True, help='Path to the metadata file (phenotypes.tsv)')
parser.add_argument('--fdr_threshold', type=float, default=0.1, help='FDR threshold for significance')
parser.add_argument('--min_abs_coeff', type=float, default=0.3, help='Minimum absolute coefficient of change')
parser.add_argument('--poly_degree', type=int, default=0, help='Degree of polynomial for non-linear regression')
parser.add_argument('--output_dir', type=str, required=True, help='Directory to save all results')
parser.add_argument('--min_fold_change', type=float, default=1.5, help='Minimum fold change between time points')
parser.add_argument('--min_variance', type=float, default=0.1, help='Minimum variance across samples')
parser.add_argument('--additional_filtering', action='store_true', help='Enable additional filtering by fold change and variance')
parser.add_argument('--aggregate_samples', action='store_true', help='Aggregate samples by time point')
parser.add_argument('--smoothing_window', type=int, default=0, help='Window size for smoothing (0 to disable)')
parser.add_argument('--smoothing_method', choices=['savgol', 'ma'], default='savgol', help='Smoothing method: savgol (Savitzky-Golay) or ma (moving average)')
parser.add_argument('--enable_umap', action='store_true', help='Generate UMAP projections')
parser.add_argument('--min_samples', type=int, default=5, help='Minimum number of samples in a neighborhood for DBSCAN')
parser.add_argument('--eps_range', type=str, default='0.1,0.9,0.1', help='Range of eps values for DBSCAN auto-selection (start, stop, step)')
parser.add_argument('--auto_eps', action='store_true', help='Auto-selection of eps')
parser.add_argument('--enable_blast', action='store_true', help='Enable fetching BLAST names')
parser.add_argument('--eps', type=float, default=0.7, help='Eps value manual')
args = parser.parse_args()

# Use arguments in code
FDR_THRESHOLD = args.fdr_threshold
MIN_ABS_COEFF = args.min_abs_coeff
POLY_DEGREE = args.poly_degree
output_dir = args.output_dir
MIN_FOLD_CHANGE = args.min_fold_change
MIN_VARIANCE = args.min_variance
additional_filtering = args.additional_filtering
aggregate_samples = args.aggregate_samples
min_samples = args.min_samples
eps_start, eps_stop, eps_step = map(float, args.eps_range.split(','))

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Get a list of all available databases
all_libraries = get_library_name()
libraries_df = pd.DataFrame(all_libraries, columns=['Library_Name'])
text_file_path = os.path.join(output_dir, 'available_libraries.txt')
with open(text_file_path, 'w') as f:
    for library in all_libraries:
        f.write(f"{library}\n")

# Saving the list of databases to an Excel file
excel_file_path = os.path.join(output_dir, 'available_libraries.xlsx')
libraries_df.to_excel(excel_file_path, index=False)

print(f"Available libraries list saved to {text_file_path} and {excel_file_path}")

# Define paths
cluster_analysis_file = os.path.join(output_dir, 'cluster_analysis_filtered.xlsx')
summary_file = os.path.join(output_dir, 'summary.xlsx')

# 1. Data Loading
counts = pd.read_csv(args.counts_file, sep='\t', index_col=0).T
metadata = pd.read_csv(args.metadata_file, sep='\t', index_col=0)

# 2. Data Preprocessing
metadata['Time'] = pd.to_numeric(metadata['Time'], errors='coerce')
common_samples = counts.index.intersection(metadata.index)
counts = counts.loc[common_samples]
metadata = metadata.loc[common_samples].sort_values('Time')
counts = counts.loc[metadata.index]
time = metadata['Time'].values

# Smoothing function
def apply_smoothing(data, window_size, method):
    if window_size < 3:
        return data
    smoothed = pd.DataFrame(index=data.index, columns=data.columns)
    for gene in data.columns:
        ts = data[gene].values
        if method == 'savgol':
            smoothed[gene] = savgol_filter(ts, window_size, 3)
        elif method == 'ma':
            smoothed[gene] = pd.Series(ts).rolling(window_size, center=True).mean().values
    return smoothed.dropna()

if args.smoothing_window > 0:
    counts = apply_smoothing(counts, args.smoothing_window, args.smoothing_method)

# Aggregate samples
if aggregate_samples:
    def aggregate_by_time(counts, metadata):
        aggregated_counts = pd.DataFrame(index=counts.index)
        for time_point, group in metadata.groupby('Time'):
            sample_ids = group.index
            aggregated_counts[f'Time_{time_point}'] = counts[sample_ids].mean(axis=1)
        return aggregated_counts
    counts = aggregate_by_time(counts, metadata)
    metadata = pd.DataFrame({'Time': sorted(set(metadata['Time']))}, index=[f'Time_{t}' for t in sorted(set(metadata['Time']))])
    time = metadata['Time'].values

# Regression analysis
significant_genes = set(counts.columns.tolist())
print(f"Initial number of genes: {len(significant_genes)}")

if args.poly_degree > 0:
    def fit_regression_model(gene_expression, degree=POLY_DEGREE):
        poly = PolynomialFeatures(degree=degree)
        X_poly = poly.fit_transform(time.reshape(-1, 1))
        X_poly = sm.add_constant(X_poly)
        model = sm.OLS(gene_expression, X_poly)
        return model.fit().params[1:], model.fit().pvalues[1:]

    coeffs_list, pvals_list = [], []
    for gene in counts.columns:
        params, pvalues = fit_regression_model(counts[gene].values)
        coeffs_list.append(params)
        pvals_list.append(pvalues)

    time_coeffs = pd.DataFrame(coeffs_list, index=counts.columns, 
                               columns=[f'Time_Coeff_{i}' for i in range(POLY_DEGREE)])
    pvals = pd.DataFrame(pvals_list, index=counts.columns,
                         columns=[f'Pvalue_{i}' for i in range(POLY_DEGREE)])

    fdr_df = pd.DataFrame({f'FDR_{i}': multipletests(pvals.iloc[:, i], method='fdr_bh')[1] 
                           for i in range(len(pvals.columns))}, index=counts.columns)

    significant_genes_after_regression = set()
    for i in range(len(time_coeffs.columns)):
        mask = (fdr_df[f'FDR_{i}'] < FDR_THRESHOLD) & (np.abs(time_coeffs.iloc[:, i]) > MIN_ABS_COEFF)
        significant_genes_after_regression.update(counts.columns[mask])

    significant_genes = significant_genes.intersection(significant_genes_after_regression)
    print(f"Number of genes after regression: {len(significant_genes)}")

# Additional filtering
if additional_filtering:
    variances = counts.var(axis=0)
    fold_changes = counts.max(axis=0) / counts.min(axis=0)
    fold_change_mask = (fold_changes >= MIN_FOLD_CHANGE) & (fold_changes != np.inf)
    variance_mask = (variances >= MIN_VARIANCE) & (variances != 0)
    significant_genes = significant_genes.intersection(
        counts.columns[fold_change_mask & variance_mask]
    )
    print(f"Number of genes after additional filtering: {len(significant_genes)}")

significant_genes = list(significant_genes)

# Normalization
scaler = StandardScaler()
normalized_significant = pd.DataFrame(
    scaler.fit_transform(counts[significant_genes]),
    index=counts.index,
    columns=significant_genes
)

# Clustering metrics
def calculate_r_squared(cluster_genes, mean_profile):
    if len(cluster_genes) == 0:  
        return 0
    r_squared_sum = 0
    for gene in cluster_genes:
        y_true = normalized_significant[gene].values
        y_pred = mean_profile
        slope, intercept, r_value, p_value, std_err = stats.linregress(y_true, y_pred)
        r_squared_sum += r_value**2
    return r_squared_sum / len(cluster_genes)

def calculate_dtw(cluster_genes, mean_profile):
    return np.mean([dtw(normalized_significant[gene].values, mean_profile) 
                    for gene in cluster_genes]) if len(cluster_genes) > 0 else 0

def calculate_autocorrelation(cluster_genes, lag=1):
    return np.mean([np.corrcoef(ts[:-lag], ts[lag:])[0,1] 
                    for ts in normalized_significant[cluster_genes].T.values]) if len(cluster_genes) > 0 else 0

# Function to calculate the silhouette
def calculate_silhouette(clusters, data):
    from sklearn.metrics import silhouette_score
    unique_labels = set(clusters)
    if -1 in unique_labels:
        unique_labels.remove(-1)
    if len(unique_labels) < 2:
        return -1
    return silhouette_score(data, clusters)
    
# Auto-select DBSCAN parameters
def find_best_eps(data, min_samples, eps_range):
    best_eps = None
    best_silhouette = -1
    for eps in np.arange(*eps_range):
        db = DBSCAN(eps=eps, min_samples=min_samples)
        labels = db.fit_predict(data)
        unique_labels = set(labels)
        if -1 in unique_labels:
            unique_labels.remove(-1)
        if len(unique_labels) < 2:
            continue
        silhouette = calculate_silhouette(labels, data)
        if silhouette > best_silhouette:
            best_silhouette = silhouette
            best_eps = eps
    return best_eps

# Clustering using DBSCAN
def evaluate_clustering_dbscan(min_samples, eps_range=(0.1, 0.9, 0.1)):
    best_eps = None
    best_silhouette = -1
    best_clusters = None
    best_cluster_metrics = None
    for eps in np.arange(*eps_range):
        db = DBSCAN(eps=eps, min_samples=min_samples)
        clusters = db.fit_predict(normalized_significant.T.values)
        unique_labels = set(clusters)
        if -1 in unique_labels:
            unique_labels.remove(-1)
        if len(unique_labels) < 2:
            continue
        silhouette = calculate_silhouette(clusters, normalized_significant.T.values)
        if silhouette > best_silhouette:
            best_silhouette = silhouette
            best_eps = eps
            best_clusters = clusters
    print(f"Best eps: {best_eps}")
    db = DBSCAN(eps=best_eps, min_samples=min_samples)
    clusters = db.fit_predict(normalized_significant.T.values)
    cluster_metrics = {}
    for cluster_id in set(clusters):
        if cluster_id == -1:
            continue
        cluster_genes = normalized_significant.columns[clusters == cluster_id]
        if len(cluster_genes) == 0: 
            print(f"Cluster {cluster_id} is empty, skipping metrics calculation.")
            continue
        mean_profile = normalized_significant[cluster_genes].mean(axis=1).values
        r_squared = calculate_r_squared(cluster_genes, mean_profile)
        avg_dtw = calculate_dtw(cluster_genes, mean_profile)
        avg_autocorr = calculate_autocorrelation(cluster_genes)
        composite = (r_squared + (1 - avg_dtw/(avg_dtw+1)) + avg_autocorr) / 3
        cluster_metrics[cluster_id] = {
            'r_squared': r_squared,
            'dtw': avg_dtw,
            'autocorr': avg_autocorr,
            'composite': composite
        }
    return best_silhouette, clusters, cluster_metrics


if args.auto_eps:
    eps = find_best_eps(normalized_significant.T.values, min_samples, (eps_start, eps_stop, eps_step))
else:
    eps = args.eps
print(f"Selected eps: {eps}")
db = DBSCAN(eps=eps, min_samples=min_samples)
clusters = db.fit_predict(normalized_significant.T.values)
unique_labels = set(clusters)
if -1 in unique_labels:
    unique_labels.remove(-1)
n_clusters = len(unique_labels)
cluster_metrics = {}
for cluster_id in unique_labels:
    cluster_genes = normalized_significant.columns[clusters == cluster_id]
    mean_profile = normalized_significant[cluster_genes].mean(axis=1).values
    r_squared = calculate_r_squared(cluster_genes, mean_profile)
    avg_dtw = calculate_dtw(cluster_genes, mean_profile)
    avg_autocorr = calculate_autocorrelation(cluster_genes)
    composite = (r_squared + (1 - avg_dtw/(avg_dtw+1)) + avg_autocorr) / 3
    cluster_metrics[cluster_id] = {
        'r_squared': r_squared,
        'dtw': avg_dtw,
        'autocorr': avg_autocorr,
        'composite': composite
    }


# Visualization functions
def plot_3d_umap_projection(embedding, clusters):
    fig = go.Figure(data=[go.Scatter3d(
        x=embedding[:, 0],
        y=embedding[:, 1],
        z=embedding[:, 2],
        mode='markers',
        marker=dict(color=clusters, size=5, colorscale='Viridis', showscale=True),
        text=[f"Cluster {c}" for c in clusters]
    )])
    fig.update_layout(scene=dict(xaxis_title='UMAP 1', yaxis_title='UMAP 2', zaxis_title='UMAP 3'),
                      title='3D UMAP Projection with Clusters')
    fig.write_html(os.path.join(output_dir, '3d_umap_projection.html'))

def plot_heatmap(cluster_genes, cluster_id):
    plt.figure(figsize=(10, 15))
    sns.heatmap(normalized_significant[cluster_genes].T, cmap='viridis', yticklabels=False)
    plt.title(f'Cluster {cluster_id} Heatmap')
    plt.savefig(os.path.join(output_dir, f'cluster_{cluster_id}_heatmap.png'), dpi=300)
    plt.close()

def plot_umap_projection(clusters):
    reducer = umap.UMAP(n_components=2, random_state=42,n_neighbors=15, min_dist=0.9)
    embedding = reducer.fit_transform(normalized_significant.T.values)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=embedding[:,0], y=embedding[:,1], hue=clusters, palette='Spectral', alpha=0.7, size=clusters)
    plt.title('UMAP Projection')
    plt.savefig(os.path.join(output_dir, 'umap_projection.png'), dpi=300)
    plt.close()

def plot_cluster_profiles(time, normalized_significant, clusters, cluster_metrics):
    all_plots = []
    for cluster_id in unique_labels:
        if cluster_id == -1 or len(normalized_significant.columns[clusters == cluster_id]) == 0:
            continue
        cluster_genes = normalized_significant.columns[clusters == cluster_id]
        mean_profile = normalized_significant[cluster_genes].mean(axis=1).values
        
        # Fetch metrics for the current cluster
        r_squared = cluster_metrics[cluster_id]['r_squared']
        dtw_value = cluster_metrics[cluster_id]['dtw']
        autocorr = cluster_metrics[cluster_id]['autocorr']
        num_genes = len(cluster_genes) 
        
        plt.figure(figsize=(10, 6))
        plt.plot(time, mean_profile, 'k-', lw=2, label=f'Cluster {cluster_id} Mean (Genes: {num_genes})')
        
        # Add individual gene profiles with alpha transparency
        for gene in np.random.choice(cluster_genes, size=min(20, len(cluster_genes)), replace=False):
            plt.plot(time, normalized_significant[gene], alpha=0.3)
        
        # Create a title that includes all three metrics and number of genes
        plt.title(f'Cluster {cluster_id} Profiles\n'
                  f'R²: {r_squared:.2f}, DTW: {dtw_value:.2f}\n'
                  f'Number of Genes: {num_genes}')
        
        # Save the figure
        plt_path = os.path.join(output_dir, f'cluster_{cluster_id}_profiles.png')
        plt.savefig(plt_path, dpi=600)
        plt.close()
        all_plots.append(plt.imread(plt_path))
    
    return all_plots

# Plotting cluster profiles and saving plots
all_plots = plot_cluster_profiles(time, normalized_significant, clusters, cluster_metrics)

# Combine all plots into one image
if all_plots:
    fig, axes = plt.subplots(1, len(all_plots), figsize=(len(all_plots) * 5, 5))
    for ax, img in zip(axes, all_plots):
        ax.imshow(img)
        ax.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'all_cluster_profiles.png'), dpi=600)
    plt.close()

# Generate visualizations
writer = pd.ExcelWriter(cluster_analysis_file, engine='xlsxwriter')
gene_stats = pd.DataFrame({
    'Gene': counts.columns,
    **(time_coeffs if args.poly_degree > 0 else {}),
    **(pvals if args.poly_degree > 0 else {}),
    **(fdr_df if args.poly_degree > 0 else {})
}).set_index('Gene')
gene_stats.to_excel(writer, sheet_name='All_genes_stats')

if args.enable_umap:
    plot_umap_projection(clusters)
    reducer = umap.UMAP(n_components=3, random_state=42)
    embedding = reducer.fit_transform(normalized_significant.T.values)
    plot_3d_umap_projection(embedding, clusters)

# Save the information about the number of genes at each filtering step in Excel

filter_info = pd.DataFrame({
    'Stage': ['Initial', 'After Additional Filtering'],
    'Number of Genes': [len(counts.columns), len(significant_genes)]
})
filter_info.to_excel(writer, sheet_name='Filter_Information')

# Plotting cluster profiles
plot_cluster_profiles(time, normalized_significant, clusters, cluster_metrics)

def shorten_sheet_name(base_name, cluster_id):
    base_name = base_name.replace(' ', '_')
    return f'Cluster_{cluster_id}_{base_name[:15]}'

def create_dotplot(top_terms, gene_set, cluster_id):
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=-np.log10(top_terms['Adjusted P-value']), y=top_terms['Term'], size=top_terms['Overlap'], hue=top_terms['Adjusted P-value'], palette='viridis')
    plt.title(f'Top Terms for Cluster {cluster_id} ({gene_set})')
    plt.legend(title='Adjusted P-value', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'cluster_{cluster_id}_{gene_set.replace(" ", "_")}_dotplot.png'), dpi=300)
    plt.close()

from PIL import Image

def create_dotplot(top_terms, gene_set, cluster_id):
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=-np.log10(top_terms['Adjusted P-value']), y=top_terms['Term'], size=top_terms['Overlap'], hue=top_terms['Adjusted P-value'], palette='viridis')
    plt.title(f'Top Terms for Cluster {cluster_id} ({gene_set})')
    plt.legend(title='Adjusted P-value', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f'cluster_{cluster_id}_{gene_set.replace(" ", "_")}_dotplot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    return plot_path

# Function to create a summary dotplot for one database ( not working)
def create_summary_dotplot_for_gene_set(gene_set, clusters, output_dir):
    image_files = []
    for cluster_id in clusters:
        if cluster_id == -1 or len(normalized_significant.columns.values[clusters == cluster_id]) == 0:
            continue
        plot_path = os.path.join(output_dir, f'cluster_{cluster_id}_{gene_set.replace(" ", "_")}_dotplot.png')
        if os.path.exists(plot_path):
            image_files.append(plot_path)
    
    if not image_files:
        print(f"No dotplot images found for gene set {gene_set}. Skipping summary dotplot.")
        return
    

    num_images = len(image_files)
    images_per_row = 3
    rows = (num_images + images_per_row - 1) // images_per_row
    
    images = [Image.open(image_file) for image_file in image_files]
    widths, heights = zip(*(i.size for i in images))

    max_width = max(widths)
    total_height = sum(heights[:images_per_row]) * rows

    new_im = Image.new('RGB', (max_width * images_per_row, total_height))

    x_offset = 0
    y_offset = 0
    for i, im in enumerate(images):
        new_im.paste(im, (x_offset, y_offset))
        x_offset += max_width
        if (i + 1) % images_per_row == 0:
            y_offset += im.size[1]
            x_offset = 0
    summary_plot_path = os.path.join(output_dir, f'{gene_set.replace(" ", "_")}_summary_dotplots.png')
    new_im.save(summary_plot_path)
    print(f"Summary dotplot for {gene_set} saved to {summary_plot_path}")


# GO analysis with multiple bases
all_gene_sets = [
    'GO_Biological_Process_2023', 
    'KEGG_2021_Human', 
    'Reactome_2022', 
    'WikiPathways_2024_Human',
    'GO_Cellular_Component_2023',
    'GO_Molecular_Function_2023',
    'DisGeNET',
    'Human_Gene_Atlas',
    'Reactome_Pathways_2024',
]
mg = mygene.MyGeneInfo()

for gene_set in all_gene_sets:
    # Create a separate Excel writer for non-GO annotations
    if gene_set != 'GO_Biological_Process_2023':
        non_go_writer = pd.ExcelWriter(os.path.join(output_dir, f'{gene_set}_annotations.xlsx'), engine='xlsxwriter')
    
    for cluster_id in unique_labels:
        if cluster_id == -1 or len(normalized_significant.columns[clusters == cluster_id]) == 0:
            continue
        cluster_genes = normalized_significant.columns[clusters == cluster_id]
        try:
            gene_symbols = mg.querymany(cluster_genes, scopes='ensembl.gene', fields='symbol', species='human')
            symbols = [gene['symbol'] for gene in gene_symbols if 'symbol' in gene]
            if len(symbols) == 0:
                raise ValueError("No symbols found")
            
            enr = enrichr(gene_list=symbols, gene_sets=[gene_set], organism='human')
            if enr.results is not None:
                top_terms = enr.results.sort_values('Adjusted P-value').head(10)
                
                # Barplot with adjusted formatting
                plt.figure(figsize=(10, 8))
                sns.barplot(x=-np.log10(top_terms['Adjusted P-value']), y='Term', data=top_terms, dodge=False)
                plt.xticks(rotation=45, ha='right')
                plt.title(f'Top Terms for Cluster {cluster_id} ({gene_set})')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, f'cluster_{cluster_id}_{gene_set.replace(" ", "_")}_barplot.png'), dpi=300)
                plt.close()
                
                # Dotplot
                create_dotplot(top_terms, gene_set, cluster_id)
                
                if gene_set == 'GO_Biological_Process_2023':
                    sheet_name = shorten_sheet_name(gene_set, cluster_id)
                    enr.results.head(1000).to_excel(writer, sheet_name=sheet_name)
                    
                    # Scatterplot for GO annotation
                    plt.figure(figsize=(10, 8))
                    sns.scatterplot(x=-np.log10(top_terms['Adjusted P-value']), y=range(len(top_terms)), hue=top_terms['Term'])
                    plt.title(f'Top Terms for Cluster {cluster_id} ({gene_set})')
                    plt.savefig(os.path.join(output_dir, f'cluster_{cluster_id}_{gene_set.replace(" ", "_")}_scatterplot.png'), dpi=300)
                    plt.close()
                
                else:
                    # Save non-GO annotations to separate files
                    sheet_name = shorten_sheet_name(gene_set, cluster_id)
                    enr.results.head(1000).to_excel(non_go_writer, sheet_name=sheet_name)
                

        except Exception as e:
            print(f"Error in cluster {cluster_id}: {str(e)}")
            continue
    
    if gene_set != 'GO_Biological_Process_2023':
        non_go_writer.close()

    create_summary_dotplot_for_gene_set(gene_set, unique_labels, output_dir)

print("GO annotation and visualization completed.")

for cluster_id in unique_labels:
    if cluster_id == -1 or len(normalized_significant.columns[clusters == cluster_id]) == 0:
        continue
    cluster_genes = normalized_significant.columns[clusters == cluster_id]
    

    cluster_genes_df = pd.DataFrame({
        'Ensembl_ID': cluster_genes,
        'Cluster_ID': [cluster_id] * len(cluster_genes)
    })
    
    cluster_genes_df.to_excel(writer, sheet_name=f'Cluster_{cluster_id}_genes', index=False)

if args.enable_blast:
    from mygene import MyGeneInfo

if args.enable_blast:
    from mygene import MyGeneInfo

    def get_local_blast_names(genes):
        mg = MyGeneInfo()
        blast_names = {}
        
        gene_info = mg.querymany(genes, scopes='ensembl.gene', fields='symbol,name,refseq', species='human')
        
        for gene in gene_info:
            if 'symbol' in gene and gene['symbol']:
                blast_names[gene['query']] = gene['symbol']
            elif 'name' in gene and gene['name']:
                blast_names[gene['query']] = gene['name']
            elif 'refseq' in gene and gene['refseq'] and 'rna' in gene['refseq']:
                blast_names[gene['query']] = gene['refseq']['rna']
            else:
                blast_names[gene['query']] = "No annotation found"
        
        return blast_names

    for cluster_id in unique_labels:
        if cluster_id == -1 or len(normalized_significant.columns[clusters == cluster_id]) == 0:
            continue
        
        cluster_genes = normalized_significant.columns[clusters == cluster_id]
        
        blast_names = get_local_blast_names(cluster_genes)
        
        blast_df = pd.DataFrame(list(blast_names.items()), columns=['Ensembl_ID', 'Annotation'])
        
        sheet_name = f'Cluster_{cluster_id}_BLAST_Names'
        blast_df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    print("BLAST names saved to the Excel file.")

metrics_df = pd.DataFrame.from_dict({k: v for k, v in cluster_metrics.items()}, orient='index')
metrics_df.to_excel(writer, sheet_name='Cluster_Metrics')

print("\nChanges in cluster metrics:")
print(metrics_df)

# Summary file
summary_writer = pd.ExcelWriter(os.path.join(output_dir, 'summary.xlsx'), engine='xlsxwriter')
summary_df = pd.DataFrame(columns=['Cluster_ID', 'Ensembl_ID', 'R_Squared', 'DTW', 'Autocorrelation', 'Composite'])

for cluster_id in unique_labels:
    if cluster_id == -1:
        continue
    cluster_genes = normalized_significant.columns[clusters == cluster_id]
    if len(cluster_genes) == 0:
        continue
    temp_df = pd.DataFrame({
        'Cluster_ID': [cluster_id] * len(cluster_genes),
        'Ensembl_ID': cluster_genes,
        'R_Squared': [cluster_metrics[cluster_id]['r_squared']] * len(cluster_genes),
        'DTW': [cluster_metrics[cluster_id]['dtw']] * len(cluster_genes),
        'Autocorrelation': [cluster_metrics[cluster_id]['autocorr']] * len(cluster_genes),
        'Composite': [cluster_metrics[cluster_id]['composite']] * len(cluster_genes)
    })
    summary_df = pd.concat([summary_df, temp_df], ignore_index=True)

summary_df.to_excel(summary_writer, sheet_name='Summary', index=False)
summary_writer.close()

writer.close()
