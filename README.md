TiSeClAN
=======
Overview 
=======
This script performs comprehensive analysis of gene expression time series data, including data preprocessing, differential expression analysis, clustering, functional enrichment analysis, and various visualizations. The pipeline is designed to identify significant temporal expression patterns and provide biological insights through multiple analysis steps. 

Requirements 
=======
    Python 3.8+
    Required libraries:
        pandas
        numpy
        matplotlib
        plotly
        scikit-learn
        tslearn
        gseapy
        mygene
        xlsxwriter
        scipy
        statsmodels
        umap-learn
        seaborn
        BioPython
 
 
Input Files 
==============
Counts file  (sorted.norm.CPMs.tsv): 
        Tab-separated file with gene counts
        Rows: Samples
        Columns: Genes(ENSEMBL only!)
         
  Metadata file  (phenotypes.tsv): 
        Tab-separated file with sample information
        Must contain 'Time' column with numeric time points

Usage 
==============
    python dynamics3.py \
    --counts_file path/to/counts_file.tsv \
    --metadata_file path/to/metadata_file.tsv \
    --output_dir path/to/output_directory \
    --fdr_threshold 0.1 \
    --min_abs_coeff 0.3 \
    --poly_degree 3 \
    --min_fold_change 1.5 \
    --min_variance 0.1 \
    --smoothing_window 5 \
    --smoothing_method savgol \
    --enable_umap \
    --auto_eps \
    --enable_blast
 
 
Parameters 
==============
Main Parameters 
**********
    --counts_file: Path to normalized counts file
    --metadata_file: Path to metadata file
    --output_dir: Directory to save results
    --fdr_threshold: FDR threshold for significance (default: 0.1)
    --min_abs_coeff: Minimum absolute coefficient of change (default: 0.3)
    --poly_degree: Degree of polynomial for regression (default: 0, 0 to disable)
     

Optional Parameters 
******
    --min_fold_change: Minimum fold change between time points (default: 1.5)
    --min_variance: Minimum variance across samples (default: 0.3)
    --additional_filtering: Enable additional filtering by fold change and variance
    --aggregate_samples: Aggregate samples by time point
    --smoothing_window: Window size for smoothing (0 to disable)
    --smoothing_method: Smoothing method ('savgol' or 'ma')
    --enable_umap: Generate UMAP projections
    --auto_eps: Auto-selection of DBSCAN eps parameter
    --enable_blast: Enable fetching BLAST names
     

Output Files 
==============
The script generates the following outputs in the specified output directory: 

Main Results 
************
    cluster_analysis_filtered.xlsx: Contains: 
        All genes statistics
        Filter information
        Cluster metrics
        Gene lists for each cluster
        GO enrichment results
         

    summary.xlsx: Summary table with cluster metrics 
     

Visualizations 
************
    Cluster profile plots (cluster_X_profiles.png)
    Heatmaps (cluster_X_heatmap.png)
    UMAP projections (umap_projection.png, 3d_umap_projection.html)
    Enrichment visualizations:
        Barplots (cluster_X_GO_Biological_Process_2023_barplot.png)
        Dotplots (cluster_X_GO_Biological_Process_2023_dotplot.png)
        Scatterplots (cluster_X_GO_Biological_Process_2023_scatterplot.png)
         
     

Additional Files 
***********
    available_libraries.txt: List of available enrichment databases
     

Analysis Workflow 
==============
Data Preprocessing  
    Normalization using StandardScaler
    Optional smoothing (Savitzky-Golay or moving average)
    ***********
 Differential Expression Analysis(optional)  
    Polynomial regression modeling
    Multiple testing correction (Benjamini-Hochberg)
    FDR threshold filter
    Minimum coefficient change filter
   Additional Filtering by:
    Fold change
    Variance
             
Clustering Analysis  
     DBSCAN clustering
    Cluster quality assessment using:
        R²
        Dynamic Time Warping (DTW)
        Autocorrelation
        Composite score
             
         

Functional Enrichment Analysis  
    Multiple databases:
            GO Biological Process
            KEGG
            Reactome
            WikiPathways
            DisGeNET
            Human Gene Atlas
             
Visualization:
    Barplots
    Dotplots
    Scatterplots
               
Visualization  
    UMAP projections (2D and 3D)
    Cluster profile plots
        Heatmaps
         

Optional Features  
    BLAST annotation retrieval
    Summary dotplots for enrichment results
         
     

Notes 
==============
    The script automatically creates necessary directories and handles file naming conflicts
    All visualizations are saved in high resolution (300-600 dpi)
    Results are organized in a structured Excel file with separate sheets for different analyses
    The script includes error handling for network-dependent operations (e.g., Enrichr queries)
     


Contact 

For questions or support, please contact [vsevolodvolodin@yandex.ru] 
License 

Version 

1.0 
Change Log 

    Initial release
     
