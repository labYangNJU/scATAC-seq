# scATAC-seq analysis

#All codes are run on R(3.5.1) and Python(2.7.10).

#The raw scATAC-seq reads to generate files in "data" are available in the GEO under accession number GSE166547.

The codes in this repository cover the following sections:

#Section 1: Barcode checking of raw reads against 10x genomics barcode whitelist.

code: 

      01_sc_atac_barcode_check_NoEdit.py


#Section 2: Remove duplicates on a cell-by-cell basis on filtered bam files.

code: 

      02_sc_atac_dedup.py


#Section 3: Identify barcodes representing genuine cells.

code: 
      
      03_Determine_ReadDepth_TssRatio.R   #Filter cells based on number of total fragments count (ReadDepth) and the ratio of fragments in TSS region (TssRatio).
      
      04_get_insert_size_distribution_per_cell.py   #Get insert size distribution per cell.
      
      05_Determine_insert_sizes_banding_score.R   #Calculate the periodicity in the frequency of insert sizes and filter out cells based on banding score.
      
      
#Section 4: Clustering analysis for scATAC-seq.

code: 

      06_snapatac_twobatch_hg38.R   #Used with SnapATAC(1.0.0).


#Section 5: Constructing peak-cell-matrix.

code: 

      07_Construct_SparseMatrix.py   #Need input files including genuine cell barcode list, peak list and deduplicated bam files.


#Section 6: Idenyify high quality, high variable and specific peaks.

code: 

      08_HighQuality_and_Variable_Peaks.R   #Idenyify high quality and high variable peaks. Need input file: peak-by-cluster proportion matrix.
      
      09_SpecificPeak_calculator.R   #Identify specific peaks. Need input files including peak-by-cell-type proportion matrix and cell-type medeian read depth file.
      
      
#Section 7: Generate normalized peak-by-cell matrix as pRCC.

code: 

      10_Generate_normalized_peak-by-cell_matrix_as_pRCC.R   #Need input files including deduplicated bam files per cluster and pan-cancer ATAC-seq peak list.


#Section 8: Optimized cicero calling and calculating gene acivities: 

code: 

      11_Optimized_cicero_calling.R   #Input files including peak-cell-matrix and cell grouping information.  Based on Cicero(1.0.15).
      

#Section 9: Identify origin-derived features for pRCC subtypes: 

code:

      12_Identify_origin_derived_features.R   #Input files including gene activities matrices from origin cells and pRCC. Based on Seurat(3.0.0) and edgeR(3.24.3).
     
     
#Section 10: Trajectory analysis between pRCC and normal samples.

code:

      13_Trajectory_analysis.R   #Input file including gene expression matrix for pRCC and normal samples and group information for all the samples. Based on Monocle(2.10.1).
      

# Citation

# Reference
https://github.com/shendurelab/fly-atac
https://github.com/shendurelab/mouse-atac
