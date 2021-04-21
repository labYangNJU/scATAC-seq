# scATAC-seq analysis

#All codes are run on R(3.5.1) and Python(2.7.10).
#The raw scATAC-seq reads to generate files in "data" are available in the GEO under accession number GSE166547.

The codes in this repository cover the following sections:

#Section 1: Barcode checking of raw reads against 10x genomics barcode whitelist.
code: 01_sc_atac_barcode_check_NoEdit.py

#Section 2: Remove duplicates on a cell-by-cell basis on filtered bam files.
code: 02_sc_atac_dedup.py

#Section 3: Identify barcodes representing genuine cells.
code: 03_Determine_ReadDepth_TssRatio.R   #Filter cells based on number of total fragments count (ReadDepth) and the ratio of fragments in TSS region (TssRatio).
      04_get_insert_size_distribution_per_cell.py   #Get insert size distribution per cell.
      05_Determine_insert_sizes_banding_score.R   #Calculate the periodicity in the frequency of insert sizes and filter out cells based on banding score.
      
#Section 4: Clustering analysis for scATAC-seq.
code: 06_snapatac_twobatch_hg38.R   #Used with SnapATAC(1.0.0).

#Section 5: Constructing peak-cell-matrix.
code: 07_Construct_SparseMatrix.py   #need input files including genuine cell barcode list, peak list and deduplicated bam files.

Optimized cicero: 
input files including peak-cell-matrix and cell grouping information.  Based on Cicero(1.0.15).

#Citation
https://github.com/shendurelab/fly-atac
