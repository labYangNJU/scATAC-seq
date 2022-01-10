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


#Section 4: Optimized cicero calling and calculating gene acivities: 

code: 

      04_Optimized_cicero_calling.R   #Input files including peak-cell-matrix and cell grouping information.  Based on Cicero(1.0.15).

# Citation

# Reference
https://github.com/shendurelab/fly-atac

https://github.com/shendurelab/mouse-atac
