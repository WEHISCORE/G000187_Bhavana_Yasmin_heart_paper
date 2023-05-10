# G000187_Bhavana_Yasmin_heart_paper

## Context

The `G000187_Bhavana_Yasmin` project contains data from 2 papers, the 'heart paper' and another, which were originally processed together as one.
Bhavana and Yasmin requested that only the 'heart paper' data be made public at this time, so this repository contains code that reproduces the data used in 'heart paper' and masks or omits the phenotypes of the other samples in the broader dataset.

## Data processing

Code used to process mini-bulk RNA-seq data for `G000187_Bhavana_Yasmin` 'heart paper' using [**scPipe**](	https://bioconductor.org/packages/scPipe/); see [`code/scPipe.R`](code/scPipe.R).
In brief:

1. FASTQ files are reformatted with `scPipe::sc_trim_barcode()` so that the barcode and UMI sequences are moved from the sequence of R1 into the read name.
2. Reads are aligned to the mouse reference genome (GRCm38.p6) using `Rsubread::align()`.
3. Aligned reads are assigned to the annotated exons of the GENCODE comprehensive gene annotation on the primary assembly ([`gencode.vM25.primary_assembly.annotation.gff3`](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gff3.gz)) with `scPipe::sc_exon_mapping()`.
4. Aligned reads that are mapped to annotated exons are demultiplexed by the cell barcode using `scPipe::sc_demultiplex()`. Only the first 7 bp of the 8 bp cell barcode are used because the quality of the last base in the cell barcode is often poor.
5. Generate the gene count matrices (UMI-deduplicated and non-deduplicated versions)
  a. UMI-deduplicated: Use `scPipe::sc_gene_counting()` to perform simple correction of the UMIs to merge UMIs with and edit distance of 1.
  b. Non-deduplicated: Use a custom function, `geneCountingNoUMIDedup()`, to create a gene count matrix while ignoring the UMI data (i.e. no deduplication based on UMIs).
6. Create [*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/) objects of each count matrix and incorporating the sample metadata.
  a. UMI-deduplicated: [`data/SCEs/G000187_Bhavana_Yasmin_heart_paper.UMI_deduped.scPipe.SCE.rds`](data/SCEs/G000187_Bhavana_Yasmin_heart_paper.UMI_deduped.scPipe.SCE.rds).
  b. Non-deduplicated: [`data/SCEs/G000187_Bhavana_Yasmin_heart_paper.not_UMI_deduped.scPipe.SCE.rds`](data/SCEs/G000187_Bhavana_Yasmin_heart_paper.not_UMI_deduped.scPipe.SCE.rds).

The [*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/) objects then undergo some additional processing; see [`code/G000187_Bhavana_Yasmin_heart_paper.preprocess.R`](code/G000187_Bhavana_Yasmin_heart_paper.preprocess.R).
In brief:

1. Dropping non-'heart paper' samples.
2. Tidying up of sample metadata.
3. Combining UMI-deduplicated and non-deduplicated data into one [*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/) object.
4. Re-ordering of genes to match the order in the original `G000187_Bhavana_Yasmin` project.
5. Incorporating gene-based annotation.
6. Calculation of some QC metrics.
7. Saving the processed [*SingleCellExperiment*](https://bioconductor.org/packages/SingleCellExperiment/) object as [`data/SCEs/G000187_Bhavana_Yasmin_heart_paper.preprocessed.SCE.rds`](data/SCEs/G000187_Bhavana_Yasmin_heart_paper.preprocessed.SCE.rds).
8. Generating various CSV files of the raw and pseudobulked UMI and gene (i.e. non-deduplicated) count matrices that are suitable for use with [degust](https://degust.erc.monash.edu/), as requested by Bhanava and Yasmin.
