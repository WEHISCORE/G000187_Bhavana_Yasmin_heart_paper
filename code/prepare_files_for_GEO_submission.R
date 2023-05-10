# Prepare (some of) G000187_Bhavana_Yasmin data for GEO submission
# Peter Hickey
# 2023-05-10

library(here)
library(SingleCellExperiment)

outdir <- here("GEO")
dir.create(outdir, recursive = TRUE)

# FASTQs -----------------------------------------------------------------------

# NOTE: It is not feasible to split the FASTQ files into 'heart paper' and
#       'bone, muscle, skin paper' because that is not how the samples were
#       pooled.
#       The samples were pooled into 'Macrophage' (RPI-12; Macrophage) and
#       'Treg' (RPI-6; Treg).

dir.create(file.path(outdir, "FASTQ"))
# NOTE: These plate-level FASTQ files are created by code/scPipe.R.
file.copy(
  from = here("extdata/scPipe/Macrophage/Macrophage.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ", "Macrophage.R1.fastq.gz"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/scPipe/Macrophage/Macrophage.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ", "Macrophage.R2.fastq.gz"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/scPipe/Treg/Treg.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/scPipe/Treg/Treg.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)

# Clean up SCE -----------------------------------------------------------------

dir.create(file.path(outdir, "SCE"))
# NOTE: Starting from the SCE used for the DE analysis.
sce <- readRDS(
  here(
    "data",
    "SCEs",
    "G000187_Bhavana_Yasmin_heart_paper.preprocessed.SCE.rds"))
# NOTE: Revert some of the changes to the rowData.
rownames(sce) <- rowData(sce)$ENSEMBL.GENEID
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))

# Treg 'heart paper' samples ---------------------------------------------------

# NOTE: Bhavana list of Treg samples to be used in the 'heart paper' (email
#       2023-02-24).
#       We can't simply rely on the sample metadata (e.g.,`sce$Tissue == Heart`)
#       because the sample metadata has been such a mess from the beginning of
#       this project.
#       I've done my best to convert Bhavana's instructions into something that
#       allows me to filter down to just the relevant samples.
treg_heart_samples <- data.frame(
  sample_descriptor_tissue = c(
    # 'Spleen'
    "Sp 3,4", "Sp 5,6", "Sp 7,8", "Sp 9,10",
    # 'Heart'
    "HT 2", "HT 3", "HT 4",
    # 'Heart Lymph node'
    "LN 1", "LN 2", "LN 3",
    # 'Heart Spleen')
    "Spl 1", "Spl 2", "Spl 3"),
  renamed_to = c(
    "D0 SPL1", "D0 SPL2", "D0 SPL3", "D0 SPL4",
    "HRT1", "HRT2", "HRT3",
    "HRT LN1", "HRT LN2", "HRT LN3",
    "HRT SPL1", "HRT SPL2", "HRT SPL3"))
# Check for typos.
stopifnot(
  all(
    treg_heart_samples[["sample_descriptor_tissue"]] %in%
      sce$sample_descriptor_tissue))

dir.create(file.path(outdir, "SCE", "Treg_heart_paper"))
sce_treg_heart <- sce[
  ,
  sce$sample_descriptor_tissue %in%
    treg_heart_samples[["sample_descriptor_tissue"]]]

# Gene counts (UMIs)
write.csv(
  x = as.data.frame(as.matrix(assay(sce_treg_heart, "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Treg_heart_paper",
      "Treg_heart_paper.gene_UMI_counts.csv.gz")),
  row.names = TRUE)
# ERCC counts (UMIs)
write.csv(
  x = as.data.frame(
    as.matrix(assay(altExp(sce_treg_heart, "ERCC"), "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Treg_heart_paper",
      "Treg_heart_paper.ERCC_UMI_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce_treg_heart)),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Treg_heart_paper",
      "Treg_heart_paper.sample_sheet.csv.gz")),
  row.names = TRUE)

# NOTE: Pseudobulk UMI counts because that's what Bhavana says they analysed
#       (email 2023-02-24).
library(scuttle)
se_treg_heart <- applySCE(
  sce_treg_heart,
  aggregateAcrossCells,
  ids = sce_treg_heart$sample_descriptor_tissue,
  use.assay.type = "UMI_counts")
# NOTE: Renaming samples as requested by Bhavana (2023-02-24).
colnames(se_treg_heart) <- treg_heart_samples[["renamed_to"]][
  match(colnames(se_treg_heart), treg_heart_samples[["sample_descriptor_tissue"]])]
# NOTE: Drop colData because it's a mess. The only relevant sample metadata is
#       already available in the 'new' colnames.
colData(se_treg_heart) <- NULL

# Gene counts (UMIs)
write.csv(
  x = as.data.frame(as.matrix(assay(se_treg_heart, "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Treg_heart_paper",
      "Treg_heart_paper.pseudobulked_gene_UMI_counts.csv.gz")),
  row.names = TRUE)
# ERCC counts (UMIs)
write.csv(
  x = as.data.frame(
    as.matrix(assay(altExp(se_treg_heart, "ERCC"), "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Treg_heart_paper",
      "Treg_heart_paper.pseudobulked_ERCC_UMI_counts.csv.gz")),
  row.names = TRUE)

# Macrophage 'heart paper' -----------------------------------------------------

# NOTE: Bhavana list of Macrophage samples to be used in the 'heart paper'
#       (email 2023-03-14; although she accidentally swapped the 'heart paper'
#       and 'bone, muscle, skin paper' samples, which is corrected here).
#       We can't simply rely on the sample metadata (e.g.,`sce$Tissue == Heart`)
#       because the sample metadata has been such a mess from the beginning of
#       this project.
#       I've done my best to convert Bhavana's instructions into something that
#       allows me to filter down to just the relevant samples.

macrophage_heart_samples <- c(
  "HRT D4 P1", "HRT D4 P2", "HRT D4 P3", "HRT D4 T1", "HRT D4 T2", "HRT D4 T3",
  "HRT D7 P4", "HRT D7 P5", "HRT D7 P6", "HRT D7 T1", "HRT D7 T3", "HRT D7 T4")
stopifnot(
  all(macrophage_heart_samples %in% sce$sample_descriptor_tissue))

dir.create(file.path(outdir, "SCE", "Macrophage_heart_paper"))
sce_macrophage_heart <- sce[
  ,
  sce$sample_descriptor_tissue %in% macrophage_heart_samples]

# NOTE: Files uploaded to GEO still use `LCE725` rather than `Macrophage` in
#       the column names and as `plate_number`.
colnames(sce_macrophage_heart) <- sub(
  "Macrophage",
  "LCE725",
  colnames(sce_macrophage_heart))
sce_macrophage_heart$plate_number <- "LCE725"

# Gene counts (UMIs)
write.csv(
  x = as.data.frame(
    as.matrix(assay(sce_macrophage_heart, "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Macrophage_heart_paper",
      "Macrophage_heart_paper.gene_UMI_counts.csv.gz")),
  row.names = TRUE)
# ERCC counts (UMIs)
write.csv(
  x = as.data.frame(
    as.matrix(assay(altExp(sce_macrophage_heart, "ERCC"), "UMI_counts"))),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Macrophage_heart_paper",
      "Macrophage_heart_paper.ERCC_UMI_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce_macrophage_heart)),
  file = gzfile(
    file.path(
      outdir,
      "SCE",
      "Macrophage_heart_paper",
      "Macrophage_heart_paper.sample_sheet.csv.gz")),
  row.names = TRUE)
