# Preprocess the 'heart paper' samples from G000187_Bhavana_Yasmin.
# The aim is to replicate they key steps of the 'preprocessing.Rmd' for just
# the 'heart paper' samples.
# Peter Hickey
# 2023-05-10

library(SingleCellExperiment)
library(here)

# Setting up the data ----------------------------------------------------------

sce_deduped <- readRDS(
  here(
    "data",
    "SCEs",
    "G000187_Bhavana_Yasmin_heart_paper.UMI_deduped.scPipe.SCE.rds"))
sce_not_deduped <- readRDS(
  here(
    "data",
    "SCEs",
    "G000187_Bhavana_Yasmin_heart_paper.not_UMI_deduped.scPipe.SCE.rds"))

# Drop non-'heart paper' samples.
sce_deduped <- sce_deduped[, !is.na(sce_deduped$sample_descriptor_tissue)]
sce_not_deduped <- sce_not_deduped[
  , !is.na(sce_not_deduped$sample_descriptor_tissue)]

stopifnot(
  identical(
    sort(rownames(altExp(sce_deduped, "ERCC"))),
    sort(rownames(altExp(sce_not_deduped, "ERCC")))))
altExp(sce_deduped, "ERCC") <-
  altExp(sce_deduped, "ERCC")[rownames(altExp(sce_not_deduped), "ERCC"), ]
colnames(sce_not_deduped) <- sub(
  "\\.not_UMI_deduped",
  "",
  colnames(sce_not_deduped))
colnames(altExp(sce_not_deduped, "ERCC")) <- sub(
  "\\.not_UMI_deduped",
  "",
  colnames(altExp(sce_not_deduped, "ERCC")))
stopifnot(
  identical(rownames(sce_deduped), rownames(sce_not_deduped)),
  identical(colnames(sce_deduped), colnames(sce_not_deduped)),
  identical(
    rownames(altExp(sce_deduped, "ERCC")),
    rownames(altExp(sce_not_deduped, "ERCC"))),
  identical(
    colnames(altExp(sce_deduped, "ERCC")),
    colnames(altExp(sce_not_deduped, "ERCC"))))

# Combine UMI and read counts a single SCE.
sce <- sce_deduped
assay(sce, "UMI_counts") <- assay(sce_deduped, "counts")
assay(sce, "read_counts") <- assay(sce_not_deduped, "counts")
assay(altExp(sce, "ERCC"), "UMI_counts") <-
  assay(altExp(sce_deduped, "ERCC"), "counts")
assay(altExp(sce, "ERCC"), "read_counts") <-
  assay(altExp(sce_not_deduped, "ERCC"), "counts")
# NOTE: Nullify some now-unrequired data.
assay(sce, "counts") <- NULL
assay(altExp(sce, "ERCC"), "counts") <- NULL
sce$UMI_deduped <- NULL

# Initial tidying of sample metadata
sce$experiment <- factor(
  ifelse(sce$plate_number == "Treg", "Treg", "Macrophage"),
  levels = c("Macrophage", "Treg"))
sce$well_position <- factor(
  sce$well_position,
  unlist(lapply(LETTERS[1:16], function(x) paste0(x, 1:24))))
# NOTE: This info is not recorded for all samples.
sce$original_plate_well <- factor(
  sce$original_plate_well,
  c(unlist(lapply(LETTERS[1:8], function(x) paste0(x, 1:12)))))

# Load and attach sample metadata
# NOTE: NAs are being parsed as the character vector "NA"
library(readxl)
sample_metadata_df <- read_excel(
  here(
    "data",
    "sample_sheets",
    "G000187_Bhavana_Yasmin_heart_paper.sample_metadata.xlsx"))

# NOTE: Bhavana modified the `Sample ID` value for some rows. I have now
#       corrected the XLSX file on google sheets and GitHub so the
#       modifications are not needed.
stopifnot(
  all(sce$sample_descriptor_tissue %in% sample_metadata_df$`Sample ID`) &
    all(sample_metadata_df$`Sample ID` %in% sce$sample_descriptor_tissue))

colData(sce) <- DataFrame(
  dplyr::left_join(
    as.data.frame(colData(sce)),
    sample_metadata_df,
    by = c(
      "sample_descriptor_tissue" = "Sample ID",
      "experiment" = "Experiment")),
  row.names = colnames(sce),
  check.names = FALSE)

# NOTE: `Condition` field is redundant for `Treg` experiment. Re-name to
#       `Treatment`(?) or similar and set to NA for Treg samples.
sce$Condition[sce$experiment == "Treg"] <- "NA"

sce$Tissue <- factor(
  sce$Tissue,
  levels = c("Bone", "Heart", "Lymph node", "Muscle", "Skin", "Spleen"))

# NOTE: `Day` is recorded as a character rather than integer in the data
#       uploaded to GEO.
sce$Day <- as.character(sce$Day)

sce$group <- interaction(
  sce$Tissue,
  sce$Day,
  sce$Condition,
  drop = TRUE,
  lex.order = TRUE)

# Re-order rows to match that of the original data processing.
rn <- readLines(here("data/original_row_order.txt"))
stopifnot(
  identical(setdiff(rn, rownames(sce)), character(0)),
  identical(setdiff(rownames(sce), rn), character(0)))
aern <- readLines(here("data/original_row_order.ERCC.txt"))
stopifnot(
  identical(setdiff(aern, rownames(altExp(sce))), character(0)),
  identical(setdiff(rownames(altExp(sce)), aern), character(0)))
sce <- sce[rn, ]
altExp(sce) <- altExp(sce)[aern, ]

# Incorporating gene-based annotation ------------------------------------------

# Extract rownames (Ensembl IDs) to use as key in database lookups.
ensembl <- rownames(sce)

# Pull out useful gene-based annotations from the Ensembl-based database.
library(EnsDb.Mmusculus.v79)
library(ensembldb)
# NOTE: These columns were customised for this project.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(EnsDb.Mmusculus.v79)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = EnsDb.Mmusculus.v79,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
# NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
ensdb_df$ENSEMBL.GENEID <- ensembl

# NOTE: Mus.musculus combines org.Mm.eg.db and
#       TxDb.Mmusculus.UCSC.mm10.knownGene (as well as others) and therefore
#       uses entrez gene and RefSeq based data.
library(Mus.musculus)
# NOTE: These columns were customised for this project.
ncbi_columns <- c(
  # From TxDB: None required
  # From OrgDB
  "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Mus.musculus)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Mus.musculus,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)

rowData(sce) <- cbind(ensdb_df, ncbi_df)

# Replace the row names of the SCE by the gene symbols (where available).
library(scuttle)
rownames(sce) <- uniquifyFeatureNames(
  ID = rownames(sce),
  # NOTE: An Ensembl ID may map to 0, 1, 2, 3, ... gene symbols.
  #       When there are multiple matches only the 1st match is used.
  names = sapply(rowData(sce)$ENSEMBL.SYMBOL, function(x) {
    if (length(x)) {
      x[[1]]
    } else {
      NA_character_
    }
  }))

# Have we sequenced enough? ----------------------------------------------------

saturation <- 1 -
  (colSums(assay(sce, "UMI_counts")) / colSums(assay(sce, "read_counts")))
sce$saturation <- saturation

# Concluding remarks -----------------------------------------------------------

saveRDS(
  sce,
  here(
    "data",
    "SCEs",
    "G000187_Bhavana_Yasmin_heart_paper.preprocessed.SCE.rds"),
  compress = "xz")

# Degust outputs ---------------------------------------------------------------

# Create pseudobulked SEs
pb_se <- aggregateAcrossCells(
  sce,
  ids = DataFrame(group = sce$group, Mouse = sce$Mouse),
  use.assay.type = assayNames(sce),
  use.altexps = NULL)

SCE2DegustCSV <- function(sce, assay.type) {
  genes_df <- data.frame(
    `Gene ID` = rowData(sce)$ENSEMBL.GENEID,
    name = rownames(sce),
    check.names = FALSE)
  counts_df <- as.data.frame(as.matrix(assay(sce, assay.type)))
  colnames(counts_df) <- sce$sample_descriptor_tissue
  cbind(genes_df, counts_df)
}

degust_dir <- here("output", "Degust")
dir.create(degust_dir)

# Write Macrophage CSVs
write.csv(
  SCE2DegustCSV(sce[, sce$experiment == "Macrophage"], "read_counts"),
  here(degust_dir, "Macrophage.heart_paper.read_counts.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(sce[, sce$experiment == "Macrophage"], "UMI_counts"),
  here(degust_dir, "Macrophage.heart_paper.UMI_counts.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(pb_se[, pb_se$experiment == "Macrophage"], "read_counts"),
  here(degust_dir, "Macrophage.heart_paper.read_counts.pseudobulked.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(pb_se[, pb_se$experiment == "Macrophage"], "UMI_counts"),
  here(degust_dir, "Macrophage.heart_paper.UMI_counts.pseudobulked.csv"),
  quote = TRUE,
  row.names = FALSE)

# Write Treg CSVs
write.csv(
  SCE2DegustCSV(sce[, sce$experiment == "Treg"], "read_counts"),
  here(degust_dir, "Treg.heart_paper.read_counts.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(sce[, sce$experiment == "Treg"], "UMI_counts"),
  here(degust_dir, "Treg.heart_paper.UMI_counts.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(pb_se[, pb_se$experiment == "Treg"], "read_counts"),
  here(degust_dir, "Treg.heart_paper.read_counts.pseudobulked.csv"),
  quote = TRUE,
  row.names = FALSE)
write.csv(
  SCE2DegustCSV(pb_se[, pb_se$experiment == "Treg"], "UMI_counts"),
  here(degust_dir, "Treg.heart_paper.UMI_counts.pseudobulked.csv"),
  quote = TRUE,
  row.names = FALSE)
