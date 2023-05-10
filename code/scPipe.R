# Process 'heart paper' data from G000187_Bhavana_Yasmin with scPipe
# Peter Hickey
# 2023-05-09

# Setup ------------------------------------------------------------------------

library(here)
library(scPipe)
library(Rsubread)

# Load sample sheet ------------------------------------------------------------

sample_sheet <- as(
    read.csv(
      here("data/sample_sheets/S000306_S000326.heart_paper.sample_sheet.csv"),
      row.names = 1),
  "DataFrame")

# Key variables ----------------------------------------------------------------

plates <- unique(sample_sheet$plate_number)
names(plates) <- plates
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$plate_number,
  unique)
outdir <- here("data", "SCEs")
dir.create(outdir, recursive = TRUE)
extdir <- file.path(
  "/vast/scratch/users/hickey/G000187_Bhavana_Yasmin_heart_paper/scPipe",
  plates)
names(extdir) <- plates
sapply(extdir, dir.create, recursive = TRUE)
# NOTE: Only using first 7 nt of barcode.
read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
# NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
#       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
organism <- "mmusculus_gene_ensembl"
# NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
#       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
gene_id_type <- "ensembl_gene_id"

# Input files ------------------------------------------------------------------

# FASTQ files
# NOTE: Only 1 pair of FASTQ files are relevant (the others are from 'empty'
#       barcodes).
r1_fq <- list.files(
  here("extdata", "FASTQ"),
  pattern = glob2rx("*R1*fastq.gz"),
  full.names = TRUE)
r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))
tx_fq <- file.path(extdir, paste0(plates, ".R2.fastq.gz"))
names(tx_fq) <- plates
barcode_fq <- gsub("R2", "R1", tx_fq)

lapply(plates, function(plate) {
  message(plate)
  cmd <- paste0(
    "cat ",
    paste(grep(plate, r1_fq, value = TRUE), collapse = " "),
    " > ",
    barcode_fq[[plate]],
    "\n",
    "cat ",
    paste(grep(plate, r2_fq, value = TRUE), collapse = " "),
    " > ",
    tx_fq[[plate]])
  system(cmd)
})

# Genome index
# NOTE: These samples didn't include ERCCs, but it doesn't hurt to have them
#       in the reference genome.
genome_index <- here("extdata", "GRCm38.p6", "GRCm38_with_ERCC")

# Genome annotation(s)
annofn <- c(
  here("extdata", "GRCm38.p6", "gencode.vM25.primary_assembly.annotation.gff3"),
  system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))

# Cell barcodes
bc_anno <- file.path(extdir, paste0(plates, ".barcode_annotation.csv"))
names(bc_anno) <- plates

for (plate in plates) {
  message(plate)
  tmp <- sample_sheet[sample_sheet$plate_number == plate, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    # NOTE: For some reason the primer name and sequence columns have been
    #       reversed in this sample sheet.
    # NOTE: Only using first 7 nt of barcode.
    barcode = strtrim(tmp$c_rt1_primer_name, 7))
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[plate]],
    quote = FALSE,
    row.names = FALSE)
}

# Output files -----------------------------------------------------------------

combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
names(combined_fq) <- names(tx_fq)
subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
exon_bam <- gsub("subread", "exon", subread_bam)

# FASTQ reformatting -----------------------------------------------------------

filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
lapply(seq_along(tx_fq), function(i) {
  message(combined_fq[i])
  sc_trim_barcode(
    outfq = combined_fq[i],
    r1 = tx_fq[i],
    r2 = barcode_fq[i],
    read_structure = read_structure,
    filter_settings = filter_settings)
})

# Aligning reads to a reference genome -----------------------------------------

align(
  index = genome_index,
  readfile1 = combined_fq,
  output_file = subread_bam,
  nthreads = 12)

# Assigning reads to annotated exons -------------------------------------------

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
lapply(seq_along(subread_bam), function(i) {
  message(i)
  sc_exon_mapping(
    inbam = subread_bam[i],
    outbam = exon_bam[i],
    annofn = annofn,
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
})

# De-multiplexing data ---------------------------------------------------------

max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
lapply(seq_along(exon_bam), function(i) {
  message(i)
  sc_demultiplex(
    inbam = exon_bam[i],
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)
})

# Gene counting deduped data ---------------------------------------------------

UMI_cor <- 1
gene_fl <- FALSE
lapply(seq_along(bc_anno), function(i) {
  message(i)
  sc_gene_counting(
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
})

# Create and save deduped SingleCellExperiment ---------------------------------

list_of_sce <- lapply(plates, function(plate) {
  message(plate)
  create_sce_by_dir(
    datadir = extdir[[plate]],
    organism = organism,
    gene_id_type = gene_id_type,
    pheno_data = sample_sheet[sample_sheet$plate_number == plate, ],
    # NOTE: Create the report separately for more fine-grained control.
    report = FALSE)
})
sce <- do.call(combineCols, c(unname(list_of_sce), fill = 0L, delayed = FALSE))
sce$UMI_deduped <- TRUE
# NOTE: Need to tidy up metadata to support use with scPipe::QC_metrics() on
#       the returned object. This assumes that all elements of `list_of_sce`
#       have the same metadata.
metadata(sce) <- list(
  scPipe = list(
    version = metadata(sce)[[1]][["version"]],
    QC_cols = metadata(sce)[[1]][["QC_cols"]],
    demultiplex_info = data.frame(
      status = metadata(sce)[[1]][["demultiplex_info"]][["status"]],
      count = rowSums(
        do.call(
          cbind,
          lapply(
            metadata(sce)[seq(1, length(metadata(sce)), 2)],
            function(x) {
              x[["demultiplex_info"]][["count"]]
            })))),
    UMI_dup_info = data.frame(
      duplication.number =
        metadata(sce)[[1]][["UMI_dup_info"]][["duplication.number"]],
      count = rowSums(
        do.call(
          cbind,
          lapply(
            metadata(sce)[seq(1, length(metadata(sce)), 2)],
            function(x) {
              x[["UMI_dup_info"]][["count"]]
            }))))),
  Biomart = metadata(sce)[[2]])

assay(sce, withDimnames = FALSE) <- as(
  assay(sce, withDimnames = FALSE),
  "dgCMatrix")
sce <- splitAltExps(
  sce,
  ifelse(grepl("^ERCC", rownames(sce)), "ERCC", "Endogenous"))

saveRDS(
  sce,
  file.path(
    outdir,
    "G000187_Bhavana_Yasmin_heart_paper.UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")

# Create QC report of deduped data----------------------------------------------

library(readr)
library(plotly)
library(DT)
library(scran)
library(Rtsne)
dir.create(here("output", "scPipe"), recursive = TRUE)
# NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
lapply(plates, function(plate) {
  message(plate)
  try(
    create_report(
      sample_name = plate,
      outdir = extdir[[plate]],
      r1 = tx_fq[[plate]],
      r2 = barcode_fq[[plate]],
      outfq = combined_fq[[plate]],
      read_structure = read_structure,
      filter_settings = filter_settings,
      align_bam = subread_bam[[plate]],
      genome_index = genome_index,
      map_bam = exon_bam[[plate]],
      exon_anno = annofn,
      stnd = stnd,
      fix_chr = fix_chr,
      barcode_anno = bc_anno[[plate]],
      max_mis = max_mis,
      UMI_cor = UMI_cor,
      gene_fl = gene_fl,
      organism = organism,
      gene_id_type = gene_id_type))

  # NOTE: Workaround bug in create_report() and stop output after 'Data summary'
  #       section.
  tmp <- readLines(file.path(extdir[[plate]], "report.Rmd"))
  # NOTE: Hotfix for https://github.com/LuyiTian/scPipe/issues/146
  tmp[109] <- "The organism is \"`r getVal(params$organism, \"unknown\")`\", and gene id type is \"`r getVal(params$gene_id_type,"
  # NOTE: Another hot-fix another long-standing bug in the scPipe QC report
  #       generation functionality.
  tmp <- c(tmp[1:168], "knitr::knit_exit()", tmp[169:length(tmp)])
  writeLines(tmp, file.path(extdir[[plate]], "report.Rmd"))
  knitr::wrap_rmd(
    file = file.path(extdir[[plate]], "report.Rmd"),
    width = 120,
    backup = NULL)
  rmarkdown::render(
    input = file.path(extdir[[plate]], "report.Rmd"),
    output_file = file.path(extdir[[plate]], "report.html"),
    knit_root_dir = ".")

  # NOTE: Copy the QC report to the repository.
  file.copy(
    from = file.path(extdir[[plate]], "report.nb.html"),
    to = here(
      "output",
      "scPipe",
      paste0(plate, ".scPipe_QC_report.nb.html")),
    overwrite = TRUE)
})

# Gene counting non-deduped data -----------------------------------------------

# NOTE: Need to use my own gene counting function because not using UMI
#       deduplication.
geneCountingNoUMIDedup <- function(outdir, bc_anno) {
  files <- list.files(file.path(outdir, "count"), full.names = TRUE)
  names(files) <- sub("\\.csv", "", basename(files))
  counts <- lapply(files, function(file) {
    message(basename(file))
    data.table::fread(file, select = 1)[, table(gene_id)]
  })
  genes <- Reduce(union, lapply(counts, names))
  x <- matrix(
    0L,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(genes, names(counts)))
  for (j in names(counts)) {
    xx <- counts[[j]]
    x[names(xx), j] <- xx
  }
  z <- cbind(
    data.frame(gene_id = rownames(x)),
    as.data.frame(x))
  data.table::fwrite(
    x = z,
    file = file.path(paste0(outdir, "_no_dedup"), "gene_count.csv"),
    row.names = FALSE,
    nThread = 1)
}

sapply(paste0(extdir, "_no_dedup"), dir.create)
lapply(names(bc_anno), function(n) {
  message(n)
  geneCountingNoUMIDedup(
    outdir = extdir[[n]],
    bc_anno = bc_anno[[n]])
})

# Create and save non-deduped SingleCellExperiment -----------------------------

list_of_no_dedup_sce <- lapply(plates, function(plate) {
  x <- data.table::fread(
    file.path(paste0(extdir[[plate]], "_no_dedup"), "gene_count.csv"))
  counts <- as.matrix(x[, -1])
  sce <- SingleCellExperiment(
    list(counts = counts),
    colData = sample_sheet[colnames(counts), ])
  rownames(sce) <- x$gene_id
  sce
})

no_dedup_sce <- do.call(
  combineCols,
  c(unname(list_of_no_dedup_sce), fill = 0L, delayed = FALSE))
no_dedup_sce$UMI_deduped <- FALSE
colnames(no_dedup_sce) <- paste0(colnames(no_dedup_sce), ".not_UMI_deduped")
assay(no_dedup_sce, withDimnames = FALSE) <- unname(
  as(assay(no_dedup_sce, withDimnames = FALSE), "dgCMatrix"))
no_dedup_sce <- splitAltExps(
  no_dedup_sce,
  ifelse(grepl("^ERCC", rownames(no_dedup_sce)), "ERCC", "Endogenous"))
# NOTE: Order genes and samples as in `sce`.
no_dedup_sce <- no_dedup_sce[rownames(sce),
                             paste0(colnames(sce), ".not_UMI_deduped")]

saveRDS(
  no_dedup_sce,
  file.path(
    outdir,
    "G000187_Bhavana_Yasmin_heart_paper.not_UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")

# Copy outputs -----------------------------------------------------------------

# NOTE: This is very fragile and will likely require modification for each
#       project.
cmd <- paste0(
  "rsync --verbose --human-readable --recursive --progress --archive \ ",
  unique(dirname(extdir)), "\ ",
  here("extdata"))
system(cmd)
