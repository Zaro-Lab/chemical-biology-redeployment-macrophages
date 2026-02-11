# Limma for co-culture UE vs SW (Uneaten SW620 vs naive SW620)
library(readxl)
library(limma)

base_dir <- "C:/Users/vishn/Developer/experiments/zaro/volk-ms-1/Volk et al MS Data Analysis for VRT"
out_dir <- file.path(base_dir, "analysis_outputs")

inputs <- list(
  list(file="SurfaceEnrichment_Tables3.xlsx", sheet="5A_CoCulture_SE_Light"),
  list(file="SurfaceEnrichment_Tables3.xlsx", sheet="5B_CoCulture_SE_Heavy")
)

parse_donors <- function(sample_names) {
  donors <- sub("^(?:Area|NArea) (RV\\d+)_.*", "\\1", sample_names)
  donors[!grepl("^RV\\d+$", donors)] <- NA
  return(donors)
}

summary_rows <- list()

for (inp in inputs) {
  file_path <- file.path(base_dir, inp$file)
  sheet <- inp$sheet
  message("Processing UE vs SW: ", inp$file, " / ", sheet)

  df <- read_excel(file_path, sheet = sheet)
  area_cols <- grep("^Area ", names(df), value=TRUE)
  if (length(area_cols) == 0) {
    message("  No Area columns; skipping")
    next
  }

  # Keep only UE and SW columns
  area_cols <- area_cols[grepl("_UE_", area_cols) | grepl("_SW_", area_cols)]
  if (length(area_cols) == 0) {
    message("  No UE/SW columns; skipping")
    next
  }

  mat <- as.matrix(df[, area_cols])
  samples <- data.frame(sample=colnames(mat), stringsAsFactors=FALSE)
  samples$condition <- ifelse(grepl("_UE_", samples$sample), "UE",
                        ifelse(grepl("_SW_", samples$sample), "SW", NA))
  samples$donor <- parse_donors(samples$sample)

  keep_samples <- !is.na(samples$condition)
  samples <- samples[keep_samples, ]
  mat <- mat[, keep_samples, drop=FALSE]

  # log2 transform with LOD imputation
  mat[mat <= 0] <- NA
  min_pos <- suppressWarnings(min(mat[mat > 0], na.rm=TRUE))
  LOD <- ifelse(is.finite(min_pos), min_pos/2, 1)
  mat_imp <- mat
  mat_imp[is.na(mat_imp)] <- LOD
  mat_log2 <- log2(mat_imp)

  # detection filter: >=2 detections in UE before impute
  ue_detect <- rowSums(!is.na(mat[, samples$condition=="UE", drop=FALSE]))
  keep_rows <- ue_detect >= 2
  mat_use <- mat_log2[keep_rows, , drop=FALSE]

  # design
  samples$condition <- factor(samples$condition, levels=c("SW","UE"))
  design <- model.matrix(~0+condition, data=samples)
  colnames(design) <- c("SW","UE")

  corfit <- duplicateCorrelation(mat_use, design, block=samples$donor)
  fit <- lmFit(mat_use, design, block=samples$donor, correlation=corfit$consensus)
  contrast <- makeContrasts(UE-SW, levels=design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)

  tt <- topTable(fit2, number=Inf, sort.by="none")
  out <- cbind(df[keep_rows, c("Protein Group","Accession","Description")], tt)

  out_file <- file.path(out_dir, paste0("limma_", sheet, "_UE_vs_SW_donorblock_LOD_log2.csv"))
  write.csv(out, out_file, row.names=FALSE)

  sig <- sum(out$adj.P.Val < 0.05, na.rm=TRUE)
  up <- sum(out$adj.P.Val < 0.05 & out$logFC > 0, na.rm=TRUE)
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    sheet=sheet, tested=nrow(out), FDR_lt_0_05=sig, up_in_UE=up, LOD=LOD,
    detection_filter="UE>=2"
  )
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(out_dir, "limma_coculture_UE_SW_summary.csv"), row.names=FALSE)
}
