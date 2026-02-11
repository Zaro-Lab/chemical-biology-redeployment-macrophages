# Run limma across all relevant sheets using the same donor-blocked + LOD impute log2 approach
library(readxl)
library(limma)

base_dir <- "C:/Users/vishn/Developer/experiments/zaro/volk-ms-1/Volk et al MS Data Analysis for VRT"
out_dir <- file.path(base_dir, "analysis_outputs")

inputs <- list(
  list(file="SurfaceEnrichment_Tables2.xlsx", sheet="4A_EatersNon_SE_Light", conds=c("E","NE")),
  list(file="SurfaceEnrichment_Tables2.xlsx", sheet="4B_EatersNon_SE_Heavy", conds=c("E","NE")),
  list(file="SurfaceEnrichment_Tables3.xlsx", sheet="5A_CoCulture_SE_Light", conds=c("SW","Mac")),
  list(file="SurfaceEnrichment_Tables3.xlsx", sheet="5B_CoCulture_SE_Heavy", conds=c("SW","Mac")),
  list(file="WholeProtein_Tables1.xlsx", sheet="3A_EatersNon_WP_Light", conds=c("E","NE")),
  list(file="WholeProtein_Tables1.xlsx", sheet="3B_EatersNon_WP_Heavy", conds=c("E","NE")),
  list(file="WholeProtein_Tables1.xlsx", sheet="3C_EatersNon_WP_Light+Heavy", conds=c("E","NE"))
)

parse_conditions <- function(sample_names) {
  cond <- ifelse(grepl("_E_", sample_names), "E",
           ifelse(grepl("_NE_", sample_names), "NE",
           ifelse(grepl("_SW_", sample_names), "SW",
           ifelse(grepl("_Mac_", sample_names), "Mac", NA))))
  return(cond)
}

parse_donors <- function(sample_names) {
  donors <- sub("^(?:Area|NArea) (RV\\d+)_.*", "\\1", sample_names)
  donors[!grepl("^RV\\d+$", donors)] <- NA
  # fallback: try anywhere
  miss <- is.na(donors)
  if (any(miss)) {
    donors[miss] <- sub(".*(RV\\d+).*", "\\1", sample_names[miss])
    donors[!grepl("^RV\\d+$", donors)] <- NA
  }
  return(donors)
}

safe_write <- function(df, path) {
  write.csv(df, path, row.names=FALSE)
}

summary_rows <- list()

for (inp in inputs) {
  file_path <- file.path(base_dir, inp$file)
  sheet <- inp$sheet
  message("Processing: ", inp$file, " / ", sheet)

  df <- read_excel(file_path, sheet = sheet)
  # choose Area or NArea columns
  area_cols <- grep("^Area ", names(df), value=TRUE)
  if (length(area_cols) == 0) {
    area_cols <- grep("^NArea ", names(df), value=TRUE)
  }
  if (length(area_cols) == 0) {
    message("  No Area/NArea columns; skipping")
    next
  }
  # filter columns to requested condition set (e.g., E/NE or SW/Mac)
  if (!is.null(inp$conds)) {
    cond_patterns <- paste0("_", inp$conds, "_")
    area_cols <- area_cols[sapply(area_cols, function(cn) any(sapply(cond_patterns, function(p) grepl(p, cn))))]
  }
  if (length(area_cols) == 0) {
    message("  No columns match requested conditions; skipping")
    next
  }

  mat <- as.matrix(df[, area_cols])

  samples <- data.frame(sample=colnames(mat), stringsAsFactors=FALSE)
  samples$condition <- parse_conditions(samples$sample)
  samples$donor <- parse_donors(samples$sample)

  # keep only samples with recognized conditions
  keep_samples <- !is.na(samples$condition)
  if (!all(keep_samples)) {
    samples <- samples[keep_samples, ]
    mat <- mat[, keep_samples, drop=FALSE]
  }

  conds <- unique(samples$condition)
  if (length(conds) != 2) {
    message("  Not exactly 2 conditions (", paste(conds, collapse=","), "); skipping")
    next
  }

  # set baseline/case order
  if (all(c("NE","E") %in% conds)) {
    baseline <- "NE"; case <- "E"
  } else if (all(c("Mac","SW") %in% conds)) {
    baseline <- "Mac"; case <- "SW"
  } else {
    # fallback: alphabetical
    baseline <- sort(conds)[1]
    case <- sort(conds)[2]
  }

  # log2 transform with LOD imputation
  mat[mat <= 0] <- NA
  min_pos <- suppressWarnings(min(mat[mat > 0], na.rm=TRUE))
  LOD <- ifelse(is.finite(min_pos), min_pos/2, 1)
  mat_imp <- mat
  mat_imp[is.na(mat_imp)] <- LOD
  mat_log2 <- log2(mat_imp)

  # detection filter: >=2 detections in case group before impute
  case_idx <- samples$condition == case
  det <- rowSums(!is.na(mat[, case_idx, drop=FALSE]))
  keep_rows <- det >= 2
  mat_use <- mat_log2[keep_rows, , drop=FALSE]

  # design
  samples$condition <- factor(samples$condition, levels=c(baseline, case))
  design <- model.matrix(~0+condition, data=samples)
  colnames(design) <- c(baseline, case)

  # donor blocking if any donor has >1 sample
  use_block <- any(table(samples$donor) > 1, na.rm=TRUE)
  if (use_block) {
    corfit <- duplicateCorrelation(mat_use, design, block=samples$donor)
    fit <- lmFit(mat_use, design, block=samples$donor, correlation=corfit$consensus)
  } else {
    fit <- lmFit(mat_use, design)
  }
  contrast <- makeContrasts(contrasts=paste0(case, "-", baseline), levels=design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)

  tt <- topTable(fit2, number=Inf, sort.by="none")
  out <- cbind(df[keep_rows, c("Protein Group","Accession","Description")], tt)

  out_file <- file.path(out_dir, paste0("limma_", sheet, "_donorblock_LOD_log2.csv"))
  safe_write(out, out_file)

  sig <- sum(out$adj.P.Val < 0.05, na.rm=TRUE)
  up <- sum(out$adj.P.Val < 0.05 & out$logFC > 0, na.rm=TRUE)
  summary_rows[[length(summary_rows)+1]] <- data.frame(
    sheet=sheet, baseline=baseline, case=case,
    tested=nrow(out), FDR_lt_0_05=sig, up_in_case=up,
    LOD=LOD, detection_filter=paste0("det>=2 in ", case),
    blocked=use_block
  )
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(out_dir, "limma_all_summary.csv"), row.names=FALSE)
}
