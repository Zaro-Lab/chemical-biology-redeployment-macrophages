library(readxl)
library(limma)

file <- "C:/Users/vishn/Developer/experiments/zaro/volk-ms-1/Volk et al MS Data Analysis for VRT/SurfaceEnrichment_Tables2.xlsx"

df <- read_excel(file, sheet = "4B_EatersNon_SE_Heavy")
cols <- grep("^Area", names(df), value=TRUE)
mat <- as.matrix(df[, cols])

# sample info
samples <- data.frame(sample=colnames(mat))
samples$condition <- ifelse(grepl("_E_", samples$sample), "E", "NE")
samples$donor <- sub("Area (RV\\d+)_.*", "\\1", samples$sample)

# log2 transform with LOD imputation
mat[mat <= 0] <- NA
min_pos <- suppressWarnings(min(mat[mat > 0], na.rm=TRUE))
LOD <- ifelse(is.finite(min_pos), min_pos/2, 1)
mat_imp <- mat
mat_imp[is.na(mat_imp)] <- LOD
mat_log2 <- log2(mat_imp)

# detection filter: >=2 E replicates detected (non-NA before impute)
e_detect <- rowSums(!is.na(mat[, samples$condition=="E", drop=FALSE]))
keep <- e_detect >= 2

mat_use <- mat_log2[keep, , drop=FALSE]

# design
samples$condition <- factor(samples$condition, levels=c("NE","E"))
design <- model.matrix(~0+condition, data=samples)
colnames(design) <- c("NE","E")

# duplicateCorrelation for donor
corfit <- duplicateCorrelation(mat_use, design, block=samples$donor)
fit <- lmFit(mat_use, design, block=samples$donor, correlation=corfit$consensus)
contrast <- makeContrasts(E-NE, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)

tt <- topTable(fit2, number=Inf, sort.by="none")

out <- cbind(df[keep, c("Protein Group","Accession","Description")], tt)

out_dir <- "C:/Users/vishn/Developer/experiments/zaro/volk-ms-1/Volk et al MS Data Analysis for VRT/analysis_outputs"
write.csv(out, file.path(out_dir, "limma_4B_heavy_donorblock_LOD_log2.csv"), row.names=FALSE)

sig <- sum(out$adj.P.Val < 0.05)
up <- sum(out$adj.P.Val < 0.05 & out$logFC > 0)

cat("LOD:", LOD, "\n")
cat("tested:", nrow(out), "\n")
cat("FDR<0.05:", sig, "\n")
cat("Up in Eaters:", up, "\n")
