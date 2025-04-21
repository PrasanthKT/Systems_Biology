library(affy)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133a.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

#setwd("/Users/prasanthkumar/Desktop/Course_Work/Spring_2025/SB/Assignment-4/GSE11352_RAW")

raw_data <- ReadAffy()
norm_data <- rma(raw_data)
expr_matrix <- exprs(norm_data)

group <- factor(c(rep("Control", 6), rep("E2_12h", 6), rep("E2_24h", 6)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(E2_12h - Control, E2_24h - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


de_12h <- topTable(fit2, coef = "E2_12h - Control", number = Inf)
symbol_ids_12h <- mapIds(hgu133a.db, keys = rownames(de_12h),
                         column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
gene_list_12h <- de_12h$logFC
names(gene_list_12h) <- symbol_ids_12h
gene_list_12h <- sort(na.omit(gene_list_12h), decreasing = TRUE)

de_24h <- topTable(fit2, coef = "E2_24h - Control", number = Inf)
symbol_ids_24h <- mapIds(hgu133a.db, keys = rownames(de_24h),
                         column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
gene_list_24h <- de_24h$logFC
names(gene_list_24h) <- symbol_ids_24h
gene_list_24h <- sort(na.omit(gene_list_24h), decreasing = TRUE)


gmt_file <- "/Users/prasanthkumar/Downloads/Human_GO_AllPathways_noPFOCR_with_GO_iea_March_01_2025_symbol.gmt.txt"
gene_sets <- read.gmt(gmt_file)


gsea_12h <- GSEA(geneList = gene_list_12h,
                 TERM2GENE = gene_sets,
                 pvalueCutoff = 0.05,
                 verbose = TRUE)

gsea_24h <- GSEA(geneList = gene_list_24h,
                 TERM2GENE = gene_sets,
                 pvalueCutoff = 0.05,
                 verbose = TRUE)


gsea_12h@result$Description <- gsub("%.*", "", gsea_12h@result$Description)
dotplot(gsea_12h, showCategory = 20, title = "GSEA Dotplot – 12h Treated vs Control")

gsea_24h@result$Description <- gsub("%.*", "", gsea_24h@result$Description)
dotplot(gsea_24h, showCategory = 20, title = "GSEA Dotplot – 24h Treated vs Control")

ridgeplot(gsea_12h, showCategory = 20) +
  ggtitle("GSEA Ridgeplot – 12h Treated vs Control")

ridgeplot(gsea_24h, showCategory = 20) +
  ggtitle("GSEA Ridgeplot – 24h Treated vs Control")

gsea_12h <- pairwise_termsim(gsea_12h)
gsea_24h <- pairwise_termsim(gsea_24h)

emapplot(gsea_12h, showCategory = 25, title = "Enrichment Map – 12h Treated vs Control")
emapplot(gsea_24h, showCategory = 25, title = "Enrichment Map – 24h Treated vs Control")

write.table(data.frame(NAME=names(gene_list_12h), RANK=gene_list_12h),
            file="GSEA_12h.rnk", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(NAME=names(gene_list_24h), RANK=gene_list_24h),
            file="GSEA_24h.rnk", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(gsea_12h@result, file = "GSEA_12h_Results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_24h@result, file = "GSEA_24h_Results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gsea_12h_df <- gsea_12h@result[, c("ID", "Description", "NES", "p.adjust")]
gsea_24h_df <- gsea_24h@result[, c("ID", "Description", "NES", "p.adjust")]

colnames(gsea_12h_df)[3] <- "NES_12h"
colnames(gsea_24h_df)[3] <- "NES_24h"

merged_df <- merge(gsea_12h_df, gsea_24h_df, by = "Description", all = TRUE)

nes_matrix <- as.matrix(merged_df[, c("NES_12h", "NES_24h")])
rownames(nes_matrix) <- merged_df$Description

nes_matrix <- na.omit(nes_matrix)

common_terms <- intersect(gsea_12h@result$ID, gsea_24h@result$ID)

nes_12h <- gsea_12h@result %>%
  filter(ID %in% common_terms) %>%
  arrange(ID) %>%
  pull(NES)

nes_24h <- gsea_24h@result %>%
  filter(ID %in% common_terms) %>%
  arrange(ID) %>%
  pull(NES)

descriptions <- gsea_12h@result %>%
  filter(ID %in% common_terms) %>%
  arrange(ID) %>%
  pull(Description)

unique_desc <- make.unique(descriptions)

nes_matrix <- data.frame(`12h` = nes_12h, `24h` = nes_24h)
rownames(nes_matrix) <- unique_desc

top_nes_matrix <- nes_matrix[1:30, ]  


heatmap(as.matrix(top_nes_matrix),
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Top 30 NES Heatmap: 12h vs 24h",
        xlab = "Timepoint",
        margins = c(6, 15))  
