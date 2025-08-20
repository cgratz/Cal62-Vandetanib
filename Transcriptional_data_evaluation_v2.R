################################################################################
# Title: RNA-Seq Analysis Pipeline for Differential Expression and Biomarker Discovery
# Description: DESeq2, sPLS-DA, and Over-Representation Analysis (ORA) for publication.
# Author: Christian Grätz
# Version: 1.0
# Last updated: 2025-05-15
#
# This script performs:
# - Normalization of count data with DESeq2
# - Differential expression analysis
# - Visualization (volcano plots, PCA, heatmaps)
# - Multivariate modeling (sPLS-DA)
# - Functional enrichment (GO, KEGG, Reactome, WikiPathways, MSigDB)
################################################################################

# Core Analysis
library(DESeq2)
library(BiocParallel)

# Data Wrangling
library(dplyr)
library(tidyverse)

# Visualization
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggVennDiagram)

# Multivariate Modeling
library(mixOmics)

# Functional Enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(ReactomePA)
library(msigdbr)

# Plot Formatting
library(RColorBrewer)
library(grid)

# Use 4 CPU cores
bpparam <- SnowParam(workers = 4, type = "SOCK")

###############################################################################
## == Perform DESeq2 normalization ==
###############################################################################
# read gene counts and save gene name identifiers
counts <- read.table("salmon.merged.gene_counts.tsv", header = TRUE, row.names = 1, sep = "\t")
gene_names <- as.data.frame(rownames(counts))
gene_names[,2] <- as.data.frame(counts[,1])
colnames(gene_names) <- c("gene_id", "gene_name")
counts <- counts[,c(2:13)]

# create metadata table
Condition = c(rep("DMSO", 3), rep("EC10", 3), rep("EC50", 3), rep("EC80", 3))
coldata <- data.frame(
  row.names = colnames(counts),
  Condition = Condition
)

comparisonList <- as.data.frame(cbind(c("EC10", "EC50", "EC80", "EC50", "EC80", "EC80"),    # create table of comparisons to be made
                                      c("DMSO", "DMSO", "DMSO", "EC10", "EC10", "EC50")))
counts <- round(counts)   # round salmon values to integers
dds <- DESeqDataSetFromMatrix(countData = counts,      # create DESeqDataSet
                              colData = coldata,
                              design = ~ Condition)   
dds_backup <- dds

# run DESeq2 analysis on dds object
dds <- DESeq(dds)
plotDispEsts(dds)
DESeq_counts <- as.data.frame(DESeq2::counts(dds, normalized = T))
DESeq_counts$gene_id <- rownames(DESeq_counts)
DESeq_counts <- merge(DESeq_counts, gene_names, by = "gene_id")
write.csv(DESeq_counts, file = "DESeq2_cellular_counts_normalized.csv", row.names = F)    # save normalized DESeq2 counts

# extract results for all comparisons
fullResults <- as.data.frame(NULL)
volcanoData <- as.data.frame(NULL)
for (i in 1:nrow(comparisonList)) {
  res <- as.data.frame(
    results(dds, alpha = 0.05, contrast = c("Condition", 
                                            comparisonList[i,1],
                                            comparisonList[i,2]))
  )
  res$gene_id <- rownames(res)
  res2 <- merge(gene_names, res, by = "gene_id")    # add gene names
  rownames(res2) <- rownames(res)
  res2$comparison <- rep(paste(comparisonList[i,1], # add comparision name
                               "vs",
                               comparisonList[i,2],
                               sep = "_"), 
                         nrow(res2))
  volcanoData <- rbind(volcanoData, res2)           # save results for later volcano plot
  res2 <- res2[res2$padj < 0.05 & abs(res2$log2FoldChange) > 1 & !is.na(res2$padj) , ]  # remove failed padj calculations 
                                                                                        # & log2FC < 1 / > -1
  res2 <- res2[order(res2$padj),]                   # sort results by padj value
  fullResults <- rbind(fullResults, res2)           # save results in fullResults table
}
upReg <- fullResults[fullResults$log2FoldChange > 0,]       # extract and save only upregulated genes
downReg <- fullResults[fullResults$log2FoldChange < 0,]     # extract and save only downregulated genes

###############################################################################
## == Create Venn Diagram ==
###############################################################################

up_gene_lists <- list(                                                  #create list of upregulated genes
  EC10 = upReg[grep("EC10_vs_DMSO",upReg$comparison), "gene_name"],
  EC50 = upReg[grep("EC50_vs_DMSO",upReg$comparison), "gene_name"],
  EC80 = upReg[grep("EC80_vs_DMSO",upReg$comparison), "gene_name"]
)
down_gene_lists <- list(                                                # create list of downregulated genes
  EC10 = downReg[grep("EC10_vs_DMSO",downReg$comparison), "gene_name"],
  EC50 = downReg[grep("EC50_vs_DMSO",downReg$comparison), "gene_name"],
  EC80 = downReg[grep("EC80_vs_DMSO",downReg$comparison), "gene_name"]
)
upReg$dir <- "up"
downReg$dir <- "down"
up_gene_list <- list(
  EC10 = upReg[grep("EC10_vs_DMSO",upReg$comparison), colnames(upReg) %in% c("gene_name", "dir") ],
  EC50 = upReg[grep("EC50_vs_DMSO",upReg$comparison),  colnames(upReg) %in% c("gene_name", "dir")],
  EC80 = upReg[grep("EC80_vs_DMSO",upReg$comparison),  colnames(upReg) %in% c("gene_name", "dir")]
)
down_gene_list <- list(
  EC10 = downReg[grep("EC10_vs_DMSO",downReg$comparison), colnames(downReg) %in% c("gene_name", "dir") ],
  EC50 = downReg[grep("EC50_vs_DMSO",downReg$comparison),  colnames(downReg) %in% c("gene_name", "dir")],
  EC80 = downReg[grep("EC80_vs_DMSO",downReg$comparison),  colnames(downReg) %in% c("gene_name", "dir")]
)
reg_list <- list(
  EC10 = rbind(up_gene_list$EC10, down_gene_list$EC10),
  EC50 = rbind(up_gene_list$EC50, down_gene_list$EC50),
  EC80 = rbind(up_gene_list$EC80, down_gene_list$EC80)
)
reg_list$EC10$compare <- paste(reg_list$EC10$gene_name, reg_list$EC10$dir, sep = "_")
reg_list$EC50$compare <- paste(reg_list$EC50$gene_name, reg_list$EC50$dir, sep = "_")
reg_list$EC80$compare <- paste(reg_list$EC80$gene_name, reg_list$EC80$dir, sep = "_")
venn_list  <- list(
  EC10 = reg_list$EC10$compare,
  EC50 = reg_list$EC50$compare,
  EC80 = reg_list$EC80$compare
)

# Venn Diagram
ggVennDiagram(venn_list, label = "count", label_alpha = 0.5, set_size = 8, label_size = 10)+
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  ggtitle("Venn Diagram: Differentially expressed genes identified by RNA-Seq (Treatment vs DMSO)") +
  theme(plot.title = element_text(hjust = 0, vjust = 1, size = 11, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) 
ggsave("VennDiagram_RNA_cellular_color.png", width = 4000, height = 4000, dpi = 600, units = "px", bg = "white")

###############################################################################
## == Create Volcano Plots ==
###############################################################################
volcanoData$significant <- "not significant"
volcanoData$significant[volcanoData$padj < 0.05 & abs(volcanoData$log2FoldChange) >1] <- "significant"
data <- volcanoData[volcanoData$comparison == unique(volcanoData$comparison)[1], ]  #adjust for each comparison
data$gene_category <- ifelse(
  data$significant == "significant" & data$log2FoldChange > 1, "Upregulated",
  ifelse(data$significant == "significant" & data$log2FoldChange < -1, "Downregulated", "Not Significant")
)
data$padj[data$padj == 0] <- 1e-323
ggplot(data, 
       aes(x = log2FoldChange, 
           y = -log10(padj), 
           color = gene_category,
           shape = gene_category)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "firebrick3", "Downregulated" = "dodgerblue3")) +
  scale_shape_manual( values = c("Not Significant" = 16, "Upregulated" = 17, "Downregulated" = 15)) +
  theme_minimal() +
  labs(
    title = paste("Volcano Plot for", unique(volcanoData$comparison)[1]),   #adjust for each comparison
    x = "log2 Fold Change",
    y = "-log10(padj)"
  ) +
  theme (
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 32),
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18),  
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),  
    legend.text = element_text(size = 16)  
  ) +
  geom_vline(xintercept = 1, color = "grey", size = 1.5) +
  geom_vline(xintercept = -1, color = "grey", size = 1.5) +
  geom_hline(yintercept = -log10(0.05), color = "grey", size = 1.5) +
  geom_text_repel(
    aes(label = ifelse(significant == "significant", gene_name, NA)),
    max.overlaps = 30,   #adjust for each comparison
    box.padding = 0.2,
    point.padding = 0.2,
    size = 3,
    segment.size = 0.2,
    force = 2     #adjust for each comparison
  ) +
  ylim(0, 200)    #adjust for each comparison
ggsave("volcano_plot_cell_EC10_vs_DMSO_col.png", height = 5000, width = 5000, unit = "px", dpi = 600, bg = "white")

###############################################################################
## == Transform Data for subsequent visualizations ==
###############################################################################
vsd <- vst(dds, blind= F)

###############################################################################
## == Principal Component Analysis==
###############################################################################

# filter out 0 variance and transpose data
vsd_data <- assay(vsd)
PCA_data <- t(vsd_data)[ ,apply(t(vsd_data), 2, var, na.rm=TRUE) != 0]

# choose top 500 genes with highest variance
PCA_data_select <- PCA_data[, order(colVars(PCA_data),decreasing=TRUE)[1:500]]
PCA_data_select <- PCA_data_select[,!is.na(colSums(PCA_data_select))]

# create design to store the group information in
PCA_design <- colData(dds)$Condition
PCA_design <- factor(PCA_design, levels = c("DMSO", "EC10", "EC50", "EC80"))
pch_values <- c("DMSO" = 0, "EC10" = 1, "EC50" = 2, "EC80" = 5)
PCA_model <- pca(PCA_data_select, ncomp = 2)

png("PCA_cellular_col.png", width = 4000, height = 4000, res = 600)
plotIndiv(
  PCA_model, 
  group = PCA_design,
  comp = c(1,2), 
  legend = T, 
  ellipse = T,
  title ="PCA plot",
  legend.title = "Treatment",
  point.lwd = 1.5,
  pch = pch_values[PCA_design],
  cex = 4,
  size.title = 8,
  size.xlabel = 22,
  size.ylabel = 22,
  size.legend = 22,
  size.legend.title = 24,
  col.per.group = c("black", "goldenrod2", "darkorange3", "firebrick4")
)
dev.off()

#Heatmap
expr_matrix <- as.data.frame(assay(vsd))    #extract expression data from vsd object
top_genes <- head(rownames(fullResults[grep("EC80_vs_DMSO", fullResults$comparison),]), 50)   # filter for 500 top padj genes
top_genes <- intersect(intersect(reg_list[[1]]$gene_name, reg_list[[2]]$gene_name), reg_list[[3]]$gene_name)    # filter for genes that are differentially expressed in all doses
top_genes <- DESeq_counts$gene_id[make.unique(DESeq_counts$gene_name) %in% top_genes]
expr_matrix <- as.matrix(expr_matrix[rownames(expr_matrix) %in% top_genes,])

# put gene names as row names
gene_id <- as.data.frame(rownames(expr_matrix))
colnames(gene_id) <- "gene_id"
gene_id <- merge(gene_id, gene_names, by = "gene_id")
gene_id$gene_name <- make.unique(gene_id$gene_name)
rownames(expr_matrix) <- gene_id$gene_name

# create heatmap
sample_annotation <- as.data.frame(colData(dds)["Condition"])
colnames(sample_annotation) <- "Treatment"
sample_ha <- HeatmapAnnotation(
  Treatment = sample_annotation$Treatment,
  col = list(Treatment = c("DMSO" = "black", "EC10" = "goldenrod2", 
                           "EC50" = "darkorange3", "EC80" = "firebrick4")),
  annotation_name_gp = gpar(fontsize = 8)
)
YlGnBu_colors <- rev(brewer.pal(3, "YlGnBu"))
col_fun <- colorRamp2(
  c(min(expr_matrix), quantile(expr_matrix, prob = 0.25), median(expr_matrix), quantile(expr_matrix, prob = 0.75), max(expr_matrix)),
  c("azure1", "lightblue1", "lightgreen", "goldenrod3", "darkred")
)
ht <-  Heatmap(
  expr_matrix,
  name = "Expression",
  top_annotation = sample_ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 7),
  col = col_fun
)
png( "heatmap_DE_genes_color.png", height = 6000, width = 3000, res = 600)
draw(ht)
dev.off()

###############################################################################
## == sparse Partial Least Squares Regression ==
###############################################################################
expr_matrix <- assay(vsd)

# filter out rows with zero variance (all identical values)
expr_matrix <- expr_matrix[MatrixGenerics::rowVars(expr_matrix) > 0, ]
temp <- as.data.frame(rep("a", times = nrow(expr_matrix)))
colnames(temp) <- "gene_id"
temp$gene_id <- rownames(expr_matrix)
temp <- merge(temp, gene_names)
rownames(expr_matrix) <- make.unique(temp$gene_name)
rm(temp)
sample_group <- factor(rep(c("DMSO", "EC10", "EC50", "EC80"), each = 3))

#create sPLS-DA model for fine-tuning
splsda_model <- splsda(
  X = t(expr_matrix),
  Y = sample_group,
  ncomp =3,
  keepX = c(100, 100, 100)  #start with 100 factors for each dimension
)

# cross-validation of the model -> takes very long to compute, perform only once, in future load the saved results
# Use cross-validation to choose optimal ncomp (number of components)
cv_result <- perf(splsda_model, validation = "Mfold", folds = 3, progressBar = T, auc = TRUE, nrepeat = 100, 
                  BPPARAM = bpparam)
save(cv_result, file = "cv_result.RData")
# load("cv_result.RData")
plot(cv_result, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
cv_result$choice.ncomp
cv_result$auc

# ROC curve and AUC
auc_values <- auroc(splsda_model)
predicted_scores <- predict(splsda_model, t(expr_matrix), dist ="centroids.dist")$predict
normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
normalized_scores <- apply(predicted_scores, 2, normalize)
summary(predicted_scores)
list.keepX <- c(1:100)
tune.splsda_model <- tune.splsda(X = t(expr_matrix), Y = sample_group, ncomp = 3, 
                                 validation = 'Mfold', folds = 3,        #keep folds at 3, because only 3 samples per group
                                 progressBar = TRUE, dist = 'max.dist',  
                                 measure = "overall", test.keepX = list.keepX, nrepeat = 50, cpus = 12
)
save(tune.splsda_model, file = "tune_splsda_model.RData")
# load("tune_splsda_model.RData")
tune.splsda_model$choice.keepX
tune.splsda_model$choice.ncomp$ncomp
tune.splsda_model$features.selection.frequency

# final sPLS-DA model with fine-tuned parameters
splsda_model <- splsda(
  X = t(expr_matrix),
  Y = sample_group,
  ncomp =3,
  keepX = c(50, 40, 25)  #50, 40, 25 was identified as ideal parameters
)
selectVar(splsda_model , comp = 1)$name
perf_res <- perf(splsda_model , validation = "Mfold", folds = 3, nrepeat = 50, progressBar = TRUE)
print(perf_res$error.rate)

# save sPLS-DA plot as png
png("cellular_sPLS-DA_Components_1-2_col.png", width = 4000, height = 4000, res = 600)
plotIndiv(
  splsda_model, 
  group = sample_group,
  comp = c(1,2), 
  legend = T, 
  ellipse = T,
  title ="sPLS-DA cellular RNA Comp. 1 & 2 - only genes with differential expression in *both* EC10 and EC80 compared to DMSO",
  legend.title = "Treatment",
  point.lwd = 1.5,
  pch = c(0, 1, 2, 5),
  cex = 5,
  size.title = 6,
  size.xlabel = 22,
  size.ylabel = 22,
  size.legend = 22,
  size.legend.title = 24,
  col.per.group = c("black", "goldenrod2", "darkorange3", "firebrick4")
)
dev.off()
png("cellular_sPLS-DA_Components_2-3_col.png", width = 4000, height = 4000, res = 600)
plotIndiv(
  splsda_model, 
  group = sample_group,
  comp = c(2,3), 
  legend = T, 
  ellipse = T,
  title ="sPLS-DA cellular RNA Comp. 2 & 3",
  legend.title = "Treatment",
  point.lwd = 1.5,
  pch = c(0, 1, 2, 5),
  cex = 5,
  size.title = 10,
  size.xlabel = 22,
  size.ylabel = 22,
  size.legend = 22,
  size.legend.title = 24,
  col.per.group = c("black", "goldenrod2", "darkorange3", "firebrick4"))
dev.off()

# read out and export the transcripts that contribute to the components
loadings_matrix <- loadings(splsda_model)
important_features_comp1 <- as.data.frame(loadings_matrix$X[,1])
important_features_comp1$abs <- abs(important_features_comp1[,1])
important_features_comp1 <- head(important_features_comp1[order(important_features_comp1$abs, decreasing = T),],50)
important_features_comp1 <- important_features_comp1[important_features_comp1[,1] != 0,]
loadings_list_comp1 <- as.data.frame(important_features_comp1[,1])
row.names(loadings_list_comp1) <- row.names(important_features_comp1)
important_features_comp2 <- as.data.frame(loadings_matrix$X[,2])
important_features_comp2$abs <- abs(important_features_comp2[,1])
important_features_comp2 <- head(important_features_comp2[order(important_features_comp2$abs, decreasing = T),],50)
important_features_comp2 <- important_features_comp2[important_features_comp2[,1] != 0,]
loadings_list_comp2 <- as.data.frame(important_features_comp2[,1])
row.names(loadings_list_comp2) <- row.names(important_features_comp2)
important_features_comp3 <- as.data.frame(loadings_matrix$X[,3])
important_features_comp3$abs <- abs(important_features_comp3[,1])
important_features_comp3 <- head(important_features_comp3[order(important_features_comp3$abs, decreasing = T),],50)
important_features_comp3 <- important_features_comp3[important_features_comp3[,1] != 0,]
loadings_list_comp3 <- as.data.frame(important_features_comp3[,1])
row.names(loadings_list_comp3) <- row.names(important_features_comp3)

# put all genes in one table
loadings_list_full <- as.data.frame(rep(NA, nrow(loadings_list_comp1)))
loadings_list_full$factors_comp1 <- rownames(loadings_list_comp1)
loadings_list_full$contr_comp1 <- loadings_list_comp1[,1]
loadings_list_full$factors_comp2 <- c(rownames(loadings_list_comp2), 
                                      rep(NA,nrow(loadings_list_full)-nrow(loadings_list_comp2)))
loadings_list_full$contr_comp2 <- c(loadings_list_comp2[,1],
                                    rep(NA, nrow(loadings_list_full)-nrow(loadings_list_comp2)))
loadings_list_full$factors_comp3 <-  c(rownames(loadings_list_comp3), 
                                       rep(NA,nrow(loadings_list_full)-nrow(loadings_list_comp3)))
loadings_list_full$contr_comp3 <- c(loadings_list_comp3[,1],
                                    rep(NA, nrow(loadings_list_full)-nrow(loadings_list_comp3)))
loadings_list_full <- loadings_list_full[,-1]
sPLSDA_genes_cell <- unique(c(
  loadings_list_full$factors_comp1, 
  loadings_list_full$factors_comp2, 
  loadings_list_full$factors_comp3
))

# check which genes from sPLS-DA are differentially expressed
sPLSDA_DE_data <- fullResults[fullResults$gene_name %in% sPLSDA_genes_cell, ]
sPLSDA_DE <- unique(sPLSDA_DE_data$gene_name)
sPLSDA_DE_data$comp1 <- sPLSDA_DE_data$gene_name %in% rownames(loadings_list_comp1)
sPLSDA_DE_data$comp2 <- sPLSDA_DE_data$gene_name %in% rownames(loadings_list_comp2)
sPLSDA_DE_data$comp3 <- sPLSDA_DE_data$gene_name %in% rownames(loadings_list_comp3)
write.csv(sPLSDA_DE_data, file = "cell_sPLSDA_DE.csv", row.names = F)

# Heatmap with expression data from genes of sPLS-DA comp1, comp2 and comp3
heat_matrix <- as.data.frame(
  expr_matrix[rownames(expr_matrix) %in% 
                c(rownames(loadings_list_comp1), 
                  rownames(loadings_list_comp2),
                  rownames(loadings_list_comp3)),]
)
DEvsDMSO <- sPLSDA_DE_data[sPLSDA_DE_data$comparison %in% c("EC10_vs_DMSO","EC50_vs_DMSO","EC80_vs_DMSO"),]
heat_matrix <- heat_matrix[rownames(heat_matrix) %in% DEvsDMSO$gene_name, ]

# create heatmap
sample_annotation <- colData(dds)["Condition"]
colnames(sample_annotation) <- "Treatment"
sample_ha <- HeatmapAnnotation(
  Treatment = sample_annotation$Treatment,
  col = list(Treatment = c("DMSO" = "black", "EC10" = "royalblue", "EC50" = "darkcyan", "EC80" = "darkorange")),
  annotation_name_gp = gpar(fontsize = 8)
)

# Heatmap with clustering
heat_matrix <- as.matrix(heat_matrix)
col_fun <- colorRamp2(
  c(min(heat_matrix), quantile(heat_matrix, prob = 0.25), median(heat_matrix), quantile(heat_matrix, prob = 0.75), 
    max(heat_matrix)),
  c("azure1", "lightblue1", "lightgreen", "goldenrod3", "darkred")
)
ht <- Heatmap(
  heat_matrix,
  name = "Expression (vst)",
  top_annotation = sample_ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  col = col_fun
)
png("cellular_sPLS-DA_heatmap_comp1+2+3_expression_col2_poster.png", width = 4000, height = 6000, res = 600)
draw(
  ht,
  column_title =  "                                                            Variance-stabilized cellular expression values for the significant genes from sPLS-DA components 1 + 2",
  column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  padding = unit(c(2,2,2,2), "mm")
)
dev.off()

###############################################################################
## == Perform Over-Representation Analysis ==
###############################################################################
biomarkers  <- data.frame(c("NDRG1", "ISG20", "HDAC5", "GPRC5C", "RDH10", 
                            "GSN", "MTHFR", "WNT9A", "AHNAK2", "ITGB3", 
                            "BIRC5", "TPX2", "LAMC3", "TOP2A", "NCAPD2", 
                            "FOXM1", "MYBL2", "CDC25A", "THBD", "MVP", "KLF9"))

DD_DE_genes_ENSEMBL <- DESeq_counts$gene_id[DESeq_counts$gene_name %in% biomarkers[,1]]

# convert to ENTREZID
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ensembl_to_entrez <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),  
  filters = "ensembl_gene_id",                        
  values = DD_DE_genes_ENSEMBL,                         
  mart = ensembl
)
DD_DE_genes_ENTREZ <- unique(                                           # keep only unique ENTREZIDs
  ensembl_to_entrez$entrezgene_id[
    ensembl_to_entrez$ensembl_gene_id %in% DD_DE_genes_ENSEMBL
  ]
)
DD_DE_genes_ENTREZ <- DD_DE_genes_ENTREZ[!is.na(DD_DE_genes_ENTREZ)]    # remove NAs
background_ENSEMBL <- rownames(counts[rowSums(counts) >= 50, ])         # define Background as all genes that were detected with at least 50 counts total
ensembl_to_entrez <- getBM(                                             # convert to ENTREZIDs, keep only unique IDs and remove NAs
  attributes = c("ensembl_gene_id", "entrezgene_id"),  
  filters = "ensembl_gene_id",                        
  values = background_ENSEMBL,                         
  mart = ensembl
)
paste(
  round(sum(is.na(ensembl_to_entrez$entrezgene_id)) / nrow(ensembl_to_entrez) *100,2), 
  "% of gene ids could not be matched!", 
  sep = ""
)
background_ENTREZ <- unique(
  ensembl_to_entrez$entrezgene_id[
    ensembl_to_entrez$ensembl_gene_id %in% background_ENSEMBL
  ]
)
background_ENTREZ <- as.character(background_ENTREZ[!is.na(background_ENTREZ)])

# define genes and universe, so it can be switched quickly if needed
genes <- DD_DE_genes_ENTREZ
universe <- background_ENTREZ

# Read in pathways data (check for updates of the gmt file!)
wikipathways_gmt <- read.gmt("wikipathways-20250210-gmt-Homo_sapiens.gmt")  #gmt file from WikiPathways Website

#ORA for all databases 
ORA_KEGG <- enrichKEGG(    # KEGG ORA
  gene = genes,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)
ORA_GOBP <- enrichGO(    # Gene Ontolog BP ORA
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",               #biological process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)
ORA_GOMF <- enrichGO(    # Gene Ontology MF ORA
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",               #molecular function
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)
ORA_Reactome <- enrichPathway(    # Reactome ORA
  gene = genes,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)
ORA_MSigDB <- enricher(    # MSigDB ORA
  gene = genes,
  TERM2GENE = msigdbr(species = "Homo sapiens", category = "H")[, c("gs_name", "entrez_gene")],
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)
ORA_WP <- enricher(    # WikiPathways ORA
  gene = genes,
  TERM2GENE = wikipathways_gmt,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe = universe
)

# remove the gmt file description from the pathways results:
ORA_WP@result$ID <- gsub(".*(WP[0-9]+).*", "\\1", ORA_WP@result$ID)
ORA_WP@result$Description <- gsub("%WikiPathways_.*Homo sapiens", "", ORA_WP@result$Description)
rownames(ORA_WP@result) <- gsub("%WikiPathways_20250210%", "_", rownames(ORA_WP@result))
rownames(ORA_WP@result) <- gsub("%Homo sapiens", "", rownames(ORA_WP@result))

#convert ENTREZ IDs back to gene symbols for all datasets: -> function that can be applied to all results
convert_to_symbol <- function(entrez_list) {
  symbols <- unname(entrez_map[unlist(strsplit(entrez_list, "/"))])
  symbols <- symbols[!is.na(symbols)]
  paste(symbols, collapse ="/")
}

# create mapping tables for all results
all_entrez_ids <- unique(unlist(c(
  strsplit(ORA_WP@result$geneID, "/"),
  strsplit(ORA_GOBP@result$geneID, "/"),
  strsplit(ORA_GOMF@result$geneID, "/"),
  strsplit(ORA_KEGG@result$geneID, "/"),
  strsplit(ORA_MSigDB@result$geneID, "/"),
  strsplit(ORA_Reactome@result$geneID, "/")
)))
entrez_to_symbol <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = all_entrez_ids,
  mart = ensembl
) %>%      
  filter(hgnc_symbol != "") #filter out empty symbols

# check for duplicates
if (any(duplicated(entrez_to_symbol$entrezgene_id))) {
  print("WARNING: Duplicate mappings for one or more Entrez IDs!")
} else {
  print("All Entrez IDs mapped uniquely")
}

# create mapping vector
entrez_map <- setNames(entrez_to_symbol$hgnc_symbol, entrez_to_symbol$entrezgene_id)

# apply the mapping function to all results:
ORA_WP@result$geneSymbol <- sapply(ORA_WP@result$geneID, convert_to_symbol)
ORA_KEGG@result$geneSymbol <- sapply(ORA_KEGG@result$geneID, convert_to_symbol)
ORA_GOMF@result$geneSymbol <- sapply(ORA_GOMF@result$geneID, convert_to_symbol)
ORA_GOBP@result$geneSymbol <- sapply(ORA_GOBP@result$geneID, convert_to_symbol)
ORA_MSigDB@result$geneSymbol <- sapply(ORA_MSigDB@result$geneID, convert_to_symbol)
ORA_Reactome@result$geneSymbol <- sapply(ORA_Reactome@result$geneID, convert_to_symbol)

# save results
save(ORA_WP, file = "ORA_WP.RData")
save(ORA_KEGG, file = "ORA_KEGG.RData")
save(ORA_GOBP, file = "ORA_GOBP.RData")
save(ORA_GOMF, file = "ORA_GOMF.RData")
save(ORA_Reactome, file = "ORA_Reactome.RData")
save(ORA_MSigDB, file = "ORA_MSigDB.RData")
# load("ORA_WP.RData")
# load("ORA_KEGG.RData")
# load("ORA_GOBP.RData")
# load("ORA_GOMF.RData")
# load("ORA_Reactome.RData")
# load("ORA_MSigDB.RData")
# combine all results into one table
ORA_results <- rbind(ORA_WP@result, ORA_KEGG@result[,-(1:2)], ORA_GOBP@result, 
                     ORA_GOMF@result, ORA_Reactome@result, ORA_MSigDB@result)
ORA_results$category <- c(rep(NA, nrow(ORA_WP@result)), ORA_KEGG@result$category, 
                          rep(NA, (nrow(ORA_results) - nrow(ORA_WP@result) - nrow(ORA_KEGG@result))))
ORA_results$subcategory <- c(rep(NA, nrow(ORA_WP@result)), ORA_KEGG@result$subcategory, 
                             rep(NA, (nrow(ORA_results) - nrow(ORA_WP@result) - nrow(ORA_KEGG@result))))
ORA_results$DB <- c(rep("WikiPathways", nrow(ORA_WP@result)),
                    rep("KEGG", nrow(ORA_KEGG@result)),
                    rep("Gene Ontology: Biological Process", nrow(ORA_GOBP@result)),
                    rep("Gene Ontology: Molecular Function", nrow(ORA_GOMF@result)),
                    rep("Reactome", nrow(ORA_Reactome@result)),
                    rep("MSigDB", nrow(ORA_MSigDB@result))
)
write.csv(ORA_results, file = "ORA_results_all_21_biomarkers.csv", row.names = F)    

# filter for padj < 0.1
ORA_results_sig <- ORA_results[ORA_results$p.adjust<0.1,]

# plot as barchart; insert break after 50 characters
ORA_results_sig <- ORA_results[ORA_results$p.adjust<0.1,]
ORA_results_sig$Description <- gsub("\\s{2,}", " ", ORA_results_sig$Description)
ORA_results_sig$Description <- gsub("(.{1,80})(\\s)", "\\1\n", ORA_results_sig$Description, perl = TRUE)
ORA_results_sig$Description <- str_wrap(ORA_results_sig$Description, width = 50)
ggplot(ORA_results_sig, aes(x = reorder(Description, FoldEnrichment), y = Count, fill = DB)) +
  geom_bar(stat = "identity", show.legend = T, color = "black") +
  geom_text(aes(label = geneSymbol), hjust = -0.1, size = 2.2) +  # Add total count labels next to bars
  coord_flip(clip = "off") +  # Flip the axes for better readability
  labs(
    x = "",
    y = "Number of contributing transcripts",
    fill = "Database",
    title = "ORA results for 21 biomarker candidates with padj < 0.1, ordered by Fold Enrichment"
  ) +
  scale_fill_manual(values = c(
    
    "Gene Ontology: Biological Process" = "ivory1",  
    "Gene Ontology: Molecular Function" = "darkseagreen2", 
    "MSigDB" = "steelblue1",
    "Reactome" = "darkorange3",
    "WikiPathways" = "firebrick4"
  )) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7)) +
  theme(plot.title = element_text(hjust = 0.8)) +
  theme(legend.position = c(1.12, 0.8)) + 
  theme(plot.margin = margin(t = 10, r = 120, b = 10, l = 0, unit = "pt"))
ggsave("ORA_barchart_col.png", height = 6000, width = 5000, dpi = 600, units = "px", bg = "white")

###############################################################################
## == Partial Least Squares Regression for the 8 biomarker genes Cq values==
###############################################################################
deltaCq <- read.table("biomarker_deltaCq.txt", head = T, sep = "\t", row.names = 1)   #read in table containing all deltaCq values
sample_group <- factor(rep(c("DMSO", "EC10", "EC50", "EC80"), each = 3))
plsda_model <- plsda(deltaCq, sample_group, ncomp = 2)
png("PLS-DA_RT-qPCR_biomarkers_col.png", width = 4000, height = 4000, res = 600)
plotIndiv(
  plsda_model, 
  group = sample_group,
  comp = c(1,2), 
  legend = T, 
  ellipse = T,
  title ="PLS-DA RT-qPCR results - deltaCq values",
  legend.title = "Treatment",
  cex = 4,
  point.lwd = 1.5,
  pch = c(0,1,2,5),
  size.title = 14,
  size.xlabel = 22,
  size.ylabel = 22,
  size.legend = 22,
  size.legend.title = 24,
  col.per.group = c("black", "goldenrod2", "darkorange3", "firebrick4")
)
dev.off()

#ROC and AUC analysis
auc_values <- auroc(plsda_model)
# Extract the ROC data from auc_values$graph.Comp2$data
roc_data <- auc_values$graph.Comp2$data
roc_data <- roc_data %>%
  mutate(Outcome = str_replace_all(Outcome, 
                                   "(\\d+\\.\\d+)", 
                                   function(x) round(as.numeric(x), 2) %>% as.character()))
# Plot the ROC curve using ggplot2
#plot for publication - works in greyscale
ggplot(roc_data, aes(x = Specificity, y = Sensitivity, color = Outcome)) +
  geom_line(linewidth = 2) +  # Adjust line thickness with linewidth
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 2) +  # Add reference line
  labs(title = "ROC Curve for PLS Regression Model",
       x = "100 - Specificity % (False Positive Rate)",
       y = "Sensitivity % (True Positive Rate)") +
  scale_color_manual(values = c("black", "goldenrod2", "darkorange3", "firebrick4")) +
  theme_minimal() +
  theme(legend.title = element_text(size = 24), 
        legend.text = element_text(size = 22),
        axis.title.x = element_text(siz = 22),
        axis.title.y = element_text(size = 22),
        plot.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 20)
  )
ggsave("PLS-DA_RT-qPCR_ROC_AUC_col.png", width = 6000, height = 4000, units = "px", dpi = 600, bg= "white")