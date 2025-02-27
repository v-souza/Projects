###########################################################
###########################################################
###########################################################
########            0 hour              ###################
########                                ###################
###########################################################
###########################################################
###########################################################


library(tibble)
library(ggrepel)

eDat <- read.table("Count.txt", header=TRUE, sep="\t", row.names = 1)
gDat <- read.table("Group.txt", header=TRUE, sep="\t")


## Counts
## filter only samples CTL and Exe0h
suppressMessages({library(dplyr)})

# Select only 'CTL' ou 'Exe0h'
selected_columns <- eDat[, grepl("^CTL|^Exe0h", names(eDat))]

# new table
new_table1 <- cbind(rownames = rownames(eDat), selected_columns)

## remove the first colllumn
new_table1$rownames <- NULL

## Treat
new_group1 <- gDat[1:10, ]

## load EdgeR ##
library(edgeR)

# ATENTION!!
new_group1$treat

# Control treatment needs to appear as the first level for the Treat Vs. Control comparison
new_group1$treat <-factor(new_group1$treat, levels=c("CTL", "Exe0h"))

# Input data
y_adj <- DGEList(counts=new_table1, samples=new_group1, group=new_group1$treat)
design_adj <- model.matrix(~ treat - 1, data = new_group1)

# Filter lowly expressed genes
keep_adj <- rowSums(cpm(y_adj)>1) >= 2
sum(keep_adj)
y_adj <- y_adj[keep_adj, , keep.lib.sizes=FALSE]

## Explanation: ">= 2" is a condition to maintain the gene, which means that the gene needs to have at least two replicas with CPM greater than 1 to be maintained.

## Explanation: after filtering, 16131 genes remained that have sufficient expression in at least two replicates with CPM greater than 1. Genes that do not meet these criteria are removed from the dataset.


# Normalization for RNA composition
y_adj <- calcNormFactors(y_adj)
y_adj$samples
norm.expr_adj <- y_adj$samples

# norm counts
logcpm_adj <- cpm(y_adj, log=TRUE)


##########################################
##      Data analyses with edgeR        ##
##########################################


# Estimate dispersion with Cox-Reid profile-adjusted likelihood (CR) method
y_adj <- estimateDisp(y_adj, design_adj)
y_adj$common.dispersion


## Explanation: A low dispersion value indicates that the expression levels of genes are relatively consistent across replicates, suggesting that there is little variability beyond what is expected due to random sampling. Conversely, a high dispersion value indicates greater variability in gene expression levels among replicates.


## plot
plotBCV(y_adj)

# Fit Generalized Linear Model (GLM) 
fit <- glmFit(y_adj, design_adj)
colnames(fit)

# Differential gene expression analysis
lrt <- glmLRT(fit, contrast=c(-1,1))

## Coefficients used in the contrast for the likelihood ratio test
topTags(lrt)

## Explanation: It indicates that the comparison is between the coefficient for treatCTL (the reference group) and the coefficient for treatExe0h (the group of interest).

## Number of genes identified as differentially expressed (Disregarding parameters)
is.de <- decideTestsDGE(lrt)
summary(decideTestsDGE(lrt))

## plot MD
plotMD(lrt, status=is.de)

# Number of up and down-regulated. FDR < 5% and |log2 fold change| ≥ 1
summary(de <- decideTestsDGE(lrt, lfc=0))

## Explanation: The absolute value of 1 on this scale indicates that there was a doubling (2-fold) change in gene expression between the two conditions.

# MA-plot
detags <- rownames(y_adj)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-0, 0), col="blue")

# To export the results
result <- as.data.frame(topTags(lrt, n=1000000))
result <- add_column(result, symbol= rownames(result), .before = "logFC")

## Check all the genes Downregulated and Upregulated

# Remove unnecessary columns
result$logCPM <- NULL
result$LR <- NULL

# Remove rows with FDR < 0.05 and logFC >= 0
result$expression <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 0,
                            ifelse(result$logFC >= 0, 'Up', 'Down'),
                            'Stable')

# Select only genes that are upregulated or downregulated
result_final <- subset(result, expression %in% c("Up", "Down"))

# Remove row names
rownames(result_final) <- NULL


##########################################
############      PLOTS        ###########
##########################################

## volcano plot

# Extract top tags from lrt
result_edgeR <- as.data.frame(topTags(lrt, n=1000000, adjust.method = "BH"))

# Extract log-transformed counts per million
count_CPM <- as.data.frame(logcpm_adj)

## Convert row names into first column
result_edgeR <- add_column(result_edgeR, gene = rownames(result_edgeR), .before = "logFC")

count_CPM <- add_column(count_CPM, gene = rownames(count_CPM), .before = "CTL_1")

## Order by gene
result_edgeR <- result_edgeR[order(result_edgeR$gene),]
count_CPM <- count_CPM[order(count_CPM$gene),]

## Combine result_edgeR columns with count_CPM
count_CPM <- add_column(count_CPM, logFC = result_edgeR$logFC, .before = "CTL_1")
count_CPM <- add_column(count_CPM, FDR = result_edgeR$FDR, .before = "CTL_1")
count_CPM <- add_column(count_CPM, gene2 = result_edgeR$gene, .before = "CTL_1")

# Check if gene and gene2 are the same and in the same order
if(all(count_CPM$gene == count_CPM$gene2)) {
  print("The values in 'gene' and 'gene2' are the same.")
} else {
  print("The values in 'gene' and 'gene2' are different.")
  different_indices <- which(count_CPM$gene != count_CPM$gene2)
  print(paste("The value in 'gene' is different at indices:", different_indices))
}

## remove gene2
count_CPM$gene2 <- NULL

## Insert -log10 FDR
count_CPM <- add_column(count_CPM, logFDR = -log10(result_edgeR$FDR), .before = "CTL_1")

## Insert -log10 p-value
count_CPM <- add_column(count_CPM, logpvalue = -log10(result_edgeR$PValue), .before = "CTL_1")


## PLOT

# Assign count_CPM to expressao_genica
expressao_genica <- count_CPM

# Determine expression based on FDR and logFC thresholds
expressao_genica$expression <- ifelse(expressao_genica$FDR < 0.05 & abs(expressao_genica$logFC) >= 0,
                                      ifelse(expressao_genica$logFC >= 0, 'Up', 'Down'),
                                      'Stable')

# Count the number of genes categorized as Up and Down
sum(expressao_genica$expression == "Up")
sum(expressao_genica$expression == "Down")

# Create volcano plot

# Filter the data to include only "Down" and "Up" genes
expressao_genica_filtered <- expressao_genica[expressao_genica$expression %in% c("Down", "Up"), ]

# Order the filtered dataframe by logFC
expressao_genica_filtered <- expressao_genica_filtered[order(expressao_genica_filtered$logFC, decreasing = TRUE), ]

# Select the top 5 upregulated genes and the top 5 downregulated genes
top_genes_up <- head(subset(expressao_genica_filtered, expression == "Up"), 5)
top_genes_down <- head(subset(expressao_genica_filtered, expression == "Down"), 5)

# Combine top_genes_up and top_genes_down into a single dataframe
top_genes <- rbind(top_genes_up, top_genes_down)

## see
print(top_genes)

# volcano plot
plot1 <- ggplot(data = expressao_genica,
                aes(x = logFC,
                    y = logFDR,
                    colour = expression)) +
  geom_point(size = 3) +
  geom_label_repel(data = top_genes, aes(label = gene), size = 5, force = 5, 
                   fill = "white", color = "black", box.padding = unit(0.5, "lines")) +  # Labeling top 10 genes with a rectangle
  scale_color_manual(values = c("red", "grey", "#009933"),
                     name = "",
                     breaks = c("Up", "Stable", "Down")) +
  geom_vline(xintercept = c(-0, 0), lty = 4, col = "black", lwd = 0.8) +
  labs(x = "log2FC",
       y = "-log10 (adj.P.Val)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  ggtitle("CTL vs 0 hour")

# To display the volcano plot
plot1

## save as image
tiff("volcano_CTL_vs_0h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
plot1 
dev.off()


## heatmap

## Only genes Down and Upregulated
head(expressao_genica_filtered)

#remove the fisrt five cols
expressao_genica_filtered <- expressao_genica_filtered[, -(1:5)]

## remove a col name expression
expressao_genica_filtered$expression <- NULL

# Convert dataframe to matrix
samp.with.rownames <- as.matrix(expressao_genica_filtered)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2

library(RColorBrewer)
cols <- brewer.pal(9, "RdYlBu")
pal <- colorRampPalette(cols)

column_annotation <- as.data.frame((colnames(samp.with.rownames)))
rownames(column_annotation) <- column_annotation$`(colnames(samp.with.rownames))`
column_annotation$Group = c(rep("Control", 5), rep("0 hour", 5))
column_annotation$`(colnames(samp.with.rownames))` <- NULL

# Convertendo a coluna "Group" em um fator com os níveis corretos
column_annotation$Group <- factor(column_annotation$Group, levels = c("Control", "0 hour"))

# Verificando os níveis da coluna "Group" após a conversão
levels(column_annotation$Group)


cc <- pheatmap(samp.with.rownames, 
               annotation_col = column_annotation, 
               scale = "row",
               fontsize=12,
               show_rownames = FALSE, 
               annotation_names_row = FALSE,
               annotation_names_col = FALSE, 
               annotation_colors = list(Group = c(Control = "#74ADD1", `0 hour` = "#A6D96A")),
               show_colnames = FALSE, 
               col = rev(pal(50)))


## print
cc

## save as image
tiff("Heatmap_CTL_vs_0h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
cc 
dev.off()


## pathway

library("clusterProfiler")
library("org.Mm.eg.db")
library("pathview")
library("ggplot2")
library("cowplot")


### Control vs 0 hours ###

## import DEGs
DEGs <- read.delim("~/VanessaPD/differential_expression/DEGs_CTL_vs_0h_NOFC.txt")

## select only genes with FDR <0.05 
significant_genes <- DEGs[DEGs$FDR < 0.05, ]

# For KEGG pathway enrichment using the gseKEGG() function, we need to convert id types. We can use the bitr function for this (included in clusterProfiler). It is normal for this call to produce some messages / warnings.
# In the bitr function, the param fromType should be the same as keyType from the gseGO function above (the annotation source). This param is used again in the next two steps: creating dedup_ids and df2.
# toType: in the bitr function has to be one of the available options from keyTypes(org.Dm.eg.db) and must map to one of ‘kegg’, ‘ncbi-geneid’, ‘ncib-proteinid’ or ‘uniprot’ because gseKEGG() only accepts one of these 4 options as it’s keytype parameter.

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids = bitr(significant_genes$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(ids) = c("symbol", "EntrezID")
df2 = merge(significant_genes, ids, by = "symbol")

# we want the log2 fold change 
original_gene_list = df2$logFC

# name the vector
names(original_gene_list) <- df2$EntrezID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

## KEGG

kegg_ctl_vs_0h <- enrichKEGG(
  names(gene_list),
  organism = "mmu",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

# table
t_kegg_ctl_vs_0h <- as.data.frame(kegg_ctl_vs_0h)

## convert gene ID to Symbol
edox <- setReadable(kegg_ctl_vs_0h, 'org.Mm.eg.db', 'ENTREZID')

# table
t_kegg_ctl_vs_0h <- as.data.frame(edox)

## remove (mus musculos name)
edox@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", edox@result$Description, fixed = T)
edox@result$Description <- gsub(pattern = " , secretion and action", replacement = "", edox@result$Description, fixed = T)
edox@result$Description


## table
t_kegg_ctl_vs_0h <- as.data.frame(edox)


## plot
NCG_barplot <- barplot(edox, showCategory = c("Osteoclast differentiation",
                                              "Parathyroid hormone synthesis, secretion and action",
                                              "IL-17 signaling pathway",
                                              "Circadian rhythm",
                                              "Adipocytokine signaling pathway",
                                              "FoxO signaling pathway",
                                              "TNF signaling pathway",
                                              "MAPK signaling pathway",
                                              "Hippo signaling pathway",
                                              "Efferocytosis",
                                              "Apoptosis",
                                              "Signaling pathways regulating pluripotency of stem cells")) +
  ggtitle("KEGG: Kyoto Encyclopedia of Genes and Genomes (0 hour)") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

# print
NCG_barplot

## save as image
tiff("barplor_CTL_vs_0h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
NCG_barplot 
dev.off()


## CNETPLOT

# Category: Environmental Information Processing
cnetplot1_0h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("FoxO signaling pathway",
                                          "TNF signaling pathway",
                                          "MAPK signaling pathway",
                                          "Hippo signaling pathway")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (0 hour)\nCategory: Environmental Information Processing",
    subtitle = NULL
  )

# print
cnetplot1_0h

# Category: Cellular Processes
cnetplot2_0h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("Efferocytosis",
                                          "Apoptosis")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (0 hour) \nCategory: Cellular Processes",
    subtitle = NULL
  )

# pint
cnetplot2_0h


# Category: Organismal Systems
cnetplot3_0h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("Circadian rhythm",
                                          "Osteoclast differentiation",
                                          "Parathyroid hormone synthesis",
                                          "IL-17 signaling pathway",
                                          "Adipocytokine signaling pathway")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (0 hour) \nCategory: Organismal Systems",
    subtitle = NULL
  )

# print
cnetplot3_0h


# All categories

cnetplot4_0h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("Circadian rhythm",
                                          "Osteoclast differentiation",
                                          "Parathyroid hormone synthesis",
                                          "IL-17 signaling pathway",
                                          "Adipocytokine signaling pathway",
                                          "FoxO signaling pathway",
                                          "TNF signaling pathway",
                                          "MAPK signaling pathway",
                                          "Hippo signaling pathway",
                                          "Efferocytosis",
                                          "Apoptosis")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (0 hour)",
    subtitle = NULL
  )

# Print
cnetplot4_0h


## save as image
tiff("cnetplot_CTL_vs_0h_all_categories_April_19_2024.tiff", units="in", width=15, height=15, res=600)
cnetplot4_0h 
dev.off()


###########################################################
###########################################################
###########################################################
########            2 hour              ###################
########                                ###################
###########################################################
###########################################################
###########################################################


## Import dataset
eDat <- read.table("Count.txt", header=TRUE, sep="\t", row.names = 1)
gDat <- read.table("Group.txt", header=TRUE, sep="\t")


## Counts
## filter only samples CTL and Exe2h

# Select only 'CTL' ou 'Exe2h'
selected_columns <- eDat[, grepl("^CTL|^Exe2h", names(eDat))]

# new table
new_table1 <- cbind(rownames = rownames(eDat), selected_columns)

## remove the first colllumn
new_table1$rownames <- NULL

## Treat
new_group1 <- subset(gDat, treat == "CTL" | treat == "Exe2h")

# ATENTION!!
new_group1$treat

# Control treatment needs to appear as the first level for the Treat Vs. Control comparison
new_group1$treat <-factor(new_group1$treat, levels=c("CTL", "Exe2h"))

# Input data
y_adj <- DGEList(counts=new_table1, samples=new_group1, group=new_group1$treat)
design_adj <- model.matrix(~ treat - 1, data = new_group1)

# Filter lowly expressed genes
keep_adj <- rowSums(cpm(y_adj)>1) >= 2
sum(keep_adj)
y_adj <- y_adj[keep_adj, , keep.lib.sizes=FALSE]


## Explanation: ">= 2" is a condition to maintain the gene, which means that the gene needs to have at least two replicas with CPM greater than 1 to be maintained.

## Explanation: after filtering, 16010 genes remained that have sufficient expression in at least two replicates with CPM greater than 1. Genes that do not meet these criteria are removed from the dataset.

# Normalization for RNA composition
y_adj <- calcNormFactors(y_adj)
y_adj$samples
norm.expr_adj <- y_adj$samples

# norm counts
logcpm_adj <- cpm(y_adj, log=TRUE)

##########################################
##      Data analyses with edgeR        ##
##########################################


# Estimate dispersion with Cox-Reid profile-adjusted likelihood (CR) method
y_adj <- estimateDisp(y_adj, design_adj)
y_adj$common.dispersion


## Explanation: A low dispersion value indicates that the expression levels of genes are relatively consistent across replicates, suggesting that there is little variability beyond what is expected due to random sampling. Conversely, a high dispersion value indicates greater variability in gene expression levels among replicates.

## plot
plotBCV(y_adj)

# Fit Generalized Linear Model (GLM) 
fit <- glmFit(y_adj, design_adj)
colnames(fit)

# Differential gene expression analysis
lrt <- glmLRT(fit, contrast=c(-1,1))

## Coefficients used in the contrast for the likelihood ratio test
topTags(lrt)

## Explanation: It indicates that the comparison is between the coefficient for treatCTL (the reference group) and the coefficient for treatExe2h (the group of interest).

## Number of genes identified as differentially expressed (Disregarding parameters)
is.de <- decideTestsDGE(lrt)
summary(decideTestsDGE(lrt))

## plot MD
plotMD(lrt, status=is.de)

# Number of up and down-regulated. FDR < 5% and |log2 fold change| ≥ 0
summary(de <- decideTestsDGE(lrt, lfc=0))

# MA-plot
detags <- rownames(y_adj)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

# To export the results
result <- as.data.frame(topTags(lrt, n=1000000))
result <- add_column(result, symbol= rownames(result), .before = "logFC")

## Check all the genes Downregulated and Upregulated

# Remove unnecessary columns
result$logCPM <- NULL
result$LR <- NULL

# Remove rows with FDR < 0.05 and logFC >= 0
result$expression <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 0,
                            ifelse(result$logFC >= 0, 'Up', 'Down'),
                            'Stable')

# Select only genes that are upregulated or downregulated
result_final <- subset(result, expression %in% c("Up", "Down"))

# Remove row names
rownames(result_final) <- NULL

##########################################
############      PLOTS        ###########
##########################################

## volcano plot

# Extract top tags from lrt
result_edgeR <- as.data.frame(topTags(lrt, n=1000000, adjust.method = "BH"))

# Extract log-transformed counts per million
count_CPM <- as.data.frame(logcpm_adj)

## Convert row names into first column
result_edgeR <- add_column(result_edgeR, gene = rownames(result_edgeR), .before = "logFC")

count_CPM <- add_column(count_CPM, gene = rownames(count_CPM), .before = "CTL_1")

## Order by gene
result_edgeR <- result_edgeR[order(result_edgeR$gene),]
count_CPM <- count_CPM[order(count_CPM$gene),]

## Combine result_edgeR columns with count_CPM
count_CPM <- add_column(count_CPM, logFC = result_edgeR$logFC, .before = "CTL_1")
count_CPM <- add_column(count_CPM, FDR = result_edgeR$FDR, .before = "CTL_1")
count_CPM <- add_column(count_CPM, gene2 = result_edgeR$gene, .before = "CTL_1")

# Check if gene and gene2 are the same and in the same order
if(all(count_CPM$gene == count_CPM$gene2)) {
  print("The values in 'gene' and 'gene2' are the same.")
} else {
  print("The values in 'gene' and 'gene2' are different.")
  different_indices <- which(count_CPM$gene != count_CPM$gene2)
  print(paste("The value in 'gene' is different at indices:", different_indices))
}

## remove gene2
count_CPM$gene2 <- NULL

## Insert -log10 FDR
count_CPM <- add_column(count_CPM, logFDR = -log10(result_edgeR$FDR), .before = "CTL_1")

## Insert -log10 p-value
count_CPM <- add_column(count_CPM, logpvalue = -log10(result_edgeR$PValue), .before = "CTL_1")


## PLOT

# Assign count_CPM to expressao_genica
expressao_genica <- count_CPM

# Determine expression based on FDR and logFC thresholds
expressao_genica$expression <- ifelse(expressao_genica$FDR < 0.05 & abs(expressao_genica$logFC) >= 0,
                                      ifelse(expressao_genica$logFC >= 0, 'Up', 'Down'),
                                      'Stable')

# Count the number of genes categorized as Up and Down
sum(expressao_genica$expression == "Up")
sum(expressao_genica$expression == "Down")

# Create volcano plot

# Filter the data to include only "Down" and "Up" genes
expressao_genica_filtered <- expressao_genica[expressao_genica$expression %in% c("Down", "Up"), ]

# Order the filtered dataframe by logFC
expressao_genica_filtered <- expressao_genica_filtered[order(expressao_genica_filtered$logFC, decreasing = TRUE), ]

# Select the top 5 upregulated genes and the top 5 downregulated genes
top_genes_up <- head(subset(expressao_genica_filtered, expression == "Up"), 5)
top_genes_down <- head(subset(expressao_genica_filtered, expression == "Down"), 5)

# Combine top_genes_up and top_genes_down into a single dataframe
top_genes <- rbind(top_genes_up, top_genes_down)

## see
print(top_genes)

# volcano plot
library(ggrepel)

plot2 <- ggplot(data = expressao_genica,
                aes(x = logFC,
                    y = logFDR,
                    colour = expression)) +
  geom_point(size = 3) +
  geom_label_repel(data = top_genes, aes(label = gene), size = 5, force = 5, 
                   fill = "white", color = "black", box.padding = unit(0.5, "lines")) +  # Labeling top 10 genes with a rectangle
  scale_color_manual(values = c("red", "grey", "#009933"),
                     name = "",
                     breaks = c("Up", "Stable", "Down")) +
  geom_vline(xintercept = c(-0, 0), lty = 4, col = "black", lwd = 0.8) +
  labs(x = "log2FC",
       y = "-log10 (adj.P.Val)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  ggtitle("CTL vs 2 hours")


# To display the volcano plot
dev.off()
plot(rnorm(50), rnorm(50))

plot2

## save as image
tiff("volcano_CTL_vs_2h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
plot2 
dev.off()


## heatmap

## Only genes Down and Upregulated
head(expressao_genica_filtered)

#remove the fisrt five cols
expressao_genica_filtered <- expressao_genica_filtered[, -(1:5)]

## remove a col name expression
expressao_genica_filtered$expression <- NULL

# Convert dataframe to matrix
samp.with.rownames <- as.matrix(expressao_genica_filtered)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2

library(RColorBrewer)
cols <- brewer.pal(9, "RdYlBu")
pal <- colorRampPalette(cols)

column_annotation <- as.data.frame((colnames(samp.with.rownames)))
rownames(column_annotation) <- column_annotation$`(colnames(samp.with.rownames))`
column_annotation$Group = c(rep("Control", 5), rep("2 hours", 5))
column_annotation$`(colnames(samp.with.rownames))` <- NULL

# Convertendo a coluna "Group" em um fator com os níveis corretos
column_annotation$Group <- factor(column_annotation$Group, levels = c("Control", "2 hours"))

# Verificando os níveis da coluna "Group" após a conversão
levels(column_annotation$Group)


cc <- pheatmap(samp.with.rownames, 
               annotation_col = column_annotation, 
               scale = "row",
               fontsize=12,
               show_rownames = FALSE, 
               annotation_names_row = FALSE,
               annotation_names_col = FALSE,
               annotation_colors = list(Group = c(Control = "#74ADD1", `2 hours` = "purple")),
               show_colnames = FALSE, 
               col = rev(pal(50)))


## Print
cc

## save as image
tiff("Heatmap_CTL_vs_2h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
cc 
dev.off()


## pathway

library("clusterProfiler")
library("org.Mm.eg.db")
library("pathview")
library("ggplot2")
library("cowplot")


### Control vs 2 hours ###

## import DEGs
DEGs <- read.delim("~/VanessaPD/differential_expression/DEGs_CTL_vs_2h_NOFC.txt")

## select only genes with FDR <0.05 
significant_genes <- DEGs[DEGs$FDR < 0.05, ]

# For KEGG pathway enrichment using the gseKEGG() function, we need to convert id types. We can use the bitr function for this (included in clusterProfiler). It is normal for this call to produce some messages / warnings.
# In the bitr function, the param fromType should be the same as keyType from the gseGO function above (the annotation source). This param is used again in the next two steps: creating dedup_ids and df2.
# toType: in the bitr function has to be one of the available options from keyTypes(org.Dm.eg.db) and must map to one of ‘kegg’, ‘ncbi-geneid’, ‘ncib-proteinid’ or ‘uniprot’ because gseKEGG() only accepts one of these 4 options as it’s keytype parameter.

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids = bitr(significant_genes$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(ids) = c("symbol", "EntrezID")
df2 = merge(significant_genes, ids, by = "symbol")

# we want the log2 fold change 
original_gene_list = df2$logFC

# name the vector
names(original_gene_list) <- df2$EntrezID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

## KEGG

kegg_ctl_vs_2h <- enrichKEGG(
  names(gene_list),
  organism = "mmu",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

# table
t_kegg_ctl_vs_2h <- as.data.frame(kegg_ctl_vs_2h)

## convert gene ID to Symbol
edox <- setReadable(kegg_ctl_vs_2h, 'org.Mm.eg.db', 'ENTREZID')

# table
t_kegg_ctl_vs_2h <- as.data.frame(edox)

## remove (mus musculos name)
edox@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", edox@result$Description, fixed = T)
edox@result$Description

## table
t_kegg_ctl_vs_2h <- as.data.frame(edox)


## plot
NCG_barplot <- barplot(edox, showCategory = c("Protein processing in endoplasmic reticulum", 
                                              "Ribosome")) +
  ggtitle("KEGG: Kyoto Encyclopedia of Genes and Genomes (2 hours)") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

# print
NCG_barplot


## save as image
tiff("barplor_CTL_vs_2h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
NCG_barplot 
dev.off()

## CNETPLOT
cnetplot1_2h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("Protein processing in endoplasmic reticulum", 
                                          "Ribosome")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (2 hours)",
    subtitle = NULL
  )

# Print
cnetplot1_2h

## save
tiff("cnetplot_CTL_vs_2h_April_18_2024.tiff", units="in", width=14, height=12, res=600)
cnetplot1_2h 
dev.off()



###########################################################
###########################################################
###########################################################
########            4 hours             ###################
########                                ###################
###########################################################
###########################################################
###########################################################



## Import dataset
eDat <- read.table("Count.txt", header=TRUE, sep="\t", row.names = 1)
gDat <- read.table("Group.txt", header=TRUE, sep="\t")


## Counts
## filter only samples CTL and Exe2h

# Select only 'CTL' ou 'Exe2h'
selected_columns <- eDat[, grepl("^CTL|^Exe4h", names(eDat))]

# new table
new_table1 <- cbind(rownames = rownames(eDat), selected_columns)

## remove the first colllumn
new_table1$rownames <- NULL

## Treat
new_group1 <- subset(gDat, treat == "CTL" | treat == "Exe4h")

# ATENTION!!
new_group1$treat

# Control treatment needs to appear as the first level for the Treat Vs. Control comparison
new_group1$treat <-factor(new_group1$treat, levels=c("CTL", "Exe4h"))

# Input data
y_adj <- DGEList(counts=new_table1, samples=new_group1, group=new_group1$treat)
design_adj <- model.matrix(~ treat - 1, data = new_group1)

# Filter lowly expressed genes
keep_adj <- rowSums(cpm(y_adj)>1) >= 2
sum(keep_adj)
y_adj <- y_adj[keep_adj, , keep.lib.sizes=FALSE]


## Explanation: ">= 2" is a condition to maintain the gene, which means that the gene needs to have at least two replicas with CPM greater than 1 to be maintained.

## Explanation: after filtering, 16010 genes remained that have sufficient expression in at least two replicates with CPM greater than 1. Genes that do not meet these criteria are removed from the dataset.

# Normalization for RNA composition
y_adj <- calcNormFactors(y_adj)
y_adj$samples
norm.expr_adj <- y_adj$samples

# norm counts
logcpm_adj <- cpm(y_adj, log=TRUE)

##########################################
##      Data analyses with edgeR        ##
##########################################


# Estimate dispersion with Cox-Reid profile-adjusted likelihood (CR) method
y_adj <- estimateDisp(y_adj, design_adj)
y_adj$common.dispersion


## Explanation: A low dispersion value indicates that the expression levels of genes are relatively consistent across replicates, suggesting that there is little variability beyond what is expected due to random sampling. Conversely, a high dispersion value indicates greater variability in gene expression levels among replicates.

## plot
plotBCV(y_adj)

# Fit Generalized Linear Model (GLM) 
fit <- glmFit(y_adj, design_adj)
colnames(fit)

# Differential gene expression analysis
lrt <- glmLRT(fit, contrast=c(-1,1))

## Coefficients used in the contrast for the likelihood ratio test
topTags(lrt)

## Explanation: It indicates that the comparison is between the coefficient for treatCTL (the reference group) and the coefficient for treatExe2h (the group of interest).

## Number of genes identified as differentially expressed (Disregarding parameters)
is.de <- decideTestsDGE(lrt)
summary(decideTestsDGE(lrt))

## plot MD
plotMD(lrt, status=is.de)

# Number of up and down-regulated. FDR < 5% and |log2 fold change| ≥ 0
summary(de <- decideTestsDGE(lrt, lfc=0))

# MA-plot
detags <- rownames(y_adj)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

# To export the results
result <- as.data.frame(topTags(lrt, n=1000000))
result <- add_column(result, symbol= rownames(result), .before = "logFC")

## Check all the genes Downregulated and Upregulated

# Remove unnecessary columns
result$logCPM <- NULL
result$LR <- NULL

# Remove rows with FDR < 0.05 and logFC >= 0
result$expression <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 0,
                            ifelse(result$logFC >= 0, 'Up', 'Down'),
                            'Stable')

# Select only genes that are upregulated or downregulated
result_final <- subset(result, expression %in% c("Up", "Down"))

# Remove row names
rownames(result_final) <- NULL

##########################################
############      PLOTS        ###########
##########################################

## volcano plot

# Extract top tags from lrt
result_edgeR <- as.data.frame(topTags(lrt, n=1000000, adjust.method = "BH"))

# Extract log-transformed counts per million
count_CPM <- as.data.frame(logcpm_adj)

## Convert row names into first column
result_edgeR <- add_column(result_edgeR, gene = rownames(result_edgeR), .before = "logFC")

count_CPM <- add_column(count_CPM, gene = rownames(count_CPM), .before = "CTL_1")

## Order by gene
result_edgeR <- result_edgeR[order(result_edgeR$gene),]
count_CPM <- count_CPM[order(count_CPM$gene),]

## Combine result_edgeR columns with count_CPM
count_CPM <- add_column(count_CPM, logFC = result_edgeR$logFC, .before = "CTL_1")
count_CPM <- add_column(count_CPM, FDR = result_edgeR$FDR, .before = "CTL_1")
count_CPM <- add_column(count_CPM, gene2 = result_edgeR$gene, .before = "CTL_1")

# Check if gene and gene2 are the same and in the same order
if(all(count_CPM$gene == count_CPM$gene2)) {
  print("The values in 'gene' and 'gene2' are the same.")
} else {
  print("The values in 'gene' and 'gene2' are different.")
  different_indices <- which(count_CPM$gene != count_CPM$gene2)
  print(paste("The value in 'gene' is different at indices:", different_indices))
}

## remove gene2
count_CPM$gene2 <- NULL

## Insert -log10 FDR
count_CPM <- add_column(count_CPM, logFDR = -log10(result_edgeR$FDR), .before = "CTL_1")

## Insert -log10 p-value
count_CPM <- add_column(count_CPM, logpvalue = -log10(result_edgeR$PValue), .before = "CTL_1")


## PLOT

# Assign count_CPM to expressao_genica
expressao_genica <- count_CPM

# Determine expression based on FDR and logFC thresholds
expressao_genica$expression <- ifelse(expressao_genica$FDR < 0.05 & abs(expressao_genica$logFC) >= 0,
                                      ifelse(expressao_genica$logFC >= 0, 'Up', 'Down'),
                                      'Stable')

# Count the number of genes categorized as Up and Down
sum(expressao_genica$expression == "Up")
sum(expressao_genica$expression == "Down")

# Create volcano plot

# Filter the data to include only "Down" and "Up" genes
expressao_genica_filtered <- expressao_genica[expressao_genica$expression %in% c("Down", "Up"), ]

# Order the filtered dataframe by logFC
expressao_genica_filtered <- expressao_genica_filtered[order(expressao_genica_filtered$logFC, decreasing = TRUE), ]

# Select the top 5 upregulated genes and the top 5 downregulated genes
top_genes_up <- head(subset(expressao_genica_filtered, expression == "Up"), 5)
top_genes_down <- head(subset(expressao_genica_filtered, expression == "Down"), 5)

# Combine top_genes_up and top_genes_down into a single dataframe
top_genes <- rbind(top_genes_up, top_genes_down)

## see
print(top_genes)

# volcano plot
library(ggrepel)

plot2 <- ggplot(data = expressao_genica,
                aes(x = logFC,
                    y = logFDR,
                    colour = expression)) +
  geom_point(size = 3) +
  geom_label_repel(data = top_genes, aes(label = gene), size = 5, force = 5, 
                   fill = "white", color = "black", box.padding = unit(0.5, "lines")) +  # Labeling top 10 genes with a rectangle
  scale_color_manual(values = c("red", "grey", "#009933"),
                     name = "",
                     breaks = c("Up", "Stable", "Down")) +
  geom_vline(xintercept = c(-0, 0), lty = 4, col = "black", lwd = 0.8) +
  labs(x = "log2FC",
       y = "-log10 (adj.P.Val)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        text = element_text(size = 12)) +
  ggtitle("CTL vs 4 hours")


# To display the volcano plot
dev.off()
plot(rnorm(50), rnorm(50))

plot2

## save as image
tiff("volcano_CTL_vs_4h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
plot2 
dev.off()


## heatmap

## Only genes Down and Upregulated
head(expressao_genica_filtered)

#remove the fisrt five cols
expressao_genica_filtered <- expressao_genica_filtered[, -(1:5)]

## remove a col name expression
expressao_genica_filtered$expression <- NULL

# Convert dataframe to matrix
samp.with.rownames <- as.matrix(expressao_genica_filtered)

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2

library(RColorBrewer)
cols <- brewer.pal(9, "RdYlBu")
pal <- colorRampPalette(cols)

column_annotation <- as.data.frame((colnames(samp.with.rownames)))
rownames(column_annotation) <- column_annotation$`(colnames(samp.with.rownames))`
column_annotation$Group = c(rep("Control", 5), rep("4 hours", 5))
column_annotation$`(colnames(samp.with.rownames))` <- NULL

# Convertendo a coluna "Group" em um fator com os níveis corretos
column_annotation$Group <- factor(column_annotation$Group, levels = c("Control", "4 hours"))

# Verificando os níveis da coluna "Group" após a conversão
levels(column_annotation$Group)


cc <- pheatmap(samp.with.rownames, 
               annotation_col = column_annotation, 
               scale = "row",
               fontsize=12,
               show_rownames = FALSE, 
               annotation_names_row = FALSE,
               annotation_names_col = FALSE,
               annotation_colors = list(Group = c(Control = "#74ADD1", `4 hours` = "#fd61bf")),
               show_colnames = FALSE, 
               col = rev(pal(50)))


## Print 
cc

## save as image
tiff("Heatmap_CTL_vs_4h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
cc 
dev.off()



## pathway

library("clusterProfiler")
library("org.Mm.eg.db")
library("pathview")
library("ggplot2")
library("cowplot")


### Control vs 4 hours ###

## import DEGs
DEGs <- read.delim("~/VanessaPD/differential_expression/DEGs_CTL_vs_4h_NOFC.txt")

## select only genes with FDR <0.05 
significant_genes <- DEGs[DEGs$FDR < 0.05, ]

# For KEGG pathway enrichment using the gseKEGG() function, we need to convert id types. We can use the bitr function for this (included in clusterProfiler). It is normal for this call to produce some messages / warnings.
# In the bitr function, the param fromType should be the same as keyType from the gseGO function above (the annotation source). This param is used again in the next two steps: creating dedup_ids and df2.
# toType: in the bitr function has to be one of the available options from keyTypes(org.Dm.eg.db) and must map to one of ‘kegg’, ‘ncbi-geneid’, ‘ncib-proteinid’ or ‘uniprot’ because gseKEGG() only accepts one of these 4 options as it’s keytype parameter.

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids = bitr(significant_genes$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Create a new dataframe df2 which has the respective entrez IDs for the gene symbols.
colnames(ids) = c("symbol", "EntrezID")
df2 = merge(significant_genes, ids, by = "symbol")

# we want the log2 fold change 
original_gene_list = df2$logFC

# name the vector
names(original_gene_list) <- df2$EntrezID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

## KEGG

kegg_ctl_vs_4h <- enrichKEGG(
  names(gene_list),
  organism = "mmu",
  keyType = "ncbi-geneid",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

# table
t_kegg_ctl_vs_4h <- as.data.frame(kegg_ctl_vs_4h)

## convert gene ID to Symbol
edox <- setReadable(kegg_ctl_vs_4h, 'org.Mm.eg.db', 'ENTREZID')

# table
t_kegg_ctl_vs_4h <- as.data.frame(edox)

## remove (mus musculos name)
edox@result$Description <- gsub(pattern = " - Mus musculus (house mouse)", replacement = "", edox@result$Description, fixed = T)
edox@result$Description

## table
t_kegg_ctl_vs_4h <- as.data.frame(edox)


## plot
NCG_barplot <- barplot(edox, showCategory = c("AGE-RAGE signaling pathway in diabetic complications", 
                                              "Arginine biosynthesis", 
                                              "Terpenoid backbone biosynthesis", 
                                              "Alanine, aspartate and glutamate metabolism", 
                                              "Bile secretion", 
                                              "Aldosterone synthesis and secretion", 
                                              "Parathyroid hormone synthesis")) +
  ggtitle("KEGG: Kyoto Encyclopedia of Genes and Genomes (4 hours)") +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

# Print
NCG_barplot


## save as image
tiff("barplor_CTL_vs_4h_April_16_2024.tiff", units="in", width=10, height=10, res=600)
NCG_barplot 
dev.off()

## CNETPLOT

cnetplot1_4h <- cnetplot(edox, color.params = list(foldChange = gene_list), circular = T, colorEdge = TRUE,
                         showCategory = c("AGE-RAGE signaling pathway in diabetic complications", 
                                          "Arginine biosynthesis", 
                                          "Terpenoid backbone biosynthesis", 
                                          "Alanine, aspartate and glutamate metabolism", 
                                          "Bile secretion", "Aldosterone synthesis and secretion", 
                                          "Parathyroid hormone synthesis")) + 
  scale_colour_gradient2(name = "fold change", low = "darkblue", mid = "white", high = "red") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(
    title = "KEGG (4 hours)",
    subtitle = NULL
  )

# Print
cnetplot1_4h

## Save
tiff("cnetplot_CTL_vs_4h_April_18_2024.tiff", units="in", width=14, height=12, res=600)
cnetplot1_4h 
dev.off()
