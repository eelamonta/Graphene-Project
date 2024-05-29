# Bulk_RNAseq_Analysis.R

install.packages("BiocManager")
install.packages("pheatmap")
install.packages("Rtsne")
instal.packages("ggridges")
install.packages("ggrepel")
install.packages("ggVennDiagram")
BiocManager::install("tximport")
BiocManager::install("ensembldb")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("tidyverse")
BiocManager::install("biomaRt")
BiocManager::install("sva")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ggVennDiagram)
library(DOSE)
library(tximport)
library(edgeR)
library(tidyr)
library(org.Hs.eg.db)
library(dplyr)
library(Rtsne)

# Set working directory to folder with kallisto output files
setwd("~/Desktop/RNAseq/Alignment/kallisto_aligned")

################## Create a transcript to gene ID matrix ################## 

# Get the transcript to Gene IDs for Ensemble
martGRCh38.108 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                                   dataset = "hsapiens_gene_ensembl", 
                                   host = 'https://jan2020.archive.ensembl.org', #might not be the right version!
                                   path = "/biomart/martservice")
GRCh38.108t2g<- biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                              "ensembl_gene_id"),
                                              mart = martGRCh38.108)

################## Import the kallisto output files with tximport ################## 

# Set the base directory containing your files
base_dir <- "~/Desktop/RNAseq/Alignment/kallisto_aligned"
# Get the samples from a txt file
samples <- read.table(file.path(base_dir, "Sample_info_norgocombined.txt"), header = TRUE, stringsAsFactors=FALSE)
# Remove rows where cell_line is "WT83" and batch "BM1"
samples_H9 <- samples[samples$cell_line != "WT83", ]
samples <- samples_H9[samples_H9$batch != "BM1", ]
# For Kallisto, describe the path to find the quant.sf files
files <- file.path(base_dir, samples$sample_name, "abundance.h5")
# Apply the sample names to "files"
names(files) <- paste0(c(samples$sample_name))
# Check if all files exist
all(file.exists(files))

# Import the abundance/counts measurements using tximport
txi_lsTPM <- tximport(files, 
                      type = "kallisto", 
                      tx2gene = GRCh38.108t2g, 
                      countsFromAbundance = "lengthScaledTPM")

# Get the samples from a txt file
samples1 <- samples[samples$label %in% c("rgo_light", "rgo_nolight"), ]
samples2 <- samples[samples$label %in% c("rgo_light", "norgo"), ]
samples3 <- samples[samples$label %in% c("rgo_nolight", "norgo"), ]

# Create a sample table
sampleTable <- data.frame(condition=factor(rep(c(samples$condition))))
sampleTable1 <- data.frame(condition=factor(rep(c(samples1$condition))))
sampleTable2 <- data.frame(condition=factor(rep(c(samples2$condition))))
sampleTable3 <- data.frame(condition=factor(rep(c(samples3$condition))))

################## Create expression data ################## 

# Create DGEList, y$counts are TPM normalized counts
y <- DGEList(txi_lsTPM$counts,
             lib.size = colSums(txi_lsTPM$counts), # here library size is the sum of counts
             norm.factors = calcNormFactors(txi_lsTPM$counts),
             samples = samples$sample_name,
             group = samples$condition) 

#Create a Homo Sapiens annotation
Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(y$counts),
                                columns=c("ENTREZID","SYMBOL"),
                                keytype="ENSEMBL",
                                multiVals="first")

# Remove duplicated terms
Hs_ann <- Hs_ann[!duplicated(Hs_ann[,1]),]
# Filter out ribosomal RNA (RPL and RPS genes) from the annotation
Hs_ann_filtered <- Hs_ann[!grepl("^RPL|^RPS", Hs_ann$SYMBOL), ]

# Apply the annotation to your limma object "y"
y$genes <- Hs_ann_filtered
y$counts <- y$counts[match(y$genes$ENSEMBL, rownames(y$counts)),]

# View the library size for each sample
y$samples

# Define labels for plot graphics
# Define colors for each condition
condition_colors <- c(
  "0rgo_nolight" = "#FF0000",
  "0.01rgo_nolight" = "#FFDB58",
  "0.1rgo_nolight" = "#FFDB58",
  "0rgo_light" = "#FF0000",
  "0.01rgo_light" = "#00FFFF",
  "0.1rgo_light" = "#00FFFF")
# Define shapes for each batch/biological replicate
batch_shapes <- c("B33" = 16, "B38" = 17)
# Create the group factor and specify levels
group <- factor(samples$condition, levels = c("0rgo_light", "0rgo_nolight", "0.01rgo_light", "0.01rgo_nolight", "0.1rgo_light", "0.1rgo_nolight"))
group1 <- factor(samples$batch, levels = c("B33", "B38"))
group2 <- factor(samples$label, levels = c("norgo_light", "norgo_nolight", "rgo_light", "rgo_nolight"))

################## Plot Unfiltered Data ################## 

# Plot the density of unfiltered gene expression for all samples within groups
unfilteredExpr<-y$counts
plotDensities(unfilteredExpr, group=samples$group, col=condition_colors,legend = FALSE)
boxplot(unfilteredExpr, las=2,col=condition_colors, main="")

# Plot MDS of unfiltered gene expression
MDS <- plotMDS(log(unfilteredExpr+1), pch = batch_shapes, col = condition_colors, labels = NULL)
legend('topright', col=condition_colors, legend=levels(group), pch = 16, cex = 0.7)
legend('topleft', legend=levels(group1), pch = batch_shapes, cex = 0.7)

# Plot PCA of unfiltered gene expression
pca_result <- prcomp(t((log(unfilteredExpr+1))))
pca_scores <- as.data.frame(pca_result$x)
pca_scores$group <- group
ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = group)) +
  scale_color_manual(values = condition_colors, guide=FALSE) +
  geom_point() +
  theme(legend.position = "none") +
  theme_minimal()

# Plot t-SNE of unfiltered gene expression
tsne_result <- Rtsne(as.matrix(t(log(unfilteredExpr+1))),perplexity = 2)
tsne_df <- as.data.frame(cbind(tsne_result$Y,samples))
ggplot(tsne_df, aes(x = `1`, y = `2`,color=condition,shape=batch)) +
  scale_color_manual(values = condition_colors) +
  geom_point() + 
  labs(x = "t-SNE Component 1", y = "t-SNE Component 2")+
  theme_minimal()

################## Filter Method #1 - edgeR and TMM ##################

# Filter lowly expressed genes via edgeR
#keep = filterByExpr(y)

# Normalize  with TMM
# Adjusts library sizes based on the assumption that most genes are not differentially expressed
#y <- calcNormFactors(y, method = "TMM")

# Plot the density of filtered gene expression for all samples within groups
#filteredExpr <- cpm(y, log=T)
#plotDensities(filteredExpr, group=samples$group, col=condition_colors,legend = FALSE)
#boxplot(filteredExpr, las=2, main="", col=condition_colors)

################## Filter Method #2 - CPM and ComBat ##################

# Calculate logCPM values
dge_logcpm <- cpm(y, log = TRUE)

# Create a design matrix with an intercept term
modcombat <- model.matrix(~1, data = samples)

# Option #1: Apply ComBat_seq for batch correction
#correctedExpr_seq <-sva::ComBat_seq(counts = y$counts, batch = samples$batch)
#correctedExpr <- correctedExpr_seq

# Option #2: Apply ComBat for batch correction
correctedExpr <-sva::ComBat(dat = dge_logcpm, batch = samples$batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)


################## Plot Filtered Expression Data ##################

# Plot the density of filtered gene expression for all samples within groups
plotDensities(correctedExpr, group=samples$group, col=condition_colors,legend = FALSE)
boxplot(correctedExpr, las=2,col=condition_colors, main="")

# Plot MDS of filtered gene expression
MDS <- plotMDS(correctedExpr, pch = batch_shapes, col = condition_colors, labels = NULL)
legend('topright', col=condition_colors, legend=levels(group), pch = 16, cex = 0.7)
legend('topleft', legend=levels(group1), pch = batch_shapes, cex = 0.7)

# Plot PCA of filtered gene expression
pca_result <- prcomp(t(correctedExpr))
pca_scores <- as.data.frame(pca_result$x)
pca_scores$group <- group
ggplot(data = pca_scores, aes(x = PC1, y = PC2, color = group)) +
  scale_color_manual(values = condition_colors, guide=FALSE) +
  geom_point() +
  theme(legend.position = "none") +
  theme_minimal()

# Plot t-SNE of filtered gene expression
tsne_result <- Rtsne(as.matrix(t(correctedExpr)),perplexity = 2)
tsne_df <- as.data.frame(cbind(tsne_result$Y,samples))
ggplot(tsne_df, aes(x = `1`, y = `2`,color=condition,shape=batch)) +
  scale_color_manual(values = condition_colors) +
  geom_point() + 
  labs(x = "t-SNE Component 1", y = "t-SNE Component 2")+
  theme_minimal()

# Optional: Visualize batch corrected data
boxplot(as.data.frame(unfilteredExpr),main="Original")
boxplot(as.data.frame(correctedExpr),main="Batch corrected")

# Optional: Save correctedExpr_Seq as a CSV
#write.csv(correctedExpr_seq, file = "correctedExpr_seq.csv")

################## Make Design Matrix ##################

# Create the group factor and specify levels
label <- factor(samples$label, levels = c("norgo", "rgo_nolight", "rgo_light"))
label1 <- factor(samples1$label, levels = c("rgo_light", "rgo_nolight"))
label2 <- factor(samples2$label, levels = c("rgo_light", "norgo"))
label3 <- factor(samples3$label, levels = c("rgo_nolight", "norgo"))

# Create the design model matrix with an intercept term
design <- model.matrix(~0+label, data = sampleTable)
design1 <- model.matrix(~0+label1, data = sampleTable1)
design2 <- model.matrix(~0+label2, data = sampleTable2)
design3 <- model.matrix(~0+label3, data = sampleTable3)

# Rename the colnames using levels of 'group'
colnames(design) <- levels(label)
colnames(design1) <- levels(label1)
colnames(design2) <- levels(label2)
colnames(design3) <- levels(label3)

# Create the contrast matrices
contr.matrix <- makeContrasts(comp1 = rgo_light - rgo_nolight,
                              comp2 = rgo_light - norgo,
                              comp3 = rgo_nolight - norgo,
                              levels = colnames(design))
contr.matrix1 <- makeContrasts(comp1 = rgo_light - rgo_nolight, levels = colnames(design1))
contr.matrix2 <- makeContrasts(comp2 = rgo_light - norgo, levels = colnames(design2))
contr.matrix3 <- makeContrasts(comp3 = rgo_nolight - norgo, levels = colnames(design3))

################## Perform Differential Expression Analysis ##################

# Remove NA counts from y
#y <- y[!is.na(y$counts),]

# Run voom
v <- voom(y, design=design, plot=TRUE)
v1 <- voom(y[,grepl("0_01rgo-Light|0_1rgo-Light|0_01rgo-NoLight|0_1rgo-NoLight",colnames(y))], design=design1, plot=TRUE)
v2 <- voom(y[,grepl("0_01rgo-Light|0_1rgo-Light|0rgo",colnames(y))], design=design2, plot=TRUE)
v3 <- voom(y[,grepl("0_01rgo-NoLight|0_1rgo-NoLight|0rgo",colnames(y))], design=design3, plot=TRUE)

# Fit the linear model
fit <- lmFit(v, design)
fit1 <- lmFit(v1, design1)
fit2 <- lmFit(v2, design2)
fit3 <- lmFit(v3, design3)

# Apply your contrasts
cfit <- contrasts.fit(fit, contrasts=contr.matrix)
cfit1 <- contrasts.fit(fit1, contrasts=contr.matrix1)
cfit2 <- contrasts.fit(fit2, contrasts=contr.matrix2)
cfit3 <- contrasts.fit(fit3, contrasts=contr.matrix3)

# eBayes method and plot mean-variance trend
efit <- eBayes(cfit)
plotSA(efit)
efit1 <- eBayes(cfit1)
plotSA(efit1)
efit2 <- eBayes(cfit2)
plotSA(efit2)
efit3 <- eBayes(cfit3)
plotSA(efit3)

# See how many genes are differentially expressed
summary(decideTests(efit))
summary(decideTests(efit1))
summary(decideTests(efit2))
summary(decideTests(efit3))

# Look only at DEGs that have at least an absolute logFC of 0.6. This recalculates t and FDR at a particular fold change cutoff.
tfit <- treat(efit,lfc=log2(0.6))
tfit1 <- treat(efit1,lfc=log2(0.6))
tfit2 <- treat(efit2,lfc=log2(0.6))
tfit3 <- treat(efit3,lfc=log2(0.6))

# See how many DEGs in each condition
dt <- decideTests(tfit)
dt1 <- decideTests(tfit1)
dt2 <- decideTests(tfit2)
dt3 <- decideTests(tfit3)

summary(dt)
summary(dt1)
summary(dt2)
summary(dt3)

# Create a topTable of all the limma results for a particular comparison
reslimma <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
reslimma1 <- topTable(efit1, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
reslimma2 <- topTable(efit2, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
reslimma3 <- topTable(efit3, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

################## Plots of t-test Analysis ##################

expression_matrix<-t(correctedExpr)

t_test_function<-function(x,y) {
  res <- try(t.test(x, y)$p.value)
  if(grepl(pattern = "Error", x = res)){ 
    return(1)
  } else {
    return(res)
  }
}

# Select one of these class indices for further t-test analysis!!!
# comparison1
class1_indices<-grep("0_01rgo-Light|0_1rgo-Light",rownames(expression_matrix))
class2_indices<-grep("0_01rgo-NoLight|0_1rgo-NoLight",rownames(expression_matrix))
# comparison2
class1_indices<-grep("0_01rgo-Light|0_1rgo-Light",rownames(expression_matrix))
class2_indices<-grep("0rgo",rownames(expression_matrix))
# comparison3
class1_indices<-grep("0_01rgo-NoLight|0_1rgo-NoLight",rownames(expression_matrix))
class2_indices<-grep("0rgo",rownames(expression_matrix))
# comparison4
class1_indices<-grep("0_01rgo-Light|0_1rgo-Light",rownames(expression_matrix))
class2_indices<-grep("0rgo|0_01rgo-NoLight|0_1rgo-NoLight",rownames(expression_matrix))

# Perform t_tests with indices and t-test function 
t_test_results<-parallel::mclapply(1:ncol(expression_matrix),function(x) t_test_function(expression_matrix[class1_indices,x],expression_matrix[class2_indices,x]))
t_test_results<-unlist(t_test_results)
names(t_test_results)<-colnames(expression_matrix)
# Create log fold change of the t-test results
log_fold_change<-parallel::mclapply(1:ncol(expression_matrix),function(x) log2(mean(expression_matrix[class1_indices,x])/mean(expression_matrix[class2_indices,x])),mc.cores=10)
log_fold_change<-unlist(log_fold_change)
names(log_fold_change)<-colnames(expression_matrix)
t_test_results<-t_test_results[!is.na(t_test_results)]
sum(t_test_results<.05)
# Optional: BH correction
#t_test_results_bh<-p.adjust(t_test_results,method="BH")
#t_test_results_bh[t_test_results_bh<.05]

## Heatmap of top 100 DEGs
top_100_genes<-t_test_results[order(t_test_results)][1:100]
for_heatmap<-expression_matrix[c(class1_indices,class2_indices),names(top_100_genes)]
pheatmap(scale(for_heatmap),cluster_rows=T,cluster_cols=T,show_rownames=T,show_colnames=FALSE)

################## Plots of Limma DEG Comparisons ##################

## Volcano Plot
# Change which table is set to input!!!
to_plot <- reslimma1

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
to_plot$diffexpressed <- "NO"
# If log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
to_plot$diffexpressed[to_plot$logFC > 0.6 & to_plot$P.Value < 0.05] <- "UP"
# If log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
to_plot$diffexpressed[to_plot$logFC < -0.6 & to_plot$P.Value < 0.05] <- "DOWN"
to_plot<-to_plot[to_plot$logFC!="NaN",]
to_plot<-to_plot[to_plot$logFC!="Inf",]
to_plot<-to_plot[to_plot$logFC!="-Inf",]
important_genes <- to_plot[to_plot$P.Value < 0.05 & abs(to_plot$logFC) > 0.6,]
important_genes <- important_genes[important_genes$SYMBOL != "NA",]

# Select for genes of interest
sig_genes <- to_plot %>% filter(SYMBOL =="ND1"|SYMBOL == "RARRES2" | SYMBOL == "GNAT1" |SYMBOL == "RGR" | SYMBOL == "ALDH3A1"|SYMBOL == "FJX1"|SYMBOL == "FJX1"|SYMBOL == "RBL1"|SYMBOL == "PRPH2"|ENSEMBL == "ENSG00000273944"|ENSEMBL == "ENSG00000166923"|ENSEMBL == "ENSG00000197921" |ENSEMBL == "ENSG00000262102" |SYMBOL == "FABP5" |SYMBOL == "GTF2IP14"  |ENSEMBL == "ENSG00000230678" |SYMBOL == "ND2" |SYMBOL == "RAB26" |SYMBOL == "ACVR1" |SYMBOL == "SLITRK2" |ENSEMBL == "ENSG00000282283" |SYMBOL == "RSPH3" |ENSEMBL == "ENSG00000262102" |SYMBOL == "ACVR1"  |ENSEMBL == "ENSG00000137312" )

# Create volcano plot with ggplot
options(ggrepel.max.overlaps = Inf)
ggplot(data = to_plot, aes(x = logFC, y = -log10(P.Value), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 5.5), xlim = c(-8, 8)) + 
  labs(color = "", x = expression("log"[2]*"FC"), y = expression("-log"[10]*"Pvalue")) + 
  geom_label_repel(data = sig_genes, # Add labels last to appear as the top layer  
                   aes(label = SYMBOL),
                   segment.color="black", 
                   force = 2,
                   nudge_y = 1) +
  geom_point(data = sig_genes,
             shape = 21,
             size = 2, 
             fill = "yellow", 
             colour = "black") + 
  ggtitle("") +
  theme_minimal() + theme(legend.position="none")

## MA plot
# Exclude if pvalue < 0.05
to_plot2 <- to_plot[to_plot$P.Value < 0.05,]
to_plot2$diffcounts <- "NO"
# If log2Foldchange > 0.6 and average expression > 3, set as "UP"
to_plot2$diffcounts[to_plot2$logFC > 0.6 & to_plot2$AveExpr > 3 ] <- "UP"
# If log2Foldchange < 0.6 and average expression > 3, set as "DOWN"
to_plot2$diffcounts[to_plot2$logFC < -0.6 & to_plot2$AveExpr > 3 ] <- "DOWN"
to_plot<-to_plot[to_plot$AveExpr!="NaN",]
to_plot<-to_plot[to_plot$AveExpr!="Inf",]
to_plot<-to_plot[to_plot$AveExpr!="-Inf",]

# Create MA plot with ggplot
ggplot(data = to_plot2, aes(x = AveExpr, y = logFC, color = diffcounts)) +
  geom_vline(xintercept = c(3), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  labs(color = "", x = expression("log"["e"]*"Mean Expression"), y = expression("log"[2]*"FC")) + 
  geom_label_repel(data = sig_genes,  
                   aes(label = SYMBOL),
                   segment.color="black", 
                   force = 2,
                   nudge_y = 1) +
  geom_point(data = sig_genes,
             shape = 21,
             size = 2, 
             fill = "yellow", 
             colour = "black") + 
  ggtitle("") +
  theme_minimal() + theme(legend.position="none")

## Venn Diagram

# Calculate number of DEGs and create lists of DEGs for each comparison
to_plot <- reslimma1
to_plot_sort <- to_plot[order(to_plot$logFC, decreasing = TRUE), ]
top <- to_plot_sort[to_plot_sort$P.Value < 0.05 & abs(to_plot_sort$logFC) > 0.6, "ENSEMBL"]
genes_to_test1 <- unique(top)
length(genes_to_test1)
to_plot <- reslimma2
to_plot_sort <- to_plot[order(to_plot$logFC, decreasing = TRUE), ]
top <- to_plot_sort[to_plot_sort$P.Value < 0.05 & abs(to_plot_sort$logFC) > 0.6, "ENSEMBL"]
genes_to_test2 <- unique(top)
length(genes_to_test2)
to_plot <- reslimma3
to_plot_sort <- to_plot[order(to_plot$logFC, decreasing = TRUE), ]
top <- to_plot_sort[to_plot_sort$P.Value < 0.05 & abs(to_plot_sort$logFC) > 0.6, "ENSEMBL"]
genes_to_test3 <- unique(top)
length(genes_to_test3)

ven<-list(comp1 = genes_to_test1, comp2 = genes_to_test2, comp3  = genes_to_test3)

# Plot venn diagram with ggVennDiagram
ggVennDiagram(ven, label_alpha = 0) + 
  ggplot2::scale_fill_gradient(low="white",high = "darkgrey")

# Look more specifically at unique DEGs in each comparison
DEGs_of_interest1 <- list()
for (i in genes_to_test1){
  if (!(i %in% genes_to_test2)){
    if (!(i %in% genes_to_test3)){
      DEGs_of_interest1<- c(DEGs_of_interest1,i)
    }
  }
}  
DEGs_of_interest2 <- list()
for (i in genes_to_test2){
  if (!(i %in% genes_to_test1)){
    if (!(i %in% genes_to_test3)){
      DEGs_of_interest2<- c(DEGs_of_interest2,i)
    }
  }
} 
DEGs_of_interest3 <- list()
for (i in genes_to_test3){
  if (!(i %in% genes_to_test2)){
    if (!(i %in% genes_to_test1)){
      DEGs_of_interest3<- c(DEGs_of_interest3,i)
    }
  }
} 
DEGs_of_interest4 <- list()
for (i in genes_to_test1){
  if (i %in% genes_to_test2){
    if (!(i %in% genes_to_test3)){
      DEGs_of_interest4<- c(DEGs_of_interest4,i)
    }
  }
}
DEGs_of_interest5 <- list()
for (i in genes_to_test1){
  if (i %in% genes_to_test2){
    if (i %in% genes_to_test3){
      DEGs_of_interest5<- c(DEGs_of_interest5,i)
    }
  }
}

################## Pathway Analysis Plots of DEGs ##################

## GO Terms
# Change which table is set to input!!!
input <- reslimma1
# only look at genes with pvalue < 0.05 and log2foldchange > 0.6
genes_to_test <- rownames(input[input$logFC>0.6 & input$P.Value<.05,])
GO_results <- enrichGO(gene =genes_to_test,
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENSEMBL", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       readable = TRUE)
as.data.frame(GO_results)

# Create GO term plots
plot(barplot(GO_results, showCategory = 20))
plot(dotplot(GO_results, showCategory=20))
edo <- pairwise_termsim(GO_results)
plot(emapplot(edo, cex_label_category= 0.75,cex_category = 1, showCategory=25))
# create table of GO term results
cluster_summary <- data.frame(GO_results@result$Description, GO_results@result$geneID, GO_results@result$Count, GO_results@result$pvalue)

## 2-column GO Term Dot Plot
# Look closer at genes that are specifically upregulated or downregulated
genes_up <- rownames(input[input$logFC>0.6 & input$P.Value<.05,])
genes_down <- rownames(input[input$logFC<c(-0.6) & input$P.Value<.05,])
GO_results_up <- enrichGO(gene =genes_up,
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENSEMBL", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       readable = TRUE)
as.data.frame(GO_results_up)

GO_results_down <- enrichGO(gene = genes_down,
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENSEMBL", 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)
as.data.frame(GO_results_down)

# Create table of most upregulated and downregulated GO terms
GO_data_up <- data.frame(GO_results_up@result$Description, GO_results_up@result$ID, GO_results_up@result$GeneRatio, GO_results_up@result$Count, GO_results_up@result$pvalue, GO_results_up@result$p.adjust)
colnames(GO_data_up) <- c("Description", "GO.ID", "GeneRatio", "Count", "pvalue", "adj.p")
GO_data_up$label <- "rgo_light"
GO_data_up<- GO_data_up[order(GO_data_up$pvalue),]
GO_data_down <- data.frame(GO_results_down@result$Description, GO_results_down@result$ID, GO_results_down@result$GeneRatio, GO_results_down@result$Count, GO_results_down@result$pvalue, GO_results_down@result$p.adjust)
colnames(GO_data_down) <- c("Description", "GO.ID", "GeneRatio", "Count", "pvalue", "adj.p")
GO_data_down$label <- "rgo_nolight"
GO_data_down<- GO_data_down[order(GO_data_down$pvalue),]
GO_data <- rbind(GO_data_up, GO_data_down)

# Select a collection of GO terms to plot
GO_data_inflam <- GO_data[grep("innate|inflamma|interfer|amyloid", GO_data$Description, ignore.case = TRUE), ]
GO_data_migration <- GO_data[grep("migration|formation", GO_data$Description, ignore.case = TRUE), ]
GO_data_synapse <- GO_data[grep("assembly|maturation|synap|axon|dentritic|dendrite|neurotrans", GO_data$Description, ignore.case = TRUE), ]
GO_data_retina <- GO_data[grep("light|photo|retina", GO_data$Description, ignore.case = TRUE), ]
GO_data_immune <- GO_data[grep("inflamma|leuko|innate|adaptive|toxic", GO_data$Description, ignore.case = TRUE), ]
GO_data_development <- GO_data[grep("develop|different|formation", GO_data$Description, ignore.case = TRUE), ]
GO_data_oxphos <- GO_data[grep("oxidation|ATP|respiration|electron", GO_data$Description, ignore.case = TRUE), ]
GO_data_ion <- GO_data[grep("channel|calcium|messenger|signaling", GO_data$Description, ignore.case = TRUE), ]

# Enter GO term selection here!!!
GO_data_selection <- GO_data_development

# Select top 100 DEG terms
GO_data_sort <- subset(GO_data_selection, pvalue < 0.05)
GO_data_sort <- GO_data_sort[order(GO_data_sort$pvalue),]
GO_data_filtered <- GO_data_sort[1:100,]

# Add comparison labels to DEG lists
label1 <- "rgo_light" 
label2 <- "rgo_nolight"
final_df<-GO_data_filtered
for (i in 1:nrow(GO_data_filtered)){
  subset<-GO_data_filtered[i,]
  if (nrow(subset)<2 & subset$label == label1){
    new_row <- data.frame(Description = GO_data_filtered$Description[i], GO.ID = GO_data_filtered$GO.ID[i], GeneRatio = 0, Count = 0, pvalue = 1, adj.p = 1, label = label2)
    final_df <- rbind(final_df, new_row)
  }
  else if (nrow(subset)<2 & subset$label == label2){
    new_row <- data.frame(Description = GO_data_filtered$Description[i], GO.ID = GO_data_filtered$GO.ID[i], GeneRatio = 0, Count = 0, pvalue = 1, adj.p = 1, label = label1)
    final_df <- rbind(final_df, new_row)
  }
}

# Optional: Select the top 20 most significant GO terms
#final_df <- final_df[!is.na(final_df$Description),]
#final_df <- final_df[order(final_df$adj.p, increasing = TRUE), ]
#final_df <- final_df[1:20,]

# Create dotplot with ggplot
ggplot(data = final_df, aes(x = label, y = Description, color = -log10(pvalue), size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 10)) +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")

## GSEA Terms
# Change which table is set to input!!!
input <- reslimma1
original_gene_list <- input$logFC
names(original_gene_list) <- input$ENSEMBL
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

# Create GSEA term plots
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
emapplot(gse, showCategory = 10)
ridgeplot(gse) + labs(x = "enrichment distribution")

################## Boxplots of Gene Experssion ##################

# Make new table for gene expression boxplots
for_boxplot<-t(expression_matrix)
row_names <- rownames(for_boxplot)
for_boxplot <- data.frame(ENSEMBL = row_names, for_boxplot, row.names = NULL)
for_boxplot_new <- for_boxplot %>% gather(key = "sample_name", value = "value", -ENSEMBL)
for_boxplot_new$sample_name <- lapply(for_boxplot_new$sample_name, function(x) gsub("\\.", "-", x))
for_boxplot_new$condition <- samples$condition[match(for_boxplot_new$sample_name, samples$sample_name)]
for_boxplot_new$rgo <- samples$rgo[match(for_boxplot_new$sample_name, samples$sample_name)]
for_boxplot_new$light <-samples$light[match(for_boxplot_new$sample_name, samples$sample_name)]
for_boxplot_new$symbol <- y$genes$SYMBOL[match(for_boxplot_new$ENSEMBL, y$genes$ENSEMBL)]
for_boxplot_new$label <- paste(for_boxplot_new$rgo, for_boxplot_new$light, sep = "_")
for_boxplot_new$old_label <-samples$label[match(for_boxplot_new$sample_name, samples$sample_name)]

# Specifications for boxplot significance
label_order <- c("norgo_nolight", "norgo_light", "rgo_nolight", "rgo_light")
for_boxplot_new$label <- factor(for_boxplot_new$label, levels = label_order)
df_filtered <- for_boxplot_new %>%
  group_by(label) %>%
  filter(n() > 1) %>%
  ungroup()
results_pairs <- combn(unique(as.character(for_boxplot_new$label)), 2, simplify = F)

# Select for genes of interest
sig_genes <- for_boxplot_new %>% filter(symbol =="ND1"|symbol == "RARRES2" | symbol == "GNAT1" |symbol == "RGR" | symbol == "ALDH3A1"|symbol == "FJX1"|symbol == "FJX1"|symbol == "RBL1"|symbol == "PRPH2"|ENSEMBL == "ENSG00000273944"|ENSEMBL == "ENSG00000166923"|ENSEMBL == "ENSG00000197921" |ENSEMBL == "ENSG00000262102" |symbol == "FABP5" |symbol == "GTF2IP14"  |ENSEMBL == "ENSG00000230678" |symbol == "ND2" |symbol == "RAB26" |symbol == "ACVR1" |symbol == "SLITRK2" |ENSEMBL == "ENSG00000282283" |symbol == "RSPH3" |ENSEMBL == "ENSG00000262102" |symbol == "ACVR1"  |ENSEMBL == "ENSG00000137312" )

# Create boxplots of gene expression using ggplot
ggplot(data = sig_genes, aes(x = label, y = value, fill = label)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values=c("#FF0000","#FFDB58","#93C572","#16E2F5"))+
  ggtitle(" ") +
  facet_wrap(~ symbol, scales = "free_y") +
  ylab("Expression") +
  xlab("Condition") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
#stat_compare_means(comparisons=results_pairs, label = "p.signif", method="wilcox.test")
#stat_compare_means(comparisons=results_pairs, method="kruskal.test")
