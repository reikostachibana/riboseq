library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

get_coldata <- function(counts_df){
  coldata <- data.frame(colnames(counts_df), colnames(counts_df))
  colnames(coldata) <- c("SampleNames", "Condition")
  coldata$Condition <- "Vehicle"
  coldata$Condition[coldata$SampleNames %like% "Ven"] <- "Ven"
  coldata$Population <- "HSC"
  coldata$Population[coldata$SampleNames %like% "GMP"] <- "GMP"
  coldata$LibraryType <- "Input"
  coldata$LibraryType[coldata$SampleNames %like% "RIBO"] <- "Ribo"
  coldata$Population <- factor(coldata$Population, levels = c("GMP", "HSC"))
  coldata$Condition <- factor(coldata$Condition, levels = c("Vehicle", "Ven"))
  coldata$LibraryType <- factor(coldata$LibraryType, levels = c("Input", "Ribo"))
  
  return (coldata)
}

comparisons <- list(
  HSC_Ven_vs_HSC_vehicle = list(
    RIBO = c("RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B"),
    RNA = c("RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_HSC_Ven_A", "RNA_HSC_Ven_B")
  ),
  GMP_Ven_vs_GMP_vehicle = list(
    RIBO = c("RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B"),
    RNA = c("RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")
  ),
  HSC_Ven_vs_GMP_Ven = list(
    RIBO = c("RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B"),
    RNA = c("RNA_HSC_Ven_A", "RNA_HSC_Ven_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")
  ),
  HSC_vehicle_vs_GMP_vehicle = list(
    RIBO = c("RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B"),
    RNA = c("RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B")
  ),
  HSC_Ven_vs_GMP_vehicle = list(
    RIBO = c("RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B", "RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B"),
    RNA = c("RNA_HSC_Ven_A", "RNA_HSC_Ven_B", "RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B")
  )
)

################################################################################
# PCA plot
################################################################################

# ribo_file <- "//wsl$/Ubuntu/home/reiko/riboseq/ribo_counts.txt"
# rna_file <- "//wsl$/Ubuntu/home/reiko/riboseq/rna_counts.txt"
ribo_file <- "/Users/reikotachibana/Documents/Chung Lab/riboseq/ribo_counts.txt"
rna_file <- "/Users/reikotachibana/Documents/Chung Lab/riboseq/rna_counts.txt"
comparison <- c(comparisons$HSC_Ven_vs_HSC_vehicle$RNA,
                comparisons$GMP_Ven_vs_GMP_vehicle$RNA)
# comparisons$HSC_Ven_vs_HSC_vehicle$RNA,
# comparisons$GMP_Ven_vs_GMP_vehicle$RNA)

ribo <- read.delim(ribo_file)
rownames(ribo) <- ribo$gene
ribo <- ribo[, -1]

rna <- read.delim(rna_file)
rownames(rna) <- rna$gene
rna <- rna[, -1]

merge <- cbind(ribo, rna)
subset <- merge[, comparison]

colData <- get_coldata(subset) 
dds <- DESeqDataSetFromMatrix(countData = subset,
                              colData = colData,
                              design = ~ Population + Condition + Population:Condition)

dds <- DESeq(dds)
res <- results(dds)

vsd <- vst(dds, blind=FALSE)
pca_plot <- plotPCA(vsd, intgroup=c("Population", "LibraryType", "Condition"),
                    returnData=TRUE)
pca_plot$Group <- interaction(pca_plot$Population, pca_plot$LibraryType, pca_plot$Condition)

ggplot(pca_plot, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ", round(attr(pca_plot, "percentVar")[1], 2), " variance")) +
  ylab(paste0("PC2: ", round(attr(pca_plot, "percentVar")[2], 2), " variance")) +
  ggtitle("RNA Samples") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size=20),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))