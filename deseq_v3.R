library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(topGO)
library(org.Mm.eg.db)

get_coldata <- function(counts_df){
  coldata <- data.frame(colnames(counts_df), colnames(counts_df))
  colnames(coldata) <- c("SampleNames", "Condition")
  coldata$Condition <- "Vehicle"
  coldata$Condition[coldata$SampleNames %like% "Ven"] <- "Ven"
  coldata$Population <- "HSC"
  coldata$Population[coldata$SampleNames %like% "GMP"] <- "GMP"
  coldata$LibraryType <- "Input"
  coldata$LibraryType[coldata$SampleNames %like% "RIBO"] <- "Ribo"
  # coldata$Length <- 28
  # coldata$Length[coldata$SampleNames %like% "1"] <- 21
  coldata$Population <- factor(coldata$Population, levels = c("GMP", "HSC"))
  coldata$Condition <- factor(coldata$Condition, levels = c("Vehicle", "Ven"))
  coldata$LibraryType <- factor(coldata$LibraryType, levels = c("Input", "Ribo"))
  # coldata$Length <- factor(coldata$Length, levels = c(28, 21))
  
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
  )
)

comparison <- c(comparisons$HSC_vehicle_vs_GMP_vehicle$RIBO, 
                comparisons$HSC_vehicle_vs_GMP_vehicle$RNA)
plot_title <- "HSC_vehicle_vs_GMP_vehicle"

################################################################################
# Run TE DESeq
################################################################################

# Get RPF and RNA counts
ribo_file <- "/Users/reikotachibana/Documents/ChungLab/riboseq/ribo_counts.txt"
rna_file <- "/Users/reikotachibana/Documents/ChungLab/riboseq/rna_counts.txt"
# ribo_file <- "//wsl$/Ubuntu/home/reiko/riboseq/ribo_counts.txt"
# rna_file <- "//wsl$/Ubuntu/home/reiko/riboseq/rna_counts.txt"

ribo <- read.delim(ribo_file)
rownames(ribo) <- ribo$transcript
ribo <- ribo[, -1]
  
rna <- read.delim(rna_file)
rownames(rna) <- rna$gene
rna <- rna[, -1]

# Filter RPF and RNA counts
ribo <- ribo[rowSums(ribo >= 10) >= 4, ] # 8400 genes
rna <- rna[rowSums(rna >= 10) >= 2, ]
common_genes <- intersect(rownames(ribo), rownames(rna))
ribo <- ribo[common_genes, ] # 8395 genes 
rna <- rna[common_genes, ]

merge <- cbind(ribo, rna)
subset <- merge[, comparison]

colData <- get_coldata(subset)
if (length(unique(as.character(colData$Population))) > 1 &
    length(unique(as.character(colData$LibraryType))) > 1) {
  dds <- DESeqDataSetFromMatrix(countData = subset,
                                colData = colData,
                                design = ~ Population + LibraryType + Population:LibraryType)
} else if (length(unique(as.character(colData$Condition))) > 1 &
      length(unique(as.character(colData$LibraryType))) > 1){
  dds <- DESeqDataSetFromMatrix(countData = subset,
                                colData = colData,
                                design = ~ Condition + LibraryType + Condition:LibraryType)
} else if (length(unique(as.character(colData$Length))) > 1 &
           length(unique(as.character(colData$LibraryType))) > 1){
  dds <- DESeqDataSetFromMatrix(countData = subset,
                                colData = colData,
                                design = ~ Length + LibraryType + Length:LibraryType)
}

dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef=4,res=res,type="apeglm")
summary(res)
# write.csv(res,
#           paste0("/Users/reikotachibana/Documents/Chung Lab/riboseq/output/DESeq_TE_", plot_title, ".csv"),
#           row.names = TRUE)

# EnhancedVolcano(res,
#                 lab = res$gene,
#                 x = 'log2FoldChange',
#                 y = 'pvalue')

################################################################################
# logFC RPF vs mRNA
################################################################################

# Ribo
ribo_comparison <- comparison[comparison %like% "RIBO"]
ribo <- ribo[, ribo_comparison]
# ribo <- ribo[rownames(ribo) %in% rownames(filtered), ]
coldata_ribo <- get_coldata(ribo)

if (length(unique(as.character(coldata_ribo$Population))) > 1) {
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo,
                                        design = ~ Population)
} else if (length(unique(as.character(coldata_ribo$Condition))) > 1){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo,
                                        design = ~ Condition)
} else if (length(unique(as.character(coldata_ribo$Condition))) > 1 &
           length(unique(as.character(coldata_ribo$Population))) > 1){
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = ribo,
                                       colData = coldata_ribo, 
                                       design = ~ Condition + Population + Condition:Population)
} else if (length(unique(as.character(coldata_ribo$Length))) > 1){
  ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                        colData = coldata_ribo,
                                        design = ~ Length)
}

ddsMat_ribo <- DESeq(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo)
res_ribo <- lfcShrink(ddsMat_ribo, coef=2,res=res_ribo,type="apeglm") 
# write.csv(res_ribo,
#           paste0("/Users/reikotachibana/Documents/Chung Lab/riboseq/output/DESeq_ribo_", plot_title, ".csv"),
#           row.names = TRUE)

# Get
# ribo_df <- data.frame(Ribo_log_baseMean = log10(res_ribo$baseMean),
#                       Ribo_logFC = res_ribo$log2FoldChange)
# ggplot(ribo_df, aes(x=Ribo_log_baseMean, y=Ribo_logFC)) +
#   geom_point() +
#   labs(title = paste0(plot_title))
#   ylim(c(-.1, .1))

# res_ribo$gene <- genes
# EnhancedVolcano(res_ribo,
#                 lab = res_ribo$gene,
#                 x = 'log2FoldChange',
#                 y = 'pvalue')

# RNA 
rna_comparison <- comparison[comparison %like% "RNA"]
rna <- rna[, rna_comparison]
# rna <- rna[rownames(rna) %in% rownames(filtered), ]
coldata_rna <- get_coldata(rna)

if (length(unique(as.character(coldata_rna$Population))) > 1) {
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, 
                                       design = ~ Population)
} else if (length(unique(as.character(coldata_rna$Condition))) > 1) {
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, 
                                       design = ~ Condition)
} else if (length(unique(as.character(coldata_rna$Condition))) > 1 &
           length(unique(as.character(coldata_rna$Population))) > 1){
  ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                       colData = coldata_rna, 
                                       design = ~ Condition + Population + Condition:Population)
}

ddsMat_rna <- DESeq(ddsMat_rna)
res_rna <- results(ddsMat_rna)
res_rna <- lfcShrink(ddsMat_rna, coef=2,res=res_rna,type="apeglm")
# write.csv(res_rna,
#           paste0("/Users/reikotachibana/Documents/Chung Lab/riboseq/output/DESeq_rna_", plot_title, ".csv"),
#           row.names = TRUE)

# input_df <- data.frame(Input_baseMean = res_rna$baseMean,
#                        Input_logFC = res_rna$log2FoldChange)
# ggplot(input_df, aes(x=Input_baseMean, y=Input_logFC)) + 
#   geom_point()

logFC <- data.frame(Gene = rownames(ribo),
                    Ribo_logFC = res_ribo$log2FoldChange,
                    Ribo_padj = res_ribo$padj, 
                    Input_logFC = res_rna$log2FoldChange,
                    Input_padj = res_rna$padj,
                    TE_logFC = res$log2FoldChange,
                    TE_padj = res$padj,
                    Group = "No change")
rownames(logFC) <- rownames(res_rna)
logFC$Group <- ifelse(res$padj > 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05, 
                      "Forwarded", logFC$Group)
logFC$Group <- ifelse(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj > 0.05,
                      "Exclusive", logFC$Group)
logFC$Group <- ifelse(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05 & 
                        res$log2FoldChange*res_rna$log2FoldChange > 0,
                      "Intensified", logFC$Group)
logFC$Group <- ifelse(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05 &
                        res$log2FoldChange*res_rna$log2FoldChange < 0,
                      "Buffered", logFC$Group)
logFC$Group <- ifelse(res$padj < 0.05 & res_ribo$padj > 0.05 & res_rna$padj < 0.05,
                      "Buffered (Special)", logFC$Group)

logFC <- na.omit(logFC)

# write.csv(logFC, paste0("/Users/reikotachibana/Documents/Chung Lab/riboseq/output/", plot_title, "_results.csv"),
#           row.names = FALSE)

ggplot(logFC, aes(x=Input_logFC, y=Ribo_logFC, color = Group)) +
  geom_point() +
  geom_text(data = logFC[(logFC$Group == "Forwarded" | logFC$Group == "Exclusive" |
                            logFC$Group == "Intensified") | logFC$Group == "Buffered" |
                           logFC$Group == "Buffered (Special)", ],
            aes(label = Gene), vjust = -0.5, check_overlap = TRUE) +
  scale_color_manual(values = c("Forwarded" = "blue",
                                "Exclusive" = "darkgreen",
                                "Intensified" = "red",
                                "Buffered" = "magenta",
                                "Buffered (Special)" = "orange",
                                "No change" = "grey")) +
  labs(title = plot_title) +
  theme(legend.position = "right",
        plot.title = element_text(size=20),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
#   ylim(-15, 4) +
#   xlim(-15, 4)

# res[res$gene %like% "Atf4", ]
# res[res$gene %like% "G3bp1", ]
# res[res$gene %like% "Mecom", ]
# res[res$gene %like% "Tgfb1", ]
# res[res$gene %like% "Six", ]
# res[res$gene %like% "Erg", ]
# res[res$gene %like% "Adam9", ]
# res[res$gene %like% "Akt3", ]
# res_temp <- results(dds, name = "PopulationHSC.LibraryTypeRibo")
# res_temp[res_temp$gene %like% "Alox5", ]

posForwarded <- logFC[logFC$Group == "Forwarded" & logFC$Input_logFC > 0, "Gene"]
writeLines(posForwarded, paste0("/Users/reikotachibana/Documents/ChungLab/riboseq/output/genes_list/", plot_title, "_posForwarded.txt"))
negForwarded <- logFC[logFC$Group == "Forwarded" & logFC$Input_logFC < 0, "Gene"]
writeLines(negForwarded, paste0("/Users/reikotachibana/Documents/ChungLab/riboseq/output/genes_list/", plot_title, "_negForwarded.txt"))

table(logFC$Group)

################################################################################
# GO Analysis
################################################################################
  
group <- "Forwarded"
gene_list <- ifelse(logFC$Group == group & logFC$Ribo_logFC < 0,
                    1, 0)
names(gene_list) <- logFC$Gene

GOdata <- new("topGOdata",
              ontology = "BP",  
              allGenes = gene_list,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher)
allRes

sum(gene_list)

# par(cex=0.5)
# showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 3, useInfo = "all")
# title(main = paste0(plot_title, " - ", group, " - logFC > 0"),
#       cex.main = 3,
#       line = -2)