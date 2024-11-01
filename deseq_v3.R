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

ribo_file <- "data/ribo_counts.txt"
rna_file <- "data/rna_counts.txt"

ribo <- read.delim(ribo_file) 
rna <- read.delim(rna_file)

ribo$transcript <- sub("-.*", "", ribo$transcript)
genes <- ribo[, "transcript"]
merge <- cbind(ribo, rna)

# HSC vehicle vs HSC Ven
# comparison <- c("RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B")
# comparison <- c("RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_HSC_Ven_A", "RNA_HSC_Ven_B")
# comparison <- c("RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B",
#                 "RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_HSC_Ven_A", "RNA_HSC_Ven_B")

# GMP vehicle vs GMP Ven
# comparison <- c("RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B")
# comparison <- c("RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")
# comparison <- c("RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B",
#                 "RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")

# HSC Ven vs GMP Ven
# comparison <- c("RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B")
# comparison <- c("RNA_HSC_Ven_A", "RNA_HSC_Ven_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")
comparison <- c("RIBO_HSC_Ven_A", "RIBO_HSC_Ven_B", "RIBO_GMP_Ven_A", "RIBO_GMP_Ven_B",
                "RNA_HSC_Ven_A", "RNA_HSC_Ven_B", "RNA_GMP_Ven_A", "RNA_GMP_Ven_B")

# HSC vehicle vs GMP vehicle
# comparison <- c("RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B")
# comparison <- c("RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B")
# comparison <- c("RNA_HSC_vehicle_A", "RNA_HSC_vehicle_B", "RNA_GMP_vehicle_A", "RNA_GMP_vehicle_B",
#                 "RIBO_HSC_vehicle_A", "RIBO_HSC_vehicle_B", "RIBO_GMP_vehicle_A", "RIBO_GMP_vehicle_B")
  
subset <- merge[, comparison]

coldata <- data.frame(colnames(subset), colnames(subset))
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

dds <- DESeqDataSetFromMatrix(countData = subset,
                              colData = coldata,
                              # design = ~ Condition + LibraryType + Condition:LibraryType)
                              design = ~ Population + LibraryType + Population:LibraryType)
                              # design = ~ Condition)
                              # design = ~ Population)
# keep <- rowSums(counts(dds) >= 10) >= 2
# dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef=4,res=res,type="apeglm")
summary(res)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Population", "LibraryType", "Condition"))

res$gene <- merge$transcript
EnhancedVolcano(res,
                lab = res$gene,
                x = 'log2FoldChange',
                y = 'pvalue')


# plotMA(res)


# Ribo
ribo_comparison <- comparison[comparison %like% "RIBO"]
ribo <- ribo[, ribo_comparison]
coldata_ribo <- get_coldata(ribo)
ddsMat_ribo <- DESeqDataSetFromMatrix(countData = ribo,
                                      colData = coldata_ribo, 
                                      # design =~ Condition)
                                      design =~ Population)
ddsMat_ribo <- DESeq(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo)
res_ribo <- lfcShrink(ddsMat_ribo, coef=2,res=res_ribo,type="apeglm") 

# RNA 
rna_comparison <- comparison[comparison %like% "RNA"]
rna <- rna[, rna_comparison]
coldata_rna <- get_coldata(rna)
ddsMat_rna <- DESeqDataSetFromMatrix(countData = rna,
                                     colData = coldata_rna, 
                                     # design =~ Condition)
                                     design =~ Population)
ddsMat_rna <- DESeq(ddsMat_rna)
res_rna <- results(ddsMat_rna)
res_rna <- lfcShrink(ddsMat_rna, coef=2,res=res_rna,type="apeglm")

logFC <- data.frame(Ribo = res_ribo$log2FoldChange,
                    Input = res_rna$log2FoldChange)
logFC <- na.omit(logFC)

ggplot(logFC, aes(x=Input, y=Ribo)) + 
  geom_point()







# max_val = max(res_ribo[,2],res_rna[,2],na.rm = T)
# plot(y=res_ribo[,2],x=res_rna[,2], xlab="RNA-seq log2 fold change", 
#      ylab = "Ribo-seq log2 fold change",asp=1,pch=16,col=rgb(128/255,128/255,128/255,0.1),ylim=c(-max_val,max_val),xlim=c(-max_val,max_val),cex=0.4)
# abline(a=0,b=1,col="gray")
# abline(h=0,v=0,col="gray")
# points(y=res_ribo[forwarded,2],x=res_rna[forwarded,2],pch=16,col=rgb(0,0,1,1))
# points(y=res_ribo[exclusive,2],x=res_rna[exclusive,2],pch=16,col=rgb(1,0,0,1))
# points(y=res_ribo[intensified,2],x=res_rna[intensified,2],pch=16,col=rgb(1,0,1,1))
# points(y=res_ribo[buffered,2],x=res_rna[buffered,2],pch=16,col=rgb(1,0,1,1))
