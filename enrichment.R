library(ggplot2)
library(dplyr)
library(stringr)

enrichment_file <- paste0("/Users/reikotachibana/Documents/ChungLab/riboseq/", 
                          "DESeq_TE_HSC_Ven_vs_GMP_Ven_KEGG-RA-WP-BP.csv")

enrichment_data <- read.csv(enrichment_file)

up_pathways <- enrichment_data %>%
  filter(ES > 0 & str_detect(pathway, "RA_")) %>%
  arrange(padj, desc(ES), pval) %>%
  head(10)

down_pathways <- enrichment_data %>%
  filter(ES < 0 & str_detect(pathway, "RA_")) %>%
  arrange(desc(padj), desc(ES), desc(pval)) %>%
  tail(10)


pathways <- rbind(up_pathways, down_pathways)
pathways$Regulation <- ifelse(pathways$ES > 0, "Upregulated", "Downregulated")
pathways$pathway <- gsub("%R-MMU-\\d+", "", pathways$pathway)

ggplot(pathways, aes(x = reorder(pathway, ES), y = ES, fill = Regulation)) +
  geom_bar(stat = "identity") + 
  coord_flip() + 
  xlab ("RA Pathway") + 
  ylab("Enrichment Score (ES)") + 
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3")) + 
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))  
