library("DESeq2")
library("readr")
library("ggplot2")

# --- Loading the sample sheet ---

samples = read.table("sample_sheet.csv", header = TRUE, sep = ",")

# --- Loading the processed final count matrix ---

count_data = read.table("deseq2ready_gene_count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")

# --- Creating DESeq2 dataset ---

deseq2_dataset = DESeqDataSetFromMatrix(countData = count_data, colData = samples, design = ~ condition)

# --- Running the differential expression analysis and ordering the result by adjusted p-value---

deseq2_dataset = DESeq(deseq2_dataset)
degs = results(deseq2_dataset)
qordered_degs = degs[order(degs$padj),]

# --- Saving the results to a CSV file ---
write.csv(as.data.frame(qordered_degs), file = "differential_expression_results.csv")

# --- Visualizing the results ---
# --- Running PCA ---

vsd = varianceStabilizingTransformation(deseq2_dataset, blind = FALSE)
pcaData = plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

# --- Creating and saving PCA plot ---

PCA_plot = ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title = "PCA Plot") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("PCA_plot.png", plot = PCA_plot)

# --- Creating and saving MA plot ---

png("MA_plot.png")
plotMA(degs, main = "MA Plot", ylim = c(-5,5))
dev.off()

# --- Creating and saving Volcano plot ---

degs$log2FoldChange[is.na(degs$log2FoldChange)] = 0
degs$padj[is.na(degs$padj)] = 1
volcano_data = data.frame(log2FoldChange = degs$log2FoldChange, 
                           negLogPval = -log10(degs$padj),
                           Gene = rownames(degs))

ggplot(volcano_data, aes(x = log2FoldChange, y = negLogPval)) +
  geom_point(aes(color = negLogPval > -log10(0.05) & abs(log2FoldChange) > 1), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave("Volcano_plot.png")
