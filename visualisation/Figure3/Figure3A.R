# Load required libraries
library(NMF)
library(bigmemory)

# Read gene expression data (genes in rows, samples in columns)
data <- read.table("aml_up_tpm.tsv", header=T, sep='\t', row.names=1)

# Feature selection: Select top 5000 genes with highest variance
row_variances <- apply(data, 1, var)
top_indices <- order(row_variances, decreasing = TRUE)[1:5000]
data <- data[top_indices, ]

# Data preprocessing: Apply log2 transformation (add pseudocount to avoid zeros)
data <- log2(data + 1.0001)

# Run NMF with rank range 2-15 for model selection
rank <- 2:15
nmf_res <- nmf(data, rank=rank, nrun=10000, method="brunet", seed=10073, .options='P40')

# Save workspace and plot cophenetic correlation
save.image(file="mads_rank2_15_nrun10000_normal.Rdata")
pdf(file="cophenetic_2_15_run10000_normal.pdf", width=8, height=7, onefile=F)
plot(nmf_res)
dev.off()

# Final model with selected rank (8)
rank <- 8
nmf_res <- nmf(data, rank=rank, nrun=10000, method="brunet", seed=123456, .options='P40')

# Save clustering results
Cluster <- predict(nmf_res)
Cluster <- as.data.frame(Cluster)
Cluster$Cluster <- paste0("C", Cluster$Cluster)
clusterOut <- rbind(ID=colnames(Cluster), Cluster)
write.table(clusterOut, file="cluster_rank8-normal-5k.txt", sep="\t", quote=F, col.names=F)

# Generate consensus matrix heatmap
pdf(file="heatmap_rank8.pdf", width=6, height=6, onefile=F)
consensusmap(nmf_res,
             annRow=NA,
             annCol=NA,
             main="Consensus matrix",
             info=FALSE)
dev.off()