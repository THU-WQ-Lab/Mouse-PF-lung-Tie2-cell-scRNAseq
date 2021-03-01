library(Seurat)

# only need one-time run
data.mmECEC = read.table('valid_counts_anno.txt', row.names=1, header=T, check.names = F)
gene.names = unique(colnames(data.mmECEC))
data.mmECEC.dedupe = data.frame(data.mmECEC[,FALSE])


for (gene in gene.names)
{
  occurence = colnames(data.mmECEC) == gene
  freq = sum(occurence)
  if (freq > 1)
  {
    data.mmECEC.dedupe[,gene] = rowSums(data.mmECEC[,occurence])
  }
  else
  {
    data.mmECEC.dedupe[,gene] = data.mmECEC[,gene]
  }
}

data.mmECEC.dedupe = data.mmECEC.dedupe[order(as.numeric(rownames(data.mmECEC.dedupe))), ]
data.mmECEC = t(data.mmECEC.dedupe)
write.table(data.mmECEC, file = 'valid_T.txt', quote = F, sep = '\t')

# seurat beginnings
data.mmECEC = read.table('valid_T.txt', row.names=1, header=T, check.names = F)
data.BLM = data.mmECEC[, as.numeric(colnames(data.mmECEC)) > 400]
data.NORM = data.mmECEC[, as.numeric(colnames(data.mmECEC)) <= 400]

seu.ECs = CreateSeuratObject(t(t(data.mmECEC)), min.cells = 10)

seu.BLM = CreateSeuratObject(t(t(data.BLM)), min.cells = 5)
seu.NORM = CreateSeuratObject(t(t(data.NORM)), min.cells = 5, is.expr = 1)

# seurat 3 read data
library(SingleCellExperiment)
data.mmECEC = read.table('valid_T.txt', row.names=1, header=T, check.names = F)

sce.ECs = SingleCellExperiment(list(counts = as.matrix(data.mmECEC)))

counts <- assay(sce.ECs, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(sce.ECs) <- log2(t(t(counts)/size.factors) + 1)
assayNames(sce.ECs)

seu.ECs = as.Seurat(sce.ECs)

# norm & scale
seu.ECs <- NormalizeData(seu.ECs)
seu.ECs <- FindVariableFeatures(seu.ECs, selection.method = "vst")

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(seu.ECs), 10)

# Plot variable features with and without labels
# plot1 = VariableFeaturePlot(seu.ECs)
# LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2

seu.ECs <- ScaleData(seu.ECs, features = rownames(seu.ECs))

seu.ECs <- RunPCA(seu.ECs, features = VariableFeatures(object = seu.ECs))

# Examine and visualize PCA results a few different ways
#print(seu.ECs[["pca"]], dims = 1:5, nfeatures = 5)

# VizDimLoadings(seu.ECs, dims = 1:2, reduction = "pca")
# 
# DimHeatmap(seu.ECs, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(seu.ECs, dims = 1:15, cells = 500, balanced = TRUE)

# Identify significant PCs
seu.ECs <- JackStraw(seu.ECs, num.replicate = 100)
seu.ECs <- ScoreJackStraw(seu.ECs, dims = 1:20)
JackStrawPlot(seu.ECs, dims = 1:20)
ElbowPlot(seu.ECs)

# Find group markers
seu.ECs$orig.ident = "Control"
seu.ECs$orig.ident[as.numeric(colnames(seu.ECs)) > 400] = "BLM"

Idents(seu.ECs) <- "orig.ident"
grp_markers = FindMarkers(seu.ECs, ident.1 = "Bleomycin", logfc.threshold = 0.58)
write.csv(grp_markers, file = "ECs_Blm_to_Norm.csv", quote = F)

# UMAP/TSNE
seu.ECs <- RunUMAP(seu.ECs, dims = 1:9)
seu.ECs <- RunTSNE(seu.ECs, dims = 1:9)

# Cluster the cells
seu.ECs <- FindNeighbors(seu.ECs, reduction = "pca", dims = 1:9)
seu.ECs <- FindClusters(seu.ECs, resolution = 0.5)

DimPlot(seu.ECs, reduction = "umap", pt.size = 1)
DimPlot(seu.ECs, reduction = "umap", group.by = "orig.ident")

DimPlot(seu.ECs, reduction = "tsne", pt.size = 1)
DimPlot(seu.ECs, reduction = "tsne", group.by = "orig.ident")

VlnPlot(seu.ECs, features = c("Lgals3")) + NoLegend()

# all markers
all.mkrs = FindAllMarkers(seu.ECs, logfc.threshold = 0.58)
all.mkrs = all.mkrs %>% subset(subset = avg_logFC > 0 & p_val_adj < 0.05)
write.table(all.mkrs, file = "all_cluster_markers.csv", sep = ",", quote = F, row.names = F)

# S/L
saveRDS(seu.ECs, file = "ECs_seurat.rds")
seu.ECs = readRDS("ECs_seurat.rds")

# reactome pathways
library(ReactomePA)
library(ggplot2)

cluster_1_mkrs = FindMarkers(seu.ECs, ident.1 = 1, logfc.threshold = 0.58)
cluster_1_mkrs$gene = rownames(cluster_1_mkrs)
# cluster_1_mkrs$avgFC = 2** cluster_1_mkrs$avg_logFC
gene_Entrez_ID = clusterProfiler::bitr(cluster_1_mkrs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
colnames(gene_Entrez_ID)[1] = "gene"
cluster_1_mkrs = merge(cluster_1_mkrs, gene_Entrez_ID, by.x = "gene")
geneList <- cluster_1_mkrs[,"avg_logFC"]
names(geneList) = as.character(cluster_1_mkrs$ENTREZID)

pathway = enrichPathway(gene = cluster_1_mkrs$ENTREZID, pvalueCutoff = 0.05, readable = T, organism = "mouse")

pw_for_plot = pathway@result %>% subset(p.adjust < 0.05)
n = nrow(pw_for_plot)
n_genes = pw_for_plot$geneID %>% strsplit(split = "/") %>% unlist() %>% unique() %>% length()
n_max_char_in_pw = pw_for_plot$Description %>% nchar() %>% max()

png(filename = "./ECs_reactome_barplot.png", width = floor(6.4*n_max_char_in_pw) + 1000, height = 18*n + 150, res = 120)
print(barplot(pathway, showCategory = n) + ylab("Number of genes matched"))
dev.off()

png(filename = "./ECs_reactome_cnetplot.png", width = floor(6.4*n_max_char_in_pw) + 2000, height = 50*n + 800, res = 100)
try(print(cnetplot(pathway, categorySize = "p.adjust", foldChange = geneList, showCategory = n) +
            scale_color_gradient2(low = "blue", high = "red", mid = "lightgrey", midpoint = 0, guide = guide_colorbar(title="log2FC"))))
dev.off()

png(filename = "./ECs_reactome_heatplot.png", width = floor(5*n_max_char_in_pw) + 20*n_genes + 120, height = 20*n + 220, res = 120)
try(print(heatplot(pathway, foldChange = geneList, showCategory = n) +
            scale_fill_gradient2(low = "blue", high = "red", mid = "lightgrey", midpoint = 0, guide = guide_colorbar(title="log2FC"))))
dev.off()

# int_genes = names(table(msc_to_EC_interaction$gene_to))
int_pw = pathway@result
# int_pw$has_interaction_genes = sapply(int_pw$geneID, function(x){return(length(intersect(int_genes, strsplit(x, "/")[[1]])) > 0)})
# int_pw = subset(int_pw, subset = has_interaction_genes == TRUE)
# int_pw$interaction_geneID = sapply(int_pw$geneID, function(x){return(paste(intersect(int_genes, strsplit(x, "/")[[1]]), collapse = "; "))})

if (nrow(int_pw) > 0) {
  write.table(int_pw %>% subset(p.adjust < 0.05), file = "./ECs_reactome_all_pathways.tsv", sep = '\t', quote = F, row.names = F)
}

DotPlot(seu.ECs, features = c("", ""))
