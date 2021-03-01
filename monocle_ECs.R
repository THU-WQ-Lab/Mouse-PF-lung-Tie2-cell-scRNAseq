# preparations
seu.ECs = readRDS("ECs_seurat.rds")
seu.BLM = seu.ECs[, as.numeric(colnames(seu.ECs)) > 400]

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

cds = SeuratWrappers::as.cell_data_set(seu.ECs)
cds <- estimate_size_factors(cds)
#cds@colData@rownames = paste0("cell_", cds@colData@rownames)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(seu.ECs[["RNA"]])

cds <- preprocess_cds(cds, num_dim = 9)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

cds@int_colData@listData$reducedDims@listData$UMAP = seu.ECs@reductions$umap@cell.embeddings
cds <- cluster_cells(cds, resolution = 0.1)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "orig.ident", cell_size = 1, group_label_size = 4, labels_per_group = 0)

plot_cells(cds, genes = c("Lgals3"))

## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds, group_cells_by = "cluster")

EndMT_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- EndMT_pr_test_res %>% subset(q_value < 0.05) %>% .[order(.$q_value), ] %>% row.names()
plot_cells(cds, genes=pr_deg_ids[13:16],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

plot_cells(cds, genes=c("Plin2", "Ctss","Nfkbia","Hif1a",
                        "C3ar1", "Lgals3", "S100a4", "Fn1",
                        "Col4a1", "Icam1",
                        "Tgfbi","Apoe"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) +
  theme(text = element_text(size = 8))

plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1)


library(ggpubr)
vlnplots = list()
i = 1
for (gene in c("C3ar1", "Lgals3", "S100a4", "Fn1",
               #"Plin2", "Icam1",
               "Col4a1",
               #"Ctss",
               "Tgfbi","Apoe","Hif1a"
               #,"Nfkbia"
               )) {
  p = VlnPlot(seu.ECs, features = gene, pt.size = 0.3) + NoLegend() +
    theme(axis.title = element_blank(), text = element_text(size = 10), plot.margin = unit(c(0,0,0,0), units = "line"), 
          axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
  vlnplots[[i]] = p
  i = i + 1
}
p = ggarrange(plotlist = vlnplots, ncol = 4, nrow = 2)
annotate_figure(p, bottom = text_grob("Cluster", size = 14, vjust = -0.3), left = text_grob("Expression Level", size = 14, rot = 90))

saveRDS(cds, file = "mon3_ECs.rds")
cds = readRDS("mon3_ECs.rds")
