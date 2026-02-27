library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
options(stringsAsFactors = F)
options(future.seed = TRUE)
set.seed(123)

# Color scheme
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
Mycol_2 <-c('#FF7F0EFF', '#1F77B4FF')

# 10 additional colors
new_10_colors <- c('#7D6B9C', '#849E5D', '#C45B4D', '#546E89', '#D9A645', 
                   '#689F88', '#B56C8C', '#4C5D4F', '#9B815F', '#5E8CBD')

# Combine to generate my46colors
my46colors <- c(my36colors, new_10_colors)

# 1. Load preprocessed data
cat("Step 1: Loading preprocessed data\n")
seurat_obj <- readRDS("./02_Results/01_Seurat_Obj_Preprocessed.rds")
# seurat_obj <- readRDS("./02_Results/01_Seurat_Obj_Preprocessed_SCT.rds")

# 2. Clustering analysis
cat("\nStep 2: Clustering analysis\n")

seurat_obj <- RunTSNE(
  seurat_obj, 
  reduction = "harmony", 
  dims = 1:15,
  perplexity = 40,  # Add parameters to ensure stability
  seed.use = 123    # Explicitly set seed for reproducible TSNE results
)

# Run FindNeighbors once (adjacency matrix only needs to be calculated once)
seurat_obj <- FindNeighbors(
  seurat_obj, 
  reduction = "harmony",
  dims = 1:15,
  verbose = FALSE
)

# Try different resolutions
resolutions <- seq(1, 4.6, by=0.3)
cat(paste("Exploring resolution range:", paste(resolutions, collapse = ", "), "\n"))

for (res in resolutions) {
  cat(paste("Clustering with resolution", res, "...\n"))
  seurat_obj <- FindClusters(
    seurat_obj, 
    resolution = res, 
    algorithm = 1,
    cluster.name = paste0("RNA_snn_res.", res),  # Explicitly name the cluster column to avoid overwriting
    verbose = FALSE
  )
  
  # Draw tsne plot (colored by current resolution)
  tsne_plot <- DimPlot(
    seurat_obj, 
    reduction = "tsne", 
    label = TRUE, 
    cols = my46colors,
    group.by = paste0("RNA_snn_res.", res)  # Specify cluster column for current resolution
  ) + ggtitle(paste("Clustering Resolution =", res))
  
  # Fix filename (avoid spaces/special characters)
  ggsave(
    paste0("./03_Figures/02_tsne_Cluster_Res_", res, ".pdf"), 
    plot = tsne_plot, 
    width = 12, 
    height = 10
  )
}

# 4. Draw clustree plot to show clustering relationships at different resolutions
cat("Drawing clustree plot to show resolution selection...\n")
library(clustree)
clustree_plot <- clustree(
  seurat_obj, 
  prefix = "RNA_snn_res.",
  theme = theme_minimal()  # Optimize plot style
) + theme(legend.position = "bottom")
ggsave("./03_Figures/02_Clustree_Resolution.pdf", plot = clustree_plot, width = 10, height = 12)

# 5. Select optimal resolution (1.4), only run FindClusters (no need to repeat FindNeighbors/RunTSNE)
cat("\nStep 3: Select optimal resolution and determine final clusters\n")
seurat_obj <- FindClusters(
  seurat_obj, 
  resolution = 4.3, 
  algorithm = 1,
  cluster.name = "seurat_clusters",  # Name as default column for convenient downstream analysis
  verbose = FALSE
)

# 3. Cell type annotation
cat("\nStep 3: Cell type annotation\n")
# seurat_obj <- JoinLayers(seurat_obj)

# 3.1 SingleR automatic annotation
cat("Step 3.1: SingleR automatic annotation...\n")
library(SingleR)
library(celldex)
library(BiocParallel)
register(SnowParam(workers = 4))  # Adjust according to computer core count

# Download and build human PBMC reference (skip if already downloaded)
ref <- celldex::HumanPrimaryCellAtlasData()  # Or MonacoImmuneData / BlueprintEncodeData

# Get normalized expression matrix
norm_counts <- GetAssayData(seurat_obj, layer = "data")

# Run SingleR (annotate by cluster)
cat("Running SingleR cluster annotation...\n")

# Extract cluster information
clusters <- seurat_obj$seurat_clusters
unique_clusters <- unique(clusters)
cat("Unique clusters:\n")
print(unique_clusters)

# Run SingleR (annotate by cluster)
singler_res <- SingleR(test = norm_counts,
                       ref = ref,
                       labels = ref$label.fine,  # Use fine labels
                       clusters = clusters)

# Check SingleR result structure
print("SingleR result rownames (clusters):\n")
print(rownames(singler_res))
print("SingleR assigned labels:\n")
print(singler_res$labels)

# Map SingleR results to each cell
# Ensure cluster ID format is consistent
cat("Ensuring cluster ID format is consistent...\n")
cluster_ids <- as.character(unique_clusters)
singler_cluster_ids <- rownames(singler_res)

cat("Cluster IDs in data:\n")
print(cluster_ids)
cat("Cluster IDs in SingleR results:\n")
print(singler_cluster_ids)

# Create mapping from cluster to label
cluster_to_label <- singler_res$labels
names(cluster_to_label) <- singler_cluster_ids

# Verify mapping
cat("Cluster to label mapping:\n")
print(cluster_to_label)

# Assign corresponding label to each cell
cat("Assigning labels to each cell...\n")
cell_cluster_ids <- as.character(clusters)

# Get label values (without names)
labels_values <- unname(cluster_to_label[cell_cluster_ids])
cat("Label values length:", length(labels_values), "\n")
cat("Top 10 label values:\n")
print(head(labels_values))

# Direct assignment (independent of name matching)
seurat_obj$singler_label <- labels_values

# Verify label assignment
cat("Label assignment for top 10 cells:\n")
print(head(seurat_obj$singler_label, 10))

# Assign corresponding score to each cell (take highest score)
cat("Assigning scores to each cell...\n")
cluster_to_score <- apply(singler_res$scores, 1, max)
names(cluster_to_score) <- singler_cluster_ids

# Get score values (without names)
scores_values <- unname(cluster_to_score[cell_cluster_ids])
cat("Score values length:", length(scores_values), "\n")
cat("Top 10 score values:\n")
print(head(scores_values))

# Direct assignment (independent of name matching)
seurat_obj$singler_score <- scores_values

# Verify score assignment
cat("Score assignment for top 10 cells:\n")
print(head(seurat_obj$singler_score, 10))

# Visualize SingleR annotation
cat("Visualizing SingleR annotation...\n")
p_singler <- DimPlot(seurat_obj, group.by = "singler_label", reduction = "tsne", label = TRUE, cols = my36colors) +
  ggtitle("SingleR Annotation")
ggsave("./03_Figures/02_tsne_SingleR.pdf", plot = p_singler, width = 12, height = 10)

# Print SingleR annotation results
cat("\nSingleR annotation results:\n")
print(singler_res$labels)

# 3.2 Manual Marker annotation (compared with SingleR)
cat("\nStep 3.2: Manual Marker annotation (compared with SingleR)...\n")
# Calculate significant differentially expressed genes for each cluster
cat("Calculating significant differentially expressed genes for each cluster...\n")
# Find marker genes for each cluster
cat("Finding marker genes for each cluster...\n")
cluster_markers <- FindAllMarkers(seurat_obj, 
                                  group.by = "seurat_clusters", 
                                  logfc.threshold = 0.25, 
                                  min.pct = 0.25, 
                                  test.use = "wilcox")
write.csv(cluster_markers, "./04_Tables/02_Cluster_Markers.csv", row.names = FALSE)

# # Calculate DEGs of cluster 4 and 9 relative to all other clusters
# # First combine cluster 4 and 9 into a new group
# Idents(seurat_obj) <- "seurat_clusters"
# seurat_obj$group <- ifelse(Idents(seurat_obj) %in% c("7","18"), "cluster_7_18", "others")
# 
# # Calculate DEGs
# markers_7_18_vs_others <- FindMarkers(
#   seurat_obj,
#   group.by = "group", 
#   ident.1 = "cluster_7_18",
#   ident.2 = "others",
#   min.pct = 0.25,   # Expression proportion threshold, adjustable
#   logfc.threshold = 0.25 # logFC threshold, adjustable
# )
# write.csv(markers_7_18_vs_others, "./04_Tables/02_Cluster_Markers_markers_7_18_vs_others.csv", row.names = T)
# 
# 
# Extract top 10 DEGs for each cluster
top20_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
write.csv(top20_markers, "./04_Tables/02_Cluster_top20_markers.csv", row.names = FALSE)

# This plot shows: does a certain cell population only appear in Responder?
tsne_split_clusters <- DimPlot(
  seurat_obj, 
  reduction = "tsne", 
  group.by = "seurat_clusters", # Key: color by [cell type/cluster]
  split.by = "Re_or_NoR",       # Key: split by [group]
  label = TRUE,
  repel = TRUE
) 
ggsave("./03_Figures/02_tsne_Response_ByCluster.pdf", plot = tsne_split_clusters, width = 16, height = 8)

# Define typical Marker genes
markers <- list(
  "CD8+ T" = c("CD3E", "CD8A"),
  "CD4+ T" = c("CD3E", "CD4"),
  "Treg" = c("CD3E", "FOXP3"),
  "B cells" = c("CD19", "MS4A1"),
  "Plasma cells" = c("MZB1", "SDC1"),
  "NK cells" = c("NCR1", "FCGR3A"),
  "Monocytes" = c("CD14", "LYZ"),
  "DCs" = c("CD1C", "CLEC9A"),
  "pDCs" = c("CLEC4C", "IL3RA"),
  "gd T cells" = c("TRDV1", "TRDC"),
  "Macrophages" = c("MARCO", "MERTK")
)
markers <- list(
  "CD8+ T"              = c("CD3E", "CD8A"),
  "CD4+ T"              = c("CD3E", "CD4"),
  "Tregs"               = c("FOXP3", "IL2RA"), 
  "NK cells"            = c("GNLY", "FCGR3A"), 
  # "Clonal T cells"      = c("TRBV6-6", "TRBV6-5"),
  # "Proliferating T cells" = c("TYMS", "DUT"), #"BIRC5", "MKI67"
  "B cells"             = c("CD79A", "MS4A1"),
  "Monocytes"           = c("CD14", "S100A8"),
  "Gamma-Delta T"       = c("TRDV1", "TRDC"),
  "Plasma cells"        = c("MZB1", "SDC1"), 
  "Macrophages"         = c("C1QA", "TREM2"), 
  "pDCs"                = c("CLEC4C", "IL3RA")#,
  # "cDCs"                = c("CLEC10A", "FCER1A")
)
# Draw Marker DotPlot
all_markers <- unique(unlist(markers))
dotplot <- DotPlot(seurat_obj, group.by = "seurat_clusters", features = all_markers) +
  coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./03_Figures/02_Marker_DotPlot.pdf", plot = dotplot, width = 20, height = 16)

# # Locate cluster with high TRDV1 expression (key gene for Gamma-Delta T cells)
# trdv1_expr <- AverageExpression(seurat_obj, features = "TRDV1", assays = "RNA", group.by = "seurat_clusters")
# trdv1_cluster <- which.max(trdv1_expr$RNA["TRDV1",])
# cat(paste("Cluster with high TRDV1 expression:", trdv1_cluster - 1, "\n"))

# Draw FeaturePlot for each Marker (for manual annotation reference)
library(patchwork)
marker_plots <- lapply(all_markers, function(gene) {
  FeaturePlot(seurat_obj, features = gene, reduction = "tsne", cols = c("lightgrey", "blue")) +
    ggtitle(gene)
})
marker_combined <- wrap_plots(marker_plots, ncol = 5)
ggsave("./03_Figures/02_Marker_Expression_tsnes.pdf", plot = marker_combined, width = 30, height = 20)

# Draw Violin plot for each Marker
library(scCustomize)
marker_VlnPlot <- Stacked_VlnPlot(seurat_object = seurat_obj, 
                                  features = all_markers, 
                                  group.by = "seurat_clusters",
                                  x_lab_rotate = TRUE,
                                  colors_use = my46colors)
ggsave("./03_Figures/02_Marker_Expression_VlnPlot.pdf", plot = marker_VlnPlot, width = 8, height = 10)


# 3.3 Comprehensive annotation: combine SingleR and Marker results
cat("\nStep 3.3: Comprehensive annotation (SingleR + Marker)...\n")
# Provide tsne with cluster numbers for easy manual checking
tsne_ref <- DimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "tsne", label = TRUE, cols = my46colors) +
  ggtitle("tsne with Cluster IDs (for manual annotation)")
ggsave("./03_Figures/02_tsne_ClusterIDs_Reference.pdf", plot = tsne_ref, width = 12, height = 10)

tsne_ref <- DimPlot(seurat_obj, group.by = "seurat_clusters", reduction = "tsne", label = TRUE, cols = my46colors) +
  ggtitle("tSNE with Cluster IDs (for manual annotation)")
ggsave("./03_Figures/02_tSNE_ClusterIDs_Reference.pdf", plot = tsne_ref, width = 12, height = 10)

# Example: Fill in after checking against SingleR and Marker plots
# Users need to modify according to actual expression
cell_type_annotations <- list(
  "0" = "CD8 T cells", 
  "1" = "CD4 T cells",  
  "2" = "CD8 T cells",
  "3" = "Regulatory T cells",  
  "4" = "CD8 T cells", 
  "5" = "NK cells", 
  "6" = "CD8 T cells", 
  "7" = "Regulatory T cells", 
  "8" = "Mon/Mac",   
  "9" = "B cells",  
  "10" = "CD8 T cells", 
  "11" = "NK cells",  
  "12" = "CD8 T cells",  
  "13" = "CD8 T cells",  
  "14" = "CD4 T cells",  
  "15" = "CD8 T cells",  
  "16" = "Mon/Mac",  
  "17" = "CD8 T cells",  
  "18" = "CD8 T cells",  
  "19" = "Plasma cells", 
  "20" = "CD4 T cells",  
  "21" = "gd T cells",  
  "22" = "CD8 T cells",  
  "23" = "pDCs",   
  "24" = "CD4 T cells", 
  "25" = "CD8 T cells",   
  "26" = "CD4 T cells", 
  "27" = "CD8 T cells",   
  "28" = "CD8 T cells",  
  "29" = "CD4 T cells", 
  "30" = "CD8 T cells", 
  "31" = "CD4 T cells", 
  "32" = "gd T cells",  
  "33" = "CD4 T cells",  
  "34" = "Mon/Mac",  
  "35" = "CD8 T cells",  
  "36" = "CD8 T cells",  
  "37" = "Plasma cells",  
  "38" = "Regulatory T cells", 
  "39" = "Plasma cells",  
  "40" = "Plasma cells", 
  "41" = "CD8 T cells",  
  "42" = "CD8 T cells",  
  "43" = "CD8 T cells",  
  "44" = "CD8 T cells" 
)
# Apply manual comprehensive annotation
Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, cell_type_annotations)
seurat_obj$manual_cellTypes <- Idents(seurat_obj)

# 3.4 Compare consistency between SingleR and manual annotation
cat("\nStep 3.4: Comparing consistency between SingleR and manual annotation...\n")
comp_df <- data.frame(
  Cluster = seurat_obj$seurat_clusters,
  SingleR = seurat_obj$singler_label,
  Manual  = seurat_obj$manual_cellTypes
)
write.csv(comp_df, "./04_Tables/02_Compare_SingleR_Manual.csv", row.names = FALSE)

# Draw consistency heatmap
library(pheatmap)
conf_mat <- table(Manual = comp_df$Manual, SingleR = comp_df$SingleR)
pheatmap(conf_mat, filename = "./03_Figures/02_SingleR_Manual_heatmap.pdf", width = 8, height = 6)

# Output final annotation results
cat("\nFinal cell type annotation results:\n")
print(table(seurat_obj$manual_cellTypes))

# Draw final cell type tsne
final_tsne <- DimPlot(seurat_obj, group.by = "manual_cellTypes", reduction = "tsne", label = TRUE, cols = my46colors)
ggsave("./03_Figures/02_tsne_CellTypes.pdf", plot = final_tsne, width = 12, height = 10, device = cairo_pdf)

final_tsne <- DimPlot(seurat_obj, group.by = "manual_cellTypes", reduction = "tsne", label = TRUE, cols = my46colors)
ggsave("./03_Figures/02_tSNE_CellTypes.pdf", plot = final_tsne, width = 12, height = 10, device = cairo_pdf)

# Draw Violin plot for each Marker
cellTypes_marker_VlnPlot <- Stacked_VlnPlot(seurat_object = seurat_obj, 
                                  features = all_markers, 
                                  group.by = "manual_cellTypes",
                                  x_lab_rotate = TRUE,
                                  colors_use = my46colors)
ggsave("./03_Figures/02_cellTypes_Marker_Expression_VlnPlot.pdf", plot = cellTypes_marker_VlnPlot, width = 8, height = 10, device = cairo_pdf)

# Checkpoint: ensure all samples (R and NR) are displayed on the same tsne
cat("Drawing R/NR distribution tsne plot...\n")
tsne_response <- DimPlot(seurat_obj, reduction = "tsne", group.by = "Re_or_NoR", cols = Mycol_2, label = TRUE)
ggsave("./03_Figures/02_tsne_Response.pdf", plot = tsne_response, width = 12, height = 10, device = cairo_pdf)

# 4. Proportion analysis (Per-Sample)
cat("\nStep 4: Proportion analysis (Per-Sample)\n")

# 1. Extract cell types and patient IDs
Idents(seurat_obj) <- "manual_cellTypes"
cell_types <- Idents(seurat_obj)
Samples <- seurat_obj$Sample

# 2. Construct cell type x patient count matrix
# This step calculates the number of each cell type per patient
counts_table <- table(cell_types, Samples)
# Keep only T cells
# counts_table <- counts_table[grep("T", rownames(counts_table)), ]

# 3. Calculate proportion (normalize by column/patient)
# margin = 2 means calculating percentage by column, making the sum of proportions for each patient 100
props_table <- prop.table(counts_table, margin = 2) * 100

# 4. Convert to long format dataframe (suitable for ggplot)
prop_data <- as.data.frame(props_table)
colnames(prop_data) <- c("CellType", "Sample", "Proportion")

# 5. Add grouping information (Response: R vs NR)
# Create a patient ID -> group mapping table
Sample_meta <- unique(seurat_obj@meta.data[, c("Sample", "Re_or_NoR")])
# Merge grouping information into proportion data
prop_data <- merge(prop_data, Sample_meta, by = "Sample")
# Rename column to match subsequent plotting code
colnames(prop_data)[colnames(prop_data) == "Re_or_NoR"] <- "Response"

# --- Core modification section End ---

# Check data (optional)
print(head(prop_data))

# 6. Statistical test (calculate P value and print)
# Now our sample size n is the number of patients, statistical results are more reliable
cat("Using Wilcoxon test to compare proportion differences of cell types between different patient groups...\n")
all_cell_types <- unique(prop_data$CellType)

for (ct in all_cell_types) {
  # Extract data for this cell type
  sub_data <- prop_data[prop_data$CellType == ct, ]
  
  # Group data
  r_vals <- sub_data$Proportion[sub_data$Response == "Responder"]
  nr_vals <- sub_data$Proportion[sub_data$Response == "Non-responder"] # Please confirm this is also the spelling in your metadata
  
  # Only test when both groups have data
  if (length(r_vals) > 1 && length(nr_vals) > 1) {
    test_res <- wilcox.test(r_vals, nr_vals)
    cat(paste0(ct, ": p-value = ", signif(test_res$p.value, 4), "\n"))
  } else {
    cat(paste0(ct, ": Insufficient samples, skipping test\n"))
  }
}

# 7. Draw boxplot
cat("Drawing boxplot for proportion differences...\n")

# Note: It is recommended to set the x-axis to Response for more intuitive comparison;
# If you want to keep x as CellType, facet is also fine, but I fine-tuned the x parameter to meet publication standards
p_box <- ggboxplot(prop_data, 
                   x = "Response",      # Change to group as X axis
                   y = "Proportion", 
                   color = "Response",
                   palette = Mycol_2,
                   add = "jitter",      # Show points for each patient
                   facet.by = "CellType", 
                   ncol = 4,
                   scales = "free_y",   # Key: make Y-axis scales independent for each facet, as proportions vary greatly among cell types
                   ylab = "Proportion (%)", 
                   xlab = "")

# Add statistical significance markers
p_box <- p_box + stat_compare_means(method = "wilcox.test", # Use Wilcoxon test 
                                    #method = "t.test",
                                    # label = "p.signif", # Show asterisks
                                    label = "p.format", # Show numerical values
                                    method.args = list(exact = F),
                                    comparisons = list(c("Responder", "Non-responder"))) # Specify comparison groups
# Increase space for P value display
p_box <- p_box + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
# Adjust theme
p_box <- p_box + theme(legend.position = "none", # Legend can be removed since x-axis is already grouped
                       axis.text.x = element_text(angle = 45, hjust = 1)) # Tilt x-axis labels to prevent overlap

# Save
ggsave("./03_Figures/02_CellType_Proportion_Boxplot_PerSample.pdf", plot = p_box, width = 8, height = 10, device = cairo_pdf)
print("Plotting completed!")

# 5. Save results
cat("\nStep 5: Saving results\n")
saveRDS(seurat_obj, "./02_Results/02_Seurat_Obj_Annotated.rds")

# Save cell type proportion data
write.csv(prop_data, "./04_Tables/02_CellType_Proportions.csv", row.names = FALSE)

cat("\nModule 2 completed!\n")
