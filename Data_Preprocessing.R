library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(ggplot2)
library(harmony)
options(stringsAsFactors = FALSE)
set.seed(123)
options(future.globals.maxSize = 8000 * 1024^2)

my36colors <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3",
  "#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A",
  "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA",
  "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35",
  "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD",
  "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175"
)
Mycol_2 <- c("#FF7F0EFF", "#1F77B4FF")

sc_data <- read.table('./00_RawData/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt', header=T, row.names=1, skip=1)
metadata <- read.table('./00_RawData/GSE120575_patient_ID_single_cells.txt', row.name=1, header=F, sep="\t")
metadata <- metadata[,1:6]
colnames(sc_data) <- rownames(metadata)

PP <- as.data.frame(str_split_fixed(metadata$V5, "_", 2))
Patient_information <- as.data.frame(cbind(PP, metadata$V5, metadata$V6, metadata$V7))
rownames(Patient_information) <- rownames(metadata)
colnames(Patient_information) <- c("Pre_or_Post", "Patient", "Sample", "Re_or_NoR", "Medicine")

attach(Patient_information)
ifPD_1 <- Medicine == "anti-PD1"
Patient_information <- Patient_information[ifPD_1,]
sc_data <- sc_data[, ifPD_1]
metadata <- metadata[ifPD_1,]
detach(Patient_information)

ercc.index <- grep(pattern = "^ERCC", x = rownames(x = sc_data), value = FALSE)
sc_data <- sc_data[-ercc.index,]

seurat_obj <- CreateSeuratObject(counts = sc_data, meta.data = metadata, min.cells = 0, min.features = 200)

seurat_obj <- AddMetaData(object = seurat_obj, metadata = Patient_information)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rp"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL][[:digit:]]")

violin_before <- VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols = my36colors, pt.size = 0.01, ncol = 3)
ggsave("./03_Figures/01_QC_Violin_Before.pdf", plot = violin_before, width = 12, height = 6)

seurat_obj_filtered <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 2)

violin_after <- VlnPlot(seurat_obj_filtered, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), cols = my36colors, pt.size = 0.01, ncol = 3)
ggsave("./03_Figures/01_QC_Violin_After.pdf", plot = violin_after, width = 12, height = 6)

saveRDS(seurat_obj_filtered, "./02_Results/01_Seurat_Obj_Prescale.rds")

seurat_obj_filtered <- readRDS("./02_Results/01_Seurat_Obj_Prescale.rds")

seurat_obj_filtered <- NormalizeData(object = seurat_obj_filtered, normalization.method = "LogNormalize", scale.factor = 1e6)

seurat_obj_filtered <- FindVariableFeatures(object = seurat_obj_filtered, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(seurat_obj_filtered), 10)
plot1 <- VariableFeaturePlot(seurat_obj_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("./03_Figures/01_VariableFeatures.pdf", plot = plot2, width = 15, height = 10)

seurat_obj_filtered <- ScaleData(object = seurat_obj_filtered, features = VariableFeatures(seurat_obj_filtered))

seurat_obj_filtered <- RunPCA(object = seurat_obj_filtered, features = VariableFeatures(seurat_obj_filtered), verbose = FALSE)
seurat_obj_filtered <- RunHarmony(
  seurat_obj_filtered,
  group.by.vars = "Sample",
  plot_convergence = TRUE
)


saveRDS(seurat_obj_filtered, "./02_Results/01_Seurat_Obj_Preprocessed.rds")
cat("\n模块一完成！\n")
