#Creation of the Seurat object PTEN_9mo after gene invalidation treated with vehicle
data2_10x = Read10X(data.dir = "filtered_feature_bc_matrix172/")
Ech1 = CreateSeuratObject(data2_10x,
                          project = "9mo_vehicle",
                          assay = "RNA",
                          min.cells = 10,
                          min.features = 100,
                          names.field = 1,
                          names.delim = "_",
                          meta.data = NULL
)
Ech1$orig.ident="Ctl"
Ech1[["percent.mt"]] <- PercentageFeatureSet(Ech1, pattern = "^mt-")
VlnPlot(Ech1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)
ggsave("VlnPlot_vehicle.tiff")

#Creation of the Seurat object PTEN_9mo after gene invalidation treated with Gem72 for 1 week
data2_10x = Read10X(data.dir = "filtered_feature_bc_matrix173/")
Ech2 = CreateSeuratObject(data2_10x,
  project = "9mo_gem72",
  assay = "RNA",
  min.cells = 10,
  min.features = 100,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
Ech2$orig.ident="Gem72"
Ech2[["percent.mt"]] <- PercentageFeatureSet(Ech2, pattern = "^mt-")
VlnPlot(Ech2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)
ggsave("VlnPlot_Gem72.tiff")
remove(data2_10x)

#Creation of the combined Seurat object 
Ech1 <- subset(Ech1, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 20)
Ech2 <- subset(Ech2, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 20)
Ech.list <- list(Ech1, Ech2)
Gem <- lapply(X = Ech.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

Gem <- FindIntegrationAnchors(object.list = Gem, dims = 1:30)
Gem <- IntegrateData(anchorset = Gem, dims = 1:30)

#Visualisation du tSNE
DefaultAssay(Gem) <- "integrated"

# Run the standard workflow for visualization and clustering
Gem <- ScaleData(Gem, verbose = FALSE)
Gem <- RunPCA(Gem, npcs = 10, verbose = FALSE)

# t-SNE and Clustering
Gem <- FindNeighbors(Gem, reduction = "pca", dims = 1:10)
Gem <- FindClusters(Gem, resolution = 0.8)
Gem <- RunTSNE(Gem, dims.use = 1:10, reduction.use = "pca", perplexity = 30)
```
