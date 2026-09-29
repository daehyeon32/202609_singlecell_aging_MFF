library(Seurat)
library(Cairo)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)
library(glue)
library(ggsignif)
library(ggplot2)
library(ggrepel)
library(doMC)
library(fastSave)
rm(list=ls())


`%nin%` <- Negate(`%in%`)

x.paralle <- T
if(x.paralle) {
  library(future)
  options(future.globals.maxSize = 50000 * 1024^2)
  plan("multicore", workers = 10)
  plan()
}


#=====================================================================
# Input / Output

rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.RDS"

output.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/"

if(!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}

mff.rds <- readRDS(rds.file)

print(mff.rds)
print(dim(mff.rds))
print(head(mff.rds@meta.data))


#=====================================================================
# Tissue 정보 생성 및 확인

mff.rds$tissue_group <- sub("__.*", "", mff.rds$orig.ident)
mff.rds$tissue_group <- factor(mff.rds$tissue_group, levels = c("Heart", "Limb_Muscle"))

print(table(mff.rds$tissue_group, useNA = "ifany"))
print(table(mff.rds$tissue_group, mff.rds$age))
print(table(mff.rds$tissue_group, mff.rds$mouse.id))


#=====================================================================
# 분석 설정

tissues <- c("Heart", "Limb_Muscle")

# 조직별로 사용할 PC 수
# 필요하면 Elbow plot 확인 후 각각 변경
dims.by.tissue <- c("Heart" = 18, "Limb_Muscle" = 18)

x.resolution <- 0.8
x.npcs <- 50

x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))


#=====================================================================
# Tissue별 독립 분석

for(tissue.name in tissues) {

  cat("\n====================================================\n")
  cat("Tissue:", tissue.name, "\n")

  x.dims <- dims.by.tissue[tissue.name]

  #--------------------------------------------------
  # Tissue별 output directory

  output.dir2 <- paste0(output.dir, tissue.name, "/")
  plot.dir <- paste0(output.dir2, "umap.plot/")

  if(!dir.exists(output.dir2)){dir.create(output.dir2, recursive = TRUE)}
  if(!dir.exists(plot.dir)){dir.create(plot.dir, recursive = TRUE)}


  #--------------------------------------------------
  # Tissue별 subset

  x <- subset(mff.rds, subset = tissue_group == tissue.name)

  x$tissue_group <- droplevels(x$tissue_group)

  cat("# of cells:", ncol(x), "\n")
  cat("# of genes:", nrow(x), "\n")

  print(table(x$age))
  print(table(x$mouse.id))
  print(table(x$orig.ident))
  print(table(x$cell_ontology_class))


  #--------------------------------------------------
  # Normalization

  DefaultAssay(x) <- "RNA"

  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)


  #--------------------------------------------------
  # Cell-cycle scoring

  tmp.s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(x))
  tmp.g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(x))

  print(c(S.genes = length(tmp.s.genes), G2M.genes = length(tmp.g2m.genes)))

  x <- CellCycleScoring(x, s.features = tmp.s.genes, g2m.features = tmp.g2m.genes, set.ident = TRUE)


  #--------------------------------------------------
  # Variable features / Scaling

  x <- FindVariableFeatures(x, selection.method = "vst")

  x <- ScaleData(
    object = x,
    features = VariableFeatures(x),
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score", "nCount_RNA")
  )


  #--------------------------------------------------
  # Variable feature plot

  top10 <- head(VariableFeatures(x), 10)

  plot.variable <- VariableFeaturePlot(x)
  plot.variable <- LabelPoints(plot = plot.variable, points = top10, repel = TRUE)

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".FindVariableGenes.png"), width = 900, height = 900)
  print(plot.variable)
  dev.off()


  #--------------------------------------------------
  # PCA

  x <- RunPCA(x, features = VariableFeatures(x), npcs = x.npcs)

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".PCElbowPlot.png"), width = 900, height = 900)
  print(ElbowPlot(x, ndims = x.npcs) + ggtitle(paste0(tissue.name, ": PCA Elbow Plot")))
  dev.off()


  #--------------------------------------------------
  # Neighbors / Clustering / UMAP
  # 아직 Harmony는 적용하지 않음

  x <- FindNeighbors(x, reduction = "pca", dims = 1:x.dims, force.recalc = TRUE)
  x <- FindClusters(x, resolution = x.resolution, random.seed = 1234)
  x <- RunUMAP(x, reduction = "pca", dims = 1:x.dims, seed.use = 1234)


  #--------------------------------------------------
  # Color 설정

  cluster.levels <- levels(x$seurat_clusters)
  cluster.colors <- setNames(x.getPalette(length(cluster.levels)), cluster.levels)

  sample.levels <- sort(unique(as.character(x$orig.ident)))
  sample.colors <- setNames(x.getPalette(length(sample.levels)), sample.levels)

  mouse.levels <- sort(unique(as.character(x$mouse.id)))
  mouse.colors <- setNames(x.getPalette(length(mouse.levels)), mouse.levels)

  age.levels <- unique(as.character(x$age))
  age.levels <- age.levels[!is.na(age.levels)]
  age.levels <- age.levels[order(as.numeric(gsub("[^0-9.]", "", age.levels)))]
  x$age_group <- factor(as.character(x$age), levels = age.levels)
  age.colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(age.levels)), age.levels)

  celltype.levels <- sort(unique(as.character(x$cell_ontology_class)))
  celltype.colors <- setNames(x.getPalette(length(celltype.levels)), celltype.levels)


  #--------------------------------------------------
  # 1. Cluster UMAP

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".UMAP.cluster_dim", x.dims, "_res", x.resolution, ".png"), width = 1400, height = 1200)
  print(DimPlot(x, reduction = "umap", group.by = "seurat_clusters", cols = cluster.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 8) + ggtitle(paste0(tissue.name, ": Seurat Clusters")))
  dev.off()


  #--------------------------------------------------
  # 2. Cell ontology class UMAP

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".UMAP.cell_ontology_class_dim", x.dims, "_res", x.resolution, ".png"), width = 1800, height = 1400)
  print(DimPlot(x, reduction = "umap", group.by = "cell_ontology_class", cols = celltype.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 6) + ggtitle(paste0(tissue.name, ": Cell Ontology Class")))
  dev.off()


  #--------------------------------------------------
  # 3. Age UMAP

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".UMAP.age_dim", x.dims, "_res", x.resolution, ".png"), width = 1400, height = 1200)
  print(DimPlot(x, reduction = "umap", group.by = "age_group", cols = age.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 7) + ggtitle(paste0(tissue.name, ": Age")))
  dev.off()


  #--------------------------------------------------
  # 4. Mouse UMAP
  # mouse.id = biological individual

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".UMAP.mouse_id_dim", x.dims, "_res", x.resolution, ".png"), width = 1600, height = 1300)
  print(DimPlot(x, reduction = "umap", group.by = "mouse.id", cols = mouse.colors, pt.size = 0.5) + ggtitle(paste0(tissue.name, ": Mouse ID")))
  dev.off()


  #--------------------------------------------------
  # 5. orig.ident UMAP
  # orig.ident = Cell Ranger technical run/library

  CairoPNG(filename = paste0(plot.dir, tissue.name, ".UMAP.origident_dim", x.dims, "_res", x.resolution, ".png"), width = 1800, height = 1400)
  print(DimPlot(x, reduction = "umap", group.by = "orig.ident", cols = sample.colors, pt.size = 0.5) + ggtitle(paste0(tissue.name, ": Technical Run")))
  dev.off()


  #--------------------------------------------------
  # Metadata별 세포 수 저장

  cell.count.summary <- x@meta.data %>%
    count(orig.ident, mouse.id, age, sex, cell_ontology_class, name = "n_cells") %>%
    arrange(age, mouse.id, cell_ontology_class)

  fwrite(cell.count.summary, paste0(output.dir2, tissue.name, ".cell_count_summary.csv"))


  #--------------------------------------------------
  # Tissue별 최종 RDS 저장

  saveRDS.pigz(
    x,
    paste0(output.dir2, "Aging.MFF.", tissue.name, ".seurat.normalization.pca.umap.RDS"),
    n.cores = 8
  )

  cat("Saved:", paste0(output.dir2, "Aging.MFF.", tissue.name, ".seurat.normalization.pca.umap.RDS"), "\n")

  rm(x)
  gc()
}


#=====================================================================
# 완료

cat("\nAll tissue-specific analyses completed.\n")
