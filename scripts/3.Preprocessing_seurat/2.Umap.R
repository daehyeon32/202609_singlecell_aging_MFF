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
output.dir = '/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/'

if(!dir.exists(output.dir)){dir.create(output.dir)}

x.paralle <- T
if(x.paralle) {
  library(future)
  options(future.globals.maxSize = 50000 * 1024^2)
  plan("multicore", workers = 10)
  plan()
}

#===========================================================
rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.RDS"
mff.rds <- readRDS(file = rds.file)

# data normalization
mff.rds <- NormalizeData(mff.rds, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS.pigz(mff.rds,  paste0(output.dir, '/Aging.MFF.seurat.metadata.filtered.normalization.RDS' ),n.cores = 8)


tmp.s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(mff.rds))
tmp.g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(mff.rds))
print(c(S.genes = length(tmp.s.genes), G2M.genes = length(tmp.g2m.genes)))
mff.rds <- CellCycleScoring(mff.rds, s.features = tmp.s.genes, g2m.features = tmp.g2m.genes, set.ident = TRUE)
# head(mff.rds@meta.data)

mff.rds <- FindVariableFeatures(object = mff.rds, selection.method = "vst") #mff.rds[["RNA"]]@var.features
mff.rds <- ScaleData(object = mff.rds, features = VariableFeatures(object = mff.rds), vars.to.regress = c( "percent.mt", "S.Score", "G2M.Score", 'nCount_RNA'))
#Warning: Requested variables to regress not in object: percent.ribo
#mff.rds[["RNA"]]@scale.data / percent.ribo?

if(!dir.exists(paste0(output.dir, '/umap.plot/'))){dir.create(paste0(output.dir, '/umap.plot/'))}
CairoPNG(filename = paste0(output.dir, "/umap.plot/FindVariableGenes.png"), width = 900, height = 900)
#print(VariableFeaturePlot(object = mff.rds))
top10 <- head(VariableFeatures(mff.rds), 10)
plot1 <- VariableFeaturePlot(mff.rds)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1)
dev.off()

#===========================================================
x.dims <- 50
mff.rds <- RunPCA(object = mff.rds, features = VariableFeatures(object = mff.rds), npcs = x.dims)

# Determine statistically significant principal components
#mff.rds <- JackStraw(object = mff.rds, num.replicate = 100, dims = x.dims)
#mff.rds <- ScoreJackStraw(object = mff.rds, dims = 1:x.dims)

CairoPNG(filename = paste0(output.dir, "/umap.plot/PCElbowPlot.png"), width = 900, height = 900)
print(ElbowPlot(object = mff.rds, ndims = x.dims))
dev.off()
# CairoPNG(filename = paste0(output.dir, "/umap.plot/JackStrawPlot.png"), width = 1800, height = 1200)
# JackStrawPlot(object = mff.rds, dims = 1:x.dims)
# dev.off()
# CairoPNG(filename = paste0(output.dir, "/umap.plot/pca.DoHeatmap.png"), width = 900, height = 5100)
# DimHeatmap(mff.rds, dims = 1:x.dims, cells = 1000, balanced = TRUE)
# dev.off()
#===========================================================

# PC elbow plot에서 stdev가 급격히 감소하는 구간 pc로 설정
# 각 PC의 stdev를 모두 합해서 정규화 (y axis_각 pc의 분산 비율 )
# PC의 분산을 차례로 더해 누적 분산(x axis)
# PC선정기준: min(누적 분산 >90% & 분산 <5%,  분산 감소폭 0.1 이상인 PC중 최대 idx)
# Determine percent of variation associated with each PC
pct <- mff.rds[["pca"]]@stdev /sum(mff.rds[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2) 
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
plot <- ggplot(plot_df, aes(cumu, pct, label= rank, color = rank > pcs)) +
  geom_text(size = 8) +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept= min(pct[pct > 5]), color = "grey") +
  xlab("Cumulative percents for each PC") +
  ylab("percent of variation associated with each PC") +
  theme_bw(base_size = 30) +
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30))
options(bitmapType="cairo")
CairoPNG(filename = paste0(output.dir, "/umap.plot/PCElbowPlot_select_dim.png"), width = 1500, height = 1000)
print(plot)
dev.off()

saveRDS.pigz(mff.rds, paste0(output.dir, '/Aging.MFF.seurat.metadata.filtered.normalization.pca.RDS' ),n.cores = 8)
#===========================================================
#Umap
x.dims <- c(18)
x.resolution <- 0.8 #c(0.2, 0.4, 0.6, 0.8, 1.0)

rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.RDS"
mff.rds <- readRDS(file = rds.file)
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

mff.rds <- FindNeighbors(object = mff.rds, dims = 1:x.dims, force.recalc = T)
mff.rds <- FindClusters(object = mff.rds, resolution = x.resolution)
mff.rds <- RunUMAP(object = mff.rds, dims = 1:x.dims)
x.sample.color <- x.getPalette(length(unique(mff.rds@meta.data$orig.ident)))
x.cluster.color <- x.getPalette(length(unique(mff.rds@meta.data$seurat_clusters)))

#Visualization
CairoPNG(filename = paste0(output.dir, "/umap.plot/MFF.UMAP.cell_dim",as.character(x.dims),'_res',as.character(x.resolution),".png"),  width = 1300, height = 1200)
print(DimPlot(object = mff.rds, reduction = "umap", label = T, pt.size = 1, label.size = 18, cols = x.cluster.color)) #, cols = x.cluster.color
dev.off()

#Color by patient
CairoPNG(filename = paste0(output.dir, "/umap.plot/MFF.UMAP.origident_dim",as.character(x.dims),'_res',as.character(x.resolution),".png"), width = 1400, height = 1200)
print(DimPlot(object = mff.rds, reduction = "umap", group.by = "orig.ident", pt.size = 1, label.size = 18, cols = x.sample.color))
dev.off()

saveRDS.pigz(mff.rds, paste0(output.dir, '/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS' ),n.cores = 8)
#--------------------------------------------------

# mff.rds <- readRDS('/BiO2/Research/UNIST-SEV-mff-2025/0_Analysis/1.preprocessing/mff.seraut.normalization.pca.umap.RDS')
# FeaturePlot(mff.rds, 'HBD', order = TRUE)
# var_list  <- VariableFeatures(mff.rds)
# var_list[grep('^HB',var_list)]


#===========================================================
# Additional UMAP plots: tissue, age, tissue-age

# Heart / Limb Muscle
mff.rds$tissue_group <- sub("__.*", "", mff.rds$orig.ident)
mff.rds$tissue_group <- factor(mff.rds$tissue_group, levels = c("Heart", "Limb_Muscle"))

# Age를 숫자 순서로 정렬
age.values <- as.character(mff.rds$age)
age.levels <- unique(age.values[!is.na(age.values)])
age.levels <- age.levels[order(as.numeric(gsub("[^0-9.]", "", age.levels)))]
mff.rds$age_group <- factor(age.values, levels = age.levels)

# Tissue와 Age 조합
observed.combinations <- unique(paste(mff.rds$tissue_group, mff.rds$age_group, sep = "_"))
tissue.age.levels <- unlist(lapply(age.levels, function(x) paste(c("Heart", "Limb_Muscle"), x, sep = "_")))
tissue.age.levels <- tissue.age.levels[tissue.age.levels %in% observed.combinations]
mff.rds$tissue_age <- factor(paste(mff.rds$tissue_group, mff.rds$age_group, sep = "_"), levels = tissue.age.levels)

# 색상 지정
tissue.colors <- c("Heart" = "#E64B35", "Limb_Muscle" = "#4DBBD5")
age.colors <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(age.levels)), age.levels)
tissue.age.colors <- setNames(brewer.pal(12, "Paired")[seq_along(tissue.age.levels)], tissue.age.levels)

# 실제 그룹 확인
print(table(mff.rds$tissue_group))
print(table(mff.rds$age_group))
print(table(mff.rds$tissue_age))

# 1. Heart vs Limb Muscle: 2색
CairoPNG(filename = paste0(output.dir, "/umap.plot/MFF.UMAP.tissue_dim", x.dims, "_res", x.resolution, ".png"), width = 1400, height = 1200)
print(DimPlot(object = mff.rds, reduction = "umap", group.by = "tissue_group", cols = tissue.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 8) + ggtitle("Tissue"))
dev.off()

# 2. Age: 최대 6색
CairoPNG(filename = paste0(output.dir, "/umap.plot/MFF.UMAP.age_dim", x.dims, "_res", x.resolution, ".png"), width = 1400, height = 1200)
print(DimPlot(object = mff.rds, reduction = "umap", group.by = "age_group", cols = age.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 7) + ggtitle("Age"))
dev.off()

# 3. Tissue × Age: 최대 12색
CairoPNG(filename = paste0(output.dir, "/umap.plot/MFF.UMAP.tissue_age_dim", x.dims, "_res", x.resolution, ".png"), width = 1600, height = 1200)
print(DimPlot(object = mff.rds, reduction = "umap", group.by = "tissue_age", cols = tissue.age.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 6) + ggtitle("Tissue and Age"))
dev.off()











 
