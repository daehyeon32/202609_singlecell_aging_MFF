# 저장된 조직별 RDS로 그림만 다시 저장합니다.
# 정규화, PCA, clustering, UMAP 좌표를 다시 계산하거나 RDS를 저장하지 않습니다.
library(Seurat)
library(Cairo)
library(RColorBrewer)
library(ggplot2)

input.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split"

#=====================================================================
# 글씨 크기 설정: 이 부분의 숫자를 직접 바꾼 뒤 파일 전체를 다시 실행하세요.
# 아래 font.size 값은 모두 pt입니다. 내부 라벨에 필요한 mm 변환은 자동 처리합니다.
# x.text.scale: 1 = 아래 설정 그대로, 1.2 = 모든 글씨 20% 확대, 0.8 = 20% 축소

x.text.scale <- 1
font.size <- list()
font.size$base <- 24           # 별도 지정하지 않은 기본 글씨
font.size$title <- 36          # 그림 제목: Heart: Seurat Clusters 등
font.size$subtitle <- 24       # 부제목 (있는 경우)
font.size$caption <- 24        # 그림 하단 설명 (있는 경우)
font.size$tag <- 36            # 패널 번호 A, B 등 (있는 경우)

font.size$x_title <- 40        # X축 제목: umap_1
font.size$y_title <- 40        # Y축 제목: umap_2
font.size$x_tick <- 40        # X축 눈금 숫자: -15, -10, -5, 0 등
font.size$y_tick <- 40         # Y축 눈금 숫자: -10, -5, 0, 5 등

font.size$legend_title <- 26   # 범례 제목 (있는 경우)
font.size$legend_text <- 30    # 범례 항목: 0, 1, 2 또는 cell type / sample 이름
font.size$strip_x <- 24        # 가로 방향 패널 제목 (있는 경우)
font.size$strip_y <- 24        # 세로 방향 패널 제목 (있는 경우)

font.size$cluster <- 40        # Cluster UMAP 내부 숫자
font.size$celltype <- 40       # Cell type UMAP 내부 세포 유형 이름
font.size$age <- 40            # Age UMAP 내부 나이 라벨
font.size$gene <- 40           # Variable feature plot 내부 유전자 이름

#=====================================================================
# 설정 적용: 아래는 크기를 변경할 때 수정할 필요 없습니다.
# element_text는 pt, DimPlot/LabelPoints의 내부 라벨은 mm를 사용합니다.
# https://ggplot2.tidyverse.org/reference/element.html
# https://ggplot2.tidyverse.org/reference/geom_text.html

font.pt <- function(key) {font.size[[key]] * x.text.scale}
font.mm <- function(key) {font.pt(key) / ggplot2::.pt}

plot.text.theme <- theme(text = element_text(size = font.pt("base")))
plot.text.theme <- plot.text.theme + theme(plot.title = element_text(size = font.pt("title"), face = "bold", hjust = 0.5), plot.subtitle = element_text(size = font.pt("subtitle")), plot.caption = element_text(size = font.pt("caption")), plot.tag = element_text(size = font.pt("tag")))
plot.text.theme <- plot.text.theme + theme(axis.title.x = element_text(size = font.pt("x_title")), axis.title.y = element_text(size = font.pt("y_title")), axis.text.x = element_text(size = font.pt("x_tick")), axis.text.y = element_text(size = font.pt("y_tick")))
plot.text.theme <- plot.text.theme + theme(legend.title = element_text(size = font.pt("legend_title")), legend.text = element_text(size = font.pt("legend_text")), legend.key.height = grid::unit(font.pt("legend_text") * 1.1, "pt"))
plot.text.theme <- plot.text.theme + theme(strip.text.x = element_text(size = font.pt("strip_x")), strip.text.y = element_text(size = font.pt("strip_y")), plot.margin = margin(15, 15, 15, 15))
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

save.png <- function(p, filename, width, height) {
  CairoPNG(filename = filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  print(p)
  cat("Saved: ", filename, "\n", sep = "")
}

for(tissue.name in c("Heart", "Limb_Muscle")) {
  rds.file <- file.path(input.dir, tissue.name, paste0("Aging.MFF.", tissue.name, ".seurat.normalization.pca.umap.RDS"))
  if(!file.exists(rds.file)){stop("RDS 파일이 없습니다: ", rds.file)}
  cat("\nReading: ", rds.file, "\n", sep = "")
  x <- readRDS(rds.file)
  if(!"umap" %in% names(x@reductions)){stop(tissue.name, ": 저장된 UMAP 좌표가 없습니다.")}
  required.columns <- c("seurat_clusters", "cell_ontology_class", "age", "mouse.id", "orig.ident")
  if(!all(required.columns %in% colnames(x@meta.data))){stop(tissue.name, ": 필요한 metadata 열이 없습니다: ", paste(setdiff(required.columns, colnames(x@meta.data)), collapse = ", "))}

  DefaultAssay(x) <- "RNA"
  plot.dir <- file.path(input.dir, tissue.name, "umap.plot", "large_text")
  dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

  # 기존 코드와 같은 색상 구성
  cluster.levels <- levels(x$seurat_clusters)
  if(is.null(cluster.levels)){cluster.levels <- sort(unique(as.character(x$seurat_clusters)))}
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

  # 각 DimPlot에서 단일 ggplot을 받아 제목/축/범례 theme을 직접 적용합니다.
  # https://satijalab.org/seurat/reference/dimplot
  # 1. Cluster UMAP
  p <- DimPlot(x, reduction = "umap", group.by = "seurat_clusters", cols = cluster.colors, pt.size = 0.5, label = TRUE, repel = TRUE, label.size = font.mm("cluster"), combine = FALSE)[[1]] + ggtitle(paste0(tissue.name, ": Seurat Clusters")) + plot.text.theme
  save.png(p, file.path(plot.dir, paste0(tissue.name, ".UMAP.cluster.large_text.png")), width = 1400, height = 1200)

  # 2. Cell type UMAP
  p <- DimPlot(x, reduction = "umap", group.by = "cell_ontology_class", cols = celltype.colors, pt.size = 0.5, label = FALSE, repel = TRUE, label.size = font.mm("celltype"), combine = FALSE)[[1]] + ggtitle(paste0(tissue.name, ": Cell Ontology Class")) + plot.text.theme
  save.png(p, file.path(plot.dir, paste0(tissue.name, ".UMAP.cell_ontology_class.large_text.png")), width = 1800, height = 1400)

  # 3. Age UMAP
  p <- DimPlot(x, reduction = "umap", group.by = "age_group", cols = age.colors, pt.size = 0.5, label = FALSE, repel = TRUE, label.size = font.mm("age"), combine = FALSE)[[1]] + ggtitle(paste0(tissue.name, ": Age")) + plot.text.theme
  save.png(p, file.path(plot.dir, paste0(tissue.name, ".UMAP.age.large_text.png")), width = 1400, height = 1200)

  # 4. Mouse UMAP
  p <- DimPlot(x, reduction = "umap", group.by = "mouse.id", cols = mouse.colors, pt.size = 0.5, combine = FALSE)[[1]] + ggtitle(paste0(tissue.name, ": Mouse ID")) + plot.text.theme
  save.png(p, file.path(plot.dir, paste0(tissue.name, ".UMAP.mouse_id.large_text.png")), width = 1600, height = 1300)

  # 5. Technical run UMAP
  p <- DimPlot(x, reduction = "umap", group.by = "orig.ident", cols = sample.colors, pt.size = 0.5, combine = FALSE)[[1]] + ggtitle(paste0(tissue.name, ": Technical Run")) + plot.text.theme
  save.png(p, file.path(plot.dir, paste0(tissue.name, ".UMAP.origident.large_text.png")), width = 1800, height = 1400)

  # 저장된 결과로 variable feature plot과 elbow plot도 다시 그립니다.
  if(length(VariableFeatures(x)) > 0) {
    top10 <- head(VariableFeatures(x), 10)
    p <- VariableFeaturePlot(x)
    p <- LabelPoints(plot = p, points = top10, repel = TRUE, size = font.mm("gene")) + plot.text.theme
    save.png(p, file.path(plot.dir, paste0(tissue.name, ".FindVariableGenes.large_text.png")), width = 900, height = 900)
  }
  if("pca" %in% names(x@reductions)) {
    p <- ElbowPlot(x, ndims = min(50, length(x[["pca"]]@stdev))) + ggtitle(paste0(tissue.name, ": PCA Elbow Plot")) + plot.text.theme
    save.png(p, file.path(plot.dir, paste0(tissue.name, ".PCElbowPlot.large_text.png")), width = 900, height = 900)
  }

  rm(x, p)
  gc(verbose = FALSE)
}

cat("\nReplotting completed.\n")
