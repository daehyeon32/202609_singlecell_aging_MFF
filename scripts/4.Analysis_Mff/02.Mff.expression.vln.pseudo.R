#=====================================================================
# 2.Mff.pseudobulk_cell_violin.R
# 새 R 세션에서 source()로 단독 실행 (Seurat 5 이상)
#
# 01: Mouse별 pseudobulk 점그래프 (mouse당 값 1개)
# 02: Age별 pseudobulk boxplot + 개별 mouse 점
# 03: 조직 전체 cell-level violin (x = age, y = Mff expression)
# 04: Cell type별 cell-level violin (x = age, y = Mff expression)
#
# Pseudobulk = log1p(10000 * sum(Mff raw UMI) / sum(all-gene raw UMI))
# Violin = 기존 RNA data의 세포별 LogNormalize 값
# Violin은 Seurat VlnPlot 사용; add.noise = FALSE; 흰 boxplot 없음
# Mff = 0인 세포 포함; 발현 cutoff 및 20-cell 조건 적용 안 함
# PCA / UMAP / Harmony / 기존 세포별 normalization 재실행 없음
# 입력 RDS 저장/덮어쓰기 없음; Heart와 Limb_Muscle 각각 저장
# 참고: https://satijalab.org/seurat/reference/aggregateexpression
#       https://satijalab.org/seurat/reference/normalizedata

library(Seurat)
library(Matrix)
library(Cairo)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(ggplot2)


#=====================================================================
# Input / Output

gene <- "Mff"
pseudobulk.scale.factor <- 10000
age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
age.colors <- setNames(brewer.pal(8, "Dark2")[seq_along(age.order)], age.order)

rds.files <- c(
  "Heart" = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Heart/Aging.MFF.Heart.seurat.normalization.pca.umap.RDS",
  "Limb_Muscle" = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Limb_Muscle/Aging.MFF.Limb_Muscle.seurat.normalization.pca.umap.RDS"
)

output.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression/"
output.subdir <- "pseudobulk_cell_violin/"

save.plot <- function(plot, filename, width, height) {
  CairoPNG(filename = filename, width = width, height = height)
  tryCatch(print(plot), finally = dev.off())
}


#=====================================================================
# Tissue별 독립 분석

for(tissue.name in names(rds.files)) {

  cat("\n====================================================\n")
  cat("Tissue:", tissue.name, "\n")
  rds.file <- rds.files[[tissue.name]]
  if(!file.exists(rds.file)){stop(paste0("RDS file이 없습니다: ", rds.file))}

  x <- readRDS(rds.file)
  if(!"RNA" %in% Assays(x)){stop("RNA assay가 없습니다.")}
  if(ncol(x) == 0){stop("RDS에 세포가 없습니다.")}
  DefaultAssay(x) <- "RNA"

  # Seurat v5에서 나뉜 counts/data layer만 메모리에서 연결
  # JoinLayers는 값을 이어 붙이며 normalization을 다시 수행하지 않음
  if(inherits(x[["RNA"]], "Assay5")) {
    for(layer.name in c("counts", "data")) {
      layer.names <- Layers(x[["RNA"]], search = paste0("^", layer.name, "($|\\.)"))
      if(length(layer.names) == 0){stop(paste0("RNA ", layer.name, " layer가 없습니다."))}
      if(length(layer.names) > 1 || layer.names[1] != layer.name){x <- JoinLayers(x, assay = "RNA", layers = layer.name, new = layer.name)}
    }
  }

  # Metadata 확인 및 나이순 정렬
  required.columns <- c("age", "mouse.id", "cell_ontology_class")
  if(!all(required.columns %in% colnames(x@meta.data))){stop("age / mouse.id / cell_ontology_class metadata를 확인하세요.")}
  cells <- colnames(x)
  cell.data <- x@meta.data[cells, required.columns, drop = FALSE]
  if(anyNA(cell.data$age) || anyNA(cell.data$mouse.id) || any(trimws(as.character(cell.data$mouse.id)) == "")){stop("age 또는 mouse.id에 결측값이 있습니다.")}
  if(any(!as.character(cell.data$age) %in% age.order)){stop("age.order에 없는 나이가 있습니다.")}

  cell.data$age_group <- factor(as.character(cell.data$age), levels = age.order)
  cell.data$age <- NULL
  cell.data$mouse.id <- as.character(cell.data$mouse.id)
  cell.data$cell_ontology_class <- as.character(cell.data$cell_ontology_class)
  cell.data$cell_ontology_class[is.na(cell.data$cell_ontology_class) | trimws(cell.data$cell_ontology_class) == ""] <- "Unannotated"
  cell.data$cell <- cells
  cell.data$tissue <- tissue.name
  x$age_group <- cell.data$age_group

  mouse.order <- cell.data %>% distinct(mouse.id, age_group) %>% arrange(age_group, mouse.id)
  if(anyDuplicated(mouse.order$mouse.id)){stop("동일한 mouse.id에 여러 age가 연결되어 있습니다.")}
  mouse.levels <- mouse.order$mouse.id
  mouse.labels <- setNames(paste0(mouse.order$mouse.id, "\n", mouse.order$age_group), mouse.levels)
  mouse.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(mouse.levels)), mouse.levels)
  cell.data$mouse.id <- factor(cell.data$mouse.id, levels = mouse.levels)


  #-------------------------------------------------------------------
  # Pseudobulk용 raw UMI와 violin용 기존 data를 각각 추출
  # 분모는 RNA counts에 있는 모든 gene; Mff만 남긴 matrix로 계산하지 않음

  rna.counts <- LayerData(x, assay = "RNA", layer = "counts")
  if(!gene %in% rownames(rna.counts)){stop(paste0("RNA counts에 ", gene, "가 없습니다."))}
  mff.data <- LayerData(x, assay = "RNA", layer = "data", features = gene)
  if(!gene %in% rownames(mff.data)){stop(paste0("RNA data에 ", gene, "가 없습니다."))}
  if(anyDuplicated(colnames(rna.counts)) || anyDuplicated(colnames(mff.data)) || !setequal(cells, colnames(rna.counts)) || !setequal(cells, colnames(mff.data))){stop("counts / data / metadata의 전체 세포 구성이 일치하지 않습니다.")}

  cell.data$mff_raw_umi <- as.numeric(rna.counts[gene, cells])
  cell.data$total_raw_umi <- as.numeric(Matrix::colSums(rna.counts)[cells])
  cell.data$expression <- as.numeric(mff.data[gene, cells])
  if(any(!is.finite(cell.data$mff_raw_umi)) || any(!is.finite(cell.data$total_raw_umi)) || any(!is.finite(cell.data$expression))){stop("발현값에 NA / Inf가 있습니다.")}
  if(any(cell.data$mff_raw_umi < 0) || any(cell.data$total_raw_umi <= 0) || any(cell.data$mff_raw_umi > cell.data$total_raw_umi) || any(cell.data$expression < 0)){stop("RNA counts / data 값의 범위를 확인하세요.")}

  # VlnPlot에 사용할 x는 유지; 원본 RDS는 저장하지 않음
  rm(rna.counts, mff.data)
  invisible(gc())


  #-------------------------------------------------------------------
  # Mouse별 pseudobulk
  # sum(세포별 모든 gene UMI) = 모든 gene의 pseudobulk count 합계
  # Mff만 분석하므로 큰 gene x mouse matrix 대신 필요한 합계만 계산
  # mouse별/age별 mean(expression)으로 pseudobulk를 계산하지 않음

  pseudobulk.mouse <- cell.data %>%
    group_by(tissue, mouse.id, age_group) %>%
    summarise(n_cells = n(), mff_umi_sum = sum(mff_raw_umi), total_umi_sum = sum(total_raw_umi), .groups = "drop") %>%
    mutate(mff_umi_fraction = mff_umi_sum / total_umi_sum, mff_per_10000_umi = pseudobulk.scale.factor * mff_umi_fraction, pseudobulk_expression = log1p(mff_per_10000_umi)) %>%
    arrange(age_group, mouse.id)

  # Age는 mouse별 pseudobulk 값을 묶는 기준; age 전체 count를 합치지 않음
  pseudobulk.age <- pseudobulk.mouse %>%
    group_by(tissue, age_group) %>%
    summarise(n_mice = n(), n_cells = sum(n_cells), mean_pseudobulk = mean(pseudobulk_expression), median_pseudobulk = median(pseudobulk_expression), q25_pseudobulk = quantile(pseudobulk_expression, 0.25), q75_pseudobulk = quantile(pseudobulk_expression, 0.75), sd_pseudobulk = sd(pseudobulk_expression), .groups = "drop") %>%
    mutate(sem_pseudobulk = sd_pseudobulk / sqrt(n_mice)) %>%
    arrange(age_group)

  # Violin 요약은 mouse 평균이 아니라 같은 age의 모든 cell을 사용
  # 존재하지 않는 age / cell type 조합은 행을 생성하거나 0으로 채우지 않음
  cell.age.summary <- cell.data %>%
    group_by(tissue, age_group) %>%
    summarise(n_cells = n(), n_mice = n_distinct(mouse.id), detected_cells = sum(expression > 0), percent_detected = mean(expression > 0) * 100, mean_expression = mean(expression), median_expression = median(expression), min_expression = min(expression), max_expression = max(expression), .groups = "drop") %>%
    arrange(age_group)

  celltype.age.summary <- cell.data %>%
    group_by(tissue, cell_ontology_class, age_group) %>%
    summarise(n_cells = n(), n_mice = n_distinct(mouse.id), detected_cells = sum(expression > 0), percent_detected = mean(expression > 0) * 100, mean_expression = mean(expression), median_expression = median(expression), min_expression = min(expression), max_expression = max(expression), .groups = "drop") %>%
    arrange(cell_ontology_class, age_group)

  output.dir2 <- paste0(output.dir, tissue.name, "/", output.subdir)
  if(!dir.exists(output.dir2)){dir.create(output.dir2, recursive = TRUE)}
  fwrite(pseudobulk.mouse, paste0(output.dir2, "01.", gene, ".pseudobulk_mouse.csv"))
  fwrite(pseudobulk.age, paste0(output.dir2, "02.", gene, ".pseudobulk_age_summary.csv"))
  fwrite(cell.data, paste0(output.dir2, "03.", gene, ".cell_level_expression.csv"))
  fwrite(cell.age.summary, paste0(output.dir2, "04.", gene, ".cell_age_summary.csv"))
  fwrite(celltype.age.summary, paste0(output.dir2, "05.", gene, ".celltype_age_summary.csv"))
  print(pseudobulk.mouse)


  #-------------------------------------------------------------------
  # 01. Sample별 pseudobulk: mouse당 값이 하나이므로 점으로 표시

  plot.mouse <- ggplot(pseudobulk.mouse, aes(x = mouse.id, y = pseudobulk_expression, color = age_group)) +
    geom_point(size = 4) +
    scale_x_discrete(limits = mouse.levels, labels = mouse.labels) +
    scale_color_manual(values = age.colors, breaks = age.order, drop = FALSE) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.1))) +
    labs(x = "Mouse ID (ordered by age)", y = paste0(gene, " pseudobulk: ln(1 + UMI per 10,000)"), color = "Age", title = paste0(tissue.name, ": ", gene, " pseudobulk per mouse")) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  save.plot(plot.mouse, paste0(output.dir2, "01.", gene, ".pseudobulk_mouse.png"), width = max(1600, length(mouse.levels) * 100), height = 1200)


  #-------------------------------------------------------------------
  # 02. Age별 pseudobulk boxplot: 점 1개 = mouse 1마리
  # Mouse가 1마리인 age는 점만 표시; boxplot 중앙선은 median

  pseudobulk.box.data <- pseudobulk.mouse %>% group_by(age_group) %>% filter(n() >= 2) %>% ungroup()

  plot.age <- ggplot(pseudobulk.mouse, aes(x = age_group, y = pseudobulk_expression)) +
    geom_boxplot(data = pseudobulk.box.data, outlier.shape = NA, width = 0.6, fill = "grey90") +
    geom_point(aes(color = mouse.id), position = position_jitter(width = 0.12, height = 0, seed = 1234), size = 4) +
    scale_x_discrete(limits = age.order, drop = FALSE) +
    scale_color_manual(values = mouse.colors, breaks = mouse.levels) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.03, 0.1))) +
    labs(x = "Age", y = paste0(gene, " pseudobulk: ln(1 + UMI per 10,000)"), color = "Mouse ID", title = paste0(tissue.name, ": ", gene, " pseudobulk by age"), subtitle = "One point per mouse; boxes summarize mice within each age") +
    theme_bw(base_size = 18) +
    theme(plot.title = element_text(hjust = 0.5))

  save.plot(plot.age, paste0(output.dir2, "02.", gene, ".pseudobulk_age_boxplot.png"), width = 1600, height = 1200)


  #-------------------------------------------------------------------
  # 03. 조직 전체 cell-level violin
  # x축 왼쪽에서 오른쪽으로 1m -> 3m -> 18m -> 21m -> 24m -> 30m
  # VlnPlot으로 단순하게 표시; add.noise = FALSE로 그림용 노이즈 끄기
  # 세포 1개 또는 값이 모두 같은 그룹은 실제 값에 점도 표시
  # 0을 포함한 모든 세포 사용; 흰 boxplot은 겹쳐 그리지 않음

  age.constant <- cell.age.summary %>% filter(min_expression == max_expression)
  expression.y.max <- max(1, max(cell.data$expression))

  plot.cell <- VlnPlot(x, features = gene, assay = "RNA", layer = "data", pt.size = 0, group.by = "age_group", cols = age.colors, add.noise = FALSE, y.max = expression.y.max, combine = FALSE)[[1]] +
    geom_point(data = age.constant, aes(x = age_group, y = mean_expression), inherit.aes = FALSE, size = 3) +
    scale_x_discrete(limits = age.order, drop = FALSE) +
    labs(x = "Age", y = paste0(gene, " expression (RNA data; LogNormalize)"), title = paste0(tissue.name, ": cell-level ", gene, " expression"), subtitle = "All cells, including zeros; cells pooled within each age") +
    NoLegend()

  save.plot(plot.cell, paste0(output.dir2, "03.", gene, ".age_cell_violin.png"), width = 1600, height = 1200)


  #-------------------------------------------------------------------
  # 04. Cell type별 cell-level violin (cell_ontology_class 사용)
  # 발현 cutoff / mouse-cell type당 20-cell 필터 모두 적용하지 않음
  # 모든 panel의 발현 축 범위를 공유; 자료가 없는 조합은 빈칸

  x$age_group <- factor(as.character(x$age), levels = age.order)
  celltype.y.max <- max(1, max(cell.data$expression))
  celltype.levels <- sort(unique(cell.data$cell_ontology_class))
  celltype.plots <- lapply(celltype.levels, function(celltype.name) {
    celltype.cells <- cell.data$cell[cell.data$cell_ontology_class == celltype.name]
    x.celltype <- subset(x, cells = celltype.cells)
    celltype.constant <- celltype.age.summary %>% filter(cell_ontology_class == celltype.name, min_expression == max_expression)

    VlnPlot(x.celltype, features = gene, assay = "RNA", layer = "data", pt.size = 0, group.by = "age_group", cols = age.colors, add.noise = FALSE, y.max = celltype.y.max, combine = FALSE)[[1]] +
      geom_point(data = celltype.constant, aes(x = age_group, y = mean_expression), inherit.aes = FALSE, size = 2) +
      scale_x_discrete(limits = age.order, drop = FALSE) +
      labs(x = "Age", y = paste0(gene, " expression"), title = paste(strwrap(celltype.name, width = 28), collapse = "\n")) +
      NoLegend()
  })

  plot.celltype <- patchwork::wrap_plots(celltype.plots, ncol = 3) + patchwork::plot_annotation(title = paste0(tissue.name, ": cell-level ", gene, " expression by cell type"), subtitle = "All cells, including zeros; no expression cutoff or minimum-cell filter")
  n.celltypes <- length(celltype.levels)
  save.plot(plot.celltype, paste0(output.dir2, "04.", gene, ".celltype_age_cell_violin.png"), width = 2400, height = max(1400, ceiling(n.celltypes / 3) * 650))
  writeLines(capture.output(sessionInfo()), paste0(output.dir2, "sessionInfo.txt"))
  cat("Saved:", output.dir2, "\n")
  rm(x)
  invisible(gc())
}
