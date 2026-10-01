library(Seurat)
library(Cairo)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
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

gene <- "Mff"
min.cells.per.mouse.celltype <- 20
mff.expression.cutoff <- 0.5

# 발현값: RNA data의 LogNormalize 값 (0인 세포도 평균에 포함)
# 08/09 그래프: mouse별 평균을 먼저 구한 뒤 나이별 동일 가중 평균 ± SEM
# 같은 mouse의 시간 경과가 아니라 서로 다른 나이군 사이의 발현 추이
# 20-cell 조건은 cell type별 mouse 그래프와 평균 추이에만 적용
age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
age.colors <- setNames(brewer.pal(8, "Dark2")[seq_along(age.order)], age.order)
stopifnot(length(mff.expression.cutoff) == 1, is.finite(mff.expression.cutoff), mff.expression.cutoff >= 0)

rds.files <- c(
  "Heart" = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Heart/Aging.MFF.Heart.seurat.normalization.pca.umap.RDS",
  "Limb_Muscle" = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Limb_Muscle/Aging.MFF.Limb_Muscle.seurat.normalization.pca.umap.RDS"
)

output.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression/"
if(!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}


#=====================================================================
# Tissue별 Mff expression 분석

for(tissue.name in names(rds.files)) {

  cat("\n====================================================\n")
  cat("Tissue:", tissue.name, "\n")

  rds.file <- rds.files[tissue.name]
  output.dir2 <- paste0(output.dir, tissue.name, "/")
  if(!dir.exists(output.dir2)){dir.create(output.dir2, recursive = TRUE)}

  if(!file.exists(rds.file)){stop(paste0("RDS file이 없습니다: ", rds.file))}

  x <- readRDS(rds.file)
  DefaultAssay(x) <- "RNA"

  if(gene %nin% rownames(x)){stop(paste0(gene, " gene이 RDS에 없습니다."))}

  # Age 축 및 Mouse ID 축/범례를 모두 나이순으로 고정
  stopifnot(all(c("age", "mouse.id", "cell_ontology_class") %in% colnames(x@meta.data)))
  if(anyNA(x$age) || anyNA(x$mouse.id) || any(as.character(x$mouse.id) == "")){stop("age 또는 mouse.id의 결측값을 확인하세요.")}
  if(any(as.character(x$age) %nin% age.order)){stop("age.order에 없는 나이가 있습니다. age.order와 age.colors를 확인하세요.")}
  age.levels <- age.order[age.order %in% as.character(x$age)]
  x$age_group <- factor(as.character(x$age), levels = age.levels)

  mouse.order <- x@meta.data %>% transmute(mouse.id = as.character(mouse.id), age_group) %>% distinct() %>% arrange(age_group, mouse.id)
  if(anyDuplicated(mouse.order$mouse.id)){stop("동일한 mouse.id에 여러 age가 연결되어 있습니다.")}
  mouse.levels <- mouse.order$mouse.id
  x$mouse.id <- factor(as.character(x$mouse.id), levels = mouse.levels)
  mouse.colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(mouse.levels)), mouse.levels)

  print(table(x$age_group))
  print(table(x$cell_ontology_class))
  print(table(x$mouse.id, x$age_group))


  #-------------------------------------------------------------------
  # Cell 단위 expression table

  cell.data <- FetchData(x, vars = c(gene, "age_group", "mouse.id", "cell_ontology_class"), layer = "data", clean = FALSE)
  colnames(cell.data)[colnames(cell.data) == gene] <- "expression"
  if(nrow(cell.data) != ncol(x) || any(!is.finite(cell.data$expression))){stop("전체 세포의 Mff 발현값을 가져왔는지 확인하세요.")}
  cell.data$age_group <- factor(as.character(cell.data$age_group), levels = age.levels)
  cell.data$mouse.id <- factor(as.character(cell.data$mouse.id), levels = mouse.levels)
  cell.data$cell_ontology_class <- as.character(cell.data$cell_ontology_class)
  cell.data$cell_ontology_class[is.na(cell.data$cell_ontology_class) | cell.data$cell_ontology_class == ""] <- "Unannotated"
  cell.data$cell <- rownames(cell.data)
  cell.data$tissue <- tissue.name
  cell.data$detected <- cell.data$expression > 0
  cell.data$above_cutoff <- cell.data$expression > mff.expression.cutoff
  cell.data$age_months <- as.numeric(sub("m$", "", as.character(cell.data$age_group)))


  #-------------------------------------------------------------------
  # Cell type × Age 요약

  celltype.age.summary <- cell.data %>%
    group_by(tissue, cell_ontology_class, age_group) %>%
    summarise(n_cells = n(), n_mice = n_distinct(mouse.id), mean_expression = mean(expression), median_expression = median(expression), percent_expressed = mean(detected) * 100, .groups = "drop") %>%
    arrange(age_group, cell_ontology_class)

  fwrite(celltype.age.summary, paste0(output.dir2, "01.", gene, ".celltype_age_summary.csv"))


  #-------------------------------------------------------------------
  # Mouse × Cell type × Age 요약

  mouse.celltype.summary <- cell.data %>%
    group_by(tissue, mouse.id, age_group, cell_ontology_class) %>%
    summarise(n_cells = n(), mean_expression = mean(expression), median_expression = median(expression), percent_expressed = mean(detected) * 100, .groups = "drop") %>%
    mutate(analysis_eligible = n_cells >= min.cells.per.mouse.celltype, age_months = as.numeric(sub("m$", "", as.character(age_group)))) %>%
    arrange(age_group, mouse.id, cell_ontology_class)

  fwrite(mouse.celltype.summary, paste0(output.dir2, "02.", gene, ".mouse_celltype_age_summary.csv"))


  #-------------------------------------------------------------------
  # Mouse × Age 조직 전체 요약

  mouse.age.summary <- cell.data %>%
    group_by(tissue, mouse.id, age_group) %>%
    summarise(n_cells = n(), mean_expression = mean(expression), median_expression = median(expression), percent_expressed = mean(detected) * 100, .groups = "drop") %>%
    mutate(age_months = as.numeric(sub("m$", "", as.character(age_group)))) %>%
    arrange(age_group, mouse.id)

  fwrite(mouse.age.summary, paste0(output.dir2, "03.", gene, ".mouse_age_summary.csv"))


  #-------------------------------------------------------------------
  # Mouse별 Mff 검출 및 cutoff 통과 세포 요약

  mouse.detection.summary <- cell.data %>%
    group_by(tissue, mouse.id, age_group) %>%
    summarise(total_cells = n(), detected_cells = sum(detected), not_detected_cells = sum(!detected), above_cutoff_cells = sum(above_cutoff), at_or_below_cutoff_cells = sum(!above_cutoff), 
      percent_detected = mean(detected) * 100, percent_not_detected = mean(!detected) * 100, percent_above_cutoff = mean(above_cutoff) * 100, percent_above_cutoff_among_detected = ifelse(sum(detected) > 0, 
        sum(above_cutoff) / sum(detected) * 100, NA_real_), cutoff_log_normalized = mff.expression.cutoff, .groups = "drop") %>%
    arrange(age_group, mouse.id)

  fwrite(mouse.detection.summary, paste0(output.dir2, "04.", gene, ".mouse_detection_cutoff_summary.csv"))
  print(mouse.detection.summary)


  #-------------------------------------------------------------------
  # 나이별 평균 추이용 요약 (mouse마다 동일한 가중치)
  # mean_of_mouse_means: mouse별 mean_expression의 평균
  # SEM: mouse 평균 사이의 SD / sqrt(mouse 수); mouse 1개이면 NA
  # Cell-level 요약인 01 CSV와 달리 세포가 많은 mouse에 가중하지 않음

  age.mean.summary <- mouse.age.summary %>%
    group_by(tissue, age_group, age_months) %>%
    summarise(n_mice = n(), n_cells = sum(n_cells), mean_of_mouse_means = mean(mean_expression), sd_mouse_means = sd(mean_expression), .groups = "drop") %>%
    mutate(sem_mouse_means = sd_mouse_means / sqrt(n_mice)) %>%
    arrange(age_months)

  mouse.celltype.plot.data <- mouse.celltype.summary %>% filter(analysis_eligible)

  celltype.age.mean.summary <- mouse.celltype.plot.data %>%
    group_by(tissue, cell_ontology_class, age_group) %>%
    summarise(n_mice = n(), n_cells = sum(n_cells), mean_of_mouse_means = mean(mean_expression), sd_mouse_means = sd(mean_expression), .groups = "drop") %>%
    complete(tissue = tissue.name, cell_ontology_class = sort(unique(cell.data$cell_ontology_class)), age_group = factor(age.levels, levels = age.levels), fill = list(n_mice = 0L, n_cells = 0L)) %>%
    mutate(age_months = as.numeric(sub("m$", "", as.character(age_group))), sem_mouse_means = ifelse(n_mice > 1, sd_mouse_means / sqrt(n_mice), NA_real_)) %>%
    arrange(age_months, cell_ontology_class)

  fwrite(age.mean.summary, paste0(output.dir2, "05.", gene, ".age_mean_trend.csv"))
  fwrite(celltype.age.mean.summary, paste0(output.dir2, "06.", gene, ".celltype_age_mean_trend.csv"))


  #-------------------------------------------------------------------
  # 1. Mff FeaturePlot

  CairoPNG(filename = paste0(output.dir2, "01.", gene, ".FeaturePlot.png"), width = 1400, height = 1200)
  print(FeaturePlot(x, features = gene, reduction = "umap", cols = c("grey90", "red"), order = TRUE, pt.size = 0.5) + ggtitle(paste0(tissue.name, ": ", gene, " expression")))
  dev.off()


  #-------------------------------------------------------------------
  # 2. 조직 전체 나이별 Mff cell-level violin plot

  CairoPNG(filename = paste0(output.dir2, "02.", gene, ".age_cell_violin.png"), width = 1400, height = 1200)
  print(VlnPlot(x, features = gene, group.by = "age_group", cols = age.colors[age.levels], pt.size = 0) + geom_boxplot(width = 0.12, 
    outlier.shape = NA, fill = "white") + scale_x_discrete(limits = age.levels) + xlab("Age") + ggtitle(paste0(tissue.name, ": ", gene, " expression by age")))
  dev.off()


  #-------------------------------------------------------------------
  # 3. Cell type × Age DotPlot
  # Color: 평균 발현량 / Size: 발현 세포 비율

  n.celltypes <- length(unique(celltype.age.summary$cell_ontology_class))
  dotplot.height <- max(1200, n.celltypes * 75)

  plot.dot <- ggplot(celltype.age.summary, aes(x = age_group, y = cell_ontology_class)) +
    geom_point(aes(size = percent_expressed, color = mean_expression)) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    scale_size(range = c(1, 12), limits = c(0, 100)) +
    scale_x_discrete(limits = age.levels) +
    xlab("Age") +
    ylab("Cell ontology class") +
    ggtitle(paste0(tissue.name, ": ", gene, " expression by cell type and age")) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "03.", gene, ".celltype_age_dotplot.png"), width = 1600, height = dotplot.height)
  print(plot.dot)
  dev.off()


  #-------------------------------------------------------------------
  # 4. Mouse 단위 조직 전체 나이별 평균 Mff

  plot.mouse.age <- ggplot(mouse.age.summary, aes(x = age_group, y = mean_expression)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
    geom_point(aes(color = mouse.id), position = position_jitter(width = 0.12, height = 0, seed = 1234), size = 4) +
    scale_x_discrete(limits = age.levels) +
    scale_color_manual(values = mouse.colors, breaks = mouse.levels) +
    xlab("Age") +
    ylab(paste0("Mean ", gene, " expression per mouse")) +
    ggtitle(paste0(tissue.name, ": mouse-level ", gene, " expression")) +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "04.", gene, ".mouse_age_mean.png"), width = 1400, height = 1200)
  print(plot.mouse.age)
  dev.off()


  #-------------------------------------------------------------------
  # 5. Cell type별 mouse 단위 나이별 평균 Mff

  if(nrow(mouse.celltype.plot.data) > 0) {
  plot.mouse.celltype <- ggplot(mouse.celltype.plot.data, aes(x = age_group, y = mean_expression, group = age_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey90") +
    geom_point(aes(color = mouse.id), position = position_jitter(width = 0.12, height = 0, seed = 1234), size = 2) +
    scale_x_discrete(limits = age.levels) +
    scale_color_manual(values = mouse.colors, breaks = mouse.levels) +
    facet_wrap(~cell_ontology_class, scales = "free_y") +
    xlab("Age") +
    ylab(paste0("Mean ", gene, " expression per mouse")) +
    ggtitle(paste0(tissue.name, ": mouse-level ", gene, " expression by cell type")) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "05.", gene, ".mouse_celltype_age_mean.png"), width = 2000, height = 1600)
  print(plot.mouse.celltype)
  dev.off()
  }


  #-------------------------------------------------------------------
  # 6. Mouse별 Mff 검출/미검출 세포 수와 비율

  detection.plot.data <- bind_rows(
    mouse.detection.summary %>% transmute(mouse.id, age_group, status = "Mff detected", n_cells = detected_cells, percent = percent_detected),
    mouse.detection.summary %>% transmute(mouse.id, age_group, status = "Mff not detected", n_cells = not_detected_cells, percent = percent_not_detected)
  )

  detection.plot.data$status <- factor(detection.plot.data$status, levels = c("Mff not detected", "Mff detected"))

  plot.detection <- ggplot(detection.plot.data, aes(x = mouse.id, y = n_cells, fill = status)) +
    geom_col() +
    geom_text(aes(label = paste0(n_cells, "\n(", round(percent, 1), "%)")), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("Mff not detected" = "grey80", "Mff detected" = "#E64B35")) +
    scale_x_discrete(limits = mouse.levels) +
    xlab("Mouse ID (ordered by age)") +
    ylab("Number of cells") +
    ggtitle(paste0(tissue.name, ": ", gene, " detected and not detected cells")) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "06.", gene, ".mouse_detection_count.png"), width = 1800, height = 1200)
  print(plot.detection)
  dev.off()


  #-------------------------------------------------------------------
  # 7. Mouse별 cutoff 초과 세포 수와 전체 세포 대비 비율

  plot.cutoff <- ggplot(mouse.detection.summary, aes(x = mouse.id, y = percent_above_cutoff, fill = age_group)) +
    geom_col() +
    geom_text(aes(label = paste0(above_cutoff_cells, "/", total_cells, "\n(", round(percent_above_cutoff, 1), "%)")), vjust = -0.2, size = 4) +
    ylim(0, 105) +
    scale_x_discrete(limits = mouse.levels) +
    scale_fill_manual(values = age.colors, breaks = age.levels) +
    xlab("Mouse ID (ordered by age)") +
    ylab("Cells above cutoff (%)") +
    ggtitle(paste0(tissue.name, ": ", gene, " expression > ", mff.expression.cutoff)) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "07.", gene, ".mouse_above_cutoff.png"), width = 1800, height = 1200)
  print(plot.cutoff)
  dev.off()


  #-------------------------------------------------------------------
  # 6-2 / 7-2. Age별 Mff 검출 및 cutoff 통과 세포 요약
  # 같은 tissue, 같은 나이의 모든 mouse에서 세포 수를 합산
  # 비율 = 해당 조건의 세포 수 / 해당 나이의 전체 세포 수 * 100
  # mouse별 비율의 평균이 아님; 20-cell 조건은 적용하지 않음
  # 6개 나이 축을 유지하되 자료가 없는 나이는 막대를 그리지 않음

  age.detection.summary <- cell.data %>%
    group_by(tissue, age_group) %>%
    summarise(n_mice = n_distinct(mouse.id), total_cells = n(), detected_cells = sum(detected), not_detected_cells = sum(!detected), above_cutoff_cells = sum(above_cutoff),
     at_or_below_cutoff_cells = sum(!above_cutoff), percent_detected = mean(detected) * 100, percent_not_detected = mean(!detected) * 100, percent_above_cutoff = mean(above_cutoff) * 100,
      percent_above_cutoff_among_detected = ifelse(sum(detected) > 0, sum(above_cutoff) / sum(detected) * 100, NA_real_), cutoff_log_normalized = mff.expression.cutoff, .groups = "drop") %>%
    arrange(age_group)

  fwrite(age.detection.summary, paste0(output.dir2, "07.", gene, ".age_detection_cutoff_summary.csv"))


  #-------------------------------------------------------------------
  # 6-2. Age별 Mff 검출/미검출 세포 수와 비율
  # 검출: expression > 0 / 미검출: expression == 0

  age.detection.plot.data <- bind_rows(
    age.detection.summary %>% transmute(age_group, status = "Mff detected", n_cells = detected_cells, percent = percent_detected),
    age.detection.summary %>% transmute(age_group, status = "Mff not detected", n_cells = not_detected_cells, percent = percent_not_detected)
  )

  age.detection.plot.data$status <- factor(age.detection.plot.data$status, levels = c("Mff not detected", "Mff detected"))

  plot.age.detection <- ggplot(age.detection.plot.data, aes(x = age_group, y = n_cells, fill = status)) +
    geom_col() +
    geom_text(aes(label = paste0(n_cells, "\n(", round(percent, 1), "%)")), position = position_stack(vjust = 0.5), size = 4) +
    geom_text(data = age.detection.summary, aes(x = age_group, y = total_cells, label = paste0("mice=", n_mice)), inherit.aes = FALSE, vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Mff not detected" = "grey80", "Mff detected" = "#E64B35")) +
    scale_x_discrete(limits = age.order, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    xlab("Age") +
    ylab("Number of cells") +
    ggtitle(paste0(tissue.name, ": ", gene, " detected and not detected cells by age")) +
    labs(subtitle = "Pooled cells within each age; percentages use all cells in that age") +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "06_2.", gene, ".age_detection_count.png"), width = 1600, height = 1200)
  print(plot.age.detection)
  dev.off()


  #-------------------------------------------------------------------
  # 7-2. Age별 cutoff 초과 세포 수와 전체 세포 대비 비율
  # cutoff와 같은 값은 제외: expression > mff.expression.cutoff

  plot.age.cutoff <- ggplot(age.detection.summary, aes(x = age_group, y = percent_above_cutoff, fill = age_group)) +
    geom_col() +
    geom_text(aes(label = paste0(above_cutoff_cells, "/", total_cells, "\n(", round(percent_above_cutoff, 1), "%)")), vjust = -0.2, size = 4) +
    scale_x_discrete(limits = age.order, drop = FALSE) +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
    scale_fill_manual(values = age.colors, breaks = age.order) +
    xlab("Age") +
    ylab("Cells above cutoff (%)") +
    ggtitle(paste0(tissue.name, ": ", gene, " expression > ", mff.expression.cutoff, " by age")) +
    labs(subtitle = "Above-cutoff cells / all cells in each age; pooled across mice") +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "07_2.", gene, ".age_above_cutoff.png"), width = 1600, height = 1200)
  print(plot.age.cutoff)
  dev.off()


  #-------------------------------------------------------------------
  # 8. 조직 전체: 나이별 평균 Mff 발현 추이
  # 작은 점: 각 mouse / 큰 점 및 선: 나이별 평균 / 오차막대: SEM
  # x축은 실제 개월 간격; 연결선은 서로 다른 나이군 평균의 추이

  plot.age.trend <- ggplot(age.mean.summary, aes(x = age_months, y = mean_of_mouse_means)) +
    geom_point(data = mouse.age.summary, aes(x = age_months, y = mean_expression, color = age_group), inherit.aes = FALSE, position = position_jitter(width = 0.15, height = 0, seed = 1234), size = 3, alpha = 0.7) +
    geom_line(aes(group = 1), color = "grey30", linewidth = 0.9) +
    geom_errorbar(aes(ymin = mean_of_mouse_means - sem_mouse_means, ymax = mean_of_mouse_means + sem_mouse_means), width = 0.5, na.rm = TRUE) +
    geom_point(aes(color = age_group), size = 5) +
    geom_text(aes(label = paste0("n=", n_mice)), y = Inf, vjust = 1.4, size = 4) +
    scale_x_continuous(breaks = as.numeric(sub("m$", "", age.levels)), labels = age.levels) +
    scale_color_manual(values = age.colors, breaks = age.levels) +
    xlab("Age (months)") +
    ylab(paste0("Mean log-normalized ", gene, " expression")) +
    ggtitle(paste0(tissue.name, ": ", gene, " mean expression by age")) +
    labs(subtitle = "Equal-weight mouse means +/- SEM; n = mice (no SEM when n = 1)") +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "08.", gene, ".age_mean_trend.png"), width = 1500, height = 1200)
  print(plot.age.trend)
  dev.off()


  #-------------------------------------------------------------------
  # 9. Cell type별: 나이별 평균 Mff 발현 추이
  # 해당 mouse-cell type에 최소 20 cells인 그룹만 포함 (위 설정값)
  # 누락 그룹은 NA로 남겨 선을 끊음; 모든 panel은 동일한 y축

  if(nrow(mouse.celltype.plot.data) > 0) {
  plot.celltype.age.trend <- ggplot(celltype.age.mean.summary, aes(x = age_months, y = mean_of_mouse_means)) +
    geom_point(data = mouse.celltype.plot.data, aes(x = age_months, y = mean_expression, color = age_group), inherit.aes = FALSE, position = position_jitter(width = 0.15, height = 0, seed = 1234), size = 1.8, alpha = 0.6) +
    geom_line(aes(group = cell_ontology_class), color = "grey30", linewidth = 0.7, na.rm = TRUE) +
    geom_errorbar(aes(ymin = mean_of_mouse_means - sem_mouse_means, ymax = mean_of_mouse_means + sem_mouse_means), width = 0.5, na.rm = TRUE) +
    geom_point(aes(color = age_group), size = 3, na.rm = TRUE) +
    geom_text(aes(label = paste0("n=", n_mice)), y = Inf, vjust = 1.4, size = 3) +
    facet_wrap(~cell_ontology_class, ncol = 3, labeller = label_wrap_gen(width = 28)) +
    scale_x_continuous(breaks = as.numeric(sub("m$", "", age.levels)), labels = age.levels) +
    scale_color_manual(values = age.colors, breaks = age.levels) +
    xlab("Age (months)") +
    ylab(paste0("Mean log-normalized ", gene, " expression")) +
    ggtitle(paste0(tissue.name, ": ", gene, " mean expression by cell type and age")) +
    labs(subtitle = paste0("Mouse means +/- SEM; n = eligible mice; >= ", min.cells.per.mouse.celltype, " cells per mouse-cell type; n=0: no estimate")) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))

  CairoPNG(filename = paste0(output.dir2, "09.", gene, ".celltype_age_mean_trend.png"), width = 2100, height = max(1200, ceiling(n.celltypes / 3) * 420))
  print(plot.celltype.age.trend)
  dev.off()
  }

  rm(x, cell.data)
  gc()
}

cat("\nMff expression draft analysis completed.\n")
