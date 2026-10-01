# Young = 3m / Old = 18m, 21m, 24m. 1m과 30m은 분석 전에 제외합니다.
# RNA assay에 합쳐진 counts와 data(LogNormalize) layer가 있다고 가정합니다.
# 세포 검정은 세포 간 마우스 의존성을 보정하지 않습니다.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

# 1. 설정 --------------------------------------------------------------------
project.dir <- "/BiO/Live/dleogus32/202609Aging_MFF"
rds.file <- file.path(
  project.dir, "Analysis/3.preprocessing",
  "Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"
)
output.base <- file.path(project.dir, "Analysis/10.Mff_expression_new")
analysis.folder <- "young_old_sensitivity_3m_vs_18_21_24m_allow_single_mouse"
gene <- "Mff"
assay <- "RNA"
tissues <- c("Heart", "Limb_Muscle")
tissue.column <- "tissue_free_annotation"
celltype.column <- "cell_ontology_class"
mouse.column <- "mouse.id"
age.column <- "age"
include.heart.and.aorta <- TRUE
young.ages <- "3m"
old.ages <- c("18m", "21m", "24m")
group.order <- c("Young", "Old")
group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
min.cells.per.mouse.celltype <- 1
min.mice.per.group.wilcoxon <- 1  # 3번 원본처럼 그룹당 1마리도 검정합니다.
min.cells.per.group.wilcoxon <- 2
cell.fc.pseudocount <- 1
pvalue.column <- "p_value"       # violin/boxplot: p_value 또는 p_adj_BH
show.pvalue.on.violin <- TRUE
show.pvalue.on.boxplot <- TRUE
save.overviews <- TRUE
save.dotplot.pdf <- TRUE

plot.width <- 900
plot.height <- 850
plot.res <- 160
plot.base.size <- 13
plot.title.size <- 20
plot.axis.title.size <- 13
plot.axis.text.size <- 30
plot.title.wrap.width <- 40
plot.violin.width <- 0.8
plot.violin.box.width <- 0.18
plot.violin.box.alpha <- 0.85
plot.violin.box.linewidth <- 0.45
plot.pvalue.size <- 10

dotplot.width <- 10.5            # inch
dotplot.height.minimum <- 5.5
dotplot.height.per.celltype <- 0.48
dotplot.res <- 300
dotplot.title.size <- 22
dotplot.axis.title.size <- 18
dotplot.axis.number.size <- 17
dotplot.celltype.text.size <- 16
dotplot.legend.size <- 13
dotplot.point.size <- 4.5
dotplot.point.stroke <- 1.2
dotplot.p.cutoff <- 0.05         # dotplot은 항상 raw p를 사용합니다.
celltype.order <- "logFC"        # logFC 또는 alphabetical
x.limit.manual <- NULL          # NULL = 자동, 숫자 = 대칭 x축 범위
show.logFC.numbers <- FALSE
plot.logFC.number.size <- 4
direction.colors <- c(
  "Old lower" = "#D95F02", "Old higher" = "#000000", "No change" = "#777777"
)

# 2. 공통 함수 ---------------------------------------------------------------
write.table.csv <- function(d, file) {
  write.csv(d, file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

save.plot <- function(p, file, width = plot.width, height = plot.height,
                      res = plot.res, pdf = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  Cairo::CairoPNG(paste0(file, ".png"), width = width, height = height,
                  res = res, bg = "white")
  tryCatch(print(p), finally = grDevices::dev.off())
  if(pdf) ggsave(
    paste0(file, ".pdf"), p, device = grDevices::cairo_pdf,
    width = width / res, height = height / res, units = "in", limitsize = FALSE
  )
}

format.p <- function(p) {
  prefix <- if(pvalue.column == "p_value") "p" else "BH p"
  if(!is.finite(p)) return(paste0(prefix, " = NA"))
  if(p == 0) return(paste0(
    prefix, " < ", formatC(.Machine$double.xmin, format = "e", digits = 2)
  ))
  value <- formatC(p, format = if(p < 0.001) "e" else "f",
                  digits = if(p < 0.001) 2 else 3)
  paste0(prefix, " = ", value)
}

plot.theme <- function() {
  theme_bw(base_size = plot.base.size) + theme(
    legend.position = "none",
    plot.title = element_text(size = plot.title.size, face = "bold"),
    axis.title = element_text(size = plot.axis.title.size),
    axis.text = element_text(size = plot.axis.text.size),
    panel.grid.minor = element_blank()
  )
}

# 동일한 Wilcoxon 함수로 cell/mouse 검정. 검정 불가 사유만 짧게 남깁니다.
test.expression <- function(d, unit) {
  young <- d$expression[d$group == "Young"]
  old <- d$expression[d$group == "Old"]
  n1 <- length(young)
  n2 <- length(old)
  minimum <- if(unit == "cell") min.cells.per.group.wilcoxon else
    min.mice.per.group.wilcoxon
  r <- data.frame(
    test_unit = unit, n_Young = n1, n_Old = n2,
    mean_Young = if(n1) mean(young) else NA_real_,
    mean_Old = if(n2) mean(old) else NA_real_,
    median_Young = if(n1) median(young) else NA_real_,
    median_Old = if(n2) median(old) else NA_real_,
    W = NA_real_, p_value = NA_real_, exact_used = NA,
    status = "not_tested", reason = ""
  )
  if(unit == "mouse" && anyDuplicated(d$mouse.id)) stop("Repeated mouse.id.")
  if(min(n1, n2) == 0) r$reason <- "missing_group" else
    if(min(n1, n2) < minimum) r$reason <- "insufficient_observations" else
      if(length(unique(c(young, old))) < 2) r$reason <- "all_values_identical"
  if(nzchar(r$reason)) return(r)
  r$exact_used <- unit == "mouse" && max(n1, n2) < 50 &&
    !anyDuplicated(c(young, old))
  fit <- wilcox.test(
    young, old, alternative = "two.sided", paired = FALSE,
    exact = r$exact_used, correct = !r$exact_used
  )
  if(is.finite(fit$p.value)) {
    r$W <- unname(fit$statistic)
    r$p_value <- fit$p.value
    r$status <- "tested"
  } else r$reason <- "nonfinite_p"
  r
}

expression.plot <- function(d, unit, pvalue, title) {
  is.cell <- unit == "cell"
  show.p <- if(is.cell) show.pvalue.on.violin else show.pvalue.on.boxplot
  top <- max(1, d$expression)
  ymax <- top * if(show.p) 1.42 else if(is.cell) 1.08 else 1.12
  p <- ggplot(d, aes(group, expression)) + geom_blank()
  if(is.cell) {
    regular <- d %>% group_by(group) %>%
      filter(n() >= 2, n_distinct(expression) > 1) %>% ungroup()
    special <- d %>% filter(!group %in% regular$group) %>%
      group_by(group) %>% summarise(value = first(expression), n = n(), .groups = "drop")
    if(nrow(regular)) p <- p +
      geom_violin(data = regular, aes(fill = group), width = plot.violin.width,
                  trim = TRUE, scale = "width") +
      geom_boxplot(data = regular, width = plot.violin.box.width,
                   fill = "white", alpha = plot.violin.box.alpha,
                   linewidth = plot.violin.box.linewidth, outlier.shape = NA) +
      scale_fill_manual(values = group.colors)
    if(any(special$n > 1)) p <- p + geom_errorbar(
      data = special[special$n > 1, ], aes(group, ymin = value, ymax = value),
      inherit.aes = FALSE, width = 0.45, color = "grey20"
    )
    if(any(special$n == 1)) p <- p + geom_point(
      data = special[special$n == 1, ], aes(group, value),
      inherit.aes = FALSE, size = 2
    )
  } else {
    boxes <- d %>% group_by(group) %>% filter(n() >= 2) %>% ungroup()
    if(nrow(boxes)) p <- p + geom_boxplot(
      data = boxes, width = 0.55, fill = "grey92", outlier.shape = NA
    )
    if(nrow(d)) p <- p +
      geom_point(aes(color = group), size = 3,
                 position = position_jitter(width = 0.08, height = 0, seed = 1234)) +
      scale_color_manual(values = group.colors, drop = FALSE)
  }
  if(show.p) p <- p +
    annotate("segment", x = 1, xend = 2, y = top * 1.10,
             yend = top * 1.10, linewidth = 0.4) +
    annotate("segment", x = c(1, 2), xend = c(1, 2), y = top * 1.065,
             yend = top * 1.10, linewidth = 0.4) +
    annotate("text", x = 1.5, y = top * 1.26, label = format.p(pvalue),
             size = plot.pvalue.size, lineheight = 1)
  p + scale_x_discrete(limits = group.order, drop = FALSE) +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(title = title, x = NULL, y = paste0(
      gene, if(is.cell) " expression (RNA LogNormalize)" else
        " pseudobulk log2(normalized CPM + 1)"
    )) + plot.theme()
}

# 3. 데이터 준비 -------------------------------------------------------------
mff.rds <- readRDS(rds.file)
DefaultAssay(mff.rds) <- assay
meta <- mff.rds@meta.data[colnames(mff.rds), , drop = FALSE]
tissue.labels <- as.character(meta[[tissue.column]])
if(include.heart.and.aorta) {
  tissue.labels[tissue.labels %in% "Heart_and_Aorta"] <- "Heart"
}
celltype.labels <- as.character(meta[[celltype.column]])
celltype.labels[is.na(celltype.labels) | trimws(celltype.labels) == ""] <- "Unannotated"
all.cells <- data.frame(
  cell = colnames(mff.rds), tissue = tissue.labels, celltype = celltype.labels,
  mouse.id = as.character(meta[[mouse.column]]), age = as.character(meta[[age.column]])
)
dotplot.results <- list()

for(tissue.name in tissues) {
  before <- all.cells[all.cells$tissue %in% tissue.name, ]
  celltypes <- sort(unique(before$celltype))  # 제외 후 빈 세포형도 유지
  cells <- before[before$age %in% c(young.ages, old.ages), ]
  if(!nrow(cells)) stop("No selected cells: ", tissue.name)
  mouse.ages <- unique(cells[, c("mouse.id", "age")])
  if(anyNA(cells$mouse.id) || any(trimws(cells$mouse.id) == "") ||
     anyDuplicated(mouse.ages$mouse.id)) stop("Check mouse.id and age: ", tissue.name)
  cells$age <- factor(cells$age, levels = c(young.ages, old.ages))
  cells$group <- factor(
    ifelse(cells$age %in% young.ages, "Young", "Old"), levels = group.order
  )
  tissue.rds <- subset(mff.rds, cells = cells$cell)
  layers <- SeuratObject::Layers(tissue.rds[[assay]])
  if(!all(c("counts", "data") %in% layers)) stop("RNA counts/data layers required.")
  counts <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "counts")
  counts <- counts[, cells$cell, drop = FALSE]
  rna.data <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "data")
  cells$expression <- as.numeric(rna.data[gene, cells$cell])
  if(any(!is.finite(cells$expression)) || any(cells$expression < 0)) {
    stop("Invalid RNA data: ", tissue.name)
  }

  # 4. 조직 × 세포형 × 마우스 raw counts 합산 및 TMM ----------------------------
  mice <- cells %>% group_by(tissue, celltype, mouse.id, age, group) %>%
    summarise(n_cells = n(), .groups = "drop") %>% arrange(celltype, group, age, mouse.id)
  mice$pb_id <- sprintf("PB%05d", seq_len(nrow(mice)))
  mapping <- cells %>% dplyr::select(cell, celltype, mouse.id) %>%
    left_join(dplyr::select(mice, celltype, mouse.id, pb_id), by = c("celltype", "mouse.id"))
  tissue.rds$pseudobulk_id <- setNames(mapping$pb_id, mapping$cell)
  pb <- AggregateExpression(
    tissue.rds, assays = assay, features = rownames(counts),
    group.by = "pseudobulk_id", return.seurat = FALSE, verbose = FALSE
  )[[assay]]
  pb <- pb[rownames(counts), mice$pb_id, drop = FALSE]
  mice$gene <- gene
  mice$gene_raw_count <- as.numeric(pb[gene, ])
  mice$library_size <- as.numeric(Matrix::colSums(pb))
  mice$eligible <- mice$n_cells >= min.cells.per.mouse.celltype & mice$library_size > 0
  mice$norm_factor <- mice$normalized_CPM <- mice$log2_CPM_plus1 <- NA_real_
  for(ct in celltypes) {
    idx <- which(mice$celltype == ct & mice$eligible)
    if(!length(idx)) next
    mat <- as.matrix(pb[, mice$pb_id[idx], drop = FALSE])
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    y <- edgeR::DGEList(counts = mat)
    if(length(idx) >= 2) y <- edgeR::calcNormFactors(y, method = "TMM") else
      y$samples$norm.factors <- 1
    cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
    mff.cpm <- if(gene %in% rownames(cpm)) as.numeric(cpm[gene, ]) else rep(0, length(idx))
    mice$norm_factor[idx] <- y$samples$norm.factors
    mice$normalized_CPM[idx] <- mff.cpm
    mice$log2_CPM_plus1[idx] <- log2(mff.cpm + 1)
  }
  if(any(!is.finite(mice$log2_CPM_plus1[mice$eligible]))) stop("Invalid mouse CPM.")

  # 5. 통계와 그림: cell/mouse 공통 반복문 -------------------------------------
  output.dir <- file.path(output.base, tissue.name, analysis.folder)
  violin.dir <- file.path(output.dir, "cell_level_violin")
  pb.dir <- file.path(output.dir, paste0("mouse_pseudobulk_TMM_min",
                                      min.cells.per.mouse.celltype))
  for(folder in c(violin.dir, pb.dir)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }
  file.map <- data.frame(
    celltype = celltypes,
    file_stem = sprintf("%03d.%s", seq_along(celltypes),
                        substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100))
  )
  write.table.csv(cells, file.path(violin.dir, "01.cell_expression.csv"))
  write.table.csv(file.map, file.path(violin.dir, "03.celltype_file_map.csv"))
  write.table.csv(mapping, file.path(pb.dir, "01.cell_to_pseudobulk.csv"))
  saveRDS(pb, file.path(pb.dir, "02.pseudobulk_raw_counts.RDS"), compress = "gzip")
  write.table.csv(mice, file.path(pb.dir, paste0("03.", gene, ".pseudobulk_mouse.csv")))
  write.table.csv(file.map, file.path(pb.dir, "06.celltype_file_map.csv"))

  for(unit in c("cell", "mouse")) {
    d.all <- if(unit == "cell") cells else mice[mice$eligible, ]
    if(unit == "mouse") d.all$expression <- d.all$log2_CPM_plus1
    results <- lapply(celltypes, function(ct) {
      d <- d.all[d.all$celltype == ct, ]
      s <- test.expression(d, unit)
      s$tissue <- tissue.name
      s$celltype <- ct
      s$gene <- gene
      s$n_mice_Young <- n_distinct(d$mouse.id[d$group == "Young"])
      s$n_mice_Old <- n_distinct(d$mouse.id[d$group == "Old"])
      s$n_cells_Young <- if(unit == "cell") s$n_Young else
        sum(d$n_cells[d$group == "Young"])
      s$n_cells_Old <- if(unit == "cell") s$n_Old else
        sum(d$n_cells[d$group == "Old"])
      s$single_mouse_in_group <- s$n_mice_Young == 1 || s$n_mice_Old == 1
      if(unit == "cell") {
        young <- expm1(d$expression[d$group == "Young"])
        old <- expm1(d$expression[d$group == "Old"])
        s$mean_normalized_Young <- if(length(young)) mean(young) else NA_real_
        s$mean_normalized_Old <- if(length(old)) mean(old) else NA_real_
        s$fc_pseudocount <- cell.fc.pseudocount
        # 원본 FoldChange의 mean.fxn과 같은 식. p 계산 여부와 독립입니다.
        s$logFC_Old_vs_Young <- log2(s$mean_normalized_Old + cell.fc.pseudocount) -
          log2(s$mean_normalized_Young + cell.fc.pseudocount)
        s$FC_Old_vs_Young <- 2^s$logFC_Old_vs_Young
      }
      s
    })
    statistics <- bind_rows(results)
    ok <- is.finite(statistics$p_value)
    statistics$p_adj_BH <- NA_real_
    statistics$p_adj_BH[ok] <- p.adjust(statistics$p_value[ok], method = "BH")
    statistics$analysis_id <- analysis.folder
    statistics$young_ages <- paste(young.ages, collapse = ";")
    statistics$old_ages <- paste(old.ages, collapse = ";")
    folder <- if(unit == "cell") violin.dir else pb.dir
    filename <- if(unit == "cell") {
      paste0("05.", gene, ".Young_vs_Old.Wilcoxon.pvalues.csv")
    } else paste0("08.", gene, ".Young_vs_Old.pvalues.csv")
    write.table.csv(statistics, file.path(folder, filename))
    if(unit == "cell") dotplot.results[[tissue.name]] <- statistics

    plots <- list()
    kind <- if(unit == "cell") "violin" else "mouse_boxplot"
    subfolder <- if(unit == "cell") "violin" else "boxplot"
    for(i in seq_along(celltypes)) {
      ct <- celltypes[i]
      d <- d.all[d.all$celltype == ct, ]
      title <- paste(strwrap(paste(tissue.name, ct, sep = " | "),
                             width = plot.title.wrap.width), collapse = "\n")
      plots[[i]] <- expression.plot(d, unit, statistics[[pvalue.column]][i], title)
      stem <- paste0(file.map$file_stem[i], ".Young_Old.", kind)
      save.plot(plots[[i]], file.path(folder, subfolder, stem))
    }
    if(save.overviews) {
      pages <- split(seq_along(plots), ceiling(seq_along(plots) / 9))
      for(page in seq_along(pages)) {
        idx <- pages[[page]]
        nc <- min(3, length(idx))
        stem <- sprintf("%s.Young_Old.%s.page%02d", gene, kind, page)
        save.plot(
          wrap_plots(plots[idx], ncol = nc), file.path(folder, "overview", stem),
          width = plot.width * nc, height = plot.height * ceiling(length(idx) / nc)
        )
      }
    }
  }
  message(tissue.name, ": violin / boxplot / statistics saved")
  rm(tissue.rds, counts, rna.data, pb, cells, mice, d.all, d, plots)
  invisible(gc())
}

# 6. Cell-mean log2FC dotplot: 두 조직이 같은 x축을 사용 ------------------------
effects <- unlist(lapply(dotplot.results, function(d) d$logFC_Old_vs_Young))
max.effect <- max(c(0, abs(effects[is.finite(effects)])))
x.limit <- max(0.25, max.effect * if(show.logFC.numbers) 1.4 else 1.15)
if(!is.null(x.limit.manual)) {
  if(!is.finite(x.limit.manual) || x.limit.manual <= 0 || x.limit.manual < max.effect) {
    stop("x.limit.manual must be positive and include all log2FC values.")
  }
  x.limit <- x.limit.manual
}
x.breaks <- pretty(c(-x.limit, x.limit), n = 5)
x.breaks <- sort(unique(c(0, x.breaks[abs(x.breaks) <= x.limit])))
significance.labels <- c(
  "Below cutoff" = paste0("Unadjusted p < ", dotplot.p.cutoff),
  "At or above cutoff" = paste0("Unadjusted p >= ", dotplot.p.cutoff), "Not tested" = "p = NA"
)

for(tissue.name in tissues) {
  d <- dotplot.results[[tissue.name]]
  d$plot_log2FC <- d$logFC_Old_vs_Young
  d$estimate_available <- is.finite(d$plot_log2FC)
  idx <- if(celltype.order == "logFC")
    order(!d$estimate_available, d$plot_log2FC, d$celltype, na.last = TRUE) else
      order(d$celltype)
  d <- d[idx, ]
  d$celltype_axis <- factor(d$celltype, levels = rev(d$celltype))
  d$plot_direction <- ifelse(d$plot_log2FC < 0, "Old lower",
                             ifelse(d$plot_log2FC > 0, "Old higher", "No change"))
  d$plot_significance <- ifelse(!is.finite(d$p_value), "Not tested",
                                ifelse(d$p_value < dotplot.p.cutoff,
                                       "Below cutoff", "At or above cutoff"))
  estimated <- d[d$estimate_available, ]
  unavailable <- d[!d$estimate_available, ]
  p <- ggplot(d, aes(y = celltype_axis)) +
    geom_vline(xintercept = 0, color = "#999999", linewidth = 0.6, linetype = "dashed")
  if(nrow(estimated)) p <- p +
    geom_point(data = estimated,
               aes(x = plot_log2FC, color = plot_direction, shape = plot_significance),
               size = dotplot.point.size, stroke = dotplot.point.stroke) +
    scale_color_manual(values = direction.colors, guide = "none") +
    scale_shape_manual(
      values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4),
      limits = names(significance.labels), labels = unname(significance.labels),
      drop = FALSE, name = NULL
    ) + guides(shape = guide_legend(
      nrow = 1, override.aes = list(color = "#444444", size = 4)
    ))
  if(nrow(unavailable)) p <- p + geom_text(
    data = unavailable, x = x.limit * 0.96, label = "NA",
    hjust = 1, color = "#888888", size = 4.5
  )
  if(show.logFC.numbers && nrow(estimated)) p <- p + geom_text(
    data = estimated, aes(x = plot_log2FC, label = sprintf("%+.2f", plot_log2FC),
                         hjust = ifelse(plot_log2FC < 0, 1.3, -0.3)),
    size = plot.logFC.number.size, color = "#444444", show.legend = FALSE
  )
  p <- p +
    scale_x_continuous(breaks = x.breaks, limits = c(-x.limit, x.limit),
                       expand = expansion(mult = 0.02)) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = 0.7)) +
    labs(title = paste(gsub("_", " ", tissue.name), gene, sep = " | "),
         x = "Cell-mean log2 fold change (Old / Young)", y = NULL) +
    theme_classic(base_size = 14, base_family = "sans") + theme(
      plot.title = element_text(size = dotplot.title.size, face = "bold",
                                margin = margin(b = 14)),
      axis.title.x = element_text(size = dotplot.axis.title.size, margin = margin(t = 12)),
      axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"),
      axis.text.y = element_text(size = dotplot.celltype.text.size, color = "black",
                                 margin = margin(r = 10)),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.35),
      legend.position = "bottom", legend.text = element_text(size = dotplot.legend.size),
      legend.key.width = grid::unit(1, "cm"), plot.margin = margin(15, 25, 12, 12)
    )
  folder <- file.path(output.base, tissue.name, analysis.folder,
                      "cell_level_violin", "log2FC_dotplot")
  height <- max(dotplot.height.minimum, 2 + nrow(d) * dotplot.height.per.celltype)
  stem <- paste0(gene, ".Young_vs_Old.cell_mean.log2FC")
  save.plot(p, file.path(folder, stem), width = dotplot.width * dotplot.res,
            height = height * dotplot.res, res = dotplot.res, pdf = save.dotplot.pdf)
  d$celltype_axis <- NULL
  write.table.csv(d, file.path(folder, paste0(stem, ".plot_data.csv")))
  message(tissue.name, ": log2FC dotplot saved")
}
