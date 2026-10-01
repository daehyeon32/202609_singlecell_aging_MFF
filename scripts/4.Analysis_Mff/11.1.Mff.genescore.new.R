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
library(harmony)
rm(list=ls())

library(patchwork)
library(edgeR)

x.paralle <- T
if(x.paralle) {
  library(future)
  options(future.globals.maxSize = 50000 * 1024^2)
  plan("multicore", workers = 50)
  plan()
}

if(packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0") {
  stop("Seurat and SeuratObject >= 5.0.0 required.")
}

rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"
mff.rds <- readRDS(rds.file)
output.dir = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/8.genescore_instant"
if(!dir.exists(output.dir)){dir.create(output.dir)}

#------------------------------------------------1.parameter-setting---------------------------------------------------------------------------

gene <- "Mff"
assay <- "RNA"
tissues <- c("Heart", "Limb_Muscle")
tissue.column <- "tissue_free_annotation"
celltype.column <- "cell_ontology_class"
mouse.column <- "mouse.id"
age.column <- "age"
age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
young.ages <- c("1m", "3m")
old.ages <- c("18m", "21m", "24m", "30m")
make.young.old.plots <- TRUE
make.age.dotplots <- TRUE

min.matched.genes <- 1
low.coverage.warning <- 0.50
min.cells.per.mouse.celltype <- 1  # No new cell-count cutoff by default.
min.mice.per.group.wilcoxon <- 2  # Latest supplied script's test-only rule.
min.cells.per.group.wilcoxon <- 2
min.observations.correlation <- 3  # Computational minimum, not reliability guarantee.
save.scored.seurat <- TRUE  # Writes a NEW tissue-specific RDS; never the input file.
save.pdf <- FALSE  # PNG default; TRUE also exports individual vector PDFs.
save.overviews <- TRUE

plot.width <- 900L
plot.height <- 850L
plot.res <- 160L
plot.base.size <- 13
plot.title.size <- 20
plot.axis.title.size <- 13
plot.axis.text.size <- 30
plot.pvalue.size <- 10  # ggplot2 text size in mm, matching the reference.
plot.title.wrap.width <- 40
plot.legend.size <- 13
plot.violin.width <- 0.8
plot.box.width <- 0.18
cell.point.size <- 0.55
cell.point.alpha <- 0.70  # 0 = transparent, 1 = opaque; increased for visibility.
mouse.point.size <- 3.5
show.regression.line <- TRUE  # TRUE: show a line in each All/Young/Old view; FALSE: hide all lines.
regression.line.color <- "#222222"
regression.line.width <- 1.0
scatter.width <- plot.width  # One plotting area per cell type in every scatter view.
dotplot.width.in <- 10.5
dotplot.height.minimum.in <- 5.5
dotplot.height.per.celltype.in <- 0.48
dotplot.res <- 300L
dotplot.base.size <- 14
dotplot.title.size <- 22
dotplot.axis.title.size <- 18
dotplot.axis.number.size <- 17
dotplot.celltype.text.size <- 16
dotplot.legend.size <- 13
dotplot.point.size <- 4.5
dotplot.point.stroke <- 1.2
dotplot.na.text.size <- 4.5
dotplot.p.cutoff <- 0.05 



include.heart.and.aorta <- TRUE  # Matches the latest supplied Mff script.
red.when.old.higher <- c("SaulSenMayo", "CoreScence_up")


age.colors <- c( "1m" = "#1B9E77", "3m" = "#D95F02", "18m" = "#7570B3", "21m" = "#E7298A", "24m" = "#66A61E", "30m" = "#E6AB02")
group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
gene.files <- c(
  Hallmark_IFN_alpha = "HALLMARK_INTERFERON_ALPHA_RESPONSE.v2026.1.Mm.grp",
  GO_IFN_beta = "GOBP_RESPONSE_TO_INTERFERON_BETA.v2026.1.Mm.grp",
  SaulSenMayo = "SAUL_SEN_MAYO.v2026.1.Mm.grp", CoreScence_up = "CoreScence_mouse_up_22.grp",
  CoreScence_down = "CoreScence_mouse_down_16.grp"
)
score.labels <- c(
  Hallmark_IFN_alpha = "Hallmark IFN-alpha", GO_IFN_beta = "GO IFN-beta", SaulSenMayo = "SaulSenMayo",
  CoreScence_up = "CoreScence up", CoreScence_down = "CoreScence down"
)
expected.sizes <- c( Hallmark_IFN_alpha = 94, GO_IFN_beta = 77, SaulSenMayo = 117, CoreScence_up = 22, CoreScence_down = 16)

score.columns <- setNames(paste0("AMS_", names(gene.files)), names(gene.files))

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


#----------------------------------------------------2.read-geneset----------------------------------------------------------------------

geneset.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Workspace/0.Meta/geneset"

read.grp <- function(file) {
  genes <- readLines(file, warn = FALSE, encoding = "UTF-8")
  genes <- sub("^\ufeff", "", genes) #BOM문자 지우기 (눈에 안보임)
  genes <- trimws(genes) #각 줄의 앞뒤 공백 지우기
  genes <- genes[nzchar(genes)] #빈줄 제거 
  genes <- unique(genes) #주석 제거
  return(genes)
}

genesets <- lapply(file.path(geneset.dir, gene.files), read.grp)
names(genesets) <- names(gene.files)

lengths(genesets)

#----------------------------------------------------3.data-preparation-----------------------------------------------------------------

DefaultAssay(mff.rds) <- assay

tissue.labels <- as.character(mff.rds@meta.data[[tissue.column]])

if(include.heart.and.aorta) {tissue.labels[tissue.labels %in% "Heart_and_Aorta"] <- "Heart"}

tissue.list <- list()

for(tissue.name in tissues) {
  cells <- colnames(mff.rds)[tissue.labels %in% tissue.name]
  if(length(cells) == 0) {stop("No cells for tissue: ", tissue.name)}

  tissue.rds <- subset(mff.rds, cells = cells)
  tissue.list[[tissue.name]] <- tissue.rds
}

sapply(tissue.list, ncol)

#----------------------------------------------------4.Addmodulescore-------------------------------------------------------------------

for(tissue.name in tissues) {
  tissue.rds <- tissue.list[[tissue.name]]
  rna.genes <- rownames(SeuratObject::LayerData(tissue.rds[[assay]], layer = "data"))

  features <- lapply(genesets, function(g) setdiff(intersect(g, rna.genes), gene))
  pool <- setdiff(rna.genes, gene)

  if(any(lengths(features) < min.matched.genes)) {stop("Insufficient matched genes in ", tissue.name)}

  message(tissue.name)
  print(lengths(features))

  temporary.names <- paste0("AgingMFF_TMP", seq_along(features))

  tissue.rds <- AddModuleScore(
    object = tissue.rds,
    features = unname(features),
    pool = pool,
    nbin = 24,
    ctrl = 100,
    assay = assay,
    name = "AgingMFF_TMP",
    seed = 1,
    search = FALSE,
    slot = "data"
  )

  for(i in seq_along(features)) {
    values <- tissue.rds@meta.data[colnames(tissue.rds), temporary.names[i]]
    tissue.rds[[score.columns[[i]]]] <- values
    tissue.rds[[temporary.names[i]]] <- NULL
  }

  tissue.list[[tissue.name]] <- tissue.rds
}

#----------------------------------------------------5.data-preparation-for-figure------------------------------------------------------
cell.data.list <- list()
mouse.data.list <- list()

for(tissue.name in tissues) {
  tissue.rds <- tissue.list[[tissue.name]]
  meta <- tissue.rds@meta.data[colnames(tissue.rds), , drop = FALSE]

  # 5-1. 세포별 metadata, Mff 발현값, score
  cells <- data.frame(cell = colnames(tissue.rds), tissue = tissue.name)
  cells$celltype <- as.character(meta[[celltype.column]])
  cells$mouse.id <- as.character(meta[[mouse.column]])
  cells$age <- as.character(meta[[age.column]])
  cells$celltype[is.na(cells$celltype) | trimws(cells$celltype) == ""] <- "Unannotated"

  rna.data <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "data")
  cells$Mff_LogNormalize <- as.numeric(rna.data[gene, cells$cell])
  cells[, unname(score.columns)] <- meta[, unname(score.columns), drop = FALSE]

  # 5-2. 연령과 Young/Old

  cells$age_months <- as.numeric(sub("m$", "", cells$age))
  cells$group <- factor(ifelse(cells$age %in% young.ages, "Young", "Old"), levels = c("Young", "Old"))
  cells$age <- factor(cells$age, levels = age.order)

  # 5-3. 조직 × 세포형 × 마우스별 평균 score
  mice <- cells %>%
    group_by(tissue, celltype, mouse.id, age, age_months, group) %>%
    summarise(n_cells = n(), across(all_of(unname(score.columns)), mean), .groups = "drop") %>%
    arrange(celltype, age_months, mouse.id)

  # 5-4. 같은 단위로 raw count 합산 후 세포형별 TMM 정규화
  mice$pb_id <- sprintf("PB%05d", seq_len(nrow(mice)))

  mapping <- cells %>%
    dplyr::select(cell, celltype, mouse.id) %>%
    left_join(dplyr::select(mice, celltype, mouse.id, pb_id), by = c("celltype", "mouse.id"))

  tissue.rds$AgingMFF_pb_id <- setNames(mapping$pb_id, mapping$cell)
  counts <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "counts")
  counts <- counts[, cells$cell, drop = FALSE]
  if(any(Matrix::colSums(counts) <= 0)) {stop("Zero-count cells in ", tissue.name)}
  if(!gene %in% rownames(counts)) {stop("Missing ", gene, " in RNA counts of ", tissue.name)}

  pb <- AggregateExpression(tissue.rds, assays = assay, features = rownames(counts), group.by = "AgingMFF_pb_id", return.seurat = FALSE, verbose = FALSE)[[assay]]
  pb <- pb[rownames(counts), mice$pb_id, drop = FALSE]

  mice$library_size <- as.numeric(Matrix::colSums(pb))
  mice$eligible <- mice$n_cells >= min.cells.per.mouse.celltype & mice$library_size > 0
  mice$Mff_log2_TMM_CPM_plus1 <- NA_real_

  for(ct in unique(mice$celltype)) {
    idx <- which(mice$celltype == ct & mice$eligible)
    if(length(idx) == 0) {next}

    mat <- as.matrix(pb[, mice$pb_id[idx], drop = FALSE])
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    y <- edgeR::DGEList(counts = mat)

    if(length(idx) >= 2) {
      y <- edgeR::calcNormFactors(y, method = "TMM")
    } else {
      y$samples$norm.factors <- 1
    }

    cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
    mff.cpm <- if(gene %in% rownames(cpm)) as.numeric(cpm[gene, ]) else rep(0, length(idx))
    mice$Mff_log2_TMM_CPM_plus1[idx] <- log2(mff.cpm + 1)
  }

  # 5-5. 조직별 데이터 보관
  cell.data.list[[tissue.name]] <- cells
  mouse.data.list[[tissue.name]] <- mice
}
#----------------------------------------------------6.making-figures-------------------------------------------------------------------

# 6-1. 공통 표시·저장 함수
format.p <- function(p) {
  if(!is.finite(p)) return("p = NA")
  if(p == 0) return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 1)))
  paste0("p = ", if(p < 0.001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
}

plot.theme <- function(legend = FALSE) {
  theme_bw(base_size = plot.base.size) +
    theme(
      plot.title = element_text(size = plot.title.size, face = "bold"),
      axis.title = element_text(size = plot.axis.title.size),
      axis.text.x = element_text(size = plot.axis.text.size, color = "black"),
      axis.text.y = element_text(size = plot.axis.text.size, color = "black"),
      panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none",
      legend.title = element_text(size = plot.legend.size),
      legend.text = element_text(size = plot.legend.size),
      plot.margin = margin(12, 16, 12, 12)
    )
}

plot.theme.dot <- function(legend = FALSE) {
  theme_bw(base_size = dotplot.base.size) +
    theme(
      plot.title = element_text(size = dotplot.title.size, face = "bold"),
      axis.title = element_text(size = dotplot.axis.title.size),
      axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"),
      axis.text.y = element_text(size = dotplot.celltype.text.size, color = "black"),
      panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none",
      legend.title = element_text(size = dotplot.legend.size),
      legend.text = element_text(size = dotplot.legend.size),
      plot.margin = margin(25, 16, 12, 12)
    )
}

y.bounds <- function(y, top = 0.30) {
  if(!length(y)) return(c(-0.1, 0.1))
  b <- range(y)
  span <- max(diff(b), abs(mean(b)) * 0.10, 0.05)
  c(b[1] - span * 0.06, b[2] + span * top)
}

save.plot <- function(p, file, width = plot.width, height = plot.height, res = plot.res, pdf = save.pdf) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  Cairo::CairoPNG(paste0(file, ".png"), width = width, height = height, res = res, bg = "white")
  tryCatch(print(p), finally = grDevices::dev.off())
  if(pdf) ggsave(paste0(file, ".pdf"), p, device = grDevices::cairo_pdf,
                width = width / res, height = height / res, units = "in", limitsize = FALSE)
}

# 6-2. 통계: 원본의 exact·연속성 보정·최소 표본 수 조건을 유지합니다.
test.score <- function(d, type, unit) {
  n <- nrow(d)
  young <- d$score[d$group == "Young"]
  old <- d$score[d$group == "Old"]
  difference <- (if(length(old)) mean(old) else NA_real_) - (if(length(young)) mean(young) else NA_real_)
  r <- data.frame(test = type, rho = NA_real_, p_value = NA_real_, mean_difference_Old_minus_Young = difference)
  if(!n) return(r)
  if(unit == "mouse" && anyDuplicated(d$mouse.id)) stop("Repeated mouse in a mouse-level test.")

  if(type %in% c("Age_Spearman", "Mff_Spearman")) {
    x <- if(type == "Age_Spearman") d$age_months else d$mff
    if(any(!is.finite(x)) || length(unique(x)) < 2 || length(unique(d$score)) < 2) return(r)
    r$rho <- unname(cor(x, d$score, method = "spearman"))
    if(n < min.observations.correlation) return(r)
    use.exact <- n < 10 && !anyDuplicated(x) && !anyDuplicated(d$score)
    run.test <- function() cor.test(x, d$score, method = "spearman", alternative = "two.sided", exact = use.exact)
  } else if(type == "Age_Kruskal_Wallis") {
    k <- length(unique(d$age))
    if(k < 2 || n <= k || length(unique(d$score)) < 2) return(r)
    run.test <- function() kruskal.test(d$score, factor(d$age))
  } else if(type == "Young_Old_Wilcoxon") {
    minimum <- if(unit == "mouse") min.mice.per.group.wilcoxon else min.cells.per.group.wilcoxon
    if(min(length(young), length(old)) < minimum || length(unique(d$score)) < 2) return(r)
    use.exact <- unit == "mouse" && max(length(young), length(old)) < 50 && !anyDuplicated(d$score)
    run.test <- function() wilcox.test(young, old, alternative = "two.sided", paired = FALSE,
                                     exact = use.exact, correct = !use.exact)
  } else stop("Unknown test: ", type)

  p <- suppressWarnings(tryCatch(run.test()$p.value, error = function(e) NA_real_))
  if(is.finite(p)) r$p_value <- p
  r
}

# 6-3. 연령별 / Young-Old score 그림
age.plot <- function(d, unit, comparison, pvalue, title) {
  levels <- if(comparison == "age") age.order else c("Young", "Old")
  palette <- if(comparison == "age") age.colors else group.colors
  d$x <- factor(d[[comparison]], levels = levels)
  p <- ggplot(d, aes(x, score)) + geom_blank()
  if(unit == "cell") {
    regular <- d %>% group_by(x) %>% filter(n() >= 2, n_distinct(score) > 1) %>% ungroup()
    special <- unique(d[!d$x %in% regular$x, c("x", "score")])
    if(nrow(regular)) p <- p +
      geom_violin(data = regular, aes(fill = x), width = plot.violin.width, trim = TRUE, scale = "width") +
      geom_boxplot(data = regular, width = plot.box.width, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = palette)
    if(nrow(special)) p <- p + geom_point(data = special, size = 2)
  } else {
    boxes <- d %>% group_by(x) %>% filter(n() >= 2) %>% ungroup()
    if(nrow(boxes)) p <- p + geom_boxplot(data = boxes, width = 0.55, fill = "grey93", outlier.shape = NA)
    if(nrow(d)) p <- p + geom_point(aes(color = x), size = mouse.point.size,
                                  position = position_jitter(width = 0.08, height = 0, seed = 1)) +
      scale_color_manual(values = palette)
  }
  bounds <- y.bounds(d$score)
  p + scale_x_discrete(limits = levels, labels = if(comparison == "age") sub("m$", "", levels) else levels) +
    annotate("text", x = length(levels) + 0.35, y = bounds[2], label = format.p(pvalue),
             hjust = 1, vjust = 1, size = plot.pvalue.size) + coord_cartesian(ylim = bounds) +
    labs(title = title, x = if(comparison == "age") "Age (months)" else NULL,
         y = if(unit == "cell") "AddModuleScore" else "Mean AddModuleScore per mouse") + plot.theme()
}

# 6-4. Scatter: 전체(연령색/그룹색), Young만, Old만
scatter.plot <- function(d, s, unit, view, color.by, title) {
  bounds <- y.bounds(d$score, top = 0.65)
  xb <- if(nrow(d)) range(d$mff) else c(0, 1)
  if(diff(xb) == 0) xb <- xb + c(-0.05, 0.05) * max(1, abs(xb[1]))
  if(view != "All") { d <- d[d$group == view, , drop = FALSE]; title <- paste(title, view, sep = " / ") }
  s <- s[s$analysis_group == view, , drop = FALSE]
  levels <- if(color.by == "group") c("Young", "Old") else
    if(view == "All") age.order else if(view == "Young") young.ages else old.ages
  d$color <- factor(d[[color.by]], levels = levels)
  p <- ggplot(d, aes(mff, score)) +
    geom_point(aes(color = color), size = if(unit == "cell") cell.point.size else mouse.point.size,
               alpha = if(unit == "cell") cell.point.alpha else 0.95)
  if(nrow(d)) p <- p +
    scale_color_manual(values = if(color.by == "age") age.colors else group.colors, breaks = levels,
                       labels = if(color.by == "age") sub("m$", "", levels) else levels, drop = FALSE,
                       name = if(color.by == "age") "Age (months)" else "Age group") +
    guides(color = guide_legend(nrow = if(length(levels) > 3) 2 else 1, override.aes = list(alpha = 1, size = 3)))
  if(show.regression.line && nrow(d) >= 2 && length(unique(d$mff)) >= 2) p <- p +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
                color = regression.line.color, linewidth = regression.line.width, show.legend = FALSE)
  label <- paste0("rho = ", if(is.finite(s$rho)) sprintf("%.2f", s$rho) else "NA", "\n", format.p(s$p_value))
  p + annotate("text", x = Inf, y = bounds[2], label = label, hjust = 1.05, vjust = 1, size = plot.pvalue.size) +
    coord_cartesian(xlim = xb, ylim = bounds) +
    labs(title = title, x = if(unit == "cell") "Mff expression (RNA LogNormalize)" else "Mff pseudobulk log2(TMM CPM + 1)",
         y = if(unit == "cell") "AddModuleScore" else "Mean AddModuleScore per mouse") + plot.theme(nrow(d) > 0)
}

# 6-5. 전체 조합에 적용: 조직 → geneset → mouse/cell → 세포형
figure.dir <- file.path(output.dir, "06.making_figures")
statistics.list <- list()
for(tissue.name in tissues) {
  celltypes <- sort(unique(cell.data.list[[tissue.name]]$celltype))
  for(gs in names(score.columns)) {
    for(unit in c("mouse", "cell")) {
      source.data <- if(unit == "mouse") mouse.data.list[[tissue.name]] else cell.data.list[[tissue.name]]
      if(unit == "mouse") source.data <- source.data[source.data$eligible, , drop = FALSE]
      d.all <- as.data.frame(source.data[, c("tissue", "celltype", "mouse.id", "age", "age_months", "group")])
      d.all$score <- source.data[[score.columns[[gs]]]]
      d.all$mff <- if(unit == "mouse") source.data$Mff_log2_TMM_CPM_plus1 else source.data$Mff_LogNormalize
      if(any(!is.finite(d.all$score)) || any(!is.finite(d.all$mff)) || anyNA(d.all[, c("age", "group", "mouse.id")]))
        stop("Missing/nonfinite input: ", tissue.name, " / ", gs, " / ", unit)
      folder <- file.path(figure.dir, tissue.name, gs, unit)
      plots <- setNames(lapply(1:6, function(i) list()), c("age", "young_old", "scatter_age", "scatter_group", "scatter_Young", "scatter_Old"))
      results <- list()

      for(i in seq_along(celltypes)) {
        ct <- celltypes[i]
        d <- d.all[d.all$celltype == ct, , drop = FALSE]
        title <- paste(paste(strwrap(paste(gsub("_", " ", tissue.name), ct, sep = " / "),
                                      width = plot.title.wrap.width), collapse = "\n"), score.labels[[gs]], sep = "\n")
        tests <- c("Age_Kruskal_Wallis", "Young_Old_Wilcoxon", "Age_Spearman", "Mff_Spearman")
        s <- bind_rows(lapply(tests, function(type) test.score(d, type, unit)))
        s$analysis_group <- "All"
        for(view in c("Young", "Old")) {
          z <- test.score(d[d$group == view, , drop = FALSE], "Mff_Spearman", unit)
          z$test <- paste0("Mff_Spearman_", view)
          z$analysis_group <- view
          s <- bind_rows(s, z)
        }
        s$tissue <- tissue.name; s$celltype <- ct; s$geneset <- gs; s$unit <- unit
        results[[ct]] <- s
        plots$age[[ct]] <- age.plot(d, unit, "age", s$p_value[s$test == "Age_Kruskal_Wallis"], title)
        if(make.young.old.plots)
          plots$young_old[[ct]] <- age.plot(d, unit, "group", s$p_value[s$test == "Young_Old_Wilcoxon"], title)
        mff.stats <- s[grepl("^Mff_Spearman", s$test), , drop = FALSE]
        plots$scatter_age[[ct]] <- scatter.plot(d, mff.stats, unit, "All", "age", title)
        plots$scatter_group[[ct]] <- scatter.plot(d, mff.stats, unit, "All", "group", title)
        plots$scatter_Young[[ct]] <- scatter.plot(d, mff.stats, unit, "Young", "age", title)
        plots$scatter_Old[[ct]] <- scatter.plot(d, mff.stats, unit, "Old", "age", title)
        stem <- sprintf("%03d.%s", i, substr(gsub("[^A-Za-z0-9_-]+", "_", ct), 1, 90))
        for(name in names(plots)) if(!is.null(plots[[name]][[ct]]))
          save.plot(plots[[name]][[ct]], file.path(folder, name, stem),
                    width = if(grepl("^scatter", name)) scatter.width else plot.width)
      }

      statistics <- bind_rows(results)
      statistics.list[[paste(tissue.name, gs, unit, sep = ".")]] <- statistics
      if(save.overviews) for(name in names(plots)) if(length(plots[[name]])) {
        nc <- min(3, length(plots[[name]]))
        width <- if(grepl("^scatter", name)) scatter.width else plot.width
        save.plot(wrap_plots(plots[[name]], ncol = nc), file.path(folder, "overview", name),
                  width = width * nc, height = plot.height * ceiling(length(plots[[name]]) / nc), pdf = FALSE)
      }

      # 6-6. Dotplot: 연령별 평균 score / Old-Young 평균 차이 + 기존 p 재사용
      if(make.age.dotplots) {
        summary <- d.all %>% group_by(celltype, age) %>% summarise(mean_score = mean(score), .groups = "drop")
        dot <- expand.grid(celltype = celltypes, age = age.order) %>% left_join(summary, by = c("celltype", "age"))
        dot$celltype <- factor(dot$celltype, levels = rev(celltypes))
        dot$age <- factor(dot$age, levels = age.order)
        kw <- statistics[statistics$test == "Age_Kruskal_Wallis", , drop = FALSE]
        kw$celltype <- factor(kw$celltype, levels = rev(celltypes))
        kw$label <- vapply(kw$p_value, format.p, character(1))
        colors <- if(gs %in% red.when.old.higher) c("#000000", "#F3F3F3", "#D95F02") else c("#D95F02", "#F3F3F3", "#000000")
        limit <- max(abs(dot$mean_score), 0.05, na.rm = TRUE)
        title <- paste(gsub("_", " ", tissue.name), unit, score.labels[[gs]], sep = " / ")
        p.age <- ggplot(dot, aes(age, celltype)) +
          geom_point(data = dot[is.finite(dot$mean_score), ], aes(color = mean_score), size = dotplot.point.size) +
          geom_text(data = dot[!is.finite(dot$mean_score), ], label = "NA", color = "grey50", size = dotplot.na.text.size) +
          geom_text(data = kw, aes(x = "p-value", label = label), size = dotplot.na.text.size) +
          scale_color_gradient2(low = colors[1], mid = colors[2], high = colors[3], midpoint = 0,
                                limits = c(-limit, limit), name = if(unit == "mouse") "Mean score\n(equal mice)" else "Mean score\n(equal cells)") +
          scale_x_discrete(limits = c(age.order, "p-value"), labels = c(sub("m$", "", age.order), "p-value")) +
          scale_y_discrete(drop = FALSE) + labs(title = title, x = "Age (months)", y = NULL) +
          plot.theme.dot(legend = TRUE) + theme(legend.position = "right")

        delta <- statistics[statistics$test == "Young_Old_Wilcoxon", , drop = FALSE]
        delta <- delta[order(is.na(delta$mean_difference_Old_minus_Young), delta$mean_difference_Old_minus_Young, delta$celltype), ]
        delta$celltype <- factor(delta$celltype, levels = rev(delta$celltype))
        delta$difference <- delta$mean_difference_Old_minus_Young
        delta$direction <- ifelse(delta$difference > 0, "Old higher", ifelse(delta$difference < 0, "Old lower", "Equal means"))
        delta$significance <- ifelse(is.na(delta$p_value), "Not tested", ifelse(delta$p_value < dotplot.p.cutoff, "Below cutoff", "At or above cutoff"))
        delta$label <- vapply(delta$p_value, format.p, character(1))
        limit <- max(abs(delta$difference), 0.05, na.rm = TRUE) * 1.25
        p.delta <- ggplot(delta, aes(difference, celltype)) +
          geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
          geom_point(data = delta[is.finite(delta$difference), ], aes(color = direction, shape = significance),
                     size = dotplot.point.size, stroke = dotplot.point.stroke) +
          geom_text(data = delta[!is.finite(delta$difference), ], x = limit * 0.96, label = "NA", size = dotplot.na.text.size, color = "grey50") +
          geom_text(aes(label = label), x = limit * 1.90, hjust = 1, size = dotplot.na.text.size) +
          scale_color_manual(values = c("Old lower" = colors[1], "Old higher" = colors[3], "Equal means" = "grey50"), guide = "none") +
          scale_shape_manual(values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4),
                             limits = c("Below cutoff", "At or above cutoff", "Not tested"), drop = FALSE, name = NULL,
                             labels = c(paste0("p < ", dotplot.p.cutoff), paste0("p >= ", dotplot.p.cutoff), "p = NA")) +
          scale_x_continuous(limits = c(-limit, limit * 1.95), breaks = pretty(c(-limit, limit))) +
          scale_y_discrete(drop = FALSE) + labs(title = title, x = "Mean score difference (Old - Young)", y = NULL) +
          plot.theme.dot(legend = TRUE)
        width <- round(dotplot.width.in * dotplot.res)
        height <- round(max(dotplot.height.minimum.in, 2 + length(celltypes) * dotplot.height.per.celltype.in) * dotplot.res)
        save.plot(p.age, file.path(folder, "dotplot", "age_mean_score"), width, height, dotplot.res)
        save.plot(p.delta, file.path(folder, "dotplot", "Young_Old_score_difference"), width, height, dotplot.res)
        write.csv(dot, file.path(folder, "dotplot", "age_mean_score.csv"), row.names = FALSE)
        write.csv(delta, file.path(folder, "dotplot", "Young_Old_score_difference.csv"), row.names = FALSE)
      }
    }
  }
}

# 6-7. 통계표 저장: BH는 원본처럼 조직 × 분석 단위 × 검정별로 유효한 p만 보정합니다.
all.statistics <- bind_rows(statistics.list) %>% group_by(tissue, unit, test) %>%
  mutate(p_adj_BH = {
    p <- rep(NA_real_, n()); ok <- is.finite(p_value)
    p[ok] <- p.adjust(p_value[ok], method = "BH"); p
  }) %>% ungroup()
write.csv(all.statistics, file.path(figure.dir, "02.all_statistics.csv"), row.names = FALSE)
