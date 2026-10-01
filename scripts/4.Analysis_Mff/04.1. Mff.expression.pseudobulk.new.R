# 조직별 RDS → mouse × celltype raw counts 합산 → TMM → 연령별 그림/검정
# RNA assay에 합쳐진 counts layer가 있다고 가정합니다.
# 한 점 = 한 마우스, 추이 그림의 검은 점 = 평균, 오차막대 = SEM.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

# 1. 설정 --------------------------------------------------------------------
project.dir <- "/BiO/Live/dleogus32/202609Aging_MFF"
tissues <- c("Heart", "Limb_Muscle")
rds.files <- setNames(
  file.path(project.dir, "Analysis/4.tissue_split", tissues,
            paste0("Aging.MFF.", tissues, ".seurat.normalization.pca.umap.RDS")),
  tissues
)
output.base <- file.path(project.dir, "Analysis/5.Mff_expression")
gene <- "Mff"
assay <- "RNA"
min.cells.per.mouse.celltype <- 1
output.folder <- paste0("mouse_pseudobulk_TMM_min", min.cells.per.mouse.celltype)
age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
age.months <- as.numeric(sub("m$", "", age.order))
age.colors <- c(
  "1m" = "#1B9E77", "3m" = "#D95F02", "18m" = "#7570B3",
  "21m" = "#E7298A", "24m" = "#66A61E", "30m" = "#E6AB02"
)
save.overviews <- TRUE
plot.width <- 1500
plot.height <- 1100
plot.res <- 160
plot.base.size <- 13
plot.title.size <- 20
plot.axis.title.x.size <- 13
plot.axis.title.y.size <- 13
plot.axis.text.x.size <- 20
plot.axis.text.y.size <- 20
plot.title.wrap.width <- 48
plot.pvalue.size <- 10

# 2. 공통 함수 ---------------------------------------------------------------
join.values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & trimws(x) != ""])))
  if(length(x)) paste(x, collapse = ";") else NA_character_
}

save.plot <- function(p, file, width = plot.width, height = plot.height) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  Cairo::CairoPNG(paste0(file, ".png"), width = width, height = height,
                  res = plot.res, bg = "white")
  tryCatch(print(p), finally = grDevices::dev.off())
}

plot.theme <- function() {
  theme_bw(base_size = plot.base.size) + theme(
    legend.position = "none",
    plot.title = element_text(size = plot.title.size, face = "bold"),
    axis.title.x = element_text(size = plot.axis.title.x.size),
    axis.title.y = element_text(size = plot.axis.title.y.size),
    axis.text.x = element_text(size = plot.axis.text.x.size),
    axis.text.y = element_text(size = plot.axis.text.y.size),
    panel.grid.minor = element_blank()
  )
}

format.p <- function(p) {
  if(!is.finite(p)) return("p = NA")
  if(p == 0) return(paste0(
    "p < ", formatC(.Machine$double.xmin, format = "e", digits = 2)
  ))
  paste0("p = ", formatC(p, format = if(p < 0.001) "e" else "f",
                        digits = if(p < 0.001) 2 else 3))
}

# 관측 연령이 2개 이상이고, 적어도 한 연령에 마우스가 2마리 이상이어야 검정.
test.kruskal <- function(d, tissue.name, ct) {
  age.n <- table(factor(d$age, levels = age.order))
  k <- sum(age.n > 0)
  r <- data.frame(
    tissue = tissue.name, celltype = ct, gene = gene,
    n_mice = nrow(d), n_age_groups = k,
    observed_ages = paste(age.order[age.n > 0], collapse = ";"),
    Mff_total_count = sum(d$gene_raw_count), H = NA_real_, df_test = NA_real_,
    PValue = NA_real_, FDR_BH_within_tissue = NA_real_,
    status = "not_tested", reason = ""
  )
  for(a in age.order) r[[paste0("n_mice_", a)]] <- as.integer(age.n[a])
  if(anyDuplicated(d$mouse.id)) stop("Repeated mouse.id in test input.")
  if(!nrow(d)) r$reason <- "no_eligible_mice" else
    if(k < 2) r$reason <- "fewer_than_two_observed_age_groups" else
      if(nrow(d) <= k) r$reason <- "no_within_age_replication" else
        if(sum(d$gene_raw_count) == 0) r$reason <- "target_all_zero" else
          if(length(unique(d$log2_CPM_plus1)) < 2) {
            r$reason <- "all_normalized_values_identical"
          }
  if(nzchar(r$reason)) return(r)
  fit <- kruskal.test(d$log2_CPM_plus1, factor(d$age))
  if(is.finite(fit$p.value)) {
    r$H <- unname(fit$statistic)
    r$df_test <- unname(fit$parameter)
    r$PValue <- fit$p.value
    r$status <- "tested"
    r$reason <- "ok"
  } else r$reason <- "nonfinite_p"
  r
}

for(tissue.name in names(rds.files)) {
  # 3. 데이터 준비 -----------------------------------------------------------
  tissue.rds <- readRDS(rds.files[[tissue.name]])
  DefaultAssay(tissue.rds) <- assay
  if(!"counts" %in% SeuratObject::Layers(tissue.rds[[assay]])) {
    stop("RNA counts layer required: ", tissue.name)
  }
  counts <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "counts")
  counts <- counts[, colnames(tissue.rds), drop = FALSE]
  meta <- tissue.rds@meta.data[colnames(counts), , drop = FALSE]
  celltypes <- as.character(meta$cell_ontology_class)
  celltypes[is.na(celltypes) | trimws(celltypes) == ""] <- "Unannotated"
  cells <- data.frame(
    cell = colnames(counts), tissue = tissue.name, celltype = celltypes,
    mouse.id = as.character(meta$mouse.id), age = factor(meta$age, levels = age.order),
    sex = if("sex" %in% names(meta)) as.character(meta$sex) else NA_character_,
    orig.ident = if("orig.ident" %in% names(meta)) as.character(meta$orig.ident) else
      NA_character_
  )
  mouse.ages <- unique(cells[, c("mouse.id", "age")])
  if(anyNA(cells$age) || anyNA(cells$mouse.id) || any(trimws(cells$mouse.id) == "") ||
     anyDuplicated(mouse.ages$mouse.id)) stop("Check mouse.id and age: ", tissue.name)

  # 4. 마우스 × 세포형별 합산 후 세포형별 TMM 정규화 ----------------------------
  mice <- cells %>% group_by(tissue, celltype, mouse.id, age) %>%
    summarise(n_cells = n(), sex = join.values(sex),
              technical_runs = join.values(orig.ident), .groups = "drop") %>%
    arrange(celltype, age, mouse.id)
  mice$pb_id <- sprintf("PB%05d", seq_len(nrow(mice)))
  cells <- left_join(cells, dplyr::select(mice, celltype, mouse.id, pb_id),
                     by = c("celltype", "mouse.id"))
  tissue.rds$pseudobulk_id <- setNames(cells$pb_id, cells$cell)
  pb <- AggregateExpression(
    tissue.rds, assays = assay, features = rownames(counts),
    group.by = "pseudobulk_id", return.seurat = FALSE, verbose = FALSE
  )[[assay]]
  pb <- pb[rownames(counts), mice$pb_id, drop = FALSE]
  mice$gene <- gene
  mice$age_months <- as.numeric(sub("m$", "", as.character(mice$age)))
  mice$gene_raw_count <- as.numeric(pb[gene, ])
  mice$library_size <- as.numeric(Matrix::colSums(pb))
  mice$eligible <- mice$n_cells >= min.cells.per.mouse.celltype & mice$library_size > 0
  mice$norm_factor <- mice$effective_library_size <- mice$normalized_CPM <- NA_real_
  mice$log2_CPM_plus1 <- NA_real_
  celltypes <- sort(unique(mice$celltype))
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
    mice$effective_library_size[idx] <- y$samples$lib.size * y$samples$norm.factors
    mice$normalized_CPM[idx] <- mff.cpm
    mice$log2_CPM_plus1[idx] <- log2(mff.cpm + 1)
  }
  if(any(!is.finite(mice$log2_CPM_plus1[mice$eligible]))) stop("Invalid mouse CPM.")

  # 5. Kruskal-Wallis 및 연령별 평균/SEM --------------------------------------
  statistics <- bind_rows(lapply(celltypes, function(ct) {
    d <- mice[mice$celltype == ct & mice$eligible, ]
    test.kruskal(d, tissue.name, ct)
  }))
  ok <- is.finite(statistics$PValue)
  statistics$FDR_BH_within_tissue[ok] <- p.adjust(statistics$PValue[ok], method = "BH")
  statistics$plot_label <- vapply(statistics$PValue, format.p, character(1))
  age.summary <- mice %>% group_by(celltype, age) %>% summarise(
    n_mice_total = n(), n_mice = sum(eligible), n_cells_total = sum(n_cells),
    n_cells_included = sum(n_cells[eligible]),
    mean_log2_CPM_plus1 = if(any(eligible)) mean(log2_CPM_plus1[eligible]) else NA_real_,
    sd_log2_CPM_plus1 = if(sum(eligible) > 1) sd(log2_CPM_plus1[eligible]) else NA_real_,
    .groups = "drop"
  )
  grid <- expand.grid(celltype = celltypes, age = age.order)
  grid$age <- factor(grid$age, levels = age.order)
  age.summary <- left_join(grid, age.summary, by = c("celltype", "age")) %>%
    arrange(celltype, age)
  for(nm in c("n_mice_total", "n_mice", "n_cells_total", "n_cells_included")) {
    age.summary[[nm]][is.na(age.summary[[nm]])] <- 0
  }
  age.summary$tissue <- tissue.name
  age.summary$age_months <- as.numeric(sub("m$", "", as.character(age.summary$age)))
  age.summary$sem_log2_CPM_plus1 <- with(age.summary,
    ifelse(n_mice > 1, sd_log2_CPM_plus1 / sqrt(n_mice), NA_real_)
  )
  file.stems <- sprintf("%03d.%s", seq_along(celltypes),
                        substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100))
  output.dir <- file.path(output.base, tissue.name, output.folder)
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(pb, file.path(output.dir, "02.pseudobulk_raw_counts.RDS"), compress = "gzip")

  # 6. Boxplot + 마우스 점 / 실제 개월 수에 따른 평균±SEM ----------------------
  box.plots <- trend.plots <- list()
  for(i in seq_along(celltypes)) {
    ct <- celltypes[i]
    d <- mice[mice$celltype == ct & mice$eligible, ]
    s <- age.summary[age.summary$celltype == ct, ]
    title <- paste(strwrap(paste(tissue.name, ct, sep = " | "),
                           width = plot.title.wrap.width), collapse = "\n")
    common <- list(
      labs(title = title, x = "Age (months)",
           y = paste0(gene, " pseudobulk log2(normalized CPM + 1)")), plot.theme()
    )
    if(nrow(d)) common <- c(common, list(scale_color_manual(values = age.colors,
                                                           drop = FALSE)))
    ymax <- max(1, d$log2_CPM_plus1 * 1.12)
    age.labels <- setNames(paste0(s$age, "\nn=", s$n_mice), as.character(s$age))
    boxes <- d %>% group_by(age) %>% filter(n() >= 2) %>% ungroup()
    p.box <- ggplot(d, aes(age, log2_CPM_plus1)) + geom_blank()
    if(nrow(boxes)) p.box <- p.box + geom_boxplot(
      data = boxes, width = 0.55, outlier.shape = NA, fill = "grey92"
    )
    if(nrow(d)) p.box <- p.box + geom_point(
      aes(color = age), size = 3,
      position = position_jitter(width = 0.08, height = 0, seed = 1234)
    )
    p.box <- p.box + common +
      scale_x_discrete(limits = age.order, labels = age.labels, drop = FALSE) +
      coord_cartesian(ylim = c(0, ymax * 1.25)) +
      annotate("text", x = Inf, y = Inf, label = statistics$plot_label[i],
               hjust = 1.03, vjust = 1.15, size = plot.pvalue.size, lineheight = 1.05)

    p.trend <- ggplot(d, aes(age_months, log2_CPM_plus1)) + geom_blank()
    if(sum(s$n_mice > 0) >= 2) p.trend <- p.trend + geom_line(
      data = s, aes(age_months, mean_log2_CPM_plus1, group = 1),
      inherit.aes = FALSE, color = "grey45", na.rm = TRUE
    )
    if(nrow(d)) p.trend <- p.trend + geom_point(
      aes(color = age), size = 2.8,
      position = position_jitter(width = 0.15, height = 0, seed = 1234)
    )
    p.trend <- p.trend + geom_point(
      data = s[s$n_mice > 0, ], aes(age_months, mean_log2_CPM_plus1),
      inherit.aes = FALSE, shape = 18, size = 4, color = "black"
    )
    sem.data <- s[s$n_mice > 1, ]
    if(nrow(sem.data)) p.trend <- p.trend + geom_errorbar(
      data = sem.data, aes(age_months, ymin = mean_log2_CPM_plus1 - sem_log2_CPM_plus1,
                          ymax = mean_log2_CPM_plus1 + sem_log2_CPM_plus1),
      inherit.aes = FALSE, width = 0.4, color = "black"
    )
    p.trend <- p.trend + common + coord_cartesian(ylim = c(0, ymax)) +
      scale_x_continuous(breaks = age.months, labels = age.order, limits = c(0, 31))
    stem <- file.stems[i]
    save.plot(p.box, file.path(output.dir, "boxplot", paste0(stem, ".mouse_boxplot")))
    save.plot(p.trend, file.path(output.dir, "age_trend", paste0(stem, ".mouse_age_trend")))
    box.plots[[i]] <- p.box
    trend.plots[[i]] <- p.trend
  }
  if(save.overviews) {
    plots <- list(mouse_boxplot = box.plots, mouse_age_trend = trend.plots)
    pages <- split(seq_along(celltypes), ceiling(seq_along(celltypes) / 9))
    for(kind in names(plots)) for(page in seq_along(pages)) {
      idx <- pages[[page]]
      nc <- min(3, length(idx))
      stem <- sprintf("%s.all_celltypes.%s.page%02d", gene, kind, page)
      save.plot(
        wrap_plots(plots[[kind]][idx], ncol = nc), file.path(output.dir, "overview", stem),
        width = 1200 * nc, height = 1000 * ceiling(length(idx) / nc)
      )
    }
  }
  message(tissue.name, ": pseudobulk RDS / figures saved")
  rm(tissue.rds, counts, pb, cells, mice, box.plots, trend.plots)
  invisible(gc())
}
