
#=====================================================================

local({
  required <- c("Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "Cairo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing) > 0){stop("Missing R packages: ", paste(missing, collapse = ", "))}
  if(packageVersion("SeuratObject") < "5.0.0"){stop("SeuratObject >= 5.0.0 is required.")}
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))

  #-------------------------------------------------------------------
  # 1. Settings

  gene <- "Mff"
  age.order <- c("1m", "3m", "24m", "30m")
  age.colors <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"), age.order)
  rds.files <- c(Heart = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Heart/Aging.MFF.Heart.seurat.normalization.pca.umap.RDS", Limb_Muscle = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Limb_Muscle/Aging.MFF.Limb_Muscle.seurat.normalization.pca.umap.RDS")
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression"

  # 세포형별로 관측된 모든 연령군을 한 번에 비교: p-value 1개
  pvalue.column <- "p_value"   # "p_value" = raw p, "p_adj_BH" = 조직 내 세포형 간 BH 보정값
  plot.pvalue.text.size <- 10   # p-value 글자만 조절 (ggplot2의 mm 단위)
  if(!pvalue.column %in% c("p_value", "p_adj_BH")){stop("pvalue.column must be p_value or p_adj_BH.")}

  # 그림 글씨 크기 (pt): 필요하면 아래 값만 변경한 뒤 전체 코드를 다시 실행
  # 부제목과 오른쪽 아래 각주는 표시하지 않습니다.
  # 기존 output directory와 PNG 파일명을 그대로 사용하여 덮어씁니다.
  plot.base.size <- 13
  plot.title.size <- 20          # 그림 제목
  plot.axis.title.x.size <- 13   # X축 제목: Age (months)
  plot.axis.title.y.size <- 13   # Y축 제목: Mff expression (...)
  plot.axis.text.x.size <- 30    # X축 눈금: age만 표시
  plot.axis.text.y.size <- 30    # Y축 눈금 숫자
  plot.title.wrap.width <- 48    # 제목 줄바꿈 기준 글자 수
  plot.violin.width <- 0.9       # 나이 한 칸 기준 최대 바이올린 폭 (빈 나이가 있어도 고정)
  plot.violin.box.width <- 0.18  # 바이올린 위에 겹치는 박스 폭
  plot.violin.box.alpha <- 0.85  # 흰색 박스 불투명도 (0~1)
  plot.violin.box.linewidth <- 0.45  # 박스 테두리와 수염 두께

  if(any(!file.exists(rds.files))){stop("Missing input RDS: ", paste(rds.files[!file.exists(rds.files)], collapse = ", "))}

  save.png <- function(plot, filename, width = 1500, height = 1100) {
    Cairo::CairoPNG(filename = filename, width = width, height = height, res = 160, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
  }

  #-------------------------------------------------------------------
  # One omnibus test uses all cell-level values and all observed age groups.
  # Empty age groups are not filled with zeros; zero-expression CELLS remain.

  make.pvalue.table <- function(d, tissue.name, ct) {
    counts <- table(factor(as.character(d$age), levels = age.order))
    observed.ages <- age.order[counts > 0L]
    missing.ages <- setdiff(age.order, observed.ages)
    n.ages <- length(observed.ages)
    n.mice.by.age <- vapply(age.order, function(a) length(unique(d$mouse.id[as.character(d$age) == a])), integer(1))
    r <- data.frame(tissue = tissue.name, celltype = ct, gene = gene, n_cells = nrow(d), n_mice = length(unique(d$mouse.id)), n_age_groups = n.ages, observed_ages = paste(observed.ages, collapse = ";"), missing_ages = paste(missing.ages, collapse = ";"), all_six_ages_present = n.ages == length(age.order), n_cells_by_age = paste(paste0(age.order, "=", as.integer(counts)), collapse = ";"), n_mice_by_age = paste(paste0(age.order, "=", n.mice.by.age), collapse = ";"), min_cells_per_observed_age = if(n.ages) min(counts[counts > 0L]) else NA_integer_, H = NA_real_, df = if(n.ages >= 2L) n.ages - 1L else NA_integer_, p_value = NA_real_, test = "Kruskal-Wallis rank sum", comparison = "All observed age groups jointly", test_unit = "cell", pvalue_method = "Tie-corrected H; asymptotic chi-squared", mouse_dependence_adjusted = FALSE, status = "not_tested", reason = "", warnings = "", stringsAsFactors = FALSE)
    if(n.ages < 2L){r$status <- "insufficient_age_groups"; r$reason <- "Fewer than two observed age groups."; return(r)}
    if(length(unique(d$expression)) == 1L){r$status <- "all_values_identical"; r$reason <- "All pooled expression values are identical; the tie correction is zero."; return(r)}
    test.warnings <- character()
    result <- tryCatch(withCallingHandlers(stats::kruskal.test(x = d$expression, g = factor(as.character(d$age), levels = observed.ages)), warning = function(w) {test.warnings <<- c(test.warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
    r$warnings <- paste(unique(test.warnings), collapse = " | ")
    if(inherits(result, "error")){r$status <- "test_error"; r$reason <- conditionMessage(result); return(r)}
    if(length(result$p.value) != 1L || !is.finite(result$p.value)){r$status <- "nonfinite_pvalue"; r$reason <- "Kruskal-Wallis returned a nonfinite p-value."; return(r)}
    r$H <- unname(result$statistic)
    r$df <- unname(result$parameter)
    r$p_value <- result$p.value
    r$status <- "ok"
    r
  }

  format.pvalue <- function(p) {
    if(!is.finite(p)){return("NA")}
    if(p == 0){return(paste0("< ", formatC(.Machine$double.xmin, format = "e", digits = 2)))}
    if(p < 0.001){return(formatC(p, format = "e", digits = 2))}
    formatC(p, format = "f", digits = 3)
  }

  add.pvalue.annotation <- function(plot, test.row, annotation.y) {
    prefix <- if(pvalue.column == "p_value") "p" else "BH-adjusted p"
    p <- test.row[[pvalue.column]][[1]]
    p.text <- if(is.finite(p) && p == 0) paste(prefix, format.pvalue(p)) else paste(prefix, "=", format.pvalue(p))
    age.note <- if(test.row$n_age_groups[[1]] < length(age.order)) paste0(" (", test.row$n_age_groups[[1]], " observed ages)") else ""
    label <- p.text
    plot + annotate("text", x = length(age.order) + 0.35, y = annotation.y, label = label, size = plot.pvalue.text.size, hjust = 1, vjust = 1, lineheight = 1.05, colour = "black")
  }

  #-------------------------------------------------------------------
  # 2. Analyze the two tissues separately

  for(tissue.name in names(rds.files)) {
    cat("\nTissue:", tissue.name, "\n")
    x <- readRDS(rds.files[[tissue.name]])
    if(!inherits(x, "Seurat")){stop("Input is not a Seurat object: ", tissue.name)}
    if(!"RNA" %in% names(x@assays)){stop("RNA assay is missing: ", tissue.name)}
    DefaultAssay(x) <- "RNA"
    if(!gene %in% rownames(x[["RNA"]])){stop("Gene is missing: ", gene, " / ", tissue.name)}

    needed <- c("age", "mouse.id", "cell_ontology_class")
    if(!all(needed %in% colnames(x@meta.data))){stop("Missing metadata columns: ", paste(setdiff(needed, colnames(x@meta.data)), collapse = ", "))}
    if(anyNA(x$age) || any(!as.character(x$age) %in% age.order)){stop("Missing or unexpected age in ", tissue.name, "; check age.order.")}
    if(anyNA(x$mouse.id) || any(trimws(as.character(x$mouse.id)) == "")){stop("Missing mouse.id in ", tissue.name)}
    mouse.age <- unique(x@meta.data[, c("mouse.id", "age"), drop = FALSE])
    if(anyDuplicated(as.character(mouse.age$mouse.id))){stop("One mouse.id is linked to multiple ages: ", tissue.name)}

    # Require normalized data explicitly; never fall back to raw counts.
    # If data layers are split by sample, join them in memory only.
    data.layers <- SeuratObject::Layers(x[["RNA"]], search = "^data($|\\.)")
    if(length(data.layers) == 0){stop("RNA normalized data layer is missing: ", tissue.name)}
    if(length(data.layers) != 1 || data.layers != "data") {
      if(!inherits(x[["RNA"]], "Assay5")){stop("Unexpected RNA data layers: ", tissue.name)}
      x <- SeuratObject::JoinLayers(x, assay = "RNA", layers = "^data($|\\.)", new = "data")
    }
    x$plot_age <- factor(as.character(x$age), levels = age.order)
    expr <- FetchData(x, vars = gene, layer = "data", clean = FALSE)
    if(nrow(expr) != ncol(x) || !setequal(rownames(expr), colnames(x))){stop("Expression retrieval did not cover every cell: ", tissue.name)}
    if(any(!is.finite(expr[[gene]])) || any(expr[[gene]] < 0)){stop("Expected finite, nonnegative RNA LogNormalize values: ", tissue.name)}
    meta <- x@meta.data[rownames(expr), , drop = FALSE]
    celltype <- as.character(meta$cell_ontology_class)
    celltype[is.na(celltype) | trimws(celltype) == ""] <- "Unannotated"
    cell.data <- data.frame(cell = rownames(expr), tissue = tissue.name, celltype = celltype, age = factor(as.character(meta$age), levels = age.order), mouse.id = as.character(meta$mouse.id), expression = expr[[gene]], stringsAsFactors = FALSE)
    if(nrow(cell.data) == 0){stop("No cells in input: ", tissue.name)}

    output.dir <- file.path(output.base, tissue.name, "cell_level_box_violin")
    for(folder in c("boxplot", "violin", "overview")){dir.create(file.path(output.dir, folder), recursive = TRUE, showWarnings = FALSE)}
    write.csv(cell.data, file.path(output.dir, "01.cell_expression.csv"), row.names = FALSE)

    # Summary is for reporting. Both plots use individual cell values.
    summary.data <- cell.data %>% group_by(celltype, age) %>% summarise(n_cells = n(), n_mice = n_distinct(mouse.id), mean_expression = mean(expression), median_expression = median(expression), q1 = unname(quantile(expression, 0.25)), q3 = unname(quantile(expression, 0.75)), min_expression = min(expression), max_expression = max(expression), percent_detected = mean(expression > 0) * 100, .groups = "drop")
    celltypes <- sort(unique(cell.data$celltype))
    grid <- expand.grid(celltype = celltypes, age = age.order, stringsAsFactors = FALSE)
    grid$age <- factor(grid$age, levels = age.order)
    summary.data <- left_join(grid, summary.data, by = c("celltype", "age")) %>% arrange(celltype, age)
    summary.data$n_cells[is.na(summary.data$n_cells)] <- 0L
    summary.data$n_mice[is.na(summary.data$n_mice)] <- 0L
    summary.data$tissue <- tissue.name
    summary.data$display <- ifelse(summary.data$n_cells == 0, "no_cells", ifelse(summary.data$n_cells == 1, "single_cell_point", ifelse(summary.data$min_expression == summary.data$max_expression, "constant_value_line", "distribution")))
    write.csv(summary.data, file.path(output.dir, "02.celltype_age_summary.csv"), row.names = FALSE, na = "NA")

    file.map <- data.frame(celltype = celltypes, file_stem = sprintf("%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100)), stringsAsFactors = FALSE)
    write.csv(file.map, file.path(output.dir, "03.celltype_file_map.csv"), row.names = FALSE)
    # One p-value per tissue x cell type, shared by boxplot and violin.
    # BH is applied across successfully tested cell types WITHIN this tissue.
    pvalue.data <- do.call(rbind, lapply(celltypes, function(ct) make.pvalue.table(cell.data[cell.data$celltype == ct, , drop = FALSE], tissue.name, ct)))
    pvalue.data$p_adj_BH <- NA_real_
    tested <- pvalue.data$status == "ok" & is.finite(pvalue.data$p_value)
    pvalue.data$p_adj_BH[tested] <- stats::p.adjust(pvalue.data$p_value[tested], method = "BH")
    pvalue.data$BH_scope <- "All successful cell-type Mff omnibus tests within tissue"
    pvalue.data$BH_n_tests <- sum(tested)
    write.csv(pvalue.data, file.path(output.dir, paste0("05.", gene, ".cellwise_KruskalWallis_omnibus.csv")), row.names = FALSE, na = "NA")
    cat("  Cellwise Kruskal-Wallis:", sum(tested), "tested cell types;", sum(!tested), "NA cell types.\n")
    print(as.data.frame(table(pvalue.data$status)), row.names = FALSE)
    settings <- c(paste0("Gene: ", gene), paste0("Tissue: ", tissue.name), "Input: existing RNA LogNormalize data, including zeros; no normalization repeated", "Plots: individual cells; violin and internal boxplot use the same values", "Test: Kruskal-Wallis rank sum (stats::kruskal.test)", "Question: all observed age groups jointly; one p-value per tissue x cell type", "No pairwise comparisons or regression-slope tests are performed", "P-value calculation: tie-corrected H statistic; asymptotic chi-squared, df=k-1", "Unit: cell; no adjustment for dependence among cells from the same mouse", "Interpretation: these cell-level tests do not establish a mouse-level age effect", "No minimum-cell cutoff; singleton age groups remain in the test", "The chi-squared approximation may be inaccurate with very small groups", "Missing ages are omitted; observed/missing ages and group counts are saved in the test CSV", "Plots with fewer than six observed ages report the observed age-group count", paste0("Plot p-value column: ", pvalue.column), "BH scope: all successful cell-type Mff omnibus tests within this tissue", paste0("BH number of tests: ", sum(tested)), "Fewer than two ages / identical pooled values / failed tests: p=NA", "All cell types are displayed, regardless of statistical significance", "Numerical p-values of zero are displayed as below the smallest positive normal double", "p-values and BH values are shared between boxplot and violin; no duplicate tests", "X-axis: age labels only; cell and mouse counts are retained in the CSVs", "If an older Wilcoxon CSV exists in this folder, it belongs to the previous pairwise analysis; this run writes the KruskalWallis_omnibus CSV")
    writeLines(settings, file.path(output.dir, "06.analysis_settings.txt"))
    box.plots <- list()
    violin.plots <- list()

    #-----------------------------------------------------------------
    # 3. One boxplot and one Seurat violin + boxplot for EVERY cell type

    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      d <- cell.data[cell.data$celltype == ct, , drop = FALSE]
      s <- summary.data[summary.data$celltype == ct, , drop = FALSE]
      age.labels <- setNames(as.character(s$age), as.character(s$age))
      data.top <- max(1, max(d$expression))
      test.row <- pvalue.data[pvalue.data$celltype == ct, , drop = FALSE]
      stopifnot(nrow(test.row) == 1L)
      annotation.y <- data.top * 1.25
      ymax <- data.top * 1.30
      single <- s[s$display == "single_cell_point", , drop = FALSE]
      constant <- s[s$display == "constant_value_line", , drop = FALSE]

      # Boxplot uses all cells, including zeros and outliers.
      p.box <- ggplot(d, aes(x = age, y = expression, fill = age)) + geom_boxplot(width = 0.6, outlier.size = 0.7, outlier.alpha = 0.4) + scale_fill_manual(values = age.colors, drop = FALSE)

      # Density needs >= 2 cells with varying values. Handle all other
      # groups explicitly, so singleton/constant groups are still shown.
      density.ages <- as.character(s$age[s$display == "distribution"])
      if(length(density.ages) > 0) {
        density.cells <- d$cell[as.character(d$age) %in% density.ages]
        x.plot <- subset(x, cells = density.cells)
        x.plot$plot_age <- factor(as.character(x.plot$plot_age), levels = age.order[age.order %in% density.ages])
        p.violin <- VlnPlot(x.plot, features = gene, assay = "RNA", layer = "data", group.by = "plot_age", cols = age.colors[levels(x.plot$plot_age)], pt.size = 0, add.noise = FALSE, y.max = ymax, combine = FALSE)[[1]]
        # 빈 나이 칸이 있을 때 자동 폭이 커지는 것을 방지합니다.
        # Seurat 바이올린 layer의 가로 폭만 고정하며 밀도 계산 설정은 유지합니다.
        for(layer.i in seq_along(p.violin$layers)) {
          if(inherits(p.violin$layers[[layer.i]]$geom, "GeomViolin")) {
            p.violin$layers[[layer.i]]$stat_params$width <- plot.violin.width
            p.violin$layers[[layer.i]]$geom_params$width <- plot.violin.width
          }
        }
        rm(x.plot)
      } else {
        p.violin <- ggplot(d, aes(x = age, y = expression))
      }
      p.violin <- p.violin + geom_blank(data = d, aes(x = age, y = expression), inherit.aes = FALSE)
      # Seurat의 기본 x/y 매핑 대신 원본 cell-level 값을 명시한다.
      # 0 발현 세포도 박스 계산에 포함하고, outlier 점만 숨긴다.
      # 값이 일정한 그룹과 세포 1개 그룹은 아래의 선/점 표시를 사용한다.
      if(length(density.ages) > 0){p.violin <- p.violin + geom_boxplot(data = d[as.character(d$age) %in% density.ages, , drop = FALSE], mapping = aes(x = age, y = expression, group = age), inherit.aes = FALSE, position = "identity", width = plot.violin.box.width, fill = "white", colour = "black", alpha = plot.violin.box.alpha, linewidth = plot.violin.box.linewidth, outlier.shape = NA)}
      if(nrow(constant) > 0){p.violin <- p.violin + geom_errorbar(data = constant, aes(x = age, ymin = min_expression, ymax = max_expression), inherit.aes = FALSE, width = 0.55, colour = "grey20")}
      if(nrow(single) > 0) {
        single.layer <- geom_point(data = single, aes(x = age, y = mean_expression), inherit.aes = FALSE, size = 2, colour = "black")
        p.box <- p.box + single.layer
        p.violin <- p.violin + single.layer
      }

      common <- list(scale_x_discrete(limits = age.order, labels = age.labels, drop = FALSE), coord_cartesian(ylim = c(0, ymax)), labs(x = "Age (months)", y = paste0(gene, " expression (RNA LogNormalize)"), title = paste(strwrap(paste0(tissue.name, " | ", ct), width = plot.title.wrap.width), collapse = "\n"), subtitle = NULL, caption = NULL), theme_bw(base_size = plot.base.size), theme(legend.position = "none", plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.title.x = element_text(size = plot.axis.title.x.size), axis.title.y = element_text(size = plot.axis.title.y.size), axis.text.x = element_text(size = plot.axis.text.x.size), axis.text.y = element_text(size = plot.axis.text.y.size), panel.grid.minor = element_blank()))
      p.box <- p.box + common
      p.violin <- p.violin + common
      p.box <- add.pvalue.annotation(p.box, test.row, annotation.y)
      p.violin <- add.pvalue.annotation(p.violin, test.row, annotation.y)
      stem <- file.map$file_stem[[i]]
      save.png(p.box, file.path(output.dir, "boxplot", paste0(stem, ".boxplot.png")))
      save.png(p.violin, file.path(output.dir, "violin", paste0(stem, ".violin.png")))
      box.plots[[i]] <- p.box
      violin.plots[[i]] <- p.violin
      cat(sprintf("  [%d/%d] %s: %d cells\n", i, length(celltypes), ct, nrow(d)))
    }

    #-----------------------------------------------------------------
    # 4. Tissue overview pages (9 cell types per page)
    # Each cell type has its own Y range; box/violin share that range.

    pages <- split(seq_along(celltypes), ceiling(seq_along(celltypes) / 9))
    for(page in seq_along(pages)) {
      idx <- pages[[page]]
      ncols <- min(3L, length(idx))
      nrows <- ceiling(length(idx) / ncols)
      box.page <- patchwork::wrap_plots(box.plots[idx], ncol = ncols)
      violin.page <- patchwork::wrap_plots(violin.plots[idx], ncol = ncols)
      save.png(box.page, file.path(output.dir, "overview", sprintf("%s.all_celltypes.boxplot.page%02d.png", gene, page)), width = 1200 * ncols, height = 1000 * nrows)
      save.png(violin.page, file.path(output.dir, "overview", sprintf("%s.all_celltypes.violin.page%02d.png", gene, page)), width = 1200 * ncols, height = 1000 * nrows)
    }
    capture.output(sessionInfo(), file = file.path(output.dir, "04.sessionInfo.txt"))
    cat("Saved:", output.dir, "\n")
    rm(x, expr, meta, cell.data, box.plots, violin.plots, p.box, p.violin, box.page, violin.page)
    invisible(gc())
  }
  cat("\nDone: both tissues, all annotated cell types (plus Unannotated if present).\n")
})
