# Aging-MFF: five AddModuleScore signatures, four analysis families.
# Run this WHOLE file with source() in R / RStudio, or Rscript in the shell.
# Input RDS is read only; all new outputs go to Analysis/9.genescore.132430/.
# RNA data must already contain LogNormalize values. No normalization/PCA rerun.
# Source references: https://satijalab.org/seurat/reference/addmodulescore
# https://satijalab.org/seurat/reference/aggregateexpression
# Figure style matches the supplied Limb_Muscle / SaulSenMayo reference.
# Four Mff overviews: All colored by age, All colored by Young/Old, Young only, Old only.
# Each view uses one plotting area, its matching Spearman rho/p, and an optional OLS line.
# show.regression.line controls the line in all four views (default TRUE).
# Main figures: 900 x 850 px at 160 dpi; dotplots: 10.5 inches wide at 300 dpi.
# Age dotplots: signed mean score by age, and Old-minus-Young mean score difference.
# Cellwise tests are exploratory: cells from one mouse are not independent mice.
# Age/Mff associations are unadjusted for age, sex, batch or RNA detection depth.

run.aging.mff.genescores <- function() {
  # ==================================================================
  # 1. SETTINGS: paths, cohort, scoring, plotting
  # ==================================================================
  rds.file <-
    "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"
  geneset.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Workspace/0.Meta/geneset"
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/9.genescore.132430"
  gene <- "Mff"
  assay <- "RNA"
  tissues <- c("Heart", "Limb_Muscle")
  tissue.column <- "tissue_free_annotation"
  celltype.column <- "cell_ontology_class"
  mouse.column <- "mouse.id"
  age.column <- "age"
  include.heart.and.aorta <- TRUE  # Matches the latest supplied Mff script.
  age.order <- c("1m", "3m", "24m", "30m")
  young.ages <- c("1m", "3m")
  old.ages <- c("24m", "30m")
  excluded.ages <- c("18m", "21m")
  make.young.old.plots <- TRUE
  make.age.dotplots <- TRUE

  # Score once per tissue, across ALL selected ages/cell types together.
  # Do not score Young/Old or individual mice separately.
  score.seed <- 1L
  score.nbin <- 24L
  score.ctrl <- 100L
  min.matched.genes <- 1L
  low.coverage.warning <- 0.50
  min.cells.per.mouse.celltype <- 1L  # No new cell-count cutoff by default.
  min.mice.per.group.wilcoxon <- 2L  # Latest supplied script's test-only rule.
  min.cells.per.group.wilcoxon <- 2L
  min.observations.correlation <- 3L  # Computational minimum, not reliability guarantee.
  save.scored.seurat <- TRUE  # Writes a NEW tissue-specific RDS; never the input file.
  save.pdf <- FALSE  # PNG default; TRUE also exports individual vector PDFs.
  save.overviews <- TRUE

  # Plot appearance settings; long function calls are wrapped for readability.
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
  # Dotplot sizing is separate in the reference script (inches, at 300 dpi).
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
  dotplot.p.cutoff <- 0.05  # Raw Wilcoxon p, not an FDR threshold.
  # Display convention: positive scores / Old-minus-Young differences are red,
  # except CoreScence_down, which uses the reversed color direction.
  red.when.old.higher <- c("Hallmark_IFN_alpha", "GO_IFN_beta", "SaulSenMayo", "CoreScence_up")
  # Every cell is plotted and used in statistics; no Mff-positive-only selection.
  age.colors <- c(
    "1m" = "#1B9E77", "3m" = "#D95F02", "24m" = "#7570B3", "30m" = "#E6AB02"
  )
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
  expected.sizes <- c(
    Hallmark_IFN_alpha = 94L, GO_IFN_beta = 77L, SaulSenMayo = 117L, CoreScence_up = 22L,
    CoreScence_down = 16L
  )
  score.columns <- setNames(paste0("AMS_", names(gene.files)), names(gene.files))

  required <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "dplyr", "ggplot2", "patchwork", "Cairo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing)){stop("Missing R packages: ", paste(missing, collapse = ", "))}
  if(packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0") {
    stop("Seurat and SeuratObject >= 5.0.0 required.")
  }
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(!file.exists(rds.file)){stop("Missing input RDS: ", rds.file)}
  if(!setequal(c(young.ages, old.ages), age.order) || length(intersect(young.ages, old.ages))) {
    stop("Young/Old must partition age.order.")
  }
  positive.integers <- c(
    score.nbin, score.ctrl, min.matched.genes, min.cells.per.mouse.celltype, min.mice.per.group.wilcoxon,
    min.cells.per.group.wilcoxon, min.observations.correlation
  )
  if(
    any(!is.finite(positive.integers)) ||
    any(positive.integers < 1 | positive.integers != floor(positive.integers)) ||
    min.observations.correlation < 3
  ) {
    stop("Invalid integer settings.")
  }
  dir.create(output.base, recursive = TRUE, showWarnings = FALSE)

  # ==================================================================
  # 2. Helper functions
  # ==================================================================
  write.csv.safe <- function(x, file) {
    utils::write.csv(x, file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  collapse.values <- function(x) {
    paste(sort(unique(as.character(x[!is.na(x) & trimws(as.character(x)) != ""]))), collapse = ";")
  }
  tissue.key <- function(x) {gsub("[[:space:]-]+", "_", tolower(trimws(as.character(x))))}
  read.grp <- function(file) {
    if(!file.exists(file)){stop("Missing GRP: ", file)}
    x <- readLines(file, warn = FALSE, encoding = "UTF-8")
    x <- trimws(sub("^\ufeff", "", x))
    x <- x[nzchar(x) & !grepl("^#", x)]
    # MSigDB GRP starts with the set name, then a # URL comment, then genes.
    # CoreScence files contain genes only: preserve their first gene.
    msigdb.titles <- c(
      "HALLMARK_INTERFERON_ALPHA_RESPONSE", "GOBP_RESPONSE_TO_INTERFERON_BETA", "SAUL_SEN_MAYO"
    )
    if(length(x) > 0L && x[1L] %in% msigdb.titles) {
      expected.title <- sub("\\.v[0-9].*$", "", tools::file_path_sans_ext(basename(file)))
      if(x[1L] != expected.title) {
        stop("MSigDB GRP title does not match its filename: ", file, " / title=", x[1L])
      }
      message("Removed MSigDB title row: ", x[1L])
      x <- x[-1L]
    }
    if(!length(x) || any(grepl("[[:space:],;]", x))) {
      stop("GRP must contain one mouse gene symbol per line after the optional MSigDB title/comments: ", file)
    }
    if(anyDuplicated(x)){warning("Duplicate GRP entries removed: ", basename(file), call. = FALSE)}
    unique(x)
  }

 ####여기까지####


  join.layer <- function(object, layer.name) {
    pattern <- paste0("^", layer.name, "($|\\.)")
    layers <- SeuratObject::Layers(object[[assay]], search = pattern)
    if(!length(layers)){stop("Missing RNA ", layer.name, " layer. Use the normalized input RDS.")}
    cells <- unlist(
      lapply(layers, function(x) colnames(SeuratObject::LayerData(object[[assay]], layer = x))),
      use.names = FALSE
    )
    if(anyDuplicated(cells) || !setequal(cells, colnames(object))) {
      stop("Duplicated/missing cells in RNA ", layer.name, " layers.")
    }
    if(length(layers) != 1L || layers != layer.name) {
      object <- SeuratObject::JoinLayers(object, assay = assay, layers = pattern, new = layer.name)
    }
    object
  }
  format.p <- function(p) {
    if(!is.finite(p)){return("p = NA")}
    if(p == 0){return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 1)))}
    paste0(
      "p = ", if(p < 0.001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3)
    )
  }
  safe.test <- function(d, type, unit) {
    # Common output preserves estimates, unavailable-test reasons and sample counts.
    n.obs <- nrow(d)
    r <- data.frame(
      test = type, unit = unit, n_observations = n.obs, n_mice = length(unique(d$mouse.id)),
      n_age_groups = length(unique(d$age)), n_Young = sum(d$group == "Young"), n_Old = sum(d$group == "Old"),
      mean_Young = if(any(d$group == "Young")) mean(d$score[d$group == "Young"]) else NA_real_,
      mean_Old = if(any(d$group == "Old")) mean(d$score[d$group == "Old"]) else NA_real_,
      mean_difference_Old_minus_Young = NA_real_, rho = NA_real_, statistic = NA_real_, p_value = NA_real_,
      exact_used = NA, status = "not_tested", reason = "", warning = "", cell_dependence_adjusted = FALSE,
      age_adjusted = FALSE, stringsAsFactors = FALSE
    )
    r$mean_difference_Old_minus_Young <- r$mean_Old - r$mean_Young
    r$interpretation <-
      if(unit == "cell") "Exploratory: same-mouse dependence is not modeled; not mouse-level inference" else
        "One observation per mouse within tissue/cell type; unadjusted association"
    fail <- function(reason) {
      r$reason <- reason
      r
    }
    if(!n.obs){return(fail("No observations"))}
    if(any(!is.finite(d$score)) || anyNA(d$group) || anyNA(d$age)) {
      return(fail("Missing/nonfinite score or age/group"))
    }
    if(unit == "mouse" && anyDuplicated(d$mouse.id)){stop("Repeated mouse in a mouse-level test.")}
    if(type %in% c("Age_Spearman", "Mff_Spearman")) {
      x <- if(type == "Age_Spearman") d$age_months else d$mff
      if(any(!is.finite(x))){return(fail("Missing/nonfinite x values"))}
      if(length(unique(x)) < 2L || length(unique(d$score)) < 2L) {
        return(fail("Constant x or score; correlation undefined"))
      }
      r$rho <- unname(stats::cor(x, d$score, method = "spearman"))
      if(n.obs < min.observations.correlation) {
        return(fail("Too few observations for correlation test; rho retained"))
      }
      use.exact <- n.obs < 10L && !anyDuplicated(x) && !anyDuplicated(d$score)
      r$exact_used <- use.exact
      test.call <- function() stats::cor.test(
        x, d$score, method = "spearman", alternative = "two.sided", exact = use.exact
      )
    } else if(type == "Age_Kruskal_Wallis") {
      if(length(unique(d$age)) < 2L){return(fail("Only one observed age group"))}
      if(n.obs <= length(unique(d$age))){return(fail("No within-age replication for omnibus comparison"))}
      if(length(unique(d$score)) < 2L){return(fail("All scores identical"))}
      test.call <- function() stats::kruskal.test(d$score, factor(d$age))
    } else if(type == "Young_Old_Wilcoxon") {
      minimum <- if(unit == "mouse") min.mice.per.group.wilcoxon else min.cells.per.group.wilcoxon
      if(r$n_Young < minimum || r$n_Old < minimum) {
        return(
          fail(
            paste0(
              "Young and Old each need >= ", minimum, " ", unit,
              " observations for this test; observations retained in plots"
            )
          )
        )
      }
      if(length(unique(d$score)) < 2L){return(fail("All scores identical"))}
      use.exact <- unit == "mouse" && max(r$n_Young, r$n_Old) < 50L && !anyDuplicated(d$score)
      r$exact_used <- use.exact
      test.call <- function() stats::wilcox.test(
        d$score[d$group == "Young"], d$score[d$group == "Old"], alternative = "two.sided", paired = FALSE,
        exact = use.exact, correct = !use.exact
      )
    } else {stop("Unknown test type: ", type)}
    warnings <- character()
    fit <- tryCatch(
      withCallingHandlers(
        test.call(),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    r$warning <- paste(unique(warnings), collapse = " | ")
    if(inherits(fit, "error")){return(fail(conditionMessage(fit)))}
    if(length(fit$p.value) != 1L || !is.finite(fit$p.value)){return(fail("Nonfinite p-value"))}
    r$status <- "tested"
    r$reason <- "ok"
    r$p_value <- fit$p.value
    r$statistic <- unname(fit$statistic)
    r
  }
  bh.reference <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- is.finite(p)
    out[ok] <- stats::p.adjust(p[ok], method = "BH")
    out
  }
  plot.theme <- function(legend = FALSE) {
    theme_bw(base_size = plot.base.size) +
      theme(
        plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(),
        plot.caption = element_blank(), axis.title = element_text(size = plot.axis.title.size),
        axis.text = element_text(size = plot.axis.text.size, color = "black"),
        panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none",
        legend.title = element_text(size = plot.legend.size),
        legend.text = element_text(size = plot.legend.size), plot.margin = margin(12, 16, 12, 12)
      )
  }
  dotplot.theme <- function(legend = FALSE) {
    theme_bw(base_size = dotplot.base.size) +
      theme(
        plot.title = element_text(size = dotplot.title.size, face = "bold", margin = margin(b = 14)),
        plot.subtitle = element_blank(), plot.caption = element_blank(),
        axis.title = element_text(size = dotplot.axis.title.size),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"),
        axis.text.y = element_text(
          size = dotplot.celltype.text.size, color = "black", margin = margin(r = 10)
        ),
        panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none",
        legend.title = element_text(size = dotplot.legend.size),
        legend.text = element_text(size = dotplot.legend.size), legend.key.width = grid::unit(1.0, "cm"),
        plot.margin = margin(15, 25, 12, 12)
      )
  }
  plot.title <- function(tissue, ct, label) {
    paste(
      paste(
        strwrap(paste0(gsub("_", " ", tissue), " / ", ct), width = plot.title.wrap.width), collapse = "\n"
      ),
      label, sep = "\n"
    )
  }
  y.bounds <- function(y) {
    if(!length(y)){return(c(-0.1, 0.1))}
    b <- range(y)
    span <- max(diff(b), abs(mean(b)) * 0.10, 0.05)
    c(b[1] - span * 0.06, b[2] + span * 0.30)
  }
  age.plot <- function(d, unit, comparison, pvalue, title) {
    levels <- if(comparison == "age") age.order else c("Young", "Old")
    palette <- if(comparison == "age") age.colors else group.colors
    d$xgroup <- factor(if(comparison == "age") d$age else d$group, levels = levels)
    summaries <- d %>%
      group_by(xgroup) %>%
      summarise(n = n(), vmin = min(score), vmax = max(score), .groups = "drop")
    variable.groups <- as.character(summaries$xgroup[summaries$n >= 2 & summaries$vmin < summaries$vmax])
    variable <- d[as.character(d$xgroup) %in% variable.groups, , drop = FALSE]
    p <- ggplot(d, aes(x = xgroup, y = score)) + geom_blank()
    if(unit == "cell") {
      if(nrow(variable)) {
        p <- p +
          geom_violin(
            data = variable, aes(fill = xgroup), width = plot.violin.width, trim = TRUE, scale = "width",
            linewidth = 0.4
          ) +
          geom_boxplot(
            data = variable, width = plot.box.width, fill = "white", outlier.shape = NA, linewidth = 0.4
          )
      }
      special <- summaries[!as.character(summaries$xgroup) %in% variable.groups, , drop = FALSE]
      if(nrow(special)) {
        p <- p +
          geom_errorbar(
            data = special, aes(x = xgroup, ymin = vmin, ymax = vmax), inherit.aes = FALSE, width = 0.4
          ) +
          geom_point(data = special, aes(x = xgroup, y = vmin), inherit.aes = FALSE, size = 2)
      }
    } else {
      box.groups <- as.character(summaries$xgroup[summaries$n >= 2])
      if(length(box.groups)) {
        p <- p +
          geom_boxplot(
            data = d[as.character(d$xgroup) %in% box.groups, , drop = FALSE], width = 0.55, fill = "grey93",
            outlier.shape = NA
          )
      }
      p <- p +
        geom_point(
          aes(color = xgroup), position = position_jitter(width = 0.08, height = 0, seed = score.seed),
          size = mouse.point.size
        )
    }
    # Add only scales that a nonempty layer actually maps through aes().
    # Cell violins map fill; mouse points map color. Fixed white/grey boxes do not.
    if(unit == "cell" && nrow(variable) > 0L){p <- p + scale_fill_manual(values = palette)}
    if(unit == "mouse" && nrow(d) > 0L){p <- p + scale_color_manual(values = palette)}
    bounds <- y.bounds(d$score)
    p +
      scale_x_discrete(
        limits = levels, labels = if(comparison == "age") sub("m$", "", levels) else levels, drop = FALSE
      ) +
      annotate(
        "text", x = length(levels) + 0.35, y = bounds[2], label = format.p(pvalue), hjust = 1, vjust = 1,
        size = plot.pvalue.size, lineheight = 1.0
      ) +
      coord_cartesian(ylim = bounds) +
      labs(
        title = title, x = if(comparison == "age") "Age (months)" else NULL,
        y = if(unit == "cell") "AddModuleScore" else "Mean AddModuleScore per mouse"
      ) +
      plot.theme()
  }
  correlation.plot <- function(d, unit, stat, title, analysis.group = "All", color.by = "age") {
    if(length(analysis.group) != 1L || !analysis.group %in% c("All", "Young", "Old")) {
      stop("Invalid scatter group.")
    }
    if(length(color.by) != 1L || !color.by %in% c("age", "group")){stop("Invalid scatter color variable.")}
    # Use the same axes for All, Young-only and Old-only views of this cell type.
    bounds <- y.bounds(d$score)
    x.bounds <- if(nrow(d)) range(d$mff) else c(0, 1)
    if(diff(x.bounds) == 0){x.bounds <- x.bounds + c(-0.05, 0.05) * max(1, abs(x.bounds[1]))}
    if(analysis.group != "All") {
      d <- d[d$group == analysis.group, , drop = FALSE]
      title <- paste0(title, " / ", analysis.group)
    }
    stat <- stat[stat$analysis_group == analysis.group, , drop = FALSE]
    if(nrow(stat) != 1L){stop("Expected one Mff correlation result for ", analysis.group)}
    label <- paste0(
      "rho = ", if(is.finite(stat$rho)) sprintf("%.2f", stat$rho) else "NA", "\n", format.p(stat$p_value)
    )
    if(color.by == "age") {
      color.levels <-
        if(analysis.group == "All") age.order else if(analysis.group == "Young") young.ages else old.ages
      d$point_color <- factor(d$age, levels = color.levels)
      palette <- age.colors
      color.labels <- sub("m$", "", color.levels)
      legend.title <- "Age (months)"
    } else {
      color.levels <- if(analysis.group == "All") c("Young", "Old") else analysis.group
      d$point_color <- factor(d$group, levels = color.levels)
      palette <- group.colors
      color.labels <- color.levels
      legend.title <- "Age group"
    }
    # Each cell type has ONE plotting area, with no Young/Old facets.
    p <- ggplot(d, aes(x = mff, y = score)) +
      geom_point(
        aes(color = point_color), size = if(unit == "cell") cell.point.size else mouse.point.size,
        alpha = if(unit == "cell") cell.point.alpha else 0.95
      )
    if(nrow(d) > 0L) {
      p <- p +
        scale_color_manual(
          values = palette, breaks = color.levels, labels = color.labels, drop = FALSE, name = legend.title
        ) +
        guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 3)))
    }
    # Fit one OLS line to the observations shown in this view.
    # All uses pooled Young + Old; Young/Old views use their own observations.
    # Point colors never split the fit. Spearman rho/p use the same shown group.
    if(show.regression.line && nrow(d) >= 2L && length(unique(d$mff)) >= 2L) {
      p <- p +
        geom_smooth(
          aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, fullrange = FALSE,
          color = regression.line.color, linewidth = regression.line.width, show.legend = FALSE
        )
    }
    p +
      annotate(
        "text", x = Inf, y = bounds[2], label = label, hjust = 1.05, vjust = 1, size = plot.pvalue.size,
        lineheight = 1.0
      ) +
      coord_cartesian(xlim = x.bounds, ylim = bounds) +
      labs(
        title = title,
        x = if(unit == "cell") "Mff expression (RNA LogNormalize)" else "Mff pseudobulk log2(TMM CPM + 1)",
        y = if(unit == "cell") "AddModuleScore" else "Mean AddModuleScore per mouse"
      ) +
      plot.theme(legend = nrow(d) > 0L)
  }
  age.dotplot <- function(d, unit, title, geneset) {
    # Fixed-size dots: color is a signed mean score, never % score-positive cells.
    d$age_factor <- factor(d$age, levels = age.order)
    d$celltype_axis <- factor(d$celltype, levels = rev(celltypes))
    observed <- d[is.finite(d$mean_score), , drop = FALSE]
    missing <- d[!is.finite(d$mean_score), , drop = FALSE]
    score.colors <- if(geneset %in% red.when.old.higher) c(
      low = "#000000", mid = "#F3F3F3", high = "#D95F02"
    ) else c(low = "#D95F02", mid = "#F3F3F3", high = "#000000")
    p <- ggplot(d, aes(x = age_factor, y = celltype_axis)) + geom_blank()
    if(nrow(observed)) {
      limit <- max(abs(observed$mean_score), 0.05)
      p <- p +
        geom_point(data = observed, aes(color = mean_score), size = dotplot.point.size) +
        scale_color_gradient2(
          low = unname(score.colors["low"]), mid = unname(score.colors["mid"]),
          high = unname(score.colors["high"]), midpoint = 0, limits = c(-limit, limit),
          name = if(unit == "mouse") "Mean score\n(equal mice)" else "Mean score\n(equal cells)"
        )
    }
    if(nrow(missing)) {
      p <- p + geom_text(data = missing, label = "NA", color = "grey50", size = dotplot.na.text.size)
    }
    p +
      scale_x_discrete(limits = age.order, labels = sub("m$", "", age.order), drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      labs(title = title, x = "Age (months)", y = NULL) +
      dotplot.theme(legend = nrow(observed) > 0L) +
      theme(legend.position = if(nrow(observed)) "right" else "none")
  }
  difference.data <- function(statistics) {
    d <- as.data.frame(statistics[statistics$test == "Young_Old_Wilcoxon", , drop = FALSE])
    if(anyDuplicated(d$celltype)){stop("Repeated cell type in score-difference summary.")}
    d$effect_available <- is.finite(d$mean_difference_Old_minus_Young)
    d$effect_reason <- ifelse(
      d$effect_available, "estimated", "Missing Young or Old group; no score difference estimated"
    )
    d$direction <- ifelse(
      !d$effect_available, NA_character_,
      ifelse(
        d$mean_difference_Old_minus_Young > 0, "Old higher",
        ifelse(d$mean_difference_Old_minus_Young < 0, "Old lower", "Equal means")
      )
    )
    d$significance <- ifelse(
      !d$effect_available, NA_character_,
      ifelse(
        d$status != "tested" | !is.finite(d$p_value), "Not tested",
        ifelse(d$p_value < dotplot.p.cutoff, "Below cutoff", "At or above cutoff")
      )
    )
    d$effect_definition <- "Arithmetic mean score in Old minus arithmetic mean score in Young; not log2FC"
    d
  }
  difference.dotplot <- function(d, unit, title, geneset) {
    # Match the original Mff log2FC dotplot design while retaining the correct
    # AddModuleScore effect: arithmetic mean score in Old minus Young.
    d <- d[
      order(!d$effect_available, d$mean_difference_Old_minus_Young, d$celltype, na.last = TRUE), ,
      drop = FALSE
    ]
    d$celltype_axis <- factor(d$celltype, levels = rev(d$celltype))
    estimated <- d[is.finite(d$mean_difference_Old_minus_Young), , drop = FALSE]
    unavailable <- d[!is.finite(d$mean_difference_Old_minus_Young), , drop = FALSE]
    limit <- max(c(abs(estimated$mean_difference_Old_minus_Young), 0.05)) * 1.25
    x.breaks <- pretty(c(-limit, limit), n = 5)
    x.breaks <- sort(unique(c(0, x.breaks[x.breaks >= -limit & x.breaks <= limit])))
    direction.colors <- if(geneset %in% red.when.old.higher) c(
      "Old lower" = "#000000", "Old higher" = "#D95F02", "Equal means" = "grey50"
    ) else c("Old lower" = "#D95F02", "Old higher" = "#000000", "Equal means" = "grey50")
    significance.labels <- c(
      "Below cutoff" = paste0("Unadjusted p < ", dotplot.p.cutoff),
      "At or above cutoff" = paste0("Unadjusted p >= ", dotplot.p.cutoff), "Not tested" = "p = NA"
    )
    display.title <- paste0(gsub("_", " ", unique(d$tissue)[1]), " | ", score.labels[[geneset]])
    x.label <-
      if(unit == "cell") "Cell-mean AddModuleScore difference (Old - Young)" else
        "Mouse-mean AddModuleScore difference (Old - Young)"
    p <- ggplot(d, aes(y = celltype_axis)) +
      geom_vline(xintercept = 0, color = "#999999", linewidth = 0.6, linetype = "dashed") +
      geom_point(
        data = estimated, aes(x = mean_difference_Old_minus_Young, color = direction, shape = significance),
        size = dotplot.point.size, stroke = dotplot.point.stroke
      ) +
      scale_color_manual(values = direction.colors, guide = "none") +
      scale_shape_manual(
        values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4),
        limits = names(significance.labels), labels = unname(significance.labels), drop = FALSE, name = NULL
      ) +
      scale_x_continuous(breaks = x.breaks, limits = c(-limit, limit), expand = expansion(mult = 0.02)) +
      scale_y_discrete(drop = FALSE, expand = expansion(add = 0.7)) +
      labs(title = display.title, subtitle = NULL, x = x.label, y = NULL) +
      theme_classic(base_size = dotplot.base.size, base_family = "sans") +
      theme(
        plot.title = element_text(size = dotplot.title.size, face = "bold", margin = margin(b = 14)),
        axis.title.x = element_text(size = dotplot.axis.title.size, margin = margin(t = 12)),
        axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"),
        axis.text.y = element_text(
          size = dotplot.celltype.text.size, color = "black", margin = margin(r = 10)
        ),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.35), legend.position = "bottom",
        legend.text = element_text(size = dotplot.legend.size), legend.key.width = grid::unit(1.0, "cm"),
        plot.margin = margin(15, 25, 12, 12)
      ) +
      guides(shape = guide_legend(nrow = 1, override.aes = list(color = "#444444", size = 4)))
    if(nrow(unavailable)) {
      p <- p +
        geom_text(
          data = unavailable, aes(y = celltype_axis), x = limit * 0.96, label = "NA", inherit.aes = FALSE,
          hjust = 1, color = "#888888", size = dotplot.na.text.size
        )
    }
    p
  }
  save.plot <- function(
    p, stem, width = plot.width, height = plot.height, allow.pdf = save.pdf, res = plot.res
  ) {
    dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(filename = paste0(stem, ".png"), width = width, height = height, res = res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    if(allow.pdf) {
      ggplot2::ggsave(
        paste0(stem, ".pdf"), plot = p, device = grDevices::cairo_pdf, width = width / res,
        height = height / res, units = "in", limitsize = FALSE
      )
    }
  }
  save.overview <- function(plots, folder, label, panel.width = plot.width) {
    if(!save.overviews || !length(plots)){return(invisible(NULL))}
    # All cell-type panels in ONE image; 3 columns, rows expand as needed.
    ncols <- min(3L, length(plots))
    nrows <- ceiling(length(plots) / ncols)
    p <- patchwork::wrap_plots(plots, ncol = ncols)
    save.plot(
      p, file.path(folder, "overview", paste0(label, ".all_celltypes")), width = panel.width * ncols,
      height = plot.height * nrows, allow.pdf = FALSE
    )
  }

  # ==================================================================
  # 3. Read gene sets / input object and validate the selected cohort
  # ==================================================================
  genesets <- lapply(file.path(geneset.dir, gene.files), read.grp)
  names(genesets) <- names(gene.files)
  manifest <- data.frame(
    geneset = names(gene.files), file = unname(file.path(geneset.dir, gene.files)),
    md5 = unname(tools::md5sum(file.path(geneset.dir, gene.files))), n_input_genes = lengths(genesets),
    n_expected = unname(expected.sizes), stringsAsFactors = FALSE
  )
  manifest$size_matches_reference <- manifest$n_input_genes == manifest$n_expected
  write.csv.safe(manifest, file.path(output.base, "00.geneset_manifest.csv"))
  cat("\nInput GRP sizes (before matching to RNA genes):\n")
  print(manifest[, c("geneset", "n_input_genes", "n_expected", "size_matches_reference")], row.names = FALSE)
  if(any(!manifest$size_matches_reference)) {
    mismatch <- manifest[!manifest$size_matches_reference, , drop = FALSE]
    details <- paste0(
      mismatch$geneset, ": input=", mismatch$n_input_genes, ", reference=", mismatch$n_expected,
      collapse = "; "
    )
    warning(
      "GRP size differs from reference: ", details,
      ". Actual input lists are used; this is not RNA gene-matching loss. See 00.geneset_manifest.csv",
      call. = FALSE
    )
  }
  mmf.rds <- readRDS(rds.file)
  if(!inherits(mmf.rds, "Seurat") || !assay %in% names(mmf.rds@assays)) {
    stop("Expected a Seurat object with RNA assay.")
  }
  if(anyDuplicated(colnames(mmf.rds))){stop("Duplicated cell names in input.")}
  needed <- c(tissue.column, celltype.column, mouse.column, age.column)
  if(!all(needed %in% colnames(mmf.rds@meta.data))) {
    stop("Missing metadata: ", paste(setdiff(needed, colnames(mmf.rds@meta.data)), collapse = ", "))
  }
  SeuratObject::DefaultAssay(mmf.rds) <- assay
  tissue.values <- as.character(mmf.rds@meta.data[[tissue.column]])
  tissue.keys <- tissue.key(tissue.values)
  if(include.heart.and.aorta){tissue.keys[tissue.keys %in% "heart_and_aorta"] <- "heart"}
  if(anyNA(tissue.keys) || any(tissue.keys == "")){stop("Missing tissue annotation; inspect metadata.")}
  observed.ages <- as.character(mmf.rds@meta.data[[age.column]])
  # 18m와 21m은 분석 대상에서 제외한다. 이 선택은 AddModuleScore보다 먼저 적용된다.
  tissue.selected <- tissue.keys %in% tissue.key(tissues)
  if(
    anyNA(observed.ages[tissue.selected]) ||
    any(!observed.ages[tissue.selected] %in% c(age.order, excluded.ages))
  ) {
    stop(
      "Requested tissues contain missing/unexpected ages. Check age.order, excluded.ages and metadata."
    )
  }
  selected <- tissue.selected & observed.ages %in% age.order
  inventory <- data.frame(
    tissue_annotation = tissue.values, analysis_tissue_key = tissue.keys, age = observed.ages,
    selected = selected, stringsAsFactors = FALSE
  ) %>%
    count(tissue_annotation, analysis_tissue_key, age, selected, name = "n_cells")
  write.csv.safe(inventory, file.path(output.base, "01.input_cohort_inventory.csv"))
  cat("\nSelected tissue/age inventory:\n")
  print(as.data.frame(inventory), row.names = FALSE)
  all.stats <- list()
  all.coverage <- list()

  # ==================================================================
  # 4. Score each tissue ONCE, then summarize the same cells by mouse
  # ==================================================================
  for(tissue.name in tissues) {
    cat("\n========== ", tissue.name, " ==========\n", sep = "")
    selected.cells <- colnames(mmf.rds)[
      tissue.keys == tissue.key(tissue.name) & observed.ages %in% age.order
    ]
    if(!length(selected.cells)){stop("No selected cells for tissue: ", tissue.name)}
    tissue.dir <- file.path(output.base, tissue.name)
    data.dir <- file.path(tissue.dir, "00_data")
    dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)
    tissue.rds <- subset(mmf.rds, cells = selected.cells)
    tissue.rds <- join.layer(tissue.rds, "counts")
    tissue.rds <- join.layer(tissue.rds, "data")
    normalized <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "data")
    counts <- SeuratObject::LayerData(tissue.rds[[assay]], layer = "counts")
    if(
      !setequal(
        colnames(normalized), colnames(tissue.rds)
      ) || !setequal(colnames(counts), colnames(tissue.rds))
    ) {
      stop("RNA layer cell coverage mismatch.")
    }
    normalized <- normalized[, colnames(tissue.rds), drop = FALSE]
    counts <- counts[, colnames(tissue.rds), drop = FALSE]
    if(
      anyDuplicated(rownames(normalized)) || anyDuplicated(rownames(counts)) || !gene %in% intersect(
        rownames(normalized), rownames(counts)
      )
    ) {
      stop("Missing Mff or duplicate RNA feature names.")
    }
    norm.values <- if(inherits(normalized, "sparseMatrix")) normalized@x else as.vector(normalized)
    count.values <- if(inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
    if(any(!is.finite(norm.values)) || any(norm.values < 0)) {
      stop("Expected nonnegative RNA LogNormalize data, not scale.data.")
    }
    if(
      any(!is.finite(count.values)) ||
      any(count.values < 0) ||
      any(abs(count.values - round(count.values)) > 1e-8)
    ) {
      stop("Expected nonnegative integer raw UMI counts.")
    }
    rm(norm.values, count.values)
    meta <- tissue.rds@meta.data[colnames(tissue.rds), , drop = FALSE]
    ct <- as.character(meta[[celltype.column]])
    ct[is.na(ct) | trimws(ct) == ""] <- "Unannotated"
    cells <- data.frame(
      cell = colnames(tissue.rds), tissue = tissue.name,
      tissue_free_annotation = as.character(meta[[tissue.column]]), celltype = ct,
      mouse.id = as.character(meta[[mouse.column]]), age = as.character(meta[[age.column]]),
      technical_run = if("orig.ident" %in% names(meta)) as.character(meta$orig.ident) else NA_character_,
      sex = if("sex" %in% names(meta)) as.character(meta$sex) else NA_character_,
      Mff_LogNormalize = as.numeric(normalized[gene, ]), Mff_raw_UMI = as.numeric(counts[gene, ]),
      nCount_RNA = as.numeric(Matrix::colSums(counts)),
      nFeature_RNA = as.numeric(Matrix::colSums(counts > 0)),
      percent.mt = if("percent.mt" %in% names(meta)) as.numeric(meta$percent.mt) else NA_real_,
      stringsAsFactors = FALSE
    )
    if(anyNA(cells$mouse.id) || any(trimws(cells$mouse.id) == "")){stop("Missing mouse.id in ", tissue.name)}
    if(any(cells$nCount_RNA <= 0)){stop("Cells with zero total raw RNA counts in ", tissue.name)}
    cells$age_months <- as.numeric(sub("m$", "", cells$age))
    cells$group <- ifelse(cells$age %in% young.ages, "Young", "Old")
    if(
      "age_young_old" %in% names(meta) &&
      (anyNA(meta$age_young_old) || any(as.character(meta$age_young_old) != cells$group))
    ) {
      stop("Existing age_young_old metadata disagrees with age-based Young/Old labels.")
    }
    mouse.age <- unique(cells[, c("mouse.id", "age")])
    if(anyDuplicated(mouse.age$mouse.id)){stop("One mouse maps to multiple ages in ", tissue.name)}
    sex.check <- cells %>%
      filter(!is.na(sex), trimws(sex) != "") %>%
      distinct(mouse.id, sex) %>%
      count(mouse.id)
    if(any(sex.check$n > 1)){stop("One mouse maps to multiple sex labels.")}

    # Mff is excluded from both the control pool and any signature containing it.
    # GRP files on disk are unchanged. Every excluded/missing gene is recorded.
    features <- lapply(genesets, function(g) setdiff(intersect(g, rownames(normalized)), gene))
    gene.audit <- do.call(rbind, lapply(names(genesets), function(gs) {
      g <- genesets[[gs]]
      present <- g %in% rownames(normalized)
      avg <- rep(NA_real_, length(g))
      if(any(present)){avg[present] <- as.numeric(Matrix::rowMeans(normalized[g[present], , drop = FALSE]))}
      data.frame(
        tissue = tissue.name, geneset = gs, gene = g, in_RNA_data = present, mean_LogNormalize = avg,
        status = ifelse(g == gene, "excluded_Mff", ifelse(present, "used", "not_in_RNA_data")),
        stringsAsFactors = FALSE
      )
    }))
    coverage <- data.frame(
      tissue = tissue.name, geneset = names(genesets), score_column = unname(score.columns),
      n_input = lengths(genesets), n_used = lengths(features),
      n_missing_RNA = vapply(genesets, function(g) sum(!g %in% rownames(normalized)), integer(1)),
      n_excluded_Mff = vapply(genesets, function(g) sum(g == gene), integer(1)), stringsAsFactors = FALSE
    )
    coverage$fraction_used <- coverage$n_used / coverage$n_input
    coverage$n_used_with_nonzero_tissue_expression <- vapply(
      names(features),
      function(gs) sum(
        gene.audit$geneset == gs & gene.audit$status == "used" & gene.audit$mean_LogNormalize > 0,
        na.rm = TRUE
      ),
      integer(1)
    )
    write.csv.safe(gene.audit, file.path(data.dir, "01.gene_membership_audit.csv"))
    write.csv.safe(coverage, file.path(data.dir, "02.geneset_coverage.csv"))
    all.coverage[[tissue.name]] <- coverage
    print(coverage[, c("geneset", "n_input", "n_used", "n_missing_RNA", "n_excluded_Mff")], row.names = FALSE)
    if(any(coverage$n_used < min.matched.genes)) {
      stop(
        "Insufficient matched genes; see gene_membership_audit.csv. Mouse symbols must match RNA row names exactly."
      )
    }
    if(any(coverage$fraction_used < low.coverage.warning)) {
      warning("Low gene-set coverage in ", tissue.name, "; inspect 02.geneset_coverage.csv", call. = FALSE)
    }
    if(any(coverage$n_excluded_Mff > 0)) {
      warning("Mff was removed from a signature to avoid self-correlation; see gene audit.", call. = FALSE)
    }
    used.dir <- file.path(data.dir, "used_genesets")
    dir.create(used.dir, recursive = TRUE, showWarnings = FALSE)
    for(gs in names(features)){writeLines(features[[gs]], file.path(used.dir, paste0(gs, ".used.grp")))}
    pool <- setdiff(rownames(normalized), gene)
    stopifnot(!gene %in% pool, all(unlist(features, use.names = FALSE) %in% pool))
    writeLines(pool, file.path(data.dir, "03.control_candidate_pool.grp"))
    temporary.names <- paste0("AgingMFF_TMP", seq_along(features))
    if(any(temporary.names %in% colnames(tissue.rds@meta.data))) {
      stop("Temporary score columns already exist; remove/rename AgingMFF_TMP columns in a copy of the input.")
    }
    tissue.rds <- tryCatch(
      Seurat::AddModuleScore(
        object = tissue.rds, features = unname(features), pool = pool, nbin = score.nbin, ctrl = score.ctrl,
        assay = assay, name = "AgingMFF_TMP", seed = score.seed, search = FALSE, slot = "data"
      ),
      error = function(e) stop(
        "AddModuleScore failed in ", tissue.name, ": ", conditionMessage(e),
        "\nCheck matched genes and score.ctrl/score.nbin; these settings were not silently changed.",
        call. = FALSE
      )
    )
    for(i in seq_along(features)) {
      values <- tissue.rds@meta.data[colnames(tissue.rds), temporary.names[i]]
      if(length(values) != nrow(cells) || any(!is.finite(values))) {
        stop("Invalid module scores: ", names(features)[i])
      }
      tissue.rds[[score.columns[i]]] <- setNames(values, cells$cell)
      tissue.rds[[temporary.names[i]]] <- NULL
      cells[[score.columns[i]]] <- values
    }
    tissue.rds@misc$Aging_MFF_genescore <- list(
      timestamp = as.character(Sys.time()), input = rds.file, manifest = manifest, features = features,
      coverage = coverage, assay = assay, slot = "data", seed = score.seed, nbin = score.nbin,
      ctrl = score.ctrl, control_candidate_pool = pool,
      scoring_scope = paste(tissue.name, "all selected ages and cell types together"), Mff_excluded = gene,
      include_heart_and_aorta = include.heart.and.aorta,
      Seurat_version = as.character(packageVersion("Seurat"))
    )
    write.csv.safe(cells, file.path(data.dir, "04.cell_scores.csv"))
    saveRDS(cells, file.path(data.dir, "04.cell_scores.rds"))
    if(save.scored.seurat) {
      saveRDS(
        tissue.rds, file.path(data.dir, paste0("Aging.MFF.", tissue.name, ".genescores.RDS")),
        compress = "gzip"
      )
    }

    mice <- cells %>%
      group_by(tissue, celltype, mouse.id, age, age_months, group) %>%
      summarise(
        n_cells = n(), sex = collapse.values(sex), technical_runs = collapse.values(technical_run),
        tissue_annotations = collapse.values(tissue_free_annotation),
        Mff_mean_LogNormalize = mean(Mff_LogNormalize), Mff_percent_detected = mean(Mff_raw_UMI > 0) * 100,
        across(all_of(unname(score.columns)), mean), .groups = "drop"
      ) %>%
      arrange(celltype, age_months, mouse.id)
    mice$pb_id <- sprintf("PB%05d", seq_len(nrow(mice)))
    mapping <- cells %>%
      select(cell, celltype, mouse.id) %>%
      left_join(mice %>% select(celltype, mouse.id, pb_id), by = c("celltype", "mouse.id"))
    if(nrow(mapping) != nrow(cells) || anyNA(mapping$pb_id) || !identical(mapping$cell, colnames(tissue.rds))) {
      stop("Cell-to-mouse mapping failed.")
    }
    tissue.rds$AgingMFF_pb_id <- setNames(mapping$pb_id, mapping$cell)
    pb <- Seurat::AggregateExpression(
      tissue.rds, assays = assay, features = rownames(counts), group.by = "AgingMFF_pb_id",
      return.seurat = FALSE, verbose = FALSE
    )[[assay]]
    if(!setequal(rownames(pb), rownames(counts)) || !setequal(colnames(pb), mice$pb_id)) {
      stop("Pseudobulk names/coverage mismatch.")
    }
    pb <- pb[rownames(counts), mice$pb_id, drop = FALSE]
    expected <- rowsum(
      cbind(total_UMI = cells$nCount_RNA, Mff_UMI = cells$Mff_raw_UMI), group = mapping$pb_id, reorder = FALSE
    )
    expected <- expected[mice$pb_id, , drop = FALSE]
    if(
      any(abs(Matrix::rowSums(pb) - Matrix::rowSums(counts)) > 1e-8) ||
      any(abs(Matrix::colSums(pb) - expected[, "total_UMI"]) > 1e-8) ||
      any(abs(as.numeric(pb[gene, ]) - expected[, "Mff_UMI"]) > 1e-8)
    ) {
      stop("Pseudobulk sum verification failed.")
    }
    mice$library_size <- as.numeric(Matrix::colSums(pb))
    mice$Mff_raw_UMI <- as.numeric(pb[gene, ])
    mice$eligible <- mice$n_cells >= min.cells.per.mouse.celltype & mice$library_size > 0
    mice$eligibility_reason <- ifelse(
      mice$eligible, "included", ifelse(mice$library_size <= 0, "zero_library_size", "below_min_cells")
    )
    mice$TMM_factor <- NA_real_
    mice$Mff_TMM_CPM <- NA_real_
    mice$Mff_log2_TMM_CPM_plus1 <- NA_real_
    mice$normalization <- "not_normalized_excluded"
    for(ct in unique(mice$celltype)) {
      idx <- which(mice$celltype == ct & mice$eligible)
      if(!length(idx)){next}
      mat <- as.matrix(pb[, mice$pb_id[idx], drop = FALSE])
      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      y <- edgeR::DGEList(counts = mat)
      if(length(idx) >= 2L){y <- edgeR::calcNormFactors(y, method = "TMM")} else {y$samples$norm.factors <- 1}
      cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
      mff.cpm <- if(gene %in% rownames(cpm)) as.numeric(cpm[gene, ]) else rep(0, length(idx))
      if(
        any(!is.finite(mff.cpm)) ||
        any(mff.cpm < 0) ||
        any(!is.finite(y$samples$norm.factors)) ||
        any(y$samples$norm.factors <= 0)
      ) {
        stop("Invalid TMM/CPM values in ", ct)
      }
      mice$TMM_factor[idx] <- y$samples$norm.factors
      mice$Mff_TMM_CPM[idx] <- mff.cpm
      mice$Mff_log2_TMM_CPM_plus1[idx] <- log2(mff.cpm + 1)
      mice$normalization[idx] <- if(length(idx) >= 2L) "TMM" else "CPM_only_single_mouse"
    }
    write.csv.safe(mapping, file.path(data.dir, "05.cell_to_pseudobulk.csv"))
    write.csv.safe(mice, file.path(data.dir, "06.mouse_scores_Mff_pseudobulk.csv"))
    saveRDS(mice, file.path(data.dir, "06.mouse_scores_Mff_pseudobulk.rds"))
    saveRDS(pb, file.path(data.dir, "07.pseudobulk_raw_counts.RDS"))
    celltypes <- sort(unique(cells$celltype))
    file.map <- data.frame(
      celltype = celltypes,
      file_stem = sprintf(
        "%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 90)
      ),
      stringsAsFactors = FALSE
    )
    write.csv.safe(file.map, file.path(data.dir, "08.celltype_file_map.csv"))
    rm(tissue.rds, normalized, counts, meta, pb, expected)
    invisible(gc())

    # ================================================================
    # 5. Four requested analysis families, per tissue and cell type
    # ================================================================
    tissue.stats <- list()
    for(gs in names(score.columns)) {
      cat("Plots and statistics: ", tissue.name, " / ", gs, "\n", sep = "")
      score.dir <- file.path(tissue.dir, gs)
      sc <- score.columns[[gs]]
      for(unit in c("mouse", "cell")) {
        source <- if(unit == "mouse") as.data.frame(mice[mice$eligible, ]) else as.data.frame(cells)
        d.all <- source[, c("tissue", "celltype", "mouse.id", "age", "age_months", "group"), drop = FALSE]
        d.all$score <- source[[sc]]
        d.all$mff <- if(unit == "mouse") source$Mff_log2_TMM_CPM_plus1 else source$Mff_LogNormalize
        if(any(!is.finite(d.all$score)) || any(!is.finite(d.all$mff))) {
          stop("Nonfinite values in plotting/testing input.")
        }
        age.folder <- file.path(score.dir, if(unit == "mouse") "01_mouse_age_score" else "03_cell_age_score")
        mff.folder <- file.path(score.dir, if(unit == "mouse") "02_mouse_Mff_score" else "04_cell_Mff_score")
        dir.create(age.folder, recursive = TRUE, showWarnings = FALSE)
        dir.create(mff.folder, recursive = TRUE, showWarnings = FALSE)
        unit.stats <- list()
        age.plots <- list()
        yo.plots <- list()
        corr.plots <- list()
        corr.groupcolor.plots <- list()
        corr.young.plots <- list()
        corr.old.plots <- list()
        summaries <- list()
        for(i in seq_along(celltypes)) {
          ct <- celltypes[i]
          d <- d.all[d.all$celltype == ct, , drop = FALSE]
          tests <- c("Age_Kruskal_Wallis", "Age_Spearman", "Mff_Spearman", "Young_Old_Wilcoxon")
          st <- do.call(rbind, lapply(tests, function(tt) safe.test(d, tt, unit)))
          st$analysis_group <- "All"
          for(group.name in c("Young", "Old")) {
            group.stat <- safe.test(d[d$group == group.name, , drop = FALSE], "Mff_Spearman", unit)
            group.stat$test <- paste0("Mff_Spearman_", group.name)
            group.stat$analysis_group <- group.name
            st <- rbind(st, group.stat)
          }
          st$tissue <- tissue.name
          st$celltype <- ct
          st$geneset <- gs
          st$n_genes_used <- coverage$n_used[coverage$geneset == gs]
          unit.stats[[ct]] <- st
          summaries[[ct]] <- d %>%
            group_by(tissue, celltype, age, age_months, group) %>%
            summarise(
              n_observations = n(), n_mice = n_distinct(mouse.id), mean_score = mean(score),
              median_score = median(score), sd_score = if(n() > 1L) sd(score) else NA_real_,
              min_score = min(score), max_score = max(score), .groups = "drop"
            )
          title <- plot.title(tissue.name, ct, score.labels[[gs]])
          age.plots[[ct]] <- age.plot(d, unit, "age", st$p_value[st$test == "Age_Kruskal_Wallis"], title)
          mff.stats <- st[grepl("^Mff_Spearman", st$test), ]
          corr.plots[[ct]] <- correlation.plot(
            d, unit, mff.stats, title, analysis.group = "All", color.by = "age"
          )
          corr.groupcolor.plots[[ct]] <- correlation.plot(
            d, unit, mff.stats, title, analysis.group = "All", color.by = "group"
          )
          if(save.overviews) {
            corr.young.plots[[ct]] <- correlation.plot(
              d, unit, mff.stats, title, analysis.group = "Young", color.by = "age"
            )
            corr.old.plots[[ct]] <- correlation.plot(
              d, unit, mff.stats, title, analysis.group = "Old", color.by = "age"
            )
          }
          stem <- file.map$file_stem[i]
          save.plot(age.plots[[ct]], file.path(age.folder, "selected_4ages", paste0(stem, ".age_score")))
          save.plot(
            corr.plots[[ct]], file.path(mff.folder, "scatter", paste0(stem, ".Mff_score")),
            width = scatter.width
          )
          save.plot(
            corr.groupcolor.plots[[ct]],
            file.path(mff.folder, "scatter", paste0(stem, ".Mff_score.by_YoungOld")), width = scatter.width
          )
          if(make.young.old.plots) {
            yo.plots[[ct]] <- age.plot(d, unit, "group", st$p_value[st$test == "Young_Old_Wilcoxon"], title)
            save.plot(yo.plots[[ct]], file.path(age.folder, "young_old", paste0(stem, ".Young_Old_score")))
          }
        }
        statistics <- do.call(rbind, unit.stats)
        rownames(statistics) <- NULL
        tissue.stats[[paste(gs, unit, sep = ".")]] <- statistics
        # Raw p displayed. Reference BH across cell types AND all five scores is added later.
        age.summary <- bind_rows(summaries) %>% arrange(celltype, age_months)
        write.csv.safe(age.summary, file.path(age.folder, "age_score_summary.csv"))
        if(make.age.dotplots) {
          dot.folder <- file.path(age.folder, "dotplot")
          dir.create(dot.folder, recursive = TRUE, showWarnings = FALSE)
          # Complete the display grid; missing means remain NA, never zero.
          dot.data <- expand.grid(celltype = celltypes, age = age.order, stringsAsFactors = FALSE) %>%
            left_join(age.summary, by = c("celltype", "age"))
          dot.data$tissue <- tissue.name
          dot.data$geneset <- gs
          dot.data$unit <- unit
          dot.data$age_months <- as.numeric(sub("m$", "", dot.data$age))
          dot.data$group <- ifelse(dot.data$age %in% young.ages, "Young", "Old")
          dot.data$n_observations[is.na(dot.data$n_observations)] <- 0L
          dot.data$n_mice[is.na(dot.data$n_mice)] <- 0L
          dot.data$summary_weighting <-
            if(unit == "mouse") "Equal weight per eligible mouse" else
              "Equal weight per cell; not balanced by mouse"
          dot.data$estimate_status <- ifelse(is.finite(dot.data$mean_score), "estimated", "no_observations")
          dot.data <- dot.data %>% arrange(celltype, age_months)
          write.csv.safe(dot.data, file.path(dot.folder, "age_mean_score.plot_data.csv"))
          dot.title <- paste0(gsub("_", " ", tissue.name), " / ", unit, "\n", score.labels[[gs]])
          dot.width <- round(dotplot.width.in * dotplot.res)
          dot.height <- round(
            max(
              dotplot.height.minimum.in, 2.0 + length(celltypes) * dotplot.height.per.celltype.in
            ) * dotplot.res
          )
          save.plot(
            age.dotplot(dot.data, unit, dot.title, geneset = gs),
            file.path(dot.folder, "age_mean_score.all_celltypes"), width = dot.width, height = dot.height,
            res = dotplot.res
          )
          # Reuses the Young/Old Wilcoxon result; no extra test for the dotplot.
          save.plot(
            difference.dotplot(difference.data(statistics), unit, dot.title, geneset = gs),
            file.path(dot.folder, "Young_Old_score_difference.all_celltypes"), width = dot.width,
            height = dot.height, res = dotplot.res
          )
        }
        save.overview(age.plots, age.folder, "selected_4ages")
        save.overview(yo.plots, age.folder, "young_old")
        save.overview(corr.plots, mff.folder, "Mff_score", panel.width = scatter.width)
        save.overview(corr.groupcolor.plots, mff.folder, "Mff_score.by_YoungOld", panel.width = scatter.width)
        save.overview(corr.young.plots, mff.folder, "Mff_score.Young")
        save.overview(corr.old.plots, mff.folder, "Mff_score.Old")
      }
    }
    statistics <- do.call(rbind, tissue.stats)
    rownames(statistics) <- NULL
    statistics <- statistics %>%
      group_by(tissue, unit, test) %>%
      mutate(
        p_adj_BH = bh.reference(p_value), n_tests_in_BH_family = sum(is.finite(p_value)),
        p_adjust_scope = "Within tissue x unit x test, across all five signatures and all cell types"
      ) %>%
      ungroup()
    write.csv.safe(statistics, file.path(tissue.dir, "all_statistics.csv"))
    for(gs in names(score.columns)) {
      for(unit in c("mouse", "cell")) {
        s <- statistics[statistics$geneset == gs & statistics$unit == unit, ]
        age.folder <- file.path(
          tissue.dir, gs, if(unit == "mouse") "01_mouse_age_score" else "03_cell_age_score"
        )
        mff.folder <- file.path(
          tissue.dir, gs, if(unit == "mouse") "02_mouse_Mff_score" else "04_cell_Mff_score"
        )
        is.mff <- grepl("^Mff_Spearman", s$test)
        write.csv.safe(s[!is.mff, ], file.path(age.folder, "statistics.csv"))
        write.csv.safe(s[is.mff, ], file.path(mff.folder, "statistics.csv"))
        if(make.age.dotplots) {
          write.csv.safe(
            difference.data(s), file.path(age.folder, "dotplot", "Young_Old_score_difference.plot_data.csv")
          )
        }
      }
    }
    all.stats[[tissue.name]] <- statistics
    rm(
      cells, mice, source, d.all, d, age.plots, yo.plots, corr.plots, corr.groupcolor.plots, corr.young.plots,
      corr.old.plots, tissue.stats
    )
    invisible(gc())
    cat("Saved: ", tissue.dir, "\n", sep = "")
  }

  # ==================================================================
  # 6. Combined results and reproducibility information
  # ==================================================================
  combined <- bind_rows(all.stats)
  write.csv.safe(combined, file.path(output.base, "02.all_statistics.csv"))
  write.csv.safe(bind_rows(all.coverage), file.path(output.base, "03.all_geneset_coverage.csv"))
  notes <- c(
    paste0("Completed: ", Sys.time()), paste0("Input: ", rds.file), paste0("Output: ", output.base),
    paste0("Heart includes Heart_and_Aorta: ", include.heart.and.aorta),
    paste0("Ages: ", paste(age.order, collapse = ", ")),
    paste0("Excluded ages: ", paste(excluded.ages, collapse = ", ")),
    paste0("Young: ", paste(young.ages, collapse = ", "), "; Old: ", paste(old.ages, collapse = ", ")), "",
    "SCORING", "Five AddModuleScore signatures; CoreScence up and down remain separate",
    "CoreScence is not the DeepScence model; no published DeepScence cutoffs/performance are implied",
    "Higher CoreScence-down score is not interpreted as more senescence",
    "Each tissue scored once across all selected ages/cell types; same control reference within that tissue",
    "Do not compare absolute score levels between tissues or between different gene sets",
    "RNA data layer must already be LogNormalize; scale.data and regressed values are not used",
    paste0(
      "seed=", score.seed, "; nbin=", score.nbin, "; ctrl=", score.ctrl, "; Seurat=", packageVersion("Seurat")
    ),
    "Mff excluded from the control candidate pool and from any input signature containing it; exclusions recorded",
    "Presence matching is exact and case-sensitive; no HGNC lookup, capitalization conversion or imputation",
    "Row-present genes with all-zero expression are retained and recorded in the gene audit",
    "GRP lists, candidate pool, versions, seed and source checksums are retained for reproducibility",
    "The control candidate pool is not the list of sampled control genes", "", "MOUSE SUMMARIES",
    "Mff: raw RNA UMI summed per tissue x cell type x mouse; technical runs combined",
    "TMM fit jointly across eligible mice in each tissue/cell type; only all-zero genes removed",
    "Mff plotted as log2(TMM-normalized CPM + 1); a single mouse uses CPM with factor=1",
    "Score: arithmetic mean of the SAME cells' AddModuleScores within that mouse/cell type",
    "Mean module score is not a score computed from pseudobulk counts and is not TMM-normalized",
    paste0("Minimum cells per mouse/cell type: ", min.cells.per.mouse.celltype),
    "Equal weight per mouse in mouse-level tests; mean-score precision varies with the number of cells", "",
    "STATISTICS", "Selected four-age plot p: Kruskal-Wallis omnibus test across observed age groups",
    "Age_Spearman in CSV: numerical age in months vs score; different from omnibus differences",
    "Young/Old plot p: two-sided Wilcoxon rank-sum; effect = mean_Old minus mean_Young (not log fold change)",
    "Mff scatter: one plotting area per cell type; All views pool Young and Old; Young/Old-only views retain the same axes",
    "All views display pooled Spearman rho/raw p; Young-only and Old-only views display subgroup Spearman rho/raw p; these are not regression slope test results",
    "Each regression line is ordinary least squares score ~ mff on the observations shown: pooled for All, subgroup-specific for Young or Old",
    "All-data color variants (age or Young/Old) use the same observations, pooled line and pooled Spearman rho/p; point colors do not split the fit",
    "Regression lines are descriptive and are not adjusted for age, sex, batch or within-mouse cell dependence",
    "A regression line requires at least two observations and two distinct Mff values; otherwise points remain without a line",
    "Regression is drawn only over Mff values observed in the displayed group, without extrapolation or confidence bands",
    paste0("Draw regression line: ", show.regression.line),
    paste0("Cell point alpha: ", cell.point.alpha, "; point size: ", cell.point.size),
    "Mff_Spearman in CSV retains the pooled result; Mff_Spearman_Young and Mff_Spearman_Old are separate subgroup tests",
    "analysis_group identifies All, Young or Old; empty/constant/too-small groups retain NA and a reason",
    "All zero-expression Mff cells are retained; no imputation or downsampling",
    "Cellwise tests do not model dependence among cells from the same mouse; exploratory only",
    "All associations are unadjusted for age, sex, batch, RNA depth or detection rate",
    "Young/Old stratification is not continuous-age adjustment: variation within Young and within Old remains",
    "A pooled Mff-score correlation can reflect between-age differences; both All views display the same pooled association",
    "Raw p displayed; BH reference stored across all five scores/cell types within tissue x unit x test",
    "Pooled, Young and Old Mff correlations have distinct test names and separate BH families",
    "Spearman: exact only for <10 observations without ties; otherwise asymptotic approximation",
    "Small-mouse-sample asymptotic p-values and correlations can be unstable",
    "Unavailable tests retain the group/cell type and record a reason; no p-value is imputed",
    paste0("Minimum observations for correlation test: ", min.observations.correlation),
    paste0(
      "Minimum Young/Old observations: mouse=", min.mice.per.group.wilcoxon, "; cell=",
      min.cells.per.group.wilcoxon
    ),
    "Negative scores are retained and shown; score zero is not a senescence classification cutoff", "",
    "OUTPUT", "01_mouse_age_score / 02_mouse_Mff_score / 03_cell_age_score / 04_cell_Mff_score",
    "Age folders contain selected-four-age and optional Young/Old figures; 00_data contains reusable scores and QC",
    "Each overview includes ALL cell types in one image with three columns, without pagination",
    "Mff overviews: Mff_score.all_celltypes.png (All, age colors), Mff_score.by_YoungOld.all_celltypes.png (All, Young/Old colors), Mff_score.Young.all_celltypes.png (Young only), Mff_score.Old.all_celltypes.png (Old only)",
    "All four Mff overviews use the same cell-type order; axes for each cell type are derived jointly from Young and Old before selecting a group",
    "All four Mff overviews are written for both mouse and cell units under their overview directories",
    "Age dotplot: x=age, y=cell type, fixed-size dots colored by signed mean AddModuleScore",
    "Age dotplot means: mouse=mean of eligible mouse means (equal mice), cell=mean of all cells (equal cells)",
    "Age dotplot color scale is symmetric about zero within each tissue/signature/unit; not standardized across different signatures or tissues",
    "Missing celltype-age combinations are shown as NA, not assigned score zero",
    "Young/Old dotplot: x=mean_Old minus mean_Young score; positive means higher score in Old",
    "Do not exponentiate scores or compute log2FC; AddModuleScore is signed",
    "Difference dotplot: red=Old lower, black=Old higher, grey=equal means",
    paste0(
      "Difference dotplot shapes: filled raw Wilcoxon p<", dotplot.p.cutoff, "; open p>=", dotplot.p.cutoff,
      "; cross if effect is available but p=NA"
    ),
    "Difference dotplot shows NA without a point for a missing Young or Old group",
    "Wilcoxon p assesses rank distributions, not a test specifically of the plotted mean difference; no confidence intervals calculated",
    "Dotplot raw p values reuse existing Young/Old tests; reference BH is retained in plot_data CSV",
    "Cell dotplots and p values remain exploratory; pooled cells are not balanced by mouse or by age",
    "Dotplots: age folders/dotplot/*.all_celltypes.png, with matching plot_data.csv files",
    "Original RDS is not overwritten; new tissue-specific RDS files have AMS_* columns",
    "No normalization, PCA, UMAP or clustering is rerun", "", "SOURCES",
    "https://satijalab.org/seurat/reference/addmodulescore",
    "https://satijalab.org/seurat/reference/aggregateexpression",
    "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.test.html",
    "https://ggplot2.tidyverse.org/reference/geom_smooth.html", "", capture.output(sessionInfo())
  )
  notes <- c(
    notes, "", "FIGURE STYLE MATCH",
    paste0("Individual age/scatter figures: ", plot.width, " x ", plot.height, " px at ", plot.res, " dpi"),
    paste0(
      "Main text sizes: base=", plot.base.size, "; title=", plot.title.size, "; axis title=",
      plot.axis.title.size, "; axis tick=", plot.axis.text.size, "; legend=", plot.legend.size, " pt"
    ),
    paste0(
      "P-value and rho annotation size: ", plot.pvalue.size, " mm; title wrap width=", plot.title.wrap.width
    ),
    "Overview retains all cell types in three columns; canvas equals individual panel dimensions times row/column counts",
    paste0(
      "Dotplots: width=", dotplot.width.in, " inches; height=max(", dotplot.height.minimum.in,
      ", 2 + n_celltypes * ", dotplot.height.per.celltype.in, ") inches; dpi=", dotplot.res
    ),
    paste0(
      "Dotplot text sizes: base=", dotplot.base.size, "; title=", dotplot.title.size, "; axis title=",
      dotplot.axis.title.size, "; axis number=", dotplot.axis.number.size, "; cell type=",
      dotplot.celltype.text.size, "; legend=", dotplot.legend.size, " pt"
    ),
    "The dotplot reference uses 300 dpi; main age/scatter figures use 160 dpi",
    "Figure style and four scatter views follow the supplied reference; input/output paths, two tissues, five signatures and scoring/statistical settings come from the supplied full script"
  )
  notes <- notes[notes != "Difference dotplot: red=Old lower, black=Old higher, grey=equal means"]
  notes <- c(
    notes, "COLOR DIRECTION",
    "Color direction is a display convention; dot fill/open shape, not color, represents the raw-p cutoff",
    "Age dotplot: Hallmark_IFN_alpha/GO_IFN_beta/SaulSenMayo/CoreScence_up use black=negative, white=zero, red=positive; CoreScence_down uses red=negative, white=zero, black=positive",
    "Difference dotplot: Hallmark_IFN_alpha/GO_IFN_beta/SaulSenMayo/CoreScence_up use red=Old higher and black=Old lower; CoreScence_down uses red=Old lower and black=Old higher; grey=equal means"
  )
  writeLines(notes, file.path(output.base, "04.analysis_settings.txt"), useBytes = TRUE)
  cat("\nCompleted all four analysis families.\nOutput: ", output.base, "\n", sep = "")
  cat("Statistics: 02.all_statistics.csv\nGene coverage: 03.all_geneset_coverage.csv\n")
  print(as.data.frame(combined %>% count(unit, test, status)), row.names = FALSE)
  invisible(list(output = output.base, statistics = combined, coverage = bind_rows(all.coverage)))
}

# Runs automatically when you source() this complete file.
run.aging.mff.genescores()
