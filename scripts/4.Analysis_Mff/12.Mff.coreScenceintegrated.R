#=====================================================================
# Standalone NEW CELL: CoreScence net score analysis
#
# Prerequisite:
#   Run Aging_MFF_genescores_reference_style.R once so that each tissue's
#   00_data/04.cell_scores.rds and 00_data/06.mouse_scores_Mff_pseudobulk.rds
#   exist. This cell does not depend on objects/functions left in memory.
#
# Definition at the single-cell level:
#   CoreScence net = AMS_CoreScence_up - AMS_CoreScence_down
# Higher net values indicate a more CoreScence-like senescence direction.
# This is a continuous transcriptional score, not a senescent-cell cutoff.
#
# Repeats the previous analysis families for the net score:
#   1) Six-age score plots and Kruskal-Wallis / age-Spearman statistics
#   2) Young-vs-Old score plots and Wilcoxon statistics
#   3) Mff-vs-score scatterplots and Spearman statistics
#   4) Six-age and Young-vs-Old difference dotplots
# Analyses are produced at cell level and mouse-summary level.
#=====================================================================

local({
  # ==================================================================
  # 1. SETTINGS
  # ==================================================================
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/6.genescore"
  tissues <- c("Heart", "Limb_Muscle")
  score.name <- "CoreScence_net"
  score.label <- "CoreScence net (up - down)"
  up.column <- "AMS_CoreScence_up"
  down.column <- "AMS_CoreScence_down"
  net.column <- "AMS_CoreScence_net"
  age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
  young.ages <- c("1m", "3m")
  old.ages <- c("18m", "21m", "24m", "30m")
  min.observations.correlation <- 3L
  min.mice.per.group.wilcoxon <- 2L
  min.cells.per.group.wilcoxon <- 2L
  make.young.old.plots <- TRUE
  make.age.dotplots <- TRUE
  save.overviews <- TRUE
  save.pdf <- FALSE

  # Main figure style
  plot.width <- 900L
  plot.height <- 850L
  plot.res <- 160L
  plot.base.size <- 13
  plot.title.size <- 20
  plot.axis.title.size <- 13
  plot.axis.text.size <- 30
  plot.pvalue.size <- 10
  plot.title.wrap.width <- 40
  plot.legend.size <- 13
  plot.violin.width <- 0.8
  plot.box.width <- 0.18
  cell.point.size <- 0.55
  cell.point.alpha <- 0.70
  mouse.point.size <- 3.5
  show.regression.line <- TRUE  # Spearman rho/p is primary; TRUE adds a descriptive OLS line.
  regression.line.color <- "#222222"
  regression.line.width <- 1.0
  scatter.width <- plot.width

  # Dotplot style
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

  # CoreScence net direction: negative/Old lower = black; positive/Old higher = red.
  color.negative <- "#000000"
  color.positive <- "#D95F02"
  color.zero <- "#777777"
  color.mid <- "#F3F3F3"
  age.colors <- c("1m" = "#1B9E77", "3m" = "#D95F02", "18m" = "#7570B3", "21m" = "#E7298A", "24m" = "#66A61E", "30m" = "#E6AB02")
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")

  required.packages <- c("dplyr", "ggplot2", "patchwork", "Cairo")
  missing.packages <- required.packages[!vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing.packages)){stop("Missing R packages: ", paste(missing.packages, collapse = ", "))}
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(!setequal(c(young.ages, old.ages), age.order) || length(intersect(young.ages, old.ages))){stop("Young/Old must partition age.order.")}

  # ==================================================================
  # 2. HELPERS
  # ==================================================================
  write.csv.safe <- function(x, file) {dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x, file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")}
  format.p <- function(p) {
    if(!is.finite(p)){return("p = NA")}
    if(p == 0){return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 1)))}
    paste0("p = ", if(p < 0.001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3))
  }
  bh.reference <- function(p) {out <- rep(NA_real_, length(p)); ok <- is.finite(p); out[ok] <- stats::p.adjust(p[ok], method = "BH"); out}
  safe.test <- function(d, type, unit) {
    n.obs <- nrow(d)
    r <- data.frame(test = type, unit = unit, n_observations = n.obs, n_mice = length(unique(d$mouse.id)), n_age_groups = length(unique(d$age)), n_Young = sum(d$group == "Young"), n_Old = sum(d$group == "Old"), mean_Young = if(any(d$group == "Young")) mean(d$score[d$group == "Young"]) else NA_real_, mean_Old = if(any(d$group == "Old")) mean(d$score[d$group == "Old"]) else NA_real_, mean_difference_Old_minus_Young = NA_real_, rho = NA_real_, statistic = NA_real_, p_value = NA_real_, exact_used = NA, status = "not_tested", reason = "", warning = "", cell_dependence_adjusted = FALSE, age_adjusted = FALSE, stringsAsFactors = FALSE)
    r$mean_difference_Old_minus_Young <- r$mean_Old - r$mean_Young
    r$interpretation <- if(unit == "cell") "Exploratory: same-mouse dependence is not modeled; not mouse-level inference" else "One observation per mouse within tissue/cell type; unadjusted association"
    fail <- function(reason) {r$reason <- reason; r}
    if(!n.obs){return(fail("No observations"))}
    if(any(!is.finite(d$score)) || anyNA(d$group) || anyNA(d$age)){return(fail("Missing/nonfinite score or age/group"))}
    if(unit == "mouse" && anyDuplicated(d$mouse.id)){stop("Repeated mouse in a mouse-level test.")}
    if(type %in% c("Age_Spearman", "Mff_Spearman")) {
      x <- if(type == "Age_Spearman") d$age_months else d$mff
      if(any(!is.finite(x))){return(fail("Missing/nonfinite x values"))}
      if(length(unique(x)) < 2L || length(unique(d$score)) < 2L){return(fail("Constant x or score; correlation undefined"))}
      r$rho <- unname(stats::cor(x, d$score, method = "spearman"))
      if(n.obs < min.observations.correlation){return(fail("Too few observations for correlation test; rho retained"))}
      use.exact <- n.obs < 10L && !anyDuplicated(x) && !anyDuplicated(d$score)
      r$exact_used <- use.exact
      test.call <- function() stats::cor.test(x, d$score, method = "spearman", alternative = "two.sided", exact = use.exact)
    } else if(type == "Age_Kruskal_Wallis") {
      if(length(unique(d$age)) < 2L){return(fail("Only one observed age group"))}
      if(n.obs <= length(unique(d$age))){return(fail("No within-age replication for omnibus comparison"))}
      if(length(unique(d$score)) < 2L){return(fail("All scores identical"))}
      test.call <- function() stats::kruskal.test(d$score, factor(d$age))
    } else if(type == "Young_Old_Wilcoxon") {
      minimum <- if(unit == "mouse") min.mice.per.group.wilcoxon else min.cells.per.group.wilcoxon
      if(r$n_Young < minimum || r$n_Old < minimum){return(fail(paste0("Young and Old each need >= ", minimum, " ", unit, " observations for this test; observations retained in plots")))}
      if(length(unique(d$score)) < 2L){return(fail("All scores identical"))}
      use.exact <- unit == "mouse" && max(r$n_Young, r$n_Old) < 50L && !anyDuplicated(d$score)
      r$exact_used <- use.exact
      test.call <- function() stats::wilcox.test(d$score[d$group == "Young"], d$score[d$group == "Old"], alternative = "two.sided", paired = FALSE, exact = use.exact, correct = !use.exact)
    } else {
      stop("Unknown test type: ", type)
    }
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(test.call(), warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
    r$warning <- paste(unique(warnings), collapse = " | ")
    if(inherits(fit, "error")){return(fail(conditionMessage(fit)))}
    if(length(fit$p.value) != 1L || !is.finite(fit$p.value)){return(fail("Nonfinite p-value"))}
    r$status <- "tested"
    r$reason <- "ok"
    r$p_value <- fit$p.value
    r$statistic <- unname(fit$statistic)
    r
  }
  plot.theme <- function(legend = FALSE) {
    theme_bw(base_size = plot.base.size) + theme(plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.title = element_text(size = plot.axis.title.size), axis.text = element_text(size = plot.axis.text.size, color = "black"), panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none", legend.title = element_text(size = plot.legend.size), legend.text = element_text(size = plot.legend.size), plot.margin = margin(12, 16, 12, 12))
  }
  dotplot.theme <- function(legend = FALSE) {
    theme_bw(base_size = dotplot.base.size) + theme(plot.title = element_text(size = dotplot.title.size, face = "bold", margin = margin(b = 14)), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.title = element_text(size = dotplot.axis.title.size), axis.title.x = element_text(margin = margin(t = 12)), axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"), axis.text.y = element_text(size = dotplot.celltype.text.size, color = "black", margin = margin(r = 10)), panel.grid.minor = element_blank(), legend.position = if(legend) "bottom" else "none", legend.title = element_text(size = dotplot.legend.size), legend.text = element_text(size = dotplot.legend.size), legend.key.width = grid::unit(1.0, "cm"), plot.margin = margin(15, 25, 12, 12))
  }
  make.title <- function(tissue, ct) {paste(paste(strwrap(paste0(gsub("_", " ", tissue), " / ", ct), width = plot.title.wrap.width), collapse = "\n"), score.label, sep = "\n")}
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
    summaries <- d %>% group_by(xgroup) %>% summarise(n = n(), vmin = min(score), vmax = max(score), .groups = "drop")
    variable.groups <- as.character(summaries$xgroup[summaries$n >= 2 & summaries$vmin < summaries$vmax])
    variable <- d[as.character(d$xgroup) %in% variable.groups, , drop = FALSE]
    p <- ggplot(d, aes(x = xgroup, y = score)) + geom_blank()
    if(unit == "cell") {
      if(nrow(variable)){p <- p + geom_violin(data = variable, aes(fill = xgroup), width = plot.violin.width, trim = TRUE, scale = "width", linewidth = 0.4) + geom_boxplot(data = variable, width = plot.box.width, fill = "white", outlier.shape = NA, linewidth = 0.4)}
      special <- summaries[!as.character(summaries$xgroup) %in% variable.groups, , drop = FALSE]
      if(nrow(special)){p <- p + geom_errorbar(data = special, aes(x = xgroup, ymin = vmin, ymax = vmax), inherit.aes = FALSE, width = 0.4) + geom_point(data = special, aes(x = xgroup, y = vmin), inherit.aes = FALSE, size = 2)}
    } else {
      box.groups <- as.character(summaries$xgroup[summaries$n >= 2])
      if(length(box.groups)){p <- p + geom_boxplot(data = d[as.character(d$xgroup) %in% box.groups, , drop = FALSE], width = 0.55, fill = "grey93", outlier.shape = NA)}
      if(nrow(d)){p <- p + geom_point(aes(color = xgroup), position = position_jitter(width = 0.08, height = 0, seed = 1L), size = mouse.point.size)}
    }
    if(unit == "cell" && nrow(variable)){p <- p + scale_fill_manual(values = palette)}
    if(unit == "mouse" && nrow(d)){p <- p + scale_color_manual(values = palette)}
    bounds <- y.bounds(d$score)
    p + scale_x_discrete(limits = levels, labels = if(comparison == "age") sub("m$", "", levels) else levels, drop = FALSE) + annotate("text", x = length(levels) + 0.35, y = bounds[2], label = format.p(pvalue), hjust = 1, vjust = 1, size = plot.pvalue.size, lineheight = 1.0) + coord_cartesian(ylim = bounds) + labs(title = title, x = if(comparison == "age") "Age (months)" else NULL, y = if(unit == "cell") "CoreScence net score" else "Mean CoreScence net score per mouse") + plot.theme()
  }
  correlation.plot <- function(d.full, unit, stat, title, analysis.group = "All", color.by = "age") {
    if(!analysis.group %in% c("All", "Young", "Old")){stop("Invalid scatter group.")}
    bounds <- y.bounds(d.full$score)
    x.bounds <- if(nrow(d.full)) range(d.full$mff) else c(0, 1)
    if(diff(x.bounds) == 0){x.bounds <- x.bounds + c(-0.05, 0.05) * max(1, abs(x.bounds[1]))}
    d <- if(analysis.group == "All") d.full else d.full[d.full$group == analysis.group, , drop = FALSE]
    if(analysis.group != "All"){title <- paste0(title, " / ", analysis.group)}
    stat <- stat[stat$analysis_group == analysis.group, , drop = FALSE]
    if(nrow(stat) != 1L){stop("Expected one Mff correlation result for ", analysis.group)}
    label <- paste0("rho = ", if(is.finite(stat$rho)) sprintf("%.2f", stat$rho) else "NA", "\n", format.p(stat$p_value))
    if(color.by == "age") {
      color.levels <- if(analysis.group == "All") age.order else if(analysis.group == "Young") young.ages else old.ages
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
    p <- ggplot(d, aes(x = mff, y = score)) + geom_point(aes(color = point_color), size = if(unit == "cell") cell.point.size else mouse.point.size, alpha = if(unit == "cell") cell.point.alpha else 0.95)
    if(nrow(d)){p <- p + scale_color_manual(values = palette, breaks = color.levels, labels = color.labels, drop = FALSE, name = legend.title) + guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 3)))}
    if(show.regression.line && nrow(d) >= 2L && length(unique(d$mff)) >= 2L){p <- p + geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, fullrange = FALSE, color = regression.line.color, linewidth = regression.line.width, show.legend = FALSE)}
    p + annotate("text", x = Inf, y = bounds[2], label = label, hjust = 1.05, vjust = 1, size = plot.pvalue.size, lineheight = 1.0) + coord_cartesian(xlim = x.bounds, ylim = bounds) + labs(title = title, x = if(unit == "cell") "Mff expression (RNA LogNormalize)" else "Mff pseudobulk log2(TMM CPM + 1)", y = if(unit == "cell") "CoreScence net score" else "Mean CoreScence net score per mouse") + plot.theme(legend = nrow(d) > 0L)
  }
  age.dotplot <- function(d, unit, title, celltypes) {
    d$age_factor <- factor(d$age, levels = age.order)
    d$celltype_axis <- factor(d$celltype, levels = rev(celltypes))
    observed <- d[is.finite(d$mean_score), , drop = FALSE]
    missing <- d[!is.finite(d$mean_score), , drop = FALSE]
    p <- ggplot(d, aes(x = age_factor, y = celltype_axis)) + geom_blank()
    if(nrow(observed)) {
      limit <- max(abs(observed$mean_score), 0.05)
      p <- p + geom_point(data = observed, aes(color = mean_score), size = dotplot.point.size) + scale_color_gradient2(low = color.negative, mid = color.mid, high = color.positive, midpoint = 0, limits = c(-limit, limit), name = if(unit == "mouse") "Mean net score\n(equal mice)" else "Mean net score\n(equal cells)")
    }
    if(nrow(missing)){p <- p + geom_text(data = missing, label = "NA", color = "grey50", size = dotplot.na.text.size)}
    p + scale_x_discrete(limits = age.order, labels = sub("m$", "", age.order), drop = FALSE) + scale_y_discrete(drop = FALSE) + labs(title = title, x = "Age (months)", y = NULL) + dotplot.theme(legend = nrow(observed) > 0L) + theme(legend.position = if(nrow(observed)) "right" else "none")
  }
  difference.data <- function(statistics) {
    d <- as.data.frame(statistics[statistics$test == "Young_Old_Wilcoxon", , drop = FALSE])
    if(anyDuplicated(d$celltype)){stop("Repeated cell type in score-difference summary.")}
    d$effect_available <- is.finite(d$mean_difference_Old_minus_Young)
    d$effect_reason <- ifelse(d$effect_available, "estimated", "Missing Young or Old group; no score difference estimated")
    d$direction <- ifelse(!d$effect_available, NA_character_, ifelse(d$mean_difference_Old_minus_Young > 0, "Old higher", ifelse(d$mean_difference_Old_minus_Young < 0, "Old lower", "Equal means")))
    d$significance <- ifelse(!d$effect_available, NA_character_, ifelse(d$status != "tested" | !is.finite(d$p_value), "Not tested", ifelse(d$p_value < dotplot.p.cutoff, "Below cutoff", "At or above cutoff")))
    d$effect_definition <- "Arithmetic mean CoreScence net score in Old minus Young; not log2FC"
    d
  }
  difference.dotplot <- function(d, unit, title, celltypes) {
    # Match the original Mff log2FC dotplot design while retaining the correct
    # net-score effect: arithmetic mean score in Old minus Young.
    d <- d[order(!d$effect_available, d$mean_difference_Old_minus_Young, d$celltype, na.last = TRUE), , drop = FALSE]
    d$celltype_axis <- factor(d$celltype, levels = rev(d$celltype))
    estimated <- d[is.finite(d$mean_difference_Old_minus_Young), , drop = FALSE]
    unavailable <- d[!is.finite(d$mean_difference_Old_minus_Young), , drop = FALSE]
    limit <- max(c(abs(estimated$mean_difference_Old_minus_Young), 0.05)) * 1.25
    x.breaks <- pretty(c(-limit, limit), n = 5)
    x.breaks <- sort(unique(c(0, x.breaks[x.breaks >= -limit & x.breaks <= limit])))
    direction.colors <- c("Old lower" = color.negative, "Old higher" = color.positive, "Equal means" = color.zero)
    significance.labels <- c("Below cutoff" = paste0("Unadjusted p < ", dotplot.p.cutoff), "At or above cutoff" = paste0("Unadjusted p >= ", dotplot.p.cutoff), "Not tested" = "p = NA")
    display.title <- paste0(gsub("_", " ", unique(d$tissue)[1]), " | ", score.label)
    x.label <- if(unit == "cell") "Cell-mean CoreScence net-score difference (Old - Young)" else "Mouse-mean CoreScence net-score difference (Old - Young)"
    p <- ggplot(d, aes(y = celltype_axis)) + geom_vline(xintercept = 0, color = "#999999", linewidth = 0.6, linetype = "dashed") + geom_point(data = estimated, aes(x = mean_difference_Old_minus_Young, color = direction, shape = significance), size = dotplot.point.size, stroke = dotplot.point.stroke) + scale_color_manual(values = direction.colors, guide = "none") + scale_shape_manual(values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4), limits = names(significance.labels), labels = unname(significance.labels), drop = FALSE, name = NULL) + scale_x_continuous(breaks = x.breaks, limits = c(-limit, limit), expand = expansion(mult = 0.02)) + scale_y_discrete(drop = FALSE, expand = expansion(add = 0.7)) + labs(title = display.title, subtitle = NULL, x = x.label, y = NULL) + theme_classic(base_size = dotplot.base.size, base_family = "sans") + theme(plot.title = element_text(size = dotplot.title.size, face = "bold", margin = margin(b = 14)), axis.title.x = element_text(size = dotplot.axis.title.size, margin = margin(t = 12)), axis.text.x = element_text(size = dotplot.axis.number.size, color = "black"), axis.text.y = element_text(size = dotplot.celltype.text.size, color = "black", margin = margin(r = 10)), axis.ticks.y = element_blank(), panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.35), legend.position = "bottom", legend.text = element_text(size = dotplot.legend.size), legend.key.width = grid::unit(1.0, "cm"), plot.margin = margin(15, 25, 12, 12)) + guides(shape = guide_legend(nrow = 1, override.aes = list(color = "#444444", size = 4)))
    if(nrow(unavailable)){p <- p + geom_text(data = unavailable, aes(y = celltype_axis), x = limit * 0.96, label = "NA", inherit.aes = FALSE, hjust = 1, color = "#888888", size = dotplot.na.text.size)}
    p
  }
  save.plot <- function(p, stem, width = plot.width, height = plot.height, allow.pdf = save.pdf, res = plot.res) {
    dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
    Cairo::CairoPNG(filename = paste0(stem, ".png"), width = width, height = height, res = res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    if(allow.pdf){ggplot2::ggsave(paste0(stem, ".pdf"), plot = p, device = grDevices::cairo_pdf, width = width / res, height = height / res, units = "in", limitsize = FALSE)}
  }
  save.overview <- function(plots, folder, label, panel.width = plot.width) {
    if(!save.overviews || !length(plots)){return(invisible(NULL))}
    ncols <- min(3L, length(plots))
    nrows <- ceiling(length(plots) / ncols)
    p <- patchwork::wrap_plots(plots, ncol = ncols)
    save.plot(p, file.path(folder, "overview", paste0(label, ".all_celltypes")), width = panel.width * ncols, height = plot.height * nrows, allow.pdf = FALSE)
  }

  # ==================================================================
  # 3. LOAD SAVED SCORES, CREATE NET SCORE, REPEAT ANALYSES
  # ==================================================================
  all.statistics <- list()
  for(tissue.name in tissues) {
    cat("\n========== CoreScence net: ", tissue.name, " ==========\n", sep = "")
    source.data.dir <- file.path(output.base, tissue.name, "00_data")
    cell.file <- file.path(source.data.dir, "04.cell_scores.rds")
    mouse.file <- file.path(source.data.dir, "06.mouse_scores_Mff_pseudobulk.rds")
    coverage.file <- file.path(source.data.dir, "02.geneset_coverage.csv")
    missing.files <- c(cell.file, mouse.file, coverage.file)[!file.exists(c(cell.file, mouse.file, coverage.file))]
    if(length(missing.files)){stop("Run the original five-signature analysis first. Missing:\n", paste(missing.files, collapse = "\n"))}

    cells <- as.data.frame(readRDS(cell.file))
    mice <- as.data.frame(readRDS(mouse.file))
    coverage <- utils::read.csv(coverage.file, stringsAsFactors = FALSE, check.names = FALSE)
    cell.required <- c("cell", "tissue", "celltype", "mouse.id", "age", "age_months", "group", "Mff_LogNormalize", up.column, down.column)
    mouse.required <- c("tissue", "celltype", "mouse.id", "age", "age_months", "group", "eligible", "Mff_log2_TMM_CPM_plus1", up.column, down.column)
    if(!all(cell.required %in% names(cells))){stop("Missing cell columns: ", paste(setdiff(cell.required, names(cells)), collapse = ", "))}
    if(!all(mouse.required %in% names(mice))){stop("Missing mouse columns: ", paste(setdiff(mouse.required, names(mice)), collapse = ", "))}
    if(!all(c("geneset", "n_used") %in% names(coverage)) || !all(c("CoreScence_up", "CoreScence_down") %in% coverage$geneset)){stop("CoreScence coverage rows are missing: ", coverage.file)}
    if(anyNA(cells[, cell.required]) || any(!is.finite(cells[[up.column]])) || any(!is.finite(cells[[down.column]])) || any(!is.finite(cells$Mff_LogNormalize))){stop("Missing/nonfinite cell-level net-score input in ", tissue.name)}
    if(anyNA(mice[, setdiff(mouse.required, "Mff_log2_TMM_CPM_plus1")]) || any(!is.finite(mice[[up.column]])) || any(!is.finite(mice[[down.column]]))){stop("Missing/nonfinite mouse net-score input in ", tissue.name)}
    if(any(!cells$age %in% age.order) || any(!mice$age %in% age.order)){stop("Unexpected ages in saved score tables.")}

    cells[[net.column]] <- cells[[up.column]] - cells[[down.column]]
    mice[[net.column]] <- mice[[up.column]] - mice[[down.column]]
    cells$net_score_definition <- paste0(up.column, " - ", down.column)
    mice$net_score_definition <- paste0("Mean cell ", up.column, " - mean cell ", down.column)

    # Verify that each mouse score equals the mean of its own cells' net scores.
    expected.mouse <- cells %>% group_by(celltype, mouse.id) %>% summarise(expected_net = mean(.data[[net.column]]), .groups = "drop")
    check.mouse <- mice %>% select(celltype, mouse.id, all_of(net.column)) %>% left_join(expected.mouse, by = c("celltype", "mouse.id"))
    if(nrow(check.mouse) != nrow(mice) || anyNA(check.mouse$expected_net) || max(abs(check.mouse[[net.column]] - check.mouse$expected_net)) > 1e-10){stop("Mouse net score does not equal the mean of its cells in ", tissue.name)}

    score.dir <- file.path(output.base, tissue.name, score.name)
    net.data.dir <- file.path(score.dir, "00_data")
    dir.create(net.data.dir, recursive = TRUE, showWarnings = FALSE)
    write.csv.safe(cells, file.path(net.data.dir, "cell_CoreScence_net_scores.csv"))
    saveRDS(cells, file.path(net.data.dir, "cell_CoreScence_net_scores.rds"), compress = "gzip")
    write.csv.safe(mice, file.path(net.data.dir, "mouse_CoreScence_net_scores.csv"))
    saveRDS(mice, file.path(net.data.dir, "mouse_CoreScence_net_scores.rds"), compress = "gzip")

    celltypes <- sort(unique(cells$celltype))
    file.map <- data.frame(celltype = celltypes, file_stem = sprintf("%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 90)), stringsAsFactors = FALSE)
    write.csv.safe(file.map, file.path(net.data.dir, "celltype_file_map.csv"))
    tissue.stats <- list()

    for(unit in c("mouse", "cell")) {
      source <- if(unit == "mouse") mice[mice$eligible, , drop = FALSE] else cells
      d.all <- source[, c("tissue", "celltype", "mouse.id", "age", "age_months", "group"), drop = FALSE]
      d.all$score <- source[[net.column]]
      d.all$mff <- if(unit == "mouse") source$Mff_log2_TMM_CPM_plus1 else source$Mff_LogNormalize
      if(any(!is.finite(d.all$score)) || any(!is.finite(d.all$mff))){stop("Nonfinite plotting/testing values: ", tissue.name, " / ", unit)}

      age.folder <- file.path(score.dir, if(unit == "mouse") "01_mouse_age_score" else "03_cell_age_score")
      mff.folder <- file.path(score.dir, if(unit == "mouse") "02_mouse_Mff_score" else "04_cell_Mff_score")
      dir.create(age.folder, recursive = TRUE, showWarnings = FALSE)
      dir.create(mff.folder, recursive = TRUE, showWarnings = FALSE)
      unit.stats <- list()
      six.plots <- list()
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
        st$geneset <- score.name
        st$n_genes_used_up <- coverage$n_used[coverage$geneset == "CoreScence_up"][[1]]
        st$n_genes_used_down <- coverage$n_used[coverage$geneset == "CoreScence_down"][[1]]
        st$score_definition <- paste0(up.column, " - ", down.column)
        unit.stats[[ct]] <- st
        summaries[[ct]] <- d %>% group_by(tissue, celltype, age, age_months, group) %>% summarise(n_observations = n(), n_mice = n_distinct(mouse.id), mean_score = mean(score), median_score = median(score), sd_score = if(n() > 1L) sd(score) else NA_real_, min_score = min(score), max_score = max(score), .groups = "drop")

        title <- make.title(tissue.name, ct)
        six.plots[[ct]] <- age.plot(d, unit, "age", st$p_value[st$test == "Age_Kruskal_Wallis"], title)
        mff.stats <- st[grepl("^Mff_Spearman", st$test), , drop = FALSE]
        corr.plots[[ct]] <- correlation.plot(d, unit, mff.stats, title, analysis.group = "All", color.by = "age")
        corr.groupcolor.plots[[ct]] <- correlation.plot(d, unit, mff.stats, title, analysis.group = "All", color.by = "group")
        corr.young.plots[[ct]] <- correlation.plot(d, unit, mff.stats, title, analysis.group = "Young", color.by = "age")
        corr.old.plots[[ct]] <- correlation.plot(d, unit, mff.stats, title, analysis.group = "Old", color.by = "age")

        stem <- file.map$file_stem[i]
        save.plot(six.plots[[ct]], file.path(age.folder, "six_ages", paste0(stem, ".age_score")))
        save.plot(corr.plots[[ct]], file.path(mff.folder, "scatter", paste0(stem, ".Mff_score")), width = scatter.width)
        save.plot(corr.groupcolor.plots[[ct]], file.path(mff.folder, "scatter", paste0(stem, ".Mff_score.by_YoungOld")), width = scatter.width)
        if(make.young.old.plots) {
          yo.plots[[ct]] <- age.plot(d, unit, "group", st$p_value[st$test == "Young_Old_Wilcoxon"], title)
          save.plot(yo.plots[[ct]], file.path(age.folder, "young_old", paste0(stem, ".Young_Old_score")))
        }
      }

      statistics <- do.call(rbind, unit.stats)
      rownames(statistics) <- NULL
      tissue.stats[[unit]] <- statistics
      age.summary <- bind_rows(summaries) %>% arrange(celltype, age_months)
      write.csv.safe(age.summary, file.path(age.folder, "age_score_summary.csv"))

      if(make.age.dotplots) {
        dot.folder <- file.path(age.folder, "dotplot")
        dir.create(dot.folder, recursive = TRUE, showWarnings = FALSE)
        dot.data <- expand.grid(celltype = celltypes, age = age.order, stringsAsFactors = FALSE) %>% left_join(age.summary, by = c("celltype", "age"))
        dot.data$tissue <- tissue.name
        dot.data$geneset <- score.name
        dot.data$unit <- unit
        dot.data$age_months <- as.numeric(sub("m$", "", dot.data$age))
        dot.data$group <- ifelse(dot.data$age %in% young.ages, "Young", "Old")
        dot.data$n_observations[is.na(dot.data$n_observations)] <- 0L
        dot.data$n_mice[is.na(dot.data$n_mice)] <- 0L
        dot.data$summary_weighting <- if(unit == "mouse") "Equal weight per eligible mouse" else "Equal weight per cell; not balanced by mouse"
        dot.data$estimate_status <- ifelse(is.finite(dot.data$mean_score), "estimated", "no_observations")
        dot.data <- dot.data %>% arrange(celltype, age_months)
        write.csv.safe(dot.data, file.path(dot.folder, "age_mean_score.plot_data.csv"))
        dot.title <- paste0(gsub("_", " ", tissue.name), " / ", unit, "\n", score.label)
        dot.width <- round(dotplot.width.in * dotplot.res)
        dot.height <- round(max(dotplot.height.minimum.in, 2.0 + length(celltypes) * dotplot.height.per.celltype.in) * dotplot.res)
        save.plot(age.dotplot(dot.data, unit, dot.title, celltypes), file.path(dot.folder, "age_mean_score.all_celltypes"), width = dot.width, height = dot.height, res = dotplot.res)
        save.plot(difference.dotplot(difference.data(statistics), unit, dot.title, celltypes), file.path(dot.folder, "Young_Old_score_difference.all_celltypes"), width = dot.width, height = dot.height, res = dotplot.res)
      }

      save.overview(six.plots, age.folder, "six_ages")
      save.overview(yo.plots, age.folder, "young_old")
      save.overview(corr.plots, mff.folder, "Mff_score", panel.width = scatter.width)
      save.overview(corr.groupcolor.plots, mff.folder, "Mff_score.by_YoungOld", panel.width = scatter.width)
      save.overview(corr.young.plots, mff.folder, "Mff_score.Young", panel.width = scatter.width)
      save.overview(corr.old.plots, mff.folder, "Mff_score.Old", panel.width = scatter.width)
    }

    statistics <- bind_rows(tissue.stats) %>% group_by(tissue, unit, test) %>% mutate(p_adj_BH = bh.reference(p_value), n_tests_in_BH_family = sum(is.finite(p_value)), p_adjust_scope = "Within CoreScence net score x tissue x unit x test, across cell types") %>% ungroup()
    write.csv.safe(statistics, file.path(score.dir, "all_statistics.csv"))
    for(unit in c("mouse", "cell")) {
      s <- statistics[statistics$unit == unit, , drop = FALSE]
      age.folder <- file.path(score.dir, if(unit == "mouse") "01_mouse_age_score" else "03_cell_age_score")
      mff.folder <- file.path(score.dir, if(unit == "mouse") "02_mouse_Mff_score" else "04_cell_Mff_score")
      is.mff <- grepl("^Mff_Spearman", s$test)
      write.csv.safe(s[!is.mff, , drop = FALSE], file.path(age.folder, "statistics.csv"))
      write.csv.safe(s[is.mff, , drop = FALSE], file.path(mff.folder, "statistics.csv"))
      if(make.age.dotplots){write.csv.safe(difference.data(s), file.path(age.folder, "dotplot", "Young_Old_score_difference.plot_data.csv"))}
    }

    settings <- c(paste0("Completed: ", Sys.time()), paste0("Tissue: ", tissue.name), paste0("Cell input: ", cell.file), paste0("Mouse input: ", mouse.file), paste0("Output: ", score.dir), paste0("Score definition: ", net.column, " = ", up.column, " - ", down.column), "The net score is computed per cell from two previously calculated AddModuleScores", "Higher net values indicate up-program activation and/or down-program suppression", "The net score is continuous and has no senescent-cell classification cutoff", "Up and down components remain available in the saved tables for diagnostic interpretation", "Mouse net score equals the arithmetic mean of the same cells' net scores; equality was verified", "No AddModuleScore, RNA normalization, PCA, UMAP, clustering or pseudobulk count aggregation is rerun", "Six-age p: cell/mouse Kruskal-Wallis across observed age groups", "Young/Old p: two-sided unpaired Wilcoxon rank-sum", "Mff association: two-sided Spearman rho/raw p", "Cell-level tests do not adjust for same-mouse cell dependence and are exploratory", "Mouse-level tests use one eligible mouse observation per tissue/cell type", "Associations are unadjusted for age, sex, batch and RNA depth", "Raw p values are displayed; BH reference is saved across cell types within net score x tissue x unit x test", paste0("Regression line shown: ", show.regression.line), "Difference dotplot: red=Old lower, black=Old higher, grey=equal means", "Age mean dotplot: red=negative net score, white=zero, black=positive net score", paste0("Dotplot filled point: raw Wilcoxon p < ", dotplot.p.cutoff, "; open point: p >= cutoff; cross: p=NA"), "CoreScence net is not the DeepScence model and does not inherit a published DeepScence cutoff/performance claim", "", capture.output(utils::sessionInfo()))
    settings <- settings[!settings %in% c("Difference dotplot: red=Old lower, black=Old higher, grey=equal means", "Age mean dotplot: red=negative net score, white=zero, black=positive net score")]
    settings <- c(settings, "Difference dotplot: black=Old lower, red=Old higher, grey=equal means", "Age mean dotplot: black=negative net score, white=zero, red=positive net score")
    writeLines(settings, file.path(score.dir, "analysis_settings_and_sessionInfo.txt"), useBytes = TRUE)
    all.statistics[[tissue.name]] <- statistics
    cat("Saved: ", score.dir, "\n", sep = "")
    rm(cells, mice, coverage, source, d.all, d, tissue.stats, statistics, six.plots, yo.plots, corr.plots, corr.groupcolor.plots, corr.young.plots, corr.old.plots)
    invisible(gc())
  }

  combined <- bind_rows(all.statistics)
  write.csv.safe(combined, file.path(output.base, "CoreScence_net.all_statistics.csv"))
  cat("\nCompleted CoreScence net analysis for both tissues.\n")
  cat("Combined statistics: ", file.path(output.base, "CoreScence_net.all_statistics.csv"), "\n", sep = "")
  invisible(list(output = output.base, statistics = combined))
})
