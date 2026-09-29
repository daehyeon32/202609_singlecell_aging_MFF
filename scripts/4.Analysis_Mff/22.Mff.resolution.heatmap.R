# Aging-MFF: one Mff-score rho heatmap per tissue.
#
# Input statistics must be generated after excluding 18m and 21m.
# Included ages: 1m, 3m, 24m, 30m.
#
# Output from one run:
#   1) Heart.Mff_score_rho_heatmap.png
#   2) Limb_Muscle.Mff_score_rho_heatmap.png
#
# x = score, y = cell type, fill = Spearman rho.
# No numbers, p-values, or symbols are printed inside the heatmap boxes.

run.integrated.mff.score.rho.heatmap <- function() {
  # ==================================================================
  # 1. SETTINGS: adjust analysis choice, fonts, size, and colors here
  # ==================================================================
  # This folder is produced by both score scripts using ages 1m, 3m, 24m, 30m.
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/9.genescore.132430"
  individual.statistics.file <- file.path(output.base, "02.all_statistics.csv")
  net.statistics.file <- file.path(output.base, "CoreScence_net.all_statistics.csv")
  output.dir <- file.path(output.base, "05.Mff_score_rho_heatmap_two_tissues")

  analysis.unit <- "cell"  # "cell" or "mouse"
  analysis.group <- "All"  # "All", "Young", or "Old"
  tissue.order <- c("Heart", "Limb_Muscle")
  score.order <- c("Hallmark_IFN_alpha", "GO_IFN_beta", "SaulSenMayo", "CoreScence_up", "CoreScence_net")
  score.labels <- c(Hallmark_IFN_alpha = "Hallmark\nIFN-alpha", GO_IFN_beta = "GO\nIFN-beta", SaulSenMayo = "SaulSenMayo", CoreScence_up = "CoreScence\nup", CoreScence_down = "CoreScence\ndown", CoreScence_net = "CoreScence\nnet")
  tissue.labels <- c(Heart = "Heart", Limb_Muscle = "Limb Muscle")

  # Font sizes: change only these values after checking the first result.
  font.family <- "sans"
  font.base.size <- 15
  font.title.size <- 24
  font.axis.title.size <- 19
  font.score.size <- 12
  font.celltype.size <- 15
  font.legend.title.size <- 15
  font.legend.text.size <- 13

  # Figure dimensions.
  plot.width.in <- 14
  plot.height.minimum.in <- 6
  plot.height.per.celltype.in <- 0.48
  png.dpi <- 300
  save.pdf <- FALSE  # TRUE also saves one PDF per tissue.

  # rho is signed, so zero is kept white.
  # Both ends use red-family colors: dark red = negative, bright red = positive.
  rho.negative.color <- "#EF3B2C"
  rho.zero.color <- "#FFFFFF"
  rho.positive.color <- "#67000D"
  rho.na.color <- "#E5E5E5"
  legend.bar.height.cm <- 8
  legend.bar.width.cm <- 0.7

  required.packages <- c("dplyr", "ggplot2")
  missing.packages <- required.packages[!vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing.packages)){stop("Missing R packages: ", paste(missing.packages, collapse = ", "))}
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(!analysis.unit %in% c("cell", "mouse")){stop("analysis.unit must be 'cell' or 'mouse'.")}
  if(!analysis.group %in% c("All", "Young", "Old")){stop("analysis.group must be 'All', 'Young', or 'Old'.")}
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # ==================================================================
  # 2. LOAD THE FIVE INDIVIDUAL SCORES AND CORESCENCE NET
  # ==================================================================
  missing.files <- c(individual.statistics.file, net.statistics.file)[!file.exists(c(individual.statistics.file, net.statistics.file))]
  if(length(missing.files)){stop("Run both score-analysis scripts with ages 1m, 3m, 24m, 30m first. Missing:\n", paste(missing.files, collapse = "\n"))}

  individual <- utils::read.csv(individual.statistics.file, stringsAsFactors = FALSE, check.names = FALSE)
  net <- utils::read.csv(net.statistics.file, stringsAsFactors = FALSE, check.names = FALSE)
  required.columns <- c("tissue", "celltype", "geneset", "unit", "test", "analysis_group", "rho", "p_value", "status", "reason", "n_observations", "n_mice")
  if(!all(required.columns %in% names(individual))){stop("Missing columns in 02.all_statistics.csv: ", paste(setdiff(required.columns, names(individual)), collapse = ", "))}
  if(!all(required.columns %in% names(net))){stop("Missing columns in CoreScence_net.all_statistics.csv: ", paste(setdiff(required.columns, names(net)), collapse = ", "))}

  individual <- individual[individual$geneset %in% score.order[score.order != "CoreScence_net"], required.columns, drop = FALSE]
  net <- net[net$geneset == "CoreScence_net", required.columns, drop = FALSE]
  combined <- bind_rows(individual, net)
  combined <- combined[grepl("^Mff_Spearman", combined$test) & combined$unit == analysis.unit & combined$analysis_group == analysis.group & combined$tissue %in% tissue.order & combined$geneset %in% score.order, , drop = FALSE]
  if(!nrow(combined)){stop("No matching Mff_Spearman results were found for unit=", analysis.unit, " and group=", analysis.group, ".")}

  key.columns <- c("tissue", "celltype", "geneset")
  duplicated.keys <- duplicated(combined[, key.columns]) | duplicated(combined[, key.columns], fromLast = TRUE)
  if(any(duplicated.keys)){stop("Duplicated tissue x celltype x score correlation results were found.")}

  utils::write.csv(combined, file.path(output.dir, paste0("rho_heatmap_plot_data.", analysis.unit, ".", analysis.group, ".csv")), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

  # ==================================================================
  # 3. CREATE AND SAVE ONE HEATMAP FOR EACH TISSUE
  # ==================================================================
  make.tissue.heatmap <- function(data, tissue.name) {
    d <- data[data$tissue == tissue.name, , drop = FALSE]
    if(!nrow(d)){warning("No results for tissue: ", tissue.name, call. = FALSE); return(NULL)}

    complete.grid <- expand.grid(celltype = sort(unique(d$celltype)), geneset = score.order, stringsAsFactors = FALSE)
    d <- complete.grid %>% left_join(d, by = c("celltype", "geneset"))
    d$score.axis <- factor(d$geneset, levels = score.order, labels = unname(score.labels[score.order]))
    d$celltype.axis <- factor(d$celltype, levels = rev(sort(unique(d$celltype))))

    p <- ggplot(d, aes(x = score.axis, y = celltype.axis, fill = rho)) + geom_tile(color = "white", linewidth = 0.8) + scale_fill_gradient2(low = rho.negative.color, mid = rho.zero.color, high = rho.positive.color, midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), na.value = rho.na.color, name = "Spearman\nrho", guide = guide_colorbar(barheight = grid::unit(legend.bar.height.cm, "cm"), barwidth = grid::unit(legend.bar.width.cm, "cm"), frame.colour = "black", ticks.colour = "black")) + scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) + coord_fixed(ratio = 0.6) + labs(title = unname(tissue.labels[tissue.name]), x = "Score", y = "Cell type") + theme_bw(base_size = font.base.size, base_family = font.family) + theme(plot.title = element_text(size = font.title.size, face = "bold", hjust = 0.5, margin = margin(b = 15)), axis.title = element_text(size = font.axis.title.size), axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12)), axis.text.x = element_text(size = font.score.size, color = "black", angle = 35, hjust = 1, vjust = 1, lineheight = 0.9), axis.text.y = element_text(size = font.celltype.size, color = "black"), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right", legend.title = element_text(size = font.legend.title.size), legend.text = element_text(size = font.legend.text.size), plot.margin = margin(15, 20, 15, 15))
    p
  }

  saved.plots <- list()
  for(tissue.name in tissue.order) {
    p <- make.tissue.heatmap(combined, tissue.name)
    if(is.null(p)){next}
    plot.height.in <- max(plot.height.minimum.in, 2.7 + length(unique(combined$celltype[combined$tissue == tissue.name])) * plot.height.per.celltype.in)
    stem <- file.path(output.dir, paste0(tissue.name, ".Mff_score_rho_heatmap"))
    ggsave(paste0(stem, ".png"), plot = p, width = plot.width.in, height = plot.height.in, units = "in", dpi = png.dpi, bg = "white", limitsize = FALSE)
    if(save.pdf){ggsave(paste0(stem, ".pdf"), plot = p, device = grDevices::cairo_pdf, width = plot.width.in, height = plot.height.in, units = "in", bg = "white", limitsize = FALSE)}
    saved.plots[[tissue.name]] <- p
  }

  message("Saved two tissue heatmaps: ", output.dir)
  invisible(list(output = output.dir, data = combined, plots = saved.plots))
}

run.integrated.mff.score.rho.heatmap()
