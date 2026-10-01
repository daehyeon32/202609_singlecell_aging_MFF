# Aging-MFF: one Mff-score rho heatmap per tissue.
#
# Output from one run:
#   1) Heart.Mff_score_rho_heatmap.rho_m0.25_p0.25.png
#   2) Limb_Muscle.Mff_score_rho_heatmap.rho_m0.25_p0.25.png
#
# x = score, y = cell type, fill = Spearman rho.
# No numbers, p-values, or symbols are printed inside the heatmap boxes.
# Raw Spearman p_value <= 0.05: positive rho = black border; negative rho = red border.
# rho == 0, missing/nonfinite values, and untested results receive no colored border.
# geom_tile reference: https://ggplot2.tidyverse.org/reference/geom_tile.html

run.integrated.mff.score.rho.heatmap <- function() {
  # ==================================================================
  # 1. SETTINGS: adjust analysis choice, fonts, size, and colors here
  # ==================================================================
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/6.genescore"
  individual.statistics.file <- file.path(output.base, "02.all_statistics.csv")
  net.statistics.file <- file.path(output.base, "CoreScence_net.all_statistics.csv")
  output.dir <- file.path(output.base, "05.Mff_score_rho_heatmap_rho_m0.25_p0.25")

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
  font.score.size <- 10
  font.celltype.size <- 15
  font.legend.title.size <- 15
  font.legend.text.size <- 13

  # Figure dimensions.
  plot.width.in <- 12
  plot.height.minimum.in <- 6
  plot.height.per.celltype.in <- 0.48
  png.dpi <- 300
  save.pdf <- FALSE  # TRUE also saves one PDF per tissue.

  # rho is signed, so zero is kept white.
  # Keep the supplied fill colors: negative = bright red; positive = dark red.
  rho.negative.color <- "#EF3B2C"
  rho.zero.color <- "#FFFFFF"
  rho.positive.color <- "#67000D"
  rho.na.color <- "#E5E5E5"
  legend.bar.height.cm <- 8
  legend.bar.width.cm <- 0.7

  rho.minimum <- -0.25
  rho.maximum <- 0.25
  rho.breaks <- c(-0.25, -0.125, 0, 0.125, 0.25)

  # 유의한 칸의 테두리 설정: CSV의 raw p_value를 사용합니다.
  # p_adj_BH를 사용하거나 p-value를 다시 계산하지 않습니다.
  significance.p.cutoff <- 0.05
  border.positive.color <- "black"  # p <= 0.05 AND rho > 0
  border.negative.color <- "red"    # p <= 0.05 AND rho < 0
  border.linewidth <- 1.4           # 테두리 두께
  border.tile.size <- 0.92          # 칸 안쪽에 그려 이웃 테두리의 겹침 방지

  required.packages <- c("dplyr", "ggplot2", "scales")
  missing.packages <- required.packages[!vapply(required.packages, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing.packages)){stop("Missing R packages: ", paste(missing.packages, collapse = ", "))}
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(!analysis.unit %in% c("cell", "mouse")){stop("analysis.unit must be 'cell' or 'mouse'.")}
  if(!analysis.group %in% c("All", "Young", "Old")){stop("analysis.group must be 'All', 'Young', or 'Old'.")}
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # ==================================================================
  # 2. LOAD INDIVIDUAL SCORES AND CORESCENCE NET
  # ==================================================================
  missing.files <- c(individual.statistics.file, net.statistics.file)[!file.exists(c(individual.statistics.file, net.statistics.file))]
  if(length(missing.files)){stop("Run both score-analysis scripts first. Missing:\n", paste(missing.files, collapse = "\n"))}

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

    # 유효한 검정 결과 중 raw p <= 0.05인 칸만 선택합니다.
    # 표시용으로 반올림한 값이 아니라 CSV의 원래 p_value/rho를 사용합니다.
    significant <- d %>%
      filter(status == "tested", is.finite(p_value), p_value <= significance.p.cutoff, is.finite(rho))
    positive <- significant %>% filter(rho > 0)
    negative <- significant %>% filter(rho < 0)

    p <- ggplot(d, aes(x = score.axis, y = celltype.axis, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.8) +
      # fill = NA: 기존 rho 색상을 유지하면서 테두리만 추가합니다.
      geom_tile(
        data = positive,
        aes(x = score.axis, y = celltype.axis),
        inherit.aes = FALSE, fill = NA,
        color = border.positive.color, linewidth = border.linewidth,
        width = border.tile.size, height = border.tile.size,
        show.legend = FALSE
      ) +
      geom_tile(
        data = negative,
        aes(x = score.axis, y = celltype.axis),
        inherit.aes = FALSE, fill = NA,
        color = border.negative.color, linewidth = border.linewidth,
        width = border.tile.size, height = border.tile.size,
        show.legend = FALSE
      ) +
      scale_fill_gradient2(
        low = rho.negative.color, mid = rho.zero.color, high = rho.positive.color,
        midpoint = 0, limits = c(rho.minimum, rho.maximum), breaks = rho.breaks,
        oob = scales::squish, na.value = rho.na.color, name = "Spearman\nrho",
        guide = guide_colorbar(
          barheight = grid::unit(legend.bar.height.cm, "cm"),
          barwidth = grid::unit(legend.bar.width.cm, "cm"),
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      coord_fixed(ratio = 0.65) +
      labs(title = unname(tissue.labels[tissue.name]), x = "Score", y = "Cell type") +
      theme_bw(base_size = font.base.size, base_family = font.family) +
      theme(
        plot.title = element_text(size = font.title.size, face = "bold", hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(size = font.axis.title.size),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.text.x = element_text(size = font.score.size, color = "black", angle = 35, hjust = 1, vjust = 1, lineheight = 0.9),
        axis.text.y = element_text(size = font.celltype.size, color = "black"),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = font.legend.title.size),
        legend.text = element_text(size = font.legend.text.size),
        plot.margin = margin(15, 20, 15, 15)
      )
    p
  }

  saved.plots <- list()
  for(tissue.name in tissue.order) {
    p <- make.tissue.heatmap(combined, tissue.name)
    if(is.null(p)){next}
    plot.height.in <- max(plot.height.minimum.in, 2.7 + length(unique(combined$celltype[combined$tissue == tissue.name])) * plot.height.per.celltype.in)
    stem <- file.path(output.dir, paste0(tissue.name, ".Mff_score_rho_heatmap.rho_m0.25_p0.25"))
    ggsave(paste0(stem, ".png"), plot = p, width = plot.width.in, height = plot.height.in, units = "in", dpi = png.dpi, bg = "white", limitsize = FALSE)
    if(save.pdf){ggsave(paste0(stem, ".pdf"), plot = p, device = grDevices::cairo_pdf, width = plot.width.in, height = plot.height.in, units = "in", bg = "white", limitsize = FALSE)}
    saved.plots[[tissue.name]] <- p
  }

  message("Saved ", length(saved.plots), " tissue heatmaps: ", output.dir)
  invisible(list(output = output.dir, data = combined, plots = saved.plots))
}

run.integrated.mff.score.rho.heatmap()
