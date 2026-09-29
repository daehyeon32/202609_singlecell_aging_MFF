# Mff expression on saved tissue-specific UMAPs: ages 1, 3, 24, 30 months.
# Run the WHOLE file with source(). Requires outputs of 3.Umap.tissue.R.
# Existing normalization and UMAP coordinates are reused, not recomputed.
# References:
# https://satijalab.org/seurat/reference/featureplot
# https://satijalab.github.io/seurat-object/reference/Layers.html

run.mff.umap <- function() {
  # --------------------------- SETTINGS ---------------------------
  input.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split"
  output.dir <- file.path(input.dir, "Mff_expression_UMAP.132430")
  tissues <- c("Heart", "Limb_Muscle")
  gene <- "Mff"
  assay <- "RNA"
  reduction <- "umap"
  age.column <- "age"
  included.ages <- c("1m", "3m", "24m", "30m")
  group.ages <- list(
    Young_1m_3m = c("1m", "3m"), Old_24m_30m = c("24m", "30m"),
    Age_3m = "3m", Age_30m = "30m"
  )
  group.labels <- c(
    Young_1m_3m = "Young (1 + 3 months)", Old_24m_30m = "Old (24 + 30 months)",
    Age_3m = "3 months", Age_30m = "30 months"
  )
  # Both limits are computed ONCE over unique selected cells from BOTH tissues.
  # Q95 is a visualization choice, not a biological threshold or universal standard.
  robust.quantile <- 0.95
  point.size <- 0.30
  # Full opacity preserves the expression-to-color mapping.
  expression.colors <- c("#D9D9D9", "#FEE8C8", "#FDBB84", "#E34A33", "#990000")
  font.base <- 16
  font.title <- 20
  font.axis <- 16
  font.legend <- 13
  width.in <- 8
  height.in <- 7
  dpi <- 300
  save.pdf <- FALSE

  required <- c("SeuratObject", "ggplot2", "patchwork", "scales")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  if(utils::packageVersion("SeuratObject") < "5.0.0") {
    stop("SeuratObject >= 5.0.0 is required (supports Assay and Assay5 inputs).")
  }
  suppressPackageStartupMessages(library(ggplot2))
  if(length(robust.quantile) != 1L || !is.finite(robust.quantile) ||
     robust.quantile <= 0 || robust.quantile > 1) stop("Invalid robust.quantile.")
  if(!all(vapply(group.ages, function(a) all(a %in% included.ages), logical(1)))) {
    stop("Every plotted group must be inside included.ages.")
  }
  rds.files <- setNames(
    file.path(input.dir, tissues,
              paste0("Aging.MFF.", tissues, ".seurat.normalization.pca.umap.RDS")),
    tissues
  )
  if(any(!file.exists(rds.files))) {
    stop("Run 3.Umap.tissue.R first. Missing RDS:\n",
         paste(rds.files[!file.exists(rds.files)], collapse = "\n"))
  }

  # Read only the gene row from each normalized data layer.
  # Never silently use counts/scale.data, drop a layer, or replace missing values with zero.
  read.expression <- function(x, wanted.cells) {
    layers <- SeuratObject::Layers(x[[assay]], search = "^data($|\\.)")
    if(!length(layers)) stop("No normalized RNA data layer.")
    parts <- lapply(layers, function(layer) {
      m <- SeuratObject::LayerData(x[[assay]], layer = layer, features = gene)
      if(nrow(m) != 1L || !identical(rownames(m), gene)) {
        stop("Gene not available in normalized layer: ", layer)
      }
      setNames(as.numeric(m[1, ]), colnames(m))
    })
    values <- do.call(c, unname(parts))
    if(is.null(names(values)) || anyDuplicated(names(values))) {
      stop("Missing or overlapping cell names in normalized layers.")
    }
    if(!all(wanted.cells %in% names(values))) {
      stop("Selected cells are missing from normalized RNA layers.")
    }
    values <- values[wanted.cells]
    if(any(!is.finite(values)) || any(values < 0)) {
      stop("Expected finite, nonnegative LogNormalize expression values.")
    }
    unname(values)
  }

  tissue.data <- list()
  inventories <- list()
  for(tissue in tissues) {
    message("Reading: ", rds.files[[tissue]])
    x <- readRDS(rds.files[[tissue]])
    if(!inherits(x, "Seurat")) stop("Expected Seurat object: ", tissue)
    if(!assay %in% names(x@assays)) stop("Missing RNA assay: ", tissue)
    if(!reduction %in% names(x@reductions)) stop("Missing saved UMAP: ", tissue)
    md <- x@meta.data
    if(!age.column %in% names(md)) stop("Missing age metadata: ", tissue)
    cells <- colnames(x)
    if(anyDuplicated(cells) || !all(cells %in% rownames(md))) {
      stop("Invalid cell names or metadata coverage: ", tissue)
    }
    md <- md[cells, , drop = FALSE]
    ages <- as.character(md[[age.column]])
    # Exact metadata labels match the supplied analysis scripts.
    if(anyNA(ages) || any(!ages %in% c(included.ages, "18m", "21m"))) {
      stop("Unexpected age labels in ", tissue, "; expected 1m/3m/18m/21m/24m/30m.")
    }
    absent <- setdiff(included.ages, ages)
    if(length(absent)) stop("Missing requested ages in ", tissue, ": ", paste(absent, collapse = ", "))
    inventories[[tissue]] <- data.frame(
      tissue = tissue, age = names(table(ages)), n_cells = as.integer(table(ages)),
      included = names(table(ages)) %in% included.ages
    )
    keep <- ages %in% included.ages
    selected.cells <- cells[keep]
    embedding <- SeuratObject::Embeddings(x[[reduction]])
    if(ncol(embedding) < 2L || anyDuplicated(rownames(embedding)) ||
       !all(selected.cells %in% rownames(embedding))) stop("Incomplete UMAP: ", tissue)
    embedding <- embedding[selected.cells, 1:2, drop = FALSE]
    if(any(!is.finite(embedding))) stop("Nonfinite UMAP coordinates: ", tissue)
    tissue.data[[tissue]] <- data.frame(
      cell = selected.cells, tissue = tissue, age = ages[keep],
      mouse.id = if("mouse.id" %in% names(md)) as.character(md$mouse.id[keep]) else NA_character_,
      UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2],
      expression = read.expression(x, selected.cells), row.names = NULL
    )
    rm(x)
    invisible(gc())
  }
  # Groups overlap (3m also belongs to Young); pool tissue tables, not plot groups.
  pooled <- unlist(lapply(tissue.data, function(d) d$expression), use.names = FALSE)
  positive <- pooled[pooled > 0]
  actual.max <- max(pooled)
  q.limit <- if(length(positive)) {
    unname(stats::quantile(positive, probs = robust.quantile, type = 7))
  } else 0
  caps <- c(actual_max = actual.max, positive_quantile = q.limit)
  # Nondegenerate display range for the all-zero edge case; actual limits are audited.
  display.caps <- ifelse(caps > 0, caps, 1)
  if(actual.max == 0) warning("Mff is zero in all selected cells; all dots will be grey.")
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, inventories), file.path(output.dir, "input_age_inventory.csv"),
            row.names = FALSE)
  write.csv(data.frame(
    mode = names(caps), reference_upper = unname(caps),
    display_upper = unname(display.caps),
    quantile_probability = c(NA_real_, robust.quantile),
    scope = "Both tissues; unique cells aged 1m,3m,24m,30m",
    n_reference_cells = length(pooled), n_positive_cells = length(positive)
  ), file.path(output.dir, "color_scale_limits.csv"), row.names = FALSE)

  save.plot <- function(p, stem, paired = FALSE) {
    w <- width.in * if(paired) 2 else 1
    ggsave(paste0(stem, ".png"), plot = p, width = w, height = height.in,
           units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
    if(save.pdf) {
      ggsave(paste0(stem, ".pdf"), plot = p, width = w, height = height.in,
             units = "in", bg = "white", limitsize = FALSE)
    }
  }
  qc <- list()
  for(tissue in tissues) {
    d <- tissue.data[[tissue]]
    # Same axes and same saved embedding for every group and color mode in this tissue.
    xlim <- range(d$UMAP_1)
    ylim <- range(d$UMAP_2)
    if(diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)
    if(diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
    write.csv(d, file.path(output.dir, paste0(tissue, ".Mff_UMAP.plot_data.csv")),
              row.names = FALSE)
    for(mode in names(caps)) {
      cap <- display.caps[[mode]]
      folder <- file.path(output.dir, tissue, mode)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      plots <- list()
      for(group in names(group.ages)) {
        z <- d[d$age %in% group.ages[[group]], , drop = FALSE]
        # All selected cells remain; draw high-expression dots last for visibility.
        z <- z[order(z$expression, z$cell), , drop = FALSE]
        mode.label <- if(mode == "actual_max") "Shared actual maximum" else {
          paste0("Shared positive-expression Q", format(100 * robust.quantile, trim = TRUE))
        }
        p <- ggplot(z, aes(UMAP_1, UMAP_2, color = expression)) +
          geom_point(size = point.size, alpha = 1, shape = 16, stroke = 0) +
          scale_color_gradientn(
            colours = expression.colors, limits = c(0, cap), oob = scales::squish,
            breaks = c(0, cap / 2, cap),
            labels = function(v) format(signif(v, 3), trim = TRUE),
            name = paste0(gene, "\nlog1p-normalized"),
            guide = guide_colorbar(barheight = grid::unit(4, "cm"))
          ) +
          coord_fixed(xlim = xlim, ylim = ylim, expand = TRUE) +
          labs(
            title = paste(gsub("_", " ", tissue), group.labels[[group]], sep = " | "),
            subtitle = paste0(mode.label, " = ", signif(caps[[mode]], 4)),
            caption = paste0(
              "n = ", nrow(z), " cells; Mff > 0: ", sum(z$expression > 0),
              if(mode != "actual_max") "\nValues above the color limit use the maximum color." else ""
            ),
            x = "UMAP 1", y = "UMAP 2"
          ) +
          theme_classic(base_size = font.base) +
          theme(
            plot.title = element_text(size = font.title, face = "bold"),
            axis.title = element_text(size = font.axis),
            axis.text = element_text(size = font.axis),
            legend.title = element_text(size = font.legend),
            legend.text = element_text(size = font.legend),
            plot.margin = margin(12, 12, 12, 12)
          )
        plots[[group]] <- p
        save.plot(p, file.path(folder, paste0(tissue, ".Mff.", group)))
        qc[[length(qc) + 1L]] <- data.frame(
          tissue = tissue, group = group, mode = mode, n_cells = nrow(z),
          n_positive = sum(z$expression > 0), panel_max = max(z$expression),
          shared_color_upper = caps[[mode]],
          n_above_color_upper = sum(z$expression > caps[[mode]])
        )
      }
      for(pair in list(c("Young_1m_3m", "Old_24m_30m"), c("Age_3m", "Age_30m"))) {
        combined <- patchwork::wrap_plots(plots[pair], ncol = 2, guides = "collect")
        save.plot(combined, file.path(folder, paste0(tissue, ".Mff.", paste(pair, collapse = "_vs_"))),
                  paired = TRUE)
      }
    }
  }
  write.csv(do.call(rbind, qc), file.path(output.dir, "panel_summary.csv"), row.names = FALSE)
  writeLines(c(
    paste0("Generated: ", Sys.time()),
    paste0("Input: ", unname(rds.files)),
    "Expression: RNA data layer, LogNormalize(scale.factor=10000), natural-log log1p units.",
    "No normalization, scaling, PCA, clustering or UMAP is rerun; no input RDS is overwritten.",
    "Saved tissue UMAPs were fitted by 3.Umap.tissue.R using all ages present at that time.",
    "18m/21m cells are not plotted, but their contribution to existing UMAP coordinates remains.",
    "Heart and Limb Muscle have independent UMAP embeddings; positions across tissues are not aligned.",
    "Young = 1m + 3m; Old = 24m + 30m. Additional views: 3m and 30m separately.",
    "Only each group's own cells are plotted; zero-expression cells remain grey.",
    "Within each tissue, all views share axes; both tissues and all groups share each color limit.",
    "Actual maximum uses all unique selected cells across both tissues.",
    paste0("Robust limit = quantile of positive Mff expression, p=", robust.quantile, ", R type=7."),
    "This quantile is a visualization choice, not a universal standard or biological cutoff.",
    "Robust clipping affects colors only; plot_data CSV retains original expression.",
    "All-zero fallback: reference upper=0, display upper=1; all dots remain grey.",
    "16 individual PNGs + 8 paired PNGs; optional matching PDFs.",
    "Figures describe individual cells; no mouse balancing or statistical testing is performed.",
    "https://satijalab.org/seurat/reference/featureplot",
    "", capture.output(sessionInfo())
  ), file.path(output.dir, "analysis_settings.txt"))
  message("Saved 16 individual and 8 paired PNGs to: ", output.dir)
  invisible(list(output.dir = output.dir, color.limits = caps))
}

run.mff.umap()
