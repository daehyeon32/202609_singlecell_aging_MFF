# Seurat FeaturePlot cell-type crops of saved tissue UMAPs; no embedding refit.
# Run the WHOLE file with source(). Requires outputs of 3.Umap.tissue.R.
# Existing normalization and UMAP coordinates are reused, not recomputed.
# References:
# https://satijalab.org/seurat/reference/featureplot
# https://satijalab.github.io/seurat-object/reference/Layers.html

run.mff.celltype.featureplot <- function() {
  # --------------------------- SETTINGS ---------------------------
  input.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split"
  output.dir <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/11.Umap_Mff_expression/FeaturePlot_celltype_crops.132430"
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  tissues <- c("Heart", "Limb_Muscle")
  gene <- "Mff"
  assay <- "RNA"
  reduction <- "umap"
  age.column <- "age"
  celltype.column <- "cell_ontology_class"
  crop.padding <- 0.08  # 8% of the square crop span on EACH side
  min.crop.fraction <- 0.02  # Minimum span relative to tissue UMAP for tiny groups
  save.overview <- TRUE  # Additional 2 x 2 view: 3m, 30m, Young, Old
  title.wrap.width <- 40
  included.ages <- c("1m", "3m", "24m", "30m")
  group.ages <- list(
    Age_3m = "3m", Age_30m = "30m",
    Young_1m_3m = c("1m", "3m"), Old_24m_30m = c("24m", "30m")
  )
  group.labels <- c(
    Young_1m_3m = "Young (1 + 3 months)", Old_24m_30m = "Old (24 + 30 months)",
    Age_3m = "3 months", Age_30m = "30 months"
  )
  # Both limits are computed ONCE over unique selected cells from BOTH tissues.
  # Q95 is a visualization choice, not a biological threshold or universal standard.
  robust.quantile <- 0.95
  point.size <- 0.50  # FeaturePlot pt.size, matching script 7
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

  required <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "scales")
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
  if(length(crop.padding) != 1L || !is.finite(crop.padding) || crop.padding < 0) {
    stop("crop.padding must be a finite nonnegative number.")
  }
  if(length(min.crop.fraction) != 1L || !is.finite(min.crop.fraction) ||
     min.crop.fraction <= 0) stop("min.crop.fraction must be positive.")
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
    if(!all(c(age.column, celltype.column) %in% names(md))) {
      stop("Missing age or cell-type metadata: ", tissue)
    }
    cells <- colnames(x)
    if(anyDuplicated(cells) || !all(cells %in% rownames(md))) {
      stop("Invalid cell names or metadata coverage: ", tissue)
    }
    md <- md[cells, , drop = FALSE]
    ages <- as.character(md[[age.column]])
    celltypes <- as.character(md[[celltype.column]])
    celltypes[is.na(celltypes) | trimws(celltypes) == ""] <- "Unannotated"
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
      cell = selected.cells, tissue = tissue, age = ages[keep], celltype = celltypes[keep],
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

  # Shared square coordinate window; panel aspect is flexible as in script 7.
  # Bounds pool ALL four selected ages within a cell type, BEFORE group selection.
  crop.bounds <- function(d, tissue.span) {
    xr <- range(d$UMAP_1)
    yr <- range(d$UMAP_2)
    span <- max(diff(xr), diff(yr), tissue.span * min.crop.fraction, 1e-6)
    half <- span * (0.5 + crop.padding)
    c(xmin = mean(xr) - half, xmax = mean(xr) + half,
      ymin = mean(yr) - half, ymax = mean(yr) + half)
  }
  save.plot <- function(p, stem, overview = FALSE) {
    multiplier <- if(overview) 2 else 1
    ggsave(paste0(stem, ".png"), plot = p,
           width = width.in * multiplier, height = height.in * multiplier,
           units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
    if(save.pdf) {
      ggsave(paste0(stem, ".pdf"), plot = p,
             width = width.in * multiplier, height = height.in * multiplier,
             units = "in", bg = "white", limitsize = FALSE)
    }
  }

  qc <- list()
  crops <- list()
  for(tissue in tissues) {
    d <- tissue.data[[tissue]]
    # Reload one tissue at a time after computing the shared color limits.
    object <- readRDS(rds.files[[tissue]])
    SeuratObject::DefaultAssay(object) <- assay
    object <- subset(object, cells = d$cell)
    data.layers <- SeuratObject::Layers(object[[assay]], search = "^data($|\\.)")
    if(length(data.layers) != 1L || data.layers != "data") {
      object <- SeuratObject::JoinLayers(
        object, assay = assay, layers = "^data($|\\.)", new = "data"
      )
    }
    # Confirm that joining layers preserves the expression used to set the limits.
    if(!isTRUE(all.equal(read.expression(object, d$cell), d$expression,
                        tolerance = 1e-10, check.attributes = FALSE))) {
      stop("Expression changed between scale calculation and plotting: ", tissue)
    }
    saved.xy <- SeuratObject::Embeddings(object[[reduction]])[d$cell, 1:2, drop = FALSE]
    if(!isTRUE(all.equal(unname(saved.xy), unname(as.matrix(d[, c("UMAP_1", "UMAP_2")])),
                        check.attributes = FALSE))) {
      stop("UMAP coordinates changed between reads: ", tissue)
    }
    tissue.span <- max(diff(range(d$UMAP_1)), diff(range(d$UMAP_2)), 1e-6)
    celltypes <- sort(unique(d$celltype))
    # The numeric prefix prevents collisions after replacing special characters.
    file.map <- data.frame(
      celltype = celltypes,
      file_stem = sprintf("%03d.%s", seq_along(celltypes),
                          substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 80))
    )
    tissue.dir <- file.path(output.dir, tissue)
    dir.create(tissue.dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(file.map, file.path(tissue.dir, "celltype_file_map.csv"), row.names = FALSE)
    # Raw embedding coordinates and expression: no translation, rotation, or rescaling.
    write.csv(d, file.path(tissue.dir, "Mff_UMAP.plot_data.csv"), row.names = FALSE)

    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      dc <- d[d$celltype == ct, , drop = FALSE]
      bounds <- crop.bounds(dc, tissue.span)
      if(any(dc$UMAP_1 < bounds["xmin"] | dc$UMAP_1 > bounds["xmax"] |
             dc$UMAP_2 < bounds["ymin"] | dc$UMAP_2 > bounds["ymax"])) {
        stop("Crop unexpectedly excludes cells: ", tissue, " / ", ct)
      }
      crops[[length(crops) + 1L]] <- data.frame(
        tissue = tissue, celltype = ct, n_cells = nrow(dc),
        xmin = unname(bounds["xmin"]), xmax = unname(bounds["xmax"]),
        ymin = unname(bounds["ymin"]), ymax = unname(bounds["ymax"]),
        padding_fraction_per_side = crop.padding
      )
      for(mode in names(caps)) {
        cap <- display.caps[[mode]]
        folder <- file.path(tissue.dir, file.map$file_stem[[i]], mode)
        dir.create(folder, recursive = TRUE, showWarnings = FALSE)
        plots <- list()
        for(group in names(group.ages)) {
          z <- dc[dc$age %in% group.ages[[group]], , drop = FALSE]
          z <- z[order(z$expression, z$cell), , drop = FALSE]
          mode.label <- if(mode == "actual_max") "Shared actual maximum" else {
            paste0("Shared positive-expression Q", format(100 * robust.quantile, trim = TRUE))
          }
          title <- paste(
            paste(strwrap(paste(gsub("_", " ", tissue), ct, sep = " | "),
                          width = title.wrap.width), collapse = "\n"),
            group.labels[[group]], sep = "\n"
          )
          if(nrow(z)) {
            p <- Seurat::FeaturePlot(
              object = object, features = gene, cells = z$cell,
              reduction = reduction, slot = "data", order = TRUE,
              min.cutoff = 0, max.cutoff = cap, keep.scale = NULL,
              cols = c(expression.colors[1], tail(expression.colors, 1)),
              pt.size = point.size, label = FALSE, raster = FALSE,
              coord.fixed = FALSE, combine = FALSE
            )[[1]]
          } else {
            # FeaturePlot cannot plot an empty cell selection. Copy its theme
            # from this cell type, but use a genuinely empty data/layer panel.
            template <- Seurat::FeaturePlot(
              object = object, features = gene, cells = dc$cell,
              reduction = reduction, slot = "data", order = TRUE,
              min.cutoff = 0, max.cutoff = cap, keep.scale = NULL,
              cols = c(expression.colors[1], tail(expression.colors, 1)),
              pt.size = point.size, label = FALSE, raster = FALSE,
              coord.fixed = FALSE, combine = FALSE
            )[[1]]
            p <- ggplot(z, aes(UMAP_1, UMAP_2, color = expression)) +
              geom_blank() + template$theme
          }
          # Explicit continuous scale keeps both color modes comparable across groups.
          # Replacing FeaturePlot's color scale may print an informational message.
          p <- p +
            scale_color_gradientn(
              colours = expression.colors, limits = c(0, cap), oob = scales::squish,
              breaks = c(0, cap / 2, cap),
              labels = function(v) format(signif(v, 3), trim = TRUE),
              name = paste0(gene, "\nlog1p-normalized"),
              guide = guide_colorbar(barheight = grid::unit(4, "cm"))
            ) +
            coord_cartesian(
              xlim = unname(bounds[c("xmin", "xmax")]),
              ylim = unname(bounds[c("ymin", "ymax")]), expand = FALSE
            ) +
            labs(
              title = title,
              subtitle = paste0(mode.label, " = ", signif(caps[[mode]], 4)),
              caption = paste0(
                "n = ", nrow(z), " cells; Mff > 0: ", sum(z$expression > 0),
                if(mode != "actual_max") "\nValues above the color limit use the maximum color." else ""
              ),
              x = "umap_1", y = "umap_2"
            ) +
            theme(
              text = element_text(size = font.base),
              plot.title = element_text(size = font.title, face = "bold", hjust = 0.5),
              axis.title = element_text(size = font.axis),
              axis.text = element_text(size = font.axis),
              legend.title = element_text(size = font.legend),
              legend.text = element_text(size = font.legend),
              plot.margin = margin(12, 12, 12, 12)
            )
          if(!nrow(z)) {
            # An absent cell type/group is not zero expression: keep an empty labeled panel.
            p <- p + annotate(
              "text", x = mean(bounds[c("xmin", "xmax")]),
              y = mean(bounds[c("ymin", "ymax")]),
              label = "No cells in this group", color = "#666666", size = 5
            )
          }
          plots[[group]] <- p
          save.plot(p, file.path(folder, paste0("Mff.", group)))
          qc[[length(qc) + 1L]] <- data.frame(
            tissue = tissue, celltype = ct, group = group, mode = mode,
            ages = paste(group.ages[[group]], collapse = ";"),
            n_cells = nrow(z), n_positive = sum(z$expression > 0),
            panel_max = if(nrow(z)) max(z$expression) else NA_real_,
            shared_color_upper = caps[[mode]],
            n_above_color_upper = sum(z$expression > caps[[mode]]),
            status = if(nrow(z)) "plotted" else "no_cells"
          )
        }
        if(save.overview) {
          # 3m | 30m on the top row, Young | Old on the bottom row.
          p <- patchwork::wrap_plots(plots[names(group.ages)], ncol = 2, guides = "collect")
          save.plot(p, file.path(folder, "Mff.four_groups"), overview = TRUE)
        }
      }
      message("Saved: ", tissue, " / ", ct, " [", i, "/", length(celltypes), "]")
    }
    rm(object)
    invisible(gc())
  }
  crop.table <- do.call(rbind, crops)
  write.csv(crop.table, file.path(output.dir, "crop_bounds.csv"), row.names = FALSE)
  write.csv(do.call(rbind, qc), file.path(output.dir, "panel_summary.csv"), row.names = FALSE)
  writeLines(c(
    paste0("Generated: ", Sys.time()), paste0("Input: ", unname(rds.files)),
    paste0("Cell-type metadata: ", celltype.column),
    "Uses the saved tissue-specific UMAP embedding and RNA normalized data.",
    "No NormalizeData, ScaleData, PCA, RunUMAP, clustering, or new RDS output.",
    "Only the requested cell type and group are drawn; other cell types are omitted.",
    "3m; 30m; Young=1m+3m; Old=24m+30m. 18m/21m cells are not plotted.",
    "The saved UMAP may have been fitted with 18m/21m; its coordinates remain unchanged.",
    "Crop bounds pool the cell type's unique cells aged 1m/3m/24m/30m.",
    "Each cell type shares exactly one square crop across its four groups and both color modes.",
    paste0("Additional padding on EACH side = ", crop.padding, " times the unpadded square span."),
    paste0("Minimum unpadded span = ", min.crop.fraction, " times the selected tissue UMAP span."),
    "Original numeric coordinates and orientation are preserved; no embedding is refitted.",
    "FeaturePlot theme and flexible panel aspect match script 7; no 1:1 display aspect lock.",
    "Joined RNA data layers in memory; expression and embedding checked before plotting.",
    "Different cell types have different zoom levels; compare positions only within a cell type/tissue.",
    "All islands and outliers remain inside the crop; no quantile cropping or cell removal.",
    "An elongated or disconnected cell type can retain internal blank space.",
    "A group with no cells gets an empty panel labeled No cells in this group, not zero expression.",
    "Unannotated includes missing/blank cell-type labels; no selected cells are discarded.",
    "Color limits pool unique selected cells across BOTH tissues and ALL cell types, as in script 5.",
    "Actual maximum and positive-expression quantile are shared across all figures.",
    paste0("Positive-expression quantile probability = ", robust.quantile, "; R type=7."),
    "Quantile cap is a display choice; it clips colors only. Original values remain in CSV.",
    "Zero expression remains grey. All-zero data use display upper=1, reference upper=0.",
    "Within each color mode: four individual PNGs per cell type.",
    paste0("Additional 2x2 overview per cell type/color mode: ", save.overview),
    paste0("Point size: ", point.size, "; individual canvas: ", width.in, " x ", height.in,
           " inches; dpi=", dpi, "; PDF=", save.pdf),
    "No statistical tests or mouse balancing are performed.",
    "", capture.output(sessionInfo())
  ), file.path(output.dir, "analysis_settings.txt"))
  message("Completed ", nrow(crop.table), " tissue-celltype crops. Output: ", output.dir)
  invisible(list(output.dir = output.dir, crop.bounds = crop.table, color.limits = caps))
}

run.mff.celltype.featureplot()
