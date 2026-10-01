
local({
  # 1. 사용자 설정: 이 부분에서 경로와 그림 크기를 조절합니다.
  gene <- "Mff"
  rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression"
  tissues <- c("Heart", "Limb_Muscle")
  tissue.column <- "tissue_free_annotation"
  include.heart.and.aorta <- TRUE  # 첨부한 최신 Young/Old 코드와 동일한 조직 선택
  age.order <- c(1, 3, 18, 21, 24, 30)
  analysis.ages <- list(full_6ages = c(1, 3, 18, 21, 24, 30), sensitivity_4ages = c(3, 18, 21, 24))
  analysis.labels <- c(full_6ages = "All ages: 1, 3, 18, 21, 24, 30 months", sensitivity_4ages = "Sensitivity: 3, 18, 21, 24 months")
  line.colors <- c(full_6ages = "#3C5488", sensitivity_4ages = "#C85946")
  age.colors <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"), as.character(age.order))
  min.cells.per.mouse.celltype <- 1L
  fdr.threshold <- 0.05
  show.age.means <- FALSE  # TRUE: 마우스별 log2(CPM+1)의 연령별 평균을 마름모로 추가
  save.pdf <- TRUE
  save.raw.pseudobulk <- TRUE
  plot.width <- 1200
  plot.height <- 900
  plot.y.min <- 4  # Y축 표시 최솟값; 회귀 계산에는 이보다 작은 값도 포함
  plot.res <- 160
  plot.base.size <- 13
  plot.title.size <- 18
  plot.axis.title.size <- 16
  plot.axis.text.size <- 30
  plot.annotation.size <- 8
  plot.title.wrap.width <- 45
  point.size <- 3
  point.jitter.months <- 0.12  # 화면에서만 좌우로 이동; 회귀는 원래 개월 수 사용

  required <- c("Seurat", "SeuratObject", "edgeR", "statmod", "Matrix", "dplyr", "ggplot2", "patchwork", "Cairo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing) > 0){stop("Missing R packages: ", paste(missing, collapse = ", "))}
  if(packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0"){stop("Seurat and SeuratObject >= 5.0.0 are required.")}
  if(packageVersion("edgeR") < "4.0.0"){stop("edgeR >= 4.0.0 is required for legacy=FALSE; no fallback to an older test.")}
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(!file.exists(rds.file)){stop("Missing input RDS: ", rds.file)}
  if(length(min.cells.per.mouse.celltype) != 1L || !is.finite(min.cells.per.mouse.celltype) || min.cells.per.mouse.celltype < 1 || min.cells.per.mouse.celltype != floor(min.cells.per.mouse.celltype)){stop("min.cells.per.mouse.celltype must be an integer >= 1.")}
  if(!all(vapply(analysis.ages, function(a) all(a %in% age.order) && !anyDuplicated(a) && length(a) >= 2, logical(1)))){stop("Invalid analysis.ages.")}
  if(!all(analysis.ages$sensitivity_4ages %in% analysis.ages$full_6ages)){stop("Sensitivity ages must be a subset of the full analysis ages.")}

  # 2. 공통 함수
  canonical.tissue <- function(z) {gsub("[[:space:]-]+", "_", tolower(trimws(as.character(z))))}
  join.values <- function(z) {paste(sort(unique(as.character(z[!is.na(z) & trimws(as.character(z)) != ""]))), collapse = ";")}
  gene.value <- function(x, i) {if(is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x[if(length(x) == 1L) 1L else i])}
  log2.cpm.plus1 <- function(log.cpm) {(pmax(log.cpm, 0) + log1p(exp(-abs(log.cpm)))) / log(2)}
  capture.call <- function(expr) {
    notes <- character()
    out <- withCallingHandlers(tryCatch(list(value = force(expr), error = NULL), error = function(e) list(value = NULL, error = conditionMessage(e))), warning = function(w) {notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")})
    out$warnings <- join.values(notes)
    out
  }
  save.plot <- function(p, stem, width = plot.width, height = plot.height) {
    Cairo::CairoPNG(filename = paste0(stem, ".png"), width = width, height = height, res = plot.res, bg = "white")
    tryCatch(print(p), finally = grDevices::dev.off())
    if(save.pdf) {
      Cairo::CairoPDF(file = paste0(stem, ".pdf"), width = width / plot.res, height = height / plot.res, bg = "white")
      tryCatch(print(p), finally = grDevices::dev.off())
    }
  }
  join.count.layer <- function(object) {
    layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts($|\\.)")
    if(length(layers) == 0L){stop("RNA raw-count layer is missing.")}
    cells <- unlist(lapply(layers, function(z) colnames(SeuratObject::LayerData(object[["RNA"]], layer = z))), use.names = FALSE)
    if(anyDuplicated(cells) || !setequal(cells, colnames(object))){stop("Duplicated or missing cells across RNA count layers.")}
    if(length(layers) != 1L || layers != "counts") {
      if(!inherits(object[["RNA"]], "Assay5")){stop("Unexpected split RNA counts in a non-Assay5 object.")}
      object <- SeuratObject::JoinLayers(object, assay = "RNA", layers = "^counts($|\\.)", new = "counts")
    }
    object
  }
  empty.curve <- function() {data.frame(tissue = character(), analysis = character(), celltype = character(), age_months = numeric(), predicted_CPM = numeric(), predicted_log2_CPM_plus1 = numeric(), stringsAsFactors = FALSE)}

  # 각 celltype의 모든 유전자로 정규화하고 QL 변동성을 추정합니다.
  # Mff 한 유전자만 뽑아서 모형을 적합하면 안 됩니다.
  fit.one <- function(pb.counts, pb.meta, ct, tissue.name, mode) {
    d <- pb.meta[pb.meta$celltype == ct & pb.meta$age_months %in% analysis.ages[[mode]], , drop = FALSE]
    d <- d[d$eligible, , drop = FALSE]
    d <- as.data.frame(d[order(d$age_months, d$mouse.id), , drop = FALSE])
    if(anyDuplicated(d$mouse.id)){stop("Repeated mouse in model: ", tissue.name, " / ", ct)}
    d$analysis <- rep(mode, nrow(d))
    for(nm in c("norm_factor", "effective_library_size", "normalized_CPM", "log2_CPM_plus1", "predicted_CPM", "predicted_log2_CPM_plus1")){d[[nm]] <- rep(NA_real_, nrow(d))}
    d$normalization_method <- rep("not_normalized", nrow(d))
    ages <- sort(unique(d$age_months))
    r <- data.frame(tissue = tissue.name, analysis = mode, celltype = ct, gene = gene, requested_ages = paste(analysis.ages[[mode]], collapse = ";"), observed_ages = paste(ages, collapse = ";"), n_mice = nrow(d), n_observed_ages = length(ages), n_mice_Mff_nonzero = sum(d$gene_raw_count > 0), n_cells = sum(d$n_cells), n_genes_nonzero = NA_integer_, residual_df_nominal = NA_real_, residual_df_adjusted_Mff = NA_real_, ql_dispersion_Mff = NA_real_, ql_prior_df_Mff = NA_real_, F_statistic = NA_real_, df_numerator = NA_real_, df_denominator = NA_real_, log2FC_per_month = NA_real_, log2FC_per_6months = NA_real_, fold_change_per_6months = NA_real_, p_value = NA_real_, p_adj_BH = NA_real_, direction = NA_character_, status = "not_tested", reason = "", normalization_warning = "", model_warning = "", scope_note = if(length(ages) == 2L) "Only two observed ages: the slope represents this age contrast." else "", stringsAsFactors = FALSE)
    curve <- empty.curve()
    finish <- function(status, reason) {r$status <- status; r$reason <- reason; list(stats = r, mice = d, curve = curve)}
    if(nrow(d) == 0L){return(finish("no_eligible_mice", "No eligible mouse in the selected age range."))}

    mat <- as.matrix(pb.counts[, d$pb_id, drop = FALSE])
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    r$n_genes_nonzero <- nrow(mat)
    normalized <- capture.call({
      y <- edgeR::DGEList(counts = mat, lib.size = d$library_size)
      if(ncol(y) >= 2L){y <- edgeR::calcNormFactors(y, method = "TMM")}
      effective <- y$samples$lib.size * y$samples$norm.factors
      if(any(!is.finite(effective)) || any(effective <= 0)){stop("Nonpositive/nonfinite effective library size.")}
      y
    })
    r$normalization_warning <- normalized$warnings
    if(!is.null(normalized$error)){return(finish("normalization_failed", normalized$error))}
    y <- normalized$value
    d$norm_factor <- y$samples$norm.factors
    d$effective_library_size <- y$samples$lib.size * d$norm_factor
    d$normalized_CPM <- d$gene_raw_count / d$effective_library_size * 1e6
    d$log2_CPM_plus1 <- log2(d$normalized_CPM + 1)
    d$normalization_method <- if(nrow(d) >= 2L) "TMM_refitted_within_analysis" else "CPM_only_single_mouse"
    if(!gene %in% rownames(y)){return(finish("Mff_all_zero", "Mff has zero counts in every included mouse; no slope test."))}
    if(length(ages) < 2L){return(finish("only_one_observed_age", "At least two distinct observed ages are needed."))}
    d$age_per6 <- (d$age_months - 3) / 6
    design <- stats::model.matrix(~ age_per6, data = d)
    rownames(design) <- d$pb_id
    rank <- qr(design)$rank
    r$residual_df_nominal <- nrow(d) - rank
    if(rank < ncol(design)){return(finish("rank_deficient_design", "Age slope is not identifiable."))}
    if(r$residual_df_nominal <= 0){return(finish("no_residual_df", "The numeric-age model needs more mice than its two coefficients."))}
    if(!identical(rownames(design), colnames(y))){stop("Design/count sample order mismatch.")}

    modeled <- capture.call({
      fit <- edgeR::glmQLFit(y, design = design, robust = TRUE, legacy = FALSE, prior.count = 0)
      gi <- match(gene, rownames(fit$coefficients))
      adjusted.df <- gene.value(fit$df.residual.adj, gi)
      if(!is.finite(adjusted.df) || adjusted.df <= 0){stop("Mff has no usable adjusted residual degrees of freedom.")}
      if(!is.null(fit$failed) && isTRUE(as.logical(gene.value(fit$failed, gi)))){stop("Mff model failed to converge.")}
      if(any(!is.finite(fit$coefficients[gi, ]))){stop("Nonfinite Mff coefficients.")}
      qlf <- edgeR::glmQLFTest(fit, coef = "age_per6")
      pv <- qlf$table[gene, "PValue"]
      if(length(pv) != 1L || !is.finite(pv) || pv < 0 || pv > 1){stop("Invalid QL trend p-value.")}
      beta <- as.numeric(fit$coefficients[gi, ])
      observed.log.cpm <- as.numeric(design %*% beta) + log(1e6)
      expected <- as.numeric(fit$fitted.values[gi, ]) / d$effective_library_size * 1e6
      if(!isTRUE(all.equal(exp(observed.log.cpm), expected, tolerance = 1e-6, check.attributes = FALSE))){stop("Prediction/offset consistency check failed.")}
      grid <- sort(unique(c(seq(min(ages), max(ages), length.out = 200L), ages)))
      new.design <- cbind(1, (grid - 3) / 6)
      predicted.log.cpm <- as.numeric(new.design %*% beta) + log(1e6)
      if(any(!is.finite(exp(predicted.log.cpm)))){stop("Predicted CPM overflows; inspect the sparse target counts.")}
      curve <- data.frame(tissue = tissue.name, analysis = mode, celltype = ct, age_months = grid, predicted_CPM = exp(predicted.log.cpm), predicted_log2_CPM_plus1 = log2.cpm.plus1(predicted.log.cpm), stringsAsFactors = FALSE)
      list(fit = fit, qlf = qlf, gi = gi, beta = beta, curve = curve, observed_log_cpm = observed.log.cpm)
    })
    r$model_warning <- modeled$warnings
    if(!is.null(modeled$error)){return(finish("model_or_test_failed", modeled$error))}
    z <- modeled$value
    r$residual_df_adjusted_Mff <- gene.value(z$fit$df.residual.adj, z$gi)
    r$ql_dispersion_Mff <- gene.value(z$fit$s2.post, z$gi)
    r$ql_prior_df_Mff <- gene.value(z$fit$df.prior, z$gi)
    r$F_statistic <- z$qlf$table[gene, "F"]
    r$df_numerator <- gene.value(z$qlf$df.test, z$gi)
    r$df_denominator <- gene.value(z$qlf$df.total, z$gi)
    r$log2FC_per_6months <- z$beta[[2]] / log(2)
    r$log2FC_per_month <- r$log2FC_per_6months / 6
    r$fold_change_per_6months <- exp(z$beta[[2]])
    r$p_value <- z$qlf$table[gene, "PValue"]
    r$direction <- if(r$log2FC_per_6months > 0) "increasing" else if(r$log2FC_per_6months < 0) "decreasing" else "flat"
    r$status <- "tested"
    d$predicted_CPM <- exp(z$observed_log_cpm)
    d$predicted_log2_CPM_plus1 <- log2.cpm.plus1(z$observed_log_cpm)
    list(stats = r, mice = d, curve = z$curve)
  }

  make.plot <- function(result, ct, tissue.name, mode, ymax) {
    d <- result$mice[result$mice$celltype == ct & is.finite(result$mice$log2_CPM_plus1), , drop = FALSE]
    curve <- result$curve[result$curve$celltype == ct, , drop = FALSE]
    r <- result$stats[result$stats$celltype == ct, , drop = FALSE]
    label <- if(r$status == "tested") paste0("p (trend) ", if(r$p_value < 1e-4) "< 0.0001" else paste0("= ", format.pval(r$p_value, digits = 3)), "\nlog2FC / 6 months = ", sprintf("%+.3f", r$log2FC_per_6months)) else paste0("p (trend) = NA\n", gsub("_", " ", r$status))
    title <- paste(strwrap(paste0(tissue.name, " | ", ct), width = plot.title.wrap.width), collapse = "\n")
    p <- ggplot(d, aes(x = age_months, y = log2_CPM_plus1))
    if(nrow(curve) > 0L){p <- p + geom_line(data = curve, aes(x = age_months, y = predicted_log2_CPM_plus1), inherit.aes = FALSE, colour = line.colors[[mode]], linewidth = 1.1)}
    p <- p + geom_point(aes(colour = factor(age_months, levels = age.order)), position = position_jitter(width = point.jitter.months, height = 0, seed = 1234), size = point.size, alpha = 0.85)
    if(show.age.means && nrow(d) > 0L) {
      means <- d %>% group_by(age_months) %>% summarise(value = mean(log2_CPM_plus1), .groups = "drop")
      p <- p + geom_point(data = means, aes(x = age_months, y = value), inherit.aes = FALSE, shape = 23, fill = "white", colour = "black", size = 3)
    }
    p + scale_colour_manual(values = age.colors, drop = FALSE) + scale_x_continuous(breaks = analysis.ages[[mode]], labels = analysis.ages[[mode]]) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + coord_cartesian(xlim = range(age.order), ylim = c(plot.y.min, ymax)) + annotate("text", x = Inf, y = Inf, label = label, hjust = 1.03, vjust = 1.25, size = plot.annotation.size) + labs(x = "Age (months)", y = paste0(gene, " log2(TMM-normalized\nCPM + 1)"), title = title, subtitle = NULL) + theme_bw(base_size = plot.base.size) + theme(legend.position = "none", plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(), axis.title = element_text(size = plot.axis.title.size), axis.text = element_text(size = plot.axis.text.size), panel.grid.minor = element_blank())
  }

  # 3. 통합 RDS를 한 번 읽고 조직별 raw pseudobulk를 만듭니다.
  mmf.rds <- readRDS(rds.file)
  if(!inherits(mmf.rds, "Seurat") || !"RNA" %in% names(mmf.rds@assays)){stop("Expected a Seurat object with RNA assay.")}
  needed <- c(tissue.column, "age", "mouse.id", "cell_ontology_class")
  if(!all(needed %in% colnames(mmf.rds@meta.data))){stop("Missing metadata: ", paste(setdiff(needed, colnames(mmf.rds@meta.data)), collapse = ", "))}
  tissue.keys <- canonical.tissue(mmf.rds[[tissue.column, drop = TRUE]])
  if(include.heart.and.aorta){tissue.keys[tissue.keys %in% "heart_and_aorta"] <- "heart"}
  if(anyNA(tissue.keys) || any(tissue.keys == "")){stop("Missing tissue annotation.")}
  cat("Analysis tissue cell counts:\n")
  print(table(tissue.keys))

  for(tissue.name in tissues) {
    cat("\nTissue:", tissue.name, "\n")
    selected <- colnames(mmf.rds)[tissue.keys == canonical.tissue(tissue.name)]
    if(length(selected) == 0L){stop("Requested tissue not found: ", tissue.name)}
    obj <- subset(mmf.rds, cells = selected)
    obj <- join.count.layer(obj)
    counts <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
    if(!inherits(counts, "dgCMatrix")){counts <- methods::as(counts, "dgCMatrix")}
    if(anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) || !setequal(colnames(counts), colnames(obj))){stop("Invalid raw-count names/coverage.")}
    counts <- counts[, colnames(obj), drop = FALSE]
    if(!gene %in% rownames(counts)){stop("Target gene missing from RNA counts: ", gene)}
    if(any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)){stop("Expected nonnegative integer RNA counts.")}
    meta <- obj@meta.data[colnames(counts), , drop = FALSE]
    age.text <- as.character(meta$age)
    if(anyNA(age.text) || any(!age.text %in% paste0(age.order, "m"))){stop("Expected age labels 1m, 3m, 18m, 21m, 24m, 30m in ", tissue.name)}
    age.months <- as.numeric(sub("m$", "", age.text))
    if(anyNA(meta$mouse.id) || any(trimws(as.character(meta$mouse.id)) == "")){stop("Missing mouse.id.")}
    mouse.age <- unique(data.frame(mouse.id = as.character(meta$mouse.id), age_months = age.months))
    if(anyDuplicated(mouse.age$mouse.id)){stop("One mouse.id maps to multiple ages.")}
    celltypes.raw <- as.character(meta$cell_ontology_class)
    celltypes.raw[is.na(celltypes.raw) | trimws(celltypes.raw) == ""] <- "Unannotated"
    cell.data <- data.frame(cell = colnames(counts), tissue = tissue.name, celltype = celltypes.raw, mouse.id = as.character(meta$mouse.id), age = age.text, age_months = age.months, technical_run = if("orig.ident" %in% names(meta)) as.character(meta$orig.ident) else NA_character_, sex = if("sex" %in% names(meta)) as.character(meta$sex) else NA_character_, tissue_annotation = as.character(meta[[tissue.column]]), stringsAsFactors = FALSE)
    sex.check <- cell.data %>% filter(!is.na(sex), trimws(sex) != "") %>% distinct(mouse.id, sex) %>% count(mouse.id)
    if(any(sex.check$n > 1L)){stop("One mouse.id maps to multiple sex labels.")}
    pb.meta <- cell.data %>% group_by(tissue, celltype, mouse.id, age, age_months) %>% summarise(n_cells = n(), technical_runs = join.values(technical_run), sex = join.values(sex), tissue_annotations = join.values(tissue_annotation), .groups = "drop") %>% arrange(celltype, age_months, mouse.id)
    pb.meta$pb_id <- sprintf("PB%05d", seq_len(nrow(pb.meta)))
    mapping <- left_join(cell.data[, c("cell", "celltype", "mouse.id")], pb.meta[, c("celltype", "mouse.id", "pb_id")], by = c("celltype", "mouse.id"))
    if(nrow(mapping) != ncol(counts) || anyNA(mapping$pb_id) || !identical(mapping$cell, colnames(counts))){stop("Cell-to-pseudobulk mapping failed.")}
    membership <- Matrix::sparseMatrix(i = seq_len(ncol(counts)), j = match(mapping$pb_id, pb.meta$pb_id), x = 1, dims = c(ncol(counts), nrow(pb.meta)), dimnames = list(colnames(counts), pb.meta$pb_id))
    pb.counts <- counts %*% membership
    expected <- rowsum(cbind(library_size = Matrix::colSums(counts), gene_count = as.numeric(counts[gene, ])), group = mapping$pb_id, reorder = FALSE)
    expected <- expected[pb.meta$pb_id, , drop = FALSE]
    if(any(abs(Matrix::rowSums(pb.counts) - Matrix::rowSums(counts)) > 1e-8) || any(abs(Matrix::colSums(pb.counts) - expected[, "library_size"]) > 1e-8) || any(abs(as.numeric(pb.counts[gene, ]) - expected[, "gene_count"]) > 1e-8)){stop("Raw-count aggregation checks failed.")}
    pb.meta$gene_raw_count <- as.numeric(pb.counts[gene, ])
    pb.meta$library_size <- as.numeric(Matrix::colSums(pb.counts))
    pb.meta$eligible <- pb.meta$n_cells >= min.cells.per.mouse.celltype & pb.meta$library_size > 0
    pb.meta$exclusion_reason <- ifelse(pb.meta$library_size <= 0, "zero_library_size", ifelse(pb.meta$n_cells < min.cells.per.mouse.celltype, "below_min_cells", "included"))
    root <- file.path(output.base, tissue.name, "age_regression")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    write.csv(pb.meta, file.path(root, "01.pseudobulk_inventory.csv"), row.names = FALSE, na = "NA")
    if(save.raw.pseudobulk){saveRDS(pb.counts, file.path(root, "02.pseudobulk_raw_counts.RDS"), compress = "gzip")}
    celltypes <- sort(unique(pb.meta$celltype))
    file.map <- data.frame(celltype = celltypes, file_stem = sprintf("%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100)), stringsAsFactors = FALSE)
    write.csv(file.map, file.path(root, "03.celltype_file_map.csv"), row.names = FALSE)
    rm(obj, counts, meta, cell.data, mapping, membership, expected)
    invisible(gc())

    # 4. 6개 연령 / 4개 연령을 각각 정규화하고 모형을 새로 적합합니다.
    results <- list()
    for(mode in names(analysis.ages)) {
      cat("  Analysis:", mode, "\n")
      out <- lapply(celltypes, function(ct) fit.one(pb.counts, pb.meta, ct, tissue.name, mode))
      stats <- bind_rows(lapply(out, function(z) z$stats))
      mice <- bind_rows(lapply(out, function(z) z$mice))
      curve <- bind_rows(lapply(out, function(z) z$curve))
      ok <- stats$status == "tested" & is.finite(stats$p_value)
      stats$p_adj_BH[ok] <- stats::p.adjust(stats$p_value[ok], method = "BH")
      stats$significant_BH <- ifelse(ok, stats$p_adj_BH < fdr.threshold, NA)
      stats$p_adjust_scope <- paste0("Successful Mff celltype trend tests within ", tissue.name, " / ", mode)
      stats$n_tests_in_BH_family <- sum(ok)
      folder <- file.path(root, mode)
      for(subfolder in c("plots", "overview")){dir.create(file.path(folder, subfolder), recursive = TRUE, showWarnings = FALSE)}
      write.csv(stats, file.path(folder, "01.Mff.age_trend_statistics.csv"), row.names = FALSE, na = "NA")
      write.csv(mice, file.path(folder, "02.Mff.mouse_values.csv"), row.names = FALSE, na = "NA")
      write.csv(curve, file.path(folder, "03.Mff.predicted_curve.csv"), row.names = FALSE, na = "NA")
      age.summary <- mice %>% group_by(celltype, age_months) %>% summarise(n_mice = n(), n_cells = sum(n_cells), n_mice_normalized = sum(is.finite(log2_CPM_plus1)), mean_log2_CPM_plus1 = if(any(is.finite(log2_CPM_plus1))) mean(log2_CPM_plus1[is.finite(log2_CPM_plus1)]) else NA_real_, sd_log2_CPM_plus1 = if(sum(is.finite(log2_CPM_plus1)) >= 2L) sd(log2_CPM_plus1[is.finite(log2_CPM_plus1)]) else NA_real_, .groups = "drop")
      grid <- expand.grid(celltype = celltypes, age_months = analysis.ages[[mode]], stringsAsFactors = FALSE)
      age.summary <- left_join(grid, age.summary, by = c("celltype", "age_months")) %>% arrange(celltype, age_months)
      for(nm in c("n_mice", "n_cells", "n_mice_normalized")){age.summary[[nm]][is.na(age.summary[[nm]])] <- 0L}
      age.summary$tissue <- tissue.name
      age.summary$analysis <- mode
      write.csv(age.summary, file.path(folder, "04.Mff.age_summary.csv"), row.names = FALSE, na = "NA")
      results[[mode]] <- list(stats = stats, mice = mice, curve = curve)
      print(as.data.frame(table(stats$status)), row.names = FALSE)
      if(any(nzchar(stats$model_warning) | nzchar(stats$normalization_warning))){cat("  Warnings captured in statistics CSV; inspect model_warning and normalization_warning.\n")}
      rm(out)
    }

    # 5. 같은 celltype의 두 분석에 동일한 X/Y축 범위를 적용합니다.
    # 숫자 X축은 3~18개월의 넓은 간격을 그대로 보여줍니다.
    compare.dir <- file.path(root, "comparison")
    dir.create(compare.dir, recursive = TRUE, showWarnings = FALSE)
    plots <- setNames(lapply(names(analysis.ages), function(z) list()), names(analysis.ages))
    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      values <- unlist(lapply(results, function(z) c(z$mice$log2_CPM_plus1[z$mice$celltype == ct], z$curve$predicted_log2_CPM_plus1[z$curve$celltype == ct])), use.names = FALSE)
      values <- values[is.finite(values)]
      ymax <- if(length(values)) max(plot.y.min + 1, max(values) * 1.32) else plot.y.min + 1
      for(mode in names(analysis.ages)) {
        p <- make.plot(results[[mode]], ct, tissue.name, mode, ymax)
        plots[[mode]][[i]] <- p
        save.plot(p, file.path(root, mode, "plots", paste0(file.map$file_stem[[i]], ".regression")))
      }
      pair <- patchwork::wrap_plots(lapply(names(analysis.ages), function(mode) plots[[mode]][[i]]), ncol = length(analysis.ages))
      save.plot(pair, file.path(compare.dir, paste0(file.map$file_stem[[i]], ".6ages_vs_4ages")), width = plot.width * length(analysis.ages))
      cat(sprintf("  Plot [%d/%d]: %s\n", i, length(celltypes), ct))
    }
    pages <- split(seq_along(celltypes), ceiling(seq_along(celltypes) / 9))
    for(mode in names(analysis.ages)) {
      for(page in seq_along(pages)) {
        idx <- pages[[page]]
        ncols <- min(3L, length(idx))
        nrows <- ceiling(length(idx) / ncols)
        p <- patchwork::wrap_plots(plots[[mode]][idx], ncol = ncols)
        save.plot(p, file.path(root, mode, "overview", sprintf("Mff.all_celltypes.page%02d", page)), width = plot.width * ncols, height = plot.height * nrows)
      }
    }
    cols <- c("celltype", "n_mice", "n_observed_ages", "observed_ages", "log2FC_per_6months", "fold_change_per_6months", "p_value", "p_adj_BH", "direction", "status", "reason")
    full <- results$full_6ages$stats[, cols, drop = FALSE]
    sens <- results$sensitivity_4ages$stats[, cols, drop = FALSE]
    comparison <- merge(full, sens, by = "celltype", all = TRUE, sort = FALSE, suffixes = c("_6ages", "_4ages"))
    comparison$tissue <- tissue.name
    both <- comparison$status_6ages == "tested" & comparison$status_4ages == "tested"
    comparison$same_direction <- ifelse(both, comparison$direction_6ages == comparison$direction_4ages, NA)
    comparison$delta_log2FC_per_6months_4ages_minus_6ages <- comparison$log2FC_per_6months_4ages - comparison$log2FC_per_6months_6ages
    write.csv(comparison, file.path(compare.dir, "00.Mff.6ages_vs_4ages.statistics.csv"), row.names = FALSE, na = "NA")

    # 설정, 해석상 주의점, 실행 환경을 함께 저장합니다.
    settings <- c(paste0("Input RDS (read only): ", rds.file), paste0("Tissue: ", tissue.name), paste0("Tissue metadata column: ", tissue.column), paste0("Include Heart_and_Aorta in Heart: ", include.heart.and.aorta), "Full analysis ages: 1,3,18,21,24,30 months", "Sensitivity ages: 3,18,21,24 months; excludes 1m and 30m", "No Young/Old grouping. Observations are different mice, not longitudinal repeated measures.", "Pseudobulk: raw UMI sum per tissue-celltype-mouse; technical runs combined.", paste0("Minimum cells per mouse-celltype: ", min.cells.per.mouse.celltype), "Only all-zero genes removed. No CPM/total-UMI gene expression filter.", "TMM is recomputed separately for each tissue-celltype-analysis; common mice may have slightly different plotted CPM values between analyses.", "Model: log(mu) = log(effective_library_size) + intercept + beta6 * (age_months - 3)/6", "edgeR v4 glmQLFit(legacy=FALSE, robust=TRUE, prior.count=0); NB dispersion estimated internally.", "All nonzero genes enter the count model; edgeR internally borrows information across genes to estimate variability.", "glmQLFTest tests age_per6 = 0; raw trend p-values shown on plots.", "log2FC_per_6months = beta6/log(2); log2FC_per_month is one sixth of this.", "The slope describes log2 expected relative expression, not the derivative of log2(CPM+1) at every age.", "Curve uses edgeR coefficients and normalized-library offsets; no ordinary lm() fit, R-squared or lm() confidence band.", "Points = log2(TMM-normalized CPM+1); curve = log2(predicted CPM+1).", "Curve spans observed ages only; gaps are model interpolation. Same X/Y limits in paired plots.", "BH family: successful Mff celltype tests separately within each tissue and each analysis.", paste0("BH significance threshold: ", fdr.threshold), "p=NA for no eligible mice, one age, zero residual df, all-zero Mff or a failed model/test. Inspect status/reason and warnings.", "The numeric-age model has two coefficients; residual df differs from a categorical-age model. A previously NA categorical-age result may be testable under the stronger numeric-age assumption.", "At least two observed ages and positive residual df are computational requirements, not guarantees of reliable inference.", "With only two observed ages, the slope summarizes that contrast; it does not establish a trajectory across several ages.", "No adjustment for sex, technical batch or other covariates. Sparse mouse-celltype pseudobulks remain included by default.", "Age-range restriction probes sensitivity; it does not remove all developmental or survivor-selection effects.", "The two analyses share mice. Agreement of signs is descriptive, not independent validation or a formal test of slope differences.", "Nonsignificant p-values do not prove absence of an age effect; interpret effect size and mouse-level scatter.", "No paper-specific Young/Old definition is claimed.", "Raw-count aggregation checks: PASSED", "", capture.output(sessionInfo()))
    settings <- c(settings, "", paste0("Plot Y-axis minimum: ", plot.y.min), "Y limits only change the displayed window. Values below the minimum remain in normalization, regression and CSVs.")
    writeLines(settings, file.path(root, "04.analysis_settings_and_sessionInfo.txt"))
    cat("Saved:", root, "\n")
    rm(pb.counts, pb.meta, results, plots)
    invisible(gc())
  }
  cat("\nDone: both tissues, full and sensitivity age regressions.\n")
})
