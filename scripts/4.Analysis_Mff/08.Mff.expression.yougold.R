# Young (1m,3m) vs Old (18m,21m,24m,30m), separately by tissue x cell type.
# Individual violin/boxplot: 900 x 850 px; p-value text size = 10 mm; p label only.
# Dotplot uses the SAME cells as the violin, including zero-expression cells.
# Linear normalized expression = expm1(existing RNA LogNormalize value), not raw UMI.
# Cell log2FC = log2((mean_Old + 1) / (mean_Young + 1)), using linear normalized means.
# A common configurable pseudocount is added AFTER each group mean.
# Seurat::FoldChange uses an explicit mean.fxn so this definition is version-independent.
# Dotplot significance reuses the violin's cellwise Wilcoxon raw p-value; no new test.
# Filled = unadjusted p < 0.05; open = unadjusted p >= 0.05; cross = test unavailable.
# Both raw p-values and BH-adjusted p-values remain in the CSV for reference.
# Effect estimates do not require successful p-values; missing-group effects remain NA.
# Dotplot output: young_old/cell_level_violin/log2FC_dotplot.
# The full-age cohort and the minimum-two-mice-per-group test rule are retained.
# Cell violin + internal box: existing RNA LogNormalize values, including zeros.
# Violin p-value: two-sided, unpaired, cell-level Wilcoxon rank-sum.
# Mouse boxplot: log2(TMM-normalized CPM+1), one point per mouse.
# Mouse p-value: two-sided, unpaired Wilcoxon rank-sum on the plotted log2(CPM + 1).
# edgeR is used only for TMM normalization and CPM, not for p-value calculation.
# Mouse Wilcoxon: exact if both groups have <50 mice and no ties; normal approximation otherwise.
# Raw UMI counts are summed by mouse x cell type, combining technical runs.
# Only all-zero genes are removed. No CPM/total-UMI gene-expression filter.
# Raw p-values displayed by default; BH across cell types per tissue AND test unit is retained in CSV only.
# Cell-level Wilcoxon does not adjust for dependence among cells from a mouse.
# Original RDS is read only. Figure subtitles and n= axis labels are absent.
# Run this complete file with source() in R.

local({
  #-------------------------------------------------------------------
  # 1. 설정

  gene <- "Mff"
  rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression"
  tissues <- c("Heart", "Limb_Muscle")
  tissue.column <- "tissue_free_annotation"  # Heart + Heart_and_Aorta를 Heart 분석 그룹으로 포함
  group.column <- "age_young_old"
  group.order <- c("Young", "Old")
  young.ages <- c("1m", "3m")
  old.ages <- c("18m", "21m", "24m", "30m")
  analysis.folder <- "young_old"
  group.colors <- c(Young = "#1B9E77", Old = "#D95F02")
  min.cells.per.mouse.celltype <- 1L
  min.mice.per.group.wilcoxon <- 2L  # 전체 연령 버전의 기존 최소 2마리 기준 유지
  min.cells.per.group.wilcoxon <- 2L  # 검정만 제한; violin/내부 box에는 모든 세포 포함
  show.pvalue.on.violin <- TRUE
  show.pvalue.on.boxplot <- TRUE
  pvalue.column <- "p_value"   # "p_value" = raw p, "p_adj_BH" = BH 보정값
  cell.fc.pseudocount <- 1        # dotplot: 로그를 푼 각 그룹 평균에 동일하게 더하는 값

  # 그림 크기와 글씨는 여기서 조절합니다 (크기 px, 글씨 pt).
  plot.width <- 900
  plot.height <- 850
  plot.res <- 160
  plot.base.size <- 13
  plot.title.size <- 20
  plot.axis.title.size <- 13
  plot.axis.text.size <- 30
  plot.title.wrap.width <- 40
  plot.violin.width <- 0.8
  plot.violin.box.width <- 0.18
  plot.violin.box.alpha <- 0.85
  plot.violin.box.linewidth <- 0.45
  plot.pvalue.size <- 10         # violin/boxplot p-value 글자 크기 (ggplot2 mm 단위)

  required <- c("Seurat", "SeuratObject", "edgeR", "Matrix", "dplyr", "ggplot2", "patchwork", "Cairo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing) > 0){stop("Missing packages: ", paste(missing, collapse = ", "))}
  if(packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0"){stop("Seurat and SeuratObject >= 5.0.0 are required.")}
  if(!pvalue.column %in% c("p_value", "p_adj_BH")){stop("pvalue.column must be p_value or p_adj_BH.")}
  if(length(cell.fc.pseudocount) != 1L || !is.finite(cell.fc.pseudocount) || cell.fc.pseudocount <= 0){stop("cell.fc.pseudocount must be a finite positive number.")}
  if(length(min.cells.per.group.wilcoxon) != 1L || !is.finite(min.cells.per.group.wilcoxon) || min.cells.per.group.wilcoxon < 2L || min.cells.per.group.wilcoxon != floor(min.cells.per.group.wilcoxon)){stop("min.cells.per.group.wilcoxon must be an integer >= 2.")}
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(length(min.cells.per.mouse.celltype) != 1L || !is.finite(min.cells.per.mouse.celltype) || min.cells.per.mouse.celltype < 1 || min.cells.per.mouse.celltype != floor(min.cells.per.mouse.celltype)){stop("min.cells.per.mouse.celltype must be an integer >= 1.")}
  if(length(min.mice.per.group.wilcoxon) != 1L || !is.finite(min.mice.per.group.wilcoxon) || min.mice.per.group.wilcoxon < 1L || min.mice.per.group.wilcoxon != floor(min.mice.per.group.wilcoxon)){stop("min.mice.per.group.wilcoxon must be an integer >= 1.")}
  if(!file.exists(rds.file)){stop("Missing input RDS: ", rds.file)}

  save.png <- function(plot, filename, width = plot.width, height = plot.height) {
    Cairo::CairoPNG(filename = filename, width = width, height = height, res = plot.res, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
  }
  join.values <- function(z) {
    z <- sort(unique(as.character(z[!is.na(z) & trimws(as.character(z)) != ""])))
    if(length(z) == 0){return(NA_character_)}
    paste(z, collapse = ";")
  }
  canonical.tissue <- function(z) {gsub("[[:space:]-]+", "_", tolower(trimws(as.character(z))))}
  plot.theme <- function() {
    theme_bw(base_size = plot.base.size) + theme(legend.position = "none", plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.title.x = element_text(size = plot.axis.title.size), axis.title.y = element_text(size = plot.axis.title.size), axis.text.x = element_text(size = plot.axis.text.size), axis.text.y = element_text(size = plot.axis.text.size), panel.grid.minor = element_blank())
  }
  join.checked.layer <- function(object, layer.name) {
    pattern <- paste0("^", layer.name, "($|\\.)")
    layers <- SeuratObject::Layers(object[["RNA"]], search = pattern)
    if(length(layers) == 0){stop("RNA ", layer.name, " layer is missing; no fallback is used.")}
    cells <- unlist(lapply(layers, function(z) colnames(SeuratObject::LayerData(object[["RNA"]], layer = z))), use.names = FALSE)
    if(anyDuplicated(cells) || !setequal(cells, colnames(object))){stop("Duplicated/missing cells in RNA ", layer.name, " layers.")}
    if(length(layers) != 1L || layers != layer.name) {
      if(!inherits(object[["RNA"]], "Assay5")){stop("Unexpected RNA layers: ", layer.name)}
      object <- SeuratObject::JoinLayers(object, assay = "RNA", layers = pattern, new = layer.name)
    }
    object
  }

  # Reference BH correction is separate for each tissue and test unit (cell or mouse).
  add.BH <- function(result, test.unit) {
    result$p_adj_BH <- NA_real_
    ok <- result$status == "tested" & is.finite(result$p_value)
    result$p_adj_BH[ok] <- stats::p.adjust(result$p_value[ok], method = "BH")
    result$p_adjust_scope <- paste0("All successful Young-vs-Old Mff celltype tests within this tissue; test unit=", test.unit)
    result$n_tests_in_BH_family <- sum(ok)
    result$analysis_id <- analysis.folder
    result$young_ages <- paste(young.ages, collapse = ";")
    result$old_ages <- paste(old.ages, collapse = ";")
    result
  }
  format.pvalue <- function(p) {
    prefix <- if(pvalue.column == "p_value") "p" else "BH p"
    if(!is.finite(p)){return(paste0(prefix, " = NA"))}
    if(p == 0){return(paste0(prefix, " < ", formatC(.Machine$double.xmin, format = "e", digits = 2)))}
    value <- if(p < 0.001) formatC(p, format = "e", digits = 2) else formatC(p, format = "f", digits = 3)
    paste0(prefix, " = ", value)
  }
  add.young.old.pvalue <- function(plot, stat.row, data.top) {
    stopifnot(nrow(stat.row) == 1L)
    label <- format.pvalue(stat.row[[pvalue.column]][[1]])
    plot + annotate("segment", x = 1, xend = 2, y = data.top * 1.10, yend = data.top * 1.10, linewidth = 0.4) + annotate("segment", x = c(1, 2), xend = c(1, 2), y = data.top * 1.065, yend = data.top * 1.10, linewidth = 0.4) + annotate("text", x = 1.5, y = data.top * 1.26, label = label, size = plot.pvalue.size, lineheight = 1.0)
  }

  # 1 cell = 1 observation for this requested exploratory cellwise test.
  # Same-mouse dependence is NOT adjusted by this Wilcoxon test.
  test.cell.young.old <- function(cell.data, celltypes, tissue.name) {
    result <- lapply(celltypes, function(ct) {
      d <- cell.data[cell.data$celltype == ct, , drop = FALSE]
      d1 <- d[d$group == "Young", , drop = FALSE]
      d2 <- d[d$group == "Old", , drop = FALSE]
      v1 <- d1$expression
      v2 <- d2$expression
      r <- data.frame(tissue = tissue.name, celltype = ct, gene = gene, group1 = "Young", group2 = "Old", n_cells_Young = length(v1), n_cells_Old = length(v2), n_mice_Young = length(unique(d1$mouse.id)), n_mice_Old = length(unique(d2$mouse.id)), mean_Young = if(length(v1)) mean(v1) else NA_real_, mean_Old = if(length(v2)) mean(v2) else NA_real_, median_Young = if(length(v1)) median(v1) else NA_real_, median_Old = if(length(v2)) median(v2) else NA_real_, W = NA_real_, probability_Old_higher_with_half_ties = NA_real_, p_value = NA_real_, method = "Wilcoxon rank-sum; two-sided; unpaired; normal approximation; tie and continuity correction", test_unit = "cell", value_tested = "Existing RNA LogNormalize expression, including zero", mouse_dependence_adjusted = FALSE, status = "not_tested", reason = "", model_warning = "", stringsAsFactors = FALSE)
      if(length(v1) == 0L || length(v2) == 0L){r$status <- "missing_group"; r$reason <- "No cells in Young or Old."; return(r)}
      if(length(v1) < min.cells.per.group.wilcoxon || length(v2) < min.cells.per.group.wilcoxon){r$status <- "insufficient_cells"; r$reason <- paste0("Each group needs at least ", min.cells.per.group.wilcoxon, " cells for this test; cells retained in plots."); return(r)}
      if(length(unique(c(v1, v2))) == 1L){r$status <- "all_values_identical"; r$reason <- "All pooled expression values are identical; rank variance is zero."; return(r)}
      warnings <- character()
      fit <- tryCatch(withCallingHandlers(stats::wilcox.test(v1, v2, alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE), warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
      r$model_warning <- paste(unique(warnings), collapse = " | ")
      if(inherits(fit, "error")){r$status <- "test_failed"; r$reason <- conditionMessage(fit); return(r)}
      if(length(fit$p.value) != 1L || !is.finite(fit$p.value) || fit$p.value < 0 || fit$p.value > 1){r$status <- "test_failed"; r$reason <- "Invalid Wilcoxon p-value."; return(r)}
      r$status <- "tested"
      r$p_value <- fit$p.value
      r$W <- unname(fit$statistic)
      r$probability_Old_higher_with_half_ties <- 1 - r$W / (as.double(length(v1)) * length(v2))
      r
    })
    add.BH(do.call(rbind, result), "cell; Wilcoxon")
  }

  # Violin과 동일한 세포로 평균 기반 log2FC를 계산한다. p-value 성공 여부와 독립.
  # mean(expm1(expression))은 로그를 푼 정규화 발현량의 평균이다.
  # Seurat 기본 mean.fxn 대신 명시적 함수를 지정하여 양쪽 평균에 같은 +1을 적용한다.
  add.cell.logfc <- function(result, cell.data) {
    result$mean_normalized_Young <- NA_real_
    result$mean_normalized_Old <- NA_real_
    result$FC_Old_vs_Young <- NA_real_
    result$logFC_Old_vs_Young <- NA_real_
    result$fc_pseudocount <- cell.fc.pseudocount
    result$fc_formula <- "log2((mean(expm1(Old)) + pseudocount) / (mean(expm1(Young)) + pseudocount))"
    result$fc_unit <- "cell"
    result$fc_weighting <- "Equal weight per cell; mice with more cells contribute more"
    result$fc_method <- "Seurat::FoldChange with explicit log2(mean(expm1(x)) + pseudocount) mean.fxn"
    result$fc_status <- "not_estimated"
    result$fc_reason <- ""
    for(j in seq_len(nrow(result))) {
      d <- cell.data[cell.data$celltype == result$celltype[j], , drop = FALSE]
      young.cells <- d$cell[d$group == "Young"]
      old.cells <- d$cell[d$group == "Old"]
      if(anyDuplicated(d$cell)){stop("Duplicated cell in fold-change input: ", result$celltype[j])}
      linear.expression <- expm1(d$expression)
      if(any(!is.finite(linear.expression)) || any(linear.expression < 0)){stop("Invalid linear normalized expression: ", result$celltype[j])}
      if(length(young.cells) > 0L){result$mean_normalized_Young[j] <- mean(linear.expression[d$group == "Young"])}
      if(length(old.cells) > 0L){result$mean_normalized_Old[j] <- mean(linear.expression[d$group == "Old"])}
      if(length(young.cells) == 0L || length(old.cells) == 0L){result$fc_status[j] <- "missing_group"; result$fc_reason[j] <- "No cells in Young or Old; no expression value imputed."; next}
      expression.matrix <- matrix(d$expression, nrow = 1L, dimnames = list(gene, d$cell))
      fc <- Seurat::FoldChange(object = expression.matrix, cells.1 = old.cells, cells.2 = young.cells, features = gene, mean.fxn = function(z) log2(rowMeans(expm1(z)) + cell.fc.pseudocount), fc.name = "logFC_Old_vs_Young")
      value <- fc[gene, "logFC_Old_vs_Young"]
      expected <- log2(result$mean_normalized_Old[j] + cell.fc.pseudocount) - log2(result$mean_normalized_Young[j] + cell.fc.pseudocount)
      if(length(value) != 1L || !is.finite(value) || !isTRUE(all.equal(as.numeric(value), expected, tolerance = 1e-10, check.attributes = FALSE))){stop("FoldChange did not match the stated linear-mean definition: ", result$celltype[j])}
      result$logFC_Old_vs_Young[j] <- as.numeric(value)
      result$FC_Old_vs_Young[j] <- 2^value
      result$fc_status[j] <- "estimated"
      result$fc_reason[j] <- "ok"
    }
    result
  }

  # One eligible independent mouse = one observation, exactly as in the boxplot.
  # Test the saved log2(TMM-normalized CPM + 1), without recomputing normalization.
  test.young.old <- function(pb.meta, celltypes, tissue.name) {
    result <- lapply(celltypes, function(ct) {
      d <- pb.meta[pb.meta$celltype == ct & pb.meta$eligible, , drop = FALSE]
      if(anyNA(d$mouse.id) || anyDuplicated(d$mouse.id)){stop("Missing/repeated mouse in statistical input: ", tissue.name, " / ", ct)}
      if(anyNA(d$group) || any(!d$group %in% group.order)){stop("Invalid mouse group: ", tissue.name, " / ", ct)}
      if(any(!is.finite(d$log2_CPM_plus1))){stop("Invalid plotted mouse expression: ", tissue.name, " / ", ct)}
      v1 <- d$log2_CPM_plus1[d$group == "Young"]
      v2 <- d$log2_CPM_plus1[d$group == "Old"]
      n1 <- length(v1)
      n2 <- length(v2)
      r <- data.frame(tissue = tissue.name, celltype = ct, gene = gene, group1 = "Young", group2 = "Old", n_mice_Young = n1, n_mice_Old = n2, n_cells_Young = sum(d$n_cells[d$group == "Young"]), n_cells_Old = sum(d$n_cells[d$group == "Old"]), mean_log2_CPM_plus1_Young = if(n1) mean(v1) else NA_real_, mean_log2_CPM_plus1_Old = if(n2) mean(v2) else NA_real_, median_log2_CPM_plus1_Young = if(n1) median(v1) else NA_real_, median_log2_CPM_plus1_Old = if(n2) median(v2) else NA_real_, W = NA_real_, p_value = NA_real_, method = "Wilcoxon rank-sum; two-sided; unpaired", test_unit = "mouse", value_tested = "Plotted log2(TMM-normalized CPM + 1), one value per eligible independent mouse", has_ties = anyDuplicated(c(v1, v2)) > 0L, exact_used = NA, continuity_correction = NA, pvalue_computation = NA_character_, status = "not_tested", reason = "", test_warning = "", stringsAsFactors = FALSE)
      r$single_mouse_in_group <- n1 == 1L || n2 == 1L
      r$inference_note <- if(r$single_mouse_in_group) "At least one group has a single independent mouse; that group's between-mouse variability cannot be checked directly. Interpret as exploratory." else ""
      if(n1 == 0L || n2 == 0L){r$status <- "missing_group"; r$reason <- "No eligible independent mouse in Young or Old."; return(r)}
      if(n1 < min.mice.per.group.wilcoxon || n2 < min.mice.per.group.wilcoxon){r$status <- "insufficient_mice"; r$reason <- paste0("Each group needs at least ", min.mice.per.group.wilcoxon, " independent mice for this analysis."); return(r)}
      if(length(unique(c(v1, v2))) == 1L){r$status <- "all_values_identical"; r$reason <- "All pooled mouse expression values are identical; rank variance is zero."; return(r)}
      # R 4.5-compatible choice: exact for small untied samples; normal approximation otherwise.
      use.exact <- n1 < 50L && n2 < 50L && !r$has_ties
      r$exact_used <- use.exact
      r$continuity_correction <- !use.exact
      r$pvalue_computation <- if(use.exact) "Exact rank-sum distribution; no ties" else "Normal approximation with tie and continuity correction"
      warnings <- character()
      fit <- tryCatch(withCallingHandlers(stats::wilcox.test(v1, v2, alternative = "two.sided", paired = FALSE, exact = use.exact, correct = !use.exact), warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}), error = function(e) e)
      r$test_warning <- paste(unique(warnings), collapse = " | ")
      if(inherits(fit, "error")){r$status <- "test_failed"; r$reason <- conditionMessage(fit); return(r)}
      if(length(fit$p.value) != 1L || !is.finite(fit$p.value) || fit$p.value < 0 || fit$p.value > 1){r$status <- "test_failed"; r$reason <- "Invalid mouse-level Wilcoxon p-value."; return(r)}
      r$status <- "tested"
      r$p_value <- fit$p.value
      r$W <- unname(fit$statistic)
      r
    })
    add.BH(do.call(rbind, result), "mouse; Wilcoxon")
  }

  #-------------------------------------------------------------------
  # 2. 최근 저장한 통합 RDS와 Young/Old 라벨 확인

  mmf.rds <- readRDS(rds.file)
  if(!inherits(mmf.rds, "Seurat") || !"RNA" %in% names(mmf.rds@assays)){stop("Expected a Seurat object with an RNA assay.")}
  needed <- unique(c("tissue", tissue.column, "age", group.column, "mouse.id", "cell_ontology_class"))
  if(!all(needed %in% colnames(mmf.rds@meta.data))){stop("Missing metadata: ", paste(setdiff(needed, colnames(mmf.rds@meta.data)), collapse = ", "))}
  DefaultAssay(mmf.rds) <- "RNA"
  tissue.values <- as.character(mmf.rds[[tissue.column, drop = TRUE]])
  tissue.keys <- canonical.tissue(tissue.values)
  # 추가한 3m Heart__10X_P7_4 (3-F-56)는 Heart_and_Aorta로 주석되어 있습니다.
  # 분석용 키에만 매핑하며 원본 metadata는 수정하지 않습니다.
  tissue.keys[tissue.keys %in% "heart_and_aorta"] <- "heart"
  cat("Input tissue and detailed tissue annotation (cell counts):\n")
  print(table(tissue = mmf.rds$tissue, annotation = tissue.values, useNA = "ifany"))
  if(anyNA(tissue.keys) || any(tissue.keys == "")){stop("Missing ", tissue.column, " metadata; inspect the annotation table above.")}
  absent <- tissues[!canonical.tissue(tissues) %in% tissue.keys]
  if(length(absent) > 0){stop("Requested tissue not found in ", tissue.column, ": ", paste(absent, collapse = ", "), "; observed values: ", paste(unique(tissue.values), collapse = ", "))}
  cat("Selected cells by analysis tissue (Heart includes Heart_and_Aorta):\n")
  print(table(factor(tissue.keys, levels = canonical.tissue(tissues), labels = tissues)))
  cat("Cells outside the selected tissue annotations:", sum(!tissue.keys %in% canonical.tissue(tissues)), "\n")

  for(tissue.name in tissues) {
    cat("\nTissue:", tissue.name, "\n")
    selected.cells <- colnames(mmf.rds)[tissue.keys == canonical.tissue(tissue.name)]
    tissue.rds <- subset(mmf.rds, cells = selected.cells)
    if(ncol(tissue.rds) == 0 || !gene %in% rownames(tissue.rds[["RNA"]])){stop("No cells or missing gene: ", tissue.name, " / ", gene)}
    tissue.rds <- join.checked.layer(tissue.rds, "counts")
    tissue.rds <- join.checked.layer(tissue.rds, "data")
    meta <- tissue.rds@meta.data[colnames(tissue.rds), , drop = FALSE]
    if(anyNA(meta$age) || any(!as.character(meta$age) %in% c(young.ages, old.ages))){stop("Missing or unexpected age: ", tissue.name)}
    expected.group <- ifelse(as.character(meta$age) %in% young.ages, "Young", "Old")
    if(anyNA(meta[[group.column]]) || any(as.character(meta[[group.column]]) != expected.group)){stop("Young/Old labels do not match the agreed age mapping: ", tissue.name)}
    if(anyNA(meta$mouse.id) || any(trimws(as.character(meta$mouse.id)) == "")){stop("Missing mouse.id: ", tissue.name)}
    mouse.age <- unique(meta[, c("mouse.id", "age"), drop = FALSE])
    if(anyDuplicated(as.character(mouse.age$mouse.id))){stop("A mouse.id maps to multiple ages: ", tissue.name)}
    tissue.rds[[group.column]] <- factor(as.character(meta[[group.column]]), levels = group.order)
    cat("Included mice by Young/Old and age (counts are cells):\n")
    included.mice <- data.frame(group = as.character(meta[[group.column]]), age = as.character(meta$age), mouse.id = as.character(meta$mouse.id), tissue_annotation = as.character(meta[[tissue.column]]), stringsAsFactors = FALSE) %>% count(group, age, mouse.id, tissue_annotation, name = "n_cells")
    print(as.data.frame(included.mice), row.names = FALSE)

    counts <- SeuratObject::LayerData(tissue.rds[["RNA"]], layer = "counts")
    if(!inherits(counts, "dgCMatrix")){counts <- methods::as(counts, "dgCMatrix")}
    if(anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) || !setequal(colnames(counts), colnames(tissue.rds))){stop("Invalid raw-count names/coverage: ", tissue.name)}
    counts <- counts[, colnames(tissue.rds), drop = FALSE]
    if(!gene %in% rownames(counts) || any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)){stop("Expected nonnegative integer RNA counts and the target gene: ", tissue.name)}
    data.genes <- rownames(SeuratObject::LayerData(tissue.rds[["RNA"]], layer = "data"))
    if(!gene %in% data.genes){stop("Target gene missing from normalized RNA data: ", tissue.name)}
    expr <- SeuratObject::FetchData(tissue.rds, vars = gene, layer = "data", clean = FALSE)
    if(!setequal(rownames(expr), colnames(tissue.rds)) || nrow(expr) != ncol(tissue.rds)){stop("Normalized expression does not cover all cells: ", tissue.name)}
    expr <- expr[colnames(tissue.rds), , drop = FALSE]
    if(any(!is.finite(expr[[gene]])) || any(expr[[gene]] < 0)){stop("Expected nonnegative RNA LogNormalize expression: ", tissue.name)}
    celltypes.raw <- as.character(meta$cell_ontology_class)
    celltypes.raw[is.na(celltypes.raw) | trimws(celltypes.raw) == ""] <- "Unannotated"
    cell.data <- data.frame(cell = colnames(tissue.rds), tissue = tissue.name, celltype = celltypes.raw, mouse.id = as.character(meta$mouse.id), age = factor(as.character(meta$age), levels = c(young.ages, old.ages)), group = factor(as.character(meta[[group.column]]), levels = group.order), sex = if("sex" %in% names(meta)) as.character(meta$sex) else NA_character_, technical_run = if("orig.ident" %in% names(meta)) as.character(meta$orig.ident) else NA_character_, gene = gene, expression = expr[[gene]], stringsAsFactors = FALSE)
    cell.data$tissue_original <- as.character(meta$tissue)
    cell.data$tissue_free_annotation <- as.character(meta[[tissue.column]])
    sex.check <- cell.data %>% filter(!is.na(sex), trimws(sex) != "") %>% distinct(mouse.id, sex) %>% count(mouse.id)
    if(any(sex.check$n > 1)){stop("A mouse.id maps to multiple sex labels: ", tissue.name)}
    celltypes <- sort(unique(cell.data$celltype))
    file.map <- data.frame(celltype = celltypes, file_stem = sprintf("%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100)), stringsAsFactors = FALSE)
    output.dir <- file.path(output.base, tissue.name, analysis.folder)
    violin.dir <- file.path(output.dir, "cell_level_violin")
    pb.dir <- file.path(output.dir, paste0("mouse_pseudobulk_TMM_min", min.cells.per.mouse.celltype))
    for(folder in c(file.path(violin.dir, c("violin", "overview")), file.path(pb.dir, c("boxplot", "overview")))){dir.create(folder, recursive = TRUE, showWarnings = FALSE)}
    write.csv(cell.data, file.path(violin.dir, "01.cell_expression.csv"), row.names = FALSE, na = "NA")
    write.csv(file.map, file.path(violin.dir, "03.celltype_file_map.csv"), row.names = FALSE)
    write.csv(file.map, file.path(pb.dir, "06.celltype_file_map.csv"), row.names = FALSE)

    #-----------------------------------------------------------------
    # 3. Cell-level summary와 Seurat violin (세포 수 cutoff 없음)

    cell.summary <- cell.data %>% group_by(celltype, group) %>% summarise(n_cells = n(), n_mice = n_distinct(mouse.id), mean_expression = mean(expression), median_expression = median(expression), q1 = unname(quantile(expression, 0.25)), q3 = unname(quantile(expression, 0.75)), min_expression = min(expression), max_expression = max(expression), percent_detected = mean(expression > 0) * 100, .groups = "drop")
    grid <- expand.grid(celltype = celltypes, group = group.order, stringsAsFactors = FALSE)
    grid$group <- factor(grid$group, levels = group.order)
    cell.summary <- left_join(grid, cell.summary, by = c("celltype", "group")) %>% arrange(celltype, group)
    cell.summary$n_cells[is.na(cell.summary$n_cells)] <- 0L
    cell.summary$n_mice[is.na(cell.summary$n_mice)] <- 0L
    cell.summary$tissue <- tissue.name
    cell.summary$display <- ifelse(cell.summary$n_cells == 0, "no_cells", ifelse(cell.summary$n_cells == 1, "single_cell_point", ifelse(cell.summary$min_expression == cell.summary$max_expression, "constant_value_line", "distribution")))
    write.csv(cell.summary, file.path(violin.dir, "02.celltype_young_old_summary.csv"), row.names = FALSE, na = "NA")
    cell.pvalue.results <- test.cell.young.old(cell.data, celltypes, tissue.name)
    cell.pvalue.results <- add.cell.logfc(cell.pvalue.results, cell.data)
    write.csv(cell.pvalue.results, file.path(violin.dir, paste0("05.", gene, ".Young_vs_Old.Wilcoxon.pvalues.csv")), row.names = FALSE, na = "NA")
    cat("Cellwise Wilcoxon status:\n")
    print(table(cell.pvalue.results$status))
    violin.plots <- list()
    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      d <- cell.data[cell.data$celltype == ct, , drop = FALSE]
      s <- cell.summary[cell.summary$celltype == ct, , drop = FALSE]
      data.top <- max(1, max(d$expression))
      ymax <- data.top * if(show.pvalue.on.violin) 1.42 else 1.08
      density.groups <- as.character(s$group[s$display == "distribution"])
      if(length(density.groups) > 0) {
        plot.rds <- subset(tissue.rds, cells = d$cell[as.character(d$group) %in% density.groups])
        plot.rds[[group.column]] <- factor(as.character(plot.rds[[group.column, drop = TRUE]]), levels = group.order[group.order %in% density.groups])
        p <- Seurat::VlnPlot(plot.rds, features = gene, assay = "RNA", layer = "data", group.by = group.column, cols = group.colors[levels(plot.rds[[group.column, drop = TRUE]])], pt.size = 0, add.noise = FALSE, y.max = ymax, combine = FALSE)[[1]]
        for(layer.i in seq_along(p$layers)) {
          if(inherits(p$layers[[layer.i]]$geom, "GeomViolin")) {
            p$layers[[layer.i]]$stat_params$width <- plot.violin.width
            p$layers[[layer.i]]$geom_params$width <- plot.violin.width
          }
        }
        rm(plot.rds)
      } else {
        p <- ggplot(d, aes(x = group, y = expression))
      }
      p <- p + geom_blank(data = d, aes(x = group, y = expression), inherit.aes = FALSE)
      # Internal boxplot also uses every LogNormalize value, including zero.
      # Constant and singleton groups are represented by the lines/points below.
      if(length(density.groups) > 0){p <- p + geom_boxplot(data = d[as.character(d$group) %in% density.groups, , drop = FALSE], mapping = aes(x = group, y = expression, group = group), inherit.aes = FALSE, position = "identity", width = plot.violin.box.width, fill = "white", colour = "black", alpha = plot.violin.box.alpha, linewidth = plot.violin.box.linewidth, outlier.shape = NA)}
      constant <- s[s$display == "constant_value_line", , drop = FALSE]
      single <- s[s$display == "single_cell_point", , drop = FALSE]
      if(nrow(constant) > 0){p <- p + geom_errorbar(data = constant, aes(x = group, ymin = min_expression, ymax = max_expression), inherit.aes = FALSE, width = 0.45, colour = "grey20")}
      if(nrow(single) > 0){p <- p + geom_point(data = single, aes(x = group, y = mean_expression), inherit.aes = FALSE, size = 2, colour = "black")}
      title <- paste(strwrap(paste0(tissue.name, " | ", ct), width = plot.title.wrap.width), collapse = "\n")
      p <- p + scale_x_discrete(limits = group.order, drop = FALSE) + coord_cartesian(ylim = c(0, ymax)) + labs(x = NULL, y = paste0(gene, " expression (RNA LogNormalize)"), title = title, subtitle = NULL, caption = NULL) + plot.theme()
      if(show.pvalue.on.violin){p <- add.young.old.pvalue(p, cell.pvalue.results[cell.pvalue.results$celltype == ct, , drop = FALSE], data.top)}
      save.png(p, file.path(violin.dir, "violin", paste0(file.map$file_stem[[i]], ".Young_Old.violin.png")))
      violin.plots[[i]] <- p
    }

    #-----------------------------------------------------------------
    # 4. Mouse x celltype별 raw-count pseudobulk

    pb.meta <- cell.data %>% group_by(tissue, celltype, mouse.id, age, group) %>% summarise(n_cells = n(), sex = join.values(sex), technical_runs = join.values(technical_run), tissue_original = join.values(tissue_original), tissue_free_annotation = join.values(tissue_free_annotation), .groups = "drop") %>% arrange(celltype, group, age, mouse.id)
    pb.meta$pb_id <- sprintf("PB%05d", seq_len(nrow(pb.meta)))
    cell.map <- left_join(cell.data[, c("cell", "tissue", "celltype", "mouse.id", "age", "group")], pb.meta[, c("celltype", "mouse.id", "pb_id")], by = c("celltype", "mouse.id"))
    if(nrow(cell.map) != ncol(counts) || anyNA(cell.map$pb_id) || !identical(cell.map$cell, colnames(counts))){stop("Cell-to-pseudobulk mapping failed: ", tissue.name)}
    tissue.rds$pseudobulk_id <- setNames(cell.map$pb_id, cell.map$cell)[colnames(tissue.rds)]
    pb.counts <- Seurat::AggregateExpression(tissue.rds, assays = "RNA", features = rownames(counts), group.by = "pseudobulk_id", return.seurat = FALSE, verbose = FALSE)[["RNA"]]
    if(!setequal(rownames(pb.counts), rownames(counts)) || !setequal(colnames(pb.counts), pb.meta$pb_id)){stop("Unexpected pseudobulk dimensions/names: ", tissue.name)}
    pb.counts <- pb.counts[rownames(counts), pb.meta$pb_id, drop = FALSE]
    expected <- rowsum(cbind(library_size = Matrix::colSums(counts), gene_count = as.numeric(counts[gene, ])), group = cell.map$pb_id, reorder = FALSE)
    expected <- expected[pb.meta$pb_id, , drop = FALSE]
    if(any(abs(Matrix::rowSums(pb.counts) - Matrix::rowSums(counts)) > 1e-8) || any(abs(Matrix::colSums(pb.counts) - expected[, "library_size"]) > 1e-8) || any(abs(as.numeric(pb.counts[gene, ]) - expected[, "gene_count"]) > 1e-8)){stop("Raw-count aggregation checks failed: ", tissue.name)}
    pb.meta$gene <- gene
    pb.meta$gene_raw_count <- as.numeric(pb.counts[gene, ])
    pb.meta$library_size <- as.numeric(Matrix::colSums(pb.counts))
    pb.meta$eligible <- pb.meta$n_cells >= min.cells.per.mouse.celltype & pb.meta$library_size > 0
    pb.meta$exclusion_reason <- ifelse(pb.meta$library_size <= 0, "zero_library_size", ifelse(pb.meta$n_cells < min.cells.per.mouse.celltype, "below_min_cells", "included"))
    pb.meta$norm_factor <- NA_real_
    pb.meta$effective_library_size <- NA_real_
    pb.meta$normalized_CPM <- NA_real_
    pb.meta$log2_CPM_plus1 <- NA_real_
    pb.meta$normalization_method <- "not_normalized_excluded"
    write.csv(cell.map, file.path(pb.dir, "01.cell_to_pseudobulk.csv"), row.names = FALSE)
    saveRDS(pb.counts, file.path(pb.dir, "02.pseudobulk_raw_counts.RDS"), compress = "gzip")
    rm(tissue.rds, counts, meta, expr, expected)
    invisible(gc())

    # 각 celltype에서 Young/Old의 모든 eligible mouse를 함께 정규화합니다.
    for(ct in celltypes) {
      idx <- which(pb.meta$celltype == ct & pb.meta$eligible)
      if(length(idx) == 0){next}
      mat <- as.matrix(pb.counts[, pb.meta$pb_id[idx], drop = FALSE])
      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      y <- edgeR::DGEList(counts = mat)
      if(length(idx) >= 2) {
        y <- edgeR::calcNormFactors(y, method = "TMM")
        method <- "TMM"
      } else {
        y$samples$norm.factors <- 1
        method <- "CPM_only_single_mouse"
      }
      factors <- y$samples$norm.factors
      effective.lib <- y$samples$lib.size * factors
      if(any(!is.finite(effective.lib)) || any(effective.lib <= 0)){stop("Invalid effective library size: ", tissue.name, " / ", ct)}
      cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
      gene.cpm <- if(gene %in% rownames(cpm)) as.numeric(cpm[gene, ]) else rep(0, length(idx))
      if(any(!is.finite(gene.cpm)) || any(gene.cpm < 0)){stop("Invalid normalized CPM: ", tissue.name, " / ", ct)}
      pb.meta$norm_factor[idx] <- factors
      pb.meta$effective_library_size[idx] <- effective.lib
      pb.meta$normalized_CPM[idx] <- gene.cpm
      pb.meta$log2_CPM_plus1[idx] <- log2(gene.cpm + 1)
      pb.meta$normalization_method[idx] <- method
      rm(mat, y, cpm)
    }
    write.csv(pb.meta, file.path(pb.dir, paste0("03.", gene, ".pseudobulk_mouse.csv")), row.names = FALSE, na = "NA")
    pb.summary <- pb.meta %>% group_by(celltype, group) %>% summarise(n_mice_total = n(), n_mice = sum(eligible), n_cells_total = sum(n_cells), n_cells_included = sum(n_cells[eligible]), mean_log2_CPM_plus1 = if(any(eligible)) mean(log2_CPM_plus1[eligible]) else NA_real_, sd_log2_CPM_plus1 = if(sum(eligible) > 1) sd(log2_CPM_plus1[eligible]) else NA_real_, .groups = "drop")
    pb.summary <- left_join(grid, pb.summary, by = c("celltype", "group")) %>% arrange(celltype, group)
    for(nm in c("n_mice_total", "n_mice", "n_cells_total", "n_cells_included")){pb.summary[[nm]][is.na(pb.summary[[nm]])] <- 0L}
    pb.summary$tissue <- tissue.name
    write.csv(pb.summary, file.path(pb.dir, "04.pseudobulk_young_old_summary.csv"), row.names = FALSE, na = "NA")
    mouse.inventory <- cell.data %>% group_by(tissue, mouse.id, age, group) %>% summarise(n_cells = n(), sex = join.values(sex), technical_runs = join.values(technical_run), tissue_original = join.values(tissue_original), tissue_free_annotation = join.values(tissue_free_annotation), .groups = "drop") %>% arrange(group, age, mouse.id)
    write.csv(mouse.inventory, file.path(pb.dir, "05.mouse_inventory.csv"), row.names = FALSE, na = "NA")

    #-----------------------------------------------------------------
    # 5. Young vs Old p-value를 별도 CSV로 저장

    pvalue.results <- test.young.old(pb.meta, celltypes, tissue.name)
    write.csv(pvalue.results, file.path(pb.dir, paste0("08.", gene, ".Young_vs_Old.pvalues.csv")), row.names = FALSE, na = "NA")
    cat("Mouse pseudobulk Wilcoxon rank-sum status:\n")
    print(table(pvalue.results$status))

    #-----------------------------------------------------------------
    # 6. Mouse-level boxplot: 모든 eligible mouse를 점으로 표시

    box.plots <- list()
    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      d <- pb.meta[pb.meta$celltype == ct & pb.meta$eligible, , drop = FALSE]
      s <- pb.summary[pb.summary$celltype == ct, , drop = FALSE]
      data.top <- if(nrow(d) > 0) max(1, max(d$log2_CPM_plus1)) else 1
      ymax <- data.top * if(show.pvalue.on.boxplot) 1.42 else 1.12
      box.groups <- as.character(s$group[s$n_mice >= 2])
      p <- ggplot(d, aes(x = group, y = log2_CPM_plus1))
      if(any(as.character(d$group) %in% box.groups)){p <- p + geom_boxplot(data = d[as.character(d$group) %in% box.groups, , drop = FALSE], width = 0.55, outlier.shape = NA, fill = "grey92")}
      p <- p + geom_point(aes(color = group), position = position_jitter(width = 0.08, height = 0, seed = 1234), size = 3) + scale_color_manual(values = group.colors, drop = FALSE) + scale_x_discrete(limits = group.order, drop = FALSE)
      stat.row <- pvalue.results[pvalue.results$celltype == ct, , drop = FALSE]
      if(show.pvalue.on.boxplot){p <- add.young.old.pvalue(p, stat.row, data.top)}
      title <- paste(strwrap(paste0(tissue.name, " | ", ct), width = plot.title.wrap.width), collapse = "\n")
      p <- p + coord_cartesian(ylim = c(0, ymax)) + labs(x = NULL, y = paste0(gene, " pseudobulk log2(normalized CPM + 1)"), title = title, subtitle = NULL, caption = NULL) + plot.theme()
      save.png(p, file.path(pb.dir, "boxplot", paste0(file.map$file_stem[[i]], ".Young_Old.mouse_boxplot.png")))
      box.plots[[i]] <- p
    }

    # Tissue마다 최대 9개 celltype을 한 페이지로 모은 overview도 저장합니다.
    pages <- split(seq_along(celltypes), ceiling(seq_along(celltypes) / 9))
    for(page in seq_along(pages)) {
      idx <- pages[[page]]
      ncols <- min(3L, length(idx))
      nrows <- ceiling(length(idx) / ncols)
      save.png(patchwork::wrap_plots(violin.plots[idx], ncol = ncols), file.path(violin.dir, "overview", sprintf("%s.Young_Old.violin.page%02d.png", gene, page)), width = plot.width * ncols, height = plot.height * nrows)
      save.png(patchwork::wrap_plots(box.plots[idx], ncol = ncols), file.path(pb.dir, "overview", sprintf("%s.Young_Old.mouse_boxplot.page%02d.png", gene, page)), width = plot.width * ncols, height = plot.height * nrows)
    }
    settings <- c(paste0("Input RDS: ", rds.file), paste0("Tissue: ", tissue.name), paste0("Tissue selection metadata: ", tissue.column), "Selection: Heart + Heart_and_Aorta; Limb_Muscle separately; original annotations retained", paste0("Gene: ", gene), paste0("Grouping metadata: ", group.column), "Young: 1m,3m; Old: 18m,21m,24m,30m", "Cell violin and internal boxplot: existing RNA data layer (LogNormalize), including zero-expression cells; no plotting cell-count cutoff", "Internal box: median and Q1-Q3, 1.5-IQR whiskers; outlier dots hidden", "Cell test: stats::wilcox.test, two-sided, unpaired, exact=FALSE, correct=TRUE; tie-corrected normal approximation", paste0("Minimum cells per group FOR WILCOXON ONLY: ", min.cells.per.group.wilcoxon), "Cell Wilcoxon does not adjust dependence among cells from the same mouse and is not mouse-level inference", "Pseudobulk: raw UMI sum per tissue-celltype-mouse; same-mouse technical runs combined", paste0("Minimum cells per mouse-celltype: ", min.cells.per.mouse.celltype), "TMM: all eligible mice jointly within each tissue-celltype; only all-zero genes removed", "No CPM/total-UMI gene expression filters", "Single eligible mouse: normalization factor=1, descriptive CPM only", "Mouse plot: log2(TMM-normalized CPM+1); one point per eligible mouse", "Mouse test: stats::wilcox.test on the plotted log2(TMM-normalized CPM+1), one observation per eligible mouse; two-sided and unpaired", "Mouse Wilcoxon: exact=TRUE when both groups have fewer than 50 mice and no pooled ties; otherwise exact=FALSE with tie and continuity correction", "Mouse Wilcoxon compares ranks; no sex/batch covariates", "edgeR is used only for TMM normalization and CPM, not statistical testing", paste0("Minimum eligible independent mice per group for mouse Wilcoxon: ", min.mice.per.group.wilcoxon), "Single-mouse groups are flagged in single_mouse_in_group and inference_note", "Minimum counts of mice/cells are computational criteria, not guarantees of reliable inference", "Reference BH correction: successful Mff celltype tests within each tissue, separately for cell-level and mouse-level Wilcoxon; plots use raw p by default", paste0("Displayed p-value column: ", pvalue.column), "Cell-test NA: missing group, too few cells, identical pooled values or failed test", "Mouse-test NA: missing group, too few eligible mice, identical pooled values, or failed Wilcoxon test", "NA values are shown on plots and explained in each statistics CSV", "Young/Old pools the listed ages; this is not a six-age omnibus test or an age-trend test", "TMM measures relative expression, not absolute transcripts per cell", paste0("Show p-values on violin: ", show.pvalue.on.violin), paste0("Show p-values on mouse boxplot: ", show.pvalue.on.boxplot), "Raw-count aggregation checks: PASSED", "Original RDS read only; no normalization/PCA/UMAP rerun", "", capture.output(sessionInfo()))
    settings <- c(settings, paste0("Individual violin/boxplot size: ", plot.width, " x ", plot.height, " px at ", plot.res, " dpi"), paste0("P-value annotation: numeric p label only; size=", plot.pvalue.size, " mm; Young-Old bracket retained"), "Cell dotplot input: the same cells and existing RNA LogNormalize values as the violin, including zeros", "Cell dotplot: linear normalized expression = expm1(LogNormalize); this does not recover raw UMI counts", paste0("Cell log2FC formula: log2((mean_Old + ", cell.fc.pseudocount, ")/(mean_Young + ", cell.fc.pseudocount, ")); means are of linear normalized expression"), "FoldChange: explicit mean.fxn adds a common pseudocount after group averaging; no version-dependent default mean function", "Cell log2FC is not a ratio/difference of mean logged violin values and does not compare medians", "Equal cell weights; groups are not balanced by mouse or old age", "Cell dotplot p-values: existing cellwise Wilcoxon results; shapes use unadjusted p_value; BH values remain in CSV for reference", "Cell log2FC available with p=NA: plot the effect as a cross; a missing group has FC=NA", "Both groups all zero: corrected FC=1 and log2FC=0; Wilcoxon p=NA", "Cell log2FC fields added to 05.Mff.Young_vs_Old.Wilcoxon.pvalues.csv", "Cell dotplots are written under cell_level_violin/log2FC_dotplot; previously saved mouse edgeR dotplots are not updated")
    writeLines(settings, file.path(pb.dir, "07.analysis_settings_and_sessionInfo.txt"))
    writeLines(settings, file.path(violin.dir, "04.analysis_settings_and_sessionInfo.txt"))
    cat("Saved:", output.dir, "\n")
    rm(pb.counts, pb.meta, cell.data, cell.map, violin.plots, box.plots)
    invisible(gc())
  }
  cat("\nDone: Young vs Old Mff analysis for both tissues.\n")

  #-------------------------------------------------------------------
  # 7. Violin과 같은 세포의 평균 log2FC 및 Wilcoxon raw p-value 점그래프
  # Young=1m,3m / Old=18m,21m,24m,30m. 위에서 저장한 cell-level CSV 사용.
  # 입력/출력 경로, 조직, 유전자, 연령, pseudocount 설정은 위 분석과 공유한다.

plot.mff.young.old.logfc <- function() {
  # ------------------------------------------------------------------
  # 1. 입력 / 출력 설정
  result.base <- output.base
  source.folder <- "cell_level_violin"
  input.filename <- paste0("05.", gene, ".Young_vs_Old.Wilcoxon.pvalues.csv")
  output.folder <- "log2FC_dotplot"
  pvalue.cutoff <- 0.05  # dotplot: 보정 전 Wilcoxon p-value 기준

  # 그림 설정: 두 조직에 같은 x축 범위를 사용합니다.
  celltype.order <- "logFC"  # "logFC": 감소가 큰 순서부터 / "alphabetical": 이름순
  x.limit.manual <- NULL  # 자동 설정. 예: 1로 바꾸면 -1 ~ +1; 범위 밖 점이 있으면 중단
  plot.width <- 10.5  # inch; 긴 세포형 이름이 있으면 늘리세요.
  plot.height.minimum <- 5.5
  plot.height.per.celltype <- 0.48
  plot.dpi <- 300
  plot.title.size <- 22
  plot.axis.title.size <- 18
  plot.axis.number.size <- 17
  plot.celltype.text.size <- 16
  plot.legend.size <- 13
  plot.point.size <- 4.5
  plot.point.stroke <- 1.2
  color.old.lower <- "#D95F02"
  color.old.higher <- "#000000"
  color.no.change <- "#777777"
  show.logFC.numbers <- FALSE  # TRUE: 각 점 옆에 log2FC 숫자 표시
  plot.logFC.number.size <- 4  # geom_text의 크기 단위는 mm

  if(!requireNamespace("ggplot2", quietly = TRUE)){stop("ggplot2 패키지를 설치해 주세요.")}
  if(!celltype.order %in% c("logFC", "alphabetical")){stop("celltype.order must be logFC or alphabetical.")}
  if(length(pvalue.cutoff) != 1L || !is.finite(pvalue.cutoff) || pvalue.cutoff <= 0 || pvalue.cutoff >= 1){stop("pvalue.cutoff must be between 0 and 1.")}
  if(length(tissues) == 0L || anyNA(tissues) || anyDuplicated(tissues)){stop("Invalid tissues setting.")}
  input.files <- setNames(file.path(result.base, tissues, analysis.folder, source.folder, input.filename), tissues)
  missing.files <- input.files[!file.exists(input.files)]
  if(length(missing.files)){stop("먼저 위의 Young-Old cell violin 및 평균 log2FC 계산 코드를 실행해 주세요. 없는 파일:\n", paste(missing.files, collapse = "\n"))}

  # ------------------------------------------------------------------
  # 2. 전체 세포형 결과 읽기 / 입력 검증
  required.columns <- c("tissue", "celltype", "gene", "group1", "group2", "logFC_Old_vs_Young", "FC_Old_vs_Young", "mean_normalized_Young", "mean_normalized_Old", "fc_pseudocount", "fc_unit", "fc_status", "fc_reason", "p_value", "p_adj_BH", "method", "test_unit", "status", "reason", "analysis_id", "young_ages", "old_ages")
  results <- setNames(vector("list", length(tissues)), tissues)
  for(tissue.name in tissues) {
    d <- utils::read.csv(input.files[[tissue.name]], stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
    missing.columns <- setdiff(required.columns, names(d))
    if(length(missing.columns)){stop("Cell violin/log2FC 결과 열이 없습니다: ", paste(missing.columns, collapse = ", "), "\n파일: ", input.files[[tissue.name]], "\n이 전체 스크립트에서 새로 저장한 cell-level 결과를 사용하세요.")}
    if(nrow(d) == 0L){stop("Empty result table: ", tissue.name)}
    if(anyNA(d$analysis_id) || any(d$analysis_id != analysis.folder) || anyNA(d$young_ages) || any(d$young_ages != paste(young.ages, collapse = ";")) || anyNA(d$old_ages) || any(d$old_ages != paste(old.ages, collapse = ";"))){stop("Full-age analysis or age selection mismatch in input CSV: ", tissue.name)}
    if(anyNA(d$celltype) || any(trimws(d$celltype) == "") || anyDuplicated(d$celltype)){stop("Missing or duplicated celltype: ", tissue.name)}
    if(anyNA(d$tissue) || any(d$tissue != tissue.name) || anyNA(d$gene) || any(d$gene != gene)){stop("Tissue or gene mismatch: ", tissue.name)}
    if(anyNA(d$group1) || anyNA(d$group2) || any(d$group1 != "Young") || any(d$group2 != "Old")){stop("Expected group1=Young and group2=Old: ", tissue.name)}
    if(anyNA(d$test_unit) || any(d$test_unit != "cell") || anyNA(d$method) || any(!grepl("Wilcoxon rank-sum", d$method, fixed = TRUE))){stop("Expected the violin's cell-level Wilcoxon results: ", tissue.name)}
    if(anyNA(d$fc_unit) || any(d$fc_unit != "cell") || anyNA(d$fc_status)){stop("Invalid cell-level fold-change metadata: ", tissue.name)}
    if(anyNA(d$status)){stop("Missing test status: ", tissue.name)}
    for(nm in c("logFC_Old_vs_Young", "FC_Old_vs_Young", "mean_normalized_Young", "mean_normalized_Old", "fc_pseudocount", "p_value", "p_adj_BH")) {
      old <- d[[nm]]
      d[[nm]] <- suppressWarnings(as.numeric(old))
      if(any(!is.na(old) & is.na(d[[nm]]))){stop("Non-numeric value in ", nm, ": ", tissue.name)}
    }
    if(any(!is.finite(d$fc_pseudocount)) || any(d$fc_pseudocount != cell.fc.pseudocount)){stop("Fold-change pseudocount mismatch: ", tissue.name)}
    tested <- d$status == "tested"
    if(any(tested & (!is.finite(d$p_value) | d$p_value < 0 | d$p_value > 1))){stop("Invalid tested p-value: ", tissue.name)}
    # p_adj_BH는 CSV에 보존하며, dotplot에서는 p_value만 판정에 사용합니다.
    d$estimate_available <- d$fc_status == "estimated" & is.finite(d$logFC_Old_vs_Young)
    if(any(d$fc_status == "estimated" & !d$estimate_available)){stop("Nonfinite estimated cell log2FC: ", tissue.name)}
    e <- d[d$estimate_available, , drop = FALSE]
    expected.fc <- log2(e$mean_normalized_Old + e$fc_pseudocount) - log2(e$mean_normalized_Young + e$fc_pseudocount)
    if(any(!is.finite(expected.fc)) || !isTRUE(all.equal(e$logFC_Old_vs_Young, expected.fc, tolerance = 1e-7, check.attributes = FALSE))){stop("Cell log2FC does not match the saved linear means: ", tissue.name)}
    d$plot_log2FC <- ifelse(d$estimate_available, d$logFC_Old_vs_Young, NA_real_)
    d$plot_direction <- ifelse(!d$estimate_available, "Not estimated", ifelse(d$plot_log2FC < 0, "Old lower", ifelse(d$plot_log2FC > 0, "Old higher", "No change")))
    d$plot_significance <- ifelse(!d$estimate_available, NA_character_, ifelse(!tested | !is.finite(d$p_value), "Not tested", ifelse(d$p_value < pvalue.cutoff, "Below cutoff", "At or above cutoff")))
    d$plot_pvalue_column <- "p_value"
    d$plot_pvalue_cutoff <- pvalue.cutoff
    d$plot_note <- ifelse(!d$estimate_available, d$fc_reason, ifelse(!tested, paste0("FC available; Wilcoxon p=NA: ", d$reason), ""))
    idx <- if(celltype.order == "logFC") order(!d$estimate_available, d$plot_log2FC, d$celltype, na.last = TRUE) else order(d$celltype)
    d <- d[idx, , drop = FALSE]
    d$plot_order_top_to_bottom <- seq_len(nrow(d))
    results[[tissue.name]] <- d
  }

  # 모든 추정값을 포함하는 대칭 범위. 조직별 CSV는 따로 두고 축 범위만 공유합니다.
  all.effects <- as.numeric(unlist(lapply(results, function(d) d$plot_log2FC[d$estimate_available]), use.names = FALSE))
  max.effect <- max(c(0, abs(all.effects)))
  x.limit <- max(0.25, max.effect * if(show.logFC.numbers) 1.4 else 1.15)
  if(!is.null(x.limit.manual)) {
    if(length(x.limit.manual) != 1L || !is.finite(x.limit.manual) || x.limit.manual <= 0){stop("x.limit.manual must be NULL or a positive number.")}
    if(x.limit.manual < max.effect){stop("x.limit.manual would hide a logFC estimate. Use at least ", signif(max.effect, 5), ".")}
    x.limit <- x.limit.manual
  }
  x.breaks <- pretty(c(-x.limit, x.limit), n = 5)
  x.breaks <- sort(unique(c(0, x.breaks[x.breaks >= -x.limit & x.breaks <= x.limit])))
  direction.colors <- c("Old lower" = color.old.lower, "Old higher" = color.old.higher, "No change" = color.no.change)
  significance.labels <- c("Below cutoff" = paste0("Unadjusted p < ", pvalue.cutoff), "At or above cutoff" = paste0("Unadjusted p >= ", pvalue.cutoff), "Not tested" = "p = NA")

  # ------------------------------------------------------------------
  # 3. 조직별 PNG / PDF / 표시 데이터 CSV 저장
  for(tissue.name in tissues) {
    d <- results[[tissue.name]]
    d$celltype_axis <- factor(d$celltype, levels = rev(d$celltype))
    estimated <- d[d$estimate_available, , drop = FALSE]
    unavailable <- d[!d$estimate_available, , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(y = celltype_axis)) + ggplot2::geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.6, linetype = "dashed") + ggplot2::geom_point(data = estimated, ggplot2::aes(x = plot_log2FC, colour = plot_direction, shape = plot_significance), size = plot.point.size, stroke = plot.point.stroke) + ggplot2::scale_colour_manual(values = direction.colors, guide = "none") + ggplot2::scale_shape_manual(values = c("Below cutoff" = 16, "At or above cutoff" = 1, "Not tested" = 4), limits = names(significance.labels), labels = unname(significance.labels), drop = FALSE, name = NULL) + ggplot2::scale_x_continuous(breaks = x.breaks, limits = c(-x.limit, x.limit), expand = ggplot2::expansion(mult = 0.02)) + ggplot2::scale_y_discrete(drop = FALSE, expand = ggplot2::expansion(add = 0.7)) + ggplot2::labs(title = paste0(gsub("_", " ", tissue.name), " | ", gene), subtitle = NULL, x = "Cell-mean log2 fold change (Old / Young)", y = NULL) + ggplot2::theme_classic(base_size = 14, base_family = "sans") + ggplot2::theme(plot.title = ggplot2::element_text(size = plot.title.size, face = "bold", margin = ggplot2::margin(b = 14)), axis.title.x = ggplot2::element_text(size = plot.axis.title.size, margin = ggplot2::margin(t = 12)), axis.text.x = ggplot2::element_text(size = plot.axis.number.size, colour = "black"), axis.text.y = ggplot2::element_text(size = plot.celltype.text.size, colour = "black", margin = ggplot2::margin(r = 10)), axis.ticks.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.35), legend.position = "bottom", legend.text = ggplot2::element_text(size = plot.legend.size), legend.key.width = grid::unit(1.0, "cm"), plot.margin = ggplot2::margin(15, 25, 12, 12)) + ggplot2::guides(shape = ggplot2::guide_legend(nrow = 1, override.aes = list(colour = "#444444", size = 4)))
    # NA 행도 세포형 이름을 유지하며, 점을 0에 배치하지 않습니다.
    if(nrow(unavailable) > 0L){p <- p + ggplot2::geom_text(data = unavailable, ggplot2::aes(y = celltype_axis), x = x.limit * 0.96, label = "NA", inherit.aes = FALSE, hjust = 1, colour = "#888888", size = 4.5)}
    if(show.logFC.numbers && nrow(estimated) > 0L){p <- p + ggplot2::geom_text(data = estimated, ggplot2::aes(x = plot_log2FC, y = celltype_axis, label = sprintf("%+.2f", plot_log2FC), hjust = ifelse(plot_log2FC < 0, 1.3, -0.3)), inherit.aes = FALSE, size = plot.logFC.number.size, colour = "#444444", show.legend = FALSE)}
    output.dir <- file.path(result.base, tissue.name, analysis.folder, source.folder, output.folder)
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    plot.height <- max(plot.height.minimum, 2.0 + nrow(d) * plot.height.per.celltype)
    stem <- paste0(gene, ".Young_vs_Old.cell_mean.log2FC")
    ggplot2::ggsave(filename = file.path(output.dir, paste0(stem, ".png")), plot = p, width = plot.width, height = plot.height, units = "in", dpi = plot.dpi, bg = "white", limitsize = FALSE)
    pdf.device <- if(isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(filename = file.path(output.dir, paste0(stem, ".pdf")), plot = p, device = pdf.device, width = plot.width, height = plot.height, units = "in", bg = "white", limitsize = FALSE)
    d$celltype_axis <- NULL
    utils::write.csv(d, file.path(output.dir, paste0(stem, ".plot_data.csv")), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    notes <- c(paste0("Source: ", input.files[[tissue.name]]), paste0("Tissue: ", tissue.name), paste0("Analysis: ", analysis.folder), paste0("Young: ", paste(young.ages, collapse = ","), "; Old: ", paste(old.ages, collapse = ",")), "Full-age analysis: Young=1m,3m; Old=18m,21m,24m,30m; no sensitivity age exclusions", "Effect: same cells as violin; linear normalized values = expm1(existing RNA LogNormalize)", paste0("Definition: log2((mean_Old + ", cell.fc.pseudocount, ")/(mean_Young + ", cell.fc.pseudocount, ")); common pseudocount added after averaging"), "Means give every cell equal weight, including zero-expression cells; not equal weight per mouse or age", "Seurat::FoldChange with explicit mean.fxn; this is not Seurat's version-dependent default pseudocount placement", "Positive: higher cell-mean normalized Mff expression in Old; negative: lower in Old", "This is not a fold change of logged means, a median ratio, or an edgeR model estimate", "Point: one tissue-celltype Mff expression contrast, not one mouse or cell", "Red: Old lower; black: Old higher; grey: estimated logFC exactly zero", paste0("Filled: cellwise Wilcoxon unadjusted p < ", pvalue.cutoff, "; open: unadjusted p >= ", pvalue.cutoff, "; cross: FC available but p=NA"), "Dotplot significance: raw p_value, without multiple-testing correction; not an FDR threshold", "For reference only: p_adj_BH across successful cellwise Mff celltype tests within each tissue is retained in the CSV", "Violin displays the selected pvalue.column (raw p by default); dotplot shapes use raw Wilcoxon p_value", "The cellwise Wilcoxon test does not adjust same-mouse cell dependence; not mouse-level inference", "FC and test availability are separate: insufficient cells for testing can still yield an FC", "All-zero expression in both groups: corrected FC=1/log2FC=0 with p=NA", "No cells in either group: FC=NA, celltype retained with NA text and no point", "FC and test reasons are saved separately in plot_data.csv", "Direction alone does not establish statistical significance; low-expression effects depend on the chosen pseudocount", "No confidence intervals are calculated", "No new normalization or tests in the dotplot step; use the results just saved by this script", paste0("Shared x-axis limits: ", -x.limit, " to ", x.limit), paste0("Celltype order: ", celltype.order), "", capture.output(utils::sessionInfo()))
    writeLines(notes, file.path(output.dir, paste0(stem, ".settings.txt")), useBytes = TRUE)
    cat("\nSaved:", output.dir, "\n")
    cat("Cell types:", nrow(d), "| plotted estimates:", nrow(estimated), "| NA:", nrow(unavailable), "| Unadjusted p <", pvalue.cutoff, ":", sum(estimated$plot_significance == "Below cutoff", na.rm = TRUE), "\n")
  }
  invisible(NULL)
}

plot.mff.young.old.logfc()
})
