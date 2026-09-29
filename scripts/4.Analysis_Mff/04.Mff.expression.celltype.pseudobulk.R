

local({
  #-------------------------------------------------------------------
  # 1. 설정

  gene <- "Mff"
  min.cells.per.mouse.celltype <- 1L  # 1로 바꾸면 모든 관측 그룹을 사용
  age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
  age.months <- as.numeric(sub("m$", "", age.order))
  age.colors <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"), age.order)
  rds.files <- c(Heart = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Heart/Aging.MFF.Heart.seurat.normalization.pca.umap.RDS", 
    Limb_Muscle = "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/4.tissue_split/Limb_Muscle/Aging.MFF.Limb_Muscle.seurat.normalization.pca.umap.RDS")
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression"
  output.folder <- paste0("mouse_pseudobulk_TMM_min", min.cells.per.mouse.celltype)

  # 그림 글씨 크기 (pt): 필요하면 아래 값만 변경한 뒤 전체 코드를 다시 실행
  # 부제목과 오른쪽 아래 각주는 표시하지 않습니다.
  # 기존 output directory와 PNG 파일명을 그대로 사용하여 덮어씁니다.
  plot.base.size <- 13
  plot.title.size <- 20          # 그림 제목
  plot.axis.title.x.size <- 13   # X축 제목: Age (months)
  plot.axis.title.y.size <- 13   # Y축 제목: Mff pseudobulk ...
  plot.axis.text.x.size <- 20    # X축 눈금: age (boxplot에는 n=mouse 수 포함)
  plot.axis.text.y.size <- 20    # Y축 눈금 숫자
  plot.title.wrap.width <- 48    # 제목 줄바꿈 기준 글자 수
  plot.pvalue.size <- 10        # boxplot 우측 상단 p-value 글씨 크기 (ggplot2 mm 단위)

  # CPM cutoff 및 UMI 합계 cutoff를 사용하지 않는다.
  # 전체 유전자는 TMM 정규화에 사용하며, 검정 입력은 mouse별 Mff log2(CPM+1)이다.

  required <- c("Seurat", "SeuratObject", "edgeR", "Matrix", "dplyr", "ggplot2", "patchwork", "Cairo")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing) > 0){stop("Missing packages: ", paste(missing, collapse = ", "), ". See installation notes at the top of this script.")}
  if(packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0"){stop("Seurat and SeuratObject >= 5.0.0 are required.")}
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  if(length(min.cells.per.mouse.celltype) != 1L || !is.finite(min.cells.per.mouse.celltype) || min.cells.per.mouse.celltype < 1 || min.cells.per.mouse.celltype != floor(min.cells.per.mouse.celltype)){stop("min.cells.per.mouse.celltype must be a positive integer.")}
  if(any(!file.exists(rds.files))){stop("Missing RDS: ", paste(rds.files[!file.exists(rds.files)], collapse = ", "))}

  save.png <- function(plot, filename, width = 1500, height = 1100) {
    Cairo::CairoPNG(filename = filename, width = width, height = height, res = 160, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
  }
  join.values <- function(z) {
    z <- sort(unique(as.character(z[!is.na(z) & trimws(as.character(z)) != ""])))
    if(length(z) == 0){return(NA_character_)}
    paste(z, collapse = ";")
  }

  #-------------------------------------------------------------------
  # 1b. tissue x cell type당 mouse-level Kruskal-Wallis 검정 1회
  # 그림에 사용한 정규화 발현값으로 관측된 모든 age를 함께 비교한다.
  # 연령은 범주형 factor로 사용하며 sex/batch 보정은 하지 않는다.

  run.kruskal.omnibus <- function(d, tissue.name, ct) {
    age.n <- table(factor(as.character(d$age), levels = age.order))
    observed.ages <- age.order[age.n > 0]
    n.age <- length(observed.ages)
    out <- data.frame(tissue = tissue.name, celltype = ct, gene = gene, method = "Kruskal-Wallis rank sum", test_unit = "mouse.id", 
      test_expression = "TMM normalized Mff log2(CPM+1)", pvalue_method = "Tie-corrected H; asymptotic chi-squared", comparison = "All observed age groups jointly",
       n_mice = nrow(d), n_age_groups = n.age, observed_ages = paste(observed.ages, collapse = ";"), missing_ages = paste(age.order[age.n == 0], collapse = ";"),
        n_singleton_age_groups = sum(age.n == 1), any_age_fewer_than_5_mice = any(age.n > 0 & age.n < 5), Mff_total_count = sum(d$gene_raw_count), H = NA_real_,
         df_test = NA_real_, PValue = NA_real_, FDR_BH_within_tissue = NA_real_, status = "not_tested", reason = "", warnings = "", stringsAsFactors = FALSE)
    for(a in age.order){out[[paste0("n_mice_", a)]] <- as.integer(age.n[a])}
    collected.warnings <- character(0)
    worker <- function() {
      r <- out
      if(nrow(d) == 0){r$reason <- "no_eligible_mice"; return(r)}
      if(anyDuplicated(d$mouse.id) || anyDuplicated(d$pb_id)){stop("Independent-mouse mapping is inconsistent.")}
      if(anyNA(d$age) || any(!as.character(d$age) %in% age.order)){stop("Missing or unexpected age in test input.")}
      if(any(!is.finite(d$log2_CPM_plus1)) || any(d$log2_CPM_plus1 < 0)){stop("Expected finite nonnegative normalized Mff values.")}
      if(n.age < 2){r$reason <- "fewer_than_two_observed_age_groups"; return(r)}
      if(nrow(d) <= n.age){r$reason <- "no_within_age_replication"; return(r)}
      if(sum(d$gene_raw_count) == 0){r$reason <- "target_all_zero"; return(r)}
      if(length(unique(d$log2_CPM_plus1)) < 2L){r$reason <- "all_normalized_values_identical"; return(r)}

      # 실제 검정: one normalized Mff value per independent mouse.
      test <- stats::kruskal.test(x = d$log2_CPM_plus1, g = factor(as.character(d$age), levels = observed.ages))
      r$H <- unname(test$statistic)
      r$df_test <- unname(test$parameter)
      r$PValue <- test$p.value
      if(!is.finite(r$PValue) || r$PValue < 0 || r$PValue > 1 || !is.finite(r$H)){r$PValue <- NA_real_; r$status <- "failed"; r$reason <- "nonfinite_or_invalid_test_result"; return(r)}
      if(!is.finite(r$df_test) || r$df_test != n.age - 1L){stop("Unexpected Kruskal-Wallis degrees of freedom.")}
      r$status <- "tested"
      r$reason <- "ok"
      r
    }
    result <- withCallingHandlers(tryCatch(worker(), error = function(e) {out$status <- "failed"; out$reason <- paste0("KruskalWallis_error: ", conditionMessage(e)); out}), 
      warning = function(w) {collected.warnings <<- c(collected.warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    if(n.age >= 2 && any(age.n > 0 & age.n < 5)){collected.warnings <- c(collected.warnings, "Small observed age groups: chi-squared approximation may be inaccurate.")}
    if(n.age >= 2 && any(age.n == 1)){collected.warnings <- c(collected.warnings, "Some observed age groups contain only one independent mouse.")}
    if(n.age > 0 && n.age < length(age.order)){collected.warnings <- c(collected.warnings, "The omnibus test covers observed ages only; missing ages were not imputed.")}
    result$warnings <- paste(unique(collected.warnings), collapse = " | ")
    result
  }

  format.test.p <- function(p) {
    if(length(p) != 1L || !is.finite(p)){return("NA")}
    if(p == 0){return(paste0("< ", formatC(.Machine$double.xmin, format = "e", digits = 2)))}
    if(p < 0.001){return(formatC(p, format = "e", digits = 2))}
    formatC(p, format = "f", digits = 3)
  }
  make.test.label <- function(r) {
    if(r$status != "tested" || !is.finite(r$PValue)){return("p = NA")}
    if(r$PValue == 0){return(paste("p", format.test.p(r$PValue)))}
    paste0("p = ", format.test.p(r$PValue))
  }

  for(tissue.name in names(rds.files)) {
    cat("\nTissue:", tissue.name, "\n")
    x <- readRDS(rds.files[[tissue.name]])
    if(!inherits(x, "Seurat") || !"RNA" %in% names(x@assays)){stop("Expected a Seurat object with an RNA assay: ", tissue.name)}
    DefaultAssay(x) <- "RNA"
    if(ncol(x) == 0 || !gene %in% rownames(x[["RNA"]])){stop("No cells or missing ", gene, ": ", tissue.name)}
    needed <- c("age", "mouse.id", "cell_ontology_class")
    if(!all(needed %in% colnames(x@meta.data))){stop("Missing metadata: ", paste(setdiff(needed, colnames(x@meta.data)), collapse = ", "))}

    #-----------------------------------------------------------------
    # 2. 원본 counts와 metadata 확인

    count.layers <- SeuratObject::Layers(x[["RNA"]], search = "^counts($|\\.)")
    if(length(count.layers) == 0){stop("RNA raw counts layer is missing: ", tissue.name)}
    if(length(count.layers) != 1L || count.layers != "counts") {
      if(!inherits(x[["RNA"]], "Assay5")){stop("Unexpected counts layers: ", tissue.name)}
      layer.cells <- unlist(lapply(count.layers, function(z) colnames(SeuratObject::LayerData(x[["RNA"]], layer = z))), use.names = FALSE)
      if(anyDuplicated(layer.cells) || !setequal(layer.cells, colnames(x))){stop("Split counts layers have duplicated or missing cells: ", tissue.name)}
      x <- SeuratObject::JoinLayers(x, assay = "RNA", layers = "^counts($|\\.)", new = "counts")
      rm(layer.cells)
    }
    counts <- SeuratObject::LayerData(x[["RNA"]], layer = "counts")
    if(!inherits(counts, "dgCMatrix")){counts <- methods::as(counts, "dgCMatrix")}
    if(anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) || !setequal(colnames(counts), colnames(x))){stop("Counts names or cell coverage are inconsistent: ", tissue.name)}
    if(any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8)){stop("Expected nonnegative integer raw UMI counts: ", tissue.name)}
    if(!gene %in% rownames(counts)){stop(gene, " is missing from raw counts: ", tissue.name)}
    meta <- x@meta.data[colnames(counts), , drop = FALSE]
    if(anyNA(meta$age) || any(!as.character(meta$age) %in% age.order)){stop("Missing or unexpected age: ", tissue.name)}
    if(anyNA(meta$mouse.id) || any(trimws(as.character(meta$mouse.id)) == "")){stop("Missing mouse.id: ", tissue.name)}
    mouse.age <- unique(meta[, c("mouse.id", "age"), drop = FALSE])
    if(anyDuplicated(as.character(mouse.age$mouse.id))){stop("A mouse.id maps to multiple ages: ", tissue.name)}
    ct <- as.character(meta$cell_ontology_class)
    ct[is.na(ct) | trimws(ct) == ""] <- "Unannotated"
    cell.meta <- data.frame(cell = colnames(counts), tissue = tissue.name, celltype = ct, mouse.id = as.character(meta$mouse.id), age = factor(as.character(meta$age), levels = age.order), sex = if("sex" %in% names(meta)) as.character(meta$sex) else NA_character_, orig.ident = if("orig.ident" %in% names(meta)) as.character(meta$orig.ident) else NA_character_, stringsAsFactors = FALSE)
    sex.check <- cell.meta %>% filter(!is.na(sex), trimws(sex) != "") %>% distinct(mouse.id, sex) %>% count(mouse.id)
    if(any(sex.check$n > 1)){stop("A mouse.id maps to multiple sex labels: ", tissue.name)}

    # 합산용 ID를 따로 만들고 원래 mouse/cell type 이름은 대응표에 보관한다.
    pb.meta <- cell.meta %>% group_by(tissue, celltype, mouse.id, age) %>% summarise(n_cells = n(), sex = join.values(sex), technical_runs = join.values(orig.ident), .groups = "drop") %>% arrange(celltype, age, mouse.id)
    pb.meta$pb_id <- sprintf("PB%05d", seq_len(nrow(pb.meta)))
    cell.meta <- left_join(cell.meta, pb.meta[, c("celltype", "mouse.id", "pb_id")], by = c("celltype", "mouse.id"))
    if(nrow(cell.meta) != ncol(counts) || anyNA(cell.meta$pb_id) || !identical(cell.meta$cell, colnames(counts))){stop("Cell-to-pseudobulk mapping failed: ", tissue.name)}
    x$pseudobulk_id <- setNames(cell.meta$pb_id, cell.meta$cell)[colnames(x)]

    #-----------------------------------------------------------------
    # 3. 모든 유전자의 raw counts를 mouse x cell type별로 합산
    # return.seurat = FALSE: 정규화되지 않은 합산 count matrix를 반환.

    pb.counts <- Seurat::AggregateExpression(x, assays = "RNA", features = rownames(counts), group.by = "pseudobulk_id", return.seurat = FALSE, verbose = FALSE)[["RNA"]]
    if(!setequal(rownames(pb.counts), rownames(counts)) || !setequal(colnames(pb.counts), pb.meta$pb_id)){stop("Unexpected pseudobulk dimensions/names: ", tissue.name)}
    pb.counts <- pb.counts[rownames(counts), pb.meta$pb_id, drop = FALSE]

    # 합산 전후에 유전자별 전체 count, 그룹별 library size, Mff 합계 확인.
    expected <- rowsum(cbind(library_size = Matrix::colSums(counts), gene_count = as.numeric(counts[gene, ])), group = cell.meta$pb_id, reorder = FALSE)
    expected <- expected[pb.meta$pb_id, , drop = FALSE]
    if(any(abs(Matrix::rowSums(pb.counts) - Matrix::rowSums(counts)) > 1e-8) || any(abs(Matrix::colSums(pb.counts) - expected[, "library_size"]) > 1e-8) || any(abs(as.numeric(pb.counts[gene, ]) - expected[, "gene_count"]) > 1e-8)){stop("Raw-count aggregation checks failed: ", tissue.name)}
    pb.meta$age_months <- as.numeric(sub("m$", "", as.character(pb.meta$age)))
    pb.meta$library_size <- as.numeric(Matrix::colSums(pb.counts))
    pb.meta$gene <- gene
    pb.meta$gene_raw_count <- as.numeric(pb.counts[gene, ])
    pb.meta$eligible <- pb.meta$n_cells >= min.cells.per.mouse.celltype & pb.meta$library_size > 0
    pb.meta$exclusion_reason <- ifelse(pb.meta$library_size <= 0, "zero_library_size", ifelse(pb.meta$n_cells < min.cells.per.mouse.celltype, "below_min_cells", "included"))
    pb.meta$norm_factor <- NA_real_
    pb.meta$effective_library_size <- NA_real_
    pb.meta$normalized_CPM <- NA_real_
    pb.meta$log2_CPM_plus1 <- NA_real_
    pb.meta$normalization_method <- "not_normalized_excluded"
    rm(x, counts, meta, expected)
    invisible(gc())

    output.dir <- file.path(output.base, tissue.name, output.folder)
    for(folder in c("boxplot", "age_trend", "overview")){dir.create(file.path(output.dir, folder), recursive = TRUE, showWarnings = FALSE)}
    write.csv(cell.meta, file.path(output.dir, "01.cell_to_pseudobulk.csv"), row.names = FALSE, na = "NA")
    saveRDS(pb.counts, file.path(output.dir, "02.pseudobulk_raw_counts.RDS"), compress = "gzip")

    #-----------------------------------------------------------------
    # 4. 각 tissue x cell type 내부에서 TMM 정규화
    # Mff 하나로 정규화하지 않는다. 전체 유전자 중 모든 샘플에서 0인 행만
    # 계산에서 제외한다. Mff가 모든 샘플에서 0이면 CPM도 그대로 0이다.
    # Normalized CPM = raw gene count / (library_size * norm_factor) * 1e6
    # Plot value = log2(normalized CPM + 1)
    # CPM에 +1을 넣으므로 edgeR::cpm(log=TRUE)의 기본 prior.count와 다르다.

    celltypes <- sort(unique(pb.meta$celltype))
    test.rows <- setNames(vector("list", length(celltypes)), celltypes)
    for(ct in celltypes) {
      idx <- which(pb.meta$celltype == ct & pb.meta$eligible)
      if(length(idx) == 0){test.rows[[ct]] <- run.kruskal.omnibus(pb.meta[idx, , drop = FALSE], tissue.name, ct); cat("  No eligible mouse groups:", ct, "\n"); next}
      mat <- as.matrix(pb.counts[, pb.meta$pb_id[idx], drop = FALSE])
      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      y <- edgeR::DGEList(counts = mat)
      if(length(idx) >= 2) {
        y <- edgeR::calcNormFactors(y, method = "TMM")
        method <- "TMM"
      } else {
        y$samples$norm.factors <- 1
        method <- "CPM_only_single_mouse"
        cat("  Only one eligible mouse; factor=1, descriptive CPM only:", ct, "\n")
      }
      factors <- y$samples$norm.factors
      effective.lib <- y$samples$lib.size * factors
      if(any(!is.finite(effective.lib)) || any(effective.lib <= 0)){stop("Invalid effective library size: ", tissue.name, " / ", ct)}
      cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
      gene.cpm <- if(gene %in% rownames(cpm)) as.numeric(cpm[gene, ]) else rep(0, length(idx))
      if(any(!is.finite(gene.cpm)) || any(gene.cpm < 0)){stop("Invalid CPM: ", tissue.name, " / ", ct)}
      pb.meta$norm_factor[idx] <- factors
      pb.meta$effective_library_size[idx] <- effective.lib
      pb.meta$normalized_CPM[idx] <- gene.cpm
      pb.meta$log2_CPM_plus1[idx] <- log2(gene.cpm + 1)
      pb.meta$normalization_method[idx] <- method
      test.rows[[ct]] <- run.kruskal.omnibus(pb.meta[idx, , drop = FALSE], tissue.name, ct)
      test.result <- test.rows[[ct]]
      cat(sprintf("  Kruskal-Wallis %s | %s: %s, p=%s; %d age groups, %d mice\n", tissue.name, ct, test.result$reason, format.test.p(test.result$PValue), test.result$n_age_groups, test.result$n_mice))
      if(nzchar(test.result$warnings)){cat("    Note: ", test.result$warnings, "\n", sep = "")}
      cat(sprintf("  %s: %d / %d mouse groups included\n", ct, length(idx), sum(pb.meta$celltype == ct)))
      rm(mat, y, cpm)
    }
    write.csv(pb.meta, file.path(output.dir, paste0("03.", gene, ".pseudobulk_mouse.csv")), row.names = FALSE, na = "NA")

    # Mff만을 미리 정한 표적 유전자로 보고, tissue 내 검정한 cell type끼리 BH 보정.
    # 같은 tissue에서 검정된 cell type들의 Mff p-value만 보정한다.
    test.summary <- dplyr::bind_rows(test.rows)
    tested <- test.summary$status == "tested" & is.finite(test.summary$PValue)
    test.summary$FDR_BH_within_tissue[tested] <- stats::p.adjust(test.summary$PValue[tested], method = "BH")
    test.summary$n_tests_BH_within_tissue <- sum(tested)
    test.summary$FDR_scope <- paste0("Mff across tested cell types within ", tissue.name)
    test.summary$plot_label <- vapply(seq_len(nrow(test.summary)), function(j) make.test.label(test.summary[j, , drop = FALSE]), character(1))
    write.csv(test.summary, file.path(output.dir, paste0("08.", gene, ".KruskalWallis_omnibus.csv")), row.names = FALSE, na = "NA")

    # 나이별 평균은 mouse의 log2(CPM+1) 값을 동일 가중치로 평균낸다.
    # 관측되지 않은 age는 n_mice=0, 발현값=NA로 남긴다.
    age.summary <- pb.meta %>% group_by(tissue, celltype, age) %>% summarise(n_mice_total = n(), n_mice = sum(eligible), n_cells_total = sum(n_cells), n_cells_included = sum(n_cells[eligible]), mean_log2_CPM_plus1 = if(any(eligible)) mean(log2_CPM_plus1[eligible]) else NA_real_, sd_log2_CPM_plus1 = if(sum(eligible) > 1) sd(log2_CPM_plus1[eligible]) else NA_real_, .groups = "drop")
    grid <- expand.grid(celltype = celltypes, age = age.order, stringsAsFactors = FALSE)
    grid$age <- factor(grid$age, levels = age.order)
    age.summary <- left_join(grid, age.summary, by = c("celltype", "age")) %>% arrange(celltype, age)
    for(nm in c("n_mice_total", "n_mice", "n_cells_total", "n_cells_included")){age.summary[[nm]][is.na(age.summary[[nm]])] <- 0L}
    age.summary$tissue <- tissue.name
    age.summary$age_months <- as.numeric(sub("m$", "", as.character(age.summary$age)))
    age.summary$sem_log2_CPM_plus1 <- ifelse(age.summary$n_mice > 1, age.summary$sd_log2_CPM_plus1 / sqrt(age.summary$n_mice), NA_real_)
    write.csv(age.summary, file.path(output.dir, paste0("04.", gene, ".pseudobulk_age_summary.csv")), row.names = FALSE, na = "NA")
    mouse.inventory <- cell.meta %>% group_by(tissue, mouse.id, age) %>% summarise(n_cells = n(), sex = join.values(sex), technical_runs = join.values(orig.ident), .groups = "drop") %>% arrange(age, mouse.id)
    write.csv(mouse.inventory, file.path(output.dir, "05.mouse_inventory.csv"), row.names = FALSE, na = "NA")
    file.map <- data.frame(celltype = celltypes, file_stem = sprintf("%03d.%s", seq_along(celltypes), substr(gsub("[^A-Za-z0-9_-]+", "_", celltypes), 1, 100)), stringsAsFactors = FALSE)
    write.csv(file.map, file.path(output.dir, "06.celltype_file_map.csv"), row.names = FALSE)

    #-----------------------------------------------------------------
    # 5. Cell type별 boxplot + mouse 점 / 실제 개월 수에 따른 추이

    box.plots <- list()
    trend.plots <- list()
    for(i in seq_along(celltypes)) {
      ct <- celltypes[[i]]
      d <- pb.meta[pb.meta$celltype == ct & pb.meta$eligible, , drop = FALSE]
      s <- age.summary[age.summary$celltype == ct, , drop = FALSE]
      age.labels <- setNames(paste0(as.character(s$age), "\nn=", s$n_mice), as.character(s$age))
      ymax <- if(nrow(d) > 0) max(1, max(d$log2_CPM_plus1) * 1.12) else 1
      box.ymax <- ymax * 1.25
      test.label <- test.summary$plot_label[match(ct, test.summary$celltype)]
      title <- paste(strwrap(paste0(tissue.name, " | ", ct), width = plot.title.wrap.width), collapse = "\n")
      common <- list(scale_color_manual(values = age.colors, drop = FALSE), labs(y = paste0(gene, " pseudobulk log2(normalized CPM + 1)"), title = title, subtitle = NULL, caption = NULL), theme_bw(base_size = plot.base.size), theme(legend.position = "none", plot.title = element_text(size = plot.title.size, face = "bold"), plot.subtitle = element_blank(), plot.caption = element_blank(), axis.title.x = element_text(size = plot.axis.title.x.size), axis.title.y = element_text(size = plot.axis.title.y.size), axis.text.x = element_text(size = plot.axis.text.x.size), axis.text.y = element_text(size = plot.axis.text.y.size), panel.grid.minor = element_blank()))

      # 같은 age의 mouse가 2마리 이상일 때 box를 그리고, 모든 mouse는 점으로 표시.
      box.ages <- as.character(s$age[s$n_mice >= 2])
      p.box <- ggplot(d, aes(x = age, y = log2_CPM_plus1))
      if(any(as.character(d$age) %in% box.ages)){p.box <- p.box + geom_boxplot(data = d[as.character(d$age) %in% box.ages, , drop = FALSE], width = 0.55, outlier.shape = NA, fill = "grey92")}
      p.box <- p.box + geom_point(aes(color = age), position = position_jitter(width = 0.08, height = 0, seed = 1234), size = 3) + scale_x_discrete(limits = age.order, labels = age.labels, drop = FALSE) + labs(x = "Age (months)") + common + coord_cartesian(ylim = c(0, box.ymax))
      p.box <- p.box + annotate("text", x = Inf, y = Inf, label = test.label, hjust = 1.03, vjust = 1.15, size = plot.pvalue.size, colour = "black", lineheight = 1.05)

      p.trend <- ggplot(d, aes(x = age_months, y = log2_CPM_plus1))
      if(sum(s$n_mice > 0) >= 2){p.trend <- p.trend + geom_line(data = s, aes(x = age_months, y = mean_log2_CPM_plus1, group = 1), inherit.aes = FALSE, colour = "grey45", na.rm = TRUE)}
      p.trend <- p.trend + geom_point(aes(color = age), position = position_jitter(width = 0.15, height = 0, seed = 1234), size = 2.8) + geom_point(data = s, aes(x = age_months, y = mean_log2_CPM_plus1), inherit.aes = FALSE, shape = 18, size = 4, colour = "black", na.rm = TRUE)
      sem.data <- s[s$n_mice > 1, , drop = FALSE]
      if(nrow(sem.data) > 0){p.trend <- p.trend + geom_errorbar(data = sem.data, aes(x = age_months, ymin = mean_log2_CPM_plus1 - sem_log2_CPM_plus1, ymax = mean_log2_CPM_plus1 + sem_log2_CPM_plus1), inherit.aes = FALSE, width = 0.4, colour = "black")}
      p.trend <- p.trend + scale_x_continuous(breaks = age.months, labels = age.order, limits = c(0, 31)) + labs(x = "Age (months)") + common + coord_cartesian(ylim = c(0, ymax))
      stem <- file.map$file_stem[[i]]
      save.png(p.box, file.path(output.dir, "boxplot", paste0(stem, ".mouse_boxplot.png")))
      save.png(p.trend, file.path(output.dir, "age_trend", paste0(stem, ".mouse_age_trend.png")))
      box.plots[[i]] <- p.box
      trend.plots[[i]] <- p.trend
    }

    # 모든 cell type의 그림을 9개씩 모아 저장한다.
    pages <- split(seq_along(celltypes), ceiling(seq_along(celltypes) / 9))
    for(page in seq_along(pages)) {
      idx <- pages[[page]]
      ncols <- min(3L, length(idx))
      nrows <- ceiling(length(idx) / ncols)
      save.png(patchwork::wrap_plots(box.plots[idx], ncol = ncols), file.path(output.dir, "overview", sprintf("%s.all_celltypes.mouse_boxplot.page%02d.png", gene, page)), width = 1200 * ncols, height = 1000 * nrows)
      save.png(patchwork::wrap_plots(trend.plots[idx], ncol = ncols), file.path(output.dir, "overview", sprintf("%s.all_celltypes.mouse_age_trend.page%02d.png", gene, page)), width = 1200 * ncols, height = 1000 * nrows)
    }
    writeLines(c(paste0("Input RDS: ", rds.files[[tissue.name]]), paste0("Gene: ", gene), paste0("Minimum cells per mouse-celltype: ", min.cells.per.mouse.celltype), "Aggregation: sum raw RNA counts over all cells in each mouse-celltype, including gene-zero cells", "Normalization: edgeR TMM separately within each tissue-celltype, using all nonzero genes from eligible mouse groups", "Single eligible mouse: normalization factor=1; descriptive CPM only", "Plot and test input: same mouse log2(normalized CPM + 1); one observation per independent mouse", "Age mean: arithmetic mean of mouse log2(CPM+1); SEM=SD/sqrt(n_mice); SEM=NA for n<2", "Missing groups: no fabricated expression values; excluded groups have normalized values NA", "Test: stats::kruskal.test; tie-corrected H statistic and asymptotic chi-squared p-value", "Comparison: all observed age groups jointly, df = n_age_groups - 1; no sex or batch adjustment; not an age trend test", "Gene expression filter: no CPM cutoff or minimum total gene count; all-zero genes are omitted only from normalization", "No edgeR QL fitting, dispersion estimation, or edgeR DEG test", "Test requirements: >=2 observed age groups, at least one age with >=2 mice, and nonidentical finite normalized target values", "Single-mouse ages are retained if another age has replication; all ages with one mouse yield NA", "Mff all zero, all normalized values identical, insufficient ages/replication, or test errors: p=NA; reasons in test CSV", "Small age groups may make the chi-squared approximation inaccurate; recorded in the warnings column", "FDR: BH over tested cell-type Mff p-values separately within each tissue; not genome-wide FDR", paste0("Boxplot annotation: raw p only in upper right; size=", plot.pvalue.size, " mm; no method name, age-group count, FDR, or NA reason"), "Age-trend plot: existing mouse values, mean and SEM; no p-value annotation", "Results: 08.Mff.KruskalWallis_omnibus.csv includes H, df, per-age mouse counts, raw p, BH FDR, status, reason and warnings", "An older 08.Mff.edgeR_omnibus.csv, if present, belongs to an earlier analysis and is not updated by this script", "TMM is a relative-expression normalization; it does not measure absolute transcripts per cell", "Raw-count aggregation checks: PASSED", "", capture.output(sessionInfo())), file.path(output.dir, "07.analysis_settings_and_sessionInfo.txt"))
    cat("Saved:", output.dir, "\n")
    rm(pb.counts, pb.meta, cell.meta, box.plots, trend.plots, test.rows, test.summary)
    invisible(gc())
  }
  cat("\nDone: mouse-level pseudobulk with Kruskal-Wallis p-values for both tissues.\n")
})
