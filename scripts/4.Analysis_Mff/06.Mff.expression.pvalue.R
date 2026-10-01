

local({
  #-------------------------------------------------------------------
  # 1. 설정

  gene <- "Mff"
  output.base <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/5.Mff_expression"
  tissues <- c("Heart", "Limb_Muscle")
  input.folder <- "mouse_pseudobulk_TMM_min1"  # 기존 03번 CSV가 저장된 폴더명
  # 기존 분석에서 min20 등을 사용했다면 해당 결과 폴더명으로 변경합니다.
  # 이 값은 읽을 폴더만 선택합니다. 세포 수 기준은 03번 CSV의 eligible을 따릅니다.

  age.order <- c("1m", "3m", "18m", "21m", "24m", "30m")
  min.mice.per.age.test <- 2L
  comparison.mode <- "all"       # "all": 모든 age 쌍, "reference": 기준 age와만 비교
  reference.age <- "1m"          # comparison.mode = "reference"일 때 사용

  comparison.mode <- match.arg(comparison.mode, c("all", "reference"))
  if(length(reference.age) != 1L || !reference.age %in% age.order){stop("reference.age must be in age.order.")}
  if(length(min.mice.per.age.test) != 1L || !is.finite(min.mice.per.age.test) || min.mice.per.age.test < 2 || min.mice.per.age.test != floor(min.mice.per.age.test)){stop("min.mice.per.age.test must be an integer >= 2.")}
  input.files <- setNames(file.path(output.base, tissues, input.folder, paste0("03.", gene, ".pseudobulk_mouse.csv")), tissues)
  if(any(!file.exists(input.files))){stop("Missing input CSV. Run the original boxplot script first or check input.folder:\n", paste(input.files[!file.exists(input.files)], collapse = "\n"))}

  #-------------------------------------------------------------------
  # 2. 입력 확인: 저장된 mouse/age/eligible/발현값을 유지합니다.

  read.pseudobulk <- function(filename, tissue.name) {
    pb.meta <- utils::read.csv(filename, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA", fileEncoding = "UTF-8-BOM")
    required <- c("tissue", "celltype", "mouse.id", "age", "gene", "eligible", "log2_CPM_plus1")
    missing.cols <- setdiff(required, names(pb.meta))
    if(length(missing.cols) > 0){stop("Missing columns in ", filename, ": ", paste(missing.cols, collapse = ", "))}
    if(nrow(pb.meta) == 0){stop("Input CSV has no mouse rows: ", filename)}
    for(nm in c("tissue", "celltype", "mouse.id", "age", "gene", "eligible")) {
      if(anyNA(pb.meta[[nm]]) || any(trimws(pb.meta[[nm]]) == "")){stop("Missing/empty ", nm, " in ", filename)}
    }
    if(any(pb.meta$tissue != tissue.name)){stop("Input tissue does not match ", tissue.name, ": ", filename)}
    if(any(pb.meta$gene != gene)){stop("Input gene does not match ", gene, ": ", filename)}
    if(any(!pb.meta$age %in% age.order)){stop("Unexpected age in ", filename, ": ", paste(setdiff(unique(pb.meta$age), age.order), collapse = ", "))}
    if(anyDuplicated(pb.meta[, c("celltype", "mouse.id"), drop = FALSE])){stop("Duplicate mouse x celltype rows in ", filename)}
    mouse.age <- unique(pb.meta[, c("mouse.id", "age"), drop = FALSE])
    if(anyDuplicated(mouse.age$mouse.id)){stop("A mouse.id maps to multiple ages in ", filename)}
    eligible.text <- toupper(trimws(pb.meta$eligible))
    if(any(!eligible.text %in% c("TRUE", "FALSE"))){stop("eligible must contain TRUE/FALSE in ", filename)}
    pb.meta$eligible <- eligible.text == "TRUE"
    pb.meta$log2_CPM_plus1 <- suppressWarnings(as.numeric(pb.meta$log2_CPM_plus1))
    included.values <- pb.meta$log2_CPM_plus1[pb.meta$eligible]
    if(any(!is.finite(included.values)) || any(included.values < 0)){stop("Eligible mice have missing/invalid log2_CPM_plus1 values in ", filename)}
    pb.meta$age <- factor(pb.meta$age, levels = age.order)
    pb.meta
  }

  #-------------------------------------------------------------------
  # 3. 이전 코드와 동일한 검정 함수

  # Welch t-test: 같은 tissue/cell type 안의 독립적인 mouse 값을 비교합니다.
  # 표의 mean difference 및 95% CI는 age2 - age1 방향입니다(CI는 보정 전).
  make.pvalue.table <- function(pb.meta, celltypes, tissue.name) {
    pairs <- utils::combn(age.order, 2, simplify = FALSE)
    if(comparison.mode == "reference"){pairs <- Filter(function(z) reference.age %in% z, pairs)}
    result <- list()
    k <- 0L
    for(ct in celltypes) {
      d <- pb.meta[pb.meta$celltype == ct & pb.meta$eligible, , drop = FALSE]
      if(anyDuplicated(d$mouse.id)){stop("Repeated mouse in statistical input: ", tissue.name, " / ", ct)}
      if(any(!is.finite(d$log2_CPM_plus1))){stop("Non-finite expression in statistical input: ", tissue.name, " / ", ct)}
      for(pair in pairs) {
        v1 <- d$log2_CPM_plus1[as.character(d$age) == pair[[1]]]
        v2 <- d$log2_CPM_plus1[as.character(d$age) == pair[[2]]]
        mean1 <- if(length(v1) > 0) mean(v1) else NA_real_
        mean2 <- if(length(v2) > 0) mean(v2) else NA_real_
        r <- data.frame(tissue = tissue.name, celltype = ct, gene = gene, age1 = pair[[1]], age2 = pair[[2]], n_mice_age1 = length(v1), n_mice_age2 = length(v2), mean_age1 = mean1, mean_age2 = mean2, mean_difference_age2_minus_age1 = mean2 - mean1, ci95_low_unadjusted = NA_real_, ci95_high_unadjusted = NA_real_, t_statistic_age1_minus_age2 = NA_real_, df = NA_real_, p_value = NA_real_, p_adj_BH = NA_real_, method = "Welch two-sample t-test; two-sided; unpaired", value_tested = "log2(normalized_CPM + 1)", status = "not_tested", reason = "", stringsAsFactors = FALSE)
        if(length(v1) < min.mice.per.age.test || length(v2) < min.mice.per.age.test) {
          r$status <- "insufficient_mice"
          r$reason <- paste0("Each age needs at least ", min.mice.per.age.test, " independent mice")
        } else if(stats::var(v1) == 0 && stats::var(v2) == 0) {
          r$status <- "zero_standard_error"
          r$reason <- "Both groups have zero variance; Welch t-test is undefined"
        } else {
          fit <- tryCatch(stats::t.test(v1, v2, alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95), error = function(e) e)
          if(inherits(fit, "error")) {
            r$status <- "test_failed"
            r$reason <- conditionMessage(fit)
          } else if(!is.finite(fit$p.value) || fit$p.value < 0 || fit$p.value > 1) {
            r$status <- "test_failed"
            r$reason <- "Test returned an invalid p-value"
          } else {
            r$status <- "tested"
            r$p_value <- fit$p.value
            r$t_statistic_age1_minus_age2 <- unname(fit$statistic)
            r$df <- unname(fit$parameter)
            r$ci95_low_unadjusted <- -fit$conf.int[[2]]
            r$ci95_high_unadjusted <- -fit$conf.int[[1]]
          }
        }
        k <- k + 1L
        result[[k]] <- r
      }
    }
    result <- do.call(rbind, result)
    ok <- result$status == "tested" & is.finite(result$p_value)
    result$p_adj_BH[ok] <- stats::p.adjust(result$p_value[ok], method = "BH")
    result$p_adjust_scope <- "All successful selected celltype-age comparisons within this tissue"
    result$n_tests_in_BH_family <- sum(ok)
    result
  }

  #-------------------------------------------------------------------
  # 4. Tissue별 계산/저장 (한 행 = cell type 하나의 age 쌍 하나)

  for(tissue.name in tissues) {
    pb.meta <- read.pseudobulk(input.files[[tissue.name]], tissue.name)
    celltypes <- sort(unique(pb.meta$celltype))
    pvalue.results <- make.pvalue.table(pb.meta, celltypes, tissue.name)
    output.file <- file.path(dirname(input.files[[tissue.name]]), paste0("08.", gene, ".age_pairwise_pvalues.csv"))
    utils::write.csv(pvalue.results, output.file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    cat("\nTissue:", tissue.name, "\n")
    cat("Successful tests:", sum(pvalue.results$status == "tested"), "/", nrow(pvalue.results), "\n")
    print(table(pvalue.results$status))
    cat("Saved:", output.file, "\n")
  }
  cat("\nDone: p-value CSVs saved for each tissue.\n")
})
