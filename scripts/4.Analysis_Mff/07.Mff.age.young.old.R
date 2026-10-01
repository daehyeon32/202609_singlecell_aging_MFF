library(Seurat)

# 1. RDS 읽기
rds.file <- "/BiO/Live/dleogus32/202609Aging_MFF/Analysis/3.preprocessing/Aging.MFF.seurat.metadata.filtered.normalization.pca.umap.RDS"

mff.rds <- readRDS(rds.file)

# 2. Young / Old 기준 설정
young.ages <- c("1m", "3m")
old.ages <- c("18m", "21m", "24m", "30m")

stopifnot(inherits(mff.rds, "Seurat"), "age" %in% colnames(mff.rds@meta.data))

age.value <- as.character(mff.rds$age)

# 누락되거나 예상하지 못한 age가 있으면 중단
stopifnot(!anyNA(age.value), all(age.value %in% c(young.ages, old.ages)))

# 3. Young / Old 라벨 추가
mff.rds$age_young_old <- factor(ifelse(age.value %in% young.ages, "Young", "Old"), levels = c("Young", "Old"))

# 4. 결과 확인: 표의 숫자는 cell 수
print(table(Age = factor(age.value, levels = c(young.ages, old.ages)), Group = mff.rds$age_young_old))

print(head(mff.rds@meta.data[, c("age", "age_young_old", "mouse.id")]))

# 5. 같은 경로에 저장
saveRDS(mff.rds, rds.file)

cat("Saved:", rds.file, "\n")
