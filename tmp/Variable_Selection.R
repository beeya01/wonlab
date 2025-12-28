############################################################
## 0. Library
############################################################
library(data.table)
library(openxlsx)
library(glmnet)
library(pROC)

############################################################
## 1. Data load
############################################################
phenotype <- fread("~/CVD_PRS/01.KoGES_phenotype/missforest/Test_CVD_pheno_T2D_noT1D.txt")

pc_result <- fread("/data5/ebjang/Project/251010_CVDprs_QC/test_PC.eigenvec")
colnames(pc_result) <- c("ID","IID", paste0("PC",1:10))

############################################################
## 2. PRS files
############################################################
cad_path <- "/data5/ebjang/Project/251027_CVD_PRS/CAD/CAD_result/test_PRS.sscore"
sbp_path <- "/data5/ebjang/Project/251027_CVD_PRS/SBP/SBP_result/test_PRS.sscore"

prs_dir   <- "/data5/ebjang/Project/251027_CVD_PRS/MultiplePRS/"
all_files <- list.files(prs_dir, pattern="\\.sscore$", full.names=TRUE)

# CAD/SBP가 prs_dir에 같이 잡히면 중복 방지
all_files <- setdiff(normalizePath(all_files), normalizePath(c(cad_path, sbp_path)))

############################################################
## 3. Helper: read .sscore -> (ID, PRS컬럼 1개)로 정리
############################################################
read_sscore_as_prs <- function(f, prs_name = NULL, id_col = "ID") {
  dt <- fread(f)
  
  # ID 컬럼 통일 (#FID / FID / ID 등)
  if ("#FID" %in% names(dt)) setnames(dt, "#FID", id_col)
  if ("FID"  %in% names(dt)) setnames(dt, "FID",  id_col)
  
  # PRS 점수 컬럼 자동 선택 (plink2 .sscore 기준)
  score_col <- if ("SCORE1_SUM" %in% names(dt)) {
    "SCORE1_SUM"
  } else if ("SCORE1_AVG" %in% names(dt)) {
    "SCORE1_AVG"
  } else {
    cand <- grep("^SCORE", names(dt), value = TRUE)
    if (length(cand) == 0) stop("No SCORE column found in: ", f)
    cand[1]
  }
  
  # PRS 변수명 만들기
  if (is.null(prs_name) || prs_name == "") {
    prs_name <- sub("\\.sscore$", "", basename(f))
    prs_name <- gsub("[^A-Za-z0-9_]+", "_", prs_name)
  }
  
  out <- dt[, .(ID = get(id_col), score = get(score_col))]
  setnames(out, "score", prs_name)
  out
}

############################################################
## 4. Build total = phenotype + PC + all PRS
############################################################
if (!("ID" %in% names(phenotype))) setnames(phenotype, names(phenotype)[1], "ID")

total <- merge(phenotype, pc_result, by = "ID", all.x = TRUE)

prs_list <- list(
  read_sscore_as_prs(cad_path, prs_name = "PRS_CAD"),
  read_sscore_as_prs(sbp_path, prs_name = "PRS_SBP")
)

if (length(all_files) > 0) prs_list <- c(prs_list, lapply(all_files, read_sscore_as_prs))

for (p in prs_list) total <- merge(total, p, by = "ID", all.x = TRUE)

############################################################
## 5. Robust recode: CVD_case -> 0/1  (⭐ 핵심 수정)
############################################################
recode_binary01 <- function(v, name="CVD_case") {
  # 원본 보존용
  vv <- v
  
  # factor/label -> character
  if (is.factor(vv)) vv <- as.character(vv)
  
  # logical
  if (is.logical(vv)) {
    out <- as.integer(vv)
    return(out)
  }
  
  # character 처리
  if (is.character(vv)) {
    x <- tolower(trimws(vv))
    x[x %in% c("", "na", "nan", "-", ".", "null")] <- NA
    
    # 흔한 라벨 매핑
    x[x %in% c("case","disease","yes","true")]    <- "1"
    x[x %in% c("control","normal","no","false")]  <- "0"
    
    suppressWarnings(num <- as.integer(x))
    
    # 0/1 또는 1/2 처리
    if (all(na.omit(num) %in% c(0L,1L))) return(as.integer(num))
    if (all(na.omit(num) %in% c(1L,2L))) return(as.integer(num - 1L))
    
    # 여기서 못 정하면 값 출력하고 중단
    print(sort(unique(na.omit(x))))
    stop(name, " cannot be recoded to 0/1. Check unique values above.")
  }
  
  # numeric/integer 처리
  if (is.numeric(vv) || is.integer(vv)) {
    num <- as.integer(vv)
    if (all(na.omit(num) %in% c(0L,1L))) return(num)
    if (all(na.omit(num) %in% c(1L,2L))) return(num - 1L)
    
    print(sort(unique(na.omit(num))))
    stop(name, " numeric values are not 0/1 or 1/2. Check unique values above.")
  }
  
  stop(name, ": unsupported type")
}

dt <- copy(total)
dt[, CVD_case := recode_binary01(CVD_case, "CVD_case")]

# sanity check
cat("CVD_case distribution:\n")
print(table(dt$CVD_case, useNA="ifany"))
stopifnot(all(na.omit(dt$CVD_case) %in% c(0L, 1L)))

############################################################
## 6. Model setting (fixed + candidate)
############################################################
# CAD/SBP PRS 컬럼명 안전장치
cad_var <- intersect(c("PRS_CAD","CAD_PRS"), names(dt))[1]
sbp_var <- intersect(c("PRS_SBP","SBP_PRS"), names(dt))[1]
stopifnot(!is.na(cad_var), !is.na(sbp_var))

fixed_vars <- c("age","sex", paste0("PC",1:10), cad_var, sbp_var)

cand_vars <- c(
  "T2D","smoke","bmi","drink","glucose","HDL","cholesterol","triglyceride",
  "BMI_PRS7","ens_PRS10","Glucose_PRS35","T2D_PRS","TC_PRS5","TG_PRS5"
)

fixed_vars <- fixed_vars[fixed_vars %in% names(dt)]
cand_vars  <- cand_vars[cand_vars %in% names(dt)]

# sex factor
if ("sex" %in% names(dt)) dt[, sex := as.factor(sex)]

# (선택) 위험요인 중 문자형인데 숫자여야 하는 것들 강제 numeric 변환
force_num <- intersect(c("age","bmi","glucose","HDL","cholesterol","triglyceride"), names(dt))
for (v in force_num) {
  if (is.character(dt[[v]]) || is.factor(dt[[v]])) {
    x <- as.character(dt[[v]])
    x <- trimws(x)
    x[x %in% c("", "NA", "NaN", "-", ".")] <- NA
    suppressWarnings(dt[[v]] <- as.numeric(x))
  }
}

use_cols <- unique(c("CVD_case", fixed_vars, cand_vars))
dt2 <- dt[complete.cases(dt[, ..use_cols])]

cat("N used (complete cases) =", nrow(dt2), "\n")
stopifnot(all(dt2$CVD_case %in% c(0L,1L)))

############################################################
## 7. LASSO selection (fixed penalty=0, candidates penalty=1)
############################################################
f_full <- as.formula(paste("CVD_case ~", paste(c(fixed_vars, cand_vars), collapse=" + ")))
x <- model.matrix(f_full, data=dt2)[, -1, drop=FALSE]
y <- dt2$CVD_case

pen <- rep(1, ncol(x)); names(pen) <- colnames(x)
for (v in fixed_vars) pen[grep(paste0("^", v), names(pen))] <- 0

set.seed(1)
cvfit <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha  = 1,
  penalty.factor = pen,
  nfolds = 10,
  type.measure = "auc"
)

# 보수적 선택: 1se (원하면 min으로 바꾸기)
lambda_use <- cvfit$lambda.1se
# lambda_use <- cvfit$lambda.min

b <- coef(cvfit, s=lambda_use)
nz <- rownames(b)[as.numeric(b) != 0]
nz <- setdiff(nz, "(Intercept)")

selected_cand <- cand_vars[sapply(cand_vars, function(v) any(grepl(paste0("^", v), nz)))]
cat("Selected candidate vars:\n", paste(selected_cand, collapse=", "), "\n\n")

############################################################
## 8. Final GLM refit + Base GLM
############################################################
final_vars <- c(fixed_vars, selected_cand)

f_final <- as.formula(paste("CVD_case ~", paste(final_vars, collapse=" + ")))
m_final <- glm(f_final, data=dt2, family=binomial())

f_base  <- as.formula(paste("CVD_case ~", paste(fixed_vars, collapse=" + ")))
m_base  <- glm(f_base, data=dt2, family=binomial())

############################################################
## 9. AUC + DeLong test
############################################################
p_final <- predict(m_final, type="response")
p_base  <- predict(m_base,  type="response")

roc_final <- roc(dt2$CVD_case, p_final, quiet=TRUE)
roc_base  <- roc(dt2$CVD_case, p_base,  quiet=TRUE)

auc_final <- as.numeric(auc(roc_final))
auc_base  <- as.numeric(auc(roc_base))

cat(sprintf("AUC (Final selected model) = %.4f\n", auc_final))
cat(sprintf("AUC (Base fixed model)     = %.4f\n", auc_base))

delong <- roc.test(roc_final, roc_base, method="delong", paired=TRUE)
print(delong)

cat("\nAUC CI (Final):\n"); print(ci.auc(roc_final))
cat("\nAUC CI (Base):\n");  print(ci.auc(roc_base))

############################################################
## 10. Save outputs (optional)
############################################################
out_dir <- '/data5/ebjang/Project/251027_CVD_PRS/Result'
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# 선택 변수 저장
fwrite(
  data.table(selected_candidate_vars = selected_cand),
  file.path(out_dir, "selected_candidate_vars.txt"),
  sep="\t"
)

# AUC 결과 저장
res <- data.table(
  model = c("final","base"),
  auc   = c(auc_final, auc_base)
)
fwrite(res, file.path(out_dir, "auc_compare_final_vs_base.txt"), sep="\t")

# 최종모형 계수 저장
coef_final <- data.table(term = names(coef(m_final)), beta = as.numeric(coef(m_final)))
fwrite(coef_final, file.path(out_dir, "final_model_coefficients.txt"), sep="\t")

cat("\nSaved to:\n",
    file.path(out_dir, "selected_candidate_vars.txt"), "\n",
    file.path(out_dir, "auc_compare_final_vs_base.txt"), "\n",
    file.path(out_dir, "final_model_coefficients.txt"), "\n")
