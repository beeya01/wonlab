library(data.table)

ss1 <- fread("/data4/ebjang/Project/251110_CVD_trans/validation/CAD/bbj_CAD.sumstats")
ss2 <- fread("/data4/ebjang/Project/251110_CVD_trans/validation/CAD/finngen_CAD.sumstats")
ss3 <- fread("/data4/ebjang/Project/251110_CVD_trans/validation/CAD/mvp_CAD.sumstats")

common_snps <- Reduce(intersect, list(
  ss1$SNP,
  ss2$SNP,
  ss3$SNP
))

ss1 <- ss1[SNP %in% common_snps]
ss2 <- ss2[SNP %in% common_snps]
ss3 <- ss3[SNP %in% common_snps]

cat("Common SNPs:", length(common_snps), "\n")

