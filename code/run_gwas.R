library(tidyverse)
library(bigsnpr)
library(future.batchtools)

cov <- readRDS("data/covar_ipsych12_15.rds")
pars <- readRDS("data/cv_subsets_ipsych12_15_info.rds")
ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")


# GWAS logreg iPSYCH

mean(has_no_na <- complete.cases(cov[1:22]))
remove <- which(!has_no_na)
cov_m <- covar_from_df(cov[1:22])

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(pars[,c(1, 3:5)], function(p, cv, ind.train, phen.train) {
  ind.keep <- !(ind.train %in% remove)
  gwas <- bigstatsr::big_univLogReg(
    ipsych$genotypes, phen.train[ind.keep], ind.train = ind.train[ind.keep],
    covar.train = cov_m[ind.train[ind.keep], ], ncores = 24)
  
  saveRDS(gwas, paste0("results/gwas/", p, "_", cv, ".rds"))
})

pmap(pars[,c(1,3)], function(p, cv) {
  file <- readRDS(paste0("results/gwas/", p, "_", cv, ".rds"))
  pval <- predict(file, log10 = FALSE)
  file_txt <- bind_cols(ipsych$map[,c(1:3,7:8)], file, tibble(pval = pval)) %>%
    rename(A1 = a1, A2 = a2, BETA = estim, SE = std.err, P = pval) %>% select(-score)
  write_delim(file_txt, paste0("results/gwas/", p, "_", cv, ".txt"))
})



# BOLT-LMM
NCORES <- 32
plan(batchtools_slurm(resources = list(
  t = "00-24:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

pars <- pars %>%
  filter(!file.exists(glue::glue("results/bolt/{p}_{cv}.pred"))) 

furrr::future_pmap(pars[,c(1,3)], function(p, cv) {
  phen <- paste0(p, "_", cv)
  system(glue::glue(
    "/home/clara/REPOS/BOLT-LMM_v2.3.4/bolt",
    " --bfile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/ipsych15/data/ipsych2015_hm3",
    " --lmm",
    " --maxModelSnps 1200000",
    " --phenoFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/data/cv_subsets_ipsych12_15.txt",
    " --phenoCol={phen}",
    " --LDscoresFile=/home/clara/REPOS/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz",
    " --geneticMapFile=/home/clara/REPOS/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz",
    " --covarFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/data/covar_ipsych12_15.txt",
    " --covarCol=sex",
    " --qCovarCol=PC{{1:20}}",
    " --qCovarCol=age",
    " --numThreads=32",
    " --statsFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.stats",
    " --predBetasFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.pred",
    " --verboseStats",
    " &>> /faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.out"
  ))
})




#Meta analysis

NCORES <- 1
plan(batchtools_slurm(resources = list(
  t = "0-01:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(pars[, c(1, 3, 8:11, 13:14)], function(p, cv, gwas_int, gwas_ext, N_int_case, N_int_control, N_ext_case, N_ext_control) {
  outfile <- paste0("results/metal/", p, "_", cv)
  run_metal <- paste0("results/metal/run_metal_", p, "_", cv, ".txt")
  write_file(glue::glue(
    "TRACKPOSITIONS ON
SCHEME STDERR
MARKER SNP
CHROMOSOME CHR
POSITION POS
ALLELE  A1 A2
PVALUE   P
STDERR SE
EFFECT BETA
DEFAULTWEIGHT {4/(1/N_int_case + 1/N_int_control)}
PROCESS {gwas_int}
POSITION BP
EFFECT log(OR)
DEFAULTWEIGHT {4/(1/N_ext_case + 1/N_ext_control)}
PROCESS {gwas_ext}
OUTFILE {outfile} .txt
ANALYZE"), run_metal)
  
  system(glue::glue("~/REPOS/METAL/build/bin/metal {run_metal}"))
})


### UKB

pars <- readRDS("data/cv_subsets_ukb_traits_info.rds")

ukb_map <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$map
gwas <- readRDS("ukbb_traits/results/gwas/all_gwas.rds")

pmap(pars[,c(1,3)], function(p, cv) {
  file <- gwas[[paste0(p, "_train_cv", cv)]]
  pval <- predict(file, log10 = FALSE)
  file_txt <- bind_cols(ukb_map[,c(1,3:6)], file, tibble(pval = pval)) %>%
    rename(CHR = chromosome, SNP = rsid, POS = physical.pos, A1 = allele1, A2 = allele2, BETA = estim, SE = std.err, P = pval) %>% select(-score)
  write_delim(file_txt, paste0("results/gwas/", p, "_", cv, ".txt"))
})

