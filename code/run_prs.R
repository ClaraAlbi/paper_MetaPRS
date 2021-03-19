library(tidyverse)
library(bigsnpr)
library(future.batchtools)
library(glue)

rsid <- snp_attach("data/ipsych_hm3_sub5k_cohort.rds")$map$SNP

pars <- readRDS("data/cv_subsets_ipsych12_15_info.rds") %>%
  rename(Ncase_ext = N_ext_case,
         Ncontrol_ext = N_ext_control,
         Ncase_int = N_int_case,
         Ncontrol_int = N_int_control) %>%
  mutate(Ncase_meta = Ncase_int + Ncase_ext,
         Ncontrol_meta = Ncontrol_int + Ncontrol_ext) %>%
  select(-cv_sub, -ind.train, -phen.train, -ind.test, -phen.test) %>%
  pivot_longer(3:12, names_to = c(".value", "type"),
               names_pattern = "(.+)_(.+)") %>%
  filter(!is.na(type)) %>%
  mutate(Neff = round(4/(1/Ncase + 1/Ncontrol), digits = 0),
         pval = case_when(type %in% c("int", "ext") ~ "P",
                          type == "meta" ~ "P-value"),
         chr = case_when(type %in% c("int", "ext") ~ "CHR",
                          type == "meta" ~ "Chromosome"),
         pos = case_when(type == "int" ~ "POS",
                         type == "ext" ~ "BP",
                         type == "meta" ~ "Position"),
         SNP = case_when(type %in% c("int", "ext") ~ "SNP",
                         type == "meta" ~ "MarkerName"),
         beta = case_when(type == "int" ~ "BETA --eff_type LOGOR",
                          type == "ext" ~ "OR --eff_type OR",
                         type == "meta" ~ "Effect --eff_type LOGOR"),
         A1 = case_when(type %in% c("int", "ext") ~ "A1",
                         type == "meta" ~ "Allele1"),
         A2 = case_when(type %in% c("int", "ext") ~ "A2",
                        type == "meta" ~ "Allele2")) %>%
  filter(case_when(type == "ext" ~ cv == 1,
                   type %in% c("int", "meta") ~ cv %in% 1:5))



#LDpred

NCORES <- 4
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(pars[, c(1:3, 6:14)], function(p, cv, type, gwas, pval, chr, pos, SNP, A1, A2, beta, Neff) {
  
  tmp <- tempfile(tmpdir = "tmp")
  file_hdf5 <- paste0(tmp, ".hdf5")
  phen <- paste0(p, "_", cv, "_", type)
  # system(glue::glue("ldpred coord --help"))
  system(glue::glue(
    "ldpred coord",
    " --gf /faststorage/jail/project/NCRR-PRS/clara/paper_prs/data/ipsych_hm3_sub5k_cohort",
    " --ssf {gwas}",
    " --only-hm3",
    " --rs {SNP} --A1 {A1} --A2 {A2} --pos {pos} --chr {chr}",
    " --pval {pval} --eff {beta}",
    " --N {Neff}",
    " --out {file_hdf5}",
    " &>> /faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/ldpred/{phen}.coord"
  ))

  # system(glue::glue("ldpred p+t --help"))
  system(glue::glue(
    "ldpred p+t",
    " --cf {file_hdf5}",
    " --ldr {round(1e6 / 3000)}",
    " --out {tmp}",
    " &>> /faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/ldpred/{phen}.pt"
  ))
  
  
  ext_pt <- c(sprintf("_P+T_r0.20_p%.4e.txt", c(10^(0:-8), 3 * 10^(-1:-5))))
  files_pt <- paste0(tmp, ext_pt)
  beta_pt <- sapply(files_pt, function(file) {
    res_pt <- bigreadr::fread2(file, select = c(3, 8))
    beta_pt <- rep(0, length(rsid))
    beta_pt[match(res_pt$sid, rsid)] <- res_pt[[2]]
    beta_pt})

  
  # system(glue::glue("ldpred gibbs --help"))
  system(glue::glue(
    "ldpred gibbs",
    " --cf {file_hdf5}",
    " --ldr {round(1e6 / 3000)}",
    " --N {Neff}",
    " --ldf {tmp}",
    " --out {tmp}",
    " &>> /faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/ldpred/{phen}.gibbs"
  ))
  

  ext_ldpred <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)))
  #files_ldpred <- list.files("tmp", pattern = paste0(str_sub(tmp, start = 5), "_LDpred"), full.names = TRUE)
  
  files_ldpred <- paste0(tmp, ext_ldpred)
  beta_ldpred <- sapply(files_ldpred, function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- rep(0, length(rsid))
    beta_ldpred[match(res_ldpred$sid, rsid)] <- res_ldpred[[2]]
    beta_ldpred
  })
  unlink(paste0(tmp, "*"))
  saveRDS(cbind(beta_pt, beta_ldpred), paste0("results/ldpred/", p, "_", cv, "_", type, ".rds"))
})


## SCT betas

G_sub <- snp_attach("data/ipsych_hm3_sub5k_cohort.rds")

ext_gwas <- pars %>% filter(type == "ext") %>%
  mutate(phen = str_sub(p, end = -5)) %>%
  select(phen, gwas) %>%
  group_by(gwas) %>%
  filter(duplicated(phen))


# Make clumping subsets
pmap(ext_gwas[1:5,], function(phen, gwas) {
  sumstats <- read_table2(gwas) %>% 
    rename(chr = CHR, rsid = SNP, pos = BP, a0 = A2, a1 = A1) %>% 
    mutate(beta = log(OR),
           z = beta/SE,
           log.pval = log(2) + pnorm(abs(z), lower.tail=FALSE, log.p = TRUE),
           log10.pval = log.pval/log(10)) %>%
    filter(rsid %in% map$rsid)
  
  lpS <- rep(0, nrow(map))
  lpS[match(sumstats$rsid, map$rsid)] <- -sumstats$log10.pval
  all_keep <- snp_grid_clumping(G_sub$genotypes,
                                G_sub$map$CHR,
                                G_sub$map$POS,
                                lpS = lpS,
                                ncores = 24,
                                exclude = which(is.na(lpS)))
  saveRDS(all_keep, glue("results/sct/{phen}_ext_clump_50k.rds"))
})

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-04:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

map <- snp_attach("data/ipsych_hm3_sub5k_cohort.rds")$map %>% rename(chr = CHR, rsid = SNP, pos = POS, a0 = a2)
G <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")$genotypes

# Make grid PRS and run lasso
furrr::future_pmap(pars[c(1:50, 52:60), c(1, 3:5, 11)], function(p, cv, ind.train, phen.train, gwas_ext) {
  phen <- str_sub(p, end = -5)
  all_keep <- readRDS(glue("results/sct/{str_sub(p, end = -5)}_ext_clump_50k.rds"))

  sumstats <- read_table2(gwas_ext) %>% 
    rename(chr = CHR, rsid = SNP, pos = BP, a0 = A2, a1 = A1) %>% 
    mutate(beta = log(OR),
           z = beta/SE,
           log.pval = log(2) + pnorm(abs(z), lower.tail=FALSE, log.p = TRUE),
           log10.pval = log.pval/log(10)) %>%
    filter(rsid %in% map$rsid)
  
  matched <- snp_match(sumstats, map, join_by_pos = FALSE)
  
  lpS <- rep(0, nrow(map)); beta <- rep(0, nrow(map))
  lpS[matched$`_NUM_ID_`] <- -matched$log10.pval
  beta[matched$`_NUM_ID_`] <- matched$beta
  
  system.time(multi_PRS <- snp_grid_PRS(G = G,
                            all_keep = all_keep,
                            betas = beta,
                            lpS = lpS,
                            ind.row = ind.train,
                            grid.lpS.thr = seq_log(0.1, 0.9999 * max(lpS, na.rm = TRUE), 50),
                            backingfile = glue("results/sct/multi_PRS_{p}_{cv}"), 
                            ncores = 24))
  

  final_mod <-  snp_grid_stacking(multi_PRS,
                                  phen.train,
                                  ncores = 24)
  print(summary(final_mod$mod))
  saveRDS(final_mod, glue("results/sct/sct_model_{p}_{cv}.rds"))
})





# Calculate PRSs from betas

#iPSYCH

G_ldpred <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")$genotypes
G_bolt <- snp_attach("/faststorage/jail/project/NCRR-PRS/clara/paper_prs/ipsych15/data/ipsych2015_hm3.rds")$genotypes

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-06:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(pars[,c(1, 3,6)], function(p, cv, ind.test) {
  
  # Bolt
  betas_bolt <- read_table2(glue("results/bolt/{p}_{cv}.pred"))$BETA
  
  pred_bolt <- big_prodVec(G_bolt, betas_bolt, ind.row = ind.test, ncores = 24)
  saveRDS(pred_bolt, glue("results/prs/{p}_{cv}_bolt.rds"))
  
  #LDpred
  betas_ldpred_ind <- readRDS(glue("results/ldpred/{p}_{cv}_int.rds"))
  betas_ldpred_meta <- readRDS(glue("results/ldpred/{p}_{cv}_meta.rds"))
  betas_ldpred_ext <- readRDS(glue("results/ldpred/{p}_1_ext.rds"))
  
  pred_ldpred <- big_prodMat(G_ldpred, as.matrix(bind_cols(betas_ldpred_ind, betas_ldpred_meta, betas_ldpred_ext)), ind.row = ind.test, ncores = 24)
  
  saveRDS(pred_ldpred[, 1:ncol(betas_ldpred_ind)], glue("results/prs/{p}_{cv}_ldpred_int.rds"))
  saveRDS(pred_ldpred[, (ncol(betas_ldpred_ind)+1):(ncol(betas_ldpred_ind) + ncol(betas_ldpred_meta))], glue("results/prs/{p}_{cv}_ldpred_meta.rds"))
  saveRDS(pred_ldpred[, (ncol(betas_ldpred_ind) + ncol(betas_ldpred_meta) + 1):(ncol(pred_ldpred))], glue("results/prs/{p}_{cv}_ldpred_ext.rds"))
  
  final_mod <- readRDS(glue("results/sct/sct_model_{p}_{cv}.rds"))
  new_beta <- final_mod$beta.G
  ind <- which(new_beta != 0)
  pred_sct <- final_mod$intercept + bigstatsr::big_prodVec(G_ldpred,
                                                       new_beta[ind],
                                                       ind.col = ind,
                                                       ind.row = ind.test)
  saveRDS(pred_sct, glue("results/prs/{p}_{cv}_sct.rds"))
})

#UKB

pars_ukb <- readRDS("data/cv_subsets_ukb_traits_info.rds") %>% filter(p == "mdd_all")

G_ukb <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$genotypes
rsid <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$map$rsid

furrr::future_pmap(pars_ukb[,c(1:2,5)], function(p, cv, ind.test) {
  
  # Bolt
  betas_bolt <- read_table2(glue("ukbb_traits/results/bolt/ipsych_{p}_train_cv{cv}.pred"))$BETA
  
  pred_bolt <- big_prodVec(G_ukb, betas_bolt, ind.row = ind.test, ncores = 24)
  saveRDS(pred_bolt, glue("results/prs/{p}_{cv}_bolt.rds"))
  
  #LDpred
  ext_ldpred <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)))
  files_ldpred_int <- glue("ukbb_traits/results/ldpred/coord_{p}_cv{cv}_int_WEIGHTS{ext_ldpred}")
  files_ldpred_meta <- glue("ukbb_traits/results/ldpred/coord_{p}_cv{cv}_meta_WEIGHTS{ext_ldpred}")
  files_ldpred_ext <- glue("ukbb_traits/results/ldpred/coord_{p}_ext_WEIGHTS{ext_ldpred}")
  files_ldpred_ext <- list.files("ukbb_traits/results/ldpred", glue("coord_mdd_ext_WEIGHTS"), full.names = T)
  
  betas_ldpred <- sapply(c(files_ldpred_int, files_ldpred_meta, files_ldpred_ext), function(file) {
    res_pt <- bigreadr::fread2(file, select = c(3, 7))
    betas_ldpred <- rep(0, length(rsid))
    betas_ldpred[match(res_pt$sid, rsid)] <- res_pt[[2]]
    betas_ldpred})
  
  pred_ldpred <- big_prodMat(G_ukb, as.matrix(betas_ldpred), ind.row = ind.test, ncores = 24)
  
  saveRDS(pred_ldpred[, 1:length(files_ldpred_int)], glue("results/prs/{p}_{cv}_ldpred_int.rds"))
  saveRDS(pred_ldpred[, (length(files_ldpred_int)+1):(length(files_ldpred_int) + length(files_ldpred_meta))], glue("results/prs/{p}_{cv}_ldpred_meta.rds"))
  saveRDS(pred_ldpred[, (length(files_ldpred_int) + length(files_ldpred_meta) + 1):(ncol(pred_ldpred))], glue("results/prs/{p}_{cv}_ldpred_ext.rds"))
  
  final_mod <- readRDS(glue("ukbb_traits/tmp-data/stacked_model_train_{p}_excl_cv{cv}.rds"))
  new_beta <- final_mod$beta.G
  ind <- which(new_beta != 0)
  pred_sct <- final_mod$intercept + bigstatsr::big_prodVec(G_ukb,
                                                           new_beta[ind],
                                                           ind.col = ind,
                                                           ind.row = ind.test)
  saveRDS(pred_sct, glue("results/prs/{p}_{cv}_sct.rds"))
})
