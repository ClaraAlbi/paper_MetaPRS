library(tidyverse)
library(bigsnpr)
library(future.batchtools)
library(glue)

rsid <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$map$rsid

sims <- readRDS("data/sims_ind.rds") %>% mutate(p = p/1000, size = round(size/1000, 0)) 

sims_ldpred <- sims  %>%
  mutate(gwas_int = glue("results/gwas/{size}_{p}_{rg}_{rep}_ind_{part_ind}.txt"),
         gwas_ext = glue("results/gwas/{size}_{p}_{rg}_{rep}_ss_{part_ss}.txt"),
         gwas_meta = glue("results/metal/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}1.txt"),
         Ncase_int = map_dbl(phen_ind, ~ sum(.x == 1)), 
         Ncontrol_int = map_dbl(phen_ind, ~ sum(.x == 0)), 
         Ncase_ext = map_dbl(phen_ss, ~ sum(.x == 1)),
         Ncontrol_ext = map_dbl(phen_ss, ~ sum(.x == 0)),
         Ncase_meta = Ncase_int + Ncase_ext,
         Ncontrol_meta = Ncontrol_int + Ncontrol_ext) %>%
  select(1:3, 9:11, 22:30) %>%
  pivot_longer(7:15, names_to = c(".value", "type"),
               names_pattern = "(.+)_(.+)") %>%
  mutate(Neff = round(4/(1/Ncase + 1/Ncontrol), digits = 0),
         pval = case_when(type %in% c("int", "ext") ~ "P",
                          type == "meta" ~ "P-value"),
         chr = case_when(type %in% c("int", "ext") ~ "CHR",
                         type == "meta" ~ "Chromosome"),
         pos = case_when(type %in% c("int", "ext") ~ "POS",
                         type == "meta" ~ "Position"),
         SNP = case_when(type %in% c("int", "ext") ~ "SNP",
                         type == "meta" ~ "MarkerName"),
         beta = case_when(type %in% c("int", "ext") ~ "BETA --eff_type LINREG",
                          type == "meta" ~ "Effect --eff_type LINREG"),
         A1 = case_when(type %in% c("int", "ext") ~ "A1",
                        type == "meta" ~ "Allele1"),
         A2 = case_when(type %in% c("int", "ext") ~ "A2",
                        type == "meta" ~ "Allele2")) %>%
  filter(type %in% c("int")) %>%
  filter(size == 50) %>%
  filter(!file.exists(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_{type}.rds"))) %>%
  mutate(pt_file = glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_{type}.pt")) %>% filter(file.exists(pt_file)) %>%
  mutate(file = map(pt_file, ~read_lines(.x)),
       rds = glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_{type}.rds"),
       tmp = map_chr(file, ~str_match(.x[str_detect(.x, "Coordinated data filename")], "name\\s*(.*?)\\s*.hdf5")[2])) %>%
  select(-file, -rds, -pt_file) 

a <- sims %>% filter(!file.exists(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_{type}.rds")))

  
#LDpred

NCORES <- 4
plan(batchtools_slurm(resources = list(
  t = "0-24:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(sims_ldpred[, c(-9, -10)], function(p, rg, rep, part_ind,  size, part_ss, type, gwas, Neff, pval, chr, pos, SNP, beta, A1, A2) {
  
  tmp <- tempfile(tmpdir = "tmp")
  file_hdf5 <- paste0(tmp, ".hdf5")
  phen <- glue("{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_{type}")
  # system(glue::glue("ldpred coord --help"))
  system(glue::glue(
    "ldpred coord",
    " --gf /faststorage/jail/project/NCRR-PRS/clara/paper_prs/UKBB/data/UKBB_imp_HM3_5k",
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
  
  
  ext_ldpred <- c(sprintf("_LDpred_p%.4e.txt", c(1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001)),
                  "_LDpred-inf.txt")
  
  files_ldpred <- paste0(tmp, ext_ldpred)

  beta_ldpred <- sapply(files_ldpred, function(file) {
    res_ldpred <- bigreadr::fread2(file, select = c(3, 7))
    beta_ldpred <- rep(0, length(rsid))
    beta_ldpred[match(res_ldpred$sid, rsid)] <- res_ldpred[[2]]
    beta_ldpred
  })
  unlink(paste0(tmp, "*"))
  saveRDS(cbind(beta_pt, beta_ldpred), paste0("results/ldpred/", phen, ".rds"))
})


# SCT

G_sub <- snp_attach("/faststorage/jail/project/NCRR-PRS/clara/paper_prs/UKBB/data/UKBB_imp_HM3_5k.rds")
ukb <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
map <- ukb$map %>% rename(chr = chromosome, pos = physical.pos, a0 = allele2, a1 = allele1) %>% select(1, 3:6)

sims <- readRDS("data/sims_ind.rds") %>% mutate(p = p/1000, size = round(size/1000, 0)) %>%
  mutate(gwas_ext = glue("results/gwas/{size}_{p}_{rg}_{rep}_ss_{part_ss}.rds")) %>% filter(size == 50) %>% 
  filter(!file.exists(glue("results/sct/sct_model_{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}.rds")))

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-04:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(sims[, c(1:3, 9:11, 16, 18, 22)], function(p, rg, rep, part_ind,  size, part_ss, train_ind, phen_ind, gwas_ext) {
  phen <- glue("{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}")
  
  sumstats <- readRDS(gwas_ext)
  lpS <- -predict(sumstats, log10 = TRUE)
  beta <- sumstats$estim
  
  all_keep <- snp_grid_clumping(G_sub$genotypes,
                                G_sub$map$chromosome,
                                G_sub$map$physical.pos,
                                lpS = lpS,
                                ncores = 24,
                                exclude = which(is.na(lpS)))
  saveRDS(all_keep, glue("results/sct/{phen}_ext_clump_50k.rds"))
  
  system.time(multi_PRS <- snp_grid_PRS(G = G,
                                        all_keep = all_keep,
                                        betas = beta,
                                        lpS = lpS,
                                        ind.row = train_ind,
                                        grid.lpS.thr = seq_log(0.1, 0.9999 * max(lpS, na.rm = TRUE), 50),
                                        backingfile = glue("results/sct/multi_PRS_{phen}"), 
                                        ncores = 24))
  
  
  final_mod <-  snp_grid_stacking(multi_PRS,
                                  phen_ind,
                                  ncores = 24)
  print(summary(final_mod$mod))
  saveRDS(final_mod, glue("results/sct/sct_model_{phen}.rds"))
})




# PRS from betas

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

G_ukb <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$genotypes
rsid <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$map$rsid

sims_prs <- sims %>% filter(size == 50) %>% 
  filter(!file.exists(glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ldpred_int.rds")))
  filter(file.exists(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_int.rds"))) %>%
  filter(file.size(glue("results/bolt/{size}_{p}_{rg}_{rep}_{part_ind}.pred")) > 0) %>%
  filter(file.exists(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ext.rds"))) %>%
  filter(file.exists(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_meta.rds"))) %>%
  filter(file.exists(glue("results/sct/sct_model_{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}.rds"))) %>%
  filter(!file.exists(glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_bolt.rds")))

furrr::future_pmap(sims_prs[,c(1:3, 9:11, 14:15)], function(p, rg, rep, part_ind, size, part_ss, test,val) {
  
  phen <- glue("{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}")
  
  ind.test <- c(val, test)
  # Bolt
  betas_bolt <- read_table2(glue("results/bolt/{size}_{p}_{rg}_{rep}_{part_ind}.pred"))$BETA

  pred_bolt <- big_prodVec(G_ukb, betas_bolt, ind.row = ind.test, ncores = 24)
  saveRDS(pred_bolt, glue("results/prs/{phen}_bolt.rds"))

  #LDpred
  betas_ldpred_meta <- readRDS(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_meta.rds"))
  betas_ldpred_ext <- readRDS(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ext.rds"))
  betas_ldpred_int <- readRDS(glue("results/ldpred/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_int.rds"))

  pred_ldpred <- big_prodMat(G_ukb, as.matrix(bind_cols(betas_ldpred_meta, betas_ldpred_ext, betas_ldpred_int)), ind.row = ind.test, ncores = 24)

  saveRDS(pred_ldpred[, 1:ncol(betas_ldpred_meta)], glue("results/prs/{phen}_ldpred_meta.rds"))
  saveRDS(pred_ldpred[, (ncol(betas_ldpred_meta) + 1):(ncol(betas_ldpred_ext))], glue("results/prs/{phen}_ldpred_ext.rds"))
  saveRDS(pred_ldpred[, (ncol(betas_ldpred_ext) + 1):(ncol(pred_ldpred))], glue("results/prs/{phen}_ldpred_int.rds"))

  final_mod <- readRDS(glue("results/sct/sct_model_{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}.rds"))
  new_beta <- final_mod$beta.G
  ind <- which(new_beta != 0)
  pred_sct <- final_mod$intercept + bigstatsr::big_prodVec(G_ukb,
                                                           new_beta[ind],
                                                           ind.col = ind,
                                                           ind.row = ind.test)
  saveRDS(pred_sct, glue("results/prs/{phen}_sct.rds"))
})
