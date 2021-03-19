library(tidyverse)
library(bigsnpr)

# Define phenotypes and cv subsets

#iPSYCH

ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")
PC <- readRDS("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/PC.rds")
log_dist <- log(bigutilsr::dist_ogk(PC))
is_homogeneous <- log_dist < 4.8

rel <- readRDS("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/rel.rds")
is_rel <- ipsych$fam$sample.ID %in% subset(rel, KINSHIP > 2^-3.5)$IID2
mean(is_rel) # 10.98%

mean(KEEP <- !is_rel & is_homogeneous) # 80.6%
ind_keep <- which(KEEP)

covariates_ipsych <- readRDS("../multi_prs/data/covariates_ipsych2015_2016j.rds") %>%
  mutate(age = 2020 - birth_year,
         is_rel = is_rel*1,
         is_homogeneous = is_homogeneous*1) %>%
  rename(sex = sexM) %>% select(-birth_year)

#saveRDS(covariates_ipsych, "data/covar_ipsych12_15.rds")
#write_delim(bind_cols(tibble(FID = ipsych$fam$family.ID, IID = ipsych$fam$sample.ID), covariates_ipsych), "data/covar_ipsych12_15.txt")

phen_file_15 <- read.csv("/project/Register/2019_06/csv/ipsych2015design_v2.csv") %>% select(pid, ends_with("I")) %>% 
  select(-skizospek2015I, -postpartum2015I) %>%
  mutate(across(c(2:7), ~ case_when(.x == 1 ~ 1,
                                    kontrol2015I == 1~ 0,
                                    TRUE ~ NA_real_))) %>%
  slice(match(ipsych$fam$family.ID, pid))
  
phen_file_12 <- read.csv("/project/Register/2016_03/design/ipsych2012design_expanded.csv") %>% select(pid, ends_with("I")) %>%
  mutate(across(c(2:7), ~ case_when(.x == 1 ~ 1,
                                    kontrol2012I == 1~ 0,
                                    TRUE ~ NA_real_))) %>%
  slice(match(ipsych$fam$family.ID, pid))

phen_file <- phen_file_15 %>% left_join(phen_file_12)

idx <- match(ipsych$fam$family.ID, phen_file$pid)

phen <- phen_file[idx,]
phen$MDD12 <- ifelse(phen$affek2012I == 1 & is.na(phen$bipol2012I), 1, ifelse(phen$affek2012I == 1 & !is.na(phen$bipol2012I), NA, phen$affek2012I))
phen$MDD15 <- ifelse(phen$affek2015I == 1 & is.na(phen$bipol2015I), 1, ifelse(phen$affek2015I == 1 & !is.na(phen$bipol2015I), NA, phen$affek2015I))


#saveRDS(phen, "data/phen_ipsych12_15.rds")

n_total <- phen %>%
  pivot_longer(2:17) %>%
  group_by(name, value) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = value, values_from = n, names_prefix = "N_")

n_eur_unrel <- phen[covariates_ipsych$is_homogeneous & !covariates_ipsych$is_rel,] %>%
  pivot_longer(2:17) %>%
  group_by(name, value) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = value, values_from = n, names_prefix = "N_")



# CV SUBSETS

# pars <- tibble(p = colnames(phen)[c(2:3, 5:7, 9:10, 12:14, 16:17)]) %>%
#   mutate(idx = map(p, ~which(!is.na(phen[[.x]]) & covariates_ipsych$is_homogeneous & !covariates_ipsych$is_rel)),
#          cv_sub = map(idx, ~sample(1:5, length(.x), replace = TRUE))) %>%
#   expand_grid(cv = 1:5) %>%
#   mutate(ind.train = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b != c]),
#          phen.train = map2(p, ind.train, ~phen[[.x]][.y]),
#          ind.test = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b == c]),
#          phen.test = map2(p, ind.test, ~phen[[.x]][.y])) %>%
#   select(-idx)

#saveRDS(pars, "data/cv_subsets_ipsych12_15.rds")

pars <- readRDS("data/cv_subsets_ipsych12_15.rds")

phen_cv <- pars %>% 
  mutate(v = list(rep(NA, nrow(cov))),
         i = pmap(list(v, ind.train, phen.train),  function(a,b,c) {
           a[b] <- c
           return(a)})) %>%
  select(p, cv, i) %>%
  unite(1:2, col = "p_cv") %>%
  pivot_wider(names_from = p_cv, values_from = i) %>%
  unnest() 

#write_delim(bind_cols(tibble(FID = ipsych$fam$family.ID, IID = ipsych$fam$sample.ID), phen_cv), "data/cv_subsets_ipsych12_15.txt")

pars_info <- readRDS("data/cv_subsets_ipsych12_15.rds") %>%
  mutate(N_int_case = map_dbl(phen.train, ~sum(.x == 1)),
         N_int_control = map_dbl(phen.train, ~sum(.x == 0))) %>%
  #select(p,cv, starts_with("N_int")) %>%
  mutate(gwas_int = paste0("results/gwas/", p, "_", cv, ".txt"),
         gwas_ext = case_when(str_detect(p, "ski") ~ "SCZ/PGC2/NoDK/daner_SCZ_noDenmark.metadaner.gz",
                              str_detect(p, "bip") ~ "BIP/PGC2/daner_PGC_BIP32b_mds7a_0416a.gz",
                              str_detect(p, "MDD") ~ "MDD/Howard/daner_howardUKB_no_iPsych.gz",
                              str_detect(p, "aut") ~ "ASD/daner_pgc_asd_euro_all_25Mar2015.txt.gz",
                              str_detect(p, "adhd") ~ "ADHD/PGC_all_euro_filtered_results.txt.gz",
                              str_detect(p, "anorek") ~ "ANOR/daner_Meta.without.danishANGI.gz"),
         gwas_ext = paste0("/faststorage/jail/project/XFiles/GWASsumstats/", gwas_ext), 
         gwas_meta = paste0("results/metal/", p, "_", cv, "1.txt"),
         N_ext_case = case_when(str_detect(p, "ski") ~ 21169,
                                str_detect(p, "bip") ~ 20040 ,
                                str_detect(p, "MDD") ~ 229897,
                                str_detect(p, "aut") ~ 5305,
                                str_detect(p, "adhd") ~ 4225,
                                str_detect(p, "anorek") ~ 11940),
         N_ext_control = case_when(str_detect(p, "ski") ~ 28117,
                                   str_detect(p, "bip") ~ 30874,
                                   str_detect(p, "MDD") ~ 544204,
                                   str_detect(p, "aut") ~ 5305,
                                   str_detect(p, "adhd") ~ 11012,
                                   str_detect(p, "anorek") ~ 33731),
         pop.prev = case_when(str_detect(p, "ski") ~ 0.01,
                              str_detect(p, "bip") ~ 0.01,
                              str_detect(p, "MDD") ~ 0.08,
                              str_detect(p, "aut") ~ 0.01,
                              str_detect(p, "adhd") ~ 0.05,
                              str_detect(p, "anorek") ~ 0.01))

#saveRDS(pars_info, "data/cv_subsets_ipsych12_15_info.rds")





# UKB 
phens <- readRDS("ukbb_traits/data/phenotypes_ukbb.rds")
subs_train <- readRDS("ukbb_traits/data/cv/ukbb_6pheno_5cv_train_subsets.rds")
subs_test <- readRDS("ukbb_traits/data/cv/ukbb_6pheno_5cv_test_subsets.rds")

pars <- tibble(p = names(phens)) %>%
  filter(p != "mdd_no_pilot") %>%
  expand_grid(cv = 1:5) %>%
  mutate(ind.train = map2(p, cv, ~which(!is.na(subs_train[[paste0(.x, "_train_cv", .y)]]))),
         phen.train = map2(p, ind.train, ~phens[[.x]][.y]),
         ind.test = map2(p, cv, ~which(!is.na(subs_test[[paste0(.x, "_train_cv", .y)]]))),
         phen.test = map2(p, ind.test, ~phens[[.x]][.y])) 

#saveRDS(pars, "data/cv_subsets_ukb_traits.rds")

pars_info <- readRDS("data/cv_subsets_ukb_traits.rds") %>%
  mutate(N_int_case = case_when(!p %in% c("height", "bmi") ~ map_dbl(phen.train, ~sum(.x == 1)),
                                p %in% c("height", "bmi") ~ map_dbl(phen.train, ~sum(!is.na(.x)))),
         N_int_control = case_when(!p %in% c("height", "bmi") ~ map_dbl(phen.train, ~sum(.x == 0)),
                                   p %in% c("height", "bmi") ~ 0)) %>%
  mutate(gwas_int = paste0("results/gwas/", p, "_", cv, ".txt"),
         gwas_ext = case_when(str_detect(p, "bc") ~ "ukbb_traits/data/sumstats/BC/29059683-GCST004988-EFO_0000305-build37.txt",
                              str_detect(p, "t2d") ~ "ukbb_traits/data/sumstats/T2D/METAANALYSIS_DIAGRAM_SE1_LDPRED_snp.txt",
                              str_detect(p, "cad") ~ "ukbb_traits/data/sumstats/CAD/26343387-GCST003116-EFO_0000378-build37.f_snp.tsv",
                              str_detect(p, "mdd_all") ~ "ukbb_traits/data/sumstats/MDD/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz",
                              str_detect(p, "height") ~ "ukbb_traits/data/sumstats/height/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_CHR_POS.txt",
                              str_detect(p, "bmi") ~ "ukbb_traits/data/sumstats/BMI/SNP_gwas_mc_merge_nogc.tbl.uniq_GIANT.txt"),
         gwas_meta = paste0("ukbb_traits/results/metal/metal_output_", p, "_cv", cv, ".txt"),
         N_ext_case = case_when(str_detect(p, "bc") ~ 12024,
                                str_detect(p, "t2d") ~ 18857 ,
                                str_detect(p, "cad") ~ 11529,
                                str_detect(p, "mdd_all") ~ 28626,
                                str_detect(p, "height") ~ 336750/2,
                                str_detect(p, "bmi") ~ 336381/2),
         N_ext_control = case_when(str_detect(p, "bc") ~ 169207,
                                   str_detect(p, "t2d") ~ 318618 ,
                                   str_detect(p, "cad") ~ 325946,
                                   str_detect(p, "mdd_all") ~ 308849,
                                   str_detect(p, "height") ~ 336750/2,
                                   str_detect(p, "bmi") ~ 336381/2),
         pop.prev = case_when(str_detect(p, "bc") ~ 0.07,
                              str_detect(p, "t2d") ~ 0.05 ,
                              str_detect(p, "cad") ~ 0.03,
                              str_detect(p, "mdd_all") ~ 0.15,
                              str_detect(p, "height") ~ NA_real_,
                              str_detect(p, "bmi") ~ NA_real_)) 


#saveRDS(pars_info, "data/cv_subsets_ukb_traits_info.rds")

