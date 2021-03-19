library(tidyverse)
library(bigsnpr)
library(glue)

cov_ipsych <- readRDS("data/covar_ipsych12_15.rds") %>% select(-is_rel, -is_homogeneous)
cov_ukb <- data.table::fread("UKBB/data/covariates_ukb.txt") %>% select(3,4,8,10:29) %>% mutate(sex = as.factor(sex))


prs <- bind_rows(readRDS("data/cv_subsets_ipsych12_15_info.rds") %>% mutate(data = "iPSYCH"),
                        readRDS("data/cv_subsets_ukb_traits_info.rds") %>% mutate(data = "UKB")) %>%
  select(-cv_sub, -ind.train, -phen.train) %>%
  rename(ind.prs = ind.test, y = phen.test) %>%
  # make subsets
  mutate(ind.val = map(ind.prs, ~sample(1:length(.x), length(.x)*0.5)),
         phen.val = map2(ind.val, y, ~.y[.x]),
         ind.test = map2(ind.prs, ind.val, ~setdiff(1:length(.x), .y)),
         phen.test = map2(ind.test, y, ~.y[.x]),
         cov = case_when(data == "iPSYCH" ~map(ind.prs, ~cov_ipsych[.x,]),
                         data == "UKB" ~map(ind.prs, ~cov_ukb[.x,])),
         cov = case_when(p == "bc" ~ map(cov, ~.x %>% select(-sex)),
                         TRUE ~ cov)) %>%
  # import prs
  mutate(prs_bolt = glue("results/prs/{p}_{cv}_bolt.rds"),
         prs_ldpred_int = glue("results/prs/{p}_{cv}_ldpred_int.rds"),
         prs_ldpred_ext = glue("results/prs/{p}_{cv}_ldpred_ext.rds"),
         prs_ldpred_meta = glue("results/prs/{p}_{cv}_ldpred_meta.rds"),
         prs_sct = glue("results/prs/{p}_{cv}_sct.rds")) %>%
  filter(file.exists(prs_bolt) & file.exists(prs_ldpred_int) & file.exists(prs_ldpred_ext) & file.exists(prs_ldpred_meta) & file.exists(prs_sct)) %>%
  mutate(across(starts_with("prs_"), ~map(.x, readRDS))) %>%
  # separate ct and ldpred scores
  mutate(across(contains("ldpred"), ~map(.x, ~.x[,15:ncol(.x)]), .names = "ld_{.col}")) %>%
  mutate(across(starts_with("prs_ldpred"), ~map(.x, ~.x[,1:14]), .names = "ct_{.col}")) %>%
  select(-prs_ldpred_meta, -prs_ldpred_ext, -prs_ldpred_int) %>%
  pivot_longer(contains("ldpred"), names_to = c("method_ss", ".value"),
               names_pattern = "(.+)_(prs_ldpred_.+)") %>%
  # select p for ldpred an
  mutate(max_ldpred_ext = pmap_dbl(list(prs_ldpred_ext, ind.val, phen.val), function(prs, idx, y) which.max(cor(y, prs[idx,])^2)),
         max_ldpred_meta = pmap_dbl(list(prs_ldpred_meta, ind.val, phen.val), function(prs, idx, y) which.max(cor(y, prs[idx,])^2)),
         max_ldpred_int = pmap_dbl(list(prs_ldpred_int, ind.val, phen.val), function(prs, idx, y) which.max(cor(y, prs[idx,])^2)),
         # re-arrange data for lm
         data_val = pmap(list(prs_bolt, prs_ldpred_ext, phen.val, ind.val, max_ldpred_ext, prs_ldpred_meta, max_ldpred_meta), function(a,b,c,d, i, m, j) tibble(y = c, prs_int = scale(a[d])[,1], prs_ext = scale(b[d,i])[,1], metagwas = scale(m[d,j])[,1])),
         data_test = pmap(list(prs_bolt, prs_ldpred_ext, phen.test, ind.test, max_ldpred_ext, prs_ldpred_meta, max_ldpred_meta), function(a,b,c,d, i, m, j) tibble(y = c, prs_int = scale(a[d])[,1], prs_ext = scale(b[d,i])[,1],  metagwas = scale(m[d,j])[,1])),
         model_lm = map(data_val, ~lm(y ~ . -metagwas, data = .x)),
         model_gwas = map(data_val, ~lm(y ~ . -prs_ext, data = .x)),
         prs_metaprs.lm = map2(data_test, model_lm, ~predict(.y, newdata = .x)),
         prs_metaprs.plus = map2(data_test, model_gwas, ~predict(.y, newdata = .x)),
         Neff_int = 4/(1/N_int_case + 1/N_int_control),
         Neff_ext = 4/(1/N_ext_case + 1/N_ext_control),
         Neff_int = case_when(p %in% c("height", "bmi") ~ N_int_case,
                              TRUE ~ Neff_int),
         Neff_ext = case_when(p %in% c("height", "bmi") ~ N_ext_case + N_ext_control,
                              TRUE ~ Neff_ext),
         prs_metaprs.neff = pmap(list(data_test, Neff_int, Neff_ext), function(d, ni, ne) (d$prs_int*sqrt(ni) + d$prs_ext*sqrt(ne)) / sqrt(ni + ne)),
         prs_ldpred_int = map2(prs_ldpred_int, max_ldpred_int, ~.x[,.y]),
         prs_ldpred_ext = map2(prs_ldpred_ext, max_ldpred_ext, ~.x[,.y]),
         prs_ldpred_meta = map2(prs_ldpred_meta, max_ldpred_meta, ~.x[,.y])) %>%
  mutate(across(.cols = c(prs_bolt, prs_ldpred_int, prs_ldpred_ext, prs_ldpred_meta, prs_sct), .fns = ~map2(.x, ind.test, function(a,b) a[b]))) %>%
  mutate(cov = map2(cov, ind.test, ~.x[.y,])) %>%
  select(p, cv, pop.prev, method_ss, phen.test, cov, starts_with("prs"), prs_metaprs.lm, prs_metaprs.plus, prs_metaprs.neff, Neff_int, Neff_ext, N_int_case, N_int_control, N_ext_case, N_ext_control) %>%
  rename(prs_metagwas = prs_ldpred_meta, prs_ldpred.int = prs_ldpred_int, prs_ldpred.ext = prs_ldpred_ext) 

#saveRDS(prs,  "results/prs/all_prs_real.rds")         

bootstrap <- function(x) {
  boots_m <- c()
  boots_sd <- c()
  if (!is.na(x)) {
    for (i in 1:1e4) {
      boots_m[i] <- mean(sample(x, replace = TRUE))
      boots_sd[i] <- sd(sample(x, replace = TRUE))}
    return(list(mean = mean(boots_m), sd = mean(boots_sd)))
  }
  }

pred_r2 <- prs %>%
  mutate(across(contains("prs"), .fns = ~pmap_dbl(list(phen.test, .x, cov), function(prs, y, cov) pcor(prs, y, cov)[1])^2, .names = "{.col}_r2")) %>%
  rename_with(~str_sub(.x, start = 5), contains("prs")) %>%
  group_by(p, pop.prev, rg, method_ss) %>%
  summarise(across(contains("r2"), list),
            Neff_int = round(mean(Neff_int), 0),
            Neff_ext = mean(Neff_ext)) %>%
  mutate(across(contains("r2"), ~map(.x, bootstrap))) %>%
  pivot_longer(cols = c(everything(), -p, -method_ss, -pop.prev, -Neff_int, -Neff_ext)) %>%
  unnest_wider(value) %>%
  separate(name, sep = "_", into = c("method", "pred")) %>%
  separate(method, c("method", "sub")) %>%
  filter(!(method == "bolt" & method_ss == "ld") & !(method == "sct" & method_ss == "ld")) %>%
  mutate(method_ss = case_when(!method %in% c("bolt","sct") ~ method_ss, TRUE ~ NA_character_)) %>%
  mutate(data = case_when(method == "bolt" | (method == "ldpred" & sub == "int") ~ "int",
                          method == "ldpred" & sub == "ext" ~ "ext",
                          TRUE ~ "ext+int")) 

pred_auc <- prs %>%
  filter(!p %in% c("bmi", "height")) %>%
  mutate(across(contains("prs") & !contains("_r2"), .fns = ~ map2_dbl(.x, phen.test, ~AUC(.x, .y)), .names = "{.col}_auc")) %>%
  rename_with(~str_sub(.x, start = 5), contains("prs")) %>%
  group_by(p, pop.prev, rg, method_ss) %>%
  summarise(across(contains("auc"), list),
            Neff_int = round(mean(Neff_int), 0),
            Neff_ext = mean(Neff_ext)) %>%
  mutate(across(contains("auc"), ~map(.x, bootstrap))) %>%
  pivot_longer(cols = c(everything(), -p, -method_ss, -pop.prev, -Neff_int, -Neff_ext)) %>%
  unnest_wider(value) %>%
  separate(name, sep = "_", into = c("method", "pred")) %>%
  separate(method, c("method", "sub")) %>%
  filter(!(method == "bolt" & method_ss == "ld") & !(method == "sct" & method_ss == "ld")) %>%
  mutate(method_ss = case_when(!method %in% c("bolt","sct") ~ method_ss, TRUE ~ NA_character_)) %>%
  mutate(data = case_when(method == "bolt" | (method == "ldpred" & sub == "int") ~ "int",
                          method == "ldpred" & sub == "ext" ~ "ext",
                          TRUE ~ "ext+int")) 

#saveRDS(pred_r2, "results/prediction_real_r2.rds")
#saveRDS(pred_auc, "results/prediction_real_auc.rds")
