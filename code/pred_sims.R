library(tidyverse)
library(bigsnpr)
library(glue)

sims <- readRDS("data/sims_ind.rds") %>% mutate(p = p/1000, size = round(size/1000, 0))

prs_sims <- sims %>%
  # make subsets
  mutate(phen.val = map2(val, phen_1, ~.y[.x]),
         phen.test = map2(test, phen_1, ~.y[.x])) %>%
  select(-c(4:8), -ind, -train, - train_ind, -train_ss, -phen_ind, -phen_ss) %>%
  filter(size == 50) %>%
  # import prs
  mutate(prs_bolt = glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_bolt.rds"),
         prs_ldpred_ext = glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ldpred_ext.rds"),
         prs_ldpred_meta = glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ldpred_meta.rds"),
         prs_ldpred_int = glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_ldpred_int.rds"),
         prs_sct = glue("results/prs/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}_sct.rds")) %>%
  filter(file.exists(prs_bolt) & file.exists(prs_ldpred_ext) & file.exists(prs_ldpred_meta) & file.exists(prs_ldpred_int) & file.exists(prs_sct)) %>%
  mutate(across(starts_with("prs_"), ~map(.x, readRDS))) %>%
  # separate ct and ldpred scores
  mutate(across(contains("ldpred"), ~map(.x, ~.x[,15:ncol(.x)]), .names = "ld_{.col}")) %>%
  mutate(across(starts_with("prs_ldpred"), ~map(.x, ~.x[,1:14]), .names = "ct_{.col}")) %>%
  select(-prs_ldpred_meta, -prs_ldpred_ext, -prs_ldpred_int)

pred <- prs_sims %>% 
  pivot_longer(contains("ldpred"), names_to = c("method_ss", ".value"),
               names_pattern = "(.+)_(prs_ldpred_.+)") %>%
  mutate(max_ldpred_ext = map2_dbl(prs_ldpred_ext, phen.val, ~which.max(cor(.y, .x[1:2500,])^2)),
         max_ldpred_meta = map2_dbl(prs_ldpred_meta, phen.val, ~which.max(cor(.y, .x[1:2500,])^2)),
         max_ldpred_int = map2_dbl(prs_ldpred_int, phen.val, ~which.max(cor(.y, .x[1:2500,])^2)),
         data_val = pmap(list(prs_bolt, prs_ldpred_ext, phen.val, max_ldpred_ext, prs_ldpred_meta, max_ldpred_meta), function(a,b,y,i,c,j) tibble(y = y, prs_int = scale(a[1:2500])[,1], prs_ext = scale(b[1:2500,i])[,1], metagwas = scale(c[1:2500,j])[,1])),
         data_test = pmap(list(prs_bolt, prs_ldpred_ext, phen.test, max_ldpred_ext, prs_ldpred_meta, max_ldpred_meta), function(a,b,y, i,c,j) tibble(y = y, prs_int = scale(a[2501:5000])[,1], prs_ext = scale(b[2501:5000,i])[,1], metagwas = scale(c[2501:5000,j])[,1])),
         model_lm = map(data_val, ~lm(y ~ . -metagwas, data = .x)),
         model_gwas = map(data_val, ~lm(y ~ . -prs_ext, data = .x)),
         prs_metaprs.lm = map2(data_test, model_lm, ~predict(.y, newdata = .x)),
         prs_metaprs.plus = map2(data_test, model_gwas, ~predict(.y, newdata = .x)),
         prs_ldpred_ext = map2(prs_ldpred_ext, max_ldpred_ext, ~.x[,.y]),
         prs_metaprs.neff = pmap(list(data_test, Neff_ind, Neff_ss), function(d, ni, ne) d$prs_int*sqrt(ni) + d$prs_ext*sqrt(ne)),
         prs_ldpred_meta = map2(prs_ldpred_meta, max_ldpred_meta, ~.x[,.y]),
         prs_ldpred_int = map2(prs_ldpred_int, max_ldpred_int, ~.x[,.y])) %>%
  mutate(across(.cols = c(prs_bolt, prs_ldpred_ext, prs_ldpred_int, prs_ldpred_meta, prs_sct), .fns = ~map(.x, ~.x[2501:5000]))) %>%
  rename(prs_metagwas = prs_ldpred_meta, prs_ldpred.ext = prs_ldpred_ext, prs_ldpred.int = prs_ldpred_int) %>%
  mutate(across(contains("prs"), .fns = ~map2_dbl(phen.test, .x, ~ summary(lm(.x ~ .y))$adj.r.squared), .names = "{.col}_r2")) %>%
  rename_with(~str_sub(.x, start = 5), contains("prs"))  %>%
  group_by(p, rg, part_ind, part_ss, method_ss) %>%
  summarise(across(contains("r2"), list), across(contains("auc"), list))  %>%
  ungroup() %>%
  mutate(across(contains("r2"), ~map(.x, bootstrap)), across(contains("auc"), ~map(.x, bootstrap))) %>%
  pivot_longer(cols = c(everything(), -p, -rg, -method_ss, -part_ind, -part_ss)) %>%
  unnest_wider(value) %>%
  separate(name, sep = "_", into = c("method", "pred")) %>%
  separate(method, c("method", "sub")) %>%
  filter(!(method == "bolt" & method_ss == "ld") & !(method == "sct" & method_ss == "ld")) %>%
  mutate(method_ss = case_when(!method %in% c("bolt","sct") ~ method_ss, TRUE ~ NA_character_))
  
#saveRDS(pred, "results/prediction_sims.rds")

