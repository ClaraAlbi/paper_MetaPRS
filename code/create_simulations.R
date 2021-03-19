library(tidyverse)
library(mvtnorm)
library(bigsnpr)
library(future.batchtools)

# Simulations
G <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$genotypes
M <- nrow(G)
h2 <- 0.5

sims <- expand_grid(p = c(1e4, 1e5),
                    rep = 1:5,
                    h2 = 0.5,
                    rg = c(0.5, 0.8, 1)) 

scaled_prod <- function(X, ind, ind.row, ind.col, beta) {
  ms <- bigstatsr::big_scale()(X, ind.row, ind.col[ind])
  bigstatsr::big_prodVec(X, beta[ind], ind.row, ind.col[ind],
              center = ms$center, scale = ms$scale)
}

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(sims, function(p, rep, h2, rg) {
  set <- sample(cols_along(G), size = p)
  effects <- sapply(1:length(set), function(x) rmvnorm(n = 1, mean = c(0,0), matrix(c(h2/M, rg*h2/M, rg*h2/M, h2/M), nrow = 2)))

  phenos <- sapply(1:2, function(t)  {
    e = effects[t,]
    gen_liab <- big_apply(G, scaled_prod, a.combine = bigparallelr::plus, 
                            ind = seq_along(set), ind.col = set, ind.row = rows_along(G), 
                            beta = e, ncores = 24)
    coeff1 <- sqrt(h2)/stats::sd(gen_liab)
    gen_liab <- gen_liab * coeff1
    env_liab <- stats::rnorm(length(gen_liab), sd = sqrt(1 - h2))
    var_env <- stats::var(env_liab)
    cov_env <- stats::cov(gen_liab, env_liab)
    coeff2 <- (sqrt(cov_env^2 + (1 - h2) * var_env) - cov_env)/var_env
    full_liab <- gen_liab + env_liab * coeff2
    (full_liab > stats::qnorm(1 - 0.2)) + 0L
    })
  saveRDS(tibble(p = p, rg = rg, rep = rep, set = list(set), effects_1 = list(effects[1,]),
                 effects_2 = list(effects[2,]),
                 phen_1 = list(phenos[,1]),
                 phen_2 = list(phenos[,2])), glue::glue("data/sims/p{p/1000}k_rg{rg}_{rep}.rds"))

})

sims_f <- list.files("data/sims", full.names = TRUE)
library(bigparallelr)
registerDoParallel(cl <- makeCluster(24))
res <- foreach(res_file = sims_f) %dopar% readRDS(res_file)
res_l <- do.call(bind_rows, res)
#saveRDS(res_l, "data/sims.rds")

# Expand sims
sims <- readRDS("data/sims.rds")

sims_ind <- sims %>%
  expand_grid(part_ind = c(0.1, 0.25, 0.5, 0.75, 0.9),
              size = c(50e3, M)) %>%
  mutate(part_ss = 1 - part_ind,
         ind = map(size, ~sample(1:M, .x)),
         train = map2(ind, size, ~sample(.x, .y*0.9)),
         test = map2(ind, train, ~sample(setdiff(.x, .y), length(.x)*0.05)),
         val = pmap(list(ind, train, test), function(a,b,c) setdiff(a, c(b, c))),
         train_ind = map2(train, part_ind, ~sample(.x, length(.x)*.y)),
         train_ss = map2(train, train_ind, ~ setdiff(.x, .y)),
         phen_ind = map2(phen_1, train_ind, ~.x[.y]),
         phen_ss = map2(phen_2, train_ss, ~.x[.y]),
         Neff_ind = map_dbl(phen_ind, ~ 4/(1/(sum(.x == 0)) + 1/(sum(.x == 1)))),
         Neff_ss = map_dbl(phen_ss, ~ 4/(1/(sum(.x == 0)) + 1/(sum(.x == 1)))))

#saveRDS(sims_ind, "data/sims_ind.rds")
 
sims_ind <- readRDS("data/sims_ind.rds")

sims_ind_bolt <- sims_ind %>% 
  mutate(v = list(rep(NA, M)),
         i = pmap(list(v, train_ind, phen_ind),  function(a,b,c) {
           a[b] <- c
           return(a)})) %>%
  select(size, p, rg, rep, part_ind, i) %>%
  mutate(p = p/1000, size = round(size/1000, 0)) %>%
  unite(1:5, col = "size_p_rg_rep_part") %>%
  pivot_wider(names_from = size_p_rg_rep_part, values_from = i) %>%
  unnest() 

#saveRDS(sims_ind_bolt, "data/sims_ind_all.rds")
fam <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")$fam
#write_delim(bind_cols(tibble(FID = fam$family.ID, IID = fam$sample.ID), sims_ind_bolt), "data/sims_ind.txt")
