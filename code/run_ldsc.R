library(tidyverse)
library(GenomicSEM)
library(future.batchtools)

pars <- bind_rows( readRDS("data/cv_subsets_ipsych12_15_info.rds") , readRDS("data/cv_subsets_ukb_traits_info.rds")) %>%
  mutate(N_meta_case = N_int_case + N_ext_case,
         N_meta_control = N_int_control + N_ext_control) 


NCORES <- 4
plan(batchtools_slurm(resources = list(
  t = "0-02:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(pars[,c(1,3, 8:17)], function(p, cv, gwas_int, gwas_ext, gwas_meta, N_int_case, N_int_control, N_ext_case, N_ext_control, N_meta_case, N_meta_control, pop.prev) {
  
  trait.names <- c(paste0(p, "_", cv, c("_ext", "_int", "_meta")))
  traits <- paste0(trait.names, ".sumstats.gz")
  
  munge(files = c(gwas_ext, gwas_int, gwas_meta),
        hm3 = "~/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/w_hm3.snplist",
        trait.names = trait.names,
        N = c(N_ext_case + N_ext_control, N_int_case + N_int_control, N_meta_case + N_meta_control))
  
  LDSCOutput <- ldsc(traits = traits, 
                     sample.prev = c(N_ext_case / (N_ext_case + N_ext_control), N_int_case / (N_int_case + N_int_control), N_meta_case / (N_meta_case + N_meta_control)),
                     population.prev = c(pop.prev, pop.prev, pop.prev),
                     ld = "~/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/",
                     wld = "~/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/",
                     trait.names = trait.names, 
                     stand = TRUE)
  saveRDS(LDSCOutput, paste0("results/ldsc/", p, "_", cv, ".rds"))
  
})

# Results

res_cor <- pars %>% 
    mutate(log = map2(p, cv, ~read_lines(paste0("results/ldsc/", .x, "_", .y, "_ext",".sumstats.gz_", .x, "_", .y, "_int", ".sumstats.gz", "_", .x, "_", .y, "_meta",".sumstats.gz", "_ldsc.log"))),
           interc_ext = map_chr(log, ~str_match(.x[str_detect(.x, "Intercept:")][1], ": (.*)")[2]),
           h2_ext = map_chr(log, ~str_match(.x[str_detect(.x, "Total Observed Scale h2")][1], ": (.*)")[2]),
           interc_int = map_chr(log, ~str_match(.x[str_detect(.x, "Intercept:")][4], ": (.*)")[2]),
           h2_int = map_chr(log, ~str_match(.x[str_detect(.x, "Total Observed Scale h2")][2], ": (.*)")[2]),
           interc_int_ext = map_chr(log, ~str_match(.x[str_detect(.x, "Cross trait Intercept")][1], ": (.*)")[2]),
           rg_int_ext = pmap_chr(list(log, p, cv), function(f, a, b) str_match(f[str_detect(f, glue("Genetic Correlation between {a}_{b}_ext and {a}_{b}_int"))], ": (.*)")[2])) %>%
  select(1,3, 8:9, 13:15, 19:24) 

#saveRDS(res_cor, "exports/export_3/ldsc_rg.rds")


