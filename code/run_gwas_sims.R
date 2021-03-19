library(tidyverse)
library(bigsnpr)
library(future.batchtools)


gdata <- snp_attach("UKBB/data/UKBB_imp_HM3.rds")
sims <- readRDS("data/sims_ind.rds") %>% mutate(p = p/1000, size = round(size/1000, 0))



# Run GWAS on part_ind and part_ss separately

sims_gwas <- sims %>%
  select(size, p, rg, rep, part_ind, part_ss, train_ind, phen_ind, train_ss, phen_ss) %>%
  pivot_longer(5:10, names_to = c(".value", "type"),
               names_pattern = "(.+)_(.+)") %>%
  filter(file.exists(glue::glue("results/gwas/{size}_{p}_{rg}_{rep}_{type}_{part}.rds"))) 

NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(sims_gwas[,1:8], function(size, p, rg, rep, type, part, train, phen) {

  gwas <- bigstatsr::big_univLinReg(
    gdata$genotypes, phen, ind.train = train, ncores = 24)
  
  saveRDS(gwas, glue::glue("results/gwas/{size}_{p}_{rg}_{rep}_{type}_{part}.rds"))
})


map <- gdata$map %>% rename(CHR = chromosome, SNP = rsid, POS = physical.pos, A1 = allele1, A2 = allele2)

pmap(sims_gwas[,1:6], function(size, p, rg, rep, type, part) {
  file <- readRDS(glue::glue("results/gwas/{size}_{p}_{rg}_{rep}_{type}_{part}.rds"))
  pval <- predict(file, log10 = FALSE)
  file_txt <- bind_cols(map[,c(1, 3:6)], file, tibble(pval = pval)) %>%
    rename(BETA = estim, SE = std.err, P = pval) %>% select(-score)
  write_delim(file_txt, glue::glue("results/gwas/{size}_{p}_{rg}_{rep}_{type}_{part}.txt"))
})



# BOLT-LMM
NCORES <- 24
plan(batchtools_slurm(resources = list(
  t = "00-24:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

sims_bolt <- sims %>% filter(file.size(glue::glue("results/bolt/{size}_{p}_{rg}_{rep}_{part_ind}.pred")) == 0 | !file.exists(glue("results/bolt/{size}_{p}_{rg}_{rep}_{part_ind}.pred"))) %>%
  filter(size == 337)

furrr::future_pmap(sims_bolt[,c(1:3, 9:10)], function(size, p, rg, rep, part_ind) {
  phen <- glue::glue("{size}_{p}_{rg}_{rep}_{part_ind}")
  system(glue::glue(
    "/home/clara/REPOS/BOLT-LMM_v2.3.4/bolt",
    " --bfile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/UKBB/data/UKBB_imp_HM3",
    " --lmm",
    " --maxModelSnps 1200000",
    " --phenoFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/data/sims_ind.txt",
    " --phenoCol={phen}",
    " --LDscoresFile=/home/clara/REPOS/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz",
    " --geneticMapFile=/home/clara/REPOS/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz",
    " --numThreads=24",
    " --statsFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.stats",
    " --predBetasFile=/faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.pred",
    " --verboseStats",
    " &>> /faststorage/jail/project/NCRR-PRS/clara/paper_prs/results/bolt/{phen}.out"
  ))
})


# Meta-analysis

NCORES <- 1
plan(batchtools_slurm(resources = list(
  t = "0-01:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

sims <- readRDS("data/sims_ind.rds") %>% mutate(p = p/1000, size = round(size/1000, 0))  %>%
  mutate(gwas_int = glue("results/gwas/{size}_{p}_{rg}_{rep}_ind_{part_ind}.txt"),
         gwas_ext = glue("results/gwas/{size}_{p}_{rg}_{rep}_ss_{part_ss}.txt"),
         N_int_case = map_dbl(phen_ind, ~ sum(.x == 1)), 
         N_int_control = map_dbl(phen_ind, ~ sum(.x == 0)), 
         N_ext_case = map_dbl(phen_ss, ~ sum(.x == 1)),
         N_ext_control = map_dbl(phen_ss, ~ sum(.x == 0))) %>% filter(size == 50)

furrr::future_pmap(sims[,c(1:3, 9:11, 22:27)], function(size, p, rg, rep, part_ind, part_ss, gwas_int, gwas_ext, N_int_case, N_int_control, N_ext_case, N_ext_control) {
  outfile <- glue("results/metal/{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}")
  run_metal <- glue("results/metal/run_metal{size}_{p}_{rg}_{rep}_{part_ind}_{part_ss}.txt")
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
PROCESS {gwas_int}
PROCESS {gwas_ext}
OUTFILE {outfile} .txt
ANALYZE"), run_metal)
  
  system(glue::glue("~/REPOS/METAL/build/bin/metal {run_metal}"))
})
