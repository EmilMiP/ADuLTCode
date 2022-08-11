# Loading packages --------------------------------------------------------
library(bigsnpr)
library(data.table)
library(tidyverse)
library(foreach)
library(SPACox)

info = tibble(
  disorder_name = c("adhd", "asd", "dep", "scz"),
  phenotype     = paste0("./phenotypes/", disorder_name, ".phen"),
  to_exclude    = paste0("./covariates/ipsych-", disorder_name, "-excluded-thirdDegree-indivs.txt")
)

# SPACox GWAS -------------------------------------------------------------


#get and save SPACox GWAS in iPSYCH
all_phenotypes = tibble(files = list.files("./gwas-results", full.names = T, "ipsych_SPACox")) %>% 
  mutate(phenotype = as.vector(str_match(files, paste0(info$disorder_name, collapse = "|")))) %>% 
  filter(str_detect(files, "noAge")) %>% 
  group_by(phenotype) %>% 
  summarise(all_parts = list(files))

lapply(1:nrow(all_phenotypes), function(ctr) {
  print(ctr)
  raw_part_nbs = sapply(str_split(all_phenotypes$all_parts[[ctr]], "part|\\.rds"), function(x) x[2]) %>% 
    as.numeric()
  comb_gwas = lapply(all_phenotypes$all_parts[[ctr]][order(raw_part_nbs)], readRDS) %>% 
    lapply(function(x) {
      as_tibble(x[[1]])
    }) %>% bind_rows()
  saveRDS(comb_gwas,
          paste0("./results/ipsych_SPACox_", all_phenotypes$phenotype[ctr], "_noAge_thirdDegree.rds"))
})



# LTM & linReg GWAS -------------------------------------------------------

all_gwas <- save_run(
  {
    all_res_files <- paste0("./gwas-results/ipsych_noAge_thirdDegree_part", rows_along(intervals), ".rds")
    all_gwas_parts <- foreach(res_file = all_res_files) %do% readRDS(res_file)
    lapply(purrr::transpose(all_gwas_parts), function(x) do.call("rbind", x))
  },
  file = "./results/ipsych_gwas_noAge_thirdDegree.rds"
)

