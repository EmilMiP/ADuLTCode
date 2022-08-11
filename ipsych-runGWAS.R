# Loading packages --------------------------------------------------------
library(bigsnpr)
library(data.table)
library(tidyverse)
library(foreach)

# Defining Functions ------------------------------------------------------

save_run <- function(code, file, timing = TRUE) {
  file <- path.expand(file)
  bigassertr::assert_ext(file, "rds")
  bigassertr::assert_dir(dirname(file))
  if (file.exists(file)) {
    readRDS(file)
  } else {
    time <- system.time(res <- code)
    if (timing) print(time)
    saveRDS(res, file)
    res
  }
}


# Loading Genetic Data ----------------------------------------------------
obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS

intervals <- bigparallelr::split_len(ncol(G), nb_split = 20)

covar = fread("covariates/ipsych-covariates.txt")
#### covar = as.matrix(covar[,c("age", "sex", "is_2012", paste0("PC", 1:20))])
covar = as.matrix(covar[,c("sex", "is_2012", paste0("PC", 1:20))])
unrel_homogen_indivs = fread("./covariates/ipsych-unrelated-homogeneous-thirdDegree-indivs.txt")

# Loading Phenotypes ------------------------------------------------------

info = tibble(
  disorder_name = c("adhd", "asd", "dep", "scz"),
  phenotype     = paste0("./phenotypes/", disorder_name, ".phen"),
  to_exclude    = paste0("./covariates/ipsych-", disorder_name, "-excluded-thirdDegree-indivs.txt")
)

pheno = lapply(1:nrow(info), function(ctr) {
  disorder = info$disorder_name[ctr]
  fread(info$phenotype[ctr]) %>% 
    as_tibble() %>% 
    select(FID, IID, status, post_gen_no_fam) %>% 
    rename(!!as.symbol(disorder) := post_gen_no_fam,
           !!as.symbol(paste0("status_", disorder)) := status) 
}) %>% purrr::reduce(., left_join) %>% 
  distinct(FID, IID, .keep_all = T) %>% 
  left_join(obj.bigsnp$fam[,1:2], ., by = c("family.ID" = "FID", "sample.ID" = "IID")) %>% 
  as_tibble() %>% 
  rename(FID = family.ID,
         IID = sample.ID)


# Performing GWAS ---------------------------------------------------------

library(future.batchtools)
NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))


res_files <- paste0("./gwas-results/ipsych_noAge_thirdDegree_part", rows_along(intervals), ".rds")
mean(already_done <- file.exists(res_files))

future.apply::future_lapply(which(!already_done), function(ic) {
  
  ind <- seq(intervals[ic, "lower"], intervals[ic, "upper"])
  all_gwas_part <- lapply(seq_along(pheno[,-c(1,2)]), function(n) {
    
    #get the raw phenotype vector
    y = pheno[[ 2 + n ]]
    #get the name of the phenotype
    cur_name = names(pheno[,-c(1,2)])[n] %>% 
      str_replace("status_", "")
    #find the controls and cases only
    to_remove = fread(paste0("./covariates/ipsych-", cur_name, "-excluded-thirdDegree-indivs.txt")) %>% pull(FID)
    to_remove_ind = pheno$FID %in% to_remove
    #perform GWAS on homogen etc. indivs.
    ind_keep <- which(!is.na(y) & !to_remove_ind)
    gwas <- bigstatsr::big_univLinReg(
      G, y[ind_keep], ind.train = ind_keep, ind.col = ind,
      covar.train = covar[ind_keep, ], ncores = NCORES)
    # if(str_detect(names(pheno)[-c(1,2)][n], "status")) table(y[ind_keep])
  })
  names(all_gwas_part) = names(pheno)[-c(1,2)]
  saveRDS(all_gwas_part, res_files[ic])
}, future.seed = T) # ~20 min for each part


