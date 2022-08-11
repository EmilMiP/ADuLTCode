# Loading packages --------------------------------------------------------
library(bigsnpr)
library(data.table)
library(tidyverse)
library(foreach)
library(SPACox)

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
##### covar = covar[,c("age", "sex", "is_2012", paste0("PC", 1:20))]
covar = covar[,c("sex", "is_2012", paste0("PC", 1:20))]
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
    mutate(event = pmin(age, AOO, na.rm = T)) %>% 
    select(FID, IID, event, status) %>% 
    rename(!!as.symbol(paste0("event_", disorder)) := event,
           !!as.symbol(paste0("status_", disorder)) := status) 
}) %>% purrr::reduce(., left_join) %>% 
  distinct(FID, IID, .keep_all = T) %>% 
  left_join(obj.bigsnp$fam[,1:2], ., by = c("family.ID" = "FID", "sample.ID" = "IID")) %>% 
  as_tibble() %>% 
  rename(FID = family.ID,
         IID = sample.ID)


# Performing GWAS ---------------------------------------------------------

library(future.batchtools)
NCORES <- 20

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))



future.apply::future_lapply(1:nrow(info), function(ctr) {
  library(survival)
  disorder = info$disorder_name[ctr]
  
  res_files <- paste0("./gwas-results/ipsych_SPACox_", disorder, "_noAge_thirdDegree_part", rows_along(intervals), ".rds")
  mean(already_done <- file.exists(res_files))
  #find the controls and cases only
  to_remove = fread(paste0("./covariates/ipsych-", disorder, "-excluded-thirdDegree-indivs.txt")) %>% pull(FID)
  
  #perform GWAS on homogen etc. indivs.
  ind_keep <- which(!(pheno$FID %in% to_remove))
  # extract individuals  
  phen = select(pheno, FID, IID, contains(disorder)) %>% 
    bind_cols(covar) %>% 
    slice(ind_keep)
  
  
  
  ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
  #SPACox null model is fitted:
  null_mod = SPACox_Null_Model(formula = as.formula(paste0("Surv(" ,"event_", disorder, ",",
                                                    "status_", disorder, ") ~ ", 
                                                    paste0("PC", 1:20, collapse = " + "), 
                                                    " + sex + is_2012")),  ####  + age
                               data = phen,
                               pIDs = phen$IID,
                               gIDs = phen$IID)
  ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
  
  library(future.apply)
  plan(tweak(multisession, workers = NCORES))
  
  future.apply::future_lapply(1:nrow(intervals), function(ic) {
      
    ind <- seq(intervals[ic, "lower"], intervals[ic, "upper"])
    
    #SPAcox association test
    scox = big_apply(X = G, a.FUN = function(X, ind, null_mod, phen, keep_indivs) {
      library(SPACox)
      ph = X[keep_indivs,ind]
      rownames(ph) = phen$IID
      colnames(ph) = phen$IID[ind]
      SPACox(null_mod, ph)
    }, ncores = 1, null_mod = null_mod, ind = ind,  phen = phen, keep_indivs = ind_keep) %>% 
      do.call("rbind",.) 
    ph.scox = list(scox) 
    names(ph.scox) = disorder
    saveRDS(ph.scox, res_files[ic])
  }, future.seed = T)
}, future.seed = T)
