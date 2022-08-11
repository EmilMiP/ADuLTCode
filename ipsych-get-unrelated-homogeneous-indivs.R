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

table(sex <- obj.bigsnp$fam$sex - 1, exclude = NULL)
sex[sex == -1] <- NA
table(is_2012 <- obj.bigsnp$fam$is_2012, exclude = NULL)

# PC ----------------------------------------------------------------------

PC <- readRDS("./genData/PC.rds")
hist(log_dist <- log(bigutilsr::dist_ogk(PC)))
is_homogeneous <- log_dist < 4.5

str(COVAR <- cbind(PC, sex, is_2012))
mean(has_no_na <- complete.cases(COVAR)) # 99.95%

# Loading Relatedness -----------------------------------------------------

rel <- readRDS("./genData/rel.rds")

is_rel <- obj.bigsnp$fam$sample.ID %in% subset(rel, KINSHIP > 2^-3.5)$IID2
mean(is_rel) # 11.0%

mean(KEEP <- !is_rel & has_no_na & is_homogeneous) # 80.2%
ind_keep <- which(KEEP)

ph = obj.bigsnp$fam[ind_keep,1:2]
colnames(ph) = c("FID", "IID")
fwrite(ph, "./covariates/ipsych-unrelated-homogeneous-indivs.txt", sep = " ", quote = F, na = "NA")

is_rel <- obj.bigsnp$fam$sample.ID %in% subset(rel, KINSHIP > 2^-2.5)$IID2
mean(is_rel) # 11.0%

mean(KEEP <- !is_rel & has_no_na & is_homogeneous) # 80.2%
ind_keep <- which(KEEP)

ph = obj.bigsnp$fam[ind_keep,1:2]
colnames(ph) = c("FID", "IID")
fwrite(ph, "./covariates/ipsych-unrelated-homogeneous-thirdDegree-indivs.txt", sep = " ", quote = F, na = "NA")

tmp = fread("./covariates/ipsych-unrelated-homogeneous-thirdDegree-indivs.txt") %>% as_tibble()
# Loading Phenotypes ------------------------------------------------------

info = tibble(
  disorder_name = c("adhd", "asd", "dep", "scz"),
  phenotype     = paste0("./phenotypes/", disorder_name, ".phen")
)

pheno = lapply(1:nrow(info), function(ctr) {
  disorder = info$disorder_name[ctr]
  fread(info$phenotype[ctr]) %>% 
    as_tibble() %>% 
    select(FID, IID, status, post_gen_no_fam) %>% 
    rename(!!as.symbol(disorder) := post_gen_no_fam,
           !!as.symbol(paste0("status_", disorder)) := status) 
}) %>% purrr::reduce(., left_join) %>% 
#  rename("family.ID" = "FID") %>% 
  distinct(FID, IID, .keep_all = T) %>% 
  left_join(obj.bigsnp$fam[,1:2], ., by = c("family.ID" = "FID", "sample.ID" = "IID")) %>% 
  as_tibble %>% 
  rename(FID = "family.ID",
         IID = "sample.ID")



controls = as_tibble(fread("~/Register/2019_04/csv/ipsych2015design.csv")) %>% 
  select(pid, kontrol2015I) %>% 
  left_join(pheno, ., by = c("FID" = "pid")) %>% 
  select(FID, IID, contains("status"), kontrol2015I) %>% 
  mutate(across(.cols = starts_with("status_"), .fns = ~ .x + kontrol2015I > 0)) %>% 
  select(contains("status_")) %>% 
  rename_with(~ str_replace(.x, "status_", ""))


# get list of individuals to EXCLUDE --------------------------------------

new_pheno = select(pheno, starts_with("status_"))
excluded = lapply(seq_along(new_pheno), function(n) {
  #get the raw phenotype vector
  y = new_pheno[[ n ]]
  #get the name of the phenotype
  cur_name = names(new_pheno)[n] %>% 
    str_replace("status_", "")
  #find the controls and cases only
  ctrl_inds = controls[[cur_name]]
  #perform GWAS on homogen etc. indivs.
  !(KEEP & ctrl_inds)
}) 

names(excluded) = names(controls)
excluded = as_tibble(excluded)

#looks like they are largely inverted, with the difference being the related or non-homogen indivs
bind_cols(excluded %>% select(adhd) %>% rename(exclude = adhd),
          controls %>% select(adhd) %>% rename(control = adhd)) %>% 
  count(exclude, control)


fam = as_tibble(obj.bigsnp$fam) %>% 
  rename(FID = family.ID,
         IID = sample.ID)

for (i in 1:ncol(excluded)) {
  cur_name = colnames(excluded)[i]
  cur_excluded = fam %>% 
    filter(excluded[[i]] | is.na(excluded[[i]])) %>% 
    select(FID, IID)
  fwrite(cur_excluded,
         paste0("./covariates/ipsych-", cur_name, "-excluded-thirdDegree-indivs.txt"),
         sep = " ",
         quote = F,
         na = "NA")
}

setdiff(select(pheno, FID, IID), cur_excluded) %>% 
  left_join(pheno) %>% 
  summarise(prev = mean(status_scz),
            n = n(),
            cases = sum(status_scz))
