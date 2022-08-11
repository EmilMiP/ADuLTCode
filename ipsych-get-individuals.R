library(dplyr)
library(bigreadr)
library(purrr)
library(lubridate)
library(bigsnpr)
library(ggplot2)
library(xgboost)

#individuals in the ipsych design:
ipsych_info = as_tibble(fread2("~/Register/2019_04/csv/ipsych2015design.csv")) #ipsych individuals
str(ipsych_info)

#link between people in the design and family members:
info = as_tibble(fread2("~/Register/2019_04/csv/stam2016h.csv")) %>%
  mutate_at(c("pid_m", "pid_f"), as.integer)
str(info)  # pids of parents and birth dates

#phenotypes for all ipsych and related individuals.
phenotype_data = as_tibble(fread2("~/Register/2021_01/phenotype2016j.csv")) # slightly updated phenotypes
str(phenotype_data)

#this file contains the pid linked to sampleIDs:
linked_data = as_tibble(fread2(
  "~/NCRR-PRS/faststorage/emil/bedfiles/sampleID_linkedIndivs/linked_indivs.txt"))
str(linked_data)


ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/dosage_ipsych2015.rds")
str(ipsych$fam)


str(cip <- bigreadr::fread2("~/NCRR-PRS/faststorage/ncrr_genetics/prevalence/RESULTS_cip.txt"))



match <- left_join(as_tibble(ipsych$fam),
                   info[c("pid", "pid_m", "pid_f")],
                   by = c("family.ID" = "pid")) %>% 
  rename("pid" = "family.ID")


# #disorder = "adhd2015I" #col name of disorder in ipsych
# disorder_code = "F9202" #col name of disorder in phenotype file
# #adhd: "F9202"
# #depression: "F3300"
# #AUT: "F8101
# #scz: F2000
# code = "adhd" #disorder code in prevalence data
# #DEP: D33
# #adhd: adhd
# #aut: asd2 (largest of the two groups)
# #scz: D20
# 
# disorder_name = "adhd"
# #depression: dep
# #adhd: adhd
# #aut: asd2 (largest of the two groups)
# #schizophrenia: scz
disorder_code_vec = c("F9202", "F3300", "F8101", "F2000")
code_vec          = c("adhd", "D33", "asd2", "D20")
disorder_name_vec = c("adhd", "dep", "asd", "scz")

for (i in seq_along(disorder_code_vec)) {
  disorder_code = disorder_code_vec[i]
  code          = code_vec[i]
  disorder_name = disorder_name_vec[i]
  
  
  cip_adhd <- cip %>%
    as_tibble() %>%
    filter(dx == code) %>%
    select(-dx, -pct)
  
  y <- cip_adhd$cip
  
  X <- bigstatsr::covar_from_df(select(cip_adhd, -cip))
  bst.params <- list(
    max_depth  = 10,
    base_score = 0,
    verbose    = 2,
    nthread    = 4,
    min_child_weight  = 10
  )
  xgb <- xgboost::xgboost(X, y, nrounds = 500, params = bst.params)
  
  all_info <- match %>% 
    select(pid, sample.ID) %>% 
    left_join(info[c("pid", "fdato", "gender")]) %>%
    mutate(fdato = dmy(fdato),
           age = interval(fdato, dmy("01/01/2016")) / years(1)) %>%
    left_join(phenotype_data[c("pid", disorder_code)]) %>%
    mutate(AOO = ifelse(!!as.symbol(disorder_code) %in% c("", NA), NA, interval(fdato, dmy(!!as.symbol(disorder_code))) / years(1))) %>%
    select(-!!as.symbol(disorder_code)) %>%
    mutate(status = !is.na(AOO),
           sexM = (gender == "M") + 0,
           birth_year = year(fdato))
  
  all_info$cip_pred <- select(all_info, age, birth_year, sexM) %>%
    as.matrix() %>%
    predict(xgb, .) %>%
    pmax(1e-3)
  all_info = all_info %>%
    group_by(birth_year, sexM) %>%
    arrange(age) %>%
    mutate(cip_pred = cummax(cip_pred)) %>%
    ungroup() %>%
    mutate(thr = qnorm(cip_pred, lower.tail = FALSE),
           lower = ifelse(is.na(AOO), -Inf, thr),
           upper = ifelse(is.na(AOO), thr,  thr))
  
  saveRDS(all_info, paste0("./phenotypes/", disorder_name, ".rds"))
}  


for (disorder_name in disorder_name_vec) {
  all_info = fread(paste0("./phenotypes/", disorder_name, ".phen")) %>%
    as_tibble() %>% 
    mutate(status = status + 0L)
  fwrite(all_info,library(dplyr)
         library(bigreadr)
         library(purrr)
         library(lubridate)
         library(bigsnpr)
         library(ggplot2)
         library(xgboost)
         
         #individuals in the ipsych design:
         ipsych_info = as_tibble(fread2("~/Register/2019_04/csv/ipsych2015design.csv")) #ipsych individuals
         str(ipsych_info)
         
         #link between people in the design and family members:
         info = as_tibble(fread2("~/Register/2019_04/csv/stam2016h.csv")) %>%
           mutate_at(c("pid_m", "pid_f"), as.integer)
         str(info)  # pids of parents and birth dates
         
         #phenotypes for all ipsych and related individuals.
         phenotype_data = as_tibble(fread2("~/Register/2021_01/phenotype2016j.csv")) # slightly updated phenotypes
         str(phenotype_data)
         
         #this file contains the pid linked to sampleIDs:
         linked_data = as_tibble(fread2(
           "~/NCRR-PRS/faststorage/emil/bedfiles/sampleID_linkedIndivs/linked_indivs.txt"))
         str(linked_data)
         
         
         ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/dosage_ipsych2015.rds")
         str(ipsych$fam)
         
         
         str(cip <- bigreadr::fread2("~/NCRR-PRS/faststorage/ncrr_genetics/prevalence/RESULTS_cip.txt"))
         
         
         
         match <- left_join(as_tibble(ipsych$fam),
                            info[c("pid", "pid_m", "pid_f")],
                            by = c("family.ID" = "pid")) %>% 
           rename("pid" = "family.ID")
         
         
         # #disorder = "adhd2015I" #col name of disorder in ipsych
         # disorder_code = "F9202" #col name of disorder in phenotype file
         # #adhd: "F9202"
         # #depression: "F3300"
         # #AUT: "F8101
         # #scz: F2000
         # code = "adhd" #disorder code in prevalence data
         # #DEP: D33
         # #adhd: adhd
         # #aut: asd2 (largest of the two groups)
         # #scz: D20
         # 
         # disorder_name = "adhd"
         # #depression: dep
         # #adhd: adhd
         # #aut: asd2 (largest of the two groups)
         # #schizophrenia: scz
         disorder_code_vec = c("F9202", "F3300", "F8101", "F2000")
         code_vec          = c("adhd", "D33", "asd2", "D20")
         disorder_name_vec = c("adhd", "dep", "asd", "scz")
         
         for (i in seq_along(disorder_code_vec)) {
           disorder_code = disorder_code_vec[i]
           code          = code_vec[i]
           disorder_name = disorder_name_vec[i]
           
           
           cip_adhd <- cip %>%
             as_tibble() %>%
             filter(dx == code) %>%
             select(-dx, -pct)
           
           y <- cip_adhd$cip
           
           X <- bigstatsr::covar_from_df(select(cip_adhd, -cip))
           bst.params <- list(
             max_depth  = 10,
             base_score = 0,
             verbose    = 2,
             nthread    = 4,
             min_child_weight  = 10
           )
           xgb <- xgboost::xgboost(X, y, nrounds = 500, params = bst.params)
           
           all_info <- match %>% 
             select(pid, sample.ID) %>% 
             left_join(info[c("pid", "fdato", "gender")]) %>%
             mutate(fdato = dmy(fdato),
                    age = interval(fdato, dmy("01/01/2016")) / years(1)) %>%
             left_join(phenotype_data[c("pid", disorder_code)]) %>%
             mutate(AOO = ifelse(!!as.symbol(disorder_code) %in% c("", NA), NA, interval(fdato, dmy(!!as.symbol(disorder_code))) / years(1))) %>%
             select(-!!as.symbol(disorder_code)) %>%
             mutate(status = !is.na(AOO),
                    sexM = (gender == "M") + 0,
                    birth_year = year(fdato))
           
           all_info$cip_pred <- select(all_info, age, birth_year, sexM) %>%
             as.matrix() %>%
             predict(xgb, .) %>%
             pmax(1e-3)
           all_info = all_info %>%
             group_by(birth_year, sexM) %>%
             arrange(age) %>%
             mutate(cip_pred = cummax(cip_pred)) %>%
             ungroup() %>%
             mutate(thr = qnorm(cip_pred, lower.tail = FALSE),
                    lower = ifelse(is.na(AOO), -Inf, thr),
                    upper = ifelse(is.na(AOO), thr,  thr))
           
           saveRDS(all_info, paste0("./phenotypes/", disorder_name, ".rds"))
         }  
         
         
         for (disorder_name in disorder_name_vec) {
           all_info = fread(paste0("./phenotypes/", disorder_name, ".phen")) %>%
             as_tibble() %>% 
             mutate(status = status + 0L)
           fwrite(all_info,
                  paste0("./phenotypes/", disorder_name, ".phen"),
                  sep = " ",
                  quote = F,
                  na = "NA")
         }
         paste0("./phenotypes/", disorder_name, ".phen"),
         sep = " ",
         quote = F,
         na = "NA")
}
