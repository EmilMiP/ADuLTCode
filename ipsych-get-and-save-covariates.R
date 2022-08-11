library(dplyr)
library(bigreadr)
library(purrr)
library(lubridate)
library(bigsnpr)
library(ggplot2)
library(xgboost)

#loading data:
PC <- readRDS("./genData/PC.rds")
#getting ids for ipsych individuals in combined data
obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")

#individuals in the ipsych design:
ipsych_info = as_tibble(fread2("~/Register/2019_04/csv/ipsych2015design.csv")) #ipsych individuals

# gotta match on family.ID
mean(obj.bigsnp$fam[,1] %in% ipsych_info$pid) #match on family.ID

#getting sex and origin:
table(sex <- obj.bigsnp$fam$sex - 1, exclude = NULL)
sex[sex == -1] <- NA
table(is_2012 <- obj.bigsnp$fam$is_2012, exclude = NULL)

#merging to get age at end of ipsych duration, i.e. 2016
covar = ipsych_info %>%
  select(pid, fdato, gender, origin2015) %>%
  mutate(new_date = dmy(fdato),
         age = interval( new_date,dmy("01/01/2016")) / dyears(1)) %>% 
  select(-new_date) %>% 
  left_join(obj.bigsnp$fam %>% as_tibble() %>% 
              select(family.ID, sample.ID),
            .,
            by = c("family.ID" = "pid")) %>% 
  rename(FID = family.ID, 
         IID = sample.ID) %>% 
  bind_cols("sex" = sex, "is_2012" = is_2012, PC)

colnames(covar)[-(1:8)] = paste0("PC", 1:20)

fwrite(covar, 
       "./covariates/ipsych-covariates.txt",
       sep = " ",
       quote = F,
       na = "NA")
