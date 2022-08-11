library(bigsnpr)
library(dplyr)
library(stringr)
library(data.table)

#G <- obj.bigsnp$genotypes
map = fread("../results/export_map.txt") %>% as_tibble()

CHR <- map$CHR
POS <- map$POS
SNP = map$SNP


ipsych_gwas = list.files("../results", "rds", full.names = T) %>% 
  str_subset("ipsych_gwas") %>% 
  str_subset("thirdDegree") %>% 
  lapply(., function(x) tibble(gwas = readRDS(x),
                               phenotype = names(gwas),
                               noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
                               method = ifelse(str_detect(phenotype, "status"), "CaseControl", "ADuLT")) %>% 
           mutate(phenotype = str_replace(phenotype, "status_", ""))) %>% 
  do.call("bind_rows", .)


for(i in 1:nrow(ipsych_gwas)) {
  print(i)
  cur_path = paste0("../results/ipsych-", ipsych_gwas$phenotype[i], if(ipsych_gwas$noAge[i] == "yes") "-noAge", "-thirdDegree-", ipsych_gwas$method[i], ".sumstats")
  fwrite(bind_cols(select(map, CHR, POS, SNP), ipsych_gwas$gwas[[i]], p_value = 10^(predict(ipsych_gwas$gwas[[i]]))),
        cur_path)
}

tmp %>% as_tibble() %>% 
  filter(info_2012_homo < 0.5 | info_2015_homo < 0.5)

tmp %>% as_tibble() %>% 
  filter(freq_2012_homo < 0.01| freq_2012_homo > 0.99)

tmp %>% as_tibble() %>% 
  filter(freq_2015_homo < 0.01| freq_2015_homo > 0.99)
