library(bigsnpr)
library(data.table)
library(dplyr)
library(stringr)
library(future.batchtools)

# names(readRDS("./results/ipsych_gwas_thirdDegree.rds"))
# outputs:
#"status_adhd" "adhd"        "status_asd"  "asd"         "status_dep"  "dep"         "status_scz"  "scz" 

phenotypes = c("status_adhd", "adhd", "status_asd", "asd", "status_dep", "dep", "status_scz", "scz")


ipsych_gwas = list.files("./results", "rds", full.names = T) %>% 
  str_subset("ipsych_gwas") %>% 
  str_subset("thirdDegree") %>% 
  lapply(., function(x) tibble(path = x,
                               phenotype = phenotypes,
                               noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
                               method = ifelse(str_detect(phenotype, "status"), "CaseControl", "ADuLT")) 
         # %>% 
         #   mutate(phenotype = str_replace(phenotype, "status_", ""))
         ) %>% 
  do.call("bind_rows", .)

spacox_gwas = list.files("./results", "rds", full.names = T) %>% 
  str_subset("ipsych_SPACox") %>% 
  lapply(., function(x) {
    tibble(path = x,
           phenotype = str_match(x, "adhd|asd|dep|scz")[,1],
           noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
           method = "SPACox")
  }) %>% do.call("bind_rows", .)

res = bind_rows(ipsych_gwas,spacox_gwas) %>% 
  relocate(method, noAge, phenotype, path) 

#no clear evidence of large differences in p value between p.value.spa and p.value.norm
tmp = readRDS(res$path[18])
tmp %>% mutate(log_diff = abs(-log10(p.value.spa) + log10(p.value.norm))) %>% 
  filter(log_diff > .001) %>% arrange(desc(log_diff))



NCORES <- 16
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))

future.apply::future_lapply(X = 1:nrow(res), function(i) {

  obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")
  
  cur_phen = str_replace(res$phenotype[i], "status_", "")
  #extract indivs
  exclude = fread(paste0("./covariates/ipsych-", cur_phen ,"-excluded-thirdDegree-indivs.txt")) %>% pull(FID)
  KEEP = which(!(obj.bigsnp$fam$family.ID %in% exclude))

  if (res$method[i] != "SPACox") {
    scores = readRDS(res$path[i])[[res$phenotype[i]]]$score #true
  } else {
    scores = readRDS(res$path[i])$z #false
  }
                  
  
#clump snps
  clumped_snps = snp_clumping(G = obj.bigsnp$genotypes,
                              infos.chr = obj.bigsnp$map$CHR,
                              thr.r2 = 0.5,
                              ind.row = KEEP,
                              S = abs(scores),
                              infos.pos = obj.bigsnp$map$POS,
                              ncores = NCORES)
  saveRDS(clumped_snps,
          paste0("./LDClumpedSNPs/singleTrait/ipsych-", cur_phen, "-",
                 res$method[i],  
                 if(res$noAge[i] == "yes") "-noAge" , "-LDclumped-snps.rds"))

}, future.seed = T)


tmp  = readRDS("./LDClumpedSNPs/singleTrait/ipsych-adhd-ADuLT-LDclumped-snps.rds")
tmp2  = readRDS("./LDClumpedSNPs/singleTrait/ipsych-adhd-CaseControl-LDclumped-snps.rds")

gwas = readRDS("./results/ipsych_gwas_thirdDegree.rds")

gwas$status_adhd %>% 
  mutate(pval = -predict(gwas$status_adhd)) %>% 
  slice(tmp2) %>% 
  filter(pval > -log10(5e-8))
