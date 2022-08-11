library(dplyr)
library(stringr)
library(bigsnpr)
library(data.table)

which_snps_are_sugg = function(gwas, method, phenotype) {
  lapply(seq_along(method[[1]]), function(n) {
    cur_method = method[[1]][n]
    if (str_detect(cur_method, "SPACox")) {
      gwas[[n]] %>% 
        mutate(sig = p.value.spa < 5e-6) %>% 
        select(sig) %>%
        rename(!!as.symbol(paste0("sig_", cur_method)) := sig)
      
    } else {
      tibble(!!as.symbol(paste0("sig_", cur_method)) := 10^(predict(gwas[[n]])) < 5e-6)
    }
    
  } ) %>% 
    do.call("bind_cols", .) %>% 
    rowSums() > 0
}

ipsych_gwas = list.files("./results", "rds", full.names = T) %>% 
  str_subset("ipsych_gwas") %>% 
  str_subset("thirdDegree") %>% 
  lapply(., function(x) tibble(gwas = readRDS(x),
                               phenotype = names(gwas),
                               noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
                               method = ifelse(str_detect(phenotype, "status"), "CaseControl", "ADuLT")) %>% 
           mutate(phenotype = str_replace(phenotype, "status_", ""))) %>% 
  do.call("bind_rows", .)

spacox_gwas = list.files("./results", "rds", full.names = T) %>% 
  str_subset("ipsych_SPACox") %>% 
  lapply(., function(x) {
    tibble(gwas = list(readRDS(x)),
           phenotype = str_match(x, "adhd|asd|dep|scz")[,1],
           noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
           method = "SPACox")
  }) %>% do.call("bind_rows", .)

res = bind_rows(ipsych_gwas,spacox_gwas) %>% 
  relocate(method, noAge, phenotype, gwas) %>% 
  group_by(noAge, phenotype) %>% 
  summarise(comb = list(gwas),
            method_ = list(method), .groups = "drop") %>% 
  mutate(sig = purrr::pmap(.l = list(comb = comb, 
                                     method = method_,
                                     phenotype = phenotype),
                           .f = ~ which_snps_are_sugg(..1, ..2, ..3))) %>% 
  select(-comb, -method_)

saveRDS(res, "./LDClumpedSNPs/ipsych-atleast-sugg-for-one.rds")

sapply(res$sig, sum)

obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS
for (i in 1:4) {
  print(i)
  #which snps do we keep?
  inds = which(res$sig[[i]])
  #extract snps
  cur_CHR = CHR[inds]
  cur_POS = POS[inds]
  
  #extract indivs
  exclude = fread(paste0("./covariates/ipsych-", res$phenotype[i] ,"-excluded-thirdDegree-indivs.txt")) %>% pull(FID)
  KEEP = !(obj.bigsnp$fam$family.ID %in% exclude)
  
  snp_subset(x = obj.bigsnp,
             ind.row = which(KEEP),
             ind.col = inds,
             backingfile = paste0("./LDClumpedSNPs/ipsych-", res$phenotype[i], "-thirdDegree-clumping"))
}




library(future.batchtools)
NCORES = 4
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * NCORES, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))


for (i in 1:nrow(res)) {
  print(i)
#future.apply::future_lapply(X = 1:nrow(res), function(i) {
  
  obj.bigsnp <- snp_attach(paste0("./LDClumpedSNPs/ipsych-", res$phenotype[i], "-thirdDegree-clumping.rds"))

  #clump snps
  clumped_snps = snp_clumping(G = obj.bigsnp$genotypes,
                              infos.chr = obj.bigsnp$map$CHR,
                              thr.r2 = 0.5,
                              infos.pos = obj.bigsnp$map$POS,
                              ncores = NCORES)
  
  #save res
  saveRDS(clumped_snps,
         paste0("./LDClumpedSNPs/ipsych-", res$phenotype[i] , ifelse(res$noAge[i] == "yes","-noAge", ""),"-thirdDegree-clumped-indices.rds"))
}#, future.seed = T)


info_tbl = tibble(indicies  = list.files("./LDClumpedSNPs", "indices", full.names = T),
                  genotypes = rep(list.files("./LDClumpedSNPs", "clumping", full.names = T) %>%
                                    str_subset("rds"), each = 2))

for (i in 1:nrow(info_tbl)) {
  print(i)
  tmp = readRDS(info_tbl$indicies[i])
  obj.bigsnp = snp_attach(info_tbl$genotypes[i])
  saveRDS(slice(obj.bigsnp$map, tmp),
          str_replace(info_tbl$genotypes[i], "clumping", "extracted-clumped-snps"))
}
