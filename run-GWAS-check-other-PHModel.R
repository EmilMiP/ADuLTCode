# library(survival)
# library(dplyr)
# library(bigsnpr)
# library(data.table)
# 
# 
# obj.snp = snp_attach("./simulatedData/genotypes-SPACox-N1000k-M20k-normal-C1000-v4.rds")
# snpInfo = fread("./simulatedData/genotypes-SPACox-N1000k-M20k-normal-C1000-v4.snpinfo") %>% as_tibble()
# phen = fread("./simulatedData/genotypes-SPACox-N1000k-M20k-normal-C1000-prev5-v4.phen") %>% as_tibble()
# keep = fread("./simulatedData/genotypes-SPACox-N1000k-M20k-normal-C1000-prev5-v4.keepers") %>% as_tibble()
# ind_causal = which(snpInfo$beta != 0)
# ind_some_null = sample(setdiff(1:nrow(snpInfo), ind_causal), size = 100) %>% sort
# 
# first_causal = obj.snp$genotypes[keep %>% pull(ids),ind_causal[1:10]]
# first_null   = obj.snp$genotypes[keep %>% pull(ids),11:20] #no causal snps here
# colnames(first_causal) = paste0("snp", 1:10)
# colnames(first_null)   = paste0("null", 1:10)
# phen2 = bind_cols(phen[keep %>% pull(ids),],as_tibble(first_causal), as_tibble(first_null))
# 
# ph_causal = coxph(Surv(event, status) ~ snp1, phen2)
# summary(ph_causal)
# ph_null   = coxph(Surv(event, status) ~ null1, phen2)
# summary(ph_null)
# 
# # tmp = summary(ph_null)
# # tmp$coefficients[5]



#   -----------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(bigsnpr)
library(survival)
library(future.batchtools)


prevs = c(0.0101, 0.025, 0.05, 0.1, 0.2, 0.5)
#constructing prevalence name to keep file names ordered
prev_text = paste0("prev",round(prevs * 100,1))
prev_text[1] = "prev01"
prev_text[2] = "prev025"


all_files = list.files("./simulatedData", "SPA", full.names = T) %>% 
  str_subset("N1000k")


phenotypes = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "prev") %>% 
    str_subset(., "true")
  tibble(full_names = ph,
         v = v) %>% 
    mutate(gen_beta = as.vector(str_match(full_names, "unif|normal")),
           prev = as.vector(str_match(full_names, paste(prev_text, "-", collapse = "|", sep = ""))),
           C = as.vector(str_match(full_names, "C250|C1000"))) %>% 
    group_by(v, gen_beta, C, prev) %>% 
    summarise(true = full_names, .groups = "drop") %>% 
    ungroup
}) %>% 
  do.call("rbind", .) 

files = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "true|phen", negate = T)
  tibble(full_names = ph,
         v = v) %>% 
    mutate(gen_beta = as.vector(str_match(full_names, "unif|normal")),
           prev = as.vector(str_match(full_names, paste(prev_text, "-", collapse = "|", sep = ""))),
           C = as.vector(str_match(full_names, "C250|C1000"))) %>% 
    group_by(v, gen_beta, C) %>% 
    summarise(grps = list(full_names), .groups = "drop") %>% 
    ungroup
}) %>% 
  do.call("rbind", .) %>% 
  left_join(., phenotypes) %>% 
  mutate(grps = purrr::map2(grps, true, ~ c(.x, .y))) %>% 
  select(-true) %>% 
  mutate(grps = purrr::map2(grps, prev, ~ .x[str_detect(.x, "keepers", negate = T) | str_detect(.x, .y)]),
         downsampling = purrr::map(grps, ~ unique(str_detect(.x, "keepers")))) %>% 
  tidyr::unnest(cols = downsampling) %>% 
  relocate(v, gen_beta, C, prev, downsampling) %>% 
  mutate(prev = str_replace(prev, "-", ""))



NCORES = 8

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

future.apply::future_lapply(1:nrow(files), function(ctr) {
  #getting identifiers
  v = files$v[[ctr]]
  gen_beta = files$gen_beta[[ctr]]
  prev = files$prev[[ctr]]
  cur_C = files$C[[ctr]]
  downsampling = files$downsampling[[ctr]]
  
  #loading data
  bigsnp_obj   = snp_attach(str_subset(files$grps[[ctr]], "rds"))
  G = bigsnp_obj$genotypes
  cur_snp_info = fread(str_subset(files$grps[[ctr]], "snpinfo")) %>% as_tibble
  cur_true     = fread(str_subset(files$grps[[ctr]], "true")) %>%
    as_tibble %>% 
    mutate(event = pmin(onset, censor))
  
  
  #is downsampling a posibility ?
  if (downsampling) {
    keepers      = fread(str_subset(files$grps[[ctr]], "keepers")) %>% as_tibble
    KEEP = which(1:nrow(G) %in% keepers[[1]])
    phen = fread(str_subset(files$grps[[ctr]], "true")) %>% 
      as_tibble() %>% 
      slice(keepers[[1]]) %>% 
      mutate(event = pmin(onset, censor))
  } else {
    KEEP = 1:nrow(G)
    phen = fread(str_subset(files$grps[[ctr]], "true")) %>% 
      as_tibble %>% 
      mutate(event = pmin(onset, censor))
  }
  
  # analysis will only be run on a subset of SNPs
  ind_causal = which(cur_snp_info$beta != 0)
  ind_some_null = sample(setdiff(1:nrow(cur_snp_info), ind_causal), size = length(ind_causal)) %>% sort
  analysed_snps = tibble(
    index = c(ind_causal, ind_some_null),
    causal = c(rep("causal", length(ind_causal)), rep("null", length(ind_some_null)))
  ) %>% arrange(index)
  
  res_files = paste0("./gwas-results/v",v, "_PHModel_N1000k_M20k_", 
                     gen_beta, "_", 
                     cur_C, "_", 
                     prev, 
                     if (downsampling) "_DS", ".rds")
  intervals = bigparallelr::split_len(nrow(analysed_snps), nb_split = 20)
  mean(already_done <- file.exists(res_files))
  
  all_gwas_parts = lapply(1:nrow(intervals), function(ic) {
    #indecies are determined
    ind = seq(intervals[ic, "lower"], intervals[ic, "upper"])
    
    big_apply(X = G, a.FUN = function(X, ind, analysed_snps, phen, ind_keep) {
      library(dplyr)
      library(survival)
  
      tmp_X = X[ind_keep, analysed_snps$index[ind]]
      colnames(tmp_X) = paste0("snp", analysed_snps$index[ind])
      lapply(as_tibble(tmp_X), function(cur_snp) {
        phen2 = bind_cols(phen, cur_snp = cur_snp)
        PHModel = coxph(Surv(event, status) ~ cur_snp, phen2)
        summary(PHModel)$coefficients[5]
      }) %>% unlist()
      
    }, ncores = NCORES, analysed_snps = analysed_snps, phen = phen, ind = ind, ind_keep = KEEP) %>% 
      do.call("c", .)

  }) %>% do.call("c", .)
  analysed_snps = analysed_snps %>% mutate(pvalue = all_gwas_parts)
  
  saveRDS(analysed_snps, res_files)
  NULL
  
  
}, future.seed = T)





