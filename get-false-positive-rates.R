library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(bigsnpr)


#SPACox, linear, and logistic regression results:
#SPACox, linear, and logistic regression results:
all_files = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part|noFH") %>% 
  str_subset("N1000k_M20k")

#organizing files
files = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         prev = as.vector(str_match(ph, paste0(rep(prev_text_vec,each = 2), c("\\.rds", "_DS"), collapse = "|"))),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = str_replace(prev, "\\.rds|_DS", ""))


phenotypes = c("linear", "SPACox")

fp = lapply(1:nrow(files), function(ctr) { #nrow(tmp)
  C = files$C[[ctr]]
  cur_gen = files$gen_beta[[ctr]]
  cur_prev = files$prev[[ctr]]
  ds = files$downsampling[[ctr]]
  
  sum_info2 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-N1000k-M20k-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_N1000k_M20k_",
                          cur_gen,"_", 
                          C,  "_", 
                          cur_prev, 
                          if (ds == "yes") "_DS", ".rds"))
    
    sum_info = lapply(seq_along(gwas), function(ctr2) {
      
      if (ctr2 < 2) {
        ph = gwas[[ctr2]] %>% 
          mutate(p_val = predict(gwas[[ctr2]], log10 = F)) %>% 
          bind_cols(snpinfo %>% select(beta)) 
        
        lapply(5*10^(-(2:8)), function(alpha) {
          ph %>% 
            filter(beta == 0) %>%
            summarise(FP = sum(p_val < alpha),
                      FP_prop   = FP / n(),
                      FP_sem    = sqrt(FP_prop * (1 - FP_prop) / n()),
                      alpha_lvl = alpha,
                      n = n())
        }) %>% 
          do.call("rbind", .) %>%
          mutate(method = phenotypes[ctr2])
        
      } else {
        ph = gwas[[ctr2]] %>% 
          bind_cols(snpinfo %>% select(beta)) 
        
        lapply(5*10^(-(2:8)), function(alpha) {
          ph %>% 
            filter(beta == 0) %>%
            summarise(FP = sum(p.value.spa < alpha),
                      FP_prop   = FP / n(),
                      FP_sem    = sqrt(FP_prop * (1 - FP_prop) / n()),
                      alpha_lvl = alpha,
                      n = n())
        }) %>% 
          do.call("bind_rows", .) %>%
          mutate(method = phenotypes[ctr2])
      }
      
    }) %>% 
      do.call("bind_rows", .) %>%
      as_tibble() %>% 
      mutate(v = v)
    
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(C = C,
           beta_gen = cur_gen,
           prev = cur_prev,
           downsampling = ds)
  
  
}) %>% 
  do.call("bind_rows", .) %>% 
  relocate(v, beta_gen, C, prev, downsampling, method)





# ADuLT -------------------------------------------------------------------


#liability based method results:
all_files2 = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part") %>% 
  str_subset(., "noFH") %>%
  str_subset("N1000k_M20k")

files2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         prev = as.vector(str_match(ph, paste(prev_text_vec, "_", collapse = "|", sep = ""))),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = str_replace(prev, "_", ""))

phenotypes2 = c("ADuLT")

fp2 = lapply(1:nrow(files2), function(ctr) { #nrow(tmp)
  C       = files2$C[[ctr]]
  cur_gen = files2$gen_beta[[ctr]]
  cur_prev = files2$prev[[ctr]]
  ds = files2$downsampling[[ctr]]
  
  sum_info2 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-N1000k-M20k-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_N1000k_M20k_", 
                          cur_gen,"_", 
                          C, "_", 
                          cur_prev, "_noFH", 
                          if (ds == "yes") "_DS", ".rds"))

    ph = gwas %>%
      mutate(p_val = predict(gwas, log10 = F)) %>% 
      bind_cols(snpinfo %>% select(beta))
    lapply(5*10^(-(2:8)), function(alpha) {
      ph %>%
        filter(beta == 0) %>%
        summarise(FP = sum(p_val < alpha),
                  FP_prop   = FP / n(),
                  FP_sem    = sqrt(FP_prop * (1 - FP_prop) / n()),
                  alpha_lvl = alpha,
                  n = n())
    }) %>% 
      do.call("bind_rows", .) %>%
      mutate(method = phenotypes2,
             v = v)
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(method = "ADuLT",
           C = C,
           beta_gen = cur_gen,
           prev = cur_prev, 
           downsampling = ds)
  
}) %>% 
  do.call("bind_rows", .) %>% 
  as_tibble() %>% 
  relocate(v, beta_gen, C, method, prev, downsampling)



# save data ---------------------------------------------------------------

saveRDS(bind_rows(fp, fp2),
        "./results/false-positive-data-SPACox.rds")
