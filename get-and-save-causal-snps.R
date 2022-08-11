library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(bigsnpr)

#SPACox, linear, and logistic regression results:
all_files = list.files("./gwas-results", full.names = T, "LTM") %>% 
  str_subset(., negate = T, pattern = "noFH|false-positive|part")  

tmp = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, 
                                    paste0("prev", rep(100*c(0.05, 0.1, 0.2, 0.5), each = 2), 
                                           c("\\.rds", "_DS.rds"), 
                                           collapse = "|"))),
         prev = str_replace_all(prev, "\\.rds|_DS\\.rds", ""),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup


res = lapply(1:nrow(tmp), function(ctr) { #nrow(tmp)
  C = tmp$C[[ctr]]
  prev_txt = tmp$prev[[ctr]]
  ds = tmp$downsampling[[ctr]]
  
  tmp3 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_LTM_", C,"_", prev_txt, if(ds == "yes") "_DS" ,".rds"))
    #getting info per method, v, maf and generative model
    sum_info = lapply(seq_along(gwas), function(n) {
      if (n < 3) {
        tmp_res = snpinfo %>% select(beta) %>%
          mutate(p_val = predict(gwas[[n]], log10 = F),
                 p_val = ifelse(p_val == 0, 5e-300, p_val),
                 chisq = qchisq(p_val, df = 1, lower.tail = F),
                 method = names(gwas)[n],
                 z_val = gwas[[n]]$score) %>% 
          filter(beta != 0)
      } else {
        tmp_res = snpinfo %>% 
          select(beta) %>%
          bind_cols(gwas[[n]]) %>% 
          select(beta, p.value.spa, z) %>% 
          rename(p_val = p.value.spa,
                 z_val = z) %>% 
          mutate(p_val = ifelse(p_val == 0, 5e-300, p_val),
                 chisq = qchisq(p_val, df = 1, lower.tail = F),
                 method = names(gwas)[n]) %>%
         filter(beta != 0)
      }
      
    }) %>% 
      do.call("bind_rows", .) %>% 
      mutate(v = v)
    
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(C = C,
           prev = prev_txt,
           downsampling = ds)
}) %>% 
  do.call("bind_rows", .) %>% 
  relocate(v, C, prev, downsampling, method)


#liability based method results:
all_files2 = list.files("./gwas-results", full.names = T, "LTM") %>% 
  str_subset(., negate = T, pattern = "part") %>% 
  str_subset(., "noFH") 

tmp2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, paste0("prev", 100*c(0.05, 0.1, 0.2, 0.5), "_", collapse = "|"))),
         prev = str_replace(prev, "_", ""),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup



res2 = lapply(1:nrow(tmp2), function(ctr) { #nrow(tmp)
  C        = tmp2$C[[ctr]]
  prev_txt = tmp2$prev[[ctr]]
  ds       = tmp2$downsampling[[ctr]]
  
  tmp4 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_LTM_", C,"_", prev_txt, "_noFH", if (ds == "yes") "_DS" ,".rds"))
    #getting info per method, v, maf and generative model
    snpinfo %>% select(beta) %>%
      mutate(p_val = predict(gwas, log10 = F),
             p_val = ifelse(p_val == 0, 5e-300, p_val),
             chisq = qchisq(p_val, df = 1, lower.tail = F),
             v = v) %>% 
        filter(beta != 0) 
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(method = "ADuLT",
           C = C,
           prev = prev_txt,
           downsampling = ds)
  
}) %>% 
  do.call("bind_rows", .) %>%   
  relocate(v, C, prev, downsampling, method)

#save the summary information for all methods
saveRDS(bind_rows(res,res2),
        "./results/ltm_causal_snps_all.rds")



# SPACox based results ----------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(bigsnpr)


#SPACox, linear, and logistic regression results:
all_files = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part|noFH") %>% 
  str_subset(., "ER")

tmp = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         event_rate = as.vector(str_match(ph, "ER[1-8]0|ER15|ER1")),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, event_rate,downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup


res = lapply(1:nrow(tmp), function(ctr) { #nrow(tmp)
  C       = tmp$C[[ctr]]
  cur_gen = tmp$gen_beta[[ctr]]
  cur_event_rate = tmp$event_rate[[ctr]]
  ds = tmp$downsampling[[ctr]]
  
  tmp3 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_", cur_gen,"_", C,  "_", cur_event_rate, if(ds == "yes") "_DS", ".rds"))
    #getting info per method, v, maf and generative model
    sum_info = lapply(seq_along(gwas), function(n) {
      if (n < 3) {
        tmp_res = snpinfo %>% select(beta) %>%
          mutate(est_beta = gwas[[n]]$estim,
                 p_val = predict(gwas[[n]], log10 = F),
                 chisq = qchisq(p_val, df = 1, lower.tail = F),
                 method = names(gwas)[n],
                 z_val = gwas[[n]]$score) %>% 
            filter(beta != 0) 
      } else {
        tmp_res = snpinfo %>% 
          select(beta) %>%
          bind_cols(gwas[[n]]) %>% 
          select(beta, p.value.spa, z) %>% 
          rename(p_val = p.value.spa,
                 z_val = z) %>% 
          mutate(chisq = qchisq(p_val, df = 1, lower.tail = F),
                 method = names(gwas)[n]) %>%
          filter(beta != 0) 
      }
      
    }) %>% 
      do.call("bind_rows", .) %>% 
      mutate(v = v)
    
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(C = C,
           beta_gen = cur_gen,
           event_rate = cur_event_rate,
           downsampling = ds)
}) %>% 
  do.call("bind_rows", .) %>% 
  relocate(v, beta_gen, C, event_rate, method, downsampling)


#liability based method results:
all_files2 = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part") %>% 
  str_subset(., "noFH") %>% 
  str_subset(., "ER")

tmp2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         event_rate = as.vector(str_match(ph, "ER[1-8]0|ER15|ER1")),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, event_rate, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup


res2 = lapply(1:nrow(tmp2), function(ctr) { #nrow(tmp)
  C       = tmp2$C[[ctr]]
  cur_gen = tmp2$gen_beta[[ctr]]
  cur_event_rate = tmp2$event_rate[[ctr]]
  ds = tmp2$downsampling[[ctr]]
  
  tmp4 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_", cur_gen,"_", C,  "_", cur_event_rate, "_noFH", if(ds == "yes") "_DS", ".rds"))
    #getting info per method, v, maf and generative model
    tmp_res = snpinfo %>% select(beta) %>%
      mutate(est_beta = gwas$estim,
             p_val = predict(gwas, log10 = F),
             chisq = qchisq(p_val, df = 1, lower.tail = F),
             v = v,
             z_val = gwas$score) %>% 
        filter(beta != 0)
  }) %>% 
    do.call("bind_rows", .) %>% 
    mutate(method = "ADuLT",
           C = C,
           beta_gen = cur_gen,
           event_rate = cur_event_rate,
           downsampling = ds)
  
}) %>% 
  do.call("bind_rows", .) %>%   
  relocate(v, beta_gen, C, method, event_rate, downsampling)

#save the summary information for all methods
saveRDS(bind_rows(res,res2),
        "./results/cox_causal_snps_all.rds")

