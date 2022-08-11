library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(bigsnpr)

prevs = c(0.0101, 0.025, 0.05, 0.1, 0.2, 0.5)
#constructing prevalence name to keep file names ordered
prev_text_vec = paste0("prev",round(prevs * 100,1))
prev_text_vec[1] = "prev01"
prev_text_vec[2] = "prev025"

#SPACox, linear, and logistic regression results:
all_files = list.files("./gwas-results", full.names = T, "LTM") %>% 
  str_subset("N1000k_M20k") %>% 
  str_subset(negate = T, pattern = "part") %>% 
  str_subset("noFH", negate = T)

tmp = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, 
                                    paste0(rep(prev_text_vec, each = 2), c("\\.rds", "_DS.rds"), 
                                           collapse = "|", sep = ""))),
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
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_LTM_N1000k_M20k_normal_", C,"_", prev_txt, if(ds == "yes") "_DS" ,".rds"))
    #getting info per method, v, maf and generative model
    sum_info = lapply(seq_along(gwas), function(n) {
      if (n < 2) {
        tmp_res = snpinfo %>% select(beta) %>%
          mutate(p_val = predict(gwas[[n]], log10 = F),
                 p_val = ifelse(p_val == 0, 5e-300, p_val),
                 chisq = qchisq(p_val, df = 1, lower.tail = F))
        #summarising data
        bind_cols(
          tmp_res %>% 
            filter(beta == 0) %>%
            summarise(FP = sum(p_val < 5e-8),
                      mean_null_chisq = mean(chisq)),
          tmp_res %>% 
            filter(beta != 0) %>%
            summarise(causal = sum(p_val < 5e-8),
                      mean_causal_chisq = mean(chisq)))
      } else {
        tmp_res = snpinfo %>% 
          select(beta) %>%
          bind_cols(gwas[[n]]) %>% 
          select(beta, p.value.spa) %>% 
          rename(p_val = p.value.spa) %>% 
          mutate(p_val = ifelse(p_val == 0, 5e-300, p_val),
                 chisq = qchisq(p_val, df = 1, lower.tail = F))
        #summarising data
        bind_cols(
          tmp_res %>% 
            filter(beta == 0) %>%
            summarise(FP = sum(p_val < 5e-8),
                      mean_null_chisq = mean(chisq)),
          tmp_res %>% 
            filter(beta != 0) %>%
            summarise(causal = sum(p_val < 5e-8),
                      mean_causal_chisq = mean(chisq)))
      }
      
    }) %>% 
      do.call("bind_rows", .) %>% 
      mutate(method = c("linear", "SPACox"),
             v = v)
    
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
  str_subset("N1000k_M20k") %>% 
  str_subset(negate = T, pattern = "part") %>% 
  str_subset("noFH")

tmp2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, paste0(prev_text_vec, "_", collapse = "|"))),
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
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_LTM_N1000k_M20k_normal_", C,"_", prev_txt, "_noFH", if(ds == "yes") "_DS" ,".rds"))
    #getting info per method, v, maf and generative model
    tmp_res = snpinfo %>% select(beta) %>%
      mutate(p_val = predict(gwas, log10 = F),
             p_val = ifelse(p_val == 0, 5e-300, p_val),
             chisq = qchisq(p_val, df = 1, lower.tail = F))
    #summarising data
    bind_cols(
      tmp_res %>%
        filter(beta == 0) %>%
        summarise(FP = sum(p_val < 5e-8),
                  mean_null_chisq = mean(chisq)),
      tmp_res %>% 
        filter(beta != 0) %>%
        summarise(causal = sum(p_val < 5e-8),
                  mean_causal_chisq = mean(chisq))) %>% 
      mutate(v = v)
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
        "./results/ltm_gwas_N1000k_M20k_results_all.rds")
