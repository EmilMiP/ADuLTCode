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
all_files = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part|noFH") %>% 
  str_subset("N1000k_M20k")

#organizing files
tmp = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         prev = as.vector(str_match(ph, paste0(rep(prev_text_vec,each = 2), c("\\.rds", "_DS"), collapse = "|"))),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = str_replace(prev, "\\.rds|_DS", ""))

#reading files and extracting relevant info
res = lapply(1:nrow(tmp), function(ctr) { #nrow(tmp)
  C       = tmp$C[[ctr]]
  cur_gen = tmp$gen_beta[[ctr]]
  prev    = tmp$prev[[ctr]]
  ds      = tmp$downsampling[[ctr]]
  
  tmp3 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-N1000k-M20k-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_N1000k_M20k_",
                          cur_gen,"_", 
                          C,  "_", 
                          prev, 
                          if (ds == "yes") "_DS", ".rds"))
    #getting info per method, v, maf and generative model
    sum_info = lapply(seq_along(gwas), function(n) {
      if (n < 2) {
        tmp_res = snpinfo %>% select(beta) %>%
          mutate(p_val = predict(gwas[[n]], log10 = F),
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
          mutate(chisq = qchisq(p_val, df = 1, lower.tail = F))
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
           beta_gen = cur_gen,
           prev = prev,
           downsampling = ds)
}) %>% 
  do.call("bind_rows", .) %>% 
  relocate(v, beta_gen, C, prev, method, downsampling)


#liability based method results:
all_files2 = list.files("./gwas-results", full.names = T, "SPACox") %>% 
  str_subset(., negate = T, pattern = "part") %>% 
  str_subset(., "noFH") %>%
  str_subset("N1000k_M20k")

tmp2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         prev = as.vector(str_match(ph, paste(prev_text_vec, "_", collapse = "|", sep = ""))),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(gen_beta, C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = str_replace(prev, "_", ""))


res2 = lapply(1:nrow(tmp2), function(ctr) { #nrow(tmp)
  C       = tmp2$C[[ctr]]
  cur_gen = tmp2$gen_beta[[ctr]]
  prev_text    = tmp2$prev[[ctr]]
  ds      = tmp2$downsampling[[ctr]]
  
  tmp4 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-SPACox-N1000k-M20k-", cur_gen, "-", C,"-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, "_SPACox_N1000k_M20k_", cur_gen,"_", C,  "_", prev_text, "_noFH", if(ds == "yes") "_DS", ".rds"))
    #getting info per method, v, maf and generative model
    tmp_res = snpinfo %>% select(beta) %>%
      mutate(p_val = predict(gwas, log10 = F),
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
           beta_gen = cur_gen,
           prev = prev_text,
           downsampling = ds)
  
}) %>% 
  do.call("bind_rows", .) %>%   
  relocate(v, beta_gen, C, method, prev, downsampling)

#save the summary information for all methods
saveRDS(bind_rows(res,res2),
        "./results/cox_gwas_N1000k_M20k_results_all.rds")

