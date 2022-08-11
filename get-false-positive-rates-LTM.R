library(data.table)
library(dplyr)
library(stringr)
library(bigstatsr)




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

files = tibble(ph = all_files) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, 
                                    paste0(rep(prev_text_vec, each = 2), c("\\.rds", "_DS.rds"), 
                                           collapse = "|", sep = ""))),
         prev = str_replace_all(prev, "\\.rds|_DS\\.rds", ""),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup

phenotypes = c("linear", "SPACox")


fp = lapply(1:nrow(files), function(ctr) { #nrow(tmp)
  C        = files$C[[ctr]]
  prev_txt = files$prev[[ctr]]
  ds       = files$downsampling[[ctr]]
  
  sum_info2 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, 
                          "_LTM_N1000k_M20k_normal_", 
                          C,"_", 
                          prev_txt, 
                          if (ds == "yes") "_DS" ,".rds"))
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
          do.call("rbind", .) %>%
          mutate(method = phenotypes[ctr2])
      }
      
    }) %>% 
      do.call("rbind", .) %>%
      as_tibble() %>% 
      mutate(v = v)
    
  }) %>% 
    do.call("rbind", .) %>% 
    mutate(C = C,
           prev = prev_txt,
           downsampling = ds)
  
 
}) %>% 
  do.call("rbind", .) %>%
  as_tibble() %>% 
  relocate(v, C, prev, downsampling)







# adult results -----------------------------------------------------------

#liability based method results:
all_files2 = list.files("./gwas-results", full.names = T, "LTM") %>% 
  str_subset("N1000k_M20k") %>% 
  str_subset(negate = T, pattern = "part") %>% 
  str_subset("noFH")

files2 = tibble(ph = all_files2) %>%   
  mutate(C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, paste0(prev_text_vec, "_", collapse = "|"))),
         prev = str_replace(prev, "_", ""),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no")) %>% 
  group_by(C, prev, downsampling) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup


phenotypes2 = c("ADuLT")

fp2 = lapply(1:nrow(files2), function(ctr) {
  C        = files2$C[[ctr]]
  prev_txt = files2$prev[[ctr]]
  ds       = files2$downsampling[[ctr]]
  
  sum_info2 = lapply(1:10, function(v) {
    snpinfo = fread(paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-", C, "-v", v, ".snpinfo")) %>% 
      as_tibble
    gwas = readRDS(paste0("./gwas-results/v", v, 
                          "_LTM_N1000k_M20k_normal_", 
                          C,"_", 
                          prev_txt, "_noFH", 
                          if(ds == "yes") "_DS" ,".rds"))
    
    #sum_info = lapply(seq_along(gwas), function(ctr2) {
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
        do.call("rbind", .) %>%
        mutate(method = phenotypes2,
               v = v)
    }) %>%
    do.call("rbind", .) %>%
    mutate(C = C,
           prev = prev_txt,
           downsampling = ds)
  }) %>% 
  do.call("rbind", .) %>%
  as_tibble() %>% 
  relocate(v, C, prev, downsampling)


# save data ---------------------------------------------------------------

saveRDS(bind_rows(fp, fp2),
        "./results/false-positive-data-LTM.rds")

