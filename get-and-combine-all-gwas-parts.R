library(data.table)
library(dplyr)
library(stringr)


prev_vec = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5)
C_vec  = c(250, 1000)
#constructing prevalence name to keep file names ordered
prev_text_vec = paste0("prev",round(prev_vec * 100,1))
prev_text_vec[1] = "prev01"
prev_text_vec[2] = "prev025"

# Saving gwas parts into one ----------------------------------------------

all_parts = tibble(ph = list.files("./gwas-results", full.names = T, pattern = "part") %>% 
                     str_subset("N1000k_M20k"),
                   v  = sapply(str_split(ph, "v|_SPACox|_LTM"), function(x) x[2])) %>%   
  mutate(gen_model = as.vector(str_match(ph, "SPACox|LTM")),
         C = as.vector(str_match(ph, "C1000|C250")),
         gen_beta = as.vector(str_match(ph, "normal")),
         prev = as.vector(str_match(ph, paste(prev_text_vec, "_", collapse = "|", sep = ""))),
         downsampling = ifelse(str_detect(ph, "_DS"),"yes", "no"),
         method = ifelse(str_detect(ph, "_noFH"), "_noFH", "")) %>% #noFH = ADuLT
  group_by(v, gen_model, gen_beta, C, prev, downsampling, method) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = str_replace(prev, "_", ""))

all_parts %>% count(gen_model, C, prev, downsampling, method) %>% as.data.frame()
all_parts %>% count(gen_model)
all_parts %>% count(C)
all_parts %>% count(prev)
all_parts %>% count(downsampling)
all_parts %>% count(method)


lapply(1:nrow(all_parts), function(ctr) {
  print(ctr)
  raw_part_nbs = sapply(str_split(all_parts$grps[[ctr]], "part|\\.rds"), function(x) x[2]) %>% as.numeric()
  if (str_detect(all_parts$method[ctr],"noFH")) {
    comb_gwas = lapply(all_parts$grps[[ctr]][order(raw_part_nbs)], readRDS) %>% 
      do.call("rbind", .)
  } else {
    tmp_storage = lapply(all_parts$grps[[ctr]][order(raw_part_nbs)], readRDS)
    names(tmp_storage[[1]])
    seq_along(tmp_storage)
    
    comb_gwas = lapply(names(tmp_storage[[1]]), function(name) {
      lapply(seq_along(tmp_storage), function(i) {
        tmp_storage[[i]][[name]]
      }) %>% do.call("rbind", .)
    })
  }
  
  saveRDS(comb_gwas, 
          paste0("./gwas-results/v", all_parts$v[ctr],"_",
                 all_parts$gen_model[ctr],
                 "_N1000k_M20k_", 
                 all_parts$gen_beta[ctr],"_", 
                 all_parts$C[ctr], "_", 
                 all_parts$prev[ctr], 
                 all_parts$method[ctr],  
                 if (all_parts$downsampling[ctr] == "yes") "_DS",
                 ".rds"))
})
