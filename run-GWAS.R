library(data.table)
library(dplyr)
library(stringr)
library(bigsnpr)
library(SPACox)
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
  phen = cur_true %>% 
    select(FID, IID, status, event) %>% 
    mutate(status = status + 0L)
  
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
  
  ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
  #SPACox null model is fitted:
  null_mod = SPACox_Null_Model(formula = Surv(event, status) ~ 1, data = phen,
                               pIDs = phen$FID, gIDs = phen$FID)
  ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
  
  intervals = bigparallelr::split_len(ncol(G), nb_split = 20)
  res_files = paste0("./gwas-results/v",v, "_SPACox_N1000k_M20k_", 
                     gen_beta, "_", 
                     cur_C, "_", 
                     prev, 
                     if (downsampling) "_DS",
                     "_part", rows_along(intervals), ".rds")
  mean(already_done <- file.exists(res_files))
  
  all_gwas_parts = lapply(1:nrow(intervals), function(ic) {
    #indecies are determined
    ind = seq(intervals[ic, "lower"], intervals[ic, "upper"])
    
    #linear and logistic regression:
    gwas     = big_univLinReg(G, phen$status, ind.col = ind, ncores = NCORES, ind.train = KEEP)
    #gwas_log = big_univLogReg(G, phen$status, ind.col = ind, ncores = NCORES)
    
    #SPAcox association test
    scox = big_apply(X = G, a.FUN = function(X, ind, null_mod, phen, ind_keep) {
      library(SPACox)
      ph = X[ind_keep, ind]
      rownames(ph) = phen$FID
      colnames(ph) = ind
      SPACox(null_mod, ph)
    }, ncores = NCORES, null_mod = null_mod, ind = ind, phen = phen, ind_keep = KEEP) %>% 
      do.call("rbind",.) %>% 
      as_tibble()
    res = list("GWAS" = gwas,
               "SPAcox" = scox)
    saveRDS(res, res_files[ic])
    NULL
  })

}, future.seed = T)




