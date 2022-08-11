library(data.table)
library(dplyr)
library(stringr)
library(bigsnpr)
library(LTFHPlus)
library(future.batchtools)

all_files = list.files("./simulatedData", "LTM", full.names = T) %>% 
  str_subset("-N1000k")

prev_vec = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5)
C_vec  = c(250, 1000)
#constructing prevalence name to keep file names ordered
prev_text_vec = paste0("prev",round(prev_vec * 100,1))
prev_text_vec[1] = "prev01"
prev_text_vec[2] = "prev025"

phenotypes = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "phen")
  tibble(full_names = ph,
         v = v) %>% 
    mutate(C = as.vector(str_match(full_names, "C1000|C250")),
           prev = as.vector(str_match(full_names, paste0(prev_text_vec, "-", collapse = "|")))) %>% 
    group_by(v, C, prev) %>% 
    summarise(true = full_names, .groups = "drop") %>% 
    ungroup
}) %>% 
  do.call("rbind", .)

files = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "phen", negate = T)
  tibble(full_names = ph,
         v = v) %>% 
    mutate(C = as.vector(str_match(full_names, "C1000|C250"))) %>% 
    group_by(v, C) %>% 
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
  relocate(v, C, prev, downsampling) %>% 
  mutate(prev = str_replace(prev, "-", "")) #multiple keeper files are present, but none are used.

NCORES = 8
#files = files %>% filter(C == "C1000", v == "7", prev == "prev5") 
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "8g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


future.apply::future_lapply(1:nrow(files), function(ctr) {
  #getting identifiers
  v = files$v[[ctr]]
  C = files$C[[ctr]]
  prev_txt = files$prev[[ctr]]
  downsampling = files$downsampling[[ctr]]
  #loading data
  bigsnp_obj   = snp_attach(str_subset(files$grps[[ctr]], "rds"))
  G = bigsnp_obj$genotypes
  
  #is downsampling a posibility ?
  if (downsampling) {
    keepers      = fread(str_subset(files$grps[[ctr]], "keepers")) %>% as_tibble
    KEEP = which(1:nrow(G) %in% keepers[[1]])
    phen = fread(str_subset(files$grps[[ctr]], "phen")) %>% 
      as_tibble() %>% 
      slice(keepers[[1]])
  } else {
    KEEP = 1:nrow(G)
    phen = fread(str_subset(files$grps[[ctr]], "phen")) %>% 
      as_tibble 
  }
  
  intervals = bigparallelr::split_len(ncol(G), nb_split = 20)
  res_files = paste0("./gwas-results/v",v, "_LTM_N1000k_M20k_normal_",
                     C, "_",
                     prev_txt, 
                     "_noFH",
                     if (downsampling) "_DS",
                     "_part", rows_along(intervals), ".rds")
  mean(already_done <- file.exists(res_files))

  all_gwas_parts = lapply(which(!already_done), function(ic) {
    #indecies are determined
    ind = seq(intervals[ic, "lower"], intervals[ic, "upper"])

    #linear regression:
    gwas = big_univLinReg(G, phen$post_gen_no_fam, ind.col = ind, ncores = NCORES, ind.train = KEEP)

    saveRDS(gwas, res_files[ic])
    NULL
  })
  NULL
}, future.seed = T)

