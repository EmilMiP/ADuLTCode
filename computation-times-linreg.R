library(data.table)
library(dplyr)
library(stringr)
library(bigsnpr)
library(LTFHPlus)
library(future.batchtools)

all_files = list.files("./simulatedData", "SPA", full.names = T)  %>% 
  str_subset("N100k-M1000k")

#using only a subset of files
time_files = tibble(v = 1:10,
                    files = list(all_files),
                    N = list(c(10e3, 50e3, 100e3)),
                    M = list(c(100e3, 250e3, 500e3, 1000e3)),
                    ncores = list(c(4, 8, 12, 16))) %>% 
  tidyr::unnest(N) %>%
  tidyr::unnest(M) %>% 
  tidyr::unnest(ncores) %>% 
  relocate(v, N, M, ncores)


NCORES = max(time_files$ncores)

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


times = future.apply::future_lapply(1:nrow(time_files), function(ctr) {
  print(ctr)
  NCORES = time_files$ncores[ctr]
  cal_N = time_files$N[ctr]
  cur_M = time_files$M[ctr]
  #loading data
  bigsnp_obj   = snp_attach(str_subset(time_files$files[[ctr]], "rds"))
  G = bigsnp_obj$genotypes
  cur_snp_info = fread(str_subset(time_files$files[[ctr]], "snpinfo")) %>% as_tibble
  
  KEEP = 1:nrow(G) %in% sample(1:100e3, size = cal_N)
  KEEP_IND = which(KEEP)
  phen = fread(str_subset(time_files$files[[ctr]], "true")) %>% 
    as_tibble()
  
  
  KEEP_M = cur_snp_info$snpid %in% sample(cur_snp_info$snpid, size = cur_M)
  KEEP_M_IND = which(KEEP_M)
  
  system.time({
    intervals = bigparallelr::split_len(cur_M, nb_split = 20)
    all_gwas_parts = lapply(1:nrow(intervals), function(ic) {
      #indecies are determined
      ind = seq(intervals[ic, "lower"], intervals[ic, "upper"])
      
      #linear regression:
      gwas     = big_univLinReg(G, phen$status[KEEP_IND], ind.col = KEEP_M_IND[ind], ncores = NCORES, ind.train = KEEP_IND)
      
    })
  })
}, future.seed = T)



time2 = tibble(time_files, times) %>% 
  select(-files) %>%
  mutate(method = "linreg") %>% 
  rename(itr = v) %>% 
  relocate(itr, ncores, N, M, method, times)
  
saveRDS(time2,
        "./results/computation-times-N100k-M1000k-linreg-raw.rds")
