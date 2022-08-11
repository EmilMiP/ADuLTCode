library(data.table)
library(dplyr)
library(stringr)
library(bigsnpr)
library(SPACox)
library(survival)
library(future.batchtools)
library(LTFHPlus)

all_files = list.files("./simulatedData", "SPACox", full.names = T) %>% 
  str_subset(., "M1000k")


#using only a subset of files
time_files = tibble(v = 1:10,
                    files = list(all_files),
                    N = list(c(10e3, 50e3, 100e3)),
                    M = list(c(100e3, 250e3, 500e3, 1000e3))) %>% 
  tidyr::unnest(N) %>%
  tidyr::unnest(M) %>% 
  relocate(v, N, M)

NCORES = 1

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


times = future.apply::future_lapply(1:nrow(time_files), function(ctr) {
  cal_N = time_files$N[ctr]
  cur_M = time_files$M[ctr]
  #loading data
  bigsnp_obj   = snp_attach(str_subset(time_files$files[[ctr]], "rds"))
  G = bigsnp_obj$genotypes
  cur_snp_info = fread(str_subset(time_files$files[[ctr]], "snpinfo")) %>% as_tibble
  phen     = fread(str_subset(time_files$files[[ctr]], "true")) %>%
    as_tibble %>% 
    mutate(event = pmin(onset, censor)) %>% 
    select(FID, IID, status, event) %>% 
    mutate(status = status + 0L)
  
  KEEP = phen$FID %in% sample(phen$FID, size = cal_N)
  KEEP_IND = which(KEEP)
  
  KEEP_M = cur_snp_info$snpid %in% sample(cur_snp_info$snpid, size = cur_M)
  KEEP_M_IND = which(KEEP_M)
  system.time({
    
    intervals = bigparallelr::split_len(cur_M, nb_split = 20)
    
    ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
    #SPACox null model is fitted:
    null_mod = SPACox_Null_Model(formula = Surv(event, status) ~ 1, data = phen[KEEP,],
                                 pIDs = phen$FID[KEEP], gIDs =  phen$FID[KEEP])
    ### DOESNT NEED TO BE RUN MORE THAN ONCE PR DATA SET
    
      all_gwas_parts = lapply(1:nrow(intervals), function(ic) {
        #indecies are determined
        ind = seq(intervals[ic, "lower"], intervals[ic, "upper"])
       
        scox = big_apply(X = G, a.FUN = function(X, ind, null_mod, KEEP, KEEP_M_IND) {
          library(SPACox)
          ph = X[KEEP,KEEP_M_IND[ind]]
          rownames(ph) = KEEP
          colnames(ph) = KEEP_M_IND[ind]
          SPACox(null_mod, ph)
        }, ncores = NCORES, null_mod = null_mod, ind = ind, KEEP = KEEP_IND, KEEP_M_IND = KEEP_M_IND)
      })
      NULL
    
  })
  
}, future.seed = T)

time_files2 = tibble(time_files, times = times) %>% 
  select(-files) %>% 
  relocate(v, N, M)
saveRDS(time_files2, 
        "./results/computation-times-N100k-M1000k-SPACox-raw.rds")
# 
# times_files = bind_cols(time_files, times = sapply(times, function(x) x[3])) %>% 
#   select(-grps)
# saveRDS(times_files,
#         "./results/computation-times-SPACox.RDS")
