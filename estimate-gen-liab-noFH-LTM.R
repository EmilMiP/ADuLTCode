library(LTFHPlus)
library(dplyr)
library(data.table)
library(stringr)
library(future.apply)
library(future.batchtools)

NCORES = 8
h2 = .5 
prev_vec = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.5)
C_vec  = c(250, 1000)
#constructing prevalence name to keep file names ordered
prev_text_vec = paste0("prev",round(prev_vec * 100,1))
prev_text_vec[1] = "prev01"
prev_text_vec[2] = "prev025"

params = expand.grid(v = 1:10,
                     C_vec = C_vec, 
                     prev_vec_ind = seq_along(prev_vec))


# list.files("./simulatedData", full.names = T, "LTM") %>% 
#   str_subset("true")


plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))



future.apply::future_lapply(1:nrow(params), function(ctr) {
  v    = params[ctr, 1]
  C    = params[ctr, 2] 
  prev = prev_vec[params[ctr, 3]]
  prev_text = prev_text_vec[params[ctr, 3]]
  
  true.file = paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-C", C, "-v", v, ".true")
  
  true = fread(true.file) %>% 
    as_tibble
  N = nrow(true)
  phen = tibble(
    ids = 1:N,
    child_gen  = true$child_gen,
    child_full = true$child_full,
    age        = runif(N, 10, 90),
    status     = (child_full > qnorm(prev, lower.tail = F)) + 0L,
    lower      = -Inf,
    upper      = Inf
  ) 
  
  #who are the cases=
  cases = which(phen$status == 1) 
  
  phen$aoo = NA
  phen$aoo[cases] = LTFHPlus::liab_to_aoo(liab = phen$child_full[cases], pop_prev = prev)
  
  phen$upper[-cases] = LTFHPlus::age_to_thres(age = phen$age[-cases], pop_prev = prev)
  phen$lower[cases] = LTFHPlus::age_to_thres(age = phen$aoo[cases], pop_prev = prev)
  
  
  phen = select(phen, ids, status, lower, upper)
  
  
  library(future.apply)
  plan(tweak(multisession, workers = NCORES))
  noFH = estimate_gen_liability_noFH(phen = phen, h2 = h2)
  
  fwrite(noFH,
         paste0("./simulatedData/genotypes-LTM-N1000k-M20k-normal-C", C, "-", prev_text, "-v", v, ".phen"),
         na = "NA",
         sep = " ",
         quote = F)
  NULL
}, future.seed = T)
