library(LTFHPlus)
library(dplyr)


itrs = expand.grid(list(itr = 1:10,
                        ncores = c(1, 2, 4, 8, 12, 16, 24, 32),
                        N = c(10e3, 50e3, 100e3)))


NCORES = max(itrs$ncores)
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = paste0(128 / 16 * 1, "g"),
  name = basename(rstudioapi::getSourceEditorContext()$path))))

tmp <- future.apply::future_lapply(1:nrow(itrs), function(i) {
  N = itrs$N[i]
  h2 = .5
  
  #calculates the sex-specific prevalences used to determine status:
  multiplier = 1
  prev = 0.05
  
  #simulate liabilities
  liabs = MASS::mvrnorm(n = N, mu = rep(0, 2), Sigma = matrix(c(h2, h2, h2, 1), ncol = 2, nrow = 2))
  
  #Age of children, parents, and siblings
  phen = tibble(
    ids = 1:N,
    child_gen  = liabs[,1],
    child_full = liabs[,2],
    age        = runif(N, 10, 100),
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
  
  #Setting up parallelization backend
  plan(tweak(multisession, workers = itrs$ncores[i]))
  #performs LT-FH++ analysis
  system.time(estimate_gen_liability_noFH(phen = phen, h2 = h2))
}, future.seed = T)

saveRDS(tmp, 
        "./results/computation_times.rds")


