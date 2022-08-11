library(bigsnpr)
library(dplyr)
library(data.table)
library(future.batchtools)

NCORES = 8
N  = 1000e3
M  = 20e3
h2 = 0.5
C_vec  = c(250, 1000)

params = expand.grid(v = 1:10,
                     C_vec = C_vec)

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "32g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

# Create and store genetic data -------------------------------------------
future.apply::future_lapply(1:nrow(params), function(ctr) {
  
  i        = params[ctr, 1]
  C        = params[ctr, 2]
  out_name = paste0("./simulatedData/genotypes-LTM-N", N/1e3,"k-M", M/1e3, "k-normal-C", C, "-v", i)
  
  
  #initialize backingfile
  G = FBM.code256(nrow = N,
                  ncol = M,
                  init = 0,
                  code = c(0,1,2,rep(NA_real_, 253)),
                  backingfile = out_name)
  mafs = runif(ncol(G),0.01, 0.49)
  #filling backingfile
  tmp = big_apply(X = G, a.FUN = function(X, ind, maf) {
    X[,ind] = t(replicate(nrow(X), rbinom(n = length(ind), size = 2, prob = maf)))
    NULL
  }, ncores = NCORES, maf = mafs)
  
  
  #making map object
  map_obj = tibble(
    chr = rep(1, M),
    snp = 1:M,
    cM = 0,
    bp = 1:M
  )
  
  #making bigsnp format
  children_info = list(genotypes = G,
                       map = map_obj,
                       fam = tibble(FID = 1:N,
                                    IID = 1:N))
  #saving bigsnp format
  snp_save(children_info)
  
  
  # Sample effect & assign status -------------------------------------------
  #get sum stats, e.g. observed maf
  sum.stats = big_colstats(X = G, ncores = NCORES) %>% 
    mutate("mean" = sum / nrow(G),
           "sd"   = sqrt(var))
  #sample causal betas & assigning them 
  betas = matrix(0, nrow = M, ncol = 1)
  beta_ind = sort(sample(x = 1:M, size = C))
  betas[beta_ind, ] = rnorm(n = C, sd = sqrt(h2/C))
  
  #combining SNP information
  snp_info = bind_cols(tibble(snpid = 1:M),
                       sum.stats,
                       tibble(beta = betas[,1]))
  fwrite(snp_info, 
         paste0(out_name, ".snpinfo"),
         sep = " ",
         quote = F,
         na = "NA")
  
  #calculating genetic contribution
  gen = big_prodMat(X = G, A.col = betas,
                    center = sum.stats$mean,
                    scale  = sum.stats$sd,
                    ncores = NCORES)
  
  
  true = tibble(
    FID = 1:N,
    IID = 1:N,
    child_gen   = gen[,1],
    child_full  = child_gen + rnorm(N, sd = sqrt(1 - h2)),
    sexM = sample(0:1, size = N, replace = T))
  
  fwrite(true,
         paste0(out_name, ".true"),
         quote = F,
         sep = " ",
         na = "NA")
  
}, future.seed = T)
