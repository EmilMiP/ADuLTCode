library(data.table)
library(dplyr)
library(stringr)
library(LTFHPlus)
library(future.batchtools)


#h2 = 0.83 #AUT
#h2 = .75 #ADHD, 2014 isabell article claims 70-80%
#h2 = .37 # DEP: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4646689/pdf/pone.0142197.pdf
#h2 = .75 # SCZ: https://pubmed.ncbi.nlm.nih.gov/28987712/

combs = tibble(
  phenotype = paste0("./phenotypes/", c("adhd.rds", "asd.rds", "dep.rds", "scz.rds")),
  h2        = c(0.75, 0.83, 0.37, 0.75)
)

NCORES = 16

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))



future.apply::future_lapply(1:nrow(combs), function(ctr) {
  phen = readRDS(combs$phenotype[ctr])
  h2   = combs$h2[ctr]
  
  
  library(future.apply)
  plan(tweak(multisession, workers = NCORES))
  
  noFH = estimate_gen_liability_noFH(phen = phen, h2 = h2) %>% 
    select(-cip_pred) %>% 
    rename(FID = pid,
           IID = sample.ID)
  
  fwrite(noFH,
         str_replace(combs$phenotype[ctr], "\\.rds", "\\.phen"),
         quote = F,
         sep = " ",
         na = "NA")
}, future.seed = T)
