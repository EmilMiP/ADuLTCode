

sumstats = readRDS("../export3/ipsych_gwas.rds")

for (i in seq_along(sumstats)) {
  cur_name = names(sumstats)[i]
  cur_phenotype = str_match(cur_name, "adhd|asd|dep|scz")[,1]
  if (str_detect(cur_name, "status")) {
    method = "CaseControl"
  } else {
    method = "ADuLT"
  }
  saveRDS(sumstats[[i]],
          paste0("../export3/ipsych_",method, "_", cur_phenotype, ".rds"))
}


sumstats_noage = readRDS("../export3/ipsych_gwas_noage.rds")

for (i in seq_along(sumstats_noage)) {
  cur_name = names(sumstats_noage)[i]
  cur_phenotype = str_match(cur_name, "adhd|asd|dep|scz")[,1]
  if (str_detect(cur_name, "status")) {
    method = "CaseControl"
  } else {
    method = "ADuLT"
  }
  saveRDS(sumstats_noage[[i]],
          paste0("../export3/ipsych_",method, "_", cur_phenotype, "_noAge.rds"))
}

