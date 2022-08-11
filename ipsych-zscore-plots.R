source("./code/plot-functions.R")

format_data = function(all_files, cur_phenotype) {
  print(cur_phenotype)
  lapply(all_files, function(path) {
    print(path)
    cur_method = str_match(path, paste0(c("adult", "CaseControl"), collapse = "|"))
    print(cur_method)
    fread(path) %>% 
      as_tibble() %>% 
      select(SNP, CHR, BP, 
             !!as.symbol(paste0("beta_", cur_method)) := BETA,
             !!as.symbol(paste0("SE_", cur_method)) := SE,
             !!as.symbol(paste0("CHISQ_BOLT_LMM_INF_", cur_method)) := CHISQ_BOLT_LMM_INF,
             !!as.symbol(paste0("P_BOLT_LMM_INF_", cur_method)) := P_BOLT_LMM_INF)
  }) %>% do.call("inner_join", .)
}



phenotypes = c("adhd", "asd", "dep", "scz")
methods = c("adult", "CaseControl")

sumstats = tibble(files = list.files("./results", "stats", full.names = T) %>% 
  str_subset("ipsych")) %>% 
  mutate(phenotype   = str_match(files, paste0(phenotypes, collapse = "|"))[,1], #apparently it returns a matrix
         method      = str_match(files, paste0(methods, collapse = "|"))[,1],
         noAge       = str_detect(files, "noAge") + 0L,
         thirdDegree = str_detect(files, "thirdDegree") + 0L ) %>% 
  filter(thirdDegree != 0) %>% 
  group_by(phenotype, noAge, thirdDegree) %>%
  summarise(all_files = list(files), .groups = "drop") %>% 
  mutate(comb = purrr::pmap(.l = list(phenotype = phenotype,
                                      all_files = all_files),
                            .f = ~ format_data(all_files = ..2,
                                               cur_phenotype = ..1)),
         zscores = purrr::map(.x = comb, 
                              .f = ~ zscore_plot(data = .x,
                               suffix2 = "adult",
                               suffix1 = "CaseControl",
                               BetaColName = "beta",
                               y_label = "Z-score (ADuLT)",
                               x_label = "Z-score (Case-Control)",
                               filter_sugg = T)),
         chisqs = purrr::map(.x = comb,
                             .f = ~ chisq_plot(data = .x,
                                               suffix2 = "adult",
                                               suffix1 = "CaseControl",
                                               y_label = expression(chi^2 ~ "- Statistic (ADuLT)"),
                                               x_label = expression(chi^2 ~ "- Statistic (Case-Control)"),
                                               filter_sugg = T)))

sumstats$zscores[[2]]
sumstats$zscores[[4]]
sumstats$zscores[[6]]
sumstats$zscores[[8]]
sumstats$chisqs[[2]]
sumstats$chisqs[[4]]
sumstats$chisqs[[6]]
sumstats$chisqs[[8]]



# Get LambdaGC ------------------------------------------------------------

get_lambdaGC = function(log) {
  cur_log = readLines(log)
  str_subset(cur_log, "lambdaGC") %>% 
    str_subset("BOLT") %>% 
    str_split("lambdaGC:") %>% sapply(., function(x) as.numeric(x[2]))
}



logs = tibble(files = list.files("./results", "log", full.names = T) %>% 
                str_subset("ipsych")) %>% 
  mutate(phenotype   = str_match(files, paste0(phenotypes, collapse = "|"))[,1], #apparently it returns a matrix
         method      = str_match(files, paste0(methods, collapse = "|"))[,1],
         noAge       = str_detect(files, "noAge") + 0L,
         thirdDegree = str_detect(files, "thirdDegree") + 0L ) %>% 
  filter(thirdDegree != 0) %>% 
  mutate(lambdaGC = purrr::map_dbl(files, ~ get_lambdaGC(.x)))



# LDclumped snps ----------------------------------------------------------


phenotypes = c("adhd", "asd", "dep", "scz")
methods = c("adult", "CaseControl")

ld_files = tibble(
  files_LD = list.files("./LDClumpedSNPs", "clumped-snps", full.names = T)
) %>% 
  mutate(phenotype = str_match(files_LD, paste0(phenotypes, collapse = "|"))[,1])


format_data_LD = function(all_files, LD_file, cur_phenotype) {
#  print(cur_phenotype)
  
  LD_snps = readRDS(LD_file) %>% as_tibble() %>% pull(SNP)
  
  lapply(all_files, function(path) {
    print(path)
    cur_method = str_match(path, paste0(c("adult", "CaseControl"), collapse = "|"))
    print(cur_method)
    fread(path) %>% 
      as_tibble() %>% 
      filter(SNP %in% LD_snps) %>% 
      select(SNP, CHR, BP, 
             !!as.symbol(paste0("beta_", cur_method)) := BETA,
             !!as.symbol(paste0("SE_", cur_method)) := SE,
             !!as.symbol(paste0("CHISQ_BOLT_LMM_INF_", cur_method)) := CHISQ_BOLT_LMM_INF,
             !!as.symbol(paste0("P_BOLT_LMM_INF_", cur_method)) := P_BOLT_LMM_INF)
  }) %>% do.call("inner_join", .)  
}

sumstats_ld = tibble(files = list.files("./results", "stats", full.names = T) %>% 
                    str_subset("ipsych")) %>% 
  mutate(phenotype   = str_match(files, paste0(phenotypes, collapse = "|"))[,1], #apparently it returns a matrix
         method      = str_match(files, paste0(methods, collapse = "|"))[,1],
         noAge       = str_detect(files, "noAge") + 0L,
         thirdDegree = str_detect(files, "thirdDegree") + 0L ) %>% 
  filter(thirdDegree != 0) %>% 
  left_join(., ld_files) %>% 
  group_by(phenotype, noAge, thirdDegree, files_LD) %>%
  summarise(all_files = list(files), .groups = "drop") %>% 
  mutate(comb = purrr::pmap(.l = list(phenotype = phenotype,
                                      ld_file = files_LD,
                                      all_files = all_files),
                            .f = ~ format_data_LD(all_files = ..3,
                                                  LD_file = ..2,
                                                  cur_phenotype = ..1)),
         zscores = purrr::map(.x = comb, 
                              .f = ~ zscore_plot(data = .x,
                                                 suffix2 = "adult",
                                                 suffix1 = "CaseControl",
                                                 BetaColName = "beta",
                                                 y_label = "Z-score (ADuLT)",
                                                 x_label = "Z-score (Case-Control)",
                                                 filter_sugg = F)),
         chisqs = purrr::map(.x = comb,
                             .f = ~ chisq_plot(data = .x,
                                               suffix2 = "adult",
                                               suffix1 = "CaseControl",
                                               y_label = expression(chi^2 ~ "- Statistic (ADuLT)"),
                                               x_label = expression(chi^2 ~ "- Statistic (Case-Control)"),
                                               filter_sugg = F)))


sumstats_ld$zscores[[2]]
sumstats_ld$zscores[[4]]
sumstats_ld$zscores[[6]]
sumstats_ld$zscores[[8]]
sumstats_ld$chisqs[[2]]
sumstats_ld$chisqs[[4]]
sumstats_ld$chisqs[[6]]
sumstats_ld$chisqs[[8]]
