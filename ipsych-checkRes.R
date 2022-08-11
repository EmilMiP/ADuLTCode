library(data.table)
library(dplyr)
library(ggplot2)
library(bigsnpr)
library(stringr)
library(purrr)

source("./function-manhattan.R")
manhattan_helper = function(gwas_path, phenotype, method, ld_snps_path, col_pval = "p_value") {
  gwas = format_help(gwas_path = gwas_path, method = method)
  ld_clumped_indx = readRDS(ld_snps_path)
  if(str_detect(method, "Case")) {
    title_string = paste0(method, " - ", toupper(phenotype))
  } else if(str_detect(method, "ADuLT")) {
    title_string = paste0(method, " - ", toupper(phenotype))
  } else {
    title_string = paste0(method, " - ", toupper(phenotype))
    col_pval = "p.value.spa"
  }
  print(title_string)
  manhattanPlot_noSave(data = gwas,
                       title = title_string,
                       pvalColName = col_pval, 
                       x_text_size = 11, 
                       basePairName = "POS",
                       ld_clumped_indx = ld_clumped_indx,
                       SNP_distance = 500)[[1]]
}

format_help = function(gwas_path, method) {
  if (str_detect(method, "SPACox")) {
    readRDS(gwas_path) %>%
      as_tibble() %>%
      bind_cols("CHR" = CHR, "POS" = POS)
  } else {
    tmp = readRDS(gwas_path) 
    tmp %>% mutate(p_value = predict(tmp,log10 = F),
                   "CHR" = CHR,
                   "POS" = POS) %>% 
      as_tibble()
  }
}

# # Loading Genetic Data ----------------------------------------------------
#obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")
#G <- obj.bigsnp$genotypes
map = fread("../results/export_map.txt", select = c("CHR", "POS"))
CHR <- map$CHR
POS <- map$POS

spacox = tibble(gwas_path = list.files("../export3/", "ipsych_SPA", full.names = T),
                phenotype = str_match(gwas_path, "adhd|asd|dep|scz")[,1],
                method    = "SPACox",
                noAge     = ifelse(str_detect(gwas_path, "noAge"), "yes", "no")) %>% 
  relocate(method, phenotype, noAge, gwas_path)



# ipsych = readRDS("./results/ipsych_gwas_thirdDegree.rds")
# ipsych = readRDS("./results/ipsych_gwas_noAge_thirdDegree.rds")
# names(ipsych) = c("status_adhd", "adhd",
#                   "status_asd", "asd",
#                   "status_dep", "dep",
#                   "status_scz", "scz")#names(pheno)[-c(1,2)]
# str(ipsych)


ipsych = tibble(gwas_path = list.files("../export3", "ipsych", full.names = T),
                phenotype = str_match(gwas_path, "adhd|asd|dep|scz")[,1],
                method    = str_match(gwas_path, "CaseControl|ADuLT")[,1],
                noAge     = ifelse(str_detect(gwas_path, "noAge"), "yes", "no")) %>% 
  filter(!str_detect(gwas_path, "SPACox|gwas|ldclump")) %>% 
  relocate(method, phenotype, noAge, gwas_path)


ld_paths = tibble(ld_snps_path = list.files("../export3/singleTraitLDs", full.names = T),
                  phenotype = str_match(ld_snps_path, "adhd|asd|dep|scz")[,1],
                  noAge = ifelse(str_detect(ld_snps_path, "noAge"), "yes", "no"),
                  method = ifelse(str_detect(ld_snps_path, "SPACox"), "SPACox", str_match(ld_snps_path,"CaseControl|ADuLT")[,1]))
# 
# sipsych = lapply(seq_along(ipsych), function(n) {
#   ipsych[[n]] %>%
#     as_tibble() %>%
#     bind_cols("CHR" = CHR, "BP" = POS) %>%
#     mutate(p_val = 10^(predict(ipsych[[n]])))
# })
# names(sipsych) = names(ipsych)

# # loading Bolt results ----------------------------------------------------
# 
# ipsych_files = list.files("./results", "stats", full.names = T) %>% 
#   str_subset("ipsych") %>% 
#   str_subset("thirdDegree")
# 
# ipsych = lapply(ipsych_files, function(x) fread(x) %>% as_tibble())
# 
# disorder = str_split(ipsych_files, "ipsych-|-adult|-Case") %>% 
#   sapply(function(x) x[2])
# method = as.vector(str_match(ipsych_files, "adult|CaseControl"))
# 
# names(ipsych) = paste0(disorder,"_", method, ifelse(str_detect(ipsych_files, "noAge"), "_noAge", "") )




# linear reg output -------------------------------------------------------
#set new y value across a list of ggplots
set_new_y = function(plots, new_y) {
  lapply(plots, function(plt) {
    pvalColName = str_subset(names(plt$data), "p.value.spa|p_value")
    plt + 
      ylim(c(min(-log10(plt$data[[pvalColName]])), new_y))
  })
}
#get the max y value from a list of ggplots
get_max_y = function(plots) {
  sapply(plots, function(plt){
    pvalColName = str_subset(names(plt$data), "p.value.spa|p_value")
    max(-log10(plt$data[[pvalColName]]))
  }) %>% max + 2
}


linRegPlots = bind_rows(ipsych, spacox) %>%
  left_join(ld_paths) %>% 
  mutate(manhattanplots = purrr::pmap(.l = list(phenotype, method, gwas_path, ld_snps_path),
                                      ~ manhattan_helper(phenotype = ..1, 
                                                         method = ..2,
                                                         gwas_path = ..3,
                                                         ld_snps_path = ..4)))




linRegCombined = linRegPlots %>% 
  group_by(phenotype, noAge) %>%
  summarise(plots = list(manhattanplots)) %>% 
  mutate(new_y_max = purrr::map_dbl(.x = plots,
                                    .f = ~ get_max_y(.x)),
         plots = purrr::map2(.x = plots,
                             .y = new_y_max,
                             .f = ~ set_new_y(.x, .y)),
         man = purrr::map(.x = plots, ~ do.call("plot_grid", c(.x, list(ncol = 1))))) %>% 
  ungroup()


linRegCombined$man[[1]]
linRegCombined$man[[2]]
linRegCombined$man[[3]]
linRegCombined$man[[4]]
linRegCombined$man[[5]]
linRegCombined$man[[6]]
linRegCombined$man[[7]]
linRegCombined$man[[8]]

for (i in 1:nrow(linRegCombined)) {
  ggsave(paste0("../plots/ipsych/ipsych-", 
                linRegCombined$phenotype[i], 
                if(linRegCombined$noAge[i] == "yes") "-noAge" , "-manhattan-all.png"),
         linRegCombined$man[[i]],
         width = 5.5,
         height = 5.5)
}


adhd_presenting_plot = plot_grid(linRegCombined$plots[[2]][[1]] + 
                                   ggtitle("ADuLT - ADHD"),
                                 linRegCombined$plots[[2]][[3]] +
                                   ggtitle("SPACox - ADHD"),
                                 linRegCombined$plots[[1]][[2]],
                                 ncol = 1)

ggsave("../plots/ipsych/adhd_presenting_plot.png",
       adhd_presenting_plot,
       width = 5.5,
       height = 5.5)

ggsave("../plots/ipsych/Figure1_manhattan.png",
       linRegCombined$plots[[2]][[1]],
       width = 5,
       height = 3)

# manhattan plots with bolt output ----------------------------------------

boltPlots = bind_rows(tibble(phenotype = names(ipsych),
                             method    = ifelse(str_detect(phenotype, "CaseControl"), "Case-Control", "ADuLT"),
                             gwas = ipsych) %>% 
                  mutate(phenotype = str_replace(phenotype, "_adult|_CaseControl", "")),
                  spacox %>% 
                    mutate(phenotype = str_replace(phenotype, "-noAge", "_noAge"))) %>% 
  mutate(manhattanplots = purrr::pmap(.l = list(phenotype, method, gwas),
                                      ~ manhattan_helper(phenotype = ..1, 
                                                         method = ..2,
                                                         .tbl = ..3,
                                                         col_pval = "P_BOLT_LMM_INF"))) %>% 
  select(-gwas)

boltPlots$manhattanplots
boltPlotsCombined = boltPlots %>% 
  group_by(phenotype) %>%
  summarise(plots = list(manhattanplots)) %>% 
  mutate(man = purrr::map(.x = plots, ~ do.call("plot_grid", c(.x, list(ncol = 1)))))

boltPlotsCombined$man[[1]]
boltPlotsCombined$man[[2]]



for (i in 1:nrow(boltPlotsCombined)) {
  ggsave(paste0("./plots/tmp-", boltPlotsCombined$phenotype[i] ,"-manhattan-all-bolt.png"),
         boltPlotsCombined$man[[i]],
         width = 4.5,
         height = 6)
}


adhd_presenting_plot_bolt = plot_grid(boltPlotsCombined$plots[[2]][[2]] +
                                        ggtitle("ADuLT - ADHD"),
                                      boltPlotsCombined$plots[[2]][[3]] +
                                        ggtitle("SPACox - ADHD"),
                                      boltPlotsCombined$plots[[1]][[1]],
                                      ncol = 1)

ggsave("./plots/ipsych/adhd_presenting_plot_bolt.png",
       adhd_presenting_plot_bolt,
       width = 4.5,
       height = 6)


# lambdaGC ----------------------------------------------------------------

get_lambdaGC = function(gwas_path) {
  if(str_detect(gwas_path, "SPACox", negate = T)) {
    gwas = fread(gwas_path) %>% as_tibble()
    median((gwas$score)^2 / qchisq(0.5, 1))
  } else {
    gwas = readRDS(gwas_path) %>% as_tibble()
    median((gwas$z)^2 / qchisq(0.5, 1))
  }
}
ipsych_gc = ipsych %>% 
  mutate(lambdaGC = purrr::map_dbl(.x = gwas_path,
                                   .f = ~ get_lambdaGC(.x)))

spacox_gc = spacox %>% 
  mutate(lambdaGC = purrr::map_dbl(.x = gwas_path,
                                   .f = ~ get_lambdaGC(.x)))
