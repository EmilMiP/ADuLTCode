library(bigsnpr)
library(dplyr)
library(stringr)
library(data.table)

#obj.bigsnp <- snp_attach("./genData/dosage_ipsych2015.rds")
#G <- obj.bigsnp$genotypes
map = fread("../results/export_map.txt") %>% as_tibble()

CHR <- map$CHR
POS <- map$POS
SNP = map$SNP


# ipsych_gwas = list.files("../results", "rds", full.names = T) %>% 
#   str_subset("ipsych_gwas") %>% 
#   str_subset("thirdDegree") %>% 
#   lapply(., function(x) tibble(gwas = readRDS(x),
#                                phenotype = names(gwas),
#                                noAge = ifelse(str_detect(x, "noAge"), "yes", "no"),
#                                method = ifelse(str_detect(phenotype, "status"), "CaseControl", "ADuLT")) %>% 
#            mutate(phenotype = str_replace(phenotype, "status_", ""))) %>% 
#   do.call("bind_rows", .)

filter_help = function(gwas_path, ld_snps_path, method) {
  ld_snps = readRDS(ld_snps_path)
  if (method != "SPACox") {
    tmp = readRDS(gwas_path) 
    tmp %>% mutate(p_value = -predict(tmp),
                   "SNP" = SNP,
                   "CHR" = CHR,
                   "POS" = POS) %>% 
      as_tibble() %>%
      slice(ld_snps) %>% 
      filter(p_value > -log10(5e-8))
  } else {
    readRDS(gwas_path) %>%
      mutate("SNP" = SNP,
             "CHR" = CHR,
             "POS" = POS) %>% 
      slice(ld_snps) %>% 
      filter(p.value.spa < 5e-8)
  }
}

ld_paths = tibble(ld_snps_path = list.files("../export3/singleTraitLDs", full.names = T),
                  phenotype = str_match(ld_snps_path, "adhd|asd|dep|scz")[,1],
                  noAge = ifelse(str_detect(ld_snps_path, "noAge"), "yes", "no"),
                  method = ifelse(str_detect(ld_snps_path, "SPACox"), "SPACox", str_match(ld_snps_path,"CaseControl|ADuLT")[,1]))


res = tibble(gwas_path = list.files("../export3", "ADuLT|CaseControl|SPACox", full.names = T),
             phenotype = str_match(gwas_path, "adhd|asd|dep|scz")[,1],
             noAge = ifelse(str_detect(gwas_path, "noAge"), "yes", "no"),
             method = ifelse(str_detect(gwas_path, "SPACox"), "SPACox", str_match(gwas_path,"CaseControl|ADuLT")[,1]))  %>% 
  relocate(method, noAge, phenotype, gwas_path) %>% 
  left_join(ld_paths) %>% 
  mutate(gwas = purrr::pmap(.l = list(gwas_path = gwas_path, ld_snps_path = ld_snps_path, method),
                            .f = ~ filter_help(..1, ..2, ..3)))

#how many sig snps does each analysis have after ldclumping ?
res %>% 
  filter((noAge == "no" & method == "CaseControl") | (noAge == "yes" & method != "CaseControl")) %>% 
  mutate(n_sig = sapply(gwas, nrow)) %>% 
  select(-gwas, -ld_snps_path, -gwas_path) %>% 
  relocate(noAge, phenotype) %>% 
  arrange(phenotype) %>% 
  as.data.frame()

#total LD clumped snps per method
res %>% 
#  filter((noAge == "no" & method == "CaseControl") | (noAge == "yes" & method != "CaseControl")) %>% 
  mutate(n_sig = sapply(gwas, nrow)) %>% 
  select(-gwas, -ld_snps_path, -gwas_path) %>% 
  group_by(noAge, method) %>% 
  summarise(tot_sig = sum(n_sig)) %>% 
  as.data.frame()



res$gwas[[9]] %>% 
  mutate("CHR:POS" = paste0(CHR, ":", POS),
         beta = paste0(round(estim, 4), "(", round(std.err, 4), ")")) %>% 
  select(SNP, "CHR:POS", beta, p_value)


# SNP        `CHR:POS`   beta            p_value
# <chr>      <chr>       <chr>             <dbl>
# rs11210887 1:44076019  0.0331(0.0049)    10.8 
# rs4660756  1:44383914  0.0284(0.0051)     7.48
# rs7563362  2:620297    -0.0361(0.0065)    7.47
# rs4916723  5:87854395  -0.0359(0.0047)   13.7 
# rs12705966 7:114248851 0.0334(0.0052)     9.86
# rs13236619 7:157827565 -0.0263(0.0046)    8.05
# rs72673548 8:93292844  -0.0423(0.0076)    7.55
# rs12346733 9:86727865  -0.0256(0.0047)    7.35
# rs57806515 11:28628549 0.0282(0.0047)     8.87
# rs704061   12:89771903 -0.0252(0.0046)    7.49
# rs4261436  14:33299482 -0.0257(0.0045)    7.86
# rs4813421  20:21258053 0.0305(0.005)      8.92


res$gwas[[2]] %>% 
  mutate("CHR:POS" = paste0(CHR, ":", POS),
         beta = paste0(round(estim, 4), "(", round(std.err, 4), ")")) %>% 
  select(SNP, "CHR:POS", beta, p_value)

# 
# # A tibble: 11 x 4
# SNP        `CHR:POS`   beta            p_value
# <chr>      <chr>       <chr>             <dbl>
# rs11210887 1:44076019  0.0201(0.003)     10.8 
# rs7563362  2:620297    -0.022(0.004)      7.47
# rs4916723  5:87854395  -0.0225(0.0029)   14.4 
# rs9969232  7:114158954 -0.019(0.0029)    10.0 
# rs13236619 7:157827565 -0.0167(0.0028)    8.70
# rs72673548 8:93292844  -0.0261(0.0047)    7.70
# rs10767730 11:28634862 0.0169(0.0028)     8.54
# rs704061   12:89771903 -0.0153(0.0028)    7.40
# rs4261436  14:33299482 -0.0151(0.0028)    7.40
# rs8085882  18:22743899 0.0156(0.0029)     7.30
# rs6035830  20:21265728 0.0184(0.0031)     8.74

res$gwas[[9]]
setdiff(res$gwas[[2]][,5:7], res$gwas[[9]][,5:7])
# # A tibble: 5 x 3
# SNP          CHR       POS
# <chr>      <int>     <int>
# 1 rs4660756      1  44383914
# 2 rs12705966     7 114248851 # In LD with snp from CC
# 3 rs12346733     9  86727865
# 4 rs57806515    11  28628549 # In LD with snp from CC
# 5 rs4813421     20  21258053 # In LD with snp from CC

setdiff(res$gwas[[9]][,5:7], res$gwas[[2]][,5:7])
# # A tibble: 4 x 3
# SNP          CHR       POS
# <chr>      <int>     <int>
# 1 rs9969232      7 114158954 # In LD with snp from ADuLT
# 2 rs10767730    11  28634862 # In LD with snp from ADuLT
# 3 rs8085882     18  22743899
# 4 rs6035830     20  21265728 # In LD with snp from ADuLT
