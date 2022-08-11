library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)


spacox = readRDS("./results/cox_causal_snps_all.rds")

base_pvals = spacox %>% 
  filter(beta_gen == "normal",
         C == "C1000",
         event_rate == "ER15",
         downsampling == "no",
         method == "ADuLT") 

spacox %>% mutate(lin_reg_p_val = rep(base_pvals)
#tidyr::pivot_longer(spacox, col = )

spacox %>% 
  filter(beta_gen == "normal",
         C == "C1000",
         event_rate == "ER15",
         downsampling == "no",
         method == "GWAS") %>%
  ggplot(aes(x = base_pvals$z_val, y = z_val)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm",
              formula = "y ~ x + 0")

spacox %>% 
  filter(beta_gen == "normal",
         C == "C1000",
         event_rate == "ER15",
         downsampling == "no",
         method == "SPACox") %>%
  ggplot(aes(x = base_pvals$z_val, y = z_val)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm",
              formula = "y ~ x + 0")
