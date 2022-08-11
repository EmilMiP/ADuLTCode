library(dplyr)
library(data.table)
library(stringr)
library(xlsx)
# helper functions --------------------------------------------------------

#get comparison to baseline phenotype:
get_diff = function(.tbl, baseline_method = "SPACox", ...) {
  baseline = filter(.tbl, str_detect(method, baseline_method)) %>% pull(causal)
  .tbl %>% mutate(diff = causal - baseline, 
                  ratio = causal / baseline)
}

# read data ---------------------------------------------------------------

#load the summary info
res = readRDS("../results/cox_gwas_N1000k_M20k_results_all.rds") %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         method = ifelse(str_detect(method, "linear"), "Case-Control", method),
         color = ifelse(str_detect(method, "ADuLT"), "red", 
                        ifelse(str_detect(method, "Case-Control"), "green", "blue")),
         prev = ifelse(prev == "prev5", "prev05", prev)) %>% 
  select(-beta_gen) %>% 
  group_by(v, C, prev, downsampling) %>% 
  group_map(.data = ., ~ get_diff(.tbl = .x), .keep = T) %>% 
  do.call("bind_rows", .) %>% 
  mutate(method = factor(method, levels = c("ADuLT", "SPACox", "Case-Control"))) 

#  ltm --------------------------------------------------------------------
#load the summary info
ltm = readRDS("../results/ltm_gwas_N1000k_M20k_results_all.rds") %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         method = ifelse(str_detect(method, "linear"), "Case-Control", method),
         color = ifelse(str_detect(method, "ADuLT"), "red", 
                        ifelse(str_detect(method, "Case-Control"), "green", "blue")),
         prev = ifelse(prev == "prev5", "prev05", prev)) %>% 
  group_by(v, C, prev, downsampling) %>% 
  group_map(.data = ., ~ get_diff(.tbl = .x, baseline_method = "SPACox"), .keep = T) %>% 
  do.call("bind_rows", .) %>% 
  mutate(method = factor(method, levels = c("ADuLT", "SPACox", "Case-Control")))



# info --------------------------------------------------------------------

combined = bind_rows(mutate(res, gen_mod = "Hazard"),
                     mutate(ltm, gen_mod = "Liability")) %>% 
  relocate(v, C, prev, method, downsampling, gen_mod)

all_summary_info = combined %>% 
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
  group_by(C, prev, method, downsampling, gen_mod) %>% 
  summarise(n_causal  = list(causal),
            obs       = list(ratio),
            power_obs = list(power),
            null_obs  = list(mean_null_chisq),
            avg_ratio    = purrr::map_dbl(obs, ~ round(mean(.x, na.rm = T),2)),
            avg_ratio_se = purrr::map_dbl(obs, ~ round(sd(.x, na.rm = T)/sqrt(length(n)),5)),
            avg_power    = purrr::map_dbl(power_obs, ~ round(mean(.x, na.rm = T), 5)),
            avg_power_se = purrr::map_dbl(power_obs, ~ round(sd(.x, na.rm = T)/sqrt(length(.x)),5)),
            avg_null     = purrr::map_dbl(null_obs, ~ round(mean(.x, na.rm = T), 2)),
            avg_null_se  = purrr::map_dbl(null_obs, ~ round(sd(.x, na.rm = T)/sqrt(length(.x)),5)),
            no_causal    = purrr::map_dbl(n_causal, ~ sum(.x == 0)),
            .groups = "drop") %>% 
  select(-n_causal, -obs, -power_obs, -null_obs) %>% 
  arrange(gen_mod, C, prev, downsampling) %>%
  relocate(gen_mod, C, prev, downsampling, method) 




# excel -------------------------------------------------------------------


#raw simulation results, per iteration
write.xlsx(as.data.frame(combined),
           "../results/summary_information.xlsx",
           sheetName = "Raw Simulation Results", 
           row.names = F)

# averaged simulation results, per parameter setup
write.xlsx(as.data.frame(all_summary_info),
           "../results/summary_information.xlsx",
           sheetName = "Averaged Simulation Results", 
           row.names = F, 
           append = T)

#average computation times per method
write.xlsx(as.data.frame(plot_mean_times),
           "../results/summary_information.xlsx",
           sheetName = "Computation Times", 
           row.names = F, 
           append = T)
