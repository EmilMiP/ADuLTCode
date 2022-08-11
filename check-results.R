library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra) 
library(stringr)

# helper functions --------------------------------------------------------

#get comparison to baseline phenotype:
get_diff = function(.tbl, baseline_method = "SPACox", ...) {
  baseline = filter(.tbl, str_detect(method, baseline_method)) %>% pull(causal)
  .tbl %>% mutate(diff = causal - baseline, 
                  ratio = causal / baseline)
}

#Get plots for a given parameter setup
#all fucntions take piped data as input and outputs one panel each
base_diff_plot = function(data, dsword, cword, legend.pos = "none", baseline_method = baseline_method) {
  data %>% 
    filter(method != baseline_method) %>% 
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = diff, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 0, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw(base_size = base_text_size) +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = right_shift, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = paste0("Difference (to ", baseline_method, ")")) +
    scale_fill_viridis_d()
}
base_ratio_plot = function(data, dsword, cword, legend.pos = "none", baseline_method = baseline_method) {
  data %>% 
    filter(method != baseline_method) %>%
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = ratio, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw(base_size = base_text_size) +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = right_shift, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = "Power Ratio") +
    scale_fill_viridis_d()
}
base_null_plot = function(data, dsword, cword, legend.pos = "none") {
  data %>% 
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = mean_null_chisq, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw(base_size = base_text_size) + 
    ylim(.95, 1.05)  +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = expression("Average null" ~ chi^2)) +
    scale_fill_viridis_d()
}
base_power_plot = function(data, dsword, cword, legend.pos = "none") {
  data %>% 
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = power, fill = method)) +
    geom_violin() + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw(base_size = base_text_size)  +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = right_shift, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = "Power") +
    scale_fill_viridis_d()
}

#combines power, ratio and null chisq into one plot
plot_helper = function(.tbl, legend = legend_plot) {
  power_plot = filter(.tbl, type == "power") %>% pull(plots)
  ratio_plot = filter(.tbl, type == "ratio") %>% pull(plots)
  nullchisq_plot = filter(.tbl, type == "nullchisq") %>% pull(plots)
  
  plot_grid(power_plot[[1]] + theme(legend.position = "none",
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks = element_blank()),
            ratio_plot[[1]] + theme(legend.position = "none",
                                    axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks = element_blank()),
            nullchisq_plot[[1]]  + theme(legend.position = "none",
                                         axis.title.x = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.ticks = element_blank()),
            legend, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1,1,1,0.2))
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


# global settings ---------------------------------------------------------
right_shift = 0
base_text_size = 12
baseline_method = "SPACox"

legend_plot = cowplot::get_legend(base_power_plot(data = res %>%
                                                    filter(downsampling == "no",
                                                           C == "C1000"), cword = "C1000", legend.pos = "bottom"))

legend_plot_nospacox = cowplot::get_legend(base_power_plot(data = res %>%
                                                             filter(downsampling == "no",
                                                                    C == "C1000",
                                                                    method != "SPACox"), cword = "C1000", legend.pos = "bottom"))

# SPACox plots ------------------------------------------------------------

diff_plots_spacox =  res %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C,
                                       baseline_method = "SPACox"),
                             .f = ~ base_diff_plot(data = ..1,
                                                   dsword = ..2,
                                                   cword = ..3,
                                                   baseline_method = ..4)),
         type = "diff")

ratio_plots_spacox = res %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C,
                                       baseline_method = "SPACox"),
                             .f = ~ base_ratio_plot(data = ..1,
                                                    dsword = ..2,
                                                    cword = ..3,
                                                    baseline_method = ..4)),
         type = "ratio")

nullchisq_plots_spacox = res %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_null_plot(data = ..1,
                                                   dsword = ..2,
                                                   cword = ..3)),
         type = "nullchisq")

power_plots_spacox = res %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_power_plot(data = ..1,
                                                    dsword = ..2,
                                                    cword = ..3)),
         type = "power")



combined_spacox = bind_rows(diff_plots_spacox, 
                            ratio_plots_spacox,
                            nullchisq_plots_spacox,
                            power_plots_spacox) %>% 
  select(-data) %>% 
  filter(type != "diff") %>% 
  tidyr::nest(data = c(plots, type)) %>% 
  mutate(combined = purrr::pmap(.l = list(.tbl = data),
                                .f = ~ plot_helper(.tbl = ..1)))


# Saving SPACox plots -----------------------------------------------------


for(i in 1:nrow(combined_spacox)) {
  
  ggsave(paste0("./plots/simulations/SPACox-",
                combined_spacox$C[i], "-",
                if(combined_spacox$downsampling[i] == "yes" ) "ds-",
                "sim-plot.png"),
         combined_spacox$combined[[i]],
         dpi = 300,
         height = 4.5,
         width = 5.5)
}


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


# ltm plots ---------------------------------------------------------------

diff_plots_ltm =  ltm %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C,
                                       baseline_method = "SPACox"),
                             .f = ~ base_diff_plot(data = ..1,
                                                   dsword = ..2,
                                                   cword = ..3,
                                                   baseline_method = ..4)),
         type = "diff")

ratio_plots_ltm = ltm %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C,
                                       baseline_method = "SPACox"),
                             .f = ~ base_ratio_plot(data = ..1,
                                                    dsword = ..2,
                                                    cword = ..3,
                                                    baseline_method = ..4)),
         type = "ratio")
nullchisq_plots_ltm = ltm %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_null_plot(data = ..1,
                                                   dsword = ..2,
                                                   cword = ..3)),
         type = "nullchisq")

power_plots_ltm = ltm %>%
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_power_plot(data = ..1,
                                                    dsword = ..2,
                                                    cword = ..3)),
         type = "power")



combined_ltm = bind_rows(diff_plots_ltm, 
                         ratio_plots_ltm,
                         nullchisq_plots_ltm,
                         power_plots_ltm) %>% 
  select(-data) %>% 
  filter(type != "diff") %>% 
  tidyr::nest(data = c(plots, type)) %>% 
  mutate(combined = purrr::pmap(.l = list(.tbl = data),
                                .f = ~ plot_helper(.tbl = ..1)))

# Saving LTM based plots --------------------------------------------------


for (i in 1:nrow(combined_ltm)) {
  ggsave(paste0("./plots/simulations/LTM-",
                combined_ltm$C[i], "-",
                if (combined_ltm$downsampling[i] == "yes" ) "ds-",
                "sim-plot.png"),
         combined_ltm$combined[[i]],
         dpi = 300,
         height = 4.5,
         width = 5.5)
}




# combined data plots -----------------------------------------------------

combined = bind_rows(mutate(res, gen_mod = "Hazard"),
                     mutate(ltm, gen_mod = "Liability")) %>% 
  relocate(v, C, prev, method, downsampling, gen_mod)


base_ratio_plot2 = function(data, dsword, cword, legend.pos = "none", baseline_method = baseline_method) {
  data %>% 
    filter(method != baseline_method) %>%
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = ratio, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(gen_mod ~ prev, scale = "free", labeller = labeller(prev = facet1_names)) +
    theme_bw(base_size = base_text_size) +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = right_shift, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = "Power Ratio") +
    scale_fill_viridis_d()
}
base_null_plot2 = function(data, dsword, cword, legend.pos = "none") {
  data %>% 
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = mean_null_chisq, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(gen_mod ~ prev, scale = "free", labeller = labeller(prev = facet1_names)) +
    theme_bw(base_size = base_text_size) + 
    ylim(.95, 1.05)  +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = expression("Average null" ~ chi^2)) +
    scale_fill_viridis_d()
}


facet1_names = c("prev01"  = "1%",
  "prev025" = "2.5%",
  "prev05"  = "5%",
  "prev10"  = "10%",
  "prev20"  = "20%",
  "prev50"  = "50%")

base_power_plot2 = function(data, dsword, cword, legend.pos = "none") {
  data %>% 
    mutate(C = cword) %>% 
    ggplot(aes(x = method, y = power, fill = method)) +
    geom_violin() + 
    facet_grid(gen_mod ~ prev, scale = "free", labeller = labeller(prev = facet1_names)) + 
    theme_bw(base_size = base_text_size)  +
    theme(legend.position = legend.pos,
          axis.title.y = element_text(margin = margin(t = 0, r = right_shift, b = 0 , l = 0))) + 
    labs(fill = "Method:",
         y = "Power") +
    scale_fill_viridis_d()
}

#combines power, ratio and null chisq into one plot
plot_helper2 = function(.tbl, plot_type, legend_list = list(legend_plot, legend_plot_nospacox)) {
  p1 = .tbl$plots[[1]]
  p2 = .tbl$plots[[2]]
  if(plot_type == "ratio") {
    plot_grid(p1 + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks = element_blank()),
              p2 + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks = element_blank()),
              legend_list[[2]], ncol = 1, align = 'v', axis = 'l', rel_heights = c(1,1,0.1),
              labels = c("A)", "B)"))
    
  } else {
    plot_grid(p1 + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks = element_blank()),
              p2 + theme(legend.position = "none",
                         axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         axis.ticks = element_blank()),
              legend_list[[1]], ncol = 1, align = 'v', axis = 'l', rel_heights = c(1,1,0.1),
              labels = c("A)", "B)"))
    
  }
}

# legend_plot = cowplot::get_legend(base_power_plot(data = res %>%
#                                                     filter(downsampling == "no",
#                                                            C == "C1000"),
#                                                   cword = "C1000",
#                                                   legend.pos = "bottom"))


ratio_plots_combined = combined %>%
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(type = "ratio",
         plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C,
                                       baseline_method = "SPACox"),
                             .f = ~ base_ratio_plot2(data = ..1,
                                                     dsword = ..2,
                                                     cword = ..3,
                                                     baseline_method = ..4)))
nullchisq_plots_combined = combined %>%
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_null_plot2(data = ..1,
                                                   dsword = ..2,
                                                   cword = ..3)),
         type = "nullchisq")

power_plots_combined = combined %>%
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
  group_by(C, downsampling) %>%
  tidyr::nest() %>%
  ungroup() %>% 
  mutate(plots = purrr::pmap(.l = list(.tbl = data,
                                       dsword = downsampling,
                                       cword = C),
                             .f = ~ base_power_plot2(data = ..1,
                                                    dsword = ..2,
                                                    cword = ..3)),
         type = "power")


combined_combined = bind_rows(ratio_plots_combined,
                         nullchisq_plots_combined,
                         power_plots_combined) %>% 
  select(-data) %>% 
  filter(type != "diff") %>% 
  tidyr::nest(data = c(plots, downsampling)) %>% 
  mutate(combined = purrr::pmap(.l = list(.tbl = data, plot_type = type),
                                .f = ~ plot_helper2(.tbl = ..1, plot_type = ..2)))


#save all combined plots
for (i in 1:nrow(combined_combined)) {
  print(i)
  ggsave(filename = paste0("../plots/simulations/combined_", 
                           combined_combined$C[i],"_", combined_combined$type[i], ".png"),
         plot = combined_combined$combined[[i]],
         device = "png",
         dpi = 300,
         width = 5,
         height = 5)
}


# power ratio table -------------------------------------------------------

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

saveRDS(all_summary_info,
        "./results/simulation-summary-table.rds")
all_summary_info = readRDS("../results/simulation-summary-table.rds")

table_data = all_summary_info %>% 
  select(gen_mod, C, prev, downsampling, method, avg_power, avg_power_se) %>% 
  filter(C == "C250",
         downsampling == "yes") %>% 
  mutate(table_val = paste0(avg_power, "(", avg_power_se, ")")) %>% 
  as.data.frame()

table_data %>% select(-C, -downsampling) %>% 
  filter(method == "ADuLT")
table_data %>% select(-C, -downsampling) %>% 
  filter(method == "SPACox")
table_data %>% select(-C, -downsampling) %>% 
  filter(method == "Case-Control")
# gen_mod    prev       method avg_power avg_power_sd       table_val
# 1    Hazard  prev01 Case-Control    0.0920      0.00492  0.092(0.00492)
# 2    Hazard prev025 Case-Control    0.0816      0.00673 0.0816(0.00673)
# 3    Hazard  prev05 Case-Control    0.0600      0.00430    0.06(0.0043)
# 4    Hazard  prev20 Case-Control    0.0316      0.00384 0.0316(0.00384)
# 5 Liability  prev01 Case-Control    0.5124      0.00799 0.5124(0.00799)
# 6 Liability prev025 Case-Control    0.4720      0.01012  0.472(0.01012)
# 7 Liability  prev05 Case-Control    0.4288      0.00756 0.4288(0.00756)
# 8 Liability  prev20 Case-Control    0.3172      0.00757 0.3172(0.00757)
