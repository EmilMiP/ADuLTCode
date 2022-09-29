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
    facet_grid(C ~ prev, scale = "free",  labeller = labeller(prev = facet1_names)) +
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
    facet_grid(C ~ prev, scale = "free",  labeller = labeller(prev = facet1_names)) +
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
    facet_grid(C ~ prev, scale = "free",  labeller = labeller(prev = facet1_names)) +
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
    facet_grid(C ~ prev, scale = "free",  labeller = labeller(prev = facet1_names)) +
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

res_phmodel = readRDS("../results/phmodel_gwas_N1000k_M20k_results_all.rds") %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         color = "black",
         prev = ifelse(prev == "prev5", "prev05", prev))
res_phIPWmodel = readRDS("../results/phmodel_IPW_gwas_N1000k_M20k_results_all.rds") %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         color = "olive",
         prev = ifelse(prev == "prev5", "prev05", prev))

#load the summary info
res = readRDS("../results/cox_gwas_N1000k_M20k_results_all.rds") %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         method = ifelse(str_detect(method, "linear"), "Linear", method),
         color = ifelse(str_detect(method, "ADuLT"), "red", 
                        ifelse(str_detect(method, "Linear"), "green", "blue")),
         prev = ifelse(prev == "prev5", "prev05", prev)) %>% 
  bind_rows(res_phmodel) %>% 
  bind_rows(res_phIPWmodel) %>% 
  select(-beta_gen) %>% 
  group_by(v, C, prev, downsampling) %>% 
  group_map(.data = ., ~ get_diff(.tbl = .x), .keep = T) %>% 
  do.call("bind_rows", .) %>% 
  mutate(method = factor(method, levels = c("ADuLT", "SPACox", "Linear", "PH", "PH_wIPW"))) 



# global settings ---------------------------------------------------------
right_shift = 0
base_text_size = 12
baseline_method = "SPACox"

facet1_names = c(
  "prev01"  = "1%",
  "prev025" = "2.5%",
  "prev05"  = "5%",
  "prev10"  = "10%",
  "prev20"  = "20%",
  "prev50"  = "50%")

legend_plot = cowplot::get_legend(base_power_plot(data = res %>%
                                                    filter(downsampling == "yes",
                                                           C == "C1000"), cword = "C1000", legend.pos = "bottom"))



# SPACox plots ------------------------------------------------------------

diff_plots_spacox =  res %>%
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
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
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
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
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
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
  filter(str_detect(prev, paste0("prev", c("01", "025", "05", "20"), collapse = "|"))) %>% 
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
combined_spacox$combined[[2]]
combined_spacox$combined[[4]]
combined_spacox$data[[4]]$plots[[3]]


for(i in which(combined_spacox$downsampling == "yes")) {
  
  ggsave(paste0("../plots/simulations/SPACox-IPW-",
                combined_spacox$C[i], "-",
                if(combined_spacox$downsampling[i] == "yes" ) "ds-",
                "sim-plot.png"),
         combined_spacox$combined[[i]],
         dpi = 300,
         height = 4.5,
         width = 5.5)
}


