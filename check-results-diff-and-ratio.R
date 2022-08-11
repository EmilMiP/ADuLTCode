library(dplyr)
library(ggplot2)

#helper function
get_diff = function(.tbl, baseline_method = "SPACox", ...) {
  baseline = filter(.tbl, str_detect(method, baseline_method)) %>% pull(causal)
  .tbl %>% mutate(diff = causal - baseline, 
                  ratio = causal / baseline)
}

#load the summary info
res = readRDS("./results/cox_gwas_results_all.rds") %>% 
  filter(str_detect(method, negate = T, "logistic")) %>% 
  mutate(power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         method = ifelse(str_detect(method, "linear"), "Linear", method),
         color = ifelse(str_detect(method, "ADuLT"), "red", 
                        ifelse(str_detect(method, "Linear"), "green", "blue"))) %>% 
  select(-beta_gen) %>% group_by(v, C, event_rate, downsampling) %>% 
  group_map(.data = ., ~ get_diff(.tbl = .x), .keep = T) %>% 
  do.call("bind_rows", .) %>% 
  mutate(method = factor(method, levels = c("ADuLT", "SPACox", "Linear")))

baseline_method = "SPACox"

# Diff --------------------------------------------------------------------
base_diff_plot = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = diff, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 0, lty = "longdash") + 
    facet_grid(C ~ event_rate, scale = "free") +
    theme_bw() +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = paste0("Difference (to ", baseline_method, ")")) +
    scale_fill_viridis_d()
}
  
  
baseline_spacox_diff_nds = base_diff_plot(data = res %>%
                                            filter(downsampling == "no",
                                                   method != baseline_method))


baseline_spacox_diff_ds = base_diff_plot(data = res %>% 
                                           filter(downsampling == "yes",
                                                  method != baseline_method))

# Ratio -------------------------------------------------------------------
base_ratio_plot = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = ratio, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ event_rate, scale = "free") +
    theme_bw() +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = paste0("Power Ratio (to ", baseline_method, ")")) +
    scale_fill_viridis_d()
}

baseline_spacox_ratio_nds = res %>% 
  filter(downsampling == "no",
         method != baseline_method) %>% 
  base_ratio_plot(data = .)

baseline_spacox_ratio_ds = res %>% 
  filter(downsampling == "yes",
         method != baseline_method) %>% 
 base_ratio_plot(data = .)


# Null Chisq --------------------------------------------------------------
base_null_plot = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = mean_null_chisq, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ event_rate, scale = "free") +
    theme_bw() + 
    ylim(.95, 1.05)  +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = expression("Average Null" ~ chi^2)) +
    scale_fill_viridis_d()
}

baseline_spacox_null_nds = res %>% 
  filter(downsampling == "no") %>% 
  base_null_plot()


baseline_spacox_null_ds = res %>% 
  filter(downsampling == "yes") %>% 
  base_null_plot()


# Power -------------------------------------------------------------------
base_power_plot = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = power, fill = method)) +
    geom_violin() + 
    facet_grid(C ~ event_rate, scale = "free") +
    theme_bw()  +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = "Power") +
    scale_fill_viridis_d()
}

baseline_spacox_power_nds = res %>% 
  filter(downsampling == "no") %>% 
  base_power_plot()


baseline_spacox_power_ds = res %>% 
  filter(downsampling == "yes") %>% 
  base_power_plot()


# Combining Plots ---------------------------------------------------------
library(cowplot)

library(grid)
library(gridExtra) 
legend = cowplot::get_legend(base_power_plot(data = res %>% 
                                               filter(downsampling == "no"), legend.pos = "bottom"))


#no downsampling
spacox_nds = plot_grid(baseline_spacox_power_nds + theme(legend.position = "none",
                                                         axis.title.x = element_blank(),
                                                         axis.text.x = element_blank(),
                                                         axis.ticks = element_blank()),
                       baseline_spacox_ratio_nds + theme(legend.position = "none",
                                                         axis.title.x = element_blank(),
                                                         axis.text.x = element_blank(),
                                                         axis.ticks = element_blank()),
                       baseline_spacox_null_nds  + theme(legend.position = "none",
                                                         axis.title.x = element_blank(),
                                                         axis.text.x = element_blank(),
                                                         axis.ticks = element_blank()),
                       legend, ncol = 1, rel_heights = c(1,1,1,0.2))
ggsave("./plots/simulations/SPACox-potential-sim-plot.png",
       spacox_nds,
       dpi = 300,
       height = 7,
       width = 5.5)

#downsampling
spacox_ds = plot_grid(baseline_spacox_power_ds + theme(legend.position = "none",
                                                       axis.title.x = element_blank(),
                                                       axis.text.x = element_blank(),
                                                       axis.ticks = element_blank()),
                      baseline_spacox_ratio_ds + theme(legend.position = "none",
                                                         axis.title.x = element_blank(),
                                                         axis.text.x = element_blank(),
                                                         axis.ticks = element_blank()),
                      baseline_spacox_null_ds  + theme(legend.position = "none",
                                                       axis.title.x = element_blank(),
                                                       axis.text.x = element_blank(),
                                                       axis.ticks = element_blank()), 
                      legend, ncol = 1, rel_heights = c(1,1,1,0.2))

ggsave("./plots/simulations/SPACox-downsampling-potential-sim-plot.png",
       spacox_ds,
       dpi = 300,
       height = 7,
       width = 5.5)


#  ltm --------------------------------------------------------------------
#load the summary info
ltm = readRDS("./results/ltm_gwas_results_all.rds") %>% 
  filter(str_detect(method, negate = T, "logistic")) %>% 
  mutate(prev = ifelse(prev == "prev5", "prev05", prev),
         power = ifelse(str_detect(C, "1000"), causal/1000, causal/250),
         method = ifelse(str_detect(method, "linear"), "Linear", method),
         color = ifelse(str_detect(method, "ADuLT"), "red", 
                        ifelse(str_detect(method, "Linear"), "green", "blue"))) %>% 
  group_by(v, C, prev, downsampling) %>% 
  group_map(.data = ., ~ get_diff(.tbl = .x, baseline_method = "SPACox"), .keep = T) %>% 
  do.call("bind_rows", .) %>% 
  mutate(method = factor(method, levels = c("ADuLT", "SPACox", "Linear")))

# diff --------------------------------------------------------------------
base_diff_plot_ltm = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = diff, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 0, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw() +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = paste0("Difference (to ", baseline_method, ")")) +
    scale_fill_viridis_d()
}

baseline_ltm_diff_nds = ltm %>% 
  filter(downsampling == "no",
         method != "SPACox") %>% 
  base_diff_plot_ltm()
 

baseline_ltm_diff_ds = ltm %>% 
  filter(downsampling == "yes",
         method != "SPACox") %>% 
  base_diff_plot_ltm()

# Ratio -------------------------------------------------------------------

base_ratio_plot_ltm = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = ratio, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw() +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = paste0("Power Ratio (to ", baseline_method, ")")) +
    scale_fill_viridis_d()
}

baseline_ltm_ratio_nds = ltm %>% 
  filter(downsampling == "no",
         method != "SPACox") %>% 
  base_ratio_plot_ltm()


baseline_ltm_ratio_ds = ltm %>% 
  filter(downsampling == "yes",
         method != "SPACox") %>% 
  base_ratio_plot_ltm()


# null Chisq --------------------------------------------------------------

base_null_plot_ltm = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = mean_null_chisq, fill = method)) +
    geom_violin() + 
    geom_hline(yintercept = 1, lty = "longdash") + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw() + 
    ylim(.95, 1.05)  +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = expression("Average Null" ~ chi^2)) +
    scale_fill_viridis_d()
}
baseline_ltm_null_nds = ltm %>% 
  filter(downsampling == "no") %>% 
  base_null_plot_ltm()

baseline_ltm_null_ds = ltm %>% 
  filter(downsampling == "yes") %>% 
  base_null_plot_ltm() 

# power -------------------------------------------------------------------
base_power_plot_ltm = function(data, legend.pos = "none") {
  data %>% 
    ggplot(aes(x = method, y = power, fill = method)) +
    geom_violin() + 
    facet_grid(C ~ prev, scale = "free") +
    theme_bw()  +
    theme(legend.position = legend.pos) + 
    labs(fill = "Method",
         y = "Power") +
    scale_fill_viridis_d()
}

baseline_ltm_power_nds = ltm %>% 
  filter(downsampling == "no") %>% 
  base_power_plot_ltm()

baseline_ltm_power_ds = ltm %>% 
  filter(downsampling == "yes") %>% 
  base_power_plot_ltm()



# Combining plots ---------------------------------------------------------

#no downsampling
ltm_nds = plot_grid(baseline_ltm_power_nds + theme(legend.position = "none",
                                                   axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks = element_blank()),
                    baseline_ltm_ratio_nds + theme(legend.position = "none",
                                                      axis.title.x = element_blank(),
                                                      axis.text.x = element_blank(),
                                                      axis.ticks = element_blank()),
                    baseline_ltm_null_nds  + theme(legend.position = "none",
                                                   axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.ticks = element_blank()),
                    legend, ncol = 1, rel_heights = c(1,1,1,0.2))

ggsave("./plots/simulations/LTM-potential-sim-plot.png",
       ltm_nds,
       dpi = 300,
       height = 7,
       width = 5.5)

#downsampling
ltm_ds = plot_grid(baseline_ltm_power_ds + theme(legend.position = "none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank(),
                                                 axis.ticks = element_blank()),
                   baseline_ltm_ratio_ds + theme(legend.position = "none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank(),
                                                 axis.ticks = element_blank()),
                   baseline_ltm_null_ds  + theme(legend.position = "none",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank(),
                                                 axis.ticks = element_blank()), 
                   legend, ncol = 1, rel_heights = c(1,1,1,0.2))

ggsave("./plots/simulations/LTM-downsampling-potential-sim-plot.png",
       ltm_ds,
       dpi = 300,
       height = 7,
       width = 5.5)

