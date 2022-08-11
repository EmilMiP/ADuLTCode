library(dplyr)
library(scales)
library(ggplot2)


# combining information ---------------------------------------------------

itrs = expand.grid(list(itr = 1:10,
                        ncores = c(1, 2, 4, 8, 12, 16, 24, 32),
                        N = c(10e3, 50e3, 100e3)))

adult = readRDS("../results/computation_times.rds") %>% 
  sapply(., function(x) x[3]) %>% 
  as_tibble() %>% 
  bind_cols(as_tibble(itrs), times = .) %>% 
  mutate(method = "ADuLT",
         M = list(c(100e3, 250e3, 500e3, 1000e3))) %>%
  tidyr::unnest(M) %>% 
  rename(times = value) %>% 
  relocate(itr, ncores, N, M, method, times) %>% 
  filter(ncores > 3, ncores < 17)

linreg = readRDS("../results/computation-times-N100k-M1000k-linreg-raw.rds") %>% 
  mutate(time_linreg = sapply(times, function(x) x[3])) %>% 
  select(-times)

spacox = readRDS("../results/computation-times-N100k-M1000k-SPACox-raw.rds") %>% 
  mutate(method = "SPACox",
         ncores = 1,
         times = sapply(times, function(x) x[3])) %>% 
  rename(itr = v) %>% 
  relocate(itr, ncores, N, M, method, times)

linreg %>% 
  group_by(ncores, N, M) %>% 
  summarise(mean = mean(time_linreg),
            median = median(time_linreg),
            se = sd(time_linreg)/sqrt(10)) %>%
  ungroup() 

combined_times = adult  %>% 
  left_join(., linreg, by = c("itr", "ncores", "N", "M")) %>% 
  mutate(times = times + time_linreg,
         method = "ADuLT+linreg") %>% 
  select(-method.x, -method.y)


mean_times = bind_rows(combined_times, spacox) %>% 
  group_by(ncores, N, M, method) %>% 
  summarise(mean_times = mean(times/60),
            se_times = sd(times/60)/sqrt(n()),
            .groups = "drop")



plot_mean_times = mean_times %>% 
  filter(ncores < 5) %>% # gets ncores == 1 for spacox and ncores == 4 for linreg + adult
  mutate(mean_times = mean_times,
         se_times   = se_times)
#labeller for M when faceting
M_labeller_list = as_labeller(c("1e+05"  = "M = 100e3",
                                "250000" = "M = 250e3",
                                "5e+05"  = "M = 500e3",
                                "1e+06"  = "M = 1000e3"))


p1 = ggplot(plot_mean_times, aes(x = N, y = mean_times, color = as.factor(method))) +
  theme_bw(base_size = 18) +
  geom_point(lwd = 1.05) +
  geom_line(lwd = 1.025) +
  labs(x = "Number of Individuals",
       y = "Time (minutes)",
       color = "Method") +
  theme(plot.title = element_text(hjust = .5, size = 20),
        panel.spacing = unit(1.1, "lines"),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = .75)) + 
  scale_x_continuous(breaks = unique(plot_mean_times$N), labels = paste0(unique(plot_mean_times$N)/1000, "k")) + 
  scale_y_continuous(breaks = seq(0, 90, by = 15)) +
  geom_errorbar(data = plot_mean_times, 
                aes(ymin = mean_times - 1.96 * se_times, ymax = mean_times + 1.96 * se_times, width = 2500)) + 
  facet_grid( ~ M, labeller = M_labeller_list)

# ggsave(filename = "./plots/computation-times.pdf",
#        plot = p1,
#        device = "pdf",
#        dpi = 300,
#        height = 4,
#        width = 8)  
ggsave(filename = "../plots/computation-times.png",
       plot = p1,
       device = "png",
       dpi = 500,
       height = 4,
       width = 9)  


saveRDS(plot_mean_times,
        "../tables/computation-times.rds")
