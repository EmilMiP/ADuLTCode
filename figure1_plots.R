# library
library(ggridges)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)


ages = fread("U:/Projects/LT-FH++/Prevalences/RESULTS_cip.txt") %>% 
  as_tibble %>% 
  mutate(sex = ifelse(sex == "K", "F", "M"))


ages %>%
  filter(dx == "adhd",
         birth_year %% 5 == 0) %>%
  group_by(sex) %>%
  ggplot(aes(x = age, y = cip, color = sex)) +
  geom_line() +
  facet_wrap(~ birth_year)



#  ridge plot -------------------------------------------------------------



dashed_line_data = ages %>% 
  filter(dx == "adhd",
         birth_year == 2000,
         age == 10 | age == 15) %>% 
  mutate(x_end_left = 0,
         y_end_right = 0,
         thr = qnorm(cip, lower.tail = F)) %>% 
  arrange(desc(cip)) 



cip_labels = paste0("CIP[" , dashed_line_data$sex ,"](", dashed_line_data$age , ") == ", round(dashed_line_data$cip*100, 1), "~'%'") %>%
  sapply( function(x) parse( text = x))

panelA = ages %>% 
  filter(dx == "adhd",
         birth_year == 2000) %>% 
  group_by(sex) %>% 
  ggplot(aes(x = age, y = cip, color = sex)) +
  geom_line() + 
  theme_bw() +
  theme(
    axis.line         = element_line(color='black'),
    plot.background   = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none") +
  ggtitle("ADHD") +
  scale_color_manual(labels = c("F", "M"), values = c("red", "blue")) +
  ylab("Cumulative Incidence Proportion") +
  geom_segment(data = dashed_line_data, aes(x = age, y = y_end_right, xend = age, yend = cip), linetype = "dotted") +
  geom_segment(data = dashed_line_data, aes(x = age, y = cip, xend = x_end_left, yend = cip), linetype = "dotted") +
  annotate(geom = "text",
           x = 0,
           hjust  = c("inward","inward", "inward", "inward"),
           y = dashed_line_data$cip + 0.002,
           label = cip_labels, size = 3, parse = T) +
  scale_y_continuous(breaks = 0:3 / 100, labels = paste0(0:3, "%"))

ggsave(plot = panelA,
       filename = "../plots/panelA.png",
       dpi = 300,
       device = "png",
       height = 3,
       width = 3
)




vals = bind_rows(tibble(x = rnorm(10e6),
              sex = "M"),
              tibble(x = rnorm(10e6),
                     sex = "F")) %>% filter(abs(x) < 5)

cip_labels = paste0("CIP[" , dashed_line_data$sex ,"](", dashed_line_data$age , ") == ", round(dashed_line_data$cip, 4)) %>% 
  sapply( function(x) parse( text = x))


thr_labels = paste0("T[" , dashed_line_data$sex ,"](", dashed_line_data$age , ") == ", round(dashed_line_data$thr, 2)) %>% 
  sapply( function(x) parse( text = x))

vals_ridge = vals %>% mutate(case10 = ifelse(sex == "M" & x > dashed_line_data %>% filter(sex == "M", age == 10) %>% pull(thr),"caseM10",
                               ifelse(sex == "F" & x > dashed_line_data %>% filter(sex == "F", age == 10) %>% pull(thr), "caseF10",
                                       "ctrl")),
                      case15 = ifelse(sex == "M" & x > dashed_line_data %>% filter(sex == "M", age == 15) %>% pull(thr),"caseM15",
                                      ifelse(sex == "F" & x > dashed_line_data %>% filter(sex == "F", age == 15) %>% pull(thr), "caseF15",
                                             "ctrl")),
                      fill_color = ifelse(sex == "F", "red", "blue"))

vals_ridge %>% count(sex, case10, case15)
ridge_p = ggplot(vals_ridge, aes(x = x, y = sex)) +
#  geom_line(size = 1) +
  geom_density_ridges(fill = "white") +
  # geom_ribbon(data = shaded_vals, aes(x = x, ymax = y, fill = sex), ymin = 0, alpha = 0.5) + 
   theme_bw() +
  scale_color_manual("white") + 
  coord_cartesian(clip = "off")
# Build ggplot and extract data
d <- ggplot_build(ridge_p)$data[[1]]

color_dat = d %>% filter((group == 1 & x > dashed_line_data %>% filter(sex == "F", age == 10) %>% pull(thr)) | (group == 2 & x > dashed_line_data %>% filter(sex == "M", age == 10) %>% pull(thr))) %>% 
  mutate(fill = ifelse(group == 1, "blue", "red"))
color_dat2 = d %>% filter((group == 1 & x > dashed_line_data %>% filter(sex == "F", age == 15) %>% pull(thr)) | (group == 2 & x > dashed_line_data %>% filter(sex == "M", age == 15) %>% pull(thr))) %>% 
  mutate(fill = ifelse(group == 1, "blue", "red"))

line_dat = color_dat %>% 
  group_by(group) %>% 
  arrange(x) %>% 
  slice(1) %>% ungroup
line_dat2 = color_dat2 %>% 
  group_by(group) %>% 
  arrange(x) %>% 
  slice(1) %>% ungroup

text_coords = bind_rows(line_dat, line_dat2) %>% arrange(x)
ridge_p_colored = ridge_p +
  geom_ribbon(
    data = transform(color_dat, sex = group, fill = fill), #transform(subset(d, x >= 2), sex = group),
    aes(x, ymin = ymin, ymax = ymax, group = group, fill = fill),
    alpha = 0.3) + 
  geom_ribbon(
    data = transform(color_dat2, sex = group, fill = fill), #transform(subset(d, x >= 2), sex = group),
    aes(x, ymin = ymin, ymax = ymax, group = group, fill = fill),
    alpha = 0.3) + 
  theme(axis.line         = element_line(color='black'),
        plot.background   = element_blank(),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.border      = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = .9),
        legend.position = "none") +
   scale_fill_manual(labels = c("F", "M"), values = c("red", "blue")) +
  labs(x = "liability",
       y = "Sex",
       title = "Age-Dependent Liability Threshold Model") +
  geom_segment(data = transform(line_dat),
               aes(x = x, xend = x, y = ymin, yend = ymax),
               color = rev(line_dat$fill),
               linetype = "dotted") +
  geom_segment(data = transform(line_dat2),
               aes(x = x, xend = x, y = ymin, yend = ymax),
               color = rev(line_dat2$fill),
               linetype = "dotted") +
  annotate(geom = "text",
           x = text_coords$x,
           y = c(text_coords$ymax[1] + 0.15, text_coords$ymin[2] - 0.15, 
                 text_coords$ymax[3] + 0.15, text_coords$ymin[4] - 0.15) ,#0.01 + c(0, dnorm(dashed_line_data$thr[2])),
           hjust = c("outward", "outward", "outward", "outward"),
           label = thr_labels)  

  
ggsave(filename = "../plots/Figure1_ridge_density.png",
       plot = ridge_p_colored,
       device = "png",
       width = 5,
       height = 3.5,
       dpi = 300)
