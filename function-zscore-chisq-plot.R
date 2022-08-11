library(dplyr)
library(ggplot2)

zscore_plot = function(data, 
                       suffix1,
                       suffix2,
                       y_label = "Z-score",
                       x_label = "Z-score",
                       pvalColName = "P_BOLT_LMM_INF",
                       SEColName = "SE",
                       BetaColName = "BETA",
                       text_size = 14
) {
  zscore_plot = data %>% 
    filter(!!as.symbol(paste0(pvalColName, suffix1)) < 5e-6 | !!as.symbol(paste0(pvalColName, suffix2)) < 5e-6) %>%
    mutate(x_sig = !!as.symbol(paste0(pvalColName, suffix1)) < 5e-8,
           y_sig = !!as.symbol(paste0(pvalColName, suffix2)) < 5e-8,
           color = x_sig + 2 * y_sig) %>%
    # mutate("z1" = !!as.symbol(paste0("BETA", suffix1)) / !!as.symbol(paste0("SE", suffix1)),
    #        "z2" = !!as.symbol(paste0("BETA", suffix2)) / !!as.symbol(paste0("SE", suffix2)))  %>%
    mutate("z1" = !!as.symbol(paste0(BetaColName, suffix1)),
           "z2" = !!as.symbol(paste0(BetaColName, suffix2))) %>% 
    ggplot(aes(x = z1, y = z2, color = as.factor(color))) +
    geom_point() +
    geom_hline(yintercept = c(-1, 1) * qnorm(5e-8, lower.tail = F), linetype = "longdash") +
    geom_vline(xintercept = c(-1, 1) * qnorm(5e-8, lower.tail = F), linetype = "longdash") +
    theme_minimal(base_size = text_size) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(aes(x = z1, y = z2),
                method = "lm", formula = "y ~ x + 0", inherit.aes = F) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = c("black", "blue", "red", "black"),
                       breaks = 0:3) +
    xlab(x_label) +
    ylab(y_label)
  return(zscore_plot)
}

chisq_plot = function(data, 
                      suffix1,
                      suffix2,
                      y_label = expression(chi^2 ~ "- Statistic"),
                      x_label = expression(chi^2 ~ "- Statistic"),
                      pvalColName = "P_BOLT_LMM_INF",
                      chisqColName = "CHISQ_BOLT_LMM_INF",
                      include_y = 30,
                      include_x = 30,
                      text_size = 14
) {
  zscore_plot = data %>% 
    filter(!!as.symbol(paste0(pvalColName, suffix1)) < 5e-6 | !!as.symbol(paste0(pvalColName, suffix2)) < 5e-6) %>%
    mutate(x_sig = !!as.symbol(paste0(pvalColName, suffix1)) < 5e-8,
           y_sig = !!as.symbol(paste0(pvalColName, suffix2)) < 5e-8,
           color = x_sig + 2 * y_sig) %>%
    ggplot(aes(x = !!as.symbol(paste0(chisqColName, suffix1)),
               y = !!as.symbol(paste0(chisqColName, suffix2)),
               color = as.factor(color))) + 
    geom_point() +
    geom_hline(yintercept = qchisq(5e-8, df = 1, lower.tail = F), linetype = "longdash") +
    geom_vline(xintercept = qchisq(5e-8, df = 1, lower.tail = F), linetype = "longdash") +
    theme_minimal(base_size = text_size) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(aes(x = !!as.symbol(paste0(chisqColName, suffix1)),
                    y = !!as.symbol(paste0(chisqColName, suffix2))), 
                method = "lm", formula = "y ~ x + 0", inherit.aes = F) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = c("black", "blue", "red", "black"),
                       breaks = 0:3) +
    xlab(x_label) +
    ylab(y_label) +
    expand_limits(y = include_y, x = include_x)
  return(zscore_plot)
}