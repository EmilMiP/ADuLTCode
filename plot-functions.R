#packages
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
library(cowplot)
library(PlotData)


zscore_plot = function(data,
                       suffix1,
                       suffix2,
                       y_label = "Z-score",
                       x_label = "Z-score",
                       pvalColName = "P_BOLT_LMM_INF",
                       SEColName = "SE",
                       BetaColName = "BETA",
                       pretty_break_x = 10,
                       pretty_break_y = 10,
                       filter_sugg = F
) {
  
  if (filter_sugg) {
    zscore_plot = data %>%
      filter(!!as.symbol(paste0(pvalColName, "_", suffix1)) < 5e-6 | !!as.symbol(paste0(pvalColName, "_", suffix2)) < 5e-6)
  } else {
    zscore_plot = data
  }
  
  zscore_plot = zscore_plot %>%
    mutate(x_sig = !!as.symbol(paste0(pvalColName, "_", suffix1)) < 5e-8,
           y_sig = !!as.symbol(paste0(pvalColName, "_", suffix2)) < 5e-8,
           color = x_sig + 2 * y_sig) %>%
    mutate("z1" = !!as.symbol(paste0(BetaColName, "_", suffix1)) / !!as.symbol(paste0(SEColName, "_", suffix1)),
           "z2" = !!as.symbol(paste0(BetaColName, "_", suffix2)) / !!as.symbol(paste0(SEColName, "_", suffix2)))  %>%
    ggplot(aes(x = z1, y = z2, color = as.factor(color))) +
    geom_point() +
    geom_hline(yintercept = c(-1, 1) * qnorm(5e-8, lower.tail = F), linetype = "longdash") +
    geom_vline(xintercept = c(-1, 1) * qnorm(5e-8, lower.tail = F), linetype = "longdash") +
    theme_minimal(base_size = 22) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(aes(x = z1, y = z2),
                method = "lm", formula = "y ~ x + 0", inherit.aes = F) +
    theme(legend.position = "none",
          axis.title = element_text(size = 20)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = pretty_break_x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = pretty_break_y)) +
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
                      pretty_break_x = 10,
                      pretty_break_y = 10,
                      filter_sugg = T
) {
  
  if (filter_sugg) {
    zscore_plot = data %>%
      filter(!!as.symbol(paste0(pvalColName, "_", suffix1)) < 5e-6 | !!as.symbol(paste0(pvalColName, "_", suffix2)) < 5e-6)
  } else {
    zscore_plot = data
  }
  
  zscore_plot = zscore_plot %>%
    mutate(x_sig = !!as.symbol(paste0(pvalColName, "_", suffix1)) < 5e-8,
           y_sig = !!as.symbol(paste0(pvalColName, "_", suffix2)) < 5e-8,
           color = x_sig + 2 * y_sig) %>%
    ggplot(aes(x = !!as.symbol(paste0(chisqColName, "_", suffix1)),
               y = !!as.symbol(paste0(chisqColName, "_", suffix2)),
               color = as.factor(color))) +
    geom_point() +
    geom_hline(yintercept = qchisq(5e-8, df = 1, lower.tail = F), linetype = "longdash") +
    geom_vline(xintercept = qchisq(5e-8, df = 1, lower.tail = F), linetype = "longdash") +
    theme_minimal(base_size = 22) +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(aes(x = !!as.symbol(paste0(chisqColName, "_", suffix1)),
                    y = !!as.symbol(paste0(chisqColName, "_", suffix2))),
                method = "lm", formula = "y ~ x + 0", inherit.aes = F) +
    theme(legend.position = "none",
          axis.title = element_text(size = 20)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = pretty_break_x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = pretty_break_y)) +
    scale_color_manual(values = c("black", "blue", "red", "black"),
                       breaks = 0:3) +
    xlab(x_label) +
    ylab(y_label)
  return(zscore_plot)
}