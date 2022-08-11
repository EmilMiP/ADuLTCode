library(ggplot2)
library(dplyr)
library(data.table)


manhattanPlot_noSave = function(data,
                                title,
                                sugg = 1e-6,
                                sigg = 5e-8,
                                pvalColName = "P_BOLT_LMM_INF",
                                basePairName = "BP",
                                chrName = "CHR",
                                SNP_distance = 250,
                                fixed_y = NA,
                                theme_base_size = 12,
                                x_text_size = theme_base_size,
                                ld_clumped_indx = NULL) {
  data[["BP"]] = as.double(data[[basePairName]]) #needed to avoid integer overflow
  data[["CHR"]] = as.integer(data[[chrName]])
  data[[pvalColName]] = as.numeric(data[[pvalColName]])
  
  data[[pvalColName]] = ifelse(is.na(data[[pvalColName]]), 1, data[[pvalColName]])
  if (!is.null(ld_clumped_indx)) {
    ld_snps = data %>% slice(ld_clumped_indx)
  }
  data = data %>% mutate(rank = rank(!!as.symbol(pvalColName))) %>% filter(rank < 50e3)
  #if(sum(data$'A1FREQ' < 0.01) > 0) {
  # ok_maf = data$'A1FREQ' >= .01
  #  data = data[ok_maf,]
  #}
  
  if (sum(data[[pvalColName]] == 1) != 0) {
    cat("Removing", sum(data[[pvalColName]] == 1), "SNPs with a p-value of 1. \n")
    data = data[!(data[[pvalColName]] == 1), ]
  }
  
  if (sum(data[[pvalColName]] == 0) != 0) {
    cat("Removing", sum(data[[pvalColName]] == 0), "SNPs with a p-value of 0. \n")
    data = data[!(data[[pvalColName]] == 0),]
  }
  
  if (any(data$CHR > 22)) {
    cat("More than 22 chromosomes detected, removing chromosome 23 and above. \n")
    data = filter(data, CHR <= 22)
  }
  #this data will form the basis for extrating the top snps and the manhattan plot
  dataplyr = data %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len = max(BP)) %>%
    dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    dplyr::left_join(data, ., by = c("CHR" = "CHR")) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum = BP + tot)
  
  
  if (sigg != 5e-8) {
    cat("Using Bonferoni corrected alpha\n")
    bonferoni_alpha = 0.05 / nrow(data)
  } else {
    cat("Using fixed significance level: 5e-8 \n")
    bonferoni_alpha = sigg
  }
  
  if (bonferoni_alpha > sugg) {
    cat("Bonferoni corrected significance level *higher* than suggested significance level!\n Bonferioni:", bonferoni_alpha, "\n suggeted:", sugg, "\n Terminating.." )
    stop()
  }
  
  cat("Working on Manhattan plot. \n")
  
  
  
  
  dataplyr[["sig_P"]]   = ifelse(dataplyr[[pvalColName]] < bonferoni_alpha, T, F)
  dataplyr[["sugg_P"]]  = ifelse((bonferoni_alpha < dataplyr[[pvalColName]]) & (dataplyr[[pvalColName]] < sugg), T, F)
  dataplyr[["non_sig"]] = ifelse((dataplyr[[pvalColName]] >= sugg), T, F)
  
  dataplyr[["aboveunif"]] = runif(n = nrow(dataplyr))
  
  dataplyr = filter(dataplyr, !!as.symbol(pvalColName) < 0.05)  #removing any snps that are not nominally significant
  
  axisdataplyr = dataplyr %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarize(center = (max(BPcum) + min(BPcum))/ 2)
  
  non_sig_data = subset(dataplyr, non_sig == T)
  sig_P_data = subset(dataplyr, sig_P == T)
  sugg_P_data = subset(dataplyr, sugg_P == T)
  
  if (!is.na(fixed_y)) {
    ylims = c(min(-log10(dataplyr[[pvalColName]])), fixed_y)
  } else {
    ymax = max(-log10(dataplyr[[pvalColName]]) + 2, -log10(bonferoni_alpha) + 2)
    ylims = c(min(-log10(dataplyr[[pvalColName]])), ymax)
  }
  
  
  column = ensym(pvalColName)
  # find top snps to mark separately
  cat("Distance used:", SNP_distance, "Kb \n")
  distance = SNP_distance * 1000
  
  column = as.name(pvalColName)
  
  sigg_store = data.frame()
  sigg_snps = dataplyr %>% dplyr::filter(bonferoni_alpha > !!column)
  
  while(nrow(sigg_snps) > 0) {
    cur_sigg_snp = sigg_snps[which.min(sigg_snps[[pvalColName]]),]
    sigg_store   = rbind(sigg_store, cur_sigg_snp)
    sigg_snps    = sigg_snps[abs(cur_sigg_snp$BPcum - sigg_snps$BPcum) > distance,]
  }
  
  
  
  manPlot = ggplot2::ggplot(dataplyr, aes(x = BPcum, y = -log10(!!column)) ) +
    ggplot2::geom_point(data = sugg_P_data, color = "orange", size = 2, alpha = .5) +
    ggplot2::geom_point(data = sig_P_data , color = hsv(0, s = .4, v = 1), size = 2, alpha = .7) +
    ggplot2::geom_point(data = non_sig_data, size = 2, alpha = .5, aes(color = as.factor(CHR))) +
    
    ggplot2::scale_color_manual(values = rep(c("gray", "darkgray"), 11)) +
    ggplot2::scale_x_continuous(labels = axisdataplyr$CHR, breaks = axisdataplyr$center, guide = guide_axis(check.overlap = TRUE)) +
    ggplot2::scale_y_continuous(limits = ylims) +
    
    ggplot2::ggtitle(title) +
    ggplot2::labs(x = "Chromosome", y = "-log10(P)") +
    
    ggplot2::geom_hline(yintercept = -log10(bonferoni_alpha)) +
    ggplot2::geom_hline(yintercept = -log10(sugg), linetype = "dashed") +
    
    ggplot2::theme_minimal(base_size = theme_base_size) +
    theme(
      plot.title = element_text(hjust = .5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = x_text_size)
    )
  if (!is.null(ld_clumped_indx)){
    cat("Assigning diamonds to genome-wide significant LD clumped SNPs. \n")
    ld_snps = left_join(ld_snps, dataplyr) %>% filter(-log10(!!column) > -log10(sigg)) 
    manPlot = manPlot + ggplot2::geom_point(data = ld_snps , fill = "red", size = 2, shape = 23)
  } else {
    if (nrow(sigg_store > 0)) {
      cat("Assigning diamonds to lowest p-value genome-wide significant SNP in window.\n")
      manPlot = manPlot + ggplot2::geom_point(data = sigg_store , fill = "red", size = 2, shape = 23)
    }
  }

    return(list("manhattanPlot" = manPlot))

}
