library(data.table)
library(dplyr)
library(stringr)

#function used in SPAcox to get time of event
get_onset = function(lambda, eta) {
  nsib = length(eta)
  U = runif(nsib)
  #onset
  lambda * sqrt( -log(U) / exp(eta) )
}

# Find Matching Event Rates -----------------------------------------------
N = 1000e3
lambda = 1
find_prev_and_ER = function(lambda, N = 1000e3, target_prev) {
  obs_prev = tibble(liab = rnorm(N)) %>% 
    mutate(onset = get_onset(lambda = lambda, liab),
           censor = rweibull(n = N, shape = 1, scale = 0.15), #censoring time || Independent
           status = onset <= censor + 0) %>% 
    pull(status) %>% mean()
  target_prev - obs_prev
}
prevs = c(0.0101, 0.025, 0.05, 0.1, 0.2, 0.5)
est = lapply(prevs, function(prev) {
  uniroot(find_prev_and_ER, c(0,3), target_prev = prev)
})


#   -----------------------------------------------------------------------

event_rates = sapply(est, function(x) x$root)


all_files = list.files("./simulatedData", full.names = T, pattern = "true") %>% 
  str_subset("normal") %>% 
  str_subset("SPACox") %>%
  str_subset("N1000k-M20")

#constructing prevalence name to keep file names ordered
prev_text = paste0("prev",round(prevs * 100,1))
prev_text[1] = "prev01"
prev_text[2] = "prev025"


lapply(seq_along(all_files), function(i) {
  cur_file = all_files[i]
  cur_true = fread(cur_file) %>% 
    as_tibble
  
  lapply(seq_along(event_rates), function(ii) {
    lambda = event_rates[ii]
    new_true = cur_true %>%
      mutate(onset = get_onset(lambda = lambda, gen + env),
             censor = rweibull(n = nrow(cur_true), shape = 1, scale = 0.15), #censoring time || Independent
             status = onset <= censor + 0) 
    
    fwrite(new_true,
           str_replace(cur_file, "-v", paste0("-", prev_text[ii], "-v")))
    
    NULL
  })
  NULL
})





all_files = list.files("./simulatedData", full.names = T, "SPACox") %>% 
  str_subset(., "prev") %>% 
  str_subset("N1000k-M20k")

tmp2 = tibble(ph = all_files) %>%   
  mutate(gen_beta = as.vector(str_match(ph, "unif|normal")),
         prev_txt = as.vector(str_match(ph, paste0(prev_text, collapse = "|", sep = "-"))),
         C = as.vector(str_match(ph, "C250|C1000"))) %>% 
  group_by(gen_beta, C, prev_txt) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup %>% 
  mutate(prev = lapply(sapply(grps, function(x) x[1]), function(cur_phen) {
  fread(cur_phen) %>%
    as_tibble() %>% 
    summarise(prev = mean(status)) %>% pull(prev)
})) %>% 
  mutate(prev = unlist(prev)) %>%
  select(-grps)

saveRDS(tmp2,
        "./results/ER-to-prev-conversion.rds")


ph = tmp2 %>% 
  mutate(ages = lapply(sapply(tmp2$grps, function(x) x[1]), function(cur_phen) {
    fread(cur_phen) %>% 
      as_tibble() %>% 
      select(age)
  }))

lapply(ph$ages, summary)

sum(is.infinite(ph$ages[[1]]$age))

ggplot(ph$ages[[1]], aes(x = age)) +
  geom_histogram()
ph$ages[[1]] %>% filter(age > 100)
# old ---------------------------------------------------------------------



## Changing names of true files, if not event rate is in file name

true_files = list.files("./simulatedData", full.names = T, "true") %>% 
  str_subset(., "ER|LTM|unif|maf", negate = T)

for (i in seq_along(true_files)) {
  cur_true = fread(true_files[i])
  
  fwrite(cur_true,
         str_replace(true_files[i], "-v", paste0("-ER15-v")))
}




# prevs = lapply(event_rates, function(lambda) {
#   cur_true %>% 
#     mutate(onset2 = get_onset(lambda = lambda, gen + env),
#            censor2 = rweibull(n = nrow(cur_true), shape = 1, scale = 0.15), #censoring time || Independent
#            status2 = onset2 <= censor2) %>% 
#     summarise(prev = mean(status2))
# }) %>% 
#   do.call("rbind", .) %>% 
#   mutate(lambda = event_rates) 
# ggplot(prevs, aes(x = lambda, prev)) +
#   geom_point()

# 
# prevs2 = lapply(seq(0,1,0.01), function(lambda) {
#   cur_true %>% 
#     mutate(onset2 = get_onset(lambda = .15, gen + env),
#            censor2 = rweibull(n = nrow(cur_true), shape = 1, scale = lambda), #censoring time || Independent
#            status2 = onset2 <= censor2) %>% 
#     summarise(prev = mean(status2))
# }) %>% 
#   do.call("rbind", .) %>% 
#   mutate(lambda = seq(0,1,0.01)) 
# 
# ggplot(prevs2, aes(x = lambda, prev)) +
#   geom_point()
