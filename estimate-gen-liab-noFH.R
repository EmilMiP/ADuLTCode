library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(future.batchtools)
library(LTFHPlus)

h2 = .5
NCORES = 16
all_files = list.files("./simulatedData", "SPA", full.names = T) %>% 
  str_subset("true") %>% 
  str_subset("N1000k-M20k") %>% 
  str_subset("prev")

plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


future.apply::future_lapply(seq_along(all_files), function(ctr) {
  cur_true = all_files[ctr]
  true = fread(cur_true) %>% 
    as_tibble %>% 
    mutate(event = pmin(onset, censor))

  
  phen = true %>% 
    arrange(event) %>% 
    mutate(prop  = pmin((cumsum(status) + 1)/n(), mean(status)),
           est   = qnorm(prop, lower.tail = F),
           age   = LTFHPlus::cir_to_age(cir = prop, pop_prev = mean(status), slope = 0.2, age_mid = 50),
           lower = ifelse(status, LTFHPlus::age_to_thres(age = age, pop_prev = mean(status), age_mid = 50, slope = 0.2), -Inf),
           upper = ifelse(status, LTFHPlus::age_to_thres(age = age, pop_prev = mean(status), age_mid = 50, slope = 0.2), #alternatively, set upper bound to Inf for cases
                          LTFHPlus::age_to_thres(age = age, pop_prev = mean(status), age_mid = 50, slope = 0.2)),
           status = status + 0L) %>%
    rename(ids = FID) %>% 
    relocate(ids, lower, upper)
  
  library(future.apply)
  plan(tweak(multisession, workers = NCORES))
  noFH = LTFHPlus::estimate_gen_liability_noFH(phen = phen, h2 = h2) %>% 
    arrange(ids)
  
  fwrite(noFH, 
         str_replace(cur_true, "true", "phen"),
         quote = F,
         na = "NA",
         sep = " ")
}, future.seed = T)

# old ---------------------------------------------------------------------

phen %>% ggplot(aes(x = event, y = prop)) + 
  geom_point() +
  geom_hline(yintercept = mean(phen$status))



phen %>% ggplot(aes(x = prop, y = est)) + 
  geom_point() 

phen %>% ggplot(aes(x = prop, y = age)) + 
  geom_point() 


#why do I get the NAs for some values, and not all?
#why only under some parameter setups, and not all?
pop_prev = .1
age_mid = 60
slope = 1/8
liab = 1.06
age_mid - log(pop_prev/stats::pnorm(liab, lower.tail = F) - 1)* 1/slope

filter(phen, is.na(age)) %>% 
  pull(est) %>% hist
qnorm(mean(phen$status), lower.tail = F)


all_files = list.files("./simulatedData", "phen", full.names = T) 


phenotypes = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "ER")
  tibble(full_names = ph,
         v = v) %>% 
    mutate(maf = as.vector(str_match(full_names, "maf[1-4]0")),
           gen_beta = as.vector(str_match(full_names, "unif|normal")),
           event_rate = as.vector(str_match(full_names, "ER[1-8]0|ER1"))) %>% 
    group_by(v, maf, gen_beta, event_rate) %>% 
    summarise(true = full_names, .groups = "drop") %>% 
    ungroup
}) %>% 
  do.call("rbind", .)

noFH %>% 
  ggplot(aes(x = post_gen_no_fam, y = gen)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm")


noFH %>% 
  ggplot(aes(x = event, y = prop)) +
  geom_point() 

phen = fread(str_replace(cur_true, "true", "phen")) %>% 
  as_tibble()
hist(phen$upper)
ggplot(phen, aes(y = gen, x = post_gen_no_fam))+
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm")

cor(phen$gen, phen$post_gen_no_fam)
cor(phen$gen, phen$status)
cor(phen$status, phen$post_gen_no_fam)


onset_cdf = ecdf(true$onset)
event_cdf = ecdf(tmp$event)

hist(rweibull(n = 100e3, shape = 1, scale = 0.15))
hist(true$onset)
plot(sort(rweibull(n = 100e3, shape = 1, scale = 0.15)), sort(true$onset))
abline(a = 0, b = 1)

tmp = phen %>% 
  left_join(true) %>%
  arrange(onset) %>% 
  mutate(prop2 = pweibull(onset, shape = 1, scale = 0.15),
         prop3 = pweibull(prop, shape = 1, scale = 0.15),
         prop4 = pweibull(event, shape = 1, scale = 0.15),
         weib_prop = est_cdf(onset))  
get_onset = function(lambda, eta) {
  nsib = length(eta)
  U = runif(nsib)
  #onset
  lambda * sqrt( -log(U) / exp(eta) )
}
true2 = true %>% 
  mutate(censor2 = rweibull(n = 100e3, shape = 1, scale = 0.15),
         onset2 = get_onset(lambda = 0.15, eta = gen + env),
         status2 = onset2 <= censor2,
         event2 = pmin(onset2, censor2)) %>% 
  arrange(event2) %>% 
  mutate(prop2 = (cumsum(status2) + 1)/n(),
         est2  = qnorm(prop2, sd = sqrt(h2),lower.tail = F),
         age2 = LTFHPlus::liab_to_aoo(liab = est2, pop_prev = mean(status2)),
         lower = ifelse(status2, LTFHPlus::age_to_thres(age = age2, pop_prev = mean(status2)), -Inf),
         upper = ifelse(status2, LTFHPlus::age_to_thres(age = age2, pop_prev = mean(status2)), #alternatively, set upper bound to Inf for cases
                        LTFHPlus::age_to_thres(age = age2, pop_prev = mean(status2)))) %>%
  rename(ids = FID) %>% 
  relocate(ids, lower, upper)

hist(true2$est2)
true2[which((is.na(true2$age2))),]
mean(true2$status2)

hist(tmp$prop4)
hist(tmp$prop3)
hist(tmp$prop2)
hist(tmp$prop)


plot(tmp$event, mean(tmp$status) * event_cdf(tmp$event), col = c("red", "black")[tmp$status + 1])

hist(mean(tmp$status) * event_cdf(tmp$onset))

tmp %>% 
  ggplot(aes(x = event, fill = as.factor(status))) +
  geom_density(alpha = .3) 
tmp %>% ggplot() +
  geom_density(aes(x = onset, fill = "red"), alpha = .3) +
  geom_density(aes(x = censor, fill = "blue"), alpha = .3)

plot(tmp$onset, tmp$prop)  

tmp %>% 
  ggplot(aes(x = onset, fill = status)) +
  geom_density(alpha = .3)

hist(tmp$prop3)

tmp %>% ggplot(aes(x = prop, y = prop2)) + 
  geom_point() +
  geom_abline()
tmp %>% ggplot(aes(x = prop, y = prop3)) + 
  geom_point() +
  geom_abline()
tmp %>% 
  ggplot(aes(x = est, fill = status)) +
  geom_density(alpha = .3)
tmp %>% 
  ggplot(aes(x = prop, fill = status)) +
  geom_histogram()
tmp %>% 
  ggplot(aes(x = est2, fill = status)) +
  geom_histogram()

tmp %>% 
  ggplot(aes(x = event, fill = status)) +
  geom_density(alpha = .3)


tmp %>% 
  ggplot(aes(x = qnorm(prop, lower.tail = F), fill = status)) +
  geom_density(alpha = .3)

hist(qnorm(tmp$prop, lower.tail = F))
tmp %>% ggplot(aes(x = gen, y = prop)) + 
  geom_point()

tmp %>% ggplot(aes(x = gen + env, fill = status)) +
  geom_density(alpha = .4)

tmp %>% 
  ggplot(aes(x = onset, y = prop, color = status)) +
  geom_point(alpha = .3)
tmp %>% 
  ggplot(aes(x = censor, y = prop, color = status)) +
  geom_point(alpha = .3)


plot()
phen = fread(str_replace(cur_true, "true", "phen")) %>%
  as_tibble

ggplot(phen, aes(x = post_gen_no_fam, y = gen)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + 0) + 
  geom_abline()

tmp2 = left_join(tmp, noFH, by = c("FID" = "ids"))

tmp2 %>% 
  ggplot(aes(x = post_gen_no_fam, y = gen)) +
  geom_point(alpha = .3) +
  geom_abline() +
  theme_bw()

library(bigsnpr)
bigsnp_obj = snp_attach("./simulatedData/genotypes-SPACox-normal-maf10-v1.rds")
G = bigsnp_obj$genotypes

gwas     = big_univLinReg(G, tmp2$post_gen_no_fam, ind.col = 1:ncol(G), ncores = 3)

gwas$p = -predict(gwas)
snpinfo = fread("./simulatedData/genotypes-SPACox-normal-maf10-v1.snpinfo") %>% 
  as_tibble
gwas = cbind(gwas, beta = snpinfo$beta) 
hist(gwas$p)



tmp3 = readRDS("./gwas-results/v1_maf10_normal.rds")
predict_pvals = function() {
  lpval <- stats::pt(xtr, df = %d, lower.tail = FALSE, log.p = TRUE)
  (log(2) + lpval) / log(10)
}

all_gwas[[1]] %>% 
  mutate(p = -predict(all_gwas[[1]]),
         beta = snpinfo$beta) %>% 
  group_by(beta != 0) %>% 
  count(p > -log10(5e-8))

all_gwas[[2]] %>% 
  mutate(p = -predict(all_gwas[[2]]),
         beta = snpinfo$beta) %>% 
  group_by(beta != 0) %>% 
  count(p > -log10(5e-8))

all_gwas[[3]] %>% 
  mutate(beta = snpinfo$beta) %>% 
  group_by(beta != 0) %>% 
  count(p.value.spa < 5e-8)

table("betas" = gwas$beta != 0, "sig" = gwas$p > -log10(5e-8))

gwas %>% 
  ggplot(aes(x = p, fill = beta != 0)) +
  geom_density(alpha = .3)

all_gwas[[2]] %>% 
  mutate(p = -predict(all_gwas[[2]]),
         beta = snpinfo$beta) %>% 
  ggplot(aes(x = p, fill = beta != 0)) +
  geom_density(alpha = .3)

all_gwas[[3]] %>% 
  mutate(beta = snpinfo$beta) %>% 
  ggplot(aes(x = -log10(p.value.spa), fill = beta != 0)) +
  geom_density(alpha = .3)
