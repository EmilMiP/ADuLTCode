library(data.table)
library(dplyr)
library(stringr)

prevs = c(0.0101, 0.025, 0.05, 0.1, 0.2, 0.5)
prev_text_vec = paste0("prev",round(prevs * 100,1))
prev_text_vec[1] = "prev01"
prev_text_vec[2] = "prev025"
#prev_text_vec[3] = "prev05"



# get prevalences first ---------------------------------------------------

all_files = list.files("./simulatedData", full.names = T, "phen") 

phenotypes = tibble(ph = all_files) %>%   
  mutate(gen_beta = as.vector(str_match(ph, "normal")),
         C = as.vector(str_match(ph, "C1000|C250")),
         prev = as.vector(str_match(ph, paste(prev_text_vec, collapse = "-|"))),
         prev = str_replace(prev, "-", ""),
         model = as.vector(str_match(ph, "SPACox|LTM"))) %>% 
  group_by(model, gen_beta, C, prev) %>% 
  summarise(grps = list(ph), .groups = "drop") %>% 
  ungroup




# Downsampling ------------------------------------------------------------
# we wont to fix the sample size to be 20k. 10k controls and 10k cases.

lapply(1:nrow(phenotypes), function(i) {
  print(i)
  lapply(phenotypes$grps[[i]], function(cur_phen) {
    phen = fread(cur_phen) %>% as_tibble()
    
    N_tot = 20e3
    n_cases = ifelse(sum(phen$status) > 10e3, 10e3, sum(phen$status))
    n_ctrls = N_tot - n_cases
    
    ds = bind_rows(
      phen %>% 
        filter(status == 0) %>% 
        slice_sample(n = n_ctrls),
      phen %>% 
        filter(status == 1) %>% 
        slice_sample(n = n_cases)
    ) %>% 
      select(ids) %>% 
      arrange(ids)
    
    fwrite(ds, 
           str_replace(cur_phen, "phen", "keepers"),
           na = "NA",
           quote = F,
           sep = " ")
  })
})
for (i in 1:nrow(phenotypes)) {
  print(i)
  cur_file = phenotypes$true[i]
  true = fread(cur_file) %>% 
    as_tibble
  prev_val = phenotypes2$prev[i]
  n_tot    = 2 * nrow(true) * prev_val
  
  tmp = true %>%
    mutate(sample_weights = ifelse(status == 1, 1, 1e-4)) %>% 
    slice_sample(n = n_tot, weight_by = sample_weights) %>%
    select(FID) %>% 
    arrange(FID)
  
  fwrite(tmp, 
         str_replace(cur_file, "true", "keepers"),
         na = "NA",
         quote = F,
         sep = " ")
}