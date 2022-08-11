library(data.table)
library(dplyr)
library(stringr)

# LTM downsampling --------------------------------------------------------


all_files = list.files("./simulatedData", "LTM", full.names = T) 


phenotypes = lapply(1:10, function(v) {
  ph = str_subset(all_files, paste0("v", v,"\\.")) %>% 
    str_subset(., "phen")
  tibble(full_names = ph,
         v = v) %>% 
    mutate(C = as.vector(str_match(full_names, "C1000|C250")),
           prev = as.vector(str_match(full_names, paste0("prev", 100*c(0.05, 0.1, 0.2, 0.5), "-", collapse = "|"))),
           prev = str_replace(prev, "-", "")) %>% 
    group_by(v, C, prev) %>% 
    summarise(true = full_names, .groups = "drop") %>% 
    ungroup
}) %>% 
  do.call("rbind", .) %>% 
  filter(prev != "prev50") # we technically do not do any downsampling for prev50.


for (i in 1:nrow(phenotypes)) {
  print(i)
  cur_file = phenotypes$true[i]
  true = fread(cur_file) %>% 
    as_tibble
  prev_val = (str_replace(phenotypes$prev[i], "prev", "") %>% as.numeric()) / 100
  n_tot    = 2 * nrow(true) * prev_val
  
  tmp = true %>%
    mutate(sample_weights = ifelse(status == 1, 1, 1e-4)) %>% 
    slice_sample(n = n_tot, weight_by = sample_weights) %>%
    select(ids) %>% 
    arrange(ids)
  
  fwrite(tmp, 
         str_replace(cur_file, "phen", "keepers"),
         na = "NA",
         quote = F,
         sep = " ")
}
