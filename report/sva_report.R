## ---- load_data -----------------------------------

rm(list = ls())

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = TRUE)

# NEED TO SAVE MORE ANALYSIS RESULTS!!!!

# params tables
params <- read_delim("~/main_project/ALSPAC_EWAS/sva_tests/sv_test_params_sims.txt", delim = "\t")

# sva vs smart sva
smart_v_normal <- read_delim("~/main_project/ALSPAC_EWAS/sva_tests/smartsva_v_sva_sims.txt", delim = "\t")
smart_v_normal_mv <- read_delim("~/main_project/ALSPAC_EWAS/sva_tests/smartsva_v_sva_simsmost_var.txt", delim = "\t")

# ncpg with 20 svs
ncpg_dat_20 <- read_delim("~/main_project/ALSPAC_EWAS/sva_tests/ncpg_comp_sims.txt", delim = "\t")
ncpg_dat_20_mv <-read_delim("~/main_project/ALSPAC_EWAS/sva_tests/ncpg_comp_simsmost_var.txt", delim = "\t")

# ncpg with 10 svs
ncpg_dat_10 <- read_delim("~/main_project/ALSPAC_EWAS/sva_tests/ncpg_comp_10_sims.txt", delim = "\t")

## ---- timings_setup -----------------------------------
cont_params <- params %>%
	dplyr::filter(dat_type == "continuous")

time_plot <- ggplot(cont_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sample))) +
	geom_line(aes(linetype = sv_type)) +
	geom_point() +
	scale_colour_discrete(name = "n_sample")

## ---- time_plot -----------------------------------
print(time_plot)

## ---- sva_vs_smart_sva_setup -----------------------------------

