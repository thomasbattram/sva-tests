## ---- load_data -----------------------------------

rm(list = ls())

pkgs <- c("tidyverse", "knitr", "captioner", "pander")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc")

# NEED TO SAVE MORE ANALYSIS RESULTS!!!!

# params tables
params <- read_delim("~/sva_tests/results/sv_test_params_sims.txt", delim = "\t")

# sva vs smart sva
smart_v_normal <- read_delim("~/sva_tests/results/smartsva_v_sva_sims.txt", delim = "\t")
smart_v_normal_mv <- read_delim("~/sva_tests/results/smartsva_v_sva_simsmost_var.txt", delim = "\t")

# ncpg with 20 svs
ncpg_dat_20 <- read_delim("~/sva_tests/results/ncpg_comp_sims.txt", delim = "\t")

# ncpg with 10 svs
ncpg_dat_10 <- read_delim("~/sva_tests/results/ncpg_comp_10_sims.txt", delim = "\t")

# mv and rand ncpg dat 
mv_vs_rand <- read_delim("~/sva_tests/results/mv_v_rand_ncpg_comp_10_sims.txt", delim = "\t")

# Captioner setup
table_nums <- captioner(prefix = "Table")
fig_nums <- captioner()

## ---- timings_setup -----------------------------------
cont_params <- params %>%
	dplyr::filter(dat_type == "continuous") %>%
	dplyr::filter(n_sv == max(n_sv))

time_plot <- ggplot(cont_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sample))) +
	geom_line(aes(linetype = sv_type)) +
	geom_point() +
	scale_colour_discrete(name = "n_sample")

fig_nums(name = "time_plot", caption = "How does time taken to perform SVA vary with sample number, SVA package and CpG number?")
time_plot_cap <- fig_nums("time_plot")

colnames(smart_v_normal) <- c("SV", "correlation", bquote("adjusted r"^"2"))

smsva_times <- params[params$sv_type == "smartsva", "time_user", drop = T]
sva_times <- params[params$sv_type == "sva", "time_user", drop = T]

table_nums(name = "smart_v_normal", caption = "What is the correlation between SVs generated using the two different packages?")
smart_v_normal_cap <- table_nums("smart_v_normal")

sv_dat_type_params <- params %>%
	dplyr::filter(n_sample == max(n_sample)) %>%
	dplyr::filter(sv_type == "smartsva")

sv_dat_type_plot <- ggplot(sv_dat_type_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sv))) +
	geom_point() +
	geom_line(aes(linetype = dat_type)) +
	scale_colour_discrete(name = "n_sv")

fig_nums(name = "time_plot2", caption = "How does time taken to perform SVA vary with number of SVs and distribution of data (continuous or binary)?")
time_plot2_cap <- fig_nums("time_plot2")

## ---- time_plot -----------------------------------
print(time_plot)

## ---- sva_vs_smart_sva_table -----------------------------------
pander(smart_v_normal)

## ---- time_plot2 -----------------------------------
print(sv_dat_type_plot)

## ---- ncpg_setup -----------------------------------

plot_res <- mv_vs_rand %>%
	gather(key = sv, value = adj_r2, -n_cpg, -cpg_subset)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

mv_vs_random_plot_res <- plot_res %>%
	dplyr::filter(n_cpg == 20000) 

mv_vs_random_plot <- ggplot(mv_vs_random_plot_res, aes(x = reorder(as.numeric(sv), sort(as.numeric(sv))), y = adj_r2, colour = as.factor(cpg_subset))) +
	geom_point() +
	geom_line(aes(group = cpg_subset)) +
	labs(x = "SV (created using all 450k CpGs)", y = bquote("Variance explained by all SVs created from a subset of CpGs (adj" ~r^2~ ")")) +
	scale_colour_discrete(name = "", 
						  breaks = c("random", "most_variable"), 
						  labels = c("random", "most variable"))

fig_nums(name = "mv_vs_random_plot", caption = "Is it better to subset the number of CpGs randomly or by most variable CpGs when running SVA?")
mv_vs_random_plot_cap <- fig_nums("mv_vs_random_plot")

plot_res <- ncpg_dat_20 %>%
	gather(key = sv, value = adj_r2, -n_cpg)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

ncpg_plot <- ggplot(plot_res, aes(x = n_cpg, y = adj_r2, colour = reorder(as.numeric(sv), sort(as.numeric(sv))))) +
	geom_line() +
	geom_point() +
	scale_colour_discrete(name = "SV")

fig_nums(name = "ncpg_plot", caption = "Variance captured of SVs made using 450k CpGs by all SVs made using random subsets of varying numbers of CpGs")
ncpg_plot_cap <- fig_nums("ncpg_plot")

## ---- mv_vs_random -----------------------------------
print(mv_vs_random_plot)

## ---- ncpg_plot -----------------------------------
print(ncpg_plot)



