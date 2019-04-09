## ---- load_data -----------------------------------

rm(list = ls())

pkgs <- c("tidyverse", "knitr", "captioner", "pander", "gridExtra")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc")

# params tables
params <- read_delim("~/sva_tests/results/sv_test_params_sims.txt", delim = "\t")
alp_params <- read_delim("~/sva_tests/results/alpha_tests_params.txt", delim = "\t")

# sva vs smart sva
smart_v_normal <- read_delim("~/sva_tests/results/smartsva_v_sva_sims.txt", delim = "\t")

# test of default smartsva parameters
alpha_tests <- read_delim("~/sva_tests/results/alpha_tests.txt", delim = "\t")

# ncpg with 20 svs
ncpg_dat_20 <- read_delim("~/sva_tests/results/ncpg_comp_sims.txt", delim = "\t")

# ncpg with 10 svs
ncpg_dat_10 <- read_delim("~/sva_tests/results/ncpg_comp_10_sims.txt", delim = "\t")

# ncpg with 20 svs using EPIC array 
ncpg_dat_850k <- read_delim("~/sva_tests/results/ncpg_comp_sims_850.txt", delim = "\t")

# mv and rand ncpg dat 
mv_vs_rand <- read_delim("~/sva_tests/results/mv_v_rand_ncpg_comp_10_sims.txt", delim = "\t")

# sv number
estimated_num <- read_delim("~/sva_tests/data/sim_estimated_sv_num.txt", delim = "\t")
load("~/sva_tests/data/cov_r2_res.RData")
nsv_dat <- read_delim("~/sva_tests/data/estimated_sv_num.txt", delim = "\t")
svs_needed <- read_delim("~/sva_tests/results/svs_needed.txt", delim = "\t")

# Captioner setup
table_nums <- captioner(prefix = "Table")
fig_nums <- captioner()
sup_table_nums <- captioner(prefix = "Supplementary Table")
sup_figure_nums <- captioner(prefix = "Supplementary Figure")

## ---- timings_setup -----------------------------------
# How does time taken vary with CpG num, SVA package and sample number
cont_params <- params %>%
	dplyr::filter(dat_type == "continuous") %>%
	dplyr::filter(n_sv == 20)

time_plot <- ggplot(cont_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sample),
									 group = as.factor(n_sample))) +
								 geom_line() +
								 geom_point() +
								 scale_colour_discrete(name = "n_sample")

fig_nums(name = "time_plot", caption = "How does time taken to perform SVA vary with sample number, SVA package and CpG number?")
time_plot_cap <- fig_nums("time_plot")

# smartsva package vs sva package correlation
colnames(smart_v_normal) <- c("SV", "correlation", bquote("adjusted r"^"2"))

smsva_times <- params[params$sv_type == "smartsva", "time_user", drop = T]
sva_times <- params[params$sv_type == "sva", "time_user", drop = T]

table_nums(name = "smart_v_normal", caption = "What is the correlation between SVs generated using the two different packages?")
smart_v_normal_cap <- table_nums("smart_v_normal")

# how does time taken vary with number of SVs and distribution of phenotype
sv_dat_type_params <- params %>%
	dplyr::filter(n_sample == max(n_sample)) %>%
	dplyr::filter(sv_type == "smartsva")

sv_dat_type_plot <- ggplot(sv_dat_type_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sv))) +
	geom_point() +
	geom_line(aes(linetype = dat_type)) +
	scale_colour_discrete(name = "n_sv")

fig_nums(name = "time_plot2", caption = "How does time taken to perform SVA vary with number of SVs and distribution of phenotype (continuous or binary)?")
time_plot2_cap <- fig_nums("time_plot2")

# comparison of smartsva parameters --> smartsva defaults vs. sva package defaults
model_param_comp <- ggplot(alp_params, aes(x = alpha, y = time_user, fill = as.factor(B))) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_grid(. ~ as.factor(n_sv))

fig_nums(name = "model_param_time_plot", caption = "How does time taken to perform SVA vary with alpha and B values?")
model_param_time_plot_cap <- fig_nums("model_param_time_plot")

## ---- time_plot -----------------------------------
print(time_plot)

## ---- sva_vs_smart_sva_table -----------------------------------
pander(smart_v_normal)

## ---- time_plot2 -----------------------------------
print(sv_dat_type_plot)

## ---- model_param_time_plot -----------------------------------
print(model_param_comp)

## ---- ncpg_setup -----------------------------------
# subsetting data to most variable CpGs or a random subset
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

# variance of SVs created using all CpGs explained by SVs created using a subset
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

# now using an EPIC array dataset! 
plot_res <- ncpg_dat_850k %>%
	gather(key = sv, value = adj_r2, -n_cpg)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

# plot it!
ncpg_plot_850k <- ggplot(plot_res, aes(x = n_cpg, y = adj_r2, colour = reorder(sv, sort(as.numeric(sv))))) +
	geom_line() +
	geom_point() +
	scale_colour_discrete(name = "SV")

fig_nums(name = "ncpg_plot_850k", caption = "Variance captured of SVs made using EPIC array CpGs by all SVs made using random subsets of varying numbers of CpGs")
ncpg_plot_850k_cap <- fig_nums("ncpg_plot_850k")

## ---- mv_vs_random -----------------------------------
print(mv_vs_random_plot)

## ---- ncpg_plot -----------------------------------
print(ncpg_plot)

## ---- ncpg_plot_850k -----------------------------------
print(ncpg_plot_850k)

## ---- nsv_setup -----------------------------------
# Comparison of rmt and num.sv() number of SVs estimated
g_estimated_num <- gather(estimated_num, key = "method", value = "nsv", -ncpg, -nsamp) %>%
	mutate(method = gsub("_nsv", "", method))

est_num_plot <- ggplot(g_estimated_num, aes(x = ncpg, y = nsv, colour = as.factor(nsamp))) +
	geom_point() +
	geom_line(aes(linetype = method)) +
	scale_colour_discrete(name = "n_sample")

fig_nums(name = "est_num_plot", caption = "Number of SVs required estimated using num.sv() and random matrix theory")
est_num_plot_cap <- fig_nums("est_num_plot")

# variance of common EWAS covariates explained by SVs
# take just two traits
sm_sva_res_list <- sva_res_list[grep("Glucose|rcont", names(sva_res_list))]
max_sv <- map_dbl(sm_sva_res_list, function(x) max(as.numeric(x$sv)))
names(max_sv) <- c("rcont", "glc")
sva_plot_res <- do.call(rbind, sm_sva_res_list) %>%
	rownames_to_column(var = "trait") %>%
	mutate(trait = gsub("\\..*", "", trait)) %>%
	mutate(trait = ifelse(grepl("Glucose", trait), "glc", trait)) %>%
	gather(key = "covariate", value = "adjusted_r2", -sv, -trait) %>%
	mutate(cov_type = case_when(covariate %in% cc ~ "cell count", 
								covariate %in% batch ~ "batch", 
								covariate == "age" ~ "age", 
								covariate %in% paste0("PC", 1:10) ~ "PC")) %>%
	mutate(sv = as.numeric(sv)) %>%
	mutate(lab = case_when(trait == "glc" & 
						   sv == max_sv["glc"] & 
						   !(covariate %in% paste0("PC", 1:10)) 
						   & !(covariate == "age") ~ covariate,
						   sv != max(sv) ~ "", 
						   trait == "rcont" & 
						   sv == max_sv["rcont"] & 
						   !(covariate %in% paste0("PC", 1:10)) & 
						   !(covariate == "age") ~ covariate))

	xlimits <- c(max(g_res$sv), max(g_res$sv) + 10)

covs_r2_p <- ggplot(sva_plot_res, aes(x = sv, y = adjusted_r2, colour = cov_type, group = covariate)) +
	geom_point() +
	geom_line() +
	geom_text_repel(aes(label = lab), xlim = xlimits, nudge_x = 1, cex = 2.5) +
	xlim(NA, max(xlimits)) +
	labs(x = "SV") +
	facet_grid(trait ~ .) 

fig_nums(name = "covs_var_exp", caption = "Variance of important covariates explained by SVs")
covs_var_exp_cap <- fig_nums("covs_var_exp")

# for the table
sim_res <- dplyr::filter(svs_needed, data == "simulated")
real_res <- dplyr::filter(svs_needed, data == "real")

sum_res <- rbind(
sim_res %>% dplyr::filter(max_sv == "not_max") %>%
	group_by(cov) %>%
	summarise(median_sv_95 = median(as.numeric(sv))) %>%
	left_join(
		sim_res %>% dplyr::filter(max_sv == "max") %>%
		group_by(cov) %>%
		summarise(median_sv_max = median(as.numeric(sv)))
		) %>%
	mutate(data = "simulated")
,
real_res %>% dplyr::filter(max_sv == "not_max") %>%
	group_by(cov) %>%
	summarise(median_sv_95 = median(as.numeric(sv))) %>%
	left_join(
		real_res %>% dplyr::filter(max_sv == "max") %>%
		group_by(cov) %>%
		summarise(median_sv_max = median(as.numeric(sv)))
		) %>%
	mutate(data = "real")
)

table_nums(name = "nsv_summary", caption = "Summary of variance of important covariates explained by SVs")
nsv_summary_cap <- table_nums("nsv_summary")

## ---- est_num_plot -----------------------------------
print(est_num_plot)

## ---- covs_var_exp -----------------------------------
print(covs_r2_p)

## ---- nsv_summary -----------------------------------
pander(sum_res)
