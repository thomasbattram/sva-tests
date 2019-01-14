# ---------------------------------------------------------------
# Analysis of results from SVA 
# ---------------------------------------------------------------
rm(list = ls())

pkgs <- c("tidyverse", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

args <- commandArgs(trailingOnly = TRUE)
most_var <- as.logical(args[1])
wd <- args[2]

# wd <- "~/sva_tests/"
# most_var <- FALSE

setwd(wd)

if (most_var) {
	mnam <- "most_var"
} else {
	mnam <- ""
}

# load most var as well for comparison at the end
load("data/sv_test_simsmost_var.RData")
sv_list_mv <- sv_list
file <- "data/sv_test_sims.RData"
params <- read_delim(paste0("results/sv_test_params_sims", mnam, ".txt"), delim = "\t")
load(file)

length(sv_list)


# ---------------------------------------------------------------
# Check timings of the results
# ---------------------------------------------------------------
# graph it!!!

all_sv_cont_params <- params %>%
	dplyr::filter(dat_type == "continuous") %>% 
	dplyr::filter(n_sv == max(n_sv))

time_plot <- ggplot(all_sv_cont_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sample))) +
	geom_line(aes(linetype = sv_type)) +
	geom_point() +
	scale_colour_discrete(name = "n_sample")

sv_dat_type_params <- params %>%
	dplyr::filter(n_sample == max(n_sample)) %>%
	dplyr::filter(sv_type == "smartsva")

sv_time_plot <- ggplot(sv_dat_type_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sv))) +
	geom_point() +
	geom_line(aes(linetype = dat_type)) +
	scale_colour_discrete(name = "n_sv")

dat_type_time_plot <- ggplot(sv_dat_type_params, aes(x = n_cpg, y = time_user, colour = as.factor(dat_type))) +
	geom_point() +
	geom_line() +
	scale_fill_discrete(name = "dat_type")

# look at data type (binary or continuous)
bin <- params[params$dat_type == "binary", "time_user"]
cont <- params[params$dat_type == "continuous", "time_user"]
cor(bin, cont) # 0.9982354


# ---------------------------------------------------------------
# Test the differences in SVs generated
# ---------------------------------------------------------------
# parameters for:
# - sv_type -- keep smartsva
# - n_cpg -- keep 450k
# - dat_type -- keep continuous
# - n_sample -- keep max sample size
# - n_sv -- keep max number of SVs

# testing sv types
n_cpg <- max(params$n_cpg)
dat_type <- "continuous"
n_sample <- max(params$n_sample)
n_sv <- max(params$n_sv)

## sv_list names follow this order: DT_NSV_SVTYPE_NCPG_NSAMP
list_nam <- with(params, paste("con",
				 max(n_sv), 
				 unique(sv_type), 
				 paste0(round(max(n_cpg)/1000), "k"), 
				 max(n_sample), 
				 sep = "_"))

str(sv_list[576])
sva_450 <- sv_list[[list_nam[2]]]
smsva_450 <- sv_list[[list_nam[1]]]

# assess the correlation between the 20 SVs produced and see how much variance one explains of the other
smart_v_normal <- data.frame(sv = 1:20, cor = NA, adj_r2 = NA)
for (j in 1:nrow(smart_v_normal)) {
	sva_sv <- sva_450$sv[, j]
	smsva_sv <- smsva_450$sv[, j]
	smart_v_normal[j, 2] <- cor(smsva_sv, sva_sv)
	fit <- lm(smsva_sv ~ sva_sv)
	smart_v_normal[j, 3] <- summary(fit)$adj.r.squared
}

write.table(smart_v_normal, file = paste0("results/smartsva_v_sva_sims", mnam, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# 450k vs lower 
list_nam <- with(params, paste("con",
				 max(n_sv), 
				 "smartsva", 
				 paste0(round(unique(n_cpg)/1000), "k"), 
				 max(n_sample), 
				 sep = "_"))

ncpg_dat <- as.data.frame(matrix(NA, nrow = length(list_nam), ncol = 21))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:21] <- paste0("sv", 1:20, "_adjr2")

smsva_450 <- sv_list[[list_nam[grep("483k", list_nam)]]]
svs_450 <- smsva_450$sv
i=1
for (i in 1:length(list_nam)) {
	temp_svs <- sv_list[[list_nam[i]]]$sv
	ncpg <- str_split(list_nam[i], "_")[[1]][4]
	ncpg <- as.numeric(gsub("k", "", ncpg)) * 1000
	ncpg_dat[i, "n_cpg"] <- ncpg
	fom <- as.formula(paste0("svs_450[, j] ~ ", paste(paste0("temp_svs[, ", 1:20, "]"), collapse = " + ")))
	for (j in 1:ncol(svs_450)) {
		fit <- lm(fom)
		ncpg_dat[i, j+1] <- summary(fit)$adj.r.squared
	}
}

plot_res <- ncpg_dat %>%
	gather(key = sv, value = adj_r2, -n_cpg)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

# plot it!
ncpg_plot <- ggplot(plot_res, aes(x = n_cpg, y = adj_r2, colour = reorder(sv, sort(as.numeric(sv))))) +
	geom_line() +
	geom_point() +
	scale_colour_discrete(name = "SV")

write.table(ncpg_dat, file = paste0("results/ncpg_comp_sims", mnam, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# Using 10 SVs, 450k vs lower 
list_nam <- with(params, paste("con",
				 max(n_sv), 
				 "smartsva", 
				 paste0(round(unique(n_cpg)/1000), "k"), 
				 max(n_sample), 
				 sep = "_"))

ncpg_dat <- as.data.frame(matrix(NA, nrow = length(list_nam), ncol = 11))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:11] <- paste0("sv", 1:10, "_adjr2")

smsva_450 <- sv_list[[list_nam[grep("483k", list_nam)]]]
svs_450 <- smsva_450$sv[,1:10]

i=1
for (i in 1:length(list_nam)) {
	temp_svs <- sv_list[[list_nam[i]]]$sv[,1:10]
	ncpg <- str_split(list_nam[i], "_")[[1]][4]
	ncpg <- as.numeric(gsub("k", "", ncpg)) * 1000
	ncpg_dat[i, "n_cpg"] <- ncpg
	fom <- as.formula(paste0("svs_450[, j] ~ ", paste(paste0("temp_svs[, ", 1:10, "]"), collapse = " + ")))
	for (j in 1:ncol(svs_450)) {
		fit <- lm(fom)
		ncpg_dat[i, j+1] <- summary(fit)$adj.r.squared
	}
}

plot_res <- ncpg_dat %>%
	gather(key = sv, value = adj_r2, -n_cpg)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

# plot it!
ncpg_plot_10svs <- ggplot(plot_res, aes(x = n_cpg, y = adj_r2, colour = reorder(sv, sort(as.numeric(sv))))) +
	geom_line() +
	geom_point() +
	scale_colour_discrete(name = "SV")

write.table(ncpg_dat, file = paste0("results/ncpg_comp_10_sims", mnam, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")

## Taking most variable CpGs 
ncpg_dat <- as.data.frame(matrix(NA, nrow = length(list_nam), ncol = 11))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:11] <- paste0("sv", 1:10, "_adjr2")

smsva_450 <- sv_list_mv[[list_nam[grep("483k", list_nam)]]]
svs_450 <- smsva_450$sv[,1:10]

i=1
for (i in 1:length(list_nam)) {
	temp_svs <- sv_list_mv[[list_nam[i]]]$sv[,1:10]
	ncpg <- str_split(list_nam[i], "_")[[1]][4]
	ncpg <- as.numeric(gsub("k", "", ncpg)) * 1000
	ncpg_dat[i, "n_cpg"] <- ncpg
	fom <- as.formula(paste0("svs_450[, j] ~ ", paste(paste0("temp_svs[, ", 1:10, "]"), collapse = " + ")))
	for (j in 1:ncol(svs_450)) {
		fit <- lm(fom)
		ncpg_dat[i, j+1] <- summary(fit)$adj.r.squared
	}
}

ncpg_dat_mv <- ncpg_dat %>%
	dplyr::mutate(cpg_subset = "most_variable")

ncpg_dat <- read_delim("results/ncpg_comp_10_sims.txt", delim = "\t")

ncpg_dat <- ncpg_dat %>%
	dplyr::mutate(cpg_subset = "random") %>%
	rbind(ncpg_dat_mv)

write.table(ncpg_dat, "results/mv_v_rand_ncpg_comp_10_sims.txt", row.names = F, col.names = T, quote = F, sep = "\t")

plot_res <- ncpg_dat %>%
	gather(key = sv, value = adj_r2, -n_cpg, -cpg_subset)
plot_res$sv <- gsub("_adjr2", "", plot_res$sv)
plot_res$sv <- gsub("sv", "", plot_res$sv)

# plot it!
rand_v_mv_plot <- ggplot(plot_res, aes(x = n_cpg, y = adj_r2, colour = reorder(sv, sort(as.numeric(sv))))) +
	geom_line(aes(linetype = cpg_subset)) +
	geom_point() +
	scale_colour_discrete(name = "SV")


# final plots to summarise it!
pdf(paste0("results/sv_plots", mnam, ".pdf"))
marrangeGrob(list(time_plot, sv_time_plot, ncpg_plot, ncpg_plot_10svs, rand_v_mv_plot), ncol = 1, nrow = 1)
dev.off()

