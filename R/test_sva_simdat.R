# ---------------------------------------------------------------
# Testing parameters for SVA
# ---------------------------------------------------------------

rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

args <- commandArgs(trailingOnly = TRUE)
most_var <- as.logical(args[1])
ms_dir <- args[2]
wd <- args[3]

# ms_dir <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/"
# wd <- "~/sva_tests/"

setwd(wd)

TP <- "FOM"

#load the samplesheet
load(paste0(ms_dir, "samplesheet/data.Robj"))
samplesheet <- dplyr::filter(samplesheet, time_point == TP) %>%
	dplyr::filter(is.na(duplicate.rm))

#load the methylation data
load(paste0(ms_dir, "betas/data.Robj"))
meth <- beta[, samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(beta)
dim(meth)

mdata <- as.matrix(meth)
mdata <- mdata[complete.cases(mdata), ]
dim(mdata)

# load in the old data if there 
if (most_var) {
	mnam <- "most_var"
} else {
	mnam <- ""
}

file <- paste0("data/sv_test_sims", mnam, ".RData")

if (file.exists(file)) {
	load(file)
	params_2 <- params
}

# ---------------------------------------------------------------
# set up parameters to be tested
# ---------------------------------------------------------------
n_sv <- 20
sv_type <- c("sva", "smartsva")
n_cpg <- c(seq(20000, 300000, 20000), nrow(mdata))
n_samp <- c(seq(100, ncol(mdata), by = 100))
dat_type <- c("binary", "continuous")
params <- expand.grid(n_sv = n_sv, sv_type = sv_type, n_cpg = n_cpg, dat_type = dat_type, n_sample = n_samp)
params$time_user <- NA
params$time_system <- NA
params$time_elapsed <- NA
set.seed(2)
bin <- sample(0:1, ncol(mdata), replace=T)
table(bin)
cont <- rnorm(ncol(mdata))
phen <- data.frame(binary = bin, continuous = cont)

if (exists("params_2")) {
	temp <- params %>%
		dplyr::select(-time_user, -time_system, -time_elapsed) %>%
		left_join(params_2) %>%
		arrange(time_user) # arranges it so that the NAs are at the bottom - good for sv_list later
	
	params <- temp
} else {
	sv_list <- list()
}

# ---------------------------------------------------------------
# run analyses
# ---------------------------------------------------------------

for (i in 1:nrow(params)) {
	print(i)
	if (!is.na(params[i, "time_user"])) {
		next
	}
	nsv <- params[i, "n_sv"]
	svtype <- params[i, "sv_type"]
	ncpg <- params[i, "n_cpg"]
	dt <- params[i, "dat_type"]
	nsamp <- params[i, "n_sample"]
	set.seed(2)
	samp <- sample(1:ncol(mdata), nsamp)

	if (ncpg != nrow(mdata)) {
		if (most_var) {
			var.idx <- order(rowVars(mdata, na.rm=T), decreasing=T)[1:ncpg]	
		} else {
			var.idx <- sample(1:nrow(mdata), ncpg)
		}
		mdat <- mdata[var.idx, , drop=F]
	} else {
		mdat <- mdata
	}

	mdat <- mdat[, samp]
	temp_phen <- phen[samp, ]
	dt_col <- ifelse(dt == "binary", "binary", "continuous")
	
	fom <- as.formula(paste0("~", dt))
	mod <- model.matrix(fom, data = temp_phen)
	mod0 <- model.matrix(~1, data = temp_phen)
	ptm <- proc.time()
	if (svtype == "sva") {
		svobj <- sva(mdat, mod, mod0=NULL, n.sv = nsv, B = 5)
	} else if (svtype == "smartsva") {
		svobj <- smartsva.cpp(mdat, mod, mod0=NULL, n.sv = nsv, B = 5)
	}

	sv_list[[paste0("p", i)]] <- svobj
	
	tim <- proc.time() - ptm
	params[i, "time_user"] <- tim[1]
	params[i, "time_system"] <- tim[2]
	params[i, "time_elapsed"] <- tim[3]
}

write.table(params, file = paste0("results/sv_test_params_sims", mnam, ".txt"), quote = F, col.names = T, row.names = F, sep = "\t")

save(params, sv_list, file = paste0("data/sv_test_sims", mnam, ".RData"))

print("FIN")

# ---------------------------------------------------------------
# Check timings of the results
# ---------------------------------------------------------------

# graph it!!!


cont_params <- params %>%
	dplyr::filter(dat_type == "continuous")

time_plot <- ggplot(cont_params, aes(x = n_cpg, y = time_user, colour = as.factor(n_sample))) +
	geom_line(aes(linetype = sv_type)) +
	geom_point() +
	scale_colour_discrete(name = "n_sample")

# ---------------------------------------------------------------
# Test the differences in SVs generated
# ---------------------------------------------------------------
# parameters for:
# - sv_type -- keep smartsva
# - n_cpg -- keep 450k
# - dat_type -- keep continuous
# - n_sample -- keep max sample size

# testing sv types

sv_type_params <- with(params, which(dat_type == "continuous" 
					& n_cpg == max(n_cpg)
					& n_sample == max(n_sample)))

sva_450 <- sv_list[[sv_type_params[1]]]
smsva_450 <- sv_list[[sv_type_params[2]]]
params[sv_type_params,]
str(sv_list[sv_type_params[2]])

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
sv_ncpg_params <- with(params, which(dat_type == "continuous" 
					& n_sample == max(n_sample) 
					& sv_type == "smartsva"))
params[sv_ncpg_params,]

ncpg_dat <- as.data.frame(matrix(NA, nrow = length(sv_ncpg_params), ncol = 21))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:21] <- paste0("sv", 1:20, "_adjr2")

smsva_450 <- sv_list[[sv_ncpg_params[length(sv_ncpg_params)]]]
svs_450 <- smsva_450$sv
i=1
for (i in 1:length(sv_ncpg_params)) {
	temp_svs <- sv_list[[sv_ncpg_params[i]]]$sv
	ncpg_dat[i, "n_cpg"] <- params[sv_ncpg_params[i], "n_cpg"]
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
sv_ncpg_params <- with(params, which(dat_type == "continuous" 
					& n_sample == max(n_sample) 
					& sv_type == "smartsva"))
params[sv_ncpg_params,]

ncpg_dat <- as.data.frame(matrix(NA, nrow = length(sv_ncpg_params), ncol = 11))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:11] <- paste0("sv", 1:10, "_adjr2")

smsva_450 <- sv_list[[sv_ncpg_params[length(sv_ncpg_params)]]]
svs_450 <- smsva_450$sv[,1:10]
i=1
for (i in 1:length(sv_ncpg_params)) {
	temp_svs <- sv_list[[sv_ncpg_params[i]]]$sv[,1:10]
	ncpg_dat[i, "n_cpg"] <- params[sv_ncpg_params[i], "n_cpg"]
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

write.table(ncpg_dat, file = paste0("results/ncpg_comp_sims", mnam, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")


# continuous vs binary
sv_dattype_params <- with(params, which(n_cpg == max(n_cpg) 
					& n_sample == max(n_sample) 
					& sv_type == "smartsva"))
params[sv_dattype_params,]

# sample size changes
sv_nsamp_params <- with(params, which(dat_type == "continuous" 
					& n_cpg == max(n_cpg) 
					& sv_type == "smartsva"))

params[sv_nsamp_params,]

# final plots to summarise it! 
pdf(paste0("results/sv_plots", mnam, ".pdf"))
marrangeGrob(list(time_plot, ncpg_plot, ncpg_plot_10svs), ncol = 1, nrow = 1)
dev.off()








