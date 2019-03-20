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
	params <- read_delim(file = paste0("results/sv_test_params_sims", mnam, ".txt"), delim = "\t")
	params_2 <- params
}

# ---------------------------------------------------------------
# set up parameters to be tested
# ---------------------------------------------------------------
n_sv <- seq(5, 60, 5)
sv_type <- c("sva", "smartsva")
sv_type <- sv_type[2]
n_cpg <- c(seq(20000, 300000, 20000), nrow(mdata))
n_samp <- c(seq(100, ncol(mdata), by = 100))
n_samp <- max(n_samp)
dat_type <- c("binary", "continuous")
dat_type <- dat_type[2]
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
i=1
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

	tim <- proc.time() - ptm
	params[i, "time_user"] <- tim[1]
	params[i, "time_system"] <- tim[2]
	params[i, "time_elapsed"] <- tim[3]

	if (nsv != max(params$n_sv) {
		next
	} else {
	## sv_list names follow this order: DT_NSV_SVTYPE_NCPG_NSAMP
	## --- need to shorted some of them
		dt <- substr(dt, 1, 3)
		ncpg <- paste0(round(ncpg/1000), "k")
		nam <- paste(dt, nsv, svtype, ncpg, nsamp, sep = "_")
		sv_list[[nam]] <- svobj	
	}
}

write.table(params, file = paste0("results/sv_test_params_sims", mnam, ".txt"), quote = F, col.names = T, row.names = F, sep = "\t")

save(sv_list, file = paste0("data/sv_test_sims", mnam, ".RData"))

print("FIN")







