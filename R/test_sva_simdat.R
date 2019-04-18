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

# most_var = FALSE
# ms_dir <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/"
# wd <- "~/sva_tests/"

setwd(wd)

TP <- "FOM"

#load the samplesheet
load(paste0(ms_dir, "samplesheet/data.Robj"))
samplesheet <- dplyr::filter(samplesheet, time_point == TP) %>%
	dplyr::filter(is.na(duplicate.rm))
samplesheet$ALN <- as.integer(samplesheet$ALN)
#load the methylation data
load(paste0(ms_dir, "betas/data.Robj"))
meth <- beta[, samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(beta)
dim(meth)

phen_dir <- "~/main_project/ALSPAC_EWAS/methyl_variance/phen/"
# load the covars
covars <- read_delim(paste0(phen_dir, "FOM/FOM.qcovar"), delim = " ", col_names = F)
colnames(covars) <- c("ALN", "Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))

pheno <- samplesheet %>%
	left_join(covars) %>%
	dplyr::select(one_of(colnames(covars)), BCD_id, BCD_plate, MSA4Plate_id, Slide) %>%
	dplyr::filter(!is.na(Bcell))

mdata <- meth[complete.cases(meth), ]
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
	params_2 <- read_delim(file = paste0("results/sv_test_params_sims", mnam, ".txt"), delim = "\t")
}

# ---------------------------------------------------------------
# set up parameters to be tested
# ---------------------------------------------------------------
n_sv <- seq(5, 60, 5)
sv_type <- c("sva", "smartsva")
n_cpg <- c(seq(20000, 300000, 20000), nrow(mdata))
n_samp <- c(seq(100, ncol(mdata), by = 100))
dat_type <- c("binary", "continuous")
covs <- c("none", "cc")
params <- expand.grid(n_sv = n_sv, sv_type = sv_type, n_cpg = n_cpg, dat_type = dat_type, n_sample = n_samp, covs = covs)
params$time_user <- NA
params$time_system <- NA
params$time_elapsed <- NA
set.seed(2)
bin <- sample(0:1, ncol(mdata), replace=T)
table(bin)
cont <- rnorm(ncol(mdata))
phen <- data.frame(binary = bin, continuous = cont)


params <- params %>%
	mutate(include = case_when(sv_type == "smartsva" &
							   dat_type == "continuous" &
							   covs == "none" &
							   n_sv < 21 ~ TRUE, 
							   sv_type == "sva" & 
							   n_cpg == max(n_cpg) &
							   n_sample == max(n_sample) &
							   covs == "none" &
							   n_sv < 21 ~ TRUE, 
							   sv_type == "smartsva" &
							   dat_type == "binary" &
							   n_cpg == max(n_cpg) &
							   n_sample == max(n_sample) & 
							   covs == "none" &
							   n_sv < 21 ~ TRUE, 
							   sv_type == "smartsva" &
							   dat_type == "continuous" &
							   covs == "none" &
							   n_sample == max(n_sample) ~ TRUE, 
							   sv_type == "smartsva" & 
							   n_sample == min(n_sample) &
							   covs == "cc" &
							   dat_type == "continuous" &
							   n_sv == 20 ~ TRUE)) %>%
	dplyr::filter(include == TRUE)


if (exists("params_2")) {
	temp <- params %>%
		full_join(params_2) %>%
		arrange(time_user) # arranges it so that the NAs are at the bottom - good for sv_list later
	
	temp2 <- temp %>%
		dplyr::select(-time_user, -time_elapsed, -time_system, -include)

	dup_rows <- duplicated(temp2)
	
	params <- temp[!dup_rows,]
} else {
	sv_list <- list()
}

# ---------------------------------------------------------------
# run analyses
# ---------------------------------------------------------------
cc_cov <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")
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
	co <- params[i, "covs"]

	if (co == "none") {
		temp_mdata <- mdata
	} else {
		temp_mdata <- mdata[, pheno$Sample_Name]
	}
	samp <- sample(1:ncol(temp_mdata), nsamp)

	if (ncpg != nrow(temp_mdata)) {
		if (most_var) {
			var.idx <- order(rowVars(temp_mdata, na.rm=T), decreasing=T)[1:ncpg]	
		} else {
			var.idx <- sample(1:nrow(temp_mdata), ncpg)
		}
		mdat <- temp_mdata[var.idx, , drop=F]

	} else {
		mdat <- temp_mdata
	}

	mdat <- mdat[, samp]
	temp_phen <- phen[samp, ]
	dt_col <- ifelse(dt == "binary", "binary", "continuous")
	
	if (co == "none") {
		fom <- as.formula(paste0("~", dt))
		mod0=NULL
	} else if (co == "cc") {
		temp_phen <- temp_phen %>%
			mutate(Sample_Name = pheno[samp, "Sample_Name"]) %>%
			left_join(pheno)
		stopifnot(all(temp_phen$Sample_Name %in% colnames(mdat)))
		fom <- as.formula(paste0("~", paste(c(dt, cc_cov), collapse = " + ")))
		fom0 <- as.formula(paste0("~", paste(cc_cov, collapse = " + ")))
		mod0 <- model.matrix(fom0, data = temp_phen)
	}
	mod <- model.matrix(fom, data = temp_phen)

	ptm <- proc.time()
	if (svtype == "sva") {
		svobj <- sva(mdat, mod, mod0=mod0, n.sv = nsv, B = 5)
	} else if (svtype == "smartsva") {
		svobj <- smartsva.cpp(mdat, mod, mod0=mod0, n.sv = nsv, B = 5)
	}

	tim <- proc.time() - ptm
	params[i, "time_user"] <- tim[1]
	params[i, "time_system"] <- tim[2]
	params[i, "time_elapsed"] <- tim[3]

	if (!nsv %in% c(max(params$n_sv),20)) {
		next
	} else {
	## sv_list names follow this order: DT_NSV_SVTYPE_NCPG_NSAMP
	## --- need to shorted some of them
		dt <- substr(dt, 1, 3)
		ncpg <- paste0(round(ncpg/1000), "k")
		nam <- paste(dt, nsv, svtype, ncpg, nsamp, sep = "_")
		if (co == "cc") nam <- paste(nam, "cc", sep = "_")
		sv_list[[nam]] <- svobj	
	}
}

write.table(params, file = paste0("results/sv_test_params_sims", mnam, ".txt"), quote = F, col.names = T, row.names = F, sep = "\t")

save(sv_list, file = paste0("data/sv_test_sims", mnam, ".RData"))

print("FIN")







