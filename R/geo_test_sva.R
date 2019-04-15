rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
geo_dat_path <- args[2]

# wd <- "~/sva_tests/"
# geo_dat_path <- "~/ewas_catalog/geo_data/"
setwd(wd)

# TP <- "FOM"

# #load the samplesheet
# load(paste0(ms_dir, "samplesheet/data.Robj"))
# samplesheet <- dplyr::filter(samplesheet, time_point == TP) %>%
# 	dplyr::filter(is.na(duplicate.rm)) %>%
# 	mutate(ALN = as.numeric(ALN))

# # sort the covariate date
# covars <- read_delim(paste0(phen_dir, "FOM/FOM.qcovar"), delim = " ", col_names = F)
# colnames(covars) <- c("ALN", "Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))

# pheno <- samplesheet %>%
# 	left_join(covars) %>%
# 	dplyr::select(one_of(colnames(covars)), BCD_id, BCD_plate, MSA4Plate_id, Slide) %>%
# 	dplyr::filter(!is.na(PC1))

# #load the methylation data
# load(paste0(ms_dir, "betas/data.Robj"))
# meth <- beta[, pheno$Sample_Name] #keep the samples that correspond to the time point you're interested in
# rm(beta)
# dim(meth)

# load the potential geo datasets
# dat_path <- "~/ewas_catalog/geo_data/"
load(paste0(geo_dat_path, "ewas-cat-cr02.rdata"))
load(paste0(geo_dat_path, "ewas-cat-cr02-get-data.rdata"))
gses <- geo.data %>%
	mutate(gse_id = tolower(rownames(.))) %>%
	dplyr::filter(ncol > 100) %>%
	dplyr::filter(nrow > 4.8e5) %>%
	dplyr::select(gse_id, everything())
	
# add the 850k dataset
gses <- rbind(gses, 
			  data.frame(gse_id = "gse116339", class = NA, nrow = NA, ncol = NA, length = NA))

# now the 850k dataset!
# x <- load(paste0(dat_path, "gse116339.rda"))
# meth_850k <- get(x)
# meth_850k <- meth_850k$data
# rm(list = x)
# rm(x)

# function for estimating number of SVs using random matrix theory
est_nsv <- function(meth, traits, df) {
	fom <- as.formula(paste0("t(meth)", " ~ ", paste(traits, collapse = " + ")))
	y.r <- t(resid(lm(fom, data = df)))
	print("Estimating number of SVs needed")
	n_sv <- EstDimRMT(y.r, FALSE)$dim + 1
	return(n_sv)
}


files <- list.files(geo_dat_path)
gse_id <- gses$gse_id[19]
nsv_dat <- lapply(gses$gse_id, function(gse_id) {
	if (!any(grepl(gse_id, files))) return(NULL)
	print("Loading file!")
	loaded_meth_dat_name <- load(paste0(geo_dat_path, gse_id, ".rda"))
	meth <- get(loaded_meth_dat_name)
	if (class(meth) == "list") meth <- meth$data
	rm(list = loaded_meth_dat_name)
	rm(loaded_meth_dat_name)
	print(dim(meth))

	# phenotype data
	phen <- data.frame(x = rnorm(ncol(meth)))

	# ----------------------------------------------------------------
	# Estimating the number of SVs 
	# ----------------------------------------------------------------
	print("Estimating number of SVs")
	# remove NAs
	meth <- na.omit(meth)
	print(dim(meth))
	if (nrow(meth) < 3e5) return(NULL)
	
	rmt_nsv <- est_nsv(as.matrix(meth), "x", phen)

	fom <- as.formula("~x")
	mod <- model.matrix(fom, data = phen)
	leek_nsv <- num.sv(as.matrix(meth), mod, method = "leek")
	nsv_out_dat <- data.frame(gse_id = gse_id, rmt = rmt_nsv, leek = leek_nsv,
							  n_sample = ncol(meth), n_cpg = nrow(meth))

	# ----------------------------------------------------------------
	# Running the SVA 
	# ----------------------------------------------------------------
	if (nrow(meth) > 5e5) {
		n_cpg <- c(seq(20000, 700000, 20000), nrow(meth))
	} else {
		n_cpg <- c(seq(20000, 300000, 20000), nrow(meth))
	}
	print("Running SVA")

	sv_list <- lapply(n_cpg, function(ncpg) {
		print(ncpg)
		set.seed(2)
		mdat <- meth[sample(1:nrow(meth), ncpg), ]
		svobj <- tryCatch({
			smartsva.cpp(mdat, mod, mod0 = NULL, n.sv = 20)
		}, error = function(e) {err_msg(e)})
		return(svobj)
	})
	# save results
	nam <- paste0("data/geo_data/", gse_id, "_sv_list.RData")
	save(sv_list, file = nam)

	# write out nsv estimation results
	return(nsv_out_dat)
})
nsv_dat <- do.call(rbind, nsv_dat)
write.table(nsv_dat, file = "results/geo_nsv_estimation.txt", quote = F, row.names = F, col.names = T, sep = "\t")
nsv_dat <- read_delim("results/geo_nsv_estimation.txt", delim = "\t")

# plot nsv_dat stuff
p <- ggplot(nsv_dat, aes(x = leek, y = rmt)) +
	geom_point()
ggsave("results/geo_nsv_estimation_plot.pdf", plot = p)


# testing the SVs! 
ncpg_dat <- as.data.frame(matrix(NA, nrow = length(sv_list), ncol = 21))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:21] <- paste0("sv", 1:20, "_adjr2")

geo_ncpg_dat <- lapply(gse$gse_id, function(gse_id) {
	nam <- paste0("data/geo_data/", gse_id, "_sv_list.RData")
	load(nam)
	if (is.na(sv_list)) return(NULL)
	svs <- sv_list[[length(sv_list)]]$sv

	for (i in seq_along(sv_list)) {
		temp_svs <- sv_list[[i]]$sv
		ncpg <- n_cpg[i]
		ncpg_dat[i, "n_cpg"] <- ncpg
		fom <- as.formula(paste0("svs[, j] ~ ", paste(paste0("temp_svs[, ", 1:20, "]"), collapse = " + ")))
		for (j in 1:ncol(temp_svs)) {
			fit <- lm(fom)
			ncpg_dat[i, j+1] <- summary(fit)$adj.r.squared
		}
	}
})

save(geo_ncpg_dat, file = "results/geo_ncpg_dat.txt")





