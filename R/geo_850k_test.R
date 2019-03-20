rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

# ms_dir <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/"
# wd <- "~/sva_tests/"
# phen_dir <- "~/main_project/ALSPAC_EWAS/methyl_variance/phen/"

setwd(wd)

TP <- "FOM"

#load the samplesheet
load(paste0(ms_dir, "samplesheet/data.Robj"))
samplesheet <- dplyr::filter(samplesheet, time_point == TP) %>%
	dplyr::filter(is.na(duplicate.rm)) %>%
	mutate(ALN = as.numeric(ALN))

# sort the covariate date
covars <- read_delim(paste0(phen_dir, "FOM/FOM.qcovar"), delim = " ", col_names = F)
colnames(covars) <- c("ALN", "Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))

pheno <- samplesheet %>%
	left_join(covars) %>%
	dplyr::select(one_of(colnames(covars)), BCD_id, BCD_plate, MSA4Plate_id, Slide) %>%
	dplyr::filter(!is.na(PC1))

#load the methylation data
load(paste0(ms_dir, "betas/data.Robj"))
meth <- beta[, pheno$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(beta)
dim(meth)
# now the 850k dataset!
load("data/gse116339.rda")
meth_850k <- gse116339$data
rm(gse116339)
dim(meth_850k)

# random phenotype
phen <- data.frame(x = rnorm(ncol(meth_850k)))

# function for estimating number of SVs using random matrix theory
est_nsv <- function(meth, traits, df) {
	fom <- as.formula(paste0("t(meth)", " ~ ", paste(traits, collapse = " + ")))
	y.r <- t(resid(lm(fom, data = df)))
	print("Estimating number of SVs needed")
	n_sv <- EstDimRMT(y.r, FALSE)$dim + 1
	return(n_sv)
}

# remove NAs!!!
meth_850k <- na.omit(meth_850k)

rmt_nsv <- est_nsv(as.matrix(meth_850k), "x", phen)

fom <- as.formula("~x")
mod <- model.matrix(fom, data=phen)
leek_nsv <- num.sv(as.matrix(meth_850k), mod, method = "leek")

# compare number of CpGs required for 20 SVs and more
n_cpg <- c(seq(20000, 700000, 20000), nrow(meth_850k))

sv_list <- lapply(n_cpg, function(ncpg) {
	print(ncpg)
	set.seed(2)
	mdat <- meth_850k[sample(1:nrow(meth_850k), ncpg), ]
	svobj <- smartsva.cpp(mdat, mod, mod0 = NULL, n.sv = 20)
	return(svobj)
})

# 850k vs lower 
ncpg_dat <- as.data.frame(matrix(NA, nrow = length(list_nam), ncol = 21))
colnames(ncpg_dat)[1] <- "n_cpg"
colnames(ncpg_dat)[2:21] <- paste0("sv", 1:20, "_adjr2")

svs_850 <- sv_list[[length(sv_list)]]$sv

for (i in seq_along(sv_list)) {
	temp_svs <- sv_list[[i]]$sv
	ncpg <- n_cpg[i]
	ncpg_dat[i, "n_cpg"] <- ncpg
	fom <- as.formula(paste0("svs_850[, j] ~ ", paste(paste0("temp_svs[, ", 1:20, "]"), collapse = " + ")))
	for (j in 1:ncol(svs)) {
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

ggsave("results/ncpg_plot_850k.pdf", plot = ncpg_plot)

write.table(ncpg_dat, file = "results/ncpg_comp_sims_850.txt", quote = F, row.names = F, col.names = T, sep = "\t")







