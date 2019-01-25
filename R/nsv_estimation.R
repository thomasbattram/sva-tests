rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

args <- commandArgs(trailingOnly = TRUE)
ms_dir <- args[1]
wd <- args[2]

# ms_dir <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/"
# wd <- "~/sva_tests/"

setwd(wd)

TP <- "FOM"

#load the samplesheet
load(paste0(ms_dir, "samplesheet/data.Robj"))
samplesheet <- dplyr::filter(samplesheet, time_point == TP) %>%
	dplyr::filter(is.na(duplicate.rm)) %>%
	mutate(ALN = as.numeric(ALN))

# sort the covariate date
covars <- read_delim("~/main_project/ALSPAC_EWAS/methyl_variance/phen/FOM/FOM.qcovar", delim = " ", col_names = F)
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


# ---------------------------------------------------------------
# sort the FOM1 phenotype data
# ---------------------------------------------------------------

# load the FOM1 data


#
# count_na <- function(dat, col_or_row = 2) {
# 	stopifnot(col_or_row %in% c(1,2))
# 	x <- apply(dat, col_or_row, function(x) {sum(is.na(x))})
# 	return(x)
# }
# # remove people
# row_na_num <- count_na(comb_dat, 1)
# summary(row_na_num)
# comb_dat <- comb_dat[row_na_num < (ncol(comb_dat)/10), ]
# dim(comb_dat)

# na_num <- count_na(comb_dat)
# summary(na_num)
# sum(na_num == 0)


# #comp_dat <- comb_dat[, na_num < 100] 
# comp_dat <- comb_dat[, na_num < (ncol(comb_dat)/10)] 
# dim(comp_dat)
# #comp_dat <- comp_dat[,1:200]
# cor_dat <- abs(cor(comp_dat, use = "pairwise.complete.obs")) # STILL LOADS OF NAs!!!!


# ---------------------------------------------------------------
# generate SVs
# ---------------------------------------------------------------

### comparing smartsva and sva package estimation of nsv
## Methylation M values (CpG by Sample)
Y <- matrix(rnorm(20*1000), 1000, 20)
df <- data.frame(pred=gl(2, 10))
## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
## Add one extra dimension to compensate potential loss of 1 degree of freedom
## in confounded scenarios (very important)
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
n.sv # 3

mod = model.matrix(~pred, data=df)
n.sv = num.sv(Y, mod, method = "leek")
n.sv # 0
n.sv = num.sv(Y, mod, method = "be")
n.sv #0

###


mdata <- as.matrix(meth)
mdata <- mdata[complete.cases(mdata), ]
dim(mdata)

## Determine the number of SVs
est_nsv <- function(meth, traits, df) {
	fom <- as.formula(paste0("t(meth)", " ~ ", paste(traits, collapse = " + ")))
	y.r <- t(resid(lm(fom, data = df)))
	print("Estimating number of SVs needed")
	n_sv <- EstDimRMT(y.r, FALSE)$dim + 1
	return(n_sv)
}

### ADD ALN here for real data!!!
df <- data.frame(rcont = rnorm(ncol(mdata)))
for (i in seq(0.05, 0.5, 0.05)) {
	col_nam <- paste0("rbin", i*10)
	df[[col_nam]] <- rbinom(ncol(mdata), 1, i)
}

traits <- colnames(df)
n_sv <- vector(mode = "numeric", length = ncol(df))
names(n_sv) <- colnames(df)

# mod <- model.matrix(~rcont, data = df)
# sva_nsv <- num.sv(mdata, mod, method = "leek")
# sva_nsv_be <- num.sv(mdata, mod, method = "be")
# mdata <- mdata[sample(1:nrow(mdata), 2000),]

for (i in colnames(df)) {
	print(i)
	n_sv[[i]] <- est_nsv(mdata, i, df)
}

df <- cbind(pheno, df)

# generate the SVs
covs <- colnames(pheno)[-c(1, 2)]
# removing slide becuase it has too many unique values (~380...)
covs <- covs[-grep("Slide", covs)]
new_df <- df

summary(pheno)

sva_list <- vector(mode = "list", length = length(traits))
names(sva_list) <- traits

for (i in 1:length(traits)) {
	trait <- traits[i]
	
	fom <- as.formula(paste0(" ~ ", trait, " + ", paste(covs, collapse = " + ")))
	fom0 <- as.formula(paste0(" ~ ", paste(covs, collapse = " + ")))
	mod <- model.matrix(fom, data = df)
	mod0 <- model.matrix(fom0, data = df)


	# mod <- model.matrix(~ rcont, data = df)
	# sva_list[[trait]] <- smartsva.cpp(mdata, mod, mod0=NULL, n.sv = nsv, VERBOSE = T)	


	nsv <- n_sv[[trait]]

	sva_list[[trait]] <- tryCatch({smartsva.cpp(mdata, mod, mod0, n.sv = nsv)},
								   error = function(e) {NULL})
	if (!is.null(svobj)) sva_list[[trait]] <- svobj
}

save(sva_list, file = "sva_list.RData")

# sva_res_list <- list()
# for (i in traits) {
# 	svs <- sva_list[[trait]]$sv
# 	sva_temp <- data.frame()
# 	for (j in seq_along(svs)) {

# 	}
# 	fit <- lm()

# 	svs
# }















