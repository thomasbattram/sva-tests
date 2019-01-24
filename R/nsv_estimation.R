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
	dplyr::filter(is.na(duplicate.rm))

#load the methylation data
load(paste0(ms_dir, "betas/data.Robj"))
meth <- beta[, samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(beta)
dim(meth)

mdata <- as.matrix(meth)
mdata <- mdata[complete.cases(mdata), ]
dim(mdata)

#
# sort the FOM1 phenotype data
#

# load the FOM1 data


#
count_na <- function(dat, col_or_row = 2) {
	stopifnot(col_or_row %in% c(1,2))
	x <- apply(dat, col_or_row, function(x) {sum(is.na(x))})
	return(x)
}
# remove people
row_na_num <- count_na(comb_dat, 1)
summary(row_na_num)
comb_dat <- comb_dat[row_na_num < (ncol(comb_dat)/10), ]
dim(comb_dat)

na_num <- count_na(comb_dat)
summary(na_num)
sum(na_num == 0)


#comp_dat <- comb_dat[, na_num < 100] 
comp_dat <- comb_dat[, na_num < (ncol(comb_dat)/10)] 
dim(comp_dat)
#comp_dat <- comp_dat[,1:200]
cor_dat <- abs(cor(comp_dat, use = "pairwise.complete.obs")) # STILL LOADS OF NAs!!!!


# ---------------------------------------------------------------
# generate SVs
# ---------------------------------------------------------------
## Determine the number of SVs
est_nsv <- function(meth, trait, df) {
	fom <- as.formula(paste0("t(", meth, ")", " ~ ", trait))
	y.r <- t(resid(lm(fom, data = get(df))))
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

n_sv <- vector(mode = "numeric", length = ncol(df))
for (i in colnames(df)) {
	print(i)
	n_sv[[i]] <- est_nsv(mdata, i, df)
}

# generate the SVs
covars <- 
traits <- 
sva_list <- list()
for (i in 1:length(traits)) {
	trait <- traits[i]
	fom <- as.formula(paste0(trait, " ~ ", paste(covars, collapse = " + ")))
	mod <- model.matrix(fom, data = df)
	mod0 <- model.matrix(~1, data = df)

	sva_list[[i]] <- smartsva.cpp(mdat, mod, mod0=NULL, n.sv = nsv, B = 5)	
}

















