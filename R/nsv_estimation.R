rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

args <- commandArgs(trailingOnly = TRUE)
ms_dir <- args[1]
wd <- args[2]
phen_dir <- args[3]

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

# ---------------------------------------------------------------
# testing estimating the number of SVs
# ---------------------------------------------------------------

nsamp <- seq(100, 1000, 100)
ncpg <- c(seq(10000, 100000, 10000), nrow(meth))
params <- expand.grid(nsamp = nsamp, ncpg = ncpg)
params$rmt_nsv <- NA
params$leek_nsv <- NA 

for (i in 1:nrow(params)) {
	print(i)
	t_nsamp <- params[i, "nsamp"]
	t_ncpg <- params[i, "ncpg"]
	## Methylation M values (CpG by Sample)
	Y <- matrix(rnorm(t_nsamp*t_ncpg), t_ncpg, t_nsamp)
	df <- data.frame(pred=gl(2, t_nsamp/2))
	## Determine the number of SVs
	Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
	## Add one extra dimension to compensate potential loss of 1 degree of freedom
	## in confounded scenarios (very important)
	params[i, "rmt_nsv"] <- EstDimRMT(Y.r, FALSE)$dim + 1

	mod = model.matrix(~pred, data=df)
	params[i, "leek_nsv"] = num.sv(Y, mod, method = "leek")
}

write.table(params, "data/sim_estimated_sv_num.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# plot it 
g_params <- gather(params, key = "method", value = "nsv", -ncpg, -nsamp) %>%
	mutate(method = gsub("_nsv", "", method))
p <- ggplot(g_params, aes(x = ncpg, y = nsv, colour = as.factor(nsamp))) +
	geom_point() +
	geom_line(aes(linetype = method)) + 
	scale_colour_discrete(name = "n_sample")

ggsave("results/estimating_nsv_methods.pdf", plot = p)

# ---------------------------------------------------------------
# sort the FOM1 phenotype data
# ---------------------------------------------------------------

# load the FOM1 data
fom_dat <- read_delim(paste0(phen_dir, "ALSPAC/FOM1_only_data_FOM.txt"), delim = "\t")


count_na <- function(dat, col_or_row = 2) {
	stopifnot(col_or_row %in% c(1,2))
	x <- apply(dat, col_or_row, function(x) {sum(is.na(x))})
	return(x)
}
# remove people
test_dat <- fom_dat %>%
	dplyr::select(-aln, -alnqlet)
row_na_num <- count_na(test_dat, 1)
summary(row_na_num)
test_dat <- test_dat[row_na_num < (ncol(test_dat)/10), ]
dim(test_dat)

na_num <- count_na(test_dat)
summary(na_num)
sum(na_num == 0)


#comp_dat <- comb_dat[, na_num < 100] 
comp_dat <- test_dat[, na_num < (ncol(test_dat)/10)] 
dim(comp_dat)
#comp_dat <- comp_dat[,1:200]
cor_dat <- abs(cor(comp_dat, use = "pairwise.complete.obs")) 
cor_dat[1:10, 1:10]

# Extract traits that correlate with many others and remove them
extract_cor <- function(dat, cutoff = 0.2) {
	x <- apply(dat, 2, function(x) {sum(x > cutoff, na.rm = T)})
	y <- x[order(x, decreasing = T)]
	return(y)
}

cor_num <- extract_cor(cor_dat)

cor_num[1]
rm_list <- vector(mode = "character")
rm_num <- 1
temp_dat <- cor_dat
# loop takes the pheno correlated with the most variables removes it then starts again
while(cor_num[1] > 1) {
	print(cor_num[1])
	rm_list[[rm_num]] <- names(cor_num[1])
	rm_num <- rm_num + 1
	temp_dat <- temp_dat[!(rownames(temp_dat) %in% names(cor_num[1])), !(colnames(temp_dat) %in% names(cor_num[1]))]
	cor_num <- extract_cor(temp_dat)
}
length(colnames(temp_dat)) # 36
# remove "blood sample room"
temp_dat <- temp_dat[-grep("Blood_sample_room", colnames(temp_dat)), -grep("Blood_sample_room", colnames(temp_dat))]

# take 10 random traits from the 35
set.seed(2)
real_dat <- fom_dat %>%
	dplyr::select(aln, one_of(sample(colnames(temp_dat), 10))) %>%
	dplyr::filter(aln %in% pheno$ALN)

colnames(real_dat) <- gsub(" ", "_", colnames(real_dat))
colnames(real_dat) <- gsub("\\%", "percent", colnames(real_dat))
colnames(real_dat) <- gsub("[[:punct:]]", "_", colnames(real_dat))

# ---------------------------------------------------------------
# generate SVs
# ---------------------------------------------------------------

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
df <- data.frame(aln = real_dat$aln, rcont = rnorm(ncol(mdata)))
for (i in seq(0.05, 0.5, 0.05)) {
	col_nam <- paste0("rbin", i*10)
	df[[col_nam]] <- rbinom(ncol(mdata), 1, i)
}
df <- df %>%
	left_join(real_dat)

traits <- colnames(df)[-grep("aln", colnames(df))]
nsv_dat <- data.frame(trait = traits, rmt = NA, leek_nsv = NA)

df <- df %>%
	left_join(pheno, by = c("aln" = "ALN"))

covs <- colnames(pheno)[-c(1, 2)]
# removing slide becuase it has too many unique values (~380...) 
# + removing BCD_plate because it doesn't work when trying to figure out number of SVs needed...
covs <- covs[-grep(c("Slide|BCD_plate"), covs)]

# temp_mdata <- mdata[sample(1:nrow(mdata), 2000), comp_df$Sample_Name]
for (i in 1:nrow(nsv_dat)) {
	print(i)
	temp_trait <- as.character(nsv_dat[i, "trait"])

	temp_df <- df %>%
		dplyr::select(Sample_Name, one_of(c(temp_trait, covs))) %>%
		na.omit()

	temp_mdata <- mdata[, temp_df$Sample_Name]

	nsv_dat[i, "rmt"] <- tryCatch({est_nsv(temp_mdata, temp_trait, temp_df)}, 
						 error = function(e) {err_msg(e)})

	fom <- as.formula(paste("~", paste(temp_trait, collapse = " + ")))
	mod <- model.matrix(fom, data=temp_df)
	nsv_dat[i, "leek_nsv"] <- tryCatch({num.sv(temp_mdata, mod, method = "leek")}, 
						  error = function(e) {err_msg(e)})

}

write.table(nsv_dat, "data/estimated_sv_num.txt", col.names = T, row.names = F, quote = F, sep = "\t")


# generate the SVs
sva_list <- vector(mode = "list", length = length(traits))
names(sva_list) <- traits
i=1
for (i in 1:length(traits)) {
	trait <- traits[i]
	
	temp_df <- df %>%
		dplyr::select(Sample_Name, one_of(c(trait, covs))) %>%
		na.omit()

	temp_mdata <- mdata[, temp_df$Sample_Name]
	fom <- as.formula(paste0(" ~ ", trait))
	# fom0 <- as.formula(paste0(" ~ ",))
	mod <- model.matrix(fom, data = temp_df)
	# mod0 <- model.matrix(fom0, data = temp_df)

	nsv <- max(unlist(nsv_dat[nsv_dat$trait == trait, c(2,3), drop = T]))

	svobj <- tryCatch({smartsva.cpp(temp_mdata, mod, mod0=NULL, n.sv = nsv)},
								   error = function(e) {err_msg(e = e, return = NULL)})
	if (!is.null(svobj)) sva_list[[trait]] <- svobj
}

# test the SVs
sva_res_list <- list()
i=traits[1]
for (i in traits) {
	print(i)
	temp_df <- df %>%
		dplyr::select(Sample_Name, one_of(c(i, covs))) %>%
		na.omit()
	
	svs <- sva_list[[i]]$sv
	if (is.null(svs)) next
	sva_temp <- as.data.frame(matrix(NA, nrow = length(covs), ncol = ncol(svs) + 1))
	colnames(sva_temp) <- c("covariate", paste0("sv", 1:ncol(svs)))
	sva_temp$covariate <- covs

	res_list <- list()
	for (j in 1:ncol(svs)) {
		print(j)
		temp_svs <- as.data.frame(svs[, 1:j, drop = F])
		
		res_list[[j]] <- apply(temp_df[, covs], 2, function(x) {summary(lm(x ~ ., data = temp_svs))$adj.r.squared})
		
	}
	res <- as.data.frame(do.call(rbind, res_list)) %>%
		mutate(sv = as.factor(1:nrow(.))) %>%
		dplyr::select(sv, everything())

	sva_res_list[[i]] <- res
}

save(sva_res_list, file = "data/cov_r2_res.RData")

# plot it all! 
plot_list <- list()
for (i in names(sva_res_list)) {
	g_res <- gather(sva_res_list[[i]], key = "covariate", value = "adjusted_r2", -sv)
	plot_list[[i]] <- ggplot(g_res, aes(x = sv, y = adjusted_r2, colour = covariate, group = covariate)) +
		geom_point() + 
		geom_line() + 
		labs(title = i)
}

pdf("results/covs_variance_explained.pdf", width = 12, height = 10)
marrangeGrob(plot_list, ncol = 1, nrow = 2)
dev.off()


# all_res <- do.call(rbind, sva_res_list) 
# rownames(all_res) <- NULL

# g_all_res <- all_res %>%
# 	gather(key = "covariate", value = "adjusted_r2", -sv)

# p <- ggplot(g_all_res, aes(x = covariate, y = adjusted_r2)) +
# 	geom_violin()

# p_violin_list <- list()
# for (i in seq(10, max(as.numeric(g_all_res$sv)), 10)) {
# 	temp_res <- g_all_res %>%
# 		dplyr::filter(sv == i)
# 	p_violin_list[[as.character(i)]] <- ggplot(temp_res, aes(x = covariate, y = adjusted_r2)) +
# 		geom_violin()	
# }

# pdf("results/covs_r2_violin.pdf", width = 15, height = 10)
# marrangeGrob(p_violin_list, nrow = 2, ncol = 1)
# dev.off()

# ggsave("results/covs_r2_violin.pdf", plot = p, width = 15, height = 10, units = "in")








