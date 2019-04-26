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
g_nsv_dat <- nsv_dat %>%
	gather(key = "method", value = "nsv", -gse_id, -n_sample, -n_cpg)

p <- ggplot(g_nsv_dat, aes(x = gse_id, y = nsv, fill = method)) +
	geom_bar(stat = "identity", position = "dodge") +
	coord_flip() +
	theme(axis.text.y = element_blank())

ggsave("results/geo_nsv_estimation_plot.pdf", plot = p)

# testing the SVs! 
gse_id <- "gse110607"
geo_ncpg_dat <- lapply(gses$gse_id, function(gse_id) {
	print(gse_id)
	# load the data
	nam <- paste0("data/geo_data/", gse_id, "_sv_list.RData")
	ifelse(file.exists(nam), load(nam), return(NULL))
	if (all(is.na(sv_list))) return(NULL)
	# make table
	ncpg_dat <- as.data.frame(matrix(NA, nrow = length(sv_list), ncol = 21))
	colnames(ncpg_dat)[1] <- "n_cpg"
	colnames(ncpg_dat)[2:21] <- paste0("sv", 1:20, "_adjr2")

	# setting ncpg
	if (nrow(ncpg_dat) == 16) {
		n_cpg <- c(seq(20000, 300000, 20000), nsv_dat[nsv_dat$gse_id == gse_id, "n_cpg", drop = T])
		# n_cpg <- c(seq(20000, 300000, 20000), 480000)
	} else {
		n_cpg <- c(seq(20000, 700000, 20000), nsv_dat[nsv_dat$gse_id == gse_id, "n_cpg"])
	}

	# extract SVs from all CpG sites
	svs <- sv_list[[length(sv_list)]]$sv

	# Assess variance of SVs, taken above, explained by SVs produced from subsets of CpG sites
	for (i in seq_along(sv_list)) {
		ncpg <- n_cpg[i]
		ncpg_dat[i, "n_cpg"] <- ncpg
		if (is.na(sv_list[[i]])) {
			ncpg_dat[i, 2:nrow(ncpg_dat)] <- NA
			next
		} else {
			temp_svs <- sv_list[[i]]$sv
		}
		fom <- as.formula(paste0("svs[, j] ~ ", paste(paste0("temp_svs[, ", 1:20, "]"), collapse = " + ")))
		for (j in 1:ncol(temp_svs)) {
			fit <- lm(fom)
			ncpg_dat[i, j+1] <- summary(fit)$adj.r.squared
		}
	}
	ncpg_dat$gse_id <- gse_id
	return(ncpg_dat)
})
fin_dat <- geo_ncpg_dat[!map_lgl(geo_ncpg_dat, is.null)] %>%
	do.call(rbind, .)
rownames(fin_dat) <- NULL

g_dat <- fin_dat %>%
	gather("sv", "adj_r2", -n_cpg, -gse_id)
g_dat$sv <- gsub("sv", "", g_dat$sv)
g_dat$sv <- gsub("_adjr2", "", g_dat$sv)

# split data into epic and non-epic arrays
g_dat_epic <- g_dat %>%
	dplyr::filter(gse_id == "gse116339")
g_dat <- g_dat %>%
	dplyr::filter(gse_id != "gse116339")
### issue with n_cpg factors!!!!
plot_dat <- function(dat, sv_n) {
	pdat <- dplyr::filter(dat, sv == sv_n) %>%
		dplyr::filter(n_cpg <= 3e5)
	p <- ggplot(pdat, aes(x = as.factor(n_cpg), y = adj_r2)) +
		geom_boxplot() +
		geom_point(aes(colour = gse_id)) +
		geom_jitter() +
		theme(legend.position = "none") +
		ggtitle(sv_n)
	return(p)
}

p <- lapply(unique(g_dat$sv), function(sv_n) {
	x <- plot_dat(g_dat, sv_n)
	return(x)
})
ggsave("results/geo_ncpg_plot.pdf", plot = marrangeGrob(p, nrow=1, ncol=1))

p_epic <- ggplot(g_dat_epic, aes(x = n_cpg, y = adj_r2,
								 colour = reorder(sv, sort(as.numeric(sv))), group = sv)) +
	geom_point() +
	geom_line() +
	scale_colour_discrete(name = "SV")

ggsave("results/geo_epic_ncpg_plot.pdf", plot = p_epic)


save(geo_ncpg_dat, file = "results/geo_ncpg_dat.RData")





