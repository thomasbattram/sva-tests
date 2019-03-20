# ---------------------------------------------------------------
# Testing parameters for SVA
# ---------------------------------------------------------------

rm(list = ls())

pkgs <- c("tidyverse", "sva", "SmartSVA", "gridExtra", "RColorBrewer")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

# ms_dir <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/"
# wd <- "~/sva_tests/"

setwd(wd)

TP="FOM"

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

n_sv <- c(5,10,20)
alpha <- c(0.25, 1)
B <- c(5, 100)
n_cpg <- 1e5
n_samp <- ncol(mdata)
dat_type <- "continuous"
params <- expand.grid(n_sv = n_sv, alpha = alpha, B = B, n_cpg = n_cpg, dat_type = dat_type, n_sample = n_samp)

params$time_user <- NA
params$time_system <- NA
params$time_elapsed <- NA
set.seed(2)
cont <- rnorm(ncol(mdata))
phen <- data.frame(continuous = cont)
i=1
sv_list <- vector(mode = "list", length = nrow(params[params$n_sv == max(params$n_sv),]))
list_num <- 0
for (i in 1:nrow(params)) {
	print(i)
	nsv <- params[i, "n_sv"]
	alph <- params[i, "alpha"]
	b <- params[i, "B"]
	ncpg <- params[i, "n_cpg"]
	dt <- params[i, "dat_type"]
	nsamp <- params[i, "n_sample"]
	set.seed(2)

	if (ncpg != nrow(mdata)) {
		var.idx <- sample(1:nrow(mdata), ncpg)
		mdat <- mdata[var.idx, , drop=F]
	} else {
		mdat <- mdata
	}
	
	fom <- as.formula(paste0("~", dt))
	mod <- model.matrix(fom, data = phen)
	mod0 <- model.matrix(~1, data = phen)
	ptm <- proc.time()

	svobj <- smartsva.cpp(mdat, mod, mod0=NULL, n.sv = nsv, B = b, alpha = alph)

	tim <- proc.time() - ptm
	params[i, "time_user"] <- tim[1]
	params[i, "time_system"] <- tim[2]
	params[i, "time_elapsed"] <- tim[3]

	if (nsv != max(params$n_sv)) {
		next
	} else {
		list_num <- list_num + 1
	## sv_list names follow this order: DT_NSV_SVTYPE_NCPG_NSAMP
	## --- need to shorted some of them
		dt <- substr(dt, 1, 3)
		ncpg <- paste0(round(ncpg/1000), "k")
		nam <- paste(dt, nsv, alph, b, ncpg, nsamp, sep = "_")
		sv_list[[list_num]] <- svobj	
		names(sv_list)[list_num] <- nam
	}
}

# correlation between outputted SVs
# end table
#     sv1 sv2 sv3 ...
# 1&2
# 1&3
# 1&4
# 2&3
# ...

# data.frames in a list
# correlate column 1 of data.frame 1 and 2
# then correlate column 1 of data.frame 1 and 3 etc.


### save it later -- strugs to do it nicely atm!!!!

x <- map(sv_list, "sv") %>%
		map(as.data.frame) 

names(x) <- c("a0.25_B5", "a1_B5", "a0.25_B100", "a1_B100")
combinations <- as.data.frame(combn(names(x), 2))
cols <- seq_along(x[[1]])
out_mat <- matrix(nrow = ncol(combinations), ncol = length(cols) + 1)
out_mat <- as.data.frame(out_mat)
colnames(out_mat) <- c("combination", paste0("sv", ncol - 1))
i=1
for (i in seq_along(combinations)) {
	num1 <- combinations[1, i]
	num2 <- combinations[2, i]
	x1 <- x[[num1]]
	x2 <- x[[num2]]
	out_mat[i, "combination"] <- paste(num1,num2, sep = ", ")
	for (j in seq_along(x1)) {
		out_mat[i, j+1] <- cor(x1[[j]], x2[[j]])
	}
}

write.table(out_mat, file = "results/alpha_tests.txt", row.names = F, col.names = T, quote = F, sep = "\t")

time_plot <- ggplot(params, aes(x = alpha, y = time_user, fill = as.factor(B))) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_grid(. ~ as.factor(n_sv))

ggsave("results/default_params_times_test.pdf", time_plot)











