# ------------------------------------------------------------------------------
# Validation against Leathers et al. (2004)
#
# This is similar to Leathers et al. (2001), but with a bottom at a finite 
# depth.
#

# Check CDF:


#
# Generate input files:
#

winpt <- function(bttd, dir, pre, pos) {
  fnm   <- paste0(dir, pre, "_b", 
    formatC(bttd, width = 5, flag = 0, digits = 2, format = "f"), "_", pos, 
      ".txt") 
  input <- c(
    "sim_nr: 10000000",
    "sim_ns: 3",
    "sim_sza(s): 10.0 20.0 40.0",
    "sim_saa(s): 0.0 0.0 0.0",
    "sim_f0: 1.0",
    "skr_resx: 1",
    "skr_resy: 1",
    "skr_fls:",
    "anc/skyrad/550_continental_clean_80_02_sza_10_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_20_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_40_sky.txt",
    "iop_na: 1.0",
    "iop_nw: 1.34",
    "iop_c: 0.6",
    "iop_nw0: 1",
    "iop_w0: 0.666667",
    "scat_tp: ff",
    "scat_trc: 1",
    "scat_mtd: itp",
    "scat_fbb: 0.0183",
    paste("btt_d:", bttd),
    "btt_tp: lambert",
    "btt_nbr: 1",
    "btt_bhr: 0.2",
    "src_tp: rad",
    "src_fov: 40",
    "src_axs: 0.0 0.0",
    "src_ref_o: 0.0 0.0 0.66",
    "src_rel_o: 0.0 0.0 0.0",
    "acc_geom: none",
    "acc_ext: 1.0",
    "acc_resx: 0.1",
    "acc_resy: 0.1",
    "str_def: 1",
    "str_ncl: 2",
    "0.0 0.0 0.0 0.0 0.0 0.0 0.150 0.12 0.0 0.0 0.0 0.0 1 1",
    "0.0 0.0 0.0 0.0 0.0 0.0 0.044 0.66 0.0 0.0 0.0 0.0 1 1",
    "str_ncn: 0",
    "str_nbx: 0"
  )
  writeLines(input, fnm)
}

bttd <- seq(1, 3, length.out = 5)
for(i in 1:length(bttd)) {
  winpt(bttd = bttd[i], dir = "run/leathers_2004/", 
    pre = "leathers_2004_Lu", pos = "")
}

#
# Run MC code:
#

dir     <- "run/leathers_2004"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  cmd <- paste("./mc_solver_shdw_v1.6.o", file.path("src_v1.6/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

ref <- matrix(NA, ncol = 5, nrow = 3)
ref[1, ] <- c(11.5,  8.2, 4.7, 1.4, 1.3)
ref[2, ] <- c( 1.6, 1.15, 0.7, 0.2, 0.2)
ref[3, ] <- c( 0.1,  0.0, 0.0, 0.0, 0.15)
rownames(ref) <- paste("SZA", c(10, 20, 40), sep = "_")
colnames(ref) <- paste("b", bttd, sep = "_")
ref <- ref/20

ref2 <- matrix(NA, ncol = 5, nrow = 3)
ref2[1, ] <- c(11.5,   4.5,  2.0,  1.4, 1.1)
ref2[2, ] <- c( 2.3,  1.15,  0.9, 0.75, 0.6)
ref2[3, ] <- c( 0.0, 0.025, 0.05,  0.1, 0.1)
rownames(ref2) <- paste("SZA", c(10, 20, 40), sep = "_")
colnames(ref2) <- paste("b", bttd, sep = "_")
ref2 <- ref2/20

res_diff <- res <- matrix(NA, ncol = 5, nrow = 3)

fls <- list.files("run/leathers_2004", "leathers_2004_Lu_b01.00_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
dif_f <- read.table(fls[1], header = TRUE, row.names = 1)[1:3,]
dif_s <- read.table(fls[3], header = TRUE, row.names = 1)[1:3,]
res[, 1] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])
res_diff[, 1] <- as.numeric(((dif_f - dif_s) / dif_f))

fls <- list.files("run/leathers_2004", "leathers_2004_Lu_b01.50_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
dif_f <- read.table(fls[1], header = TRUE, row.names = 1)[1:3,]
dif_s <- read.table(fls[3], header = TRUE, row.names = 1)[1:3,]
res[, 2] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])
res_diff[, 2] <- as.numeric(((dif_f - dif_s) / dif_f))

fls <- list.files("run/leathers_2004", "leathers_2004_Lu_b02.00_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
dif_f <- read.table(fls[1], header = TRUE, row.names = 1)[1:3,]
dif_s <- read.table(fls[3], header = TRUE, row.names = 1)[1:3,]
res[, 3] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])
res_diff[, 3] <- as.numeric(((dif_f - dif_s) / dif_f))

fls <- list.files("run/leathers_2004", "leathers_2004_Lu_b02.50_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
dif_f <- read.table(fls[1], header = TRUE, row.names = 1)[1:3,]
dif_s <- read.table(fls[3], header = TRUE, row.names = 1)[1:3,]
res[, 4] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])
res_diff[, 4] <- as.numeric(((dif_f - dif_s) / dif_f))

fls <- list.files("run/leathers_2004", "leathers_2004_Lu_b03.00_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
dif_f <- read.table(fls[1], header = TRUE, row.names = 1)[1:3,]
dif_s <- read.table(fls[3], header = TRUE, row.names = 1)[1:3,]
res[, 5] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])
res_diff[, 5] <- as.numeric(((dif_f - dif_s) / dif_f))

rownames(res) <- rownames(ref)
colnames(res) <- colnames(ref)

ref
round(res, 3)

#
# Plots:
#

library(magrittr)

tmp_lstat <- function (stats, digits = 3) {
    r2 <- substitute(expression(R^2 == r2), env = list(r2 = round(stats$r2, 
        digits))) %>% eval()
    rmsd <- substitute(expression(RMSD == rmsd ~ units), list(rmsd = round(stats$rmsd, 
        digits+1), units = stats$units)) %>% eval()
    mapd <- substitute(expression(MAPD == mapd ~ "\\%"), list(mapd = round(stats$mapd, 
        digits))) %>% eval()
    bias <- substitute(expression(Bias == bias ~ "\\%"), list(bias = round(stats$bias, 
        digits))) %>% eval()
    n <- substitute(expression(N == n), list(n = stats$n)) %>% 
        eval()
    rang <- substitute(expression(Range == a - b ~ units), list(a = round(stats$rang[1], 
        digits), b = round(stats$rang[2], digits), units = stats$units)) %>% 
        eval()
    ltext <- c(r2, rmsd, mapd, bias, n, rang)
    return(ltext)
}


library(tikzDevice)

fnm <- "val_leathers_2"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
xlab <- "Bottom depth (m)"
ylab <- "$\\epsilon$ (unitless)"
cols <- c("red", "darkgreen", "black") #rainbow(3, start = 0, end = 0.8)
par(mar = c(5, 5, 3, 2))
plot(NA, xlim = c(1, 3), ylim = c(0, 0.6), xlab = xlab, ylab = ylab, cex.lab = 2, 
  cex.axis = 2)
for(i in 1:nrow(ref2)) lines(bttd, ref2[i, ], col = cols[i], lwd = 2)
for(i in 1:nrow(ref)) points(bttd, ref[i, ], col = 1, bg = cols[i], pch = 24)
for(i in 1:nrow(res)) points(bttd, res[i, ], col = 1, bg = cols[i], pch = 21)
legend(0.037, 0.03, tmp_lstat(rho::rstat(as.numeric(ref), as.numeric(res), 
  units = "")), bty = "n", cex = 1.5, y.intersp = 0.6)
sza <- c(10, 20, 40)
legend(1.7, 0.65, sza, title = "$\\theta_{\\text{s}}$ ($^{\\circ}$):", bty = "n", 
  pch = 21, col = 1, pt.bg = cols, cex = 1.5, y.intersp = 0.6, pt.cex = 1)
legend(2.1, 0.63, c("BMC3D", "Analytical", "This study"), bty = "n", 
  pch = c(24, NA, 21), lty = c(NA, 1, NA), col = 1, pt.bg = 1, cex = 1.5, 
  x.intersp = 0.6, y.intersp = 0.6, pt.cex = 1)
legend("right", "\\textbf{D}", bty = "n", cex = 1.6)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.6:
# REF    d_1.0  d_1.5 d_2.0 d_2.5  d_3.0
# SZA_10 0.575 0.4100 0.235  0.07 0.0650
# SZA_20 0.080 0.0575 0.035  0.01 0.0100
# SZA_40 0.005 0.0000 0.000  0.00 0.0075
#
# RES    d_1.0  d_1.5 d_2.0 d_2.5  d_3.0 // BOTH BUG AND APPARETLY WRONG, CONSIDERING PREVIOUS
# SZA_10 0.729 0.3490 0.160    NA 0.0720
# SZA_20 0.349 0.1490 0.089    NA 0.0500
# SZA_40 0.010 0.0060 0.006    NA 0.0070
#
# Benchmark execution time: 6.98, 9.33, 9.66, 11.33, 13.16 min (1E7 rays)

# Version 1.5:
# REF    d_1.0  d_1.5 d_2.0 d_2.5  d_3.0
# SZA_10 0.575 0.4100 0.235  0.07 0.0650
# SZA_20 0.080 0.0575 0.035  0.01 0.0100
# SZA_40 0.005 0.0000 0.000  0.00 0.0075
#
# RES    d_1.0  d_1.5 d_2.0 d_2.5  d_3.0
# SZA_10 0.528  0.205 0.098 0.061  0.047
# SZA_20 0.135  0.062 0.043 0.033  0.028
# SZA_40 0.007  0.005 0.005 0.005  0.006
#
# Benchmark execution time: 0.40, 0.51, 0.60, 0.70, 0.72 min (1E6 rays)


