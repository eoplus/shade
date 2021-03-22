# ------------------------------------------------------------------------------
# Validation against Leathers et al. (2001)
#
# Their study provided good tabulated values for comparison. Their instrument is 
# the first model of the Hyper-TRSB which essentially is two concentric 
# cylinders. The top one is a buoy of 15 cm radius and with bottom part at 12 cm
# depth and the instrument is a cylinder of 4.4 cm with the objective at 66 cm 
# depth. No above water structures were modelled. The FOV reporting is not 
# consistent. Base on their image of the CDF, it appears that a 28 degree FOV 
# (14 degrees half FOV) is the appropriate value. Velow is the CDF of a uniform
# probability function between 0 and 14 degrees. It compares well with their 
# CDF, though their CDF emits less photons at small angles and that may result 
# in slightly larger shading estimated errors in my simulations. 
#
# The comparison here is made against his set of black sky with b = 2*a. The 
# diffuse compoenent was reportedly computed with an analytical sky radiance 
# model at 480 nm and SZA of 30 degrees. It is compared here against 550 nm at 
# 30 degrees and totally diffuse.
#
# Though they provided values for a shallow bottom, those values were wrong and
# corrected in a later paper.

# Check CDF:

par(mar = c(5, 5, 3, 2))
x1 <- seq(0, 14, length.out = 100)
x2 <- seq(14, 25, length.out = 100)
plot(c(x1, x2), c(x1/14, x2*0+1), type = "l", col = 2, xaxs = "i", yaxs = "i", 
  xlab = "Angle", ylab = "CDF")
abline(v = c(5, 10, 15, 20), lty = 2, col = "grey")
abline(h = c(0.2, 0.4, 0.6, 0.8), lty = 2, col = "grey")

#
# Generate input files:
#

winpt <- function(c, dir, pre, pos) {
  fnm   <- paste0(dir, pre, "_a", 
    formatC(c * 1 / 3, width = 5, flag = 0, digits = 2, format = "f"), "_", pos, 
      ".txt") 
  input <- c(
    "sim_nr: 10000000",
    "sim_ns: 8",
    "sim_sza(s): 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0",
    "sim_saa(s): 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0",
    "sim_f0: 1.0",
    "skr_resx: 1",
    "skr_resy: 1",
    "skr_fls:",
    "anc/skyrad/550_continental_clean_80_02_sza_00_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_10_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_20_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_30_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_40_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_50_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt",
    "anc/skyrad/550_continental_clean_80_02_sza_70_sky.txt",
    "iop_na: 1.0",
    "iop_nw: 1.34",
    paste("iop_c:", c),
    "iop_nw0: 1",
    "iop_w0: 0.666667",
    "scat_tp: ff",
    "scat_trc: 1",
    "scat_mtd: itp",
    "scat_fbb: 0.0183",
    "btt_d: INFINITY",
    "src_tp: rad",
    "src_fov: 28.0",
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

iop_a <- c(0.05, 0.1, 0.2, 0.5, 1.0)
iop_c <- iop_a + 2* iop_a 

for(i in 1:length(iop_c)) {
  winpt(c = iop_c[i], dir = "run/leathers_2001/", 
    pre = "leathers_2001_Lu", pos = "")
}

#
# Run MC code:
#

dir     <- "run/leathers_2001"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  if(file.exists(file.path("run/leathers_2001", gsub(".txt", "_out_summary.txt", fls[i])))) next
  cmd <- paste("./mc_solver_shdw_v1.6.o", file.path("src/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

ref <- matrix(NA, ncol = 5, nrow = 8)
ref[1, ] <- c(11.1, 20.8, 33.3, 49.3, 58.1)
ref[2, ] <- c(6.9, 12.3, 20.4, 35.0, 46.5)
ref[3, ] <- c(2.2, 4.1, 7.6, 17.3, 28.8)
ref[4, ] <- c(1.2, 2.2, 4.9, 11.3, 20.9)
ref[5, ] <- c(0.7, 1.8, 3.5, 8.8, 16.4)
ref[6, ] <- c(0.6, 1.1, 2.8, 6.9, 13.9)
ref[7, ] <- c(0.6, 1.0, 2.4, 6.1, 12.1)
ref[8, ] <- c(0.7, 1.1, 2.3, 5.3, 11.2)
rownames(ref) <- paste("SZA", seq(0, 70, 10), sep = "_")
colnames(ref) <- paste("a", iop_a, sep = "_")
ref <- ref/100

res <- matrix(NA, ncol = 5, nrow = 8)

fls <- list.files("run/leathers_2001", "leathers_2001_Lu_a01.00_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[, 5] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])

fls <- list.files("run/leathers_2001", "leathers_2001_Lu_a00.50_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[, 4] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])

fls <- list.files("run/leathers_2001", "leathers_2001_Lu_a00.20_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[, 3] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])

fls <- list.files("run/leathers_2001", "leathers_2001_Lu_a00.10_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[, 2] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])

fls <- list.files("run/leathers_2001", "leathers_2001_Lu_a00.05_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[, 1] <- as.numeric(((dir_f - dir_s) / dir_f)[, 1])

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

fnm <- "val_leathers_1"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
lim  <- c(0.003, 0.80)
xlab <- "Leathers et al. (2001) $\\epsilon$ (unitless)"
ylab <- "This study $\\epsilon$ (unitless)"
col  <- rep(rainbow(8, start = 0, end = 0.8), 5)
pch  <- rep(21:25, each = 8)
par(mar = c(5, 5, 3, 2))
plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, cex.lab = 2, 
  cex.axis = 2, log = "xy")
points(ref, res, pch = pch, col = 1, bg = col)
abline(a = 0, b = 1, lwd = 2)
legend(0.037, 0.03, tmp_lstat(rho::rstat(as.numeric(ref), as.numeric(res), 
  units = "")), bty = "n", cex = 1.5, y.intersp = 0.6)
sza <- seq(0, 70, 10)
legend(0.0025, 1.09, sza, title = "$\\theta_{\\text{s}}$ ($^{\\circ}$):", bty = "n", 
  pch = 21, col = 1, pt.bg = rainbow(8, start = 0, end = 0.8), cex = 1.5, 
  y.intersp = 0.6, pt.cex = 1)
vec <- iop_c
legend(0.008, 1.09, vec, title = "$c_{\\text{t}}$ (m$^{-1}$):", bty = "n", 
  pch = 21:25, col = 1, pt.bg = 1, cex = 1.5, 
  y.intersp = 0.6, pt.cex = 1)
legend("right", "\\textbf{B}", bty = "n", cex = 1.6)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.6:
# REF    a_0.05 a_0.1 a_0.2 a_0.5   a_1
# SZA_0   0.111 0.208 0.333 0.493 0.581
# SZA_10  0.069 0.123 0.204 0.350 0.465
# SZA_20  0.022 0.041 0.076 0.173 0.288
# SZA_30  0.012 0.022 0.049 0.113 0.209
# SZA_40  0.007 0.018 0.035 0.088 0.164
# SZA_50  0.006 0.011 0.028 0.069 0.139
# SZA_60  0.006 0.010 0.024 0.061 0.121
# SZA_70  0.007 0.011 0.023 0.053 0.112
#
# RES    a_0.05 a_0.1 a_0.2 a_0.5   a_1
# SZA_0   0.111 0.198 0.329 0.509 0.615
# SZA_10  0.098 0.173 0.286 0.449 0.571
# SZA_20  0.038 0.071 0.132 0.261 0.406
# SZA_30  0.012 0.025 0.050 0.128 0.239
# SZA_40  0.009 0.018 0.034 0.092 0.171
# SZA_50  0.006 0.015 0.025 0.078 0.143
# SZA_60  0.005 0.012 0.023 0.065 0.123
# SZA_70  0.005 0.011 0.021 0.061 0.117
#
# Benchmark execution time: 27.28, 27.46, 25.71, 25.38, 28.85 min (1E7 rays)

# Version 1.5:
# REF    a_0.05 a_0.1 a_0.2 a_0.5   a_1
# SZA_0   0.111 0.208 0.333 0.493 0.581
# SZA_10  0.069 0.123 0.204 0.350 0.465
# SZA_20  0.022 0.041 0.076 0.173 0.288
# SZA_30  0.012 0.022 0.049 0.113 0.209
# SZA_40  0.007 0.018 0.035 0.088 0.164
# SZA_50  0.006 0.011 0.028 0.069 0.139
# SZA_60  0.006 0.010 0.024 0.061 0.121
# SZA_70  0.007 0.011 0.023 0.053 0.112
#
# RES    a_0.05 a_0.1 a_0.2 a_0.5   a_1
# SZA_0   0.095 0.171 0.283 0.454 0.568
# SZA_10  0.068 0.122 0.209 0.351 0.480
# SZA_20  0.020 0.037 0.077 0.171 0.282
# SZA_30  0.012 0.021 0.046 0.109 0.205
# SZA_40  0.008 0.016 0.034 0.082 0.165
# SZA_50  0.007 0.013 0.029 0.073 0.135
# SZA_60  0.006 0.011 0.022 0.061 0.113
# SZA_70  0.005 0.010 0.021 0.051 0.106
#
# Benchmark execution time: 2.83, 2.75, 2.61, 2.61, 2.63 min (1E6 rays)





