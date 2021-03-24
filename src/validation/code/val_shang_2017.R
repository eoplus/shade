# ------------------------------------------------------------------------------
# Validation against Shang et al. (2017)

#
# Generate input files:
#

winpt <- function(c, dir, pre, pos) {
  fnm   <- paste0(dir, pre, "_a", 
    formatC(c * 1 / 4, width = 5, flag = 0, digits = 2, format = "f"), "_", pos, 
      ".txt") 
  input <- c(
    "sim_nr: 1000000",
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
    "iop_w0: 0.75",
    "scat_tp: ff",
    "scat_trc: 1",
    "scat_mtd: itp",
    "scat_fbb: 0.016",
    "btt_d: INFINITY",
    "src_tp: rad",
    "src_fov: 10.0",
    "src_axs: 0.0 0.0",
    "src_ref_o: 0.0 0.0 -0.025",
    "src_rel_o: 0.0 0.0 0.0",
    "acc_geom: none",
    "acc_ext: 1.0",
    "acc_resx: 0.1",
    "acc_resy: 0.1",
    "str_def: 1",
    "str_ncl: 0",
    "str_ncn: 1",
    "0.0 0.0 -0.11 0.0 0.0 0.0 18.21 0.08 0.16 0.0 0.0 0.0 0.0 1 0",
    "str_nbx: 0"
  )
  writeLines(input, fnm)
}

iop_a <- c(0.05, 0.1, 0.2, 0.5)
iop_c <- iop_a + 3 * iop_a 

for(i in 1:length(iop_c)) {
  winpt(c = iop_c[i], dir = "run/shang_2017/", 
    pre = "shang_2017_Lu", pos = "")
}

#
# Run MC code:
#

dir     <- "run/shang_2017"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  if(file.exists(file.path("run/shang_2017", gsub(".txt", "_out_summary.txt", fls[i])))) next
  cmd <- paste("./mc_solver_shdw_v1.6.o", file.path("src/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

ref <- matrix(NA, ncol = 4, nrow = 1)
ref[1, ] <- c(0.018, 0.032, 0.06, 0.12)
rownames(ref) <- paste("bb_a", c(0.05), sep = "_")
colnames(ref) <- paste("a", iop_a, sep = "_")

res <- matrix(NA, ncol = 4, nrow = 1)

fls   <- list.files("run/shang_2017", "shang_2017_Lu_a00.05_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[1, 1] <- as.numeric(((dir_f - dir_s) / dir_f)[4, 1])

fls <- list.files("run/shang_2017", "shang_2017_Lu_a00.10_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[1, 2] <- as.numeric(((dir_f - dir_s) / dir_f)[4, 1])

fls <- list.files("run/shang_2017", "shang_2017_Lu_a00.20_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[1, 3] <- as.numeric(((dir_f - dir_s) / dir_f)[4, 1])

fls <- list.files("run/shang_2017", "shang_2017_Lu_a00.50_", full.names = T)
dir_f <- read.table(fls[5], header = TRUE, row.names = 1)
dir_s <- read.table(fls[7], header = TRUE, row.names = 1)
res[1, 4] <- as.numeric(((dir_f - dir_s) / dir_f)[4, 1])


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
# SZA_0   0.094 0.168 0.283 0.458 0.566
# SZA_10  0.069 0.123 0.209 0.352 0.467
# SZA_20  0.020 0.040 0.074 0.163 0.287
# SZA_30  0.011 0.024 0.047 0.106 0.214
# SZA_40  0.008 0.017 0.037 0.084 0.160
# SZA_50  0.006 0.015 0.028 0.068 0.133
# SZA_60  0.005 0.013 0.023 0.057 0.113
# SZA_70  0.005 0.012 0.023 0.048 0.103

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





