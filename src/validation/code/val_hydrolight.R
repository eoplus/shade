# ------------------------------------------------------------------------------
# Validation against HYDROLIGHT code version 4.2
#
# HYDROLIGHT was used to solve the RTE for a single wavelength. The absorption
# coefficient was a = {0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0}, with the scattering 
# coefficient being b = 2*a (single scattering albedo of 2/3). The phase 
# function was taken to be the Fournier-Forand with a backscattering ratio of 
# 0.0183 (equivalent to average Petzold). Inelastic scattering was not included.
# The sky radiance was set as idealized with diffuse component described as 
# uniform and the diffuse fraction of 30%, with a SZA of 30 degrees.
#
# Two sets of simulations were performed: infinite bottom and bottom at 2 m 
# depth with a lambertian bi-hemispherical reflectance of 0.2. Since when the 
# bottom is at an infinite distance the reflectance at a given optical depth 
# (surface in this case) depends only on the single scattering albedo and phase 
# function, the only absorption value used was 0.1.
#
# The program was kept with its 10 degree FOV polar caps and a radiance sensor 
# with a 10 degree FOV will be used for comparison. The refractive index of air 
# of 1 and water of 1.34 were used.
#
# Results Version 1.5:
# The most striking feature of the comparison is a slight inclination on the 
# relation given by -0.0001761 + 1.0541773x, i.e., an average ~5% higher values 
# are observed in my simulations. The MAPD was 3.166 % and the bias was 0.124 %.
# In optically deep waters, the difference was 3.49 % (underestimation) going to 
# -5.30 % (overestimation) at the shallowest depth.

# ------------------------------------------------------------------------------
# water-leaving hemispherical-directional reflectance (1/sr):

#
# Generate input files:
#

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 30.0",
  "sim_saa(s): 0.0",
  "sim_f0: 1.0",
  "skr_resx: 1",
  "skr_resy: 1",
  "skr_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_30_sky.txt", 
  "iop_na: 1.0",
  "iop_nw: 1.34",
  "iop_c: 1.0",
  "iop_nw0: 9",
  "iop_w0: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9",
  "scat_tp: ff",
  "scat_trc: 1",
  "scat_mtd: itp",
  "scat_fbb: 0.0183",
  "btt_d: INFINITY",
  "src_tp: rad",
  "src_fov: 10.0",
  "src_axs: 0.0 0.0",
  "src_ref_o: 0.0 0.0 -0.0",
  "src_rel_o: 0.0 0.0 0.0",
  "acc_geom: none",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)
writeLines(input, "run/hydrolight/val_Lu_inf.txt")

dir <- "run/hydrolight/"
rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
cmd <- paste("./mc_solver_free_v1.6.o", "src_v1.6/validation/code/run/hydrolight/val_Lu_inf.txt")
system(cmd)

setwd(rd)

winpt <- function(c, dir, pre, pos) {
  fnm   <- paste0(dir, pre, "_a", 
    formatC(c * 1 / 3, width = 5, flag = 0, digits = 2, format = "f"), "_", pos, 
      ".txt") 
  input <- c(
    "sim_nr: 10000000",
    "sim_ns: 1",
    "sim_sza(s): 30.0",
    "sim_saa(s): 0.0",
    "sim_f0: 1.0",
    "skr_resx: 1",
    "skr_resy: 1",
    "skr_fls:",
    "anc/skyrad/550_continental_clean_80_02_sza_30_sky.txt", 
    "iop_na: 1.0",
    "iop_nw: 1.34",
    paste("iop_c:", c),
    "iop_nw0: 1",
    "iop_w0: 0.666667",
    "scat_tp: ff",
    "scat_trc: 1",
    "scat_mtd: itp",
    "scat_fbb: 0.0183",
    "btt_d: 2.0",
    "btt_tp: lambert",
    "btt_nbr: 1",
    "btt_bhr: 0.2",
    "src_tp: rad",
    "src_fov: 10.0",
    "src_axs: 0.0 0.0",
    "src_ref_o: 0.0 0.0 -0.0",
    "src_rel_o: 0.0 0.0 0.0",
    "acc_geom: none",
    "acc_ext: 1.0",
    "acc_resx: 0.1",
    "acc_resy: 0.1",
    "str_def: 0"
  )
  writeLines(input, fnm)
}

iop_a <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
iop_c <- iop_a + 2* iop_a 

for(i in 1:length(iop_c)) {
  winpt(c = iop_c[i], dir = "run/hydrolight/", pre = "val_Lu", pos = "btt")
}

dir     <- "run/hydrolight"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)[1:7]

dir <- "run/hydrolight/"
rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  cmd <- paste("./mc_solver_free_v1.6.o", file.path("src_v1.6/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

library(magrittr)

fls_h42    <- c(
  list.files("run/hydrolight/h42", "Mmono_btt", full.names = T),
  list.files("run/hydrolight/h42", "Mmono_inf", full.names = T)
)

h42_Lu_btt <- data.frame(
  iop_c      = iop_c, 
  bttd       = 2, 
  HL_Rrs_dif = numeric(length(iop_c)), 
  HL_Rrs_dir = numeric(length(iop_c)),
  MC_Rrs_dif = numeric(length(iop_c)), 
  MC_Rrs_dir = numeric(length(iop_c))
)
j = 1
for(i in seq(1, 14, 2)) {
  h42_Lu_btt$HL_Rrs_dif[j] <- readLines(fls_h42[i]) %>%
                              '['(24) %>%
                              strsplit(" ") %>%
                              unlist() %>%
                              as.numeric()%>%
                              '['(4)
 j = j+1
}
j = 1
for(i in seq(2, 14, 2)) {
  h42_Lu_btt$HL_Rrs_dir[j] <- readLines(fls_h42[i]) %>%
                              '['(24) %>%
                              strsplit(" ") %>%
                              unlist() %>%
                              as.numeric()%>%
                              '['(4)
 j = j+1
}

fls_mc_dif_mn <- list.files("run/hydrolight", "btt_out_dif_f_mn", full.names = T)
fls_mc_dir_mn <- list.files("run/hydrolight", "btt_out_dir_f_mn", full.names = T)

for(i in 1:length(fls_mc_dir_mn)) {
  tmp_dif <- read.table(fls_mc_dif_mn[i], header = T)[2, 1]
  tmp_dir <- read.table(fls_mc_dir_mn[i], header = T)[1, 1]
  h42_Lu_btt[i, "MC_Rrs_dif"] <- tmp_dif
  h42_Lu_btt[i, "MC_Rrs_dir"] <- tmp_dir
}

h42_Lu_btt

(h42_Lu_btt[, "HL_Rrs_dir"] - h42_Lu_btt[, "MC_Rrs_dir"]) / h42_Lu_btt[, "HL_Rrs_dir"]
(h42_Lu_btt[, "HL_Rrs_dif"] - h42_Lu_btt[, "MC_Rrs_dif"]) / h42_Lu_btt[, "HL_Rrs_dif"]



h42_Lu_inf <- data.frame(
  iop_c      = rep(1, 9),
  w0         = (10 - 1:9) / 10,  
  HL_Rrs_dif = numeric(9), 
  HL_Rrs_dir = numeric(9), 
  MC_Rrs_dif = numeric(9), 
  MC_Rrs_dir = numeric(9)
)
j = 1
for(i in seq(15, 32, 2)) {
  h42_Lu_inf$HL_Rrs_dif[j] <- readLines(fls_h42[i]) %>%
                              '['(24) %>%
                              strsplit(" ") %>%
                              unlist() %>%
                              as.numeric()%>%
                              '['(4)
 j = j+1
}
j = 1
for(i in seq(16, 32, 2)) {
  h42_Lu_inf$HL_Rrs_dir[j] <- readLines(fls_h42[i]) %>%
                              '['(24) %>%
                              strsplit(" ") %>%
                              unlist() %>%
                              as.numeric()%>%
                              '['(4)
 j = j+1
}

fls_mc_dif_mn <- list.files("run/hydrolight", "inf_out_dif_f_mn", full.names = T)
fls_mc_dir_mn <- list.files("run/hydrolight", "inf_out_dir_f_mn", full.names = T)
tmp_dif <- read.table(fls_mc_dif_mn[1], header = T)[2, ]
tmp_dir <- read.table(fls_mc_dir_mn[1], header = T)[1, ]
h42_Lu_inf[, "MC_Rrs_dif"] <- rev(as.numeric(tmp_dif))
h42_Lu_inf[, "MC_Rrs_dir"] <- rev(as.numeric(tmp_dir))

h42_Lu_inf

(h42_Lu_inf[, "HL_Rrs_dir"] - h42_Lu_inf[, "MC_Rrs_dir"]) / h42_Lu_inf[, "HL_Rrs_dir"]
(h42_Lu_inf[, "HL_Rrs_dif"] - h42_Lu_inf[, "MC_Rrs_dif"]) / h42_Lu_inf[, "HL_Rrs_dif"]


#
# Print out:
#

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

fnm <- "val_hydrolight_1"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
lim  <- range(h42_Lu_inf[, c("HL_Rrs_dir", "HL_Rrs_dif", "MC_Rrs_dif", "MC_Rrs_dir")])
xlab <- "HYDROLIGHT $R_{\\text{rs}}$ (sr$^{-1}$)"
ylab <- "This study $R_{\\text{rs}}$ (sr$^{-1}$)"
col  <- rainbow(9, start = 0, end = 0.8)
par(mar = c(5, 5, 3, 2))
plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, cex.lab = 2, 
  cex.axis = 2, log = "xy")
points(h42_Lu_inf[, "HL_Rrs_dif"], h42_Lu_inf[, "MC_Rrs_dif"], pch = 21, col = 1, 
  bg = col)
points(h42_Lu_inf[, "HL_Rrs_dir"], h42_Lu_inf[, "MC_Rrs_dir"], pch = 22, col = 1, 
  bg = col)
abline(a = 0, b = 1, lwd = 2)
legend(0.00035, 0.0005, tmp_lstat(rho::rstat(c(h42_Lu_inf[, "HL_Rrs_dif"], h42_Lu_inf[, "HL_Rrs_dir"]), 
  c(h42_Lu_inf[, "MC_Rrs_dif"], h42_Lu_inf[, "MC_Rrs_dir"]), units = "sr$^{-1}$"), digits = 4), bty = "n", cex = 1.5, y.intersp = 0.6)
vec <- seq(0.1, 0.9, 0.1)
legend(0.000072, 0.011, vec, title = "$\\omega_{0}$ (unitless):", bty = "n", 
  pch = 21, col = 1, pt.bg = rev(col), cex = 1.5, y.intersp = 0.6, pt.cex = 1)
legend(0.0003, 0.011, c("Diffuse", "Direct"), title = "Illumination:", bty = "n", 
  pch = 21:22, col = 1, pt.bg = 1, cex = 1.5, y.intersp = 0.6, pt.cex = 1)
legend("right", "\\textbf{A}", bty = "n", cex = 1.6)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

tmp_lstat(rho::rstat(c(h42_Lu_inf[, "HL_Rrs_dir"]), 
  c(h42_Lu_inf[, "MC_Rrs_dir"]), units = "sr$^{-1}$"), digits = 4)


fnm <- "val_hydrolight_2"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 7, height = 7, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
lim  <- range(h42_Lu_btt[, c("HL_Rrs_dir", "HL_Rrs_dif", "MC_Rrs_dif", "MC_Rrs_dir")])
xlab <- "HYDROLIGHT $R_{\\text{rs}}$ (sr$^{-1}$)"
ylab <- "This study $R_{\\text{rs}}$ (sr$^{-1}$)"
col  <- rainbow(7, start = 0, end = 0.8)
par(mar = c(5, 5, 3, 2))
plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, cex.lab = 2, 
  cex.axis = 2, log = "xy")
points(h42_Lu_btt[, "HL_Rrs_dif"], h42_Lu_btt[, "MC_Rrs_dif"], pch = 21, col = 1, 
  bg = col)
points(h42_Lu_btt[, "HL_Rrs_dir"], h42_Lu_btt[, "MC_Rrs_dir"], pch = 22, col = 1, 
  bg = col)
abline(a = 0, b = 1, lwd = 2)
legend(0.0036, 0.005, tmp_lstat(rho::rstat(c(h42_Lu_btt[, "HL_Rrs_dif"], h42_Lu_btt[, "HL_Rrs_dir"]), 
  c(h42_Lu_btt[, "MC_Rrs_dif"], h42_Lu_btt[, "MC_Rrs_dir"]), units = "sr$^{-1}$"), digits = 4), bty = "n", cex = 1.5, y.intersp = 0.6)
vec <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
legend(0.00135, 0.036, vec, title = "$c_{\\text{t}}$ (m$^{-1}$):", bty = "n", 
  pch = 21, col = 1, pt.bg = col, cex = 1.5, y.intersp = 0.6, pt.cex = 1)
legend(0.003, 0.036, c("Diffuse", "Direct"), title = "Illumination:", bty = "n", 
  pch = 21:22, col = 1, pt.bg = 1, cex = 1.5, y.intersp = 0.6, pt.cex = 1)
legend("right", "\\textbf{B}", bty = "n", cex = 1.6)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.6:
#
# Optically deep:
# iop_c  w0 HL_Rrs_dif HL_Rrs_dir   MC_Rrs_dif   MC_Rrs_dir
#     1 0.9 8.1159e-03 8.0848e-03 8.048902e-03 7.893042e-03
#     1 0.8 3.2858e-03 3.2910e-03 3.284856e-03 3.172578e-03
#     1 0.7 1.8294e-03 1.8464e-03 1.836338e-03 1.750806e-03
#     1 0.6 1.1462e-03 1.1640e-03 1.153310e-03 1.086132e-03
#     1 0.5 7.5220e-04 7.6746e-04 7.581793e-04 7.059887e-04
#     1 0.4 4.9626e-04 5.0811e-04 5.009382e-04 4.618779e-04
#     1 0.3 3.1669e-04 3.2516e-04 3.200971e-04 2.927246e-04
#     1 0.2 1.8375e-04 1.8906e-04 1.859340e-04 1.689169e-04
#     1 0.1 8.1317e-05 8.3822e-05 8.237726e-05 7.445078e-05
#
# Benchmark execution time: 21.70 min (1E7 rays)
#
# Finite bottom depth:
# iop_c bttd HL_Rrs_dif HL_Rrs_dir  MC_Rrs_dif  MC_Rrs_dir
#  0.15    2  0.0274010  0.0291880 0.027717750 0.031193110
#  0.30    2  0.0215350  0.0232290 0.021804490 0.024485360
#  0.60    2  0.0136540  0.0150640 0.013868030 0.015606150
#  1.50    2  0.0042614  0.0048300 0.004350255 0.004826642
#  3.00    2  0.0017667  0.0018526 0.001776245 0.001797149
#  6.00    2  0.0015534  0.0015719 0.001552703 0.001504515
# 15.00    2  0.0015522  0.0015702 0.001560659 0.001489743
#
# Benchmark execution time: 1.35, 1.98, 2.73, 3.81, 4.90, 5.71, 6.06 min (1E7 rays)


# Version 1.5:
#
# Optically deep:
# iop_c  w0 HL_Rrs_dif HL_Rrs_dir   MC_Rrs_dif   MC_Rrs_dir
#     1 0.9 8.1159e-03 8.0848e-03 8.048037e-03 7.894907e-03
#     1 0.8 3.2858e-03 3.2910e-03 3.280672e-03 3.176226e-03
#     1 0.7 1.8294e-03 1.8464e-03 1.833113e-03 1.756914e-03
#     1 0.6 1.1462e-03 1.1640e-03 1.150612e-03 1.092123e-03
#     1 0.5 7.5220e-04 7.6746e-04 7.557028e-04 7.106305e-04
#     1 0.4 4.9626e-04 5.0811e-04 4.986986e-04 4.649196e-04
#     1 0.3 3.1669e-04 3.2516e-04 3.182416e-04 2.943881e-04
#     1 0.2 1.8375e-04 1.8906e-04 1.846171e-04 1.696088e-04
#     1 0.1 8.1317e-05 8.3822e-05 8.170325e-05 7.460746e-05
#
# Benchmark execution time: 2.70 min (1E6 rays)
#
# Finite bottom depth:
# iop_c bttd HL_Rrs_dif HL_Rrs_dir  MC_Rrs_dif  MC_Rrs_dir
#  0.15    2  0.0274010  0.0291880 0.027694820 0.031169280
#  0.30    2  0.0215350  0.0232290 0.021815300 0.024500120
#  0.60    2  0.0136540  0.0150640 0.013867760 0.015598430
#  1.50    2  0.0042614  0.0048300 0.004349216 0.004863121
#  3.00    2  0.0017667  0.0018526 0.001777061 0.001794850
#  6.00    2  0.0015534  0.0015719 0.001557771 0.001495310
# 15.00    2  0.0015522  0.0015702 0.001554816 0.001496086
#
# Benchmark execution time: 0.13, 0.23, 0.35, 0.50, 0.65, 0.68, 0.76 min (1E6 rays)

