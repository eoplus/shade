# ------------------------------------------------------------------------------
# Validation against Gordon 1985
#
# The study present an evaluation of the shadowing casted by a rectangular 
# surface (bi-dimensional). The tabulated values are for a demosntration with 
# the HG phase function and semi-infinite homogenous medium with unit refractive 
# index and a beam attenuation of 0.1 (1/m). Illumination is from the zenith 
# with a black sky.
#

#
# Generate input files:
#

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 0.0",
  "sim_saa(s): 0.0",
  "sim_f0: 1.0",
  "skr_resx: 1",
  "skr_resy: 1",
  "skr_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_00_sky.txt",
  "iop_na: 1.0",
  "iop_nw: 1.0",
  "iop_c: 0.1",
  "iop_nw0: 1",
  "iop_w0: 0.8",
  "scat_tp: hg",
  "scat_g: 0.0",
  "btt_d: INFINITY",
  "src_tp: rad",
  "src_fov: 0.001",
  "src_axs: 0.0 0.0",
  "src_ref_o: 0.0 0.0 0.0",
  "src_rel_o: 0.0 0.0 0.0",
  "acc_geom: none",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 1",
  "str_ncl: 0",
  "str_ncn: 0",
  "str_ncb: 1",
  "-4.5 0.0 -0.0000005 0.0 0.0 0.0 6.55 38.4 0.000001 0.0 0.0 0.0 0.0 1 1"
)
writeLines(input, "run/gordon_1985/gordon_1985_Lu_hg000.txt")

input[grep("scat_g: ", input)] <- "scat_g: 0.75"
writeLines(input, "run/gordon_1985/gordon_1985_Lu_hg750.txt")

input[grep("scat_g: ", input)] <- "scat_g: 0.875"
writeLines(input, "run/gordon_1985/gordon_1985_Lu_hg875.txt")

input[grep("scat_g: ", input)] <- "scat_g: 0.95"
writeLines(input, "run/gordon_1985/gordon_1985_Lu_hg950.txt")

#
# Run MC code:
#

dir     <- "run/gordon_1985"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  cmd <- paste("./mc_solver_shdw_v1.6.o", file.path("src/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

ref <- c(hg_000 = 0.1251, hg_750 = 0.1233, hg_875 = 0.1164, hg_950 = 0.1133)

fls_mc_dir_fmn <- list.files(dir, "out_dir_f_mn", full.names = T)
fls_mc_dir_smn <- list.files(dir, "out_dir_s_mn", full.names = T)

tmp1f <- read.table(fls_mc_dir_fmn[1], header = T)[1, 1]
tmp1s <- read.table(fls_mc_dir_smn[1], header = T)[1, 1]
tmp1  <- (tmp1f - tmp1s) / tmp1f

tmp2f <- read.table(fls_mc_dir_fmn[2], header = T)[1, 1]
tmp2s <- read.table(fls_mc_dir_smn[2], header = T)[1, 1]
tmp2  <- (tmp2f - tmp2s) / tmp2f

tmp3f <- read.table(fls_mc_dir_fmn[3], header = T)[1, 1]
tmp3s <- read.table(fls_mc_dir_smn[3], header = T)[1, 1]
tmp3  <- (tmp3f - tmp3s) / tmp3f

tmp4f <- read.table(fls_mc_dir_fmn[4], header = T)[1, 1]
tmp4s <- read.table(fls_mc_dir_smn[4], header = T)[1, 1]
tmp4  <- (tmp4f - tmp4s) / tmp4f

res <- c(hg_000 = tmp1, hg_750 = tmp2, hg_875 = tmp3, hg_950 = tmp4)

rbind(ref, res)

# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.6:
#        hg_000    hg_750    hg_875    hg_950
# ref 0.1251000 0.1233000 0.1164000 0.1133000
# res 0.1246287 0.1233775 0.1196655 0.1196591
#
# Benchmark execution time: 2.38, 5.08, 5.86, 6.85 min (1E7 rays)

# Version 1.5:
#        hg_000    hg_750    hg_875    hg_950
# ref 0.1251000 0.1233000 0.1164000 0.1133000
# res 0.1246092 0.1234743 0.1208892 0.1260079
#
# Benchmark execution time: 0.25, 0.53, 0.68, 0.80 min (1E6 rays)




