# ------------------------------------------------------------------------------
# Validation against Mobley et al. (1993)
#
# The study present a set o canonical "problems". For all problems, the 
# refractive index of the water is 1.34. Only the Lu and E0u sensors are 
# validated here (Ed sensor not implemented). Only the sensor at 1 optical depth
# (1 m since beam attenuation = 1) is simulated here. The FOV for the radiance 
# sensor was not specified. A 20 degree FOV used here.
#
# Problem 2:
# Air-water flat, semi-infinite homogenous medium, black sky with SZA of 60 
# degres, Ed = 0.5. Petzold scattering, no inelastic scattering, single 
# scattering albedo is 0.2 and 0.9.
#
# Problem 6:
# Same as problem 2 but with Lambertian bottom at 5 meters and bi-hemispherical
# reflectance of 0.5.
#

#
# Generate input files:
#

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 60.0",
  "sim_saa(s): 0.0",
  "sim_f0: 0.5",
  "sky_resx: 1",
  "sky_resy: 1",
  "sky_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt", 
  "iop_na: 1.0",
  "iop_nw: 1.34",
  "iop_c: 1.0",
  "iop_nw0: 1",
  "iop_w0: 0.2",
  "scat_tp: rayleigh",
  "scat_trc: 1",
  "btt_d: INFINITY",
  "sns_tp: lu",
  "sns_fov: 10.0",
  "sns_pos: 0.0 0.0 1.0",
  "acc_geom: 0",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)

writeLines(input, "run/mobley_1993/mobley_1993_Lu_p1_od1_w02.txt")

input[grep("sim_nr", input)] <- "sim_nr: 1000000"
input[grep("iop_w0", input)] <- "iop_w0: 0.9"
writeLines(input, "run/mobley_1993/mobley_1993_Lu_p1_od1_w09.txt")

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 60.0",
  "sim_saa(s): 0.0",
  "sim_f0: 0.5",
  "sky_resx: 1",
  "sky_resy: 1",
  "sky_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt", 
  "iop_na: 1.0",
  "iop_nw: 1.34",
  "iop_c: 1.0",
  "iop_nw0: 1",
  "iop_w0: 0.2",
  "scat_tp: rayleigh",
  "scat_trc: 1",
  "btt_d: INFINITY",
  "sns_tp: e0u",
  "sns_pos: 0.0 0.0 1.0",
  "acc_geom: 0",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)

writeLines(input, "run/mobley_1993/mobley_1993_E0u_p1_od1_w02.txt")

input[grep("sim_nr", input)] <- "sim_nr: 1000000"
input[grep("iop_w0", input)] <- "iop_w0: 0.9"
writeLines(input, "run/mobley_1993/mobley_1993_E0u_p1_od1_w09.txt")

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 60.0",
  "sim_saa(s): 0.0",
  "sim_f0: 0.5",
  "sky_resx: 1",
  "sky_resy: 1",
  "sky_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt", 
  "iop_na: 1.0",
  "iop_nw: 1.34",
  "iop_c: 1.0",
  "iop_nw0: 1",
  "iop_w0: 0.2",
  "scat_tp: petzold",
  "scat_trc: 1",
  "btt_d: INFINITY",
  "sns_tp: lu",
  "sns_fov: 10.0",
  "sns_pos: 0.0 0.0 1.0",
  "acc_geom: 0",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)

writeLines(input, "run/mobley_1993/mobley_1993_Lu_p2_od1_w02.txt")

input[grep("sim_nr", input)] <- "sim_nr: 1000000"
input[grep("iop_w0", input)] <- "iop_w0: 0.9"
writeLines(input, "run/mobley_1993/mobley_1993_Lu_p2_od1_w09.txt")

input[grep("bttd", input)]   <- "bttd: 5.0"
input[grep("bttdhr", input)] <- "bttdhr: 0.5"
writeLines(input, "run/mobley_1993/mobley_1993_Lu_p6_od1_w02.txt")

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 60.0",
  "sim_saa(s): 0.0",
  "sim_f0: 0.5",
  "sky_resx: 1",
  "sky_resy: 1",
  "sky_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt", 
  "iop_na: 1.0",
  "iop_nw: 1.34",
  "iop_c: 1.0",
  "iop_nw0: 1",
  "iop_w0: 0.2",
  "scat_tp: petzold",
  "scat_trc: 1",
  "btt_d: INFINITY",
  "sns_tp: e0u",
  "sns_pos: 0.0 0.0 1.0",
  "acc_geom: 0",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)
input[grep("sns_tp", input)] <- "sns_tp: e0u"
writeLines(input, "run/mobley_1993/mobley_1993_E0u_p2_od1_w09.txt")

input[grep("sim_nr", input)] <- "sim_nr: 1000000"
input[grep("iop_w0", input)] <- "iop_w0: 0.2"
writeLines(input, "run/mobley_1993/mobley_1993_E0u_p2_od1_w02.txt")

input[grep("bttd", input)]   <- "bttd: 5.0"
input[grep("bttdhr", input)] <- "bttdhr: 0.5"
writeLines(input, "run/mobley_1993/mobley_1993_E0u_p6_od1_w02.txt")

#
# Run MC code:
#

dir     <- "run/mobley_1993"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

rd  <- getwd()
ex  <- "/home/alexandre/Documents/active_research/owa_shading"

setwd(ex)
for(i in 1:length(fls)) {
  print(fls[i])
  cmd <- paste("./mc_solver_free_v1.5.o", file.path("src_v1.5/validation/code", dir, fls[i]))
  system(cmd)
}

setwd(rd)

#
# Load data and compare:
#

ref_p1_od1_w09 <- c(Lu = 0.04850, E0u = 0.372)
ref_p1_od1_w02 <- c(Lu = 0.0017200, E0u = 0.0134)
ref_p2_od1_w09 <- c(Lu = 0.0069900, E0u = 0.0931)
ref_p2_od1_w02 <- c(Lu = 0.0000547, E0u = 0.000966)
ref_p6_od1_w02 <- c(Lu = 0.0000684, E0u = 0.000981)
ref <- rbind(ref_p1_od1_w09, ref_p1_od1_w02, ref_p2_od1_w09, ref_p2_od1_w02,
 ref_p6_od1_w02)

fls_mc_dir_mn <- list.files(dir, "out_dir_f_mn", full.names = T)

res_p1_od1_w02 <- c(
  Lu = read.table(fls_mc_dir_mn[6], header = T)[1, 1],
  E0u = read.table(fls_mc_dir_mn[1], header = T)[1, 1]
)

res_p1_od1_w09 <- c(
  Lu = read.table(fls_mc_dir_mn[7], header = T)[1, 1],
  E0u = read.table(fls_mc_dir_mn[2], header = T)[1, 1]
)

res_p2_od1_w02 <- c(
  Lu = read.table(fls_mc_dir_mn[8], header = T)[1, 1],
  E0u = read.table(fls_mc_dir_mn[3], header = T)[1, 1]
)

res_p2_od1_w09 <- c(
  Lu = read.table(fls_mc_dir_mn[9], header = T)[1, 1],
  E0u = read.table(fls_mc_dir_mn[4], header = T)[1, 1]
)

res_p6_od1_w02 <- c(
  Lu = read.table(fls_mc_dir_mn[10], header = T)[1, 1],
  E0u = read.table(fls_mc_dir_mn[5], header = T)[1, 1]
)
res <- rbind(res_p1_od1_w09, res_p1_od1_w02, res_p2_od1_w09, 
  res_p2_od1_w02, res_p6_od1_w02)

rbind(ref, res)

# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.5:
#                          Lu          E0u
# ref_p1_od1_w09 4.850000e-02 0.3720000000
# ref_p1_od1_w02 1.720000e-03 0.0134000000
# ref_p2_od1_w09 6.990000e-03 0.0931000000
# ref_p2_od1_w02 5.470000e-05 0.0009660000
# ref_p6_od1_w02 6.840000e-05 0.0009810000
#
# res_p1_od1_w09 4.887109e-02 0.3062101000
# res_p1_od1_w02 1.706958e-03 0.0107262100
# res_p2_od1_w09 6.534022e-03 0.0003363508
# res_p2_od1_w02 5.611341e-05 0.0003574784
# res_p6_od1_w02 6.791724e-03 0.0003527031
#
# Benchmark execution time: 0.08, 1.08, 0.13, 0.13, 0.13, 0.11, 1.11, 0.15, 1.78, 1.71  min (1E6 rays)

ref_p1_od1_w09 0.0485000000 0.3720000000
ref_p1_od1_w02 0.0017200000 0.0134000000
ref_p2_od1_w09 0.0069900000 0.0931000000
ref_p2_od1_w02 0.0000547000 0.0009660000
ref_p6_od1_w02 0.0000684000 0.0009810000

res_p1_od1_w09 0.0488266500 0.3064041000
res_p1_od1_w02 0.0017063530 0.0107182300
res_p2_od1_w09 0.0095954430 0.0003353583
res_p2_od1_w02 0.0000556493 0.0003652727
res_p6_od1_w02 0.0070082490 0.0003614099



