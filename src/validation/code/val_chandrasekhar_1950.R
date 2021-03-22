# ------------------------------------------------------------------------------
# Validation against Chandrasekhar (1950)
#
# The results are presented in Gordon (1985) for a simple system: semi-infinite
# medium with unit refractive index, isotropic scattering, single scattering 
# albedo of 0.8 and unit beam attenuation. Direct ilumination from 60 degrees.
#

#
# Generate input files:
#

input <- c(
  "sim_nr: 10000000",
  "sim_ns: 1",
  "sim_sza(s): 60.0",
  "sim_saa(s): 0.0",
  "sim_f0: 1.0",
  "skr_resx: 1",
  "skr_resy: 1",
  "skr_fls:",
  "anc/skyrad/550_continental_clean_80_02_sza_60_sky.txt",
  "iop_na: 1.0",
  "iop_nw: 1.0",
  "iop_c: 1.0",
  "iop_nw0: 1",
  "iop_w0: 0.8",
  "scat_tp: isotropic",
  "btt_d: INFINITY",
  "src_tp: rad",
  "src_fov: 0.001",
  "src_axs: 0.0 0.0",
  "src_ref_o: 0.0 0.0 -0.0",
  "src_rel_o: 0.0 0.0 0.0",
  "acc_geom: none",
  "acc_ext: 1.0",
  "acc_resx: 0.1",
  "acc_resy: 0.1",
  "str_def: 0"
)
writeLines(input, "run/chandrasekhar_1950/chandrasekhar_1950_Lu.txt")

input[grep("src_tp", input)] <- "src_tp: pir"
input <- input[-grep("src_fov", input)]
writeLines(input, "run/chandrasekhar_1950/chandrasekhar_1950_Eu.txt")

#
# Run MC code:
#

dir     <- "run/chandrasekhar_1950"
fls     <- list.files(dir, ".txt")
fls.out <- list.files(dir, "out")
fls.tmp <- list.files(dir, "temp")
fls     <- setdiff(setdiff(fls, fls.out), fls.tmp)

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

ref <- c(Lu = 0.09583, Eu = 0.36799)

fls_mc_dir_mn <- list.files(dir, "out_dir_f_mn", full.names = T)
res <- c(
  Lu = read.table(fls_mc_dir_mn[2], header = T)[1, 1],
  Eu = read.table(fls_mc_dir_mn[1], header = T)[1, 1]
)

rbind(ref, res)

# ------------------------------------------------------------------------------
# Previous validations
#

# Version 1.6:
#
#             Lu        Eu
# ref 0.09583000 0.3679900
# res 0.09585718 0.3679013
#
# Benchmark execution time: 1.80 to 2.13 min (1E7 rays)

# Version 1.5:
#
#             Lu       Eu
# ref 0.09583000 0.367990
# res 0.09585895 0.367974
#
# Benchmark execution time: 0.15 to 0.23 min (1E6 rays)

