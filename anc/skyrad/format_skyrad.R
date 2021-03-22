################################################################################
# format_skyrad.R
#
# Version 2021-01-11
#
# This R script is inteded to format the sky radiance distributions provided in 
# Castagna et al. (2019) for input to the Backward Monte Carlo code. It takes 
# the grid and calculates 1 x 1 degree cell averages, than normalizes it to 
# unity downwelling diffuse plane irradiance. Note that the values are per solid 
# angle and this is the appropriate quantity for weighting the Monte Carlo 
# simulations. 
#
# Castagna, A.; Johnson, C. B.; Voss, K.; Dierssen, H. M.; Patrick, H.; Germer, 
#   T. A.; Sabbe, K.; Vyverman, W. 2019. Uncertainty in global downwelling plane 
#   irradiance estimates from sintered polytetrafluoroethylene plaque radiance 
#   measurements. Applied Optics 58, 16, 4497. DOI: 10.1364/AO.58.004497
#
################################################################################

# Create uniform distribution:
#
# The uniform distribution is described by:
#
# integral integral L0 sin(theta) cos(theta) dthetha dphi = pi * L0,
# 
# where L0 is a constant. For an illumination of 1, L0 = 1/pi = 0.31831. 
#

tmp <- matrix((1 / pi), nrow = 91, ncol = 361)
tmp <- tmp[, 1:360] + tmp[, 2:361]
tmp <- tmp[1:90, ] + tmp[2:91, ]
tmp <- tmp / 4

theta <- (0.5:89.5) * pi / 180

scale <- 1 / sum(tmp * cos(theta) * sin(theta) * (pi / 180)^2) # ds*dy in radians
tmp   <- tmp * scale # remove small inacuracies from coarse resolution.

tmp <- formatC(tmp, digits = 6, format = "E")
write.table(tmp, file = "isotropic_sky.txt", col.names = FALSE, 
  row.names = FALSE, quote = FALSE)

# Create cardioidal distribution:
#
# The downwelling irradiance from a cardioidal distribution is described by:
#
# integral integral L0 (1 + 2 cos(theta)) sin(theta) cos(theta) dthetha dphi,
# 
# where (1 + 2 cos(theta)) describes the polar angle dependency. For an 
# irradiance of 1, is necessary to divide L0 by the integral, which gives
# L0 = 1 / 7.330383.

theta <- (0:90) * pi / 180
int   <- integrate(function(x) { 2 * pi * (1 + 2 * cos(x)) * cos(x) * sin(x)}, lower = 0, upper = pi/2)$val
tmp   <- matrix((1 / int) * (1 + 2 * cos(theta)), nrow = 91, ncol = 361)
tmp   <- tmp[, 1:360] + tmp[, 2:361]
tmp   <- tmp[1:90, ] + tmp[2:91, ]
tmp   <- tmp / 4

theta <- (0.5:89.5) * pi / 180
scale <- 1 / sum(tmp * cos(theta) * sin(theta) * (pi / 180)^2) # ds*dy in radians
tmp   <- tmp * scale # remove small inacuracies from coarse resolution.

tmp <- formatC(tmp, digits = 6, format = "E")
write.table(tmp, file = "cardioidal_sky.txt", col.names = FALSE, 
  row.names = FALSE, quote = FALSE)

# Convert data from Skyrad simulations:
#
# As before, the irradiance from a sky radiance distribution is:
#
# integral integral Lu sin(theta) cos(theta) dthetha dphi,
#
# but since Lu is not constant, it cannot be taken out of teh integral. The 
# solution is first to normalize the sky radiance by its diffuse downwelling 
# irradiance calculated with the equation above.
#
# The normalization for the Monte Carlo run will then be to divide each result
# by the average sky radiance:
#
# integral integral Lu sin(theta) cos(theta) dthetha dphi / 
#   integral integral sin(theta) cos(theta) dthetha dphi = 1/pi.
#
# So it does not matter the sky radiance distribution, 
#

fls   <- list.files("skyrad_original", pattern = "sky.txt", full.names = TRUE)
fls_i <- list.files("skyrad_original", pattern = "diffuse.txt", full.names = TRUE)
 
wave    = 550                 # any in seq(from = 350, to = 1000, by = 10)
aerosol = "maritime_clean" # any of "continental_clean", "continental_polluted", "maritime_clean", "maritime_polluted"
relumid = 95                  # any of 50, 80, 95
aot     = "02"                # any of "01", "02"

base    <- paste(aerosol, relumid, aot, sep = "_")
flaer   <- grep(base, fls, value = T)
flaer_i <- grep(base, fls_i, value = T)

theta <- (0.5:89.5) * pi / 180
for(i in 1:length(flaer)) {
  tmp <- read.table(flaer[i], header = T)
  tmp <- tmp[, paste0("l", wave)]
  dim(tmp) <- c(91, 361)
  tmp <- tmp
  tmp <- tmp[, 1:360] + tmp[, 2:361]
  tmp <- tmp[1:90, ] + tmp[2:91, ]
  tmp <- tmp / 4

  tmp_i <- read.table(flaer_i[i], header = T, row.names = 1)
  tmp   <- tmp / tmp_i[as.character(wave), ]

  scale <- 1 / sum(tmp * cos(theta) * sin(theta) * (pi / 180)^2) # ds*dy in radians
  tmp   <- tmp * scale # remove small inacuracies from coarse resolution.

  tmp <- formatC(tmp, digits = 6, format = "E")
  write.table(tmp, file = gsub("skyrad_original/", paste0(wave, "_"), flaer[i]), 
    col.names = FALSE, row.names = FALSE, quote = FALSE)
}


