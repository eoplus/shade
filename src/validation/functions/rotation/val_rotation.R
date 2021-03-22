
rad <- function (deg) { deg * pi / 180 }

rot_mat_ZYX <- function(alpha, beta, gama) {
 vec <- c(
   cos(alpha), -sin(alpha), 0,
   sin(alpha), cos(alpha), 0,
   0, 0, 1
 )
 Rz = matrix (vec, 3, 3, byrow = T)
 vec <- c(
   cos(beta), 0, sin(beta),
   0, 1, 0,
   -sin(beta), 0, cos(beta)
 )
 Ry = matrix (vec, 3, 3, byrow = T)
 vec <- c(
   1, 0, 0,
   0, cos(gama), -sin(gama),
   0, sin(gama), cos(gama)
 )
 Rx = matrix (vec, 3, 3, byrow = T)
 Rz %*% Ry %*% Rx
}

rot_ZYX_frame = function(alpha, beta, gama, main) {
  require(plot3D)
  scatter3D(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
    phi = 30, bty = "g", type = "p", ticktype = "detailed", lwd = 2, 
    main = main)
  arrows3D(x0 = -1, y0 = 0, z0 = 0, x1 = 1, y1 = 0, z1 = 0, col = 2, add = TRUE)
  arrows3D(x0 = 0, y0 = -1, z0 = 0, x1 = 0, y1 = 1, z1 = 0, col = 3, add = TRUE)
  arrows3D(x0 = 0, y0 = 0, z0 = -1, x1 = 0, y1 = 0, z1 = 1, col = 4, add = TRUE)
  R <- rot_mat_ZYX(alpha, beta, gama)
  ax <- R %*% matrix(c(-1,0,0,1,0,0), ncol = 2)
  ay <- R %*% matrix(c(0,-1,0,0,1,0), ncol = 2)
  az <- R %*% matrix(c(0,0,-1,0,0,1), ncol = 2)
  arrows3D(x0 = ax[1,1], y0 = ax[2,1], z0 = ax[3,1], x1 = ax[1,2], y1 = ax[2,2], 
    z1 = ax[3,2], col = 2, add = TRUE, lty = 2)
  arrows3D(x0 = ay[1,1], y0 = ay[2,1], z0 = ay[3,1], x1 = ay[1,2], y1 = ay[2,2], 
    z1 = ay[3,2], col = 3, add = TRUE, lty = 2)
  arrows3D(x0 = az[1,1], y0 = az[2,1], z0 = az[3,1], x1 = az[1,2], y1 = az[2,2], 
    z1 = az[3,2], col = 4, add = TRUE, lty = 2)
}

rot_mat_ZYZ <- function(phi, theta, alpha) {
 vec <- c(
   cos(phi), -sin(phi), 0,
   sin(phi),  cos(phi), 0,
   0, 0, 1
 )
 Rz = matrix (vec, 3, 3, byrow = T)
 vec <- c(
   cos(theta), 0, sin(theta),
   0, 1, 0,
   -sin(theta), 0, cos(theta)
 )
 Ry = matrix (vec, 3, 3, byrow = T)
 vec <- c(
   cos(alpha), -sin(alpha), 0,
   sin(alpha), cos(alpha), 0,
   0, 0, 1
 )
 Rz2 = matrix (vec, 3, 3, byrow = T)
 Rz %*% Ry %*% Rz2
}

rot_ZYZ_frame = function(phi, theta, alpha, main) {
  require(plot3D)
  scatter3D(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
    phi = 30, bty = "g", type = "p", ticktype = "detailed", lwd = 2, 
    main = main)
  arrows3D(x0 = -1, y0 = 0, z0 = 0, x1 = 1, y1 = 0, z1 = 0, col = 2, add = TRUE)
  arrows3D(x0 = 0, y0 = -1, z0 = 0, x1 = 0, y1 = 1, z1 = 0, col = 3, add = TRUE)
  arrows3D(x0 = 0, y0 = 0, z0 = -1, x1 = 0, y1 = 0, z1 = 1, col = 4, add = TRUE)
  R <- rot_mat_ZYZ(phi, theta, alpha)
  ax <- R %*% matrix(c(-1,0,0,1,0,0), ncol = 2)
  ay <- R %*% matrix(c(0,-1,0,0,1,0), ncol = 2)
  az <- R %*% matrix(c(0,0,-1,0,0,1), ncol = 2)
  arrows3D(x0 = ax[1,1], y0 = ax[2,1], z0 = ax[3,1], x1 = ax[1,2], y1 = ax[2,2], 
    z1 = ax[3,2], col = 2, add = TRUE, lty = 2)
  arrows3D(x0 = ay[1,1], y0 = ay[2,1], z0 = ay[3,1], x1 = ay[1,2], y1 = ay[2,2], 
    z1 = ay[3,2], col = 3, add = TRUE, lty = 2)
  arrows3D(x0 = az[1,1], y0 = az[2,1], z0 = az[3,1], x1 = az[1,2], y1 = az[2,2], 
    z1 = az[3,2], col = 4, add = TRUE, lty = 2)
}

library(plot3D)

png("val_rotation_ZYX.png", height = 2000, width = 2000, res = 230)
  par(mfcol = c(2, 2), mar = c(2, 2, 2, 2))
  scatter3D(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
    phi = 30, bty = "g", type = "p", ticktype = "detailed", lwd = 2, 
    main = "Original")
  arrows3D(x0 = -1, y0 = 0, z0 = 0, x1 = 1, y1 = 0, z1 = 0, col = 2, add = TRUE)
  arrows3D(x0 = 0, y0 = -1, z0 = 0, x1 = 0, y1 = 1, z1 = 0, col = 3, add = TRUE)
  arrows3D(x0 = 0, y0 = 0, z0 = -1, x1 = 0, y1 = 0, z1 = 1, col = 4, add = TRUE)
  rot_ZYX_frame(rad(60), rad(45), rad(0), "(2) Rotation around Y = 45º")
  rot_ZYX_frame(rad(60), rad(0), rad(0), "(1) Rotation around Z = 60º")
  rot_ZYX_frame(rad(60), rad(45), rad(90), "(3) Rotation around X = 90º")
dev.off()

png("val_rotation_ZY0.png", height = 2000, width = 2000, res = 230)
  par(mfcol = c(2, 2), mar = c(2, 2, 2, 2))
  scatter3D(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
    phi = 30, bty = "g", type = "p", ticktype = "detailed", lwd = 2, 
    main = "Original")
  arrows3D(x0 = -1, y0 = 0, z0 = 0, x1 = 1, y1 = 0, z1 = 0, col = 2, add = TRUE)
  arrows3D(x0 = 0, y0 = -1, z0 = 0, x1 = 0, y1 = 1, z1 = 0, col = 3, add = TRUE)
  arrows3D(x0 = 0, y0 = 0, z0 = -1, x1 = 0, y1 = 0, z1 = 1, col = 4, add = TRUE)
  rot_ZYX_frame(rad(60), rad(45), rad(0), "(2) Rotation around Y = 45º")
  rot_ZYX_frame(rad(60), rad(0), rad(0), "(1) Rotation around Z = 60º")
dev.off()

png("val_rotation_ZYZ.png", height = 2000, width = 2000, res = 230)
  par(mfcol = c(2, 2), mar = c(2, 2, 2, 2))
  scatter3D(0, 0, 0, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
    phi = 30, bty = "g", type = "p", ticktype = "detailed", lwd = 2, 
    main = "Original")
  arrows3D(x0 = -1, y0 = 0, z0 = 0, x1 = 1, y1 = 0, z1 = 0, col = 2, add = TRUE)
  arrows3D(x0 = 0, y0 = -1, z0 = 0, x1 = 0, y1 = 1, z1 = 0, col = 3, add = TRUE)
  arrows3D(x0 = 0, y0 = 0, z0 = -1, x1 = 0, y1 = 0, z1 = 1, col = 4, add = TRUE)
  rot_ZYZ_frame(rad(60), rad(45), rad(0), "(2) Rotation around Y = 45º")
  rot_ZYZ_frame(rad(60), rad(0), rad(0), "(1) Rotation around Z = 60º")
  rot_ZYZ_frame(rad(60), rad(45), rad(90), "(3) Rotation around Z = 90º")
dev.off()

