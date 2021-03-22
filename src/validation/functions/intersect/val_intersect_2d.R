
dir <- "plots"
fls <- list.files(dir, "out_spr", full.names = F)
library(raster)

for(i in 1:length(fls)) {
  tmp1  <- as.matrix(read.table(file.path(dir, fls[i]), header = FALSE))
  tmpr1 <- raster(tmp1[-c(1, nrow(tmp1)), -c(1, ncol(tmp1))])
  extent(tmpr1) <- extent(-0.995, 0.995, -0.995, 0.995)
  tmpr1 = flip(tmpr1, "y")

  fn <- paste0(fls[i], ".png")
  fn <- gsub("out", "fig", fn)
  png(file.path(dir, fn), height = 1200, width = 1200, res = 230)
    par(mar = c(5, 5, 2, 2))
    image(tmpr1,  zlim = c(0, 1), xlab = "X position (m)", 
      ylab = "Y position (m)")
    abline(v = 0, lty = 2, col = "grey40")
    abline(h = 0, lty = 2, col = "grey40")
    rect(-0.25, -0.25, 0.25, 0.25, border = 2)
  dev.off()
}

