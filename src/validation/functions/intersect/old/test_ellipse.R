
dir <- "plots/ellipse"
fls <- list.files(dir, "out_spr", full.names = F)
library(raster)

for(i in 1:length(fls)) {
  tmp1  <- as.matrix(read.table(file.path(dir, fls[i]), header = FALSE))
  tmpr1 <- raster(tmp1[-c(1, nrow(tmp1)), -c(1, ncol(tmp1))])
  extent(tmpr1) <- extent(-0.995, 0.995, -0.995, 0.995)

  fn <- paste0(fls[i], ".png")
  fn <- gsub("out", "fig", fn)
  png(file.path(dir, fn), height = 1200, width = 1200, res = 230)
    par(mar = c(5, 5, 2, 2))
    image(tmpr1,  zlim = c(0, 1), xlab = "X position (m)", 
      ylab = "Y position (m)")
    abline(v = 0, lty = 2, col = "grey40")
    abline(h = 0, lty = 2, col = "grey40")
    a <- seq(0, 2*pi, length.out = 360)
    x <- 0.5 * cos(a)
    y <- 0.5 * sin(a)
    lines(x,y, col = 2)
  dev.off()
}

