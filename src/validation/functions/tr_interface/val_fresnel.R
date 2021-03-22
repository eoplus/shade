
fls <- list.files("res/val_fresnel", "txt", full.names = T)

aw <- read.table(fls[1], skip = 1, header = F)

rm <- tm <- matrix(0, ncol = 16, nrow = nrow(aw))

rm[, 1]  <- aw[, 3]
rm[, 2]  <- aw[, 4]
rm[, 5]  <- aw[, 4]
rm[, 6]  <- aw[, 3]
rm[, 11] <- aw[, 5]
rm[, 12] <- -aw[, 6]
rm[, 15] <- aw[, 6]
rm[, 16] <- aw[, 5]

tm[, 1]  <- aw[, 7]
tm[, 2]  <- aw[, 8]
tm[, 5]  <- aw[, 8]
tm[, 6]  <- aw[, 7]
tm[, 11] <- aw[, 9]
tm[, 12] <- -aw[, 10]
tm[, 15] <- aw[, 10]
tm[, 16] <- aw[, 9]

na  <- complex(real = 1.000278, imaginary = -0.000E-00)
nw  <- complex(real = 1.342033, imaginary = -2.442E-09)

tb <- Re(atan(nw/na)) * 180 / pi

library(tikzDevice)

fnm <- "val_fresnel_air_water"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 14, height = 14, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mar = c(5, 5, 1, 1), mfcol = c(4, 4))
xlab = "$\\theta_{\\text{i}}$ ($^{\\circ}$)"
ylab = "T or R"
for(i in 1:16) {
  plot(NA, xlim = c(0, 90), ylim = c(-0.5, 2), xaxs = "i", yaxs = "i", 
    xlab = xlab, ylab = ylab, cex.lab = 2, cex.axis = 2)
  abline(h = 0)
  lines(aw[, 1], rm[, i], col = 2, lwd = 2)
  lines(aw[, 1], tm[, i], col = 4, lwd = 2)
  abline(v = tb, lty = 2, lwd = 2)
}
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

rm <- rm / rm[, 1]
tm <- tm / tm[, 1]

fnm <- "val_fresnel_air_water_reduced"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 14, height = 14, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mar = c(5, 5, 1, 1), mfcol = c(4, 4))
xlab = "$\\theta_{\\text{i}}$ ($^{\\circ}$)"
ylab = "T or R"
for(i in 1:16) {
  plot(NA, xlim = c(0, 90), ylim = c(-1, 1), xaxs = "i", yaxs = "i", 
    xlab = xlab, ylab = ylab, cex.lab = 2, cex.axis = 2)
  abline(h = 0)
  lines(aw[, 1], rm[, i], col = 2, lwd = 2)
  lines(aw[, 1], tm[, i], col = 4, lwd = 2)
  abline(v = tb, lty = 2, lwd = 2)
}
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


aw <- read.table(fls[2], skip = 1, header = F)

rm <- tm <- matrix(0, ncol = 16, nrow = nrow(aw))

rm[, 1]  <- aw[, 3]
rm[, 2]  <- aw[, 4]
rm[, 5]  <- aw[, 4]
rm[, 6]  <- aw[, 3]
rm[, 11] <- aw[, 5]
rm[, 12] <- -aw[, 6]
rm[, 15] <- aw[, 6]
rm[, 16] <- aw[, 5]

tm[, 1]  <- aw[, 7]
tm[, 2]  <- aw[, 8]
tm[, 5]  <- aw[, 8]
tm[, 6]  <- aw[, 7]
tm[, 11] <- aw[, 9]
tm[, 12] <- -aw[, 10]
tm[, 15] <- aw[, 10]
tm[, 16] <- aw[, 9]

tc  <- Re(asin(na / nw)) * 180 / pi

fnm <- "val_fresnel_water_air"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 14, height = 14, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mar = c(5, 5, 1, 1), mfcol = c(4, 4))
xlab = "$\\theta_{\\text{i}}$ ($^{\\circ}$)"
ylab = "T or R"
for(i in 1:16) {
  plot(NA, xlim = c(0, 90), ylim = c(-0.6, 1), xaxs = "i", yaxs = "i", 
    xlab = xlab, ylab = ylab, cex.lab = 2, cex.axis = 2)
  abline(h = 0)
  lines(aw[, 1], rm[, i], col = 2, lwd = 2)
  lines(aw[, 1], tm[, i], col = 4, lwd = 2)
  abline(v = tc, lty = 2, lwd = 2)
}
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


rm <- rm / rm[, 1]
tm <- tm / tm[, 1]

fnm <- "val_fresnel_water_air_reduced"
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")
tikz(paste0("plots/",fnm,".tex"), width = 14, height = 14, standAlone = TRUE, pointsize = 18,
  packages = c("\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mar = c(5, 5, 1, 1), mfcol = c(4, 4))
xlab = "$\\theta_{\\text{i}}$ ($^{\\circ}$)"
ylab = "T or R"
for(i in 1:16) {
  plot(NA, xlim = c(0, 90), ylim = c(-1, 1), xaxs = "i", yaxs = "i", 
    xlab = xlab, ylab = ylab, cex.lab = 2, cex.axis = 2)
  abline(h = 0)
  lines(aw[, 1], rm[, i], col = 2, lwd = 2)
  lines(aw[, 1], tm[, i], col = 4, lwd = 2)
  abline(v = tc, lty = 2, lwd = 2)
}
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

