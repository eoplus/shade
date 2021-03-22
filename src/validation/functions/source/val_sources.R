
fl <- "plots/lu_10.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/lu_10_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_lu_10"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{S}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "EF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 5), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2)
ylab = "$\\tilde{S}(\\psi,\\phi) \\sin \\psi / \\int\\limits_{0}^{10^{\\circ}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 5), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$\\int\\limits_{0}^{\\psi} \\tilde{S}(\\psi) \\sin \\psi \\, d\\psi / \\int\\limits_{0}^{10^{\\circ}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 5), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 5, 0.25)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = (0:5) * pi / 180, labels = 0:5, cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


fl <- "plots/ed_eu.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/ed_eu_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_ed_eu"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{S}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "EF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2)
ylab = "$\\tilde{S}(\\psi,\\phi) \\sin \\psi / \\int\\limits_{0}^{\\frac{\\pi}{2}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$\\int\\limits_{0}^{\\psi} \\tilde{S}(\\psi) \\sin \\psi \\, d\\psi / \\int\\limits_{0}^{\\frac{\\pi}{2}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 3)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = c(0, 20, 40, 60, 80) * pi / 180, labels = c(0, 20, 40, 60, 80), 
  cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))



fl <- "plots/ed0_eu0.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/ed0_eu0_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_ed0_eu0"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{S}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "EF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = F, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2, 
  ylim = 1/2/pi * c(0.9, 1.1))
ylab = "$\\tilde{S}(\\psi,\\phi) \\sin \\psi / \\int\\limits_{0}^{\\frac{\\pi}{2}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$\\int\\limits_{0}^{\\psi} \\tilde{S}(\\psi) \\sin \\psi \\, d\\psi / \\int\\limits_{0}^{\\frac{\\pi}{2}} \\tilde{S}(\\psi,\\phi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 90), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Emission angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 3)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = c(0, 20, 40, 60, 80) * pi / 180, labels = c(0, 20, 40, 60, 80), 
  cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


