
fl <- "plots/isotropic.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/isotropic_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "isotropic"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

fl <- "plots/rayleigh_air.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/rayleigh_air_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "rayleigh_air"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


fl <- "plots/rayleigh_water.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/rayleigh_water_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "rayleigh_water"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", cex.main = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
box(lwd = 2)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

fl <- "plots/hg.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/hg_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "hg"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", 
  log = "y", cex.main = 2)
box(lwd = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0.1, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

fl <- "plots/ff_full.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/ff_full_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "ff_full"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", 
  log = "y", cex.main = 2)
box(lwd = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0.1, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))

fl <- "plots/ff_trunc.txt"
data1 <- read.table(fl, header = T)
fl <- "plots/ff_trunc_smp.txt"
data2 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "ff_trunc"
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", 
  log = "y", cex.main = 2)
box(lwd = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
box(lwd = 2)
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "Density"
brks = seq(0, 180, 5)
hist(data2[, "psi"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", ylim = c(range(data1[, "PDF"], na.rm = T)), 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0.1, 180, 50) * pi / 180, labels = seq(0, 180, 50), cex.axis = 2)
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))






library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "all"
cols = rev(rainbow(n = 7, start = 0, end = 0.8))
tikz(paste0("plots/",fnm,".tex"), width = 9.8, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,2), mar = c(5, 5, 3, 2))
xlab = "Scattering angle ($\\psi$, $^{\\circ}$)"
ylab = "$\\tilde{\\beta}(\\psi)$ (sr$^{-1}$)"
fl <- "plots/isotropic.txt"
data1 <- read.table(fl, header = T)
plot(data1[, "Angle"], data1[, "PF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = cols[1], xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), lwd = 2, main = "\\textbf{Angular dependency}", 
  log = "y", ylim = c(1E-3,1E3), cex.main = 2)
fl <- "plots/rayleigh_air.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PF"], lwd = 2, col = cols[2])
fl <- "plots/rayleigh_water.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PF"], lwd = 2, col = cols[3])
fl <- "plots/hg.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PF"], lwd = 2, col = cols[4])
fl <- "plots/ff_full.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PF"], lwd = 2, col = cols[5])
fl <- "plots/ff_trunc.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PF"], lwd = 2, col = cols[6])
box(lwd = 2)
ylab = "$2 \\pi \\tilde{\\beta}(\\psi) \\sin \\psi$"
fl <- "plots/isotropic.txt"
data1 <- read.table(fl, header = T)
plot(data1[, "Angle"], data1[, "PDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = cols[1], xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{PDF}", cex.main = 2, lwd = 2, 
  ylim = c(0, 1.5))
fl <- "plots/rayleigh_air.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PDF"], lwd = 2, col = cols[2])
fl <- "plots/rayleigh_water.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PDF"], lwd = 2, col = cols[3])
fl <- "plots/hg.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PDF"], lwd = 2, col = cols[4])
fl <- "plots/ff_full.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PDF"], lwd = 2, col = cols[5])
fl <- "plots/ff_trunc.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "PDF"], lwd = 2, col = cols[6])
box(lwd = 2)
ylab = "$2 \\pi \\int\\limits_{0}^{\\psi} \\tilde{\\beta}(\\psi) \\sin \\psi \\, d\\psi$"
fl <- "plots/isotropic.txt"
data1 <- read.table(fl, header = T)
plot(data1[, "Angle"], data1[, "CDF"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = cols[1], xpd = T, cex.lab = 2, cex.axis = 2, 
  xlim = c(0, 180), main = "\\textbf{CDF}", cex.main = 2, lwd = 2)
fl <- "plots/rayleigh_air.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "CDF"], lwd = 2, col = cols[2])
fl <- "plots/rayleigh_water.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "CDF"], lwd = 2, col = cols[3])
fl <- "plots/hg.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "CDF"], lwd = 2, col = cols[4])
fl <- "plots/ff_full.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "CDF"], lwd = 2, col = cols[5])
fl <- "plots/ff_trunc.txt"
data1 <- read.table(fl, header = T)
lines(data1[, "Angle"], data1[, "CDF"], lwd = 2, col = cols[6])
box(lwd = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))




