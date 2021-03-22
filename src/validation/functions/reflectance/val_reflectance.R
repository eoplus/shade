
fl <- "plots/lambert_cdf_i_00_000.txt"
data2 <- read.table(fl, header = T)

fl <- "plots/lambert_brdf_i_00_000_smp.txt"
data3 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_brdf_lamb_i_00_000"
tikz(paste0("plots/",fnm,".tex"), width = 14.7, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,3), mar = c(5, 5, 3, 2))
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 1.5)
hist(data3[, "Psi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 90, 20) * pi / 180, labels = seq(0, 90, 20), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 360, 6.0)
hist(data3[, "Phi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 360, 50) * pi / 180, labels = seq(0, 360, 50), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "PDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "PDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\theta} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "CDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\phi}\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta\\, d\\phi / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "CDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


fl <- "plots/lambert_cdf_i_60_000.txt"
data2 <- read.table(fl, header = T)

fl <- "plots/lambert_brdf_i_60_000_smp.txt"
data3 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_brdf_lamb_i_60_000"
tikz(paste0("plots/",fnm,".tex"), width = 14.7, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,3), mar = c(5, 5, 3, 2))
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 1.5)
hist(data3[, "Psi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 90, 20) * pi / 180, labels = seq(0, 90, 20), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 360, 6.0)
hist(data3[, "Phi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 360, 50) * pi / 180, labels = seq(0, 360, 50), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "PDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "PDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\theta} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "CDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\phi}\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta\\, d\\phi / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "CDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))



fl <- "plots/minnaert_cdf_i_00_000.txt"
data2 <- read.table(fl, header = T)

fl <- "plots/minnaert_brdf_i_00_000_smp.txt"
data3 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_brdf_minn_i_00_000"
tikz(paste0("plots/",fnm,".tex"), width = 14.7, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,3), mar = c(5, 5, 3, 2))
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 1.5)
hist(data3[, "Psi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 90, 20) * pi / 180, labels = seq(0, 90, 20), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 360, 6.0)
hist(data3[, "Phi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 360, 50) * pi / 180, labels = seq(0, 360, 50), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "PDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "PDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\theta} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "CDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\phi}\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta\\, d\\phi / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "CDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


fl <- "plots/minnaert_cdf_i_60_000.txt"
data2 <- read.table(fl, header = T)

fl <- "plots/minnaert_brdf_i_60_000_smp.txt"
data3 <- read.table(fl, header = T)

library(tikzDevice)

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
  "\\usepackage[T1]{fontenc}", 
  "\\usetikzlibrary{calc}", "\\usepackage{amssymb}"), tikzDocumentDeclaration = 
  "\\documentclass[18pt]{article}")

fnm <- "val_brdf_minn_i_60_000"
tikz(paste0("plots/",fnm,".tex"), width = 14.7, height = 9.8, 
  standAlone = TRUE, pointsize = 18, packages = c("\\usepackage{tikz}", 
  "\\usepackage[active,tightpage,psfixbb]{preview}",
  "\\PreviewEnvironment{pgfpicture}", "\\setlength\\PreviewBorder{0pt}",
  "\\usepackage{amssymb}", "\\usepackage{amsmath}"))
par(mfcol = c(2,3), mar = c(5, 5, 3, 2))
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 90, 1.5)
hist(data3[, "Psi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 90, 20) * pi / 180, labels = seq(0, 90, 20), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "Density"
brks = seq(0, 360, 6.0)
hist(data3[, "Phi_r"]*pi/180, xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
  xpd = T, cex.lab = 2, cex.axis = 2, freq = F, breaks = brks*pi/180, 
  xaxt = "n", 
  main = "\\textbf{Samples}", cex.main = 2, lwd = 2)
axis(1, at = seq(0, 360, 50) * pi / 180, labels = seq(0, 360, 50), cex.axis = 2)
box(lwd = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "PDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "PDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{PDF}", cex.main = 2)
xlab = "Reflected $\\theta$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\theta} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta / \\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta$"
plot(data2[, "Psi_r"], data2[, "CDF_Psi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
xlab = "Reflected $\\phi$ ($^{\\circ}$)"
ylab = "$\\int\\limits_{0}^{\\phi}\\int\\limits_{0}^{\\frac{\\pi}{2}} f(\\theta_{\\text{i}}; \\phi_{\\text{i}}; \\theta_{\\text{r}}; \\phi_{\\text{r}}) \\cos \\theta \\sin \\theta \\, d\\theta\\, d\\phi / \\rho_{\\text{dh}}(\\theta_{\\text{i}}; \\phi_{\\text{i}})$"
plot(data2[, "Phi_r"], data2[, "CDF_Phi_r"], type = "l", xaxs = "i", yaxs = "i", 
  xlab = xlab, ylab = ylab, col = 2, xpd = T, cex.lab = 2, cex.axis = 2, 
  lwd = 2, main = "\\textbf{CDF}", cex.main = 2)
dev.off()
tools::texi2pdf(paste0("plots/",fnm,".tex"))
file.rename(paste0(fnm,".pdf"), paste0("plots/",fnm,".pdf"))
file.remove(paste0(fnm,".aux"), paste0(fnm,".log"))


