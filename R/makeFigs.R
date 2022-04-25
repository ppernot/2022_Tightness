figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

# Functions ####
fUpNorm = function(X, indX = 1:length(X), p = 0.95) {
  # Half range of 100p% prediction interval for normal draws
  qt((1+p)/2, df = length(indX)-1) *
    sd(X[indX]) *
    sqrt(1+1/length(indX))
}
fUp = function(X, indX = 1:length(X), p = 0.95){
  0.5*diff(ErrViewLib::vhd(X[indX],p = c(0.5*(1-p),0.5*(1+p))))
}
testBPCI = function(S,M,p,method = "wilsoncc") {
  ci = DescTools::BinomCI(S, M, conf.level = 0.95, method =method)
  as.logical(p >= ci[,2] & p <= ci[,3])
}

# Datasets ####
set.seed(123)
N   = 1000
s2   = rchisq(N, df = 4)
uE  = 0.01 * sqrt(s2 / mean(s2))
E   = rnorm(N, 0, uE)
SYNT01 = list(E = E, uE = uE)
SYNT03 = list(E = E, uE = rep(0.01, N))

E   = rnorm(N, 0, 0.01)
SYNT02 = list(E = E, uE = uE)
SYNT04 = list(E = E, uE = rep(0.01, N))

V      = seq(1, 3, length.out = 1000)
Ecal   = rnorm(V, 0, 0.01 * V)
Etest  = rnorm(V, 0, 0.01 * V)
SYNT05 = list(V = V, Ecal = Ecal, Etest = Etest)


# Fig_01 ####
png(file = paste0(figDir,'/Fig_01a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT01$uE, SYNT01$E,
  label = 1,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_01b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT01$uE, SYNT01$E,
  runQuant = TRUE,
  cumMAE = TRUE,
  label = 2,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_01c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT02$uE, SYNT02$E,
  runExt = TRUE,
  cumMAE = TRUE,
  ylim = c(-0.05,0.05),
  label = 3,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

# Fig_02 ####
png(file = paste0(figDir,'/Fig_02a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  1:N, SYNT04$E/SYNT04$uE,
  type = 'horiz',
  xlab = 'Points index',
  ylim = c(-5,5),
  label = 1,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_02b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  1:N, SYNT04$E/SYNT04$uE,
  type = 'horiz',
  runQuant = TRUE,
  xlab = 'Points index',
  ylim = c(-5,5),
  label = 2,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_02c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  1:N, SYNT03$E/SYNT03$uE,
  type = 'horiz',
  runQuant = TRUE,
  xlab = 'Points index',
  ylim = c(-5,5),
  label = 3,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

# Fig_03 ####
pSeq = c(0.5, 0.75, 0.95)

png(file = paste0(figDir,'/Fig_03.png'),
    width = gPars$reso, height = gPars$reso)

par(
  mfrow = c(1, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  lwd = gPars$lwd
)

pchs = 15:(15+length(pSeq))

for (ipr in seq_along(pSeq)) {
  p = pSeq[ipr]
  load(file = paste0('../Data/PICPpowTest_',p,'_.RData'))

  # For each case estimate the size necessary for power=0.8
  xs = ys = c()
  for(i in 1:nrow(ms)) {
    x = mSeq
    y = pw[i,]
    z = ms[i,]
    ip = which(y-0.8 >=0)[1]
    if(is.na(ip) | ip == 1) {
      xs[i] = NA
      ys[i] = NA
    } else {
      xs[i] = x[ip-1] + (x[ip]-x[ip-1])*(y[ip-1]-0.8)/(y[ip-1]-y[ip])
      ys[i] = z[ip-1] + (z[ip]-z[ip-1])*(y[ip-1]-0.8)/(y[ip-1]-y[ip])
    }
  }
  if(ipr==1) {
    plot(ys-p,xs,
         type = 'b',
         log = 'y',
         pch = pchs[ipr],
         col = cols[ipr],
         xlim = c(-0.1,0.1),
         xlab = expression(nu[p]-p),
         ylim = c(100,2000),
         ylab = 'Sample size',
         yaxs = 'i'
    )
    grid(equilogs = FALSE)
  } else {
    points(ys-p,xs,
           type = 'b',
           pch = pchs[ipr],
           col = cols[ipr])
  }
}
legend(
  'topleft', bty = 'n',
  legend = pSeq,
  title = 'p',
  cex = 0.8,
  lty = 1,
  pch = pchs,
  col = cols
)
box()
dev.off()

# Fig_04 ####
png(file = paste0(figDir,'/Fig_04a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  SYNT01$E, 1.96*SYNT01$uE,
  ordX = SYNT01$uE,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = FALSE,
  xlab = 'Prediction uncertainty, uE',
  ylim = c(0.5,1),
  label = 1,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_04b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  SYNT02$E, 1.96*SYNT02$uE,
  ordX = SYNT02$uE,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = FALSE,
  xlab = 'Prediction uncertainty, uE',
  ylim = c(0.5,1),
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_04c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLRR(
  SYNT02$E, 1.96*SYNT02$uE,
  ordX = SYNT02$uE,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = FALSE,
  ylim = c(0.4,1.6),
  xlab = 'Prediction uncertainty, uE',
  label = 3,
  gPars = gPars
)
dev.off()

# Fig_05 ####

png(file = paste0(figDir,'/Fig_05a.png'),
    width = gPars$reso, height = gPars$reso)

uE = sd(SYNT05$Ecal)
Z  = SYNT05$Etest/uE
ErrViewLib::plotEvsPU(
  V, Z,
  type = 'horiz',
  runQuant = TRUE,
  xlab = 'Predicted value, V',
  label = 1,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_05b.png'),
    width = gPars$reso, height = gPars$reso)
nMC = 1e4
p = c(0.25,0.5,0.75,0.95)
M = length(SYNT05$Etest)
U = matrix(NA,ncol=length(p),nrow=M)
V = seq(1, 3, length.out = 1000)

pOK = pSrate = pSlo = pSup = c()
for(i in seq_along(p)) {
  pi = p[i]
  # Define expanded prediction uncertainty
  U[,i] = fUp(SYNT05$Ecal, p = pi)

}

ErrViewLib::plotLCP(
  SYNT05$Etest, U,
  prob = p,
  ordX = V,
  logX = FALSE,
  mycols = c(7,5,3,2),
  slide = TRUE,
  nBin = 6,
  xlab = 'Predicted value, V',
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_05c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLRR(
  SYNT05$Etest, U[,4],
  prob = 0.95,
  ordX = V,
  logX = FALSE,
  mycols = 2,
  slide = FALSE,
  nBin = 6,
  ylim = c(0,3),
  xlab = 'Predicted value, V',
  label = 3,
  # legLoc = 'top',
  gPars = gPars
)
dev.off()

# Fig_06 ####
# Nb: these might be slow. For faster rendering,
# add "method = 'cho'"

png(file = paste0(figDir,'/Fig_06a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  SYNT01$uE, SYNT01$E/SYNT01$uE,
  nBin = 6,
  method = 'cho',
  slide = TRUE,
  xlab = 'Prediction uncertainty, uE',
  # ylim = c(0.5,1),
  label = 1,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_06b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  SYNT02$uE, SYNT02$E/SYNT02$uE,
  nBin = 6,
  method = 'cho',
  slide = TRUE,
  xlab = 'Prediction uncertainty, uE',
  ylim = c(0,10),
  label = 2,
  gPars = gPars
)
dev.off()

# Fig_07 ####
D = read.table('../Data/PRO2022_Data.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)

R  = D[,1]
C  = D[,2]
ua = D[,3] # U95
ub = D[,4] # U95
E  = R-C

logX = TRUE
xlim = range(c(ua,ub))
xlab = 'U95'

runExt = FALSE
runQuant = TRUE
logX = TRUE
xlim = range(c(ua,ub))
xlab = 'U95'
ylab = 'Error'

png(file = paste0(figDir,'/Fig_07a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  ua, E,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  xlim = xlim,
  xlab = xlab,
  ylab = ylab,
  title = 'Model a',
  label = 1,
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_07d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  ub, E,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  xlim = xlim,
  xlab = xlab,
  ylab = ylab,
  title = 'Model b',
  scalePoints = scalePoints,
  label = 4,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_07b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  E, ua,
  ordX = ua,
  logX = logX,
  mycols = 2,
  slide = TRUE,
  nBin = 2,
  xlim = xlim,
  ylim = c(0.9,1),
  xlab = xlab,
  # title = 'Model a',
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_07e.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  E, ub,
  ordX = ub,
  logX = logX,
  mycols = 2,
  slide = TRUE,
  nBin = 2,
  xlim = xlim,
  ylim = c(0.9,1),
  xlab = xlab,
  # title = 'Model b',
  label = 5,
  gPars = gPars
)
dev.off()

logX = TRUE
slide = FALSE
xlim = range(c(ua,ub))
xlab = 'U95'

png(file = paste0(figDir,'/Fig_07c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLRR(
  E, ua,
  ordX = ua,
  logX = logX,
  mycols = 2,
  slide = slide,
  nBin = 8,
  xlim = xlim,
  ylim = c(0,10),
  xlab = xlab,
  # title = 'Model a',
  label = 3,
  legLoc = 'top',
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_07f.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLRR(
  E, ub,
  ordX = ub,
  logX = logX,
  mycols = 2,
  slide = slide,
  nBin = 8,
  xlim = xlim,
  ylim = c(0,10),
  xlab = xlab,
  # title = 'Model b',
  label = 6,
  legLoc = 'top',
  gPars = gPars
)
dev.off()

# Fig_08 ####

D = read.csv('../Data/LIN2021_RBFE.csv', header = TRUE,
             check.names = FALSE, stringsAsFactors = FALSE)
systems = D[,1]
R = D[,2]
V = D[,3]
uE = uV = D[,4]
E = R-V
N = length(R)
Z = E / uE

xlab1 = expression(Calculated~paste(Delta,Delta,G)~group("[",kcal/mol,"]"))

label = 1
png(file = paste0(figDir,'/Fig_08a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  V, E, uy = uE,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  xlab = xlab1,
  ylab = 'Error [kcal/mol]',
  plotBA = TRUE,
  outLiers = FALSE,
  ylim = c(-5,5),
  labels = systems,
  topMar = 3,
  gPars = gPars,
  scalePoints = scalePoints
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()


runExt = FALSE
runQuant = TRUE
logX = TRUE
xlab = 'uE [kcal/mol]'
ylab = 'Error [kcal/mol]'

png(file = paste0(figDir,'/Fig_08b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uE, E,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  # xlim = xlim,
  xlab = xlab,
  ylim = c(-4,4),
  ylab = ylab,
  # title = 'Model a',
  scalePoints = scalePoints,
  label = 2,
  gPars = gPars
)
dev.off()


png(file = paste0(figDir,'/Fig_08c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uE, Z,
  nBin = 3,
  logX = logX,
  method = 'cho',
  slide = TRUE,
  xlab = xlab,
  # ylim = c(0,10),
  label = 3,
  gPars = gPars
)
dev.off()

# Correct linear trend
reg = lm(E~V)
E  = residuals(reg)
Z  = E / uE

label = 4
png(file = paste0(figDir,'/Fig_08d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  V, E, uy = uE,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  xlab = xlab1,
  ylab = 'Error [kcal/mol]',
  plotBA = TRUE,
  outLiers = FALSE,
  ylim = c(-5,5),
  labels = systems,
  topMar = 3,
  gPars = gPars,
  scalePoints = scalePoints
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

png(file = paste0(figDir,'/Fig_08e.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uE, E,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  xlab = xlab,
  ylim = c(-4,4),
  ylab = ylab,
  scalePoints = scalePoints,
  label = 5,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_08f.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uE, Z,
  nBin = 3,
  logX = logX,
  method = 'cho',
  slide = TRUE,
  xlab = xlab,
  label = 6,
  gPars = gPars
)
dev.off()


# Add missing variance
vMiss = var(E) - mean(uE^2)
uEc = sqrt(uE^2 + vMiss)

label = 7
png(file = paste0(figDir,'/Fig_08g.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  V, E, uy = uEc,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  xlab = xlab1,
  ylab = 'Error [kcal/mol]',
  plotBA = TRUE,
  outLiers = FALSE,
  ylim = c(-5,5),
  labels = systems,
  topMar = 3,
  gPars = gPars,
  scalePoints = scalePoints
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

png(file = paste0(figDir,'/Fig_08h.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uEc, E,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  xlab = xlab,
  ylim = c(-4,4),
  ylab = ylab,
  scalePoints = scalePoints,
  label = 8,
  gPars = gPars
)
dev.off()

Zc = E / uEc
png(file = paste0(figDir,'/Fig_08i.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uEc, Zc,
  nBin = 3,
  logX = logX,
  method = 'cho',
  slide = TRUE,
  xlab = xlab,
  label = 9,
  gPars = gPars
)
dev.off()

# Fig_09 ####

D = read.table('../Data/ZHE2022_AIQM1.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)

Ea  = D[,2]
uEa = D[,3]
Za  = Ea/uEa

D = read.table('../Data/ZHE2022_ANI-1ccx.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)

Eb  = D[,2]
uEb = D[,3] #sqrt(D[,3]^2+2)
Zb  = Eb/uEb

logX = TRUE
xlim = range(c(uEa,uEb))
ylim = range(c(Ea,Eb))
xlab = 'uE [kcal/mol]'

runExt = FALSE
runQuant = TRUE
cumMAE = FALSE
logX = TRUE
ylab = 'Error [kcal/mol]'

png(file = paste0(figDir,'/Fig_09a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uEa, Ea,
  runExt = runExt,
  runQuant = runQuant,
  cumMAE = cumMAE,
  logX = logX,
  xlim = xlim,
  xlab = xlab,
  ylim = c(-10,10),
  ylab = ylab,
  scalePoints = scalePoints,
  title = 'AIQM1',
  label = 1,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_09b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uEb, Eb,
  runExt = runExt,
  runQuant = runQuant,
  cumMAE = cumMAE,
  logX = logX,
  xlim = xlim,
  xlab = xlab,
  ylim = c(-10,10),
  ylab = ylab,
  scalePoints = scalePoints,
  title = 'ANI-1ccx',
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_09c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uEb, Zb,
  nBin = 5,
  logX = logX,
  method = 'cho',
  slide = TRUE,
  xlab = xlab,
  # ylim = c(0,10),
  title = 'ANI-1ccx',
  label = 3,
  gPars = gPars
)
dev.off()
