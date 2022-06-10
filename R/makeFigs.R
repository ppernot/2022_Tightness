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
set.seed(1234)
N   = 1000

# s2   = rchisq(N, df = 4)
# uE  = 0.01 * sqrt(s2 / mean(s2))
# E   = rnorm(N, 0, uE)
#
# SYNT01 = list(E = E, uE = uE)
# SYNT03 = list(E = E, uE = rep(0.01, N))
#
# E   = rnorm(N, 0, 0.01)
# SYNT02 = list(E = E, uE = uE)
# SYNT04 = list(E = E, uE = rep(0.01, N))


V = uE = E = rep(0,N)
for(i in 1:N) {
  V[i]   = runif(1,-2,2)
  uE[i]  = 0.01*(1 + V[i]^2)
  E[i]   = rnorm(1, 0, uE[i])
}
muE = sd(E)
SYNT01 = list(V = V, E = E, uE = uE)
SYNT03 = list(V = V, E = E, uE = rep(muE,N))

E   = rnorm(N, 0, muE)
SYNT02 = list(V = V, E = E, uE = uE)
SYNT04 = list(V = V, E = E, uE = rep(muE, N))

# Check R2 coeff for linear regression
x = SYNT01$uE
y = abs(SYNT01$E)
summary(lm(y~0+x))
summary(lm(y~x))
cor(x,y)^2


# Fig_01 ####
png(file = paste0(figDir,'/Fig_01a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT01$uE, SYNT01$E,
  label = 1,
  xlim = c(0,0.06),
  title = 'SYNT01',
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
  # logX = TRUE,
  label = 2,
  xlim = c(0,0.06),
  title = 'SYNT01',
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
  xlim = c(0,0.06),
  ylim = c(-0.15,0.15),
  label = 3,
  title = 'SYNT02',
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

# Fig_02 ####
# png(file = paste0(figDir,'/Fig_02a.png'),
#     width = gPars$reso, height = gPars$reso)
# ErrViewLib::plotEvsPU(
#   1:N, SYNT04$E/SYNT04$uE,
#   type = 'horiz',
#   xlab = 'Points index',
#   ylim = c(-5,5),
#   label = 1,
#   scalePoints = scalePoints,
#   gPars = gPars
# )
# dev.off()


png(file = paste0(figDir,'/Fig_02a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT01$V, SYNT01$E/SYNT01$uE,
  type = 'horiz',
  runQuant = TRUE,
  xlab = 'Predicted value, V',
  xlim = range(SYNT01$V),
  ylim = c(-5,5),
  label = 1,
  title = 'SYNT01',
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_02b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  SYNT03$V, SYNT03$E/SYNT03$uE,
  type = 'horiz',
  runQuant = TRUE,
  xlab = 'Predicted value, V',
  xlim = range(SYNT03$V),
  ylim = c(-5,5),
  label = 2,
  title = 'SYNT03',
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
  # xlim = range(SYNT03$V),
  ylim = c(-5,5),
  label = 3,
  title = 'SYNT03',
  scalePoints = scalePoints,
  gPars = gPars
)
dev.off()

# png(file = paste0(figDir,'/Fig_02d.png'),
#     width = gPars$reso, height = gPars$reso)
# ErrViewLib::plotEvsPU(
#   SYNT04$V, SYNT04$E/SYNT04$uE,
#   type = 'horiz',
#   runQuant = TRUE,
#   # xlab = 'Points index',
#   xlim = range(SYNT04$V),
#   ylim = c(-6,6),
#   label = 4,
#   title = 'SYNT04',
#   scalePoints = scalePoints,
#   gPars = gPars
# )
# dev.off()


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
         col = gPars$cols[ipr],
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
           col = gPars$cols[ipr])
  }
}
legend(
  'topleft', bty = 'n',
  legend = pSeq,
  title = 'p',
  cex = 0.8,
  lty = 1,
  pch = pchs,
  col = gPars$cols
)
box()
dev.off()

# Fig_04 ####
png(file = paste0(figDir,'/Fig_04a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  SYNT01$E, 1.96*SYNT01$uE,
  ordX = SYNT01$V,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = TRUE,
  ylim = c(0.7,1),
  title = 'SYNT01',
  label = 1,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_04b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLCP(
  SYNT03$E, 1.96*SYNT03$uE,
  ordX = SYNT03$V,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = TRUE,
  ylim = c(0.7,1),
  title = 'SYNT03',
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_04c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLRR(
  SYNT03$E, 1.96*SYNT03$uE,
  ordX = SYNT03$V,
  prob = 0.95,
  mycols = 2,
  nBin = 6,
  slide = FALSE,
  title = 'SYNT02',
  ylim = c(0,3),
  label = 3,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_04d.png'),
    width = gPars$reso, height = gPars$reso)
nMC = 1e4
p = c(0.25,0.5,0.75,0.95)
M = length(SYNT03$E)
U = matrix(NA,ncol=length(p),nrow=M)
pOK = pSrate = pSlo = pSup = c()
for(i in seq_along(p)) {
  # Define expanded prediction uncertainty
  U[,i] = fUp(SYNT03$E, p = p[i])
}

ErrViewLib::plotLCP(
  SYNT03$E, U,
  prob = p,
  ordX = SYNT03$V,
  logX = FALSE,
  mycols = c(7,5,3,2),
  slide = TRUE,
  nBin = 6,
  title = 'SYNT03',
  xlab = 'Predicted value, V',
  label = 4,
  gPars = gPars
)
dev.off()

# Fig_05 ####
# Nb: these might be slow. For faster rendering,
# add "method = 'cho'"

png(file = paste0(figDir,'/Fig_05a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  SYNT01$uE, SYNT01$E/SYNT01$uE,
  nBin = 6,
  method = 'cho',
  slide = TRUE,
  xlab = 'Prediction uncertainty, uE',
  col = 2,
  xlim = c(0.009,0.046),
  ylim = c(0,8),
  label = 1,
  gPars = gPars
)
ErrViewLib::plotLZV(
  SYNT02$uE, SYNT02$E/SYNT02$uE,
  nBin = 6,
  method = 'cho',
  slide = TRUE,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  'topright', bty='n',
  legend = c('SYNT01','SYNT02'),
  col = gPars$cols[c(2,5)],
  lwd = 2*gPars$lwd,
  lty = 1,
  pch = 16
)
dev.off()


png(file = paste0(figDir,'/Fig_05b.png'),
    width = gPars$reso, height = gPars$reso)
lims = c(0.002,0.03)
ErrViewLib::plotCalVar(
  SYNT01$uE, SYNT01$E,
  slide = FALSE,
  nBoot = 1500,
  logX = TRUE,
  nBin = 10,
  # xlim = lims,
  # ylim = lims,
  label = 2,
  col = 2,
  gPars = gPars
)
ErrViewLib::plotCalVar(
  SYNT02$uE, SYNT02$E,
  slide = FALSE,
  nBoot = 1500,
  nBin = 10,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  'topleft', bty='n',
  legend = c('SYNT01','SYNT02'),
  col = gPars$cols[c(2,5)],
  lwd = 2*gPars$lwd,
  lty = 1,
  pch = 16
)
dev.off()

png(file = paste0(figDir,'/Fig_05c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  SYNT01$E, SYNT01$uE,
  # legend = 'SYNT01',
  ylim = c(0,1.1),
  label = 3,
  gPars = gPars
)
ErrViewLib::plotConfidence(
  SYNT02$E, SYNT02$uE,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  'bottomleft', bty='n',
  legend = c('Oracle','SYNT01','SYNT02'),
  col = gPars$cols[c(1,2,5)],
  lwd = 2*gPars$lwd,
  lty = c(2,1,1),
  pch = NA
)
dev.off()

# Fig_06 ####
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

png(file = paste0(figDir,'/Fig_06a.png'),
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

png(file = paste0(figDir,'/Fig_06b.png'),
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
  label = 2,
  gPars = gPars
)
dev.off()


logX = TRUE
slide = FALSE
xlim = range(c(ua,ub))
xlab = 'U95'

png(file = paste0(figDir,'/Fig_06c.png'),
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
  legend = '   Model a',
  label = 3,
  legLoc = 'topright',
  legNcol = 1,
  gPars = gPars
)
ErrViewLib::plotLRR(
  E, ub,
  ordX = ub,
  logX = logX,
  mycols = 5,
  slide = slide,
  nBin = 8,
  add = TRUE,
  gPars = gPars
)
legend(
  0.47, 9.5, bty='n',
  cex = 0.8,
  legend = 'Model b',
  col = gPars$cols[5],
  lwd = gPars$lwd,
  lty = 1,
  pch = 16
)
dev.off()

png(file = paste0(figDir,'/Fig_06d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  E, ua,
  legend = 'Model a',
  type = 'p',
  label = 4,
  gPars = gPars
)
ErrViewLib::plotConfidence(
  E, ub,
  add = TRUE,
  col = 5,
  type = 'l',
  gPars = gPars
)
legend(
  6, 0.1, bty='n',
  legend = '  Model b',
  col = gPars$cols[5],
  lwd = 2*gPars$lwd,
  lty = 1,
  pch = NA
)
dev.off()

# Fig_07 ####

## BAK2021 ####
D = read.table('../Data/BAK2021_Data.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
systems = D[,1]
R  = D[,2]
C  = D[,3]
uC = D[,4] # U95
E  = R-C
uE = uC

png(file = paste0(figDir,'/Fig_07a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uE, E,
  runQuant = TRUE,
  logX  = TRUE,
  title = 'BAK2021',
  xlab = 'Prediction uncertainty, U95',
  scalePoints = scalePoints,
  label = 1,
  gPars = gPars)
dev.off()

png(file = paste0(figDir,'/Fig_07d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  E, uE,
  label = 4,
  gPars = gPars)
dev.off()

## PAN2015 ####
fileName = '../Data/PAN2015_Data.csv'
table = read.table(
  fileName,
  header = FALSE,
  stringsAsFactors = FALSE,
  encoding = 'utf8',
  colClasses = c(
    'character',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'character',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'numeric'
  )
)
# reshape table 12 columns -> 6 columns
tab = table[, 1:6]
colnames(tab) = paste0('V', 7:12)
tab = rbind(tab, table[, 7:12])
# Get rid of last line (duplicated to fill original table)
tab = tab[-nrow(tab), 1:4]
colnames(tab) = c('Species', 'Expt', 'BEEF', 'uBEEF')
ord = order(tab[, 3])
Calc  = tab[ord, 3]
uCalc = tab[ord, 4]
Ref   = tab[ord, 2]
E     = Ref - Calc
uE    = uCalc

png(file = paste0(figDir,'/Fig_07b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uE, E,
  runQuant = TRUE,
  logX = TRUE,
  title = 'PAN2015',
  label = 2,
  scalePoints = scalePoints,
  gPars = gPars)
dev.off()

png(file = paste0(figDir,'/Fig_07e.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  E, uE,
  ylim = c(0,1.2),
  label = 5,
  gPars = gPars)
dev.off()

## PAR2019 ####
file='../Data/PAR2019_Data.csv'
dat = read.csv(file)

method='mu'
Calc  = dat[[method]]
Ref   = dat[['expt.']]
uCalc = dat[['sigma']]
E   = Ref - Calc
uE  = uCalc

png(file = paste0(figDir,'/Fig_07c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotEvsPU(
  uE, E,
  runQuant = TRUE,
  logX = TRUE,
  title = 'PAR2019',
  label = 3,
  scalePoints = scalePoints,
  gPars = gPars)
dev.off()

png(file = paste0(figDir,'/Fig_07f.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  E, uE,
  ylim = c(0,1.2),
  label = 6,
  gPars = gPars)
dev.off()

# Fig_08 ####
D = read.csv('../Data/LIN2021_RBFE.csv', header = TRUE,
             check.names = FALSE, stringsAsFactors = FALSE)

N = 5 # Sampling size
systems = D[,1]
R = D[,2]
V = D[,3]
uE = uV = D[,4] / sqrt(N)
E = R-V
Z = E / uE

# Correct linear trend
reg = lm(E~V)
uEb = uE
Eb  = residuals(reg)
Zb  = Eb / uEb

# Impact of mean experimental uncertainty on Var(Z)
ErrViewLib::varZCI(Z,method = "cho")
uE_aug = sqrt(uE^2+0.4^2)
ErrViewLib::varZCI(E/uE_aug,method = "cho")

ErrViewLib::varZCI(Zb,method = "cho")
uEb_aug = sqrt(uEb^2+0.4^2)
ErrViewLib::varZCI(Eb/uEb_aug,method = "cho")



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
  xlim = c(min(uE),2),
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
ErrViewLib::plotConfidence(
  E, uE,
  legend = 'Orig. data',
  label = 3,
  gPars = gPars
)
ErrViewLib::plotConfidence(
  Eb, uEb,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  6, 0.1, bty='n',
  legend = '  Corr. data',
  col = gPars$cols[5],
  lwd = 2*gPars$lwd,
  lty = 1,
  pch = NA
)
dev.off()

label = 4
png(file = paste0(figDir,'/Fig_08d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  V, Eb, uy = uEb,
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
  uEb, Eb,
  runExt = runExt,
  runQuant = runQuant,
  logX = logX,
  xlab = xlab,
  xlim = c(min(uEb),2),
  ylim = c(-4,4),
  ylab = ylab,
  scalePoints = scalePoints,
  label = 5,
  gPars = gPars
)
dev.off()

# png(file = paste0(figDir,'/Fig_08f.png'),
#     width = gPars$reso, height = gPars$reso)
# ErrViewLib::plotCalVar(
#   uE, E,
#   slide = FALSE,
#   nBoot = 1500,
#   logX = TRUE,
#   nBin = 10,
#   label = 6,
#   col = 2,
#   # legend = 'Orig. data',
#   gPars = gPars
# )
# ErrViewLib::plotCalVar(
#   uEb, Eb,
#   slide = FALSE,
#   nBoot = 1500,
#   logX = TRUE,
#   nBin = 10,
#   add = TRUE,
#   col = 5,
#   # legend = 'Orig. data',
#   gPars = gPars
# )
# legend(
#   'bottomright', bty='n',
#   legend = c('Orig. data','Corr. data'),
#   col = gPars$cols[c(2,5)],
#   lwd = gPars$lwd,
#   lty = 3,
#   pch = 16
# )
# dev.off()

png(file = paste0(figDir,'/Fig_08f.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  V, Z,
  nBin = 10,
  method = 'cho',
  slide = TRUE,
  xlab = xlab1,
  ylim = c(0,400),
  col = 2,
  varZ = 2,
  label = 6,
  gPars = gPars
)
ErrViewLib::plotLZV(
  V, Zb,
  nBin = 10,
  method = 'cho',
  slide = TRUE,
  col = 5,
  add = TRUE,
  gPars = gPars
)
legend(
  'topleft', bty='n',
  legend = c('Orig. data','Corr. data'),
  col = gPars$cols[c(2,5)],
  lwd = gPars$lwd,
  lty = 3,
  pch = 16
)
dev.off()


# Fig_09 ####

D = read.table('../Data/ZHE2022_AIQM1.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
N = 8 # sampling size

Ea  = D[,2]
uEa = D[,3] #/ sqrt(N)
Za  = Ea/uEa

D = read.table('../Data/ZHE2022_ANI-1ccx.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)

Eb  = D[,2]
uEb = D[,3]  #/ sqrt(N)
Zb  = Eb/uEb


# Impact of experimental uncertainty
var(Za)
var(Zb)
# (N-1)/(N-3)
uExp = 0.1
var(Ea/(sqrt(uEa^2+uExp^2)))
var(Eb/(sqrt(uEb^2+uExp^2)))


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
ErrViewLib::plotConfidence(
  Ea, uEa,
  legend = 'AIQM1',
  label = 1,
  gPars = gPars
)
ErrViewLib::plotConfidence(
  Eb, uEb,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  6, 0.1, bty='n',
  legend = '  ANI-1ccx',
  col = gPars$cols[5],
  lwd = 2*gPars$lwd,
  lty = 1,
  pch = NA
)
dev.off()
png(file = paste0(figDir,'/Fig_09b.png'),
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
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_09c.png'),
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
  label = 3,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_09d.png'),
    width = gPars$reso, height = gPars$reso)

lims = c(0.04,60)
ErrViewLib::plotCalVar(
  uEa, Ea,
  slide = FALSE,
  nBoot = 1500,
  logX = TRUE,
  nBin = 10,
  xlim = lims,
  ylim = lims,
  col = 2,
  gPars = gPars
)
ErrViewLib::plotCalVar(
  uEb, Eb,
  slide = FALSE,
  nBoot = 1500,
  logX = TRUE,
  nBin = 10,
  add = TRUE,
  col = 5,
  gPars = gPars
)
legend(
  'bottomright', bty = 'n',
  legend = c('AIQM1','ANI-1ccx'),
  col = gPars$cols[c(2,5)],
  pch = 19,
  lty = 3
)
dev.off()

