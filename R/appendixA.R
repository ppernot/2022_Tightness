figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2
set.seed(123)

## Unit variance distributions ####
Normal = function(N)
  rnorm(N)
T3 = function(N,df=3)
  rt(N, df = 3) / sqrt(3)
Uniform = function(N)
  runif(N, -sqrt(3), sqrt(3))
Laplace = function(N, df = 1)
  normalp::rnormp(N, p = df) / sqrt(df^(2/df)*gamma(3/df)/gamma(1/df))
ftab = c('Uniform','Normal','Laplace','T3')
# Oracle = function(N) # Scaled-shifted Bernoulli
#   2*rbinom(N,size=1,prob=0.5) - 1
# Beta = function(N,p=1e-2)
#   2*rbeta(N,shape1=p,shape2=p) -1

nMC = 10^5
N = 5
uTrue = 1/sqrt(N)

resu = resm = rest = resz = list()
for (k in seq_along(ftab)) {
  fun = get(ftab[k])
  resu[[ftab[k]]] = resm[[ftab[k]]] = rest[[ftab[k]]] = resz[[ftab[k]]] = rep(0,nMC)
  for(j in 1:nMC) {
    S = fun(N) # Random sample
    mu  = mean(S)
    umu = sd(S)/sqrt(N)
    resm[[ftab[k]]][j] = mu
    resu[[ftab[k]]][j] = umu
    rest[[ftab[k]]][j] = mu/umu
    resz[[ftab[k]]][j] = mu/uTrue
  }
}

# Fig_A01 ####
png(file = paste0(figDir,'/Fig_A01.png'),
    width = 2*gPars$reso, height = 2*gPars$reso)
par(mfrow = c(2,2),
    mar = c(3,3,2,1),
    tcl = gPars$tcl,
    mgp = gPars$mgp,
    pty = gPars$pty,
    lwd = 2*gPars$lwd,
    cex = gPars$cex)

for (k in seq_along(ftab)) {
  D = density(rest[[ftab[k]]])
  D$y = D$y / max(D$y) * dt(0, df = N - 1)
  plot(
    D$x, D$y,
    type = 'l',
    main = ftab[k],
    xlim = c(-3,4),
    xlab = 'Score',
    xaxs = 'i',
    ylim = c(0, 0.5),
    yaxs = 'i',
    ylab = 'Density',
    col  = gPars$cols[2]
  )
  grid()
  curve(
    dt(x, df = N - 1),
    from = -4,
    to = 4,
    n = 1000,
    add = TRUE,
    lty = 2,
    col  = gPars$cols[2]
  )
  D = density(resz[[ftab[k]]])
  D$y = D$y / max(D$y) * dnorm(0)
  lines(D$x, D$y,
        col  = gPars$cols[5])
  curve(
    dnorm(x),
    from = -4,
    to = 4,
    n = 1000,
    add = TRUE,
    lty = 2,
    col  = gPars$cols[5]
  )
  legend(
    c(-3,0.45), bty = 'n', cex = 1, xjust = 0,
    title = paste0(
         'Var(T)=',signif(var(rest[[ftab[k]]]),2),'\n',
         'Var(Z)=',signif(var(resz[[ftab[k]]]),2)
       ),
    legend =c('t-score','z-score'),
    col = gPars$cols[c(2,5)],
    lwd = 2*gPars$lwd,
    pch = NA
  )
  box()
}
dev.off()

# Fig_A03 ####

png(file = paste0(figDir,'/Fig_A03a.png'),
    width = gPars$reso, height = gPars$reso)
sel = sample.int(nMC,size = 1000)
X = resu[['Normal']][sel]
Y = resm[['Normal']][sel]
ErrViewLib::plotEvsPU(
  X , Y ,
  runQuant = TRUE,
  # cumMAE = TRUE,
  scalePoints = scalePoints,
  label = 1,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A03b.png'),
    width = gPars$reso, height = gPars$reso)
uE = resu[['Normal']]
E  = resm[['Normal']]
ErrViewLib::plotConfidence(
  E, uE,
  legend = 'Noisy data',
  oracle = FALSE,
  probref = TRUE,
  conf_probref = TRUE,
  label = 2, ylim = c(0,1.1),
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A03c.png'),
    width = gPars$reso, height = gPars$reso)
uE = resu[['Normal']]
Z  = rest[['Normal']]
ErrViewLib::plotLZV(
  uE, Z,
  method = 'cho',
  xlab = 'uE',
  varZ = (N-1)/(N-3),
  label = 3,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A03d.png'),
    width = gPars$reso, height = gPars$reso)
uE = resu[['Normal']]
Z  = rest[['Normal']]
ErrViewLib::plotLZV(
  1:length(Z), Z,
  method = 'cho',
  xlab = 'Point index',
  nBin = 10,
  label = 3,
  gPars = gPars
)
abline(h=(N-1)/(N-3),lwd = gPars$lwd, col=2, lty=2)
dev.off()

png(file = paste0(figDir,'/Fig_A03e.png'),
    width = gPars$reso, height = gPars$reso)
uE = resu[['Normal']]
E  = resm[['Normal']]
ErrViewLib::plotRelDiag(
  uE, E,
  nBin = 10,
  nBoot = 1000,
  BSmethod = 'perc',
  label = 4,
  gPars = gPars
)
dev.off()


# Convergence of Var(T) ####
## Distributions
Normal = function(N)
  rnorm(N)
T3 = function(N, df = 3)
  rt(N, df = df)
Uniform = function(N)
  runif(N, -1, 1)
Exp1 = function(N, df = 1)
  normalp::rnormp(N, p = df)
Oracle = function(N)
  # Scaled-shifted Bernoulli
  2 * rbinom(N, size = 1, prob = 0.5) - 1
Beta = function(N, p = 0.5)
  2 * rbeta(N, shape1 = p, shape2 = p) - 1
Exp4 = function(N, df = 4)
  normalp::rnormp(N, p = df)

nMC = 10^5

ftab = c('Beta','Uniform','Exp4','Normal','Exp1','T3')

nSeq= c(5:14,seq(15,30,by=5))

resuVarT = resuMeanT = list()
for (k in seq_along(ftab)) {
  fun = get(ftab[k])
  resuVarT[[ftab[k]]]  = rep(0,length(nSeq))
  resuMeanT[[ftab[k]]] = rep(0,length(nSeq))
  for(i in seq_along(nSeq)) {
    N = nSeq[i]
    mu = umu = rep(0,nMC)
    for(j in 1:nMC) {
      S = fun(N) # Random sample
      mu[j]  = mean(S)
      umu[j] = sd(S)/sqrt(N)
    }
    sel = umu != 0
    t = mu[sel]/umu[sel]
    resuMeanT[[ftab[k]]][i] = mean(t)
    resuVarT[[ftab[k]]][i]  = var(t)
  }
}

for (k in seq_along(ftab)) {
  print(c(resuMeanT[[ftab[k]]][1],resuVarT[[ftab[k]]][1]))
}

# Fig_A02 ####
ftabp = c('Uniform','Exp4','Normal','Exp1','T3')
png(file = paste0(figDir,'/Fig_A02.png'),
    width = gPars$reso, height = gPars$reso)
par(mfrow = c(1,1),
    mar = c(3,3,2,1),
    tcl = gPars$tcl,
    mgp = gPars$mgp,
    pty = gPars$pty,
    lwd = 2*gPars$lwd,
    cex = gPars$cex)

for (k in seq_along(ftabp)) {
  if(k==1) {
    plot(
      nSeq,resuVarT[[ftabp[k]]],
      type = 'l', log = 'x',
      xlab = 'n',
      ylab = 'Var(T)',
      ylim = c(0.95,3),
      col=gPars$cols[k])
    grid(lwd=2)

  } else {
    lines(
      nSeq,resuVarT[[ftabp[k]]],
      col=gPars$cols[k])
  }
}
law = (nSeq-1)/(nSeq-3)
icol = which(ftabp=='Normal')
points(nSeq,law,
       pch=19, col = gPars$cols[icol])

abline(h=1, lty =2)
box()
legend(
  'topright', bty = 'n',
  legend = ftabp,
  col = gPars$cols,
  lty = 1,
  pch = NA
)
dev.off()


# Heteroscedastic case ####
set.seed(123)
nMC = 10^4
N = 5
resvH = reseH = resuH = resmH = ressH = restH = reszH = rep(0,nMC)
for(i in 1:nMC) {
  V   = runif(1,-2,2)
  uE  = 0.01*(1 + V^2)
  E   = rnorm(1, 0, uE)
  S   = rnorm(N, 0, uE)
  mu  = mean(S)
  umu = sd(S)
  reseH[i] = E
  resvH[i] = V
  resmH[i] = mu
  resuH[i] = umu / sqrt(N)
  ressH[i] = uE
  restH[i] = mu / (umu / sqrt(N))
  reszH[i] = V  / uE
}

# Fig_A04 ####
sel = sample.int(nMC,size = 1000)
png(file = paste0(figDir,'/Fig_A04a.png'),
    width = gPars$reso, height = gPars$reso)
X = resuH[sel]
Y = resmH[sel]
ErrViewLib::plotEvsPU(
  X , Y ,
  xlim = c(0,0.04),
  runQuant = TRUE,
  # cumMAE = TRUE,
  scalePoints = scalePoints,
  label = 1,
  # xlim = c(0,3),
  title = 'n = 5',
  gPars = gPars
)
dev.off()


png(file = paste0(figDir,'/Fig_A04b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  resmH, resuH,
  legend = 'Noisy data',
  oracle = FALSE,
  probref = TRUE,
  conf_probref = TRUE,
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A04c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  resuH, restH,
  method = 'cho',
  xlab = 'uE',
  nBin = 10,
  slide = FALSE,
  ylim = c(0,4),
  varZ =(N-1)/(N-3),
  label = 3,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A04d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  resvH, restH,
  method = 'cho',
  nBin = 10,
  slide = FALSE,
  ylim = c(0,4),
  varZ =(N-1)/(N-3),
  label = 4,
  gPars = gPars
)
dev.off()

set.seed(123)
nMC = 10^4
N = 10
resvH = reseH = resuH = resmH = ressH = restH = reszH = rep(0,nMC)
for(i in 1:nMC) {
  V   = runif(1,-2,2)
  uE  = 0.01*(1 + V^2)
  E   = rnorm(1, 0, uE)
  S   = rnorm(N, 0, uE)
  mu  = mean(S)
  umu = sd(S)
  reseH[i] = E
  resvH[i] = V
  resmH[i] = mu
  resuH[i] = umu / sqrt(N)
  ressH[i] = uE
  restH[i] = mu / (umu / sqrt(N))
  reszH[i] = V  / uE
}


# Fig_A05 ####
sel = sample.int(nMC,size = 1000)
png(file = paste0(figDir,'/Fig_A05a.png'),
    width = gPars$reso, height = gPars$reso)
X = resuH[sel]
Y = resmH[sel]
ErrViewLib::plotEvsPU(
  X , Y ,
  xlim = c(0,0.03),
  runQuant = TRUE,
  # cumMAE = TRUE,
  scalePoints = scalePoints,
  label = 1,
  # xlim = c(0,3),
  title = 'n = 10',
  gPars = gPars
)
dev.off()


png(file = paste0(figDir,'/Fig_A05b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotConfidence(
  resmH, resuH,
  legend = 'Noisy data',
  oracle = FALSE,
  probref = TRUE,
  conf_probref = TRUE,
  label = 2,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A05c.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  resuH, restH,
  method = 'cho',
  xlab = 'uE',
  nBin = 10,
  slide = FALSE,
  ylim = c(0,4),
  varZ =(N-1)/(N-3),
  label = 3,
  gPars = gPars
)
dev.off()

png(file = paste0(figDir,'/Fig_A05d.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  resvH, restH,
  method = 'cho',
  nBin = 10,
  slide = FALSE,
  ylim = c(0,4),
  varZ =(N-1)/(N-3),
  label = 4,
  gPars = gPars
)
dev.off()
