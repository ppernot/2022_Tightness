# Estimate power for PICP testing
#
# # Functions ####
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

# Params ####
nMC  = 1e4 # Number of Monte Carlo runs
mSeq = c(seq(100,1000,by=50),seq(1000,3000,by=100)) # Sample sizes
sSeq = seq(0.5, 1.5, by = 0.025) # Sigma of candidate distribution
pSeq = c(0.5, 0.75, 0.95) # Target probabilities

for (p in pSeq) {
  ms = pw = matrix(0.0,ncol = length(mSeq), nrow = length(sSeq))
  Up = qnorm((1+p)/2) # Half-range of probability interval
  for(k in seq_along(mSeq)) {
    M = mSeq[k]     # Total Group size
    for (j in seq_along(sSeq)) {
      OK   = 0
      Stab = rep(0, nMC)
      for (i in 1:nMC) {
        Etest   = rnorm(M, sd = sSeq[j]) # Random sample
        S       = sum(abs(Etest) <= Up)  # Number od successes
        Stab[i] = S / M                  # PICP
        OK      = OK + testBPCI(S, M, p, method = "wilsoncc")
      }
      ms[j, k] = mean(Stab)   # Mean PICP
      pw[j, k] = 1 - OK / nMC # Power
    }

  }
  save(mSeq,sSeq,p,ms,pw,
       file = paste0('../Data/PICPpowTest_',p,'_.RData'))
}
