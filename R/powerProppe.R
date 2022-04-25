# Simulate power study in PRO2022 setup

nMC = 1e4
p = 0.95 # Target probability
mSeq = seq(80, 200, by=5)

s = 0.698 # Ensures approx. 0.995 PICP

ms = pw = c()
for(k in seq_along(mSeq)) {
  M = mSeq[k]
  OK   = 0
  Stab = rep(0, nMC)
  for (i in 1:nMC) {
    Up      = qnorm((1+p)/2)
    Etest   = rnorm(M, sd = s)
    S       = sum(abs(Etest) <= Up)
    Stab[i] = S / M
    OK      = OK + testBPCI(S, M, p, method = "wilsoncc")
  }
  ms[k] = mean(Stab)
  pw[k] = 1 - OK / nMC
}

plot(mSeq, pw,
     type = 'b',
     pch = 19,
     col = cols[5],
     xlab = 'M',
     ylab = 'Power',
     main = expression(nu[0.95] == 0.995)
)
abline(h=0.8,lty=2, col = cols[2])
