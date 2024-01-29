## ----echo=TRUE----------------------------------------------------------------
library(kitagawa)

## ----echo=TRUE----------------------------------------------------------------
library(psd)
data(Tohoku)
toh_orig <- with(subset(Tohoku, epoch=='seismic'), {
  cbind(
    scale(1e3*areal, scale=FALSE), # scale strain to nanostrain, remove mean
    scale(1e2*pressure.pore, scale=FALSE) # scale hPa to Pa, remove mean
  )
})
colnames(toh_orig) <- c('input','output')
toh.dat <- window(ts(toh_orig), 100, 2400)

## ----echo=FALSE, fig.show='hold', fig.width=7., fig.height=4.5----------------
library(RColorBrewer)
Set1 <- brewer.pal(8, 'Set1')
par(mar=c(3,3,0.2,0.2))
plot(toh.dat, yax.flip = TRUE, main="Strain and Pressure: 2011 M9 Tohoku")

## ----echo=FALSE, fig.show='hold', fig.width=7., fig.height=4.5----------------
windat <- scale(window(toh.dat, 1400, 1600))
plot(windat[,'input'], lty=5, type='l', ylab='', main='Rescaled Input and Output')
lines(-windat[,'output'], col=2, lwd=1.5)
lines(as.vector(time(windat)), scale(apply(windat,1, function(x){x <- abs(x); atan2(x[1],x[2])}))/2, lty=1, col=4)
legend('topleft', c('Input strain','Output pressure','Internal angle'), col=c(1,2,4), lty=c(2,1,1), lwd=c(1,1.5,1))

## ----echo=TRUE----------------------------------------------------------------
m <- lm(output ~ input - 1, as.data.frame(toh.dat))
strain_scaling <- coef(m)
signif(strain_scaling, 3) # GPa/strain

## ----echo=FALSE, fig.show='hold', fig.width=4.5, fig.height=4.5---------------
IO <- as.matrix(toh.dat)
plot(IO[,1], IO[,2], 
     asp=1, col=NA, 
     main='Pressure-strain correlation',
     xlab="Input (strain)", 
     ylab="Output (pore pressure)")
grid()
points(IO[,1], IO[,2], pch=3)
abline(m, col=2)

## -----------------------------------------------------------------------------
k <- 2*130 # number to start out with
gam <- seq(0.001, 1, by=0.001)
gamrat <- 2 * gam / (1 - gam)
Pgam <- pf(k*gamrat, 2, 4*k)

## ----echo=FALSE, fig.show='hold', fig.width=5.5, fig.height=4.5---------------
k2 <- 100
Pgam2 <- pf(k2*gamrat, 2, 4*k2)
k3 <- 10
Pgam3 <- pf(k3*gamrat, 2, 4*k3)
x.g <- ((1 - gam)*gamrat/2)
plot(x.g, Pgam, type='l', 
     main=expression(F(2*","~4*k)), 
     xlab=expression(gamma), 
     ylab=expression(p(gamma,k)), log='x')
lines(x.g, Pgam2, lty=5)
lines(x.g, Pgam3, lty=2)
legend('bottomright', parse(text=c(sprintf("k==%s",c(k,k2,k3)))), lty=c(1,5,2))

## ----echo=TRUE----------------------------------------------------------------
#!order_matters
class(toh.dat)
toh_to_pspec <- toh.dat[,c('input','output')]
toh.cs <- psd::pspectrum(toh_to_pspec, ntap.init=k, verbose=FALSE)

## ----echo=TRUE----------------------------------------------------------------
class(toh.cs)
str(toh.cs)

## ----echo=TRUE, fig.show='hold', fig.width=5., fig.height=7-------------------

f <- as.vector(toh.cs[['freq']]) # frequency
lf <- log10(f)
p <- 1/f # period

# coherence
Coh <- toh.cs[['coh']]
# wrapped phase (in radians)
Phi <- as.vector(toh.cs[['phase']])
suppressPackageStartupMessages(can.unwrap <- require(signal))
if (can.unwrap){
  # unwrap if possible
  Phi <- signal::unwrap(Phi)
}
# Admittance or Gain
G <- Mod(toh.cs[['transfer']])
G <- Coh * G
# Tapers
K <- toh.cs[['taper']] # 'tapers' object
k <- as.numeric(K)
# Uncertainty in the admittance
G.err <- sqrt((1 - Coh) / k)

## ----echo=TRUE----------------------------------------------------------------
csd <- data.frame(f, p, lf, Coh, k, G, G.err, Phi = Phi * 180 / pi)
csd.f <- subset(csd, p <= 100)

## ----echo=FALSE, fig.show='hold', fig.width=6., fig.height=6.5----------------

layout(matrix(1:3), heights=c(2,2,1.5))
par(oma=c(1,1,3,1), cex=0.8, las=1, tcl=-0.2, mgp=c(2,0.3,0))

par(mar=c(0.1, 4, 1, 4))
plot(Coh ~ f, csd.f, 
     log='x',
     xlab='', ylab='',
     type='l', xaxt='n', 
     ylim=c(-0.6,1), 
     xaxs='i', yaxt='n',
     yaxs='i', frame=FALSE)
abline(h=c(0,1), lty=3)
mtext('Coherence', font=2, line=-2.5, adj=0.3)
mtext('No. tapers', side=1, font=3, line=-2.7, adj=0.2)
axis(2, at=seq(0,1,by=0.2))
#axis(1, col=NA, col.ticks=1, labels=FALSE)
axis(3, col=NA, col.ticks=1)
title("Cross Spectrum: Pressure from Strain (Tohoku M9)", outer=TRUE)
par(new=TRUE)
plot(K, xaxs='i', yaxs='i', axes=FALSE, ylim=c(0,1000),xlab='', ylab='')
axis(4, at=seq(0,300,by=100))
par(new=FALSE)
box()

par(mar=c(0.1, 4, 0, 4))
nsig <- 2
with(csd.f, {
  lG <- log10(G)
  Delt <- nsig*log10(exp(1))*G.err/G
  Upper <- lG + Delt
  Lower <- lG - Delt
  ylim <- range(c(lG)) + log10(2)*c(-1,2)
  plot(f, lG, type='l', 
       log='x',
       yaxt='n', xlab='', ylab='',
       xaxt='n', 
       col=NA, frame=FALSE,
       xaxs='i', yaxs='i', ylim=ylim)
  polygon(c(f, rev(f)), c(Upper, rev(Lower)), col='lightcyan', border='lightgrey')
  lmsc <- log10(abs(strain_scaling))
  abline(h=lmsc, col=2, lty=2)
  mtext("scaling\nfrom lm", side=4, at=lmsc, col=2, line=0.2, font=3)
  lines(f, lG)
})
mtext('Admittance', font=2, line=-4.3, adj=0.3)
mtext(parse(text=sprintf("%s * sigma ~ 'uncert.'", nsig)), cex=1, adj=0.3, line=-5.5,  col='cyan4')
ll <- c(1,2,5)
lbls <- c(ll/1000, ll/100, ll/10, ll, ll*10)
ats <- log10(lbls)
lbls[c(F,T,T)] <- ""
axis(2, at=ats, labels=lbls)
box()

par(mar=c(1, 4, 0, 4))
degadd <- 180
plot(Phi + degadd ~ f, csd.f, 
     log='x',
     type='l', #col='lightgrey', 
     xlab='Frequency, Hz', ylab="Degrees (<0 = lag)",
     xaxs='i', frame=FALSE,
     #yaxs='i', 
     yaxt='n')

lmphs <- 180 * sign(strain_scaling) + degadd
abline(h=lmphs, col=2, lty=2)
mtext("sign of\nlm coef.", side=4, at=lmphs, col=2, line=0.2, font=3)

mtext(sprintf('Relative Phase',ifelse(can.unwrap," (Unwrapped)","")), font=2, line=-4.0, adj=0.3)
axis(2)
axis(1, at=10**(-3:1), labels=paste0("(",c(1000,100,10,1,"1/10"),'s)'), line=1.6, lwd=0)
box()

