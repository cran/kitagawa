## ----include=FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----eval=TRUE, echo=FALSE-----------------------------------------
library(knitr)
options(width=69)
knitr::opts_chunk$set(tidy = TRUE, size="small")

## ----eval=TRUE, echo=TRUE, label=PRELIM----------------------------
library(RColorBrewer)
Set1 <- brewer.pal(8, "Set1")
library(signal, warn.conflicts=FALSE)
library(kitagawa)

## ----eval=TRUE, echo=TRUE, label=PARAMS1---------------------------
S. <- 1e-5    # Storativity [nondimensional]
T. <- 1e-4    # Transmissivity [m**2 / s]
D. <- T./S.   # Diffusivity [m**2 / s]
Ta <- 50      # Aquifer thickness [m] #100
Hw <- z <- 50 # Depth to water table [m] #10

# Using ANO1 stats from Kit Tbl 1
Rc. <- 0.075		# Radius of cased portion of well [m]
Lc. <- 570		# Length of cased portion of well [m]
Rs. <- 0.135		# Radius of screened portion of well [m]
Ls. <- 15		# Length of screened portion of well [m]
Vw. <- sensing_volume(Rc., Lc., Rs., Ls.) 	# volume of fluid [m**3]
#
# parameters assumed by well_response:
#	rho=1000		# density of rock [kg/m**3]
#	Kf=2.2e9		# Bulk modulus of fluid [Pascals]
#	grav=9.81	# gravitational acceleration [m/s**2]
rhog <- 9.81*1000
# Kitagawa Fig 7: Ku B / Kw Aw = 3 => Aw==4.8 at 40GPa
Ku. <- 40e9		# Bulk modulus [Pascals]
B. <- 0.5		# Skemptons ratio [nondimensional]

## ----eval=TRUE, echo=TRUE, label=PARAMS2---------------------------
# Frequencies
Q <- 10**seq(-5,2,by=0.05)					# [nondimensional]
lQ <- log10(Q)
omega <- omega_norm(Q, z, D., invert=TRUE)		# [Hz]

Phase <- function(Z){
	Phs. <- Arg(Z) # will wrap to -pi/pi
	uPhs. <- signal::unwrap(Phs., tol=pi/30)
	return(data.frame(Phs=Phs., uPhs=uPhs.))
}
	
# Responses converted to pressure if TRUE
asP <- FALSE
ZasP <- FALSE

## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=KITRESP----
wrsp <- well_response(omega, T.=T., S.=S., Vw.=Vw., Rs.=Rs., Ku.=Ku., B.=B., Avs=1, Aw=1, as.pressure=asP)
plot(wrsp) # uses plot.wrsp method
crsp <- wrsp[["Response"]][,2]	# Complex response
kGain <- Mod(crsp)/Ku./B.			# Amplitude (or Gain)
kP <- Phase(crsp)				# Phase

## ----eval=TRUE, echo=FALSE, fig.height=7, fig.width=5.5, label=KITRESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
plot(lQ, (kGain), type="l", #ylim=c(0, 1.15),
     xaxt="n", #yaxs="i",
     lwd=2,
     ylab=expression( Z / Epsilon * B * kappa[u] ), 
     xlab="",
     main="")
log10_ticks()
mtext("Sealed Well Response (KITAGAWA): Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain", adj=0)
par(mar=c(4,4,1.1,2))
plot(lQ, kP$Phs*180/pi, type="l", lty=3, #ylim=c(-190, -130), 
	 xaxt="n",
     ylab=expression(Z ~ "rel." ~ Epsilon), 
     xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, kP$uPhs*180/pi, type="l", lwd=2)
log10_ticks()
mtext("(b) Phase", adj=0)

## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=COOPERRESP----
wrsp <- open_well_response(omega, T.=T., S.=S., Ta=Ta, Hw=Hw, 
	model = "cooper", as.pressure=ZasP)
plot(wrsp)
crsp <- wrsp[["Response"]][,2]
cGain <- Mod(crsp)
cP <- Phase(crsp)

## ----eval=TRUE, echo=FALSE, fig.height=7, fig.width=5.5, label=COOPERRESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
plot(lQ, cGain, type="l", #ylim=c(0, 1.15),
     xaxt="n", #yaxs="i",
     lwd=2,
     ylab=expression( Z / H ),
     xlab="",
     main="")
log10_ticks()
mtext("Open Well Response (COOPER): Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain", adj=0)
par(mar=c(4,4,1.1,2))
plot(lQ, cP$Phs*180/pi, type="l", lwd=2, #ylim=c(-190, -130), 
	 xaxt="n",
     ylab=expression(Z ~ "rel." ~ H), 
     xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, cP$uPhs*180/pi, type="l", lty=3)
log10_ticks()
mtext("(b) Phase", adj=0)

## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=HSIEHRESP----
wrsp <- open_well_response(omega, T.=T., S.=S.,  Ta=Ta, Hw=Hw, model = "hsieh", as.pressure=ZasP)
plot(wrsp)
crsp <- wrsp[["Response"]][,2]
hGain <- Mod(crsp)
hP <- Phase(crsp)

## ----eval=TRUE, echo=FALSE, fig.height=7, fig.width=5.5, label=HSIEHRESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
plot(lQ, hGain, type="l", #ylim=c(0, 1.15),
     xaxt="n", #yaxs="i",
     lwd=2,
     ylab=expression( Z / H ), 
     xlab="",
     main="")
log10_ticks()
mtext("Open Well Response (HSIEH): Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain", adj=0)
par(mar=c(4,4,1.1,2))
plot(lQ, hP$Phs*180/pi, type="l", lty=3, #ylim=c(-190, -130), 
	 xaxt="n",
     ylab=expression(Z ~ "rel." ~ H), 
     xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, hP$uPhs*180/pi, type="l", lwd=2)
log10_ticks()
mtext("(b) Phase", adj=0)

## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=LIURESP----
wrsp <- open_well_response(omega, T.=T., S.=S.,  Ta=Ta, Hw=Hw, model = "liu", as.pressure=ZasP)
plot(wrsp)
crsp <- wrsp[["Response"]][,2]
lGain <- Mod(crsp)
lP <- Phase(crsp)

## ----eval=TRUE, echo=FALSE, fig.height=7, fig.width=5.5, label=LIURESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
plot(lQ, lGain, type="l", #ylim=c(0, 1.15),
     xaxt="n", #yaxs="i",
     lwd=2,
     ylab=expression( Z / H ),
     xlab="",
     main="")
log10_ticks()
mtext("Open Well Response (LIU): Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain", adj=0)
par(mar=c(4,4,1.1,2))
plot(lQ, lP$Phs*180/pi, type="l", lwd = 2, #ylim=c(-190, -130), 
	 xaxt="n",
     ylab=expression(Z ~ "rel." ~ H),
     xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, lP$uPhs*180/pi, type="l", lty=3)
log10_ticks()
mtext("(b) Phase", adj=0)

## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=WANGRESP----
wrsp <- open_well_response(omega, T.=T., S.=S., leak = 1e-8,
                           model = "wang", as.pressure=asP)
plot(wrsp)
crsp <- wrsp[["Response"]][,2]
rGain <- Mod(crsp)
rP <- Phase(crsp)

## ----eval=TRUE, echo=TRUE, label=WANGRESPDAT-----------------------
Transmiss <- c(1e0, 1e-2, 1e-4, 1e-6, 1e-8)
Storativ  <- c(1e-2, 1e-4, 1e-6, 1e-8)
omeg      <- 1.9322736 / 86400 # M2 in Hz
leak      <- 10^seq(-11, -3, 0.2)

## ----eval=TRUE, echo=FALSE, fig.height=4.7, fig.width=4, label=WANGRESPFIG----

pt.cex = 0.35
layout(matrix(c(1,2), ncol=1), heights=c(0.5,0.5))
  par(mar=c(2,4,1,1), 
      oma=c(2,0.1,1,0.1), 
      tcl=-0.3,
      mgp=c(2.5, 0.5, 0), las=1)
plot(c(), c(),  log = 'x',
     xlim = range(leak), type='l', ylim = c(-80, 100),
     ylab = 'Phase shift (deg)', xaxt='n', frame=FALSE)
lats <- seq(-11,-3,by=2)
axis(1, at=10^lats, labels=parse(text=sprintf("10^%s",lats)))
axis(3, at=10^lats, labels=FALSE)
axis(4, labels=FALSE)
box()

for (i in seq_along(Transmiss)) {
  T. <- Transmiss[i]
  for (j in seq_along(Storativ)) {
    S. <- Storativ[j]
    wang  <- open_well_response(omeg, 
                                T., 
                                S.,
                                Rs. = 0.1,
                                model = 'wang', 
                                leak = leak,
                                freq.units = 'Hz', 
                                as.pressure = FALSE)
    wang <- wang[["Response"]][,2]
    
    col <- ifelse(i <=3, 'firebrick4', i)
    if(i == 4) col <- 'dodgerblue4'
    if(i == 5) col <- 'forestgreen'
    points(x = leak, Arg(wang) * 180/pi, type = 'o', col = col,
           pch = j, cex = pt.cex)
    
  }
}
legend('bottomright', 
       col = c('firebrick4', 'dodgerblue4', 'forestgreen'),
       lty = 1,
       legend = c('1.0 >= T >= 1e-4',
                  'T = 1e-6',
                  'T = 1e-8'),
       cex = 0.5)
  
plot(c(), c(),  log = 'xy',
     xlim = range(leak), type='l', ylim = c(1e-10, 1),
     ylab = 'Amplitude ratio',
     xlab = 'Specific Leakage (1/s)', axes=FALSE)
axis(1, at=10^lats, labels=parse(text=sprintf("10^%s",lats)))
lats.y <- seq(-10,-1,by=3)
axis(3, at=10^lats, labels=FALSE)
axis(2, at=10^lats.y, labels=parse(text=sprintf("10^%s",lats.y)))
axis(4, at=10^lats.y, labels=FALSE)
box()

for (i in seq_along(Transmiss)) {
  T. <- Transmiss[i]
  for (j in seq_along(Storativ)) {
    S. <- Storativ[j]
    wang  <- open_well_response(omeg, 
                                T., 
                                S.,
                                Rs. = 0.1,
                                model = 'wang', 
                                leak = leak,
                                freq.units = 'Hz', 
                                as.pressure = FALSE)
    wang <- wang[["Response"]][,2]
    
    col <- ifelse(i <=3, 'firebrick4', i)
    if(i == 4) col <- 'dodgerblue4'
    if(i == 5) col <- 'forestgreen'
    points(x = leak, Mod(wang), type = 'o', col = col, 
           pch = j, cex = pt.cex)
    
  }
}
legend('bottomleft', 
       pch = 1:4,
       legend = c('S = 1e-2',
                  'S = 1e-4',
                  'S = 1e-6',
                  'S = 1e-8'),
       pt.cex = pt.cex,
       cex = 0.5)
mtext("Leakage Coefficient (1/s)", side=1, line=2)


## ----eval=TRUE, echo=TRUE, fig.height=4.7, fig.width=4, label=ROJRESP----
wrsp <- open_well_response(omega, T.=T., S.=S., z=z, model = "rojstaczer", as.pressure=asP)
plot(wrsp)
crsp <- wrsp[["Response"]][,2]
rGain <- Mod(crsp)
rP <- Phase(crsp)

## ----eval=TRUE, echo=FALSE, fig.height=6, fig.width=5.5, label=ROJRESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)
plot(lQ, (rGain), type="l", 
	 #ylim=c(0, 1.15),
     xaxt="n", 
     lwd=2,
     ylab=expression( Z / Epsilon ),
     xlab="", 
     main="")
log10_ticks()
mtext("Open Well Response (ROJSTACZER): Harmonic Strain", font=2, line=1.5)
mtext("(a) Gain", adj=0)
# relative to static-confined areal strain response", adj=0)
par(mar=c(4,4,1.1,2))
plot(lQ, rP$Phs*180/pi, type="l", lty=3, ylim=c(130, 180), 
	 xaxt="n",
     ylab=expression(Z ~ "rel." ~ Epsilon), 
     xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
abline(h=-180, col="grey")
lines(lQ, rP$uPhs*180/pi, type="l", lwd=2)
log10_ticks()
mtext("(b) Phase", adj=0)

## ----eval=TRUE, echo=FALSE, fig.height=6, fig.width=5.5, label=ALLRESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)

# kitagawa and rojstaczer
alim <- range(pretty(log10(c(kGain,rGain))))
plot(lQ, log10(kGain), col=Set1[1],
	 type="l", ylim=alim, 
	 #ylim=c(-2.5, 0.2),
     yaxt="n", xaxt="n", 
     lwd=2,
     ylab=expression(log[10] ~ Z / Epsilon),
     xlab="",
     main="")
lines(lQ, log10(rGain), lwd=2, col=Set1[2])
legend("bottomright", 
	c("Kitagawa et al (2011) -- sealed","Rojstaczer et al (1988) -- open"), 
	lwd=3, col=Set1[1:2], bty="n")
log10_ticks(2, major.ticks=-5:5)
log10_ticks()
mtext("Harmonic Strain Well Responses", font=2, line=1.5)
mtext("(a) Gain", adj=0)

par(mar=c(4,4,1.1,2))
plot(lQ, kP$Phs*180/pi -180, lty=3, lwd=1.5, col=Set1[1],
	type="l", ylim=90*c(0,-1), 
	yaxt="n",
	xaxt="n",
    ylab=expression(Z ~ "rel. -180" ~ Epsilon), 
    xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
lbls <- ats <- seq(-90,90,by=15)
lbls[seq_along(lbls)%%2==0] <- ""
axis(2, at=ats, labels=lbls)
lines(lQ, rP$Phs*180/pi-180, lty=3, lwd=1.5, col=Set1[2])
# unwrapped phase
lines(lQ, kP$uPhs*180/pi-180, lwd=2, col=Set1[1])
lines(lQ, rP$uPhs*180/pi-180, lwd=2, col=Set1[2])
log10_ticks()
mtext("(b) Anti-Phase", adj=0)

## ----eval=TRUE, echo=FALSE, fig.height=6, fig.width=5.5, label=ALLORESPFIG----
par(mfrow=c(2,1), 
    mar=c(2.5,4,3.1,2), 
    oma=rep(0.1,4), 
    omi=rep(0.1,4), 
    las=1)

# cooper, liu, hsieh
alim <- range(pretty(log10(c(hGain,cGain,lGain))))
plot(lQ, log10(hGain), col=Set1[3],
	 type="l", ylim=alim,
     xaxt="n", yaxt="n", 
     lwd=2,
     ylab=expression(log[10] ~ Z/H),
     xlab="",
     main="")
lines(lQ, log10(cGain), lwd=2, col=Set1[1])
lines(lQ, log10(lGain), lwd=2, col=Set1[2])
legend("bottomleft", 
	c("Cooper et al (1965)","Liu et al (1989)","Hsieh et al (1987)"), 
	lwd=3, col=Set1[1:3], bty="n")
log10_ticks(2, major.ticks=-9:2)
log10_ticks()
mtext("Harmonic Pressure-head Well Responses (Open)", font=2, line=1.5)
mtext("(a) Gain", adj=0)

par(mar=c(4,4,1.1,2))
plot(lQ, cP$uPhs*180/pi, lty=3, lwd=1.5, col=Set1[1],
	type="l", ylim=185*c(0,-1), 
	yaxt="n",
	xaxt="n",
    ylab=expression(Z ~ "rel." ~ H), 
    xlab=expression("Dimensionless frequency," ~ Q ==z^2 * omega / 2 * D  ))
lbls <- ats <- seq(-180,180,by=30)
lbls[seq_along(lbls)%%2==0] <- ""
axis(2, at=ats, labels=lbls)
lines(lQ, lP$uPhs*180/pi, lty=3, lwd=1.5, col=Set1[2])
lines(lQ, hP$uPhs*180/pi, lty=3, lwd=1.5, col=Set1[3])
# unwrapped phase
lines(lQ, cP$Phs*180/pi, lwd=2, col=Set1[1])
lines(lQ, lP$Phs*180/pi, lwd=2, col=Set1[2])
lines(lQ, hP$Phs*180/pi, lwd=2, col=Set1[3])
log10_ticks()
mtext("(b) Phase", adj=0)

