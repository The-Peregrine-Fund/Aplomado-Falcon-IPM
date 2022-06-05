# Goodness-of-fit
library (jagsUI)

load("outputs\\ipm-jags-noimm-GOF-pois-excludetransl.Rdata")
mod <- "no immigration P"
load("outputs\\ipm-jags-noimm-GOF-pois-zip.Rdata")
mod <- "no immigration ZIP"
load("outputs\\ipm-jags-imm-GOF-pois-excludetransl.Rdata")
mod <- "immigration P"
load("outputs\\ipm-jags-imm-GOF-pois-zip.Rdata")
mod <- "immigration ZIP"

load("C:\\Users\\rolek.brian\\Documents\\GitHub\\Aplomado-Falcon-IPM\\outputs\\ipm-e.Rdata")
load("C:\\Users\\rolek.brian\\Documents\\GitHub\\Aplomado-Falcon-IPM\\outputs\\ipm-ie.Rdata")


par(mfrow=c(2,3), oma=c(0,0,5,0))
plot(out$sims.list$dd.obs, out$sims.list$dd.rep,
     main="Poisson regression model\nfor Productivity",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values",
     xlim=c(0,40), ylim=c(0,40),
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=40, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$dd.rep > out$sims.list$dd.obs),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")

hist(out$sims.list$tvm.rep, nclass=50,
     xlab="variance/mean ", main=NA, axes=FALSE)
abline(v=out$mean$tvm.obs, col="red")
axis(1); axis(2)

plot(jitter(out$sims.list$dmape.obs[,1], amount=300), 
     jitter(out$sims.list$dmape.rep[,1], amount=300),
     main="State-space model\n for Breeder counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,8000), ylim=c(0,8000), 
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=8000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$dmape.rep[,1] > out$sims.list$dmape.obs[,1]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")


plot(jitter(out$sims.list$dmape.obs[,2], amount=300), 
     jitter(out$sims.list$dmape.rep[,2], amount=300),
     main="State-space model\n for non-breeder counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,20000), ylim=c(0,20000), 
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=30000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$dmape.rep[,2] > out$sims.list$dmape.obs[,2]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")


plot(jitter(out$sims.list$dmape.obs[,3], amount=200), 
     jitter(out$sims.list$dmape.rep[,3], amount=200),
     main="State-space model\n for first-year counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,15000), ylim=c(0,15000),
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=15000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$dmape.rep[,3] > out$sims.list$dmape.obs[,3]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")

mtext(mod, outer=T, side=3, xpd=NA)



par(mfrow=c(2,2))
plot(jitter(out$sims.list$tturn.obs[,1], amount=1), 
     jitter(out$sims.list$tturn.rep[,1], amount=1),
     main="State-space model\n for breeder counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,40), ylim=c(0,40),
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=6000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$tturn.rep[,1] > out$sims.list$tturn.obs[,1]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")

plot(jitter(out$sims.list$tturn.obs[,2], amount=3), 
     jitter(out$sims.list$tturn.rep[,2], amount=3),
     main="State-space model\n for non-breeder counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,20), ylim=c(0,30),
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=6000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$tturn.rep[,2] > out$sims.list$tturn.obs[,2]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")

plot(jitter(out$sims.list$tturn.obs[,3], amount=3), 
     jitter(out$sims.list$tturn.rep[,3], amount=3),
     main="State-space model\n for first-year counts",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,20), ylim=c(0,30),
     pch=16, cex=0.1, col="gray10")
curve(1*x, from=0, to=6000, add=T, lty=2, lwd=2, col="blue")
bp <- round(mean(out$sims.list$tturn.rep[,3] > out$sims.list$tturn.obs[,3]),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")


