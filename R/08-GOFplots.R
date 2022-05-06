# Goodness-of-fit
library (jagsUI)
load("outputs\\ipm-jags-imm-GOF.Rdata")
plot(out$sims.list$dmape.obs[,1], out$sims.list$dmape.rep[,1])

par(mfrow=c(1,2))
plot(out$sims.list$dd.obs, out$sims.list$dd.rep,
     main="Poisson regression model",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(200,400), ylim=c(200,400))
curve(1*x, from=0, to=400, add=T, lty=2)
bp <- round(mean(out$sims.list$dd.rep > out$sims.list$dd.obs),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")

hist(out$sims.list$tvm.rep, nclass=50,
     xlab="variance/mean ", main=NA, axes=FALSE)
abline(v=out$mean$tvm.obs, col="red")
axis(1); axis(2)


par(mfrow=c(1,2))
plot(out$sims.list$dmape.obs[,3], out$sims.list$dmape.rep[,3],
     main="State-space model",
     ylab="Discrepancy replicate values",
     xlab="Discrepancy observed values", 
     xlim=c(0,200000), ylim=c(0,200000))
curve(1*x, from=0, to=200000, add=T, lty=2)
bp <- round(mean(out$sims.list$dd.rep > out$sims.list$dd.obs),2)
loc <- ifelse(bp < 0.5, "topleft", "bottomright")
legend(loc, legend=bquote(p[B]==.(bp)), bty="n")