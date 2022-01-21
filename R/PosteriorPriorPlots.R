x <- seq(from=0, to=4, by=0.01)
df <- data.frame(x=x, prior=dnorm(x,0,1))
df2 <- data.frame(omega=out$sims.list$omega1)

# in R base
par(mar=c(5,5,2,5))
plot(df$x, df$prior, type="l", lwd=2, xlim=c(0,4),
     ylab="Prior density (black line)", main="")
par(new=TRUE)
plot(density(df2$omega), type="l", lwd=2, xlim=c(0,4),
     main="", col="blue",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = c(0,30,60))      # Add second axis
mtext("Posterior density (blue line)", side = 4, line = 3)  

##########
# plot number of immigrants each year
##########
NAd <- 
# in ggplot
library(ggplot2)
ggplot() +
  geom_density(data=df2, aes(x=omega)) +
  geom_line(data=df, aes(x=prior, y=x)) +
  scale_y_continuous(
  name = "Prior density", limits=c(0,0.4),
  sec.axis = sec_axis( trans=~.*10, name="Posterior density")
)
