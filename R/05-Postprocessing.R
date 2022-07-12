## ---- ipm postprocess 1 --------
load("outputs\\ipm-e.Rdata")
# load(".\\ipm-jags-no-imm.Rdata")
library (jagsUI)
# plots
library (ggplot2)
library (gridExtra)
library (tidybayes)
# Function to compute highest density interval. From Kruschke 2011.
source("R\\HDIofMCMC.R")

data_summary <- function(x) {
  m <- median(x)
  ymin <- HDIofMCMC(x, credMass=0.85)[[1]]
  ymax <- HDIofMCMC(x, credMass=0.85)[[2]]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

data_summary2 <- function(x) {
  m <- median(x)
  ymin <- HDIofMCMC(x, credMass=0.95)[[1]]
  ymax <- HDIofMCMC(x, credMass=0.95)[[2]]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

O.surv <- O.trans <- A.surv <- A.trans <- B.surv <- B.trans <- array(NA, dim=c(dim(out$sims.list$OSalpha1)[1],3))
B.recap <-array(NA, dim=c(dim(out$sims.list$mu.pB1)[1],2))
A.recap <-array(NA, dim=c(dim(out$sims.list$mu.pA1)[1],2))
##########################
# Juveniles
#########################
# Survival
# female - male
O.surv[,2]<- (out$sims.list$OSalpha1[,1,1]+out$sims.list$OSalpha1[,1,2])/2 - 
              (out$sims.list$OSalpha1[,2,1]+out$sims.list$OSalpha1[,2,2])/2 
# wild -hacked
O.surv[,3]<- (out$sims.list$OSalpha1[,1,1]+out$sims.list$OSalpha1[,2,1])/2 - 
  (out$sims.list$OSalpha1[,1,2]+out$sims.list$OSalpha1[,2,2])/2
# Recruitment
# female - male
O.trans[,2]<- out$sims.list$OBRalpha1[,1] - out$sims.list$OBRalpha1[,2]

##########################
# Floaters
#########################
# Recruitment
# wild - hacked
A.trans[,1]<- (out$sims.list$ABRalpha1[,1,1] + out$sims.list$ABRalpha1[,2,1])/2 -
  (out$sims.list$ABRalpha1[,1,2] + out$sims.list$ABRalpha1[,2,2])/2
# female - male
A.trans[,2]<- (out$sims.list$ABRalpha1[,1,1] + out$sims.list$ABRalpha1[,1,2])/2 -
  (out$sims.list$ABRalpha1[,2,1] + out$sims.list$ABRalpha1[,2,2])/2
# Recapture prob  
# wild - hacked
# mu.pF1 [effort, h]
A.recap[,1] <- (out$sims.list$mu.pA1[,1,1] + out$sims.list$mu.pA1[,2,1])/2 - 
              (out$sims.list$mu.pA1[,1,2] + out$sims.list$mu.pA1[,2,2])/2
# high effort - low effort
A.recap[,2] <- (out$sims.list$mu.pA1[,2,1] + out$sims.list$mu.pA1[,2,2])/2 - 
            (out$sims.list$mu.pA1[,1,1] + out$sims.list$mu.pA1[,1,2])/2
# mu.pB[effort[t]]
# high effort - low effort
B.recap[,1] <- out$sims.list$mu.pB1[,2] - out$sims.list$mu.pB1[,1] 

tab <- data.frame(stage= c("First-year", "First-year", "First-year", "Nonbreeder", "Nonbreeder", "Nonbreeder", "Nonbreeder", "Breeder"),
                  param= c("Survival", "Survival", "Transition", "Transition", "Transition", "Resight", "Resight", "Resight"),
                  comp= c("Sex", "Hacked",  "Sex", "Hacked", "Sex", "Hacked", "Effort", "Effort" ),
                  median=NA, LHDI95=NA, UHDI95=NA, P=NA, important=NA)
tab[1, 4:6]<- round(data_summary2(O.surv[,2]),3)
tab[2, 4:6]<- round(data_summary2(O.surv[,3]),3)
tab[3, 4:6]<- round(data_summary2(O.trans[,2]),3)
tab[4, 4:6]<- round(data_summary2(A.trans[,1]),3)
tab[5, 4:6]<- round(data_summary2(A.trans[,2]),3)
tab[6, 4:6]<- round(data_summary2(A.recap[,1]),3)
tab[7, 4:6]<- round(data_summary2(A.recap[,2]),3)
tab[8, 4:6]<- round(data_summary2(B.recap[,1]),3)
tab$important <- ifelse( ((tab$LHDI95 >= 0 & tab$UHDI95>=0) |
                             (tab$LHDI95 <= 0 & tab$UHDI95<=0)) &
                            (tab$LHDI95 != 0 & tab$UHDI95!=0), 
                          "yes", "no") 
tab$P[1] <- round(sum(O.surv[,2]>0)/length(O.surv[,2]),2)
tab$P[2] <- round(sum(O.surv[,3]>0)/length(O.surv[,3]),2)
tab$P[3] <- round(sum(O.trans[,2]<0)/length(O.trans[,2]),2)
tab$P[4] <- round(sum(A.trans[,1]>0)/length(A.trans[,1]),2)
tab$P[5] <- round(sum(A.trans[,2]<0)/length(A.trans[,2]),2)
tab$P[6] <- round(sum(A.recap[,1]>0)/length(A.recap[,1]),2)
tab$P[7] <- round(sum(A.recap[,2]>0)/length(A.recap[,2]),2)
tab$P[8] <- round(sum(B.recap[,1]>0)/length(B.recap[,1]),2)
print(tab)

## ---- ipm postprocess 2 --------
# combine survival data
temp.df<- data.frame(Draws=O.surv[,3], Cat="Hacked", State="First-year")
temp.df1<- data.frame(Draws=O.surv[,2], Cat="Sex", State="First-year")
temp.df2<- data.frame(Draws=0, Cat="Hacked", State="Breeder")
temp.df3<- data.frame(Draws=0, Cat="Sex", State="Breeder")
temp.df4<- data.frame(Draws=0, Cat="Hacked", State="Non-breeder")
temp.df5<- data.frame(Draws=0, Cat="Sex", State="Non-breeder")
df.surv <- rbind(temp.df, temp.df1)
df.surv$combined<- factor(paste(df.surv$State, df.surv$Cat))
df.surv$combined <- factor(df.surv$combined, levels=levels(df.surv$combined)[c(3,4,5,6,1,2)])

# combine transitions data
temp.df<- data.frame(Draws=0, Cat="Hacked", State="First-year")
temp.df1<- data.frame(Draws=0, Cat="Sex", State="First-year")
temp.df2<- data.frame(Draws=0, Cat="Hacked", State="Breeder")
temp.df3<- data.frame(Draws=0, Cat="Sex", State="Breeder")
temp.df4<- data.frame(Draws=A.trans[,1], Cat="Hacked", State="Non-breeder")
temp.df5<- data.frame(Draws=A.trans[,2], Cat="Sex", State="Non-breeder")
df.trans<- rbind(temp.df4, temp.df5)
df.trans$combined<- factor(paste(df.trans$State, df.trans$Cat))
df.trans$combined <- factor(df.trans$combined, levels=levels(df.trans$combined)[c(3,4,5,6,1,2)])

# combine resight data
temp.df2<- data.frame(Draws=0, Cat="Hacked", State="Breeder")
temp.df3<- data.frame(Draws=0, Cat="Sex", State="Breeder")
temp.df4<- data.frame(Draws=B.recap[,1], Cat="Effort", State="Breeder")
temp.df5<- data.frame(Draws=0, Cat="Hacked", State="Non-breeder")
temp.df6<- data.frame(Draws=0, Cat="Sex", State="Non-breeder")
temp.df7<- data.frame(Draws=A.recap[,2], Cat="Effort", State="Non-breeder")
df.recap<- rbind(temp.df4, temp.df7)
df.recap$combined<- factor(paste(df.recap$State, df.recap$Cat))
df.recap$combined <- factor(df.recap$combined, levels=c(levels(df.recap$combined)[c(4,5,6,1,2,3)])) #"Juvenile Hacked", "Juvenile Sex", 
# fix level order for plots
ord <- c("First-year","Non-breeder","Breeder")
df.surv$State <-  factor(df.surv$State, levels=ord, labels=ord) 
df.trans$State <-  factor(df.trans$State, levels=ord, labels=ord) 
df.recap$State <-  factor(df.recap$State, levels=ord[-1], labels=ord[-1]) 

# plots
library (ggplot2)
library (gridExtra)
source("R/HDIofMCMC.R")
data_summary <- function(x) {
  m <- median(x)
  ymin <- HDIofMCMC(x, credMass=0.85)[[1]]
  ymax <- HDIofMCMC(x, credMass=0.85)[[2]]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

data_summary2 <- function(x) {
  m <- median(x)
  ymin <- HDIofMCMC(x, credMass=0.95)[[1]]
  ymax <- HDIofMCMC(x, credMass=0.95)[[2]]
  return(c(y=m,ymin=ymin,ymax=ymax))
}

df.surv$Parameter <- "Survival"
df.trans$Parameter <- "Recruitment"
df.recap$Parameter <- "Resight"
df.all <- rbind(df.surv, df.trans, df.recap)
df.all$pc <- paste(df.all$Parameter, df.all$Cat)

txt <- 30
lwd <- 2
lwd2 <- 1
colors <- c("Survival" = "#f7f7f7", "Recruitment" = "#cccccc", "Resight" = "#969696")
colors <- c("First-year" = "#7fc97f", "Non-breeder" = "#beaed4", "Breeder" = "#fdc086")
f_labels <- data.frame(Parameter=c("Survival", "Recruitment", "Resight"),
                      label=c("Greater wild-hatched \nand females", "", ""))
pall <- ggplot(df.all, aes(x = combined, y = Draws)) +
  scale_y_continuous(breaks=c(-1.0, -0.5, 0, 0.5, 1.0),  
                     labels=c(-1.0, "", 0, "", 1.0), limits = c(-1, 1.2)) +
  geom_hline(yintercept=0, linetype="solid", size=lwd2) +
  stat_halfeye(aes(fill=State), point_interval = median_hdi,
               .width = c(0.85, 0.95), interval_size_domain=c(0.1,4),
               fatten_point=2) +
  scale_fill_manual(name='Life stage', values=colors) +
  #stat_summary(fun.data=data_summary,  geom="pointrange", size=lwd) +
  #stat_summary(fun.data=data_summary2,  geom="pointrange", size=lwd2) +
  theme_classic() + 
  theme (text = element_text(size=txt)) +  
  ylab("Difference") + xlab ("") +
  scale_x_discrete(breaks=levels(factor(df.all$combined)),
                   labels=c("Hacked","Sex", "Hacked", "Sex", "Effort", "Effort" )) + 
  facet_wrap(~Parameter, scales="free_x") 
  #geom_text(x = 0.5, y = 0.9, aes(label = label), data = f_labels)
pall 
tiff(".//figs//SurvivalDiffs_IPM_Reduced.tiff", 
    height=6, width=12, res=300, units="in")
pall 
dev.off()


## ---- ipm postprocess 3 --------
########################
# Plot mean estimates with 95% HDIs 
#######################
df4<- df3<- df2 <- df1 <- data.frame(param=NA , md=NA, lhdi85=NA, uhdi85=NA, lhdi95=NA, uhdi95=NA)
# params with 1 estimate
for (i in 1:3){
ind <- c(17,18,21)[i]
df1[i,1] <- names(out$sims.list)[ind]
df1[i,2] <- median(out$sims.list[[ind]] )
df1[i,3:4] <- HDIofMCMC(out$sims.list[[ind]], credMass=0.85)
df1[i,5:6] <- HDIofMCMC(out$sims.list[[ind]])
}

# params with 2 estimates
ind <- 19
df2[1,1] <- paste(names(out$sims.list)[ind], "_1", sep="")
df2[1,2] <- median(out$sims.list[[ind]][,1] )
df2[1,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1], credMass=0.85)
df2[1,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1], credMass=0.95)
df2[2,1] <- paste(names(out$sims.list)[ind], "_2", sep="")
df2[2,2] <- median(out$sims.list[[ind]][,2] )
df2[2,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2], credMass=0.85)
df2[2,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2], credMass=0.95)
ind <- 23
df2[3,1] <- paste(names(out$sims.list)[ind], "_1", sep="")
df2[3,2] <- median(out$sims.list[[ind]][,1] )
df2[3,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1], credMass=0.85)
df2[3,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1], credMass=0.95)
df2[4,1] <- paste(names(out$sims.list)[ind], "_2", sep="")
df2[4,2] <- median(out$sims.list[[ind]][,2] )
df2[4,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2], credMass=0.85)
df2[4,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2], credMass=0.95)

ind <- 16
df4[1,1] <- paste(names(out$sims.list)[ind], "_1_1", sep="")
df4[1,2] <- median(out$sims.list[[ind]][,1,1] )
df4[1,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.85)
df4[1,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.95)
df4[2,1] <- paste(names(out$sims.list)[ind], "_2_1", sep="")
df4[2,2] <- median(out$sims.list[[ind]][,2,1] )
df4[2,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.85)
df4[2,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.95)
df4[3,1] <- paste(names(out$sims.list)[ind], "_1_2", sep="")
df4[3,2] <- median(out$sims.list[[ind]][,1,2] )
df4[3,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.85)
df4[3,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.95)
df4[4,1] <- paste(names(out$sims.list)[ind], "_2_2", sep="")
df4[4,2] <- median(out$sims.list[[ind]][,2,2] )
df4[4,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.85)
df4[4,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.95)
ind <- 20
df4[5,1] <- paste(names(out$sims.list)[ind], "_1_1", sep="")
df4[5,2] <- median(out$sims.list[[ind]][,1,1] )
df4[5,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.85)
df4[5,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.95)
df4[6,1] <- paste(names(out$sims.list)[ind], "_2_1", sep="")
df4[6,2] <- median(out$sims.list[[ind]][,2,1] )
df4[6,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.85)
df4[6,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.95)
df4[7,1] <- paste(names(out$sims.list)[ind], "_1_2", sep="")
df4[7,2] <- median(out$sims.list[[ind]][,1,2] )
df4[7,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.85)
df4[7,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.95)
df4[8,1] <- paste(names(out$sims.list)[ind], "_2_2", sep="")
df4[8,2] <- median(out$sims.list[[ind]][,2,2] )
df4[8,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.85)
df4[8,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.95)
ind <- 22
df4[9,1] <- paste(names(out$sims.list)[ind], "_1_1", sep="")
df4[9,2] <- median(out$sims.list[[ind]][,1,1] )
df4[9,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.85)
df4[9,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,1], credMass=0.95)
df4[10,1] <- paste(names(out$sims.list)[ind], "_2_1", sep="")
df4[10,2] <- median(out$sims.list[[ind]][,2,1] )
df4[10,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.85)
df4[10,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,1], credMass=0.95)
df4[11,1] <- paste(names(out$sims.list)[ind], "_1_2", sep="")
df4[11,2] <- median(out$sims.list[[ind]][,1,2] )
df4[11,3:4] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.85)
df4[11,5:6] <- HDIofMCMC(out$sims.list[[ind]][,1,2], credMass=0.95)
df4[12,1] <- paste(names(out$sims.list)[ind], "_2_2", sep="")
df4[12,2] <- median(out$sims.list[[ind]][,2,2] )
df4[12,3:4] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.85)
df4[12,5:6] <- HDIofMCMC(out$sims.list[[ind]][,2,2], credMass=0.95)


df<- rbind(df1,df2,df4) 
df<- df[c(8,9,10,11,1,2, # survival
          4,5,12,13,14,15,3, # recruitment
          16,17,18,19,6,7),] # detection
state <- c("First-year\nfemale wild", "First-year\nmale wild", "First-year\nfemale hacked", "First-year\nmale hacked", "Nonbreeder", "Breeder", 
  "First-year to\nbreeder\nfemale", "First-year to\nbreeder\nmale", "Nonbreeder\nto breeder\nfemale wild" , "Nonbreeder\nto breeder\nmale wild", "Nonbreeder\nto breeder\nfemale hacked", "Nonbreeder\nto breeder\nmale hacked", "Breeder\nto\nnonbreeder",
  "Nonbreeder\nwild\nlow effort", "Nonbreeder\nwild\nhigh effort", "Nonbreeder\nhacked\nlow effort", "Nonbreeder\nhacked\nhigh effort","Breeder\nlow effort", "Breeder\nhigh effort")
statenum <- 1:length(state)

lwd <- 1.5
lwd2 <- 0.5
pcex <- 1
tcl=-0.2

tiff(file=".\\figs\\Survival_estimates.tiff", width=3.25, height=2.3, res=300, unit="in")
par(mfrow=c(3,1), mar=c(1.4,1,0.3,0.3), oma=c(1,1,0,0)) #bottom, left, top, right
plot(statenum[1:6], df$md[1:6], ylim=c(0,1), xlim=c(0.5,6.5), 
     xaxt="n", yaxt="n", xlab="", ylab="Probability",
     cex=pcex)
title("Survival", line=-1, cex.main=0.8, font.main=1)
arrows(x0=statenum[1:6], y0=df$lhdi85[1:6], y1=df$uhdi85[1:6], length=0, lwd=lwd)
arrows(x0=statenum[1:6], y0=df$lhdi95[1:6], y1=df$uhdi95[1:6], length=0, lwd=lwd2)
axis(1, at=statenum[1:6], labels=state[1:6], cex.axis=0.5, padj=c(-1,-1,-1,-1,-4,-4), tcl=tcl)
axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=0.5, tcl=tcl, padj=2.5)

plot(statenum[7:13], df$md[7:13], ylim=c(0,1.3), xlim=c(6.5,13.5), 
     xaxt="n", yaxt="n",xlab="", ylab="Probability",
     cex=pcex)
title("Recruitment", line=-1, cex.main=0.8, font.main=1)
arrows(x0=statenum[7:13], y0=df$lhdi85[7:13], y1=df$uhdi85[7:13], length=0, lwd=lwd)
arrows(x0=statenum[7:13], y0=df$lhdi95[7:13], y1=df$uhdi95[7:13], length=0, lwd=lwd2)
axis(1, at=statenum[7:13], labels=state[7:13], cex.axis=0.5,  tcl=tcl, padj=-0.2)
axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=0.5, tcl=tcl, padj=2.5)

plot(statenum[14:19], df$md[14:19], ylim=c(0,1), xlim=c(13.5,19.5), 
     xaxt="n", yaxt="n", xlab="Estimate", ylab="Probability", 
     cex=pcex)
title("Resight", line=-1, cex.main=0.8, font.main=1)
arrows(x0=statenum[14:19], y0=df$lhdi85[14:19], y1=df$uhdi85[14:19], length=0, lwd=lwd)
arrows(x0=statenum[14:19], y0=df$lhdi95[14:19], y1=df$uhdi95[14:19], length=0, lwd=lwd2)
axis(1, at=statenum[14:19], labels=state[14:19], cex.axis=0.5, tcl=tcl, padj=0)
axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), cex.axis=0.5, tcl=tcl, padj=2.5)
mtext("Probability", side =2, outer=T, cex=0.5)
dev.off()

###############################
# Plot time series
##############################
#######################
# Plot survival over time
######################
library (ggplot2)
library (gridExtra)
library (ggmcmc)
source(".\\R\\HDIofMCMC.R")
logit <- function(p=p){
  log(p/(1-p))
}
data_summary <- function(x) {
  m <- plogis(apply(logit(x), 2, median) )
  ymin <- plogis(apply(logit(x), 2, HDIofMCMC, credMass=0.95))[1,]
  ymax <- plogis(apply(logit(x), 2, HDIofMCMC, credMass=0.95))[2,]
  return(list(y=m,ymin=ymin,ymax=ymax))
}
data_summary2 <- function(x) {
  m <- apply(x, 2, median) 
  ymin <- apply(x, 2, HDIofMCMC, credMass=0.95)[1,]
  ymax <- apply(x, 2, HDIofMCMC, credMass=0.95)[2,]
  return(list(y=m,ymin=ymin,ymax=ymax))
}

tiff(file=".\\figs\\Survival_overTime.tiff", 
     width=3.14, height=3.1, units="in", res=300)
par(mfrow=c(3,2), mar=c(1,1,0.3,0.3), oma=c(1,3,0,0))
ca <- 1
# Survival first year wild female
dat<- data_summary(out$sims.list$eta.OSalpha[,1,1,])
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), tcl=-0.2, cex.axis=0.5,
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=ca, 
     padj=0.5, las=1, hadj=0.5)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=ca, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year\nwild female", cex.main=1, line=-2)

# Survival first-year wild male
dat<- data_summary(out$sims.list$eta.OSalpha[,2,1,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), 
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c("","",""), tcl=-0.2, cex.axis=ca, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=ca, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year\nwild male", cex.main=1, line=-2)

# Survival first-year hacked female
dat<- data_summary(out$sims.list$eta.OSalpha[,1,2,])
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), tcl=-0.2, cex.axis=0.5,
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=ca, 
     padj=0.5, las=1, hadj=0.5)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=ca, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year\nhacked female", cex.main=1, line=-2)

# Survival first-year hacked male
dat<- data_summary(out$sims.list$eta.OSalpha[,2,2,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), 
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c("","",""), tcl=-0.2, cex.axis=ca, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=ca, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year\nhacked male", cex.main=1, line=-2)

# Survival nonbreeder
dat<- data_summary(out$sims.list$eta.ASalpha[,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), 
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=ca, 
     padj=0.5, las=1, hadj=0.5)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=ca, padj=-1)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder", cex.main=1, line=-1)

# Survival Breeder
dat<- data_summary(out$sims.list$eta.BSalpha[,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1), 
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c("","",""), tcl=-0.2, cex.axis=ca, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=ca, padj=-1)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Breeder", cex.main=1, line=-1)
mtext("Survival probability", side=2, outer=T, cex=1, padj=-1.25)
dev.off()

###########
# Transitions 1
###########
tiff(file=".\\figs\\Recruit_1_overTime.tiff", 
     width=3.14, height=2.3, units="in", res=300)
par(mfrow=c(2,2), mar=c(1,1,0.3,0.3), oma=c(0,1,0,0))
# Transition Juvenile female 
dat<- data_summary(out$sims.list$eta.OBRalpha[,1,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year female to breeder", cex.main=0.5, line=-0.5)

# Transition Juvenile male
dat<- data_summary(out$sims.list$eta.OBRalpha[,2,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("First-year male to breeder", cex.main=0.5, line=-0.5)

# Transition breeder to nonbreeder
dat<- data_summary(out$sims.list$eta.BARalpha[,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Breeder to nonbreeder", cex.main=0.5, line=-0.5)
mtext("Recruitment probability", side=2, outer=T, cex=0.5)
dev.off()

##############
# Transitions 2
################
tiff(file=".\\figs\\Recruit_2_overTime.tiff", 
     width=3.14, height=2.3, units="in", res=300)
par(mfrow=c(2,2), mar=c(1,1,0.3,0.3), oma=c(0,1,0,0))
# Transition Nonbreeder female wild
dat<- data_summary(out$sims.list$eta.ABRalpha[,1,1,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder female wild to breeder", cex.main=0.5, line=-0.5)

# Transition Nonbreeder male wild
dat<- data_summary(out$sims.list$eta.ABRalpha[,2,1,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder male wild to breeder", cex.main=0.5, line=-0.5)

# Transition Nonbreeder female hacked
dat<- data_summary(out$sims.list$eta.ABRalpha[,1,2,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder female hacked to breeder", cex.main=0.5, line=-0.5)

# Transition Nonbreeder male hacked
dat<- data_summary(out$sims.list$eta.ABRalpha[,2,2,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder male hacked to breeder", cex.main=0.5, line=-0.5)
mtext("Recruitment probability", side=2, outer=T, cex=0.5)
dev.off()

################
# Resight 
#################
tiff(file=".\\figs\\Resight_overTime.tiff", 
     width=3.14, height=2.3, units="in", res=300)
par(mfrow=c(2,2), mar=c(1,1,0.3,0.3), oma=c(0,1,0,0))
# Recapture nonbreeder wild
dat<- data_summary(out$sims.list$eta.pA[,1,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder wild", cex.main=0.5, line=-0.5)

# Recapture nonbreeder hacked
dat<- data_summary(out$sims.list$eta.pA[,2,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Nonbreeder hacked", cex.main=0.5, line=-0.5)

# Recapture breeder 
dat<- data_summary(out$sims.list$eta.pB[,]) 
plot( 1994:2018, dat$y, type="n", ylab="Probability", xlab="Year", ylim=c(0,1),
      xaxt="n", yaxt="n")
axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), tcl=-0.2, cex.axis=0.5, padj=3)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.5, padj=-4)
polygon(c(1994:2018, 2018:1994), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1994:2018, dat$y)
title ("Breeder", cex.main=0.5, line=-0.5)
mtext("Resight probability", side=2, outer=T, cex=0.5)
dev.off()

###################
# Fecundity 
###########################
tiff(file=".\\figs\\Fecundity_overTime.tiff", 
     width=3.14, height=2.3, res=300, units="in")
par(mfrow=c(1,1), mar=c(2,1,0.3,0.3), oma=c(0,1,0,0))
# Fecundity breeder 
dat<- data_summary2(out$sims.list$F)  
plot( 1993:2018, dat$y, type="n", ylab="Fecundity", xlab="Year", ylim=c(0,max(dat$ymax)),
      xaxt="n", yaxt="n")
axis(2, at=c(0,1,2), labels=c(0,1,2), tcl=-0.2, cex.axis=1, hadj=-0.5, las=1)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), 
     tcl=-0.2, cex.axis=1, padj=-1.5)
polygon(c(1993:2018, 2018:1993), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1993:2018, dat$y)
abline(v=2012, lty=2)
mtext("Fecundity", side=2, cex=1, line=1)
dev.off()

#################################
# Population estimate plots
####################################
data_summary2 <- function(x) {
  m <- apply(x, 2, median) 
  ymin <- apply(x, 2, HDIofMCMC, credMass=0.95)[1,]
  ymax <- apply(x, 2, HDIofMCMC, credMass=0.95)[2,]
  return(list(y=m,ymin=ymin,ymax=ymax))
}

##############
# Abundance
###############
load(".\\data\\counts.Rdata")
load(".\\data\\data-7states.Rdata")
datl$countOM <- datl$countJM-datl$aug

tiff(file=".\\figs\\Abundance_overTime.tiff", 
     width=3.14, height=2.3, res=300, units="in")
par(mfrow=c(2,2), mar=c(0.5,1,0,0), oma=c(0.5,1,0.3,0.3))
dat<- data_summary2(out$sims.list$Ntot)  
plot( 1993:2018, dat$y, type="n", ylab="Abundance", xlab="Year", ylim=c(0,max(dat$ymax)),
      xaxt="n", yaxt="n")
axis(2, at=c(0,60,120), labels=c(0,60,120), tcl=-0.2, cex.axis=0.75, 
     padj=2,  las=0)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=0.75, padj=-4)
polygon(c(1993:2018, 2018:1993), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1993:2018, dat$y)
lines( 1993:2018, datl$countBM+datl$countFM+datl$countOM+datl$aug , lty=2)
title ("All males", cex.main=0.75, line=-0.5)

dat<- data_summary2(out$sims.list$N[,1,])  
plot( 1993:2018, dat$y, type="n", ylab="Abundance", xlab="Year", ylim=c(min(dat$ymin),max(dat$ymax)),
      xaxt="n", yaxt="n")
axis(2, at=c(0,15,30), labels=c(0,15,30), tcl=-0.2, cex.axis=0.75, padj=2)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c("","","","",""), tcl=-0.2, 
     cex.axis=0.75, padj=-4)
polygon(c(1993:2018, 2018:1993), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1993:2018, dat$y)
lines( 1993:2018, datl$countOM , lty=2)
title ("First-year males", cex.main=0.75, line=-0.5)

dat<- data_summary2(out$sims.list$NF)  
plot( 1993:2018, dat$y, type="n", ylab="Abundance", xlab="Year", ylim=c(min(dat$ymin),max(dat$ymax)),
      xaxt="n", yaxt="n")
axis(2, at=c(0,15,30), labels=c(0,15,30), tcl=-0.2, cex.axis=0.75, padj=2)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.75, padj=-2.5)
polygon(c(1993:2018, 2018:1993), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1993:2018, dat$y)
lines( 1993:2018, datl$countFM , lty=2)
title ("Nonbreeder males", cex.main=0.75, line=-0.5)

dat<- data_summary2(out$sims.list$NB)  
plot( 1993:2018, dat$y, type="n", ylab="Abundance", xlab="Year", ylim=c(min(dat$ymin),max(dat$ymax)),
      xaxt="n", yaxt="n")
axis(2, at=c(0,30,60), labels=c(0,30,60), tcl=-0.2, cex.axis=0.75, padj=2)
axis(1, at=c(1995,2000,2005,2010,2015), labels=c(1995,"",2005,"",2015), tcl=-0.2, 
     cex.axis=0.75, padj=-2.5)
polygon(c(1993:2018, 2018:1993), c(dat$ymax, rev(dat$ymin)), border=NA, col="gray")
lines( 1993:2018, dat$y )
lines( 1993:2018, datl$countBM , lty=2)
title ("Breeder males", cex.main=0.75, line=-0.5)
mtext("Abundance", side=2, outer=T, cex=1)
dev.off()

####################
## ---- ipm population growth rates -------
# Code adapted from Kery and Schaub 2012
####################
## ---- ipm population growth rates -------
# Code adapted from Kery and Schaub 2012
load("outputs\\ipm-em.Rdata")
source(".\\HDIofMCMC.R")
load(".\\data\\counts.Rdata")
aug<- df.F$male_hacked
nyears <- 26
lam <- array(NA, dim = c(nyears-1,3), dimnames=list(1:25, c("md", "lhdi", "uhdi")))
fit <- array(NA, dim = c(nyears-1,8,3), dimnames=list(1:(nyears-1), 
                                                      c("phiO_fem_wild", "phiO_male_wild", "phiO_fem_hacked", "phiO_male_hacked", "phiF", "phiB", "fec", "hacked"), c("md", "lhdi", "uhdi")))

data_summary <- function(x) {
  m <- apply(x, 2, median) 
  ymin <- apply(x, 2, HDIofMCMC, credMass=0.95)[1,]
  ymax <- apply(x, 2, HDIofMCMC, credMass=0.95)[2,]
  return(list(y=m,ymin=ymin,ymax=ymax))
}

lam.sims <- array(NA, dim=c(nrow(out$sims.list$Ntot), ncol(out$sims.list$Ntot)-1 ) )
for (t in 1:ncol(lam.sims)){    
  lam.sims[,t] <-  out$sims.list$Ntot[,t+1]/(out$sims.list$Ntot[,t])
} #t
# skip first year bc 0/0 returns NAs and screws up functions
lam[,1] <- data_summary(lam.sims)$y
lam[,2] <- data_summary(lam.sims)$ymin
lam[,3] <- data_summary(lam.sims)$ymax

# Survival FY wild female
fit[,1,1] <- data_summary(out$sims.list$eta.OSalpha[,1,1,])$y
fit[,1,2] <- data_summary(out$sims.list$eta.OSalpha[,1,1,])$ymin
fit[,1,3] <- data_summary(out$sims.list$eta.OSalpha[,1,1,])$ymax

# Survival FY wild male
fit[,2,1] <- data_summary(out$sims.list$eta.OSalpha[,2,1,])$y
fit[,2,2] <- data_summary(out$sims.list$eta.OSalpha[,2,1,])$ymin
fit[,2,3] <- data_summary(out$sims.list$eta.OSalpha[,2,1,])$ymax

# Survival FY hacked female
fit[,3,1] <- data_summary(out$sims.list$eta.OSalpha[,1,2,])$y
fit[,3,2] <- data_summary(out$sims.list$eta.OSalpha[,1,2,])$ymin
fit[,3,3] <- data_summary(out$sims.list$eta.OSalpha[,1,2,])$ymax

# Survival FY hacked male
fit[,4,1] <- data_summary(out$sims.list$eta.OSalpha[,2,2,])$y
fit[,4,2] <- data_summary(out$sims.list$eta.OSalpha[,2,2,])$ymin
fit[,4,3] <- data_summary(out$sims.list$eta.OSalpha[,2,2,])$ymax

# Survival Nonbreeder
fit[,5,1] <- data_summary(out$sims.list$eta.ASalpha[,])$y
fit[,5,2] <- data_summary(out$sims.list$eta.ASalpha[,])$ymin
fit[,5,3] <- data_summary(out$sims.list$eta.ASalpha[,])$ymax

# Survival Nonbreeder
fit[,6,1] <- data_summary(out$sims.list$eta.BSalpha[,])$y
fit[,6,2] <- data_summary(out$sims.list$eta.BSalpha[,])$ymin
fit[,6,3] <- data_summary(out$sims.list$eta.BSalpha[,])$ymax

# Fecundity 
fit[,7,1] <- data_summary(out$sims.list$F[,1:25])$y
fit[,7,2] <- data_summary(out$sims.list$F[,1:25])$ymin
fit[,7,3] <- data_summary(out$sims.list$F[,1:25])$ymax

# Number hacked 
fit[,8,1] <- aug[1:25]

# Calculate some correlation coefficients
l.aug<- ifelse(aug[1:25]==0, log(0.001), log(aug[1:25]))
n.iter <- nrow(out$sims.list$mu.F)
correl.h <- array(NA, dim = c(n.iter,8))
for (i in 1:n.iter){
  correl.h[i,1] <- cor(lam.sims[i,], out$sims.list$eta.OSalpha[i,1,1,1:25])
  correl.h[i,2] <- cor(lam.sims[i,], out$sims.list$eta.OSalpha[i,2,1,1:25])
  correl.h[i,3] <- cor(lam.sims[i,], out$sims.list$eta.OSalpha[i,1,2,1:25])
  correl.h[i,4] <- cor(lam.sims[i,], out$sims.list$eta.OSalpha[i,2,2,1:25])
  correl.h[i,5] <- cor(lam.sims[i,], out$sims.list$eta.ASalpha[i,1:25])
  correl.h[i,6] <- cor(lam.sims[i,], out$sims.list$eta.BSalpha[i,1:25])
  correl.h[i,7] <- cor(lam.sims[i,], out$sims.list$F[i,2:26])
  correl.h[i,8] <- cor(lam.sims[i,], l.aug)
}

# Credible intervals of correlation coefficients
correl.est <- array(NA, dim=c(8,3))
correl.est[1,2:3] <- HDIofMCMC(correl.h[,1], 0.90)
correl.est[2,2:3] <- HDIofMCMC(correl.h[,2], 0.90)
correl.est[3,2:3] <- HDIofMCMC(correl.h[,3], 0.90)
correl.est[4,2:3] <- HDIofMCMC(correl.h[,4], 0.90)
correl.est[5,2:3] <- HDIofMCMC(correl.h[,5], 0.90)
correl.est[6,2:3] <- HDIofMCMC(correl.h[,6], 0.90)
correl.est[7,2:3] <- HDIofMCMC(correl.h[,7], 0.90)
correl.est[8,2:3] <- HDIofMCMC(correl.h[,8], 0.90)

# Compute the posterior modes of correlation coefficients
for (j in 1:8){
  m <- density(correl.h[,j], na.rm = TRUE)
  correl.est[j,1]<- m$x[which(m$y==max(m$y))]
}
correl.est <- round(correl.est,2)

# Probability that correlation coefficients (r) > 0
n.iter <- nrow(out$sims.list$mu.F)
P <- c()
P[1] <- sum(correl.h[!is.na(correl.h[,1]),1]>0)/n.iter
P[2] <- sum(correl.h[!is.na(correl.h[,2]),2]>0)/n.iter
P[3] <- sum(correl.h[!is.na(correl.h[,3]),3]>0)/n.iter
P[4] <- sum(correl.h[!is.na(correl.h[,4]),4]>0)/n.iter
P[5] <- sum(correl.h[!is.na(correl.h[,5]),5]>0)/n.iter
P[6] <- sum(correl.h[!is.na(correl.h[,6]),6]>0)/n.iter
P[7] <- sum(correl.h[!is.na(correl.h[,7]),7]>0)/n.iter
P[8] <- sum(correl.h[!is.na(correl.h[,8]),8]>0)/n.iter
P <- round(P,2)

# Plot Fig. 9
library (viridis)
tiff(file="figs\\Corr_plots.tiff", 
     width=4.0, height=2.3, res=600, units="in")
txt <- 0.5
txt2 <- 0.6
pt.sz <- 0.5
mar<- c(4, 1, 0.5, 0.5)
lwd <- 0.7
par(mfrow = c(2, 4), mar = mar, mgp=c(1, 0.6, 0), las = 1, oma=c(0,3,0,0), cex = txt)
linecol <- c("grey50")

par(mar = mar)
plot(y = lam[,1], fit[,"phiB",1], type = "n", 
     xlim = c(0.3, 1.0), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3))
axis(1, at=c(0.3, 0.65, 1), labels=c("0.3","0.65", "1"))
mtext(text="Breeder survival", side=1, line=1.5, cex=txt)
segments(x0=fit[,"phiB",2], y0=lam[,1], x1=fit[,"phiB",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiB",1], y0=lam[,2], x1=fit[,"phiB",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiB",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0.3, y = 2.8, paste("r = ",correl.est[6,1], " (",correl.est[6,2],", ", correl.est[6,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0.3, y = 2.3, paste("P(r>0) = ", P[6], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

par(mar = mar)
plot(y = lam[,1], l.aug, type = "n", 
     xlim = c(log(0.001), log(70)), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=F)
axis(1, at=c(log(0.001), log(1), log(70)), labels=c(0, 1, 70))
mtext(text="Number hacked\n(log)", side=1, line=2.5, cex=txt)
segments(x0=l.aug, y0=lam[,2], x1=l.aug, y1=lam[,3], col = linecol, lwd=lwd)
#points(y = lam[,1], x=l.aug, pch = 19, col = viridis(length(lam[,1])), cex=pt.sz)
points(y = lam[,1], x=l.aug, pch = 19, col = "black", cex=pt.sz)
text(x = log(0.001), y = 2.8, paste("r = ",correl.est[8,1], " (",correl.est[8,2],", ", correl.est[8,3], ")", sep=""),
     pos = 4, font = 3, cex = txt2)
text(x = log(0.001), y = 2.3, paste("P(r>0) = ", P[8], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

plot(y = lam[,1], fit[,"phiO_fem_wild",1], type = "n", 
     xlim = c(0.0, 0.8), ylim = c(0, 3), yaxt="n", xaxt="n", 
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=F)
axis(1, at=c(0, 0.4, 0.8))
mtext(text="First-year survival\nwild female", side=1, line=2.5, cex=txt)
segments(x0=fit[,"phiO_fem_wild",2], y0=lam[,1], x1=fit[,"phiO_fem_wild",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiO_fem_wild",1], y0=lam[,2], x1=fit[,"phiO_fem_wild",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiO_fem_wild",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0, y = 2.8, paste("r = ",correl.est[1,1], " (",correl.est[1,2],", ", correl.est[1,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0, y = 2.3, paste("P(r>0) = ", P[1], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

par(mar = mar)
plot(y = lam[,1], fit[,"phiO_male_wild",1], type = "n", 
     xlim = c(0.0, 0.5), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=c("", "", "", ""))
axis(1, at=c(0, 0.25, 0.5), labels=c("0", "0.25", "0.5"))
mtext(text="First-year survival\nwild male", side=1, line=2.5, cex=txt)
segments(x0=fit[,"phiO_male_wild",2], y0=lam[,1], x1=fit[,"phiO_male_wild",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiO_male_wild",1], y0=lam[,2], x1=fit[,"phiO_male_wild",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiO_male_wild",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0, y = 2.8, paste("r = ",correl.est[2,1], " (",correl.est[2,2],", ", correl.est[2,3], ")", sep=""), pos = 4, font = 3, cex = txt2)
text(x = 0, y = 2.3, paste("P(r>0) = ", P[2], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

##########
# plot.new()
# legend("topleft", legend=c(1994:2018)[c(2,7,12,17,22)] , col=viridis(length(lam[,1]))[c(2,7,12,17,22)], 
#        pch=16, xpd=NA, ncol=1)

par(mar = mar)
plot(y = lam[,1], fit[,"phiO_fem_hacked",1], type = "n", 
     xlim = c(0.0, 0.5), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3))
axis(1, at=c(0, 0.25, 0.5), labels=c("0", "0.25", "0.5"))
mtext(text="First-year survival\nhacked female", side=1, line=2.5, cex=txt)
segments(x0=fit[,"phiO_fem_hacked",2], y0=lam[,1], x1=fit[,"phiO_fem_hacked",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiO_fem_hacked",1], y0=lam[,2], x1=fit[,"phiO_fem_hacked",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiO_fem_hacked",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0, y = 2.8, paste("r = ",correl.est[3,1], " (",correl.est[3,2],", ", correl.est[3,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0, y = 2.3, paste("P(r>0) = ", P[3], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

par(mar = mar)
plot(y = lam[,1], fit[,"phiO_male_hacked",1], type = "n", 
     xlim = c(0.0, 0.5), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=c("", "", "", ""))
axis(1, at=c(0, 0.25, 0.5), labels=c("0", "0.25", "0.5"))
mtext(text="First-year survival\nhacked male", side=1, line=2.5, cex=txt)
segments(x0=fit[,"phiO_male_hacked",2], y0=lam[,1], x1=fit[,"phiO_male_hacked",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiO_male_hacked",1], y0=lam[,2], x1=fit[,"phiO_male_hacked",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiO_male_hacked",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0, y = 2.8, paste("r = ",correl.est[4,1], " (",correl.est[4,2],", ", correl.est[4,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0, y = 2.3, paste("P(r>0) = ", P[4], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

par(mar = mar)
plot(y = lam[,1], fit[,"phiF",1], type = "n", 
     xlim = c(0.6, 0.9), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=c("", "", "", ""))
axis(1, at=c(0.6, 0.75, 0.9), labels=c("0.6", "0.75", "0.9"))
mtext(text="Nonbreeder\nsurvival", side=1, line=2.5, cex=txt)
segments(x0=fit[,"phiF",2], y0=lam[,1], x1=fit[,"phiF",3], y1=lam[,1], col = linecol, lwd=lwd)
segments(x0=fit[,"phiF",1], y0=lam[,2], x1=fit[,"phiF",1], y1=lam[,3], col = linecol, lwd=lwd)
points(y = lam[,1], fit[,"phiF",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0.6, y = 2.8, paste("r = ",correl.est[5,1], " (",correl.est[5,2],", ", correl.est[5,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0.6, y = 2.3, paste("P(r>0) = ", P[5], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

par(mar = mar)
plot(y = lam[2:25,1], fit[1:24,"fec",1], type = "n", 
     xlim = c(0, 2.2), ylim = c(0, 3), yaxt="n", xaxt="n",
     ylab = "", xlab = "", frame = FALSE, pch = 19)
axis(2, at=c(0, 1.0, 2, 3), labels=c("", "", "", ""))
axis(1, at=c(0, 1.1, 2.2))
mtext(text="Fecundity", side=1, line=1.5, cex=txt)
segments(x0=fit[2:25,"fec",2], y0=lam[2:25,1], x1=fit[2:25,"fec",3], y1=lam[2:25,1], col = linecol, lwd=lwd)
segments(x0=fit[2:25,"fec",1], y0=lam[2:25,2], x1=fit[2:25,"fec",1], y1=lam[2:25,3], col = linecol, lwd=lwd)
points(y = lam[2:25], fit[2:25,"fec",1], pch = 19, col = "black", cex=pt.sz)
text(x = 0, y = 2.8, paste("r = ",correl.est[7,1], " (",correl.est[7,2],", ", correl.est[7,3], ")", sep=""), 
     pos = 4, font = 3, cex = txt2)
text(x = 0, y = 2.3, paste("P(r>0) = ", P[7], sep=""), pos = 4, font = 3, cex = txt2)
abline(h=1, lty=2, lwd=1)

mtext("Population growth rate", side=2, line=1.5, las=3, outer=T, cex=0.7)

dev.off()

##############################################
# Calculate correlations 
# between demographic rates
#############################################
n.iter <- dim(out$sims.list$l.mu.F)[[1]]
correl.h <- array(NA, dim = c(n.iter,10))
for (i in 1:n.iter){
  correl.h[i,1] <- cor(out$sims.list$eta.JSalpha[i,1,2:25], out$sims.list$eta.JSalpha[i,2,2:25])
  correl.h[i,2] <- cor(out$sims.list$eta.JSalpha[i,1,2:25], out$sims.list$eta.FSalpha[i,2:25])
  correl.h[i,3] <- cor(out$sims.list$eta.JSalpha[i,1,2:25], out$sims.list$eta.BSalpha[i,2:25])
  correl.h[i,4] <- cor(out$sims.list$eta.JSalpha[i,1,2:25], out$sims.list$F[i,2:25])
  
  correl.h[i,5] <- cor(out$sims.list$eta.JSalpha[i,2,2:25], out$sims.list$eta.FSalpha[i,2:25])
  correl.h[i,6] <- cor(out$sims.list$eta.JSalpha[i,2,2:25], out$sims.list$eta.BSalpha[i,2:25])
  correl.h[i,7] <- cor(out$sims.list$eta.JSalpha[i,2,2:25], out$sims.list$F[i,2:25])
  
  correl.h[i,8] <- cor(out$sims.list$eta.FSalpha[i,2:25], out$sims.list$eta.BSalpha[i,2:25])
  correl.h[i,9] <- cor(out$sims.list$eta.FSalpha[i,2:25], out$sims.list$F[i,2:25])
  
  correl.h[i,10] <- cor(out$sims.list$eta.BSalpha[i,2:25], out$sims.list$F[i,2:25])
}

# Credible intervals of correlation coefficients
correl.est <- array(NA, dim=c(10,3))
for (k in 1:10){
  correl.est[k,2:3] <- HDIofMCMC(correl.h[,k], 0.85)
  m <- density(correl.h[,k], na.rm = TRUE)
  correl.est[k,1]<- m$x[which(m$y==max(m$y))]
}
correl.est <- round(correl.est,2)
colnames(correl.est) <- c("mode", "LHDI_85", "UHDI_85")
# Compute the posterior modes of correlation coefficients
correl.est2 <- array(NA, dim=c(10,2))
for (k in 1:10){
  correl.est2[k,1] <- mean(correl.h[,k])
  correl.est2[k,2] <- sd(correl.h[,k])
}
correl.est2 <- round(correl.est2,2)
colnames(correl.est2) <- c("mean", "sd")

# Probability that correlation coefficients (r) > 0
n.iter <- dim(out$sims.list$l.mu.F)[[1]]
P <- c()
for (k in 1:10){
  P[k] <- sum(correl.h[!is.na(correl.h[,k]),k]>0)/n.iter
}
P <- round(P,2)

corest<- cbind(correl.est, correl.est2)
rownames(corest) <- c("surv.fy.f surv.j.m", "surv.fy.f surv.nb", "surv.fy.f surv.b", "surv.fy.f F",
                      "surv.fy.m surv.nb", "surv.fy.m surv.b", "surv.fy.m F",
                      "surv.nb surv.b", "surv.nb F", "surv.b F")


######################################
# TABLES
# Overall estimates
######################################
data_summary <- function(x) {
  med <- median(x )
  mn <- mean(x) 
  sd <- sd(x)
  var <- var(x)
  return(c(med,mn,sd,var))
}

data_summary2 <- function(x) {
  med <- median(x )
  mn <- mean(x) 
  sd <- sd(x)
  var <- var(matrix(x))
  return(c(med,mn,sd,var))
}

tab <- data.frame(matrix(NA, nrow=46, ncol=6))
colnames(tab) <- c("parameter", "state sex hacked", "median", "mean", "sd", "var")

# Survival
tab[1,]<- c("Survival", "juvenile female", data_summary(out$sims.list$JSalpha1[,1]))
tab[2,]<- c("Survival", "juvenile male", data_summary(out$sims.list$JSalpha1[,2]))
tab[3,]<- c("Survival",  "nonbreeder", data_summary(out$sims.list$FSalpha1))
tab[4,]<- c("Survival",  "breeder", data_summary(out$sims.list$BSalpha1))
tab[5,]<- c("Transition", "juvenile female", data_summary(out$sims.list$JBRalpha1[,1]))
tab[6,]<- c("Transition", "juvenile male", data_summary(out$sims.list$JBRalpha1[,2]))
tab[7,]<- c("Transition", "nonbreeder female wild", data_summary(out$sims.list$FBRalpha1[,1,1]))
tab[8,]<- c("Transition", "nonbreeder male wild", data_summary(out$sims.list$FBRalpha1[,2,1]))
tab[9,]<- c("Transition", "nonbreeder female hacked", data_summary(out$sims.list$FBRalpha1[,1,2]))
tab[10,]<- c("Transition", "nonbreeder male hacked", data_summary(out$sims.list$FBRalpha1[,2,2]))
tab[11,]<- c("Transition", "breeder", data_summary(out$sims.list$BFRalpha1))
tab[12,]<- c("Resight", "nonbreeder wild, low effort", data_summary(out$sims.list$mu.pF1[,1,1]))
tab[13,]<- c("Resight", "nonbreeder hacked, low effort", data_summary(out$sims.list$mu.pF1[,1,2]))
tab[14,]<- c("Resight", "nonbreeder wild, high effort", data_summary(out$sims.list$mu.pF1[,2,1]))
tab[15,]<- c("Resight", "nonbreeder hacked, high effort", data_summary(out$sims.list$mu.pF1[,2,2]))
tab[16,]<- c("Resight", "breeder, low effort", data_summary(out$sims.list$mu.pB1[,1]))
tab[17,]<- c("Resight", "breeder, high effort", data_summary(out$sims.list$mu.pB1[,2]))
# Immigration and Fecundity
tab[18,]<- c("Fecundity <2012", "breeder", data_summary(exp(out$sims.list$l.mu.F[,1])))
tab[19,]<- c("Fecundity >=2012", "breeder", data_summary(exp(out$sims.list$l.mu.F[,2])))
#tab[20,]<- c("Fecundity all years", "breeder", data_summary(exp(out$sims.list$l.mu.F)))
tab[21,]<- c("Immigration", "nonbreeder", data_summary(exp(out$sims.list$l.omegaF)))
tab[22,]<- c("Immigration", "breeder", data_summary(exp(out$sims.list$l.omegaB)))
# Population numbers
tab[23,]<- c("Number", "juveniles (includes hacked)", data_summary2(out$sims.list$NJ))
tab[24,]<- c("Number", "juveniles (excludes hacked)", data_summary2(out$sims.list$N[,8,]))
tab[25,]<- c("Number", "nonbreeder", data_summary2(out$sims.list$NF))
tab[26,]<- c("Number", "breeder", data_summary2(out$sims.list$NB))
tab[27,]<- c("Number", "total", data_summary2(out$sims.list$Ntot))
# Temporal variance
tab[28,]<- c("Survival temporal SD", "juvenile female", data_summary(plogis(out$sims.list$sigma.JS.phi[,1])))
tab[29,]<- c("Survival temporal SD", "juvenile male", data_summary(plogis(out$sims.list$sigma.JS.phi[,2])))
tab[30,]<- c("Survival temporal SD",  "nonbreeder", data_summary(plogis(out$sims.list$sigma.FS.phi)))
tab[31,]<- c("Survival temporal SD",  "breeder", data_summary(plogis(out$sims.list$sigma.BS.phi)))
tab[32,]<- c("Transition temporal SD", "juvenile female", data_summary(plogis(out$sims.list$sigma.JBR.psi[,1])))
tab[33,]<- c("Transition temporal SD", "juvenile male", data_summary(plogis(out$sims.list$sigma.JBR.psi[,2])))
tab[34,]<- c("Transition temporal SD", "nonbreeder female wild", data_summary(plogis(out$sims.list$sigma.FBR.psi[,1,1])))
tab[35,]<- c("Transition temporal SD", "nonbreeder male wild", data_summary(plogis(out$sims.list$sigma.FBR.psi[,2,1])))
tab[36,]<- c("Transition temporal SD", "nonbreeder female hacked", data_summary(plogis(out$sims.list$sigma.FBR.psi[,1,2])))
tab[37,]<- c("Transition temporal SD", "nonbreeder male hacked", data_summary(plogis(out$sims.list$sigma.FBR.psi[,2,2])))
tab[38,]<- c("Transition temporal SD", "breeder", data_summary(plogis(out$sims.list$sigma.BFR.psi)))
tab[39,]<- c("Resight temporal SD", "nonbreeder wild", data_summary(plogis(out$sims.list$sigma.pF[,1])))
tab[40,]<- c("Resight temporal SD", "nonbreeder hacked", data_summary(plogis(out$sims.list$sigma.pF[,2])))
tab[41,]<- c("Resight temporal SD", "breeder", data_summary(plogis(out$sims.list$sigma.pB)))
# Immigration and Fecundity
tab[42,]<- c("Fecundity <2012 temporal SD", "breeder", data_summary(exp(out$sims.list$sigma.F[,1])))
tab[43,]<- c("Fecundity >=2012 temporal SD", "breeder", data_summary(exp(out$sims.list$sigma.F[,2])))
tab[44,]<- c("Fecundity territory SD", "breeder", data_summary(exp(out$sims.list$sigma.prod)))
tab[45,]<- c("Immigration temporal SD", "nonbreeder", data_summary(exp(out$sims.list$sigma.omegaF)))
tab[46,]<- c("Immigration temporal SD", "breeder", data_summary(exp(out$sims.list$sigma.omegaB)))

tab[,3] <- round(as.numeric(tab[,3]),3)
tab[,4] <- round(as.numeric(tab[,4]),3)
tab[,5] <- round(as.numeric(tab[,5]),3)
tab[,6] <- round(as.numeric(tab[,6]),3)

write.csv(tab, "Tables//ForPVA_backtransformed.csv")


tab2 <- data.frame(matrix(NA, nrow=48, ncol=6))
colnames(tab2) <- c("parameter", "state sex hacked", "median", "mean", "sd", "var")
# Survival
tab2[1,]<- c("Survival_logit", "juvenile female", data_summary(out$sims.list$JSalpha[,1]))
tab2[2,]<- c("Survival_logit", "juvenile male", data_summary(out$sims.list$JSalpha[,2]))
tab2[3,]<- c("Survival_logit",  "nonbreeder", data_summary(out$sims.list$FSalpha))
tab2[4,]<- c("Survival_logit",  "breeder", data_summary(out$sims.list$BSalpha))
tab2[5,]<- c("Transition_logit", "juvenile female", data_summary(out$sims.list$JBRalpha[,1]))
tab2[6,]<- c("Transition_logit", "juvenile male", data_summary(out$sims.list$JBRalpha[,2]))
tab2[7,]<- c("Transition_logit", "nonbreeder female wild", data_summary(out$sims.list$FBRalpha[,1,1]))
tab2[8,]<- c("Transition_logit", "nonbreeder male wild", data_summary(out$sims.list$FBRalpha[,2,1]))
tab2[9,]<- c("Transition_logit", "nonbreeder female hacked", data_summary(out$sims.list$FBRalpha[,1,2]))
tab2[10,]<- c("Transition_logit", "nonbreeder male hacked", data_summary(out$sims.list$FBRalpha[,2,2]))
tab2[11,]<- c("Transition_logit", "breeder", data_summary(out$sims.list$BFRalpha))
tab2[12,]<- c("Resight_logit", "nonbreeder wild, low effort", data_summary(out$sims.list$mu.pF[,1,1]))
tab2[13,]<- c("Resight_logit", "nonbreeder hacked, low effort", data_summary(out$sims.list$mu.pF[,1,2]))
tab2[14,]<- c("Resight_logit", "nonbreeder wild, high effort", data_summary(out$sims.list$mu.pF[,2,1]))
tab2[15,]<- c("Resight_logit", "nonbreeder hacked, high effort", data_summary(out$sims.list$mu.pF[,2,2]))
tab2[16,]<- c("Resight_logit", "breeder, low effort", data_summary(out$sims.list$mu.pB[,1]))
tab2[17,]<- c("Resight_logit", "breeder, high effort", data_summary(out$sims.list$mu.pB[,2]))
# Immigration and Fecundity
tab2[18,]<- c("Fecundity <2012_log", "breeder", data_summary((out$sims.list$l.mu.F[,1])))
tab2[19,]<- c("Fecundity >=2012_log", "breeder", data_summary((out$sims.list$l.mu.F[,2])))
tab2[20,]<- c("Immigration_log", "nonbreeder", data_summary((out$sims.list$l.omegaF)))
tab2[21,]<- c("Immigration_log", "breeder", data_summary((out$sims.list$l.omegaB)))
# Population numbers
tab2[22,]<- c("Abundance_log", "juveniles (includes hacked)", data_summary2(log(out$sims.list$NJ)))
tab2[23,]<- c("Abundance_log", "juveniles (excludes hacked)", data_summary2(ifelse(out$sims.list$N[,8,]==0, log(0.001),log(out$sims.list$N[,8,]))))
tab2[24,]<- c("Abundance_log", "nonbreeder", data_summary2(ifelse(out$sims.list$NF==0, log(0.001),log(out$sims.list$NF))))
tab2[25,]<- c("Abundance_log", "breeder", data_summary2(ifelse(out$sims.list$NB==0, log(0.001),log(out$sims.list$NB))))
tab2[26,]<- c("Abundance_log", "total", data_summary2(log(out$sims.list$Ntot)))
tab2[27,]<- c("Abundance_Obs Err_SD_log", "juveniles (includes hacked)", data_summary(out$sims.list$sigma.JM))
tab2[28,]<- c("Abundance_Obs Err_SD_log", "nonbreeder", data_summary(out$sims.list$sigma.FM))
tab2[29,]<- c("Abundance_Obs Err_SD_log", "breeder", data_summary(out$sims.list$sigma.BM))
# Temporal variance
tab2[30,]<- c("Survival temporal SD_logit", "juvenile female", data_summary((out$sims.list$sigma.JS.phi[,1])))
tab2[31,]<- c("Survival temporal SD_logit", "juvenile male", data_summary((out$sims.list$sigma.JS.phi[,2])))
tab2[32,]<- c("Survival temporal SD_logit",  "nonbreeder", data_summary((out$sims.list$sigma.FS.phi)))
tab2[33,]<- c("Survival temporal SD_logit",  "breeder", data_summary((out$sims.list$sigma.BS.phi)))
tab2[34,]<- c("Transition temporal SD_logit", "juvenile female", data_summary((out$sims.list$sigma.JBR.psi[,1])))
tab2[35,]<- c("Transition temporal SD_logit", "juvenile male", data_summary((out$sims.list$sigma.JBR.psi[,2])))
tab2[36,]<- c("Transition temporal SD_logit", "nonbreeder female wild", data_summary(plogis(out$sims.list$sigma.FBR.psi[,1,1])))
tab2[37,]<- c("Transition temporal SD_logit", "nonbreeder male wild", data_summary((out$sims.list$sigma.FBR.psi[,2,1])))
tab2[38,]<- c("Transition temporal SD_logit", "nonbreeder female hacked", data_summary((out$sims.list$sigma.FBR.psi[,1,2])))
tab2[39,]<- c("Transition temporal SD_logit", "nonbreeder male hacked", data_summary((out$sims.list$sigma.FBR.psi[,2,2])))
tab2[40,]<- c("Transition temporal SD_logit", "breeder", data_summary((out$sims.list$sigma.BFR.psi)))
tab2[41,]<- c("Resight temporal SD_logit", "nonbreeder wild", data_summary((out$sims.list$sigma.pF[,1])))
tab2[42,]<- c("Resight temporal SD_logit", "nonbreeder hacked", data_summary((out$sims.list$sigma.pF[,2])))
tab2[43,]<- c("Resight temporal SD_logit", "breeder", data_summary((out$sims.list$sigma.pB)))
# Immigration and Fecundity
tab2[44,]<- c("Fecundity <2012 temporal SD_log", "breeder", data_summary((out$sims.list$sigma.F[,1])))
tab2[45,]<- c("Fecundity >=2012 temporal SD_log", "breeder", data_summary((out$sims.list$sigma.F[,2])))
tab2[46,]<- c("Fecundity territory_log", "breeder", data_summary((out$sims.list$sigma.prod)))
tab2[47,]<- c("Immigration temporal SD_log", "nonbreeder", data_summary((out$sims.list$sigma.omegaF)))
tab2[48,]<- c("Immigration temporal SD_log", "breeder", data_summary((out$sims.list$sigma.omegaB)))

tab2[,3] <- round(as.numeric(tab2[,3]),3)
tab2[,4] <- round(as.numeric(tab2[,4]),3)
tab2[,5] <- round(as.numeric(tab2[,5]),3)
tab2[,6] <- round(as.numeric(tab2[,6]),3)
write.csv(tab2, "Tables//ForPVA_untransformed.csv")

## ---- data summaries --------
table(datl$y)
prod.yr <- tapply(datl$prod, list(datl$year.p), sum, na.rm=T)
mean(prod.yr, na.rm=T)
sd(prod.yr, na.rm=T)
range(prod.yr, na.rm=T)
mean(datl$prod, na.rm=T)
sd(datl$prod, na.rm=T)

## ---- immigration ---------
load("outputs\\ipm-ie.Rdata")
source("R/HDIofMCMC.R")
options(scipen = 100)
x <- seq(from=0, to=4, by=0.01)
df <- data.frame(x=x, prior=dnorm(x,0,1))
df2 <- data.frame(omega=out$sims.list$omega)

par(mar=c(5,5,2,5))
plot(df$x, df$prior, type="l", lwd=2, xlim=c(0,4),
     ylab="Prior density (black line)", 
     xlab=expression(paste("Immigration rate (",omega,")", sep="")),
     main="")
par(new=TRUE)
plot(density(df2$omega), type="l", lwd=2, xlim=c(0,4),
     main="", col="blue",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = c(0,30,60))      # Add second axis
mtext("Posterior density (blue line)", side = 4, line = 3)  

# model estimated number of immigrants and imm rate
mnI <- apply(out$sims.list$NI[,3:26], 1, mean)
median(mnI)
HDIofMCMC(mnI)
median(out$sims.list$omega)
HDIofMCMC(out$sims.list$omega)

hist(out$sims.list$omega, main="Histogram of Immigration", xlab="Posterior draws of immigration")
abline(v=median(out$sims.list$omega), lwd=2)
abline(v=HDIofMCMC(out$sims.list$omega, cred=0.95), lty=2, lwd=2)
median(out$sims.list$omega)
HDIofMCMC(out$sims.list$omega, cred=0.95)