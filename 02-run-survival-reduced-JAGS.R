## ---- reduced --------
# JAGS must be installed on your computer
library (jagsUI)
load("./data/data-7states.Rdata")
m<- c("surv7-jags-reduced")
modfl <- paste("./", m, ".txt", sep="")
sink(modfl)
cat("
    model{
    ####################################################
    ####################################################
    # Mark-resight-recovery data
    #   Observations (po) = y  
    #     1 seen first year (age 0)
    #     2 seen nonbreeder
    #     3 seen breeder
    #     4 recovered dead
    #     5 not seen
    #   States (ps)
    #     1 alive first year
    #     2 alive nonbreeder
    #     3 alive breeder
    #     4 dead
    #     5 dead not recovered
    #     6 emigrated alive
    #     7 emigrated dead
    #   Groups
    #     1 wild born
    #     2 hacked
    #   Sex
    #     1 female
    #     2 male
    #   Effort
    #     1 low
    #     2 high
    ###################################################
    # PARAMETERS
    #   sO: survival probability first year 
    #       (sO as in the letter O rather than zero so jags can parse)
    #   sA: survival probability nonbreeders
    #   sB: survival probability breeders
    #   psiOB: recruitment probability from first-year to breeder
    #   psiAB: recruitment probability from nonbreeders to breeder
    #   psiBA: recruitment probability from breeder to nonbreeders 
    #   pO: resight probability first-year
    #   pA: resight probability nonbreeders
    #   pB: resight probability breeder
    #   em: probability of emigration
    #   r: probability of dead recovery
    ###################################################
    # Priors and constraints
    ###################################################
    l.mu.em <- logit(mu.em)
    mu.em ~ dunif(0,1)

    sigma.AS.s ~ dunif(0,10)
    ASalpha <- logit(ASalpha1)
    ASalpha1 ~ dunif(0, 1)
    
    sigma.BS.s ~ dunif(0,10)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    sigma.BAR.psi ~ dunif(0,10)
    BARalpha <- logit(BARalpha1)
    BARalpha1 ~ dunif(0, 1)
    r ~ dunif(0, 1)
    sigma.pB ~ dunif(0,10)
    sigma.pA ~ dunif(0,10)
    
    for (k in 1:2){
    mu.pB[k]<- logit(mu.pB1[k])
    mu.pB1[k] ~ dunif(0, 1)
    mu.pA[k] <- logit(mu.pA1[k])
    mu.pA1[k] ~ dunif(0,1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.ASalpha[t]) <- ASalpha + eps.AS.s[t]
    eps.AS.s[t] ~ dnorm(0, 1/(sigma.AS.s*sigma.AS.s))
    logit(eta.BSalpha[t]) <- BSalpha + eps.BS.s[t]
    eps.BS.s[t] ~ dnorm(0, 1/(sigma.BS.s*sigma.BS.s))    
    logit(eta.BARalpha[t]) <- BARalpha + eps.BAR.psi[t]
    eps.BAR.psi[t] ~ dnorm(0, 1/(sigma.BAR.psi*sigma.BAR.psi))
    logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
    eps.pB[t] ~ dnorm(0, 1/(sigma.pB*sigma.pB))
    } #t
    
    for(s in 1:2){
    sigma.OS.s[s] ~ dunif(0,10)
    OSalpha[s] <- logit(OSalpha1[s]) 
    OSalpha1[s] ~ dunif(0, 1)
    sigma.OBR.psi[s] ~ dunif(0,10)
    OBRalpha[s] <- logit(OBRalpha1[s]) 
    OBRalpha1[s] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.OSalpha[s,t]) <- OSalpha[s] + eps.OS.s[s,t]
    eps.OS.s[s,t] ~ dnorm(0, 1/(sigma.OS.s[s]*sigma.OS.s[s]))
    logit(eta.OBRalpha[s,t]) <- OBRalpha[s] + eps.OBR.psi[s,t]
    eps.OBR.psi[s,t] ~ dnorm(0, 1/(sigma.OBR.psi[s]*sigma.OBR.psi[s]))
    } #t
    
    for (h in 1:2){
    sigma.ABR.psi[s,h] ~ dunif(0,10)
    ABRalpha[s,h] <- logit(ABRalpha1[s,h])
    ABRalpha1[s,h] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.ABRalpha[s,h,t]) <- ABRalpha[s,h] + eps.ABR.psi[s,h,t]
    eps.ABR.psi[s,h,t] ~ dnorm(0, 1/(sigma.ABR.psi[s,h]*sigma.ABR.psi[s,h]))
    } } }   #t #h #s
    
    for (t in 1:(n.yr-1)){
    logit(eta.pA[t]) <- mu.pA[effort[t]] + eps.pA[t]
    eps.pA[t] ~ dnorm(0, 1/(sigma.pA*sigma.pA)) 
    logit(em[t]) <-  l.mu.em # eps.em[t]    
} #t 
    
    ################################
    # Likelihood for survival
    ################################
    for (i in 1:nind){
    for (t in 1:(n.yr-1)){
    #Survival
    sO[i,t] <- eta.OSalpha[sex[i],t] # first year
    sA[i,t] <- eta.ASalpha[t] # nonbreeder
    sB[i,t] <- eta.BSalpha[t] # breeder
    #Recruitment
    psiOB[i,t] <- eta.OBRalpha[sex[i],t] # first year to breeder
    psiAB[i,t] <- eta.ABRalpha[sex[i], hacked[i],t] # nonbreeder to breeder
    psiBA[i,t] <- eta.BARalpha[t] # breeder to nonbreeder
    #Re-encounter
    pA[i,t] <- eta.pA[t] # resight of nonbreeders
    pB[i,t] <- eta.pB[t]  # resight of breeders
    }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in first[i]:(n.yr-1)){
    ps[1,i,t,1] <- 0
    ps[1,i,t,2] <- sO[i,t] * (1-psiOB[i,t]) * (1-em[t])
    ps[1,i,t,3] <- sO[i,t] * psiOB[i,t] * (1-em[t])
    ps[1,i,t,4] <- (1-sO[i,t]) * r * (1-em[t])
    ps[1,i,t,5] <- (1-sO[i,t]) * (1-r) * (1-em[t])
    ps[1,i,t,6] <- sO[i,t] * em[t]
    ps[1,i,t,7] <- (1-sO[i,t]) * (1-r) * em[t]
    
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- sA[i,t] * (1-psiAB[i,t]) * (1-em[t])
    ps[2,i,t,3] <- sA[i,t] * psiAB[i,t] * (1-em[t])
    ps[2,i,t,4] <- (1-sA[i,t]) * r * (1-em[t])
    ps[2,i,t,5] <- (1-sA[i,t]) * (1-r) * (1-em[t])
    ps[2,i,t,6] <- sA[i,t] * em[t]
    ps[2,i,t,7] <- (1-sA[i,t]) * (1-r) * em[t]
    
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- sB[i,t] * psiBA[i,t] * (1-em[t])
    ps[3,i,t,3] <- sB[i,t] * (1-psiBA[i,t]) * (1-em[t])
    ps[3,i,t,4] <- (1-sB[i,t]) * r * (1-em[t])
    ps[3,i,t,5] <- (1-sB[i,t]) * (1-r) * (1-em[t])
    ps[3,i,t,6] <- sB[i,t] * em[t]
    ps[3,i,t,7] <- (1-sB[i,t]) * (1-r) * em[t]
    
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 0
    ps[4,i,t,5] <- 1
    ps[4,i,t,6] <- 0
    ps[4,i,t,7] <- 0
    
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 1
    ps[5,i,t,6] <- 0
    ps[5,i,t,7] <- 0
    
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    ps[6,i,t,7] <- 0
    
    ps[7,i,t,1] <- 0
    ps[7,i,t,2] <- 0
    ps[7,i,t,3] <- 0
    ps[7,i,t,4] <- 0
    ps[7,i,t,5] <- 0
    ps[7,i,t,6] <- 0
    ps[7,i,t,7] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 1 
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pA[i,t]
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 1-pA[i,t]
    
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- pB[i,t]
    po[3,i,t,4] <- 0
    po[3,i,t,5] <- 1-pB[i,t]
    
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- 1
    po[4,i,t,5] <- 0
    
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1
    
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 1
    
    po[7,i,t,1] <- 0
    po[7,i,t,2] <- 0
    po[7,i,t,3] <- 0
    po[7,i,t,4] <- 0
    po[7,i,t,5] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]]
    for (t in (first[i]+1):n.yr){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:7])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:5])
    } #t
    } #i
    } #model
    ",fill = TRUE)
sink()  

# Initial values
get.first <- function(x) min(which(x!=5))
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=5))
l <- apply(datl$y, 1, get.last)
TFmat <- is.na(z.inits) & is.na(datl$z)
for (i in 1:dim(TFmat)[1]){  TFmat[i,1:f[i]] <- FALSE }
z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 

inits <- function(){list(z = z.inits)} 

params <- c(
  "r", "l.mu.em", "mu.em", "em",
  "OSalpha", "ASalpha", "BSalpha", "OBRalpha", "ABRalpha", "BARalpha",  "mu.pA", "mu.pB",
  "OSalpha1", "ASalpha1", "BSalpha1","OBRalpha1", "ABRalpha1", "BARalpha1","mu.pA1", "mu.pB1",
  "eps.OS.s", "eps.AS.s", "eps.BS.s", "eps.OBR.psi", "eps.ABR.psi", "eps.BAR.psi", "eps.pA", "eps.pB", 
  "eta.OSalpha", "eta.ASalpha", "eta.BSalpha", "eta.OBRalpha", "eta.ABRalpha", "eta.BARalpha", "eta.pA",  "eta.pB", 
  "sigma.OS.s", "sigma.AS.s", "sigma.BS.s", "sigma.OBR.psi", "sigma.ABR.psi", "sigma.BAR.psi", "sigma.pA", "sigma.pB"  
)

# MCMC settings
# ni <- 200000; nt <- 100; nb <- 100000; nc <- 3; na <- 10000 # actual run # takes about 1 week on a high performance computer
ni <- 10000; nt <- 1; nb <- 5000; nc <- 1; na <- 100 # time-saving run, but inadequate for convergence
out <- jags(datl, inits, params, modfl,  
            n.chains = nc, n.thin = nt, n.burnin = nb, n.adapt=na, n.iter=ni, 
            parallel=T, module=c("glm", "bugs"))
save(file=paste("./", m, ".Rdata", sep=""), list="out")
