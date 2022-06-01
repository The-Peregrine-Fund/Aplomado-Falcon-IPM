## ---- ipm2 --------
#################################
# The model
################################
library (jagsUI)
load(".\\data\\data-7states.Rdata")
datl$pp <- c(0, 0, tapply(datl$prod, datl$year.p, sum, na.rm=T), 42)
datl$pp[c(15,17)] <- NA
m<- c("ipm-jags-no-imm")
modfl <- paste(".\\", m, ".txt", sep="")
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
    #     4 not seen
    #   States (ps)
    #     1 alive first year
    #     2 alive nonbreeder
    #     3 alive breeder
    #     4 dead
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
    #       (this is the letter 'O' rather than zero so jags can parse)
    #   sA: survival probability nonbreeders
    #   sB: survival probability breeders
    #   psiOB: recruitment probability from first-year to breeder
    #   psiAB: recruitment probability from nonbreeders to breeder
    #   psiBA: recruitment probability from breeder to nonbreeders 
    #   pO: resight probability first-year
    #   pA: resight probability nonbreeders
    #   pB: resight probability breeder
    #   F: fecundity
    #   omega: immigration
    ###################################################
    # Priors and constraints
    ###################################################
    # Counts for first year are known with certainty
    N[1,1] <- 0 # Wild born first-year
    N[2,1] <- 0 # First-year to nonbreeder
    N[3,1] <- 0  # First-year to breeders
    N[4,1] <- 0  # Nonbreeder wild to nonbreeder
    N[5,1] <- 0  # Nonbreeder wild to breeders
    N[6,1] <- 0  # Breeders to nonbreeder
    N[7,1] <- 0  # Breeders to breeders
    N[8,1] <- 0  # Wild first-year
    N[9,1] <- 0  # Nonbreeder hacked to nonbreeder
    # N[10,1] <- 0  # Immigrant to breeder
    # N[11,1] <- 0  # Immigrant to nonbreeder
    F[1] <- 0 # No breeding during the first year, all hacked birds
    
    sigma.BM ~ dnorm(0, 1/(2*2) )T(0,)
    sigma.AM ~ dnorm(0, 1/(2*2) )T(0,)
    sigma.OM ~ dnorm(0, 1/(2*2) )T(0,) 
    
    for (m in 1:2){ mu.F[m] ~ dunif(0,5) } # m # limits to help run model
    sigma.F ~ dnorm(0, 1/(2*2) )T(0,)
    
    # Survival loops for demographic categories by sex, hacked, effort 
    
    sigma.AS.phi ~ dunif(0,10)
    ASalpha <- logit(ASalpha1)
    ASalpha1 ~ dunif(0, 1)
    
    sigma.BS.phi ~ dunif(0,10)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    sigma.BAR.psi ~ dunif(0,10)
    BARalpha <- logit(BARalpha1)
    BARalpha1 ~ dunif(0, 1)
    
    sigma.pB ~ dunif(0,10)
    
    for (k in 1:2){
    mu.pB[k]<- logit(mu.pB1[k])
    mu.pB1[k] ~ dunif(0, 1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.ASalpha[t]) <- ASalpha + eps.AS.phi[t]
    eps.AS.phi[t] ~ dnorm(0, 1/(sigma.AS.phi*sigma.AS.phi) )
    
    logit(eta.BSalpha[t]) <- BSalpha + eps.BS.phi[t]
    eps.BS.phi[t] ~ dnorm(0, 1/(sigma.BS.phi*sigma.BS.phi) )    
    
    logit(eta.BARalpha[t]) <- BARalpha + eps.BAR.psi[t]
    eps.BAR.psi[t] ~ dnorm(0, 1/(sigma.BAR.psi*sigma.BAR.psi))
    
    logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
    eps.pB[t] ~ dnorm(0, 1/(sigma.pB*sigma.pB) )
    } #t
    
    for(s in 1:2){
    sigma.OS.phi[s] ~ dunif(0,10)
    OSalpha[s] <- logit(OSalpha1[s]) 
    OSalpha1[s] ~ dunif(0, 1)
    
    
    sigma.OBR.psi[s] ~ dunif(0,10)
    OBRalpha[s] <- logit(OBRalpha1[s]) 
    OBRalpha1[s] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.OSalpha[s,t]) <- OSalpha[s] + eps.OS.phi[s,t]
    eps.OS.phi[s,t] ~ dnorm(0, 1/(sigma.OS.phi[s]*sigma.OS.phi[s]) )
    
    logit(eta.OBRalpha[s,t]) <- OBRalpha[s] + eps.OBR.psi[s,t]
    eps.OBR.psi[s,t] ~ dnorm(0, 1/(sigma.OBR.psi[s]*sigma.OBR.psi[s]) )
    } #t
    
    for (h in 1:2){
    sigma.ABR.psi[s,h] ~ dunif(0,10)
    ABRalpha[s,h] <- logit(ABRalpha1[s,h])
    ABRalpha1[s,h] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.ABRalpha[s,h,t]) <- ABRalpha[s,h] + eps.ABR.psi[s,h,t]
    eps.ABR.psi[s,h,t] ~ dnorm(0, 1/(sigma.ABR.psi[s,h]*sigma.ABR.psi[s,h]) )
    } #t
    } } #h #s
    
    for (h in 1:2){
    sigma.pA[h] ~ dunif(0,10)
    for (k in 1:2){
    mu.pA[k,h] <- logit(mu.pA1[k,h])
    mu.pA1[k,h] ~ dunif(0,1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.pA[h,t]) <- mu.pA[effort[t], h] + eps.pA[h,t]
    eps.pA[h,t] ~ dnorm(0, 1/(sigma.pA[h]*sigma.pA[h]) )
    }} #t #h
    
    #######################
    # Derived params
    #######################
    for (t in 3:(n.yr-1)){    
    lambda[t] <-  Ntot[t+1]/(Ntot[t])
    loglambda[t] <- log(lambda[t])
    } #t
    
    ###############################
    # Likelihood for productivity
    ###############################
    for (t in 2:(n.yr-1)){ 
      pp[t] ~ dpois( F[t]*NB[t] )
      log(F[t]) <- log(mu.F[manage[t]]) + eps.F[t]
      eps.F[t] ~ dnorm (0, 1/(sigma.F*sigma.F) )
    }
    # GOF
    for (t in 1:(n.yr-1)){ 
      J.rep[t] ~ dpois( F[t]*NB[t] )
      J.exp[t] <- F[t]*NB[t]
      d.obs[t] <- pp[t]* log((pp[t]+0.001)/(J.exp[t]+0.001)) - (pp[t]-J.exp[t])
      d.rep[t] <- J.rep[t]*log((J.rep[t]+0.001)/(J.exp[t]+0.001)) - (J.rep[t]-J.exp[t])
    } # t
    dd.obs <- sum(d.obs)
    tvm.obs <- sd(pp)^2/mean(pp)
    dd.rep <- sum(d.rep)
    tvm.rep <- sd(J.rep)^2/mean(J.rep)
    
    ################################
    # Likelihood for counts
    ################################
   for (t in 1:(n.yr-1)){ 
    # Number of wild born juvs
    N[1,t+1] ~ dpois( (NO[t]*eta.OSalpha[2,t]*eta.OBRalpha[2,t] + # first year males to breeders
    NF[t]*eta.ASalpha[t]*eta.ABRalpha[2,2,t] + # nonbreeder male hacked adults to breeder
    NF[t]*eta.ASalpha[t]*eta.ABRalpha[2,1,t] + # nonbreeder male wild adults to breeder adults
    NB[t]*eta.BSalpha[t]*(1-eta.BARalpha[t]) )* # breeders to breeders
    F[t+1]/2) 
    # first year to nonbreeder adults
    N[2,t+1] ~ dbin(eta.OSalpha[2,t]*(1-eta.OBRalpha[2,t]), NO[t] )
    # first year to breeders
    N[3,t+1] ~ dbin(eta.OSalpha[2,t]*eta.OBRalpha[2,t], NO[t] )
    # Nonbreeder male wild adults to nonbreeder adults
    N[4,t+1] ~ dbin(eta.ASalpha[t]*(1-eta.ABRalpha[2,1,t]), NF[t] )
    # Nonbreeder male wild adults to breeders
    N[5,t+1] ~ dbin(eta.ASalpha[t]*eta.ABRalpha[2,1,t], NF[t] )
    # Breeders to nonbreeder adults
    N[6,t+1] ~ dbin(eta.BSalpha[t]*(eta.BARalpha[t]), NB[t] )
    # Breeders to breeders
    N[7,t+1] ~ dbin(eta.BSalpha[t]*(1-eta.BARalpha[t]), NB[t] )    
    # Nonbreeder hacked adults to nonbreeder adults
    N[8,t+1] ~ dbin(eta.ASalpha[t]*(1-eta.ABRalpha[2,2,t]), NF[t] )
    # Nonbreeder hacked adults to breeders
    N[9,t+1] ~ dbin(eta.ASalpha[t]*eta.ABRalpha[2,2,t], NF[t] )
    } # t
    
    for (t in 1:n.yr){
    Ntot[t] <- sum(N[c(1,2,3,4,5,6,7,8,9),t]) + aug[t] # total number
    NB[t] <-  sum(N[c(3,5,7,9),t]) # number of breeders
    NF[t] <-  sum(N[c(2,4,6,8),t]) # number of nonbreeders
    NO[t] <-  N[1,t] + aug[t] # number of first years
    } # t
    
    # Observation process    
    for (t in 2:n.yr){
    countBM[t] ~ dpois(NB[t]) # breeding males
    countFM[t] ~ dpois(NF[t]) # nonbreeding adult males
    countJM[t] ~ dpois(NO[t]) # first year males
    } # t
    
    ###################
    # Assess fit of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 1:n.yr){ 
    c.expB[t] <- NB[t] + 0.001 # expected counts adult breeder
    c.expA[t] <- NF[t] + 0.001 # nonbreeder
    c.expO[t] <- NO[t] + 0.001 # first year
    c.obsB[t] <- countBM[t] + 0.001
    c.obsA[t] <- countFM[t] + 0.001
    c.obsO[t] <- countJM[t] + 0.001
    dssm.obsB[t] <- abs( ( (c.obsB[t]) - (c.expB[t]) ) / (c.obsB[t]+0.001)  )
    dssm.obsA[t] <- abs( ( (c.obsA[t]) - (c.expA[t]) ) / (c.obsA[t]+0.001)  )
    dssm.obsO[t] <- abs( ( (c.obsO[t]) - (c.expO[t]) ) / (c.obsO[t]+0.001)  )
    } # t
    dmape.obs[1] <- sum(dssm.obsB[2:n.yr])
    dmape.obs[2] <- sum(dssm.obsA[2:n.yr])
    dmape.obs[3] <- sum(dssm.obsO[2:n.yr])
    # Compute fit statistic for replicate data
    # Mean absolute error
    for (t in 1:n.yr){ 
    c.repB[t] ~ dpois(NB[t] ) # expected counts
    c.repA[t] ~ dpois(NF[t] ) 
    c.repO[t] ~ dpois(NO[t] ) 
    } # t
    for (t in 1:n.yr){ 
    dssm.repB[t] <- abs( ( (c.repB[t]) - (c.expB[t]) ) / (c.repB[t]+0.001) )
    dssm.repA[t] <- abs( ( (c.repA[t]) - (c.expA[t]) ) / (c.repA[t]+0.001) )
    dssm.repO[t] <- abs( ( (c.repO[t]) - (c.expO[t]) ) / (c.repO[t]+0.001) )
    } # t
    dmape.rep[1] <- sum(dssm.repB[2:n.yr])
    dmape.rep[2] <- sum(dssm.repA[2:n.yr])
    dmape.rep[3] <- sum(dssm.repO[2:n.yr])
    
    # Test statistic for number of turns
    for (t in 1:(n.yr-2)){
    tt1.obsB[t] <- step(countBM[t+2] - countBM[t+1])
    tt2.obsB[t] <- step(countBM[t+1] - countBM[t])
    tt3.obsB[t] <- equals(tt1.obsB[t] + tt2.obsB[t], 1)
    tt1.obsA[t] <- step(countFM[t+2] - countFM[t+1])
    tt2.obsA[t] <- step(countFM[t+1] - countFM[t])
    tt3.obsA[t] <- equals(tt1.obsA[t] + tt2.obsA[t], 1)
    tt1.obsO[t] <- step(countJM[t+2] - countJM[t+1])
    tt2.obsO[t] <- step(countJM[t+1] - countJM[t])
    tt3.obsO[t] <- equals(tt1.obsO[t] + tt2.obsO[t], 1)
    } # t
    tturn.obs[1] <- sum(tt3.obsB)
    tturn.obs[2] <- sum(tt3.obsA)
    tturn.obs[3] <- sum(tt3.obsO)
    
    for (t in 1:(n.yr-2)){
    tt1.repB[t] <- step(c.repB[t+2] - c.repB[t+1])
    tt2.repB[t] <- step(c.repB[t+1] - c.repB[t])
    tt3.repB[t] <- equals(tt1.repB[t] + tt2.repB[t], 1)
    tt1.repA[t] <- step(c.repA[t+2] - c.repA[t+1])
    tt2.repA[t] <- step(c.repA[t+1] - c.repA[t])
    tt3.repA[t] <- equals(tt1.repA[t] + tt2.repA[t], 1)
    tt1.repO[t] <- step(c.repO[t+2] - c.repO[t+1])
    tt2.repO[t] <- step(c.repO[t+1] - c.repO[t])
    tt3.repO[t] <- equals(tt1.repO[t] + tt2.repO[t], 1)
    } # t
    tturn.rep[1] <- sum(tt3.repB)
    tturn.rep[2] <- sum(tt3.repA)
    tturn.rep[3] <- sum(tt3.repO)
    
    ################################
    # Likelihood for survival
    ################################
    for (i in 1:nind){
    for (t in 1:(n.yr-1)){
    #Survival
    phiO[i,t] <- eta.OSalpha[sex[i],t] # first year
    phiA[i,t] <- eta.ASalpha[t] # nonbreeder
    phiB[i,t] <- eta.BSalpha[t] # breeder
    #Recruitment
    psiOB[i,t] <- eta.OBRalpha[sex[i],t] # first year to breeder
    psiAB[i,t] <- eta.ABRalpha[sex[i], hacked[i],t] # nonbreederto breeder
    psiBA[i,t] <- eta.BARalpha[t] # breeder to nonbreeder
    #Re-encounter
    pA[i,t] <- eta.pA[hacked[i],t] # resight of nonbreeders
    pB[i,t] <- eta.pB[t]  # resight of breeders
    }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in first[i]:(n.yr-1)){
    ps[1,i,t,1] <- 0
    ps[1,i,t,2] <- phiO[i,t] * (1-psiOB[i,t])
    ps[1,i,t,3] <- phiO[i,t] * psiOB[i,t]
    ps[1,i,t,4] <- 1-phiO[i,t]
    
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phiA[i,t] * (1-psiAB[i,t])
    ps[2,i,t,3] <- phiA[i,t] * psiAB[i,t]
    ps[2,i,t,4] <- 1-phiA[i,t]
    
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- phiB[i,t] * psiBA[i,t]
    ps[3,i,t,3] <- phiB[i,t] * (1-psiBA[i,t])
    ps[3,i,t,4] <- 1-phiB[i,t]
    
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 1 
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pA[i,t]
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 1-pA[i,t]
    
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- pB[i,t]
    po[3,i,t,4] <- 1-pB[i,t]
    
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]]
    for (t in (first[i]+1):n.yr){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:4])
    } #t
    } #i
    
    } #model
    ",fill = TRUE)
sink()

datl$y[datl$y==5] <- 4

get.first <- function(x) min(which(x!=4))
f <- apply(datl$y, 1, get.first)

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 2:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

inits <- function(){list(z = ms.init.z(datl$y, f) )}

params <- c(
  "N", "NB", "NF", "NO", "NI", "Ntot", 
  "F",  "sigma.F", "eps.F", "mu.F",   "omega",
  "dmape.obs", "dmape.rep", "tvm.obs", "tvm.rep", "dd.obs", "dd.rep",
  "tturn.obs", "tturn.rep",
  "OSalpha", "ASalpha", "BSalpha", "OBRalpha", "ABRalpha", "BARalpha",  "mu.pA", "mu.pB",
  "OSalpha1", "ASalpha1", "BSalpha1","OBRalpha1", "ABRalpha1", "BARalpha1","mu.pA1", "mu.pB1",
  "eps.OS.s", "eps.AS.s", "eps.BS.s", "eps.OBR.psi", "eps.ABR.psi", "eps.BAR.psi", "eps.pA", "eps.pB", 
  "eta.OSalpha", "eta.ASalpha", "eta.BSalpha", "eta.OBRalpha", "eta.ABRalpha", "eta.BARalpha", "eta.pA",  "eta.pB", 
  "sigma.OS.s", "sigma.AS.s", "sigma.BS.s", "sigma.OBR.psi", "sigma.ABR.psi", "sigma.BAR.psi", "sigma.pA", "sigma.pB"  
)

datl$z <- known.state.ms(datl$y, 4)

# MCMC settings
ni <- 200000; nt <- 50; nb <- 150000; nc <- 3; na <- 1000 
out <- jags(datl, inits, params, modfl,  
            n.chains = nc, n.thin = nt, n.burnin = nb, n.adapt=na, n.iter=ni, 
            parallel=T, module=c("glm", "bugs"))
#save(file=paste("./", m, ".Rdata", sep=""), list="out")
