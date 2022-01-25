## ---- ipm2 --------
#################################
# The model
################################
library (jagsUI)
load(".\\data\\data-7states.Rdata")
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
    sigma.prod ~ dunif(0,10)
    
    for (m in 1:2){ mu.F[m] ~ dunif(0,5) } # m # limits to help run model
    sigma.F ~ dunif(0,10)
    
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
    for (k in 1:K){ prod[k] ~ dpois( prod.mu[year.p[k] , ter[k]] ) }
    for (j in 1:n.ter){
    for (t in 3:(n.yr-1)){ 
    log(prod.mu[t,j]) <- log(F[t]) + eps.prod[j] # link productivity to fecundity
    } # t
    eps.prod[j] ~ dnorm(0, 1/(sigma.prod*sigma.prod))  
    } # j
    
    for (t in 1:(n.yr-1)){ 
    log(F[t+1]) <- log(mu.F[manage[t+1]]) + eps.F[t+1]
    eps.F[t+1] ~ dnorm (0, 1/(sigma.F*sigma.F) )
    
    ################################
    # Likelihood for counts
    ################################
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
    l.countBM[t] ~ dnorm(log(NB[t]), 1/(sigma.BM*sigma.BM)) # breeding males
    l.countFM[t] ~ dnorm(log(NF[t]), 1/(sigma.AM*sigma.AM)) # nonbreeding adult males
    l.countJM[t] ~ dnorm(log(NO[t]), 1/(sigma.OM*sigma.OM)) # first year males
    } # t
    
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
  "N", "NB", "NF", "NO", "NI", "Ntot", "sigma.BM", "sigma.AM", "sigma.OM",
  "F", "sigma.prod", "sigma.F", "eps.prod", "eps.F", "mu.F",
  #  "omegaA", "omegaB", "omegaA1", "omegaB1", "sigma.omegaA", "sigma.omegaB", "eps.omegaA", "eps.omegaB",
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
