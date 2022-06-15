## ---- ipm2 --------
#################################
# The model
################################
library (jagsUI)
load(".\\data\\data-7states.Rdata")
datl$pp <- c(0, 0, tapply(datl$prod, datl$year.p, sum, na.rm=T), 42)
datl$pp[c(15,17)] <- NA
datl$countOM <- datl$countJM-datl$aug

m<- c("ipm-e")
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
    #     4 seen dead
    #     5 not seen
    #   States (ps)
    #     1 alive first year
    #     2 alive nonbreeder
    #     3 alive breeder
    #     4 Recovered recently dead
    #     5 Dead not recovered/long dead
    #     6 Emigrated and alive
    #     7 Emigrated and dead
    #   Groups
    #     1 wild hatched
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
    #   pA: resight probability nonbreeders
    #   pB: resight probability breeder
    #   F: fecundity
    #   omega: immigration
    ###################################################
    # Priors and constraints
    ###################################################
    # Counts for first year are known with certainty
    N[1,1] <- 0 # Wild born first-year
    N[2,1] <- 0 # first-year hacked to nonbreeder adults
    N[3,1] <- 0 # First-year wild to nonbreeder
    N[4,1] <- 0  # First-year hacked to breeders    
    N[5,1] <- 0  # First-year wild to breeders
    N[6,1] <- 0  # Nonbreeder wild to nonbreeder
    N[7,1] <- 0  # Nonbreeder wild to breeders
    N[8,1] <- 0  # Breeders to nonbreeder
    N[9,1] <- 0  # Breeders to breeders
    N[10,1] <- 0  # Nonbreeder hacked to nonbreeder 
    N[11,1] <- 0  # Nonbreeder hacked to breeder
    N[12,1] <- 0  # Immigrant to breeder
    N[13,1] <- 0  # Immigrant to nonbreeder
    N[14,1] <- 0 # breeder to emigrant
    N[15,1] <- 0 # nonbreeder to emigrant
    N[16,1] <- 0 # first-year to emigrant
    N[17,1] <- 0 # first-year hacked to emigrant
    F[1] <- 0 # No breeding during the first year, all hacked birds
    
    rNB ~ dunif(0,50)
    # omegaB1 ~ dnorm(0, 1/(1*1) )T(0,) # upper limit near 5 imm per breeder to help run model
    # omegaA1 ~ dnorm(0, 1/(1*1) )T(0,)
    # sigma.omegaB ~ dnorm(0, 1/(2*2) )T(0,)
    # sigma.omegaA ~ dnorm(0, 1/(2*2) )T(0,)
    for (m in 1:2){ mu.F[m] ~ dunif(0,5) } # m # limits to help run model
    sigma.F ~ dnorm(0, 1/(2*2) )T(0,)
    mu.em <- logit(mu.em1)    
    mu.em1 ~ dunif(0,1)
    sigma.em ~ dnorm(0, 1/(2*2) )T(0,)
    r ~ dunif(0,1)
    
    # Survival loops for demographic categories by sex, hacked, effort 
    
    sigma.AS.s ~ dnorm(0, 1/(2*2) )T(0,)
    ASalpha <- logit(ASalpha1)
    ASalpha1 ~ dunif(0, 1)
    
    sigma.BS.s ~ dnorm(0, 1/(2*2) )T(0,)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    sigma.BAR.psi ~ dnorm(0, 1/(2*2) )T(0,)
    BARalpha <- logit(BARalpha1)
    BARalpha1 ~ dunif(0, 1)
    
    sigma.pB ~ dnorm(0, 1/(2*2) )T(0,)
    
    for (k in 1:2){
    mu.pB[k]<- logit(mu.pB1[k])
    mu.pB1[k] ~ dunif(0, 1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.ASalpha[t]) <- ASalpha + eps.AS.phi[t]
    eps.AS.phi[t] ~ dnorm(0, 1/(sigma.AS.s*sigma.AS.s) )
    
    logit(eta.BSalpha[t]) <- BSalpha + eps.BS.phi[t]
    eps.BS.phi[t] ~ dnorm(0, 1/(sigma.BS.s*sigma.BS.s) )    
    
    logit(eta.BARalpha[t]) <- BARalpha + eps.BAR.psi[t]
    eps.BAR.psi[t] ~ dnorm(0, 1/(sigma.BAR.psi*sigma.BAR.psi))
    
    logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
    eps.pB[t] ~ dnorm(0, 1/(sigma.pB*sigma.pB) )
    } #t
    
    for(s in 1:2){
    
    sigma.OBR.psi[s] ~ dnorm(0, 1/(2*2) )T(0,)
    OBRalpha[s] <- logit(OBRalpha1[s]) 
    OBRalpha1[s] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.OBRalpha[s,t]) <- OBRalpha[s] + eps.OBR.psi[s,t]
    eps.OBR.psi[s,t] ~ dnorm(0, 1/(sigma.OBR.psi[s]*sigma.OBR.psi[s]) )
    } #t
    
    for (h in 1:2){
    sigma.ABR.psi[s,h] ~ dnorm(0, 1/(2*2) )T(0,)
    ABRalpha[s,h] <- logit(ABRalpha1[s,h])
    ABRalpha1[s,h] ~ dunif(0, 1)
    
    sigma.OS.s[s,h] ~ dnorm(0, 1/(2*2) )T(0,)
    OSalpha[s,h] <- logit(OSalpha1[s,h])
    OSalpha1[s,h] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.ABRalpha[s,h,t]) <- ABRalpha[s,h] + eps.ABR.psi[s,h,t]
    eps.ABR.psi[s,h,t] ~ dnorm(0, 1/(sigma.ABR.psi[s,h]*sigma.ABR.psi[s,h]) )
    
    logit(eta.OSalpha[s,h,t]) <- OSalpha[s,h] + eps.OS.s[s,h,t]
    eps.OS.s[s,h,t] ~ dnorm(0, 1/(sigma.OS.s[s,h]*sigma.OS.s[s,h]) )
    } #t
    } } #h #s
    
    for (h in 1:2){
    sigma.pA[h] ~ dnorm(0, 1/(2*2) )T(0,)
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
    for (t in 2:n.yr){ 
    pp[t] ~ dpois( F[t]*NB[t] )
    log(F[t]) <- log(mu.F[manage[t]]) + eps.F[t]
    eps.F[t] ~ dnorm (0, 1/(sigma.F*sigma.F) )
    }
    # GOF fecundity
    for (t in 1:n.yr){ 
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
    # Likelihood for immigration 
    ###############################
    for (t in 1:(n.yr-1)){
    omegaB[t] <- 0 # set to zero, no immigration
    omegaA[t] <- 0 # set to zero, no immigration
    # log(omegaB[t]) <-  log(omegaB1) + eps.omegaB[t] # breeder
    # log(omegaA[t]) <- log(omegaA1) + eps.omegaA[t] # nonbreeder
    # eps.omegaB[t] ~ dnorm(0, 1/(sigma.omegaB*sigma.omegaB))
    # eps.omegaA[t] ~ dnorm(0, 1/(sigma.omegaA*sigma.omegaA))
    logit(em[t]) <- mu.em + eps.em[t]
    eps.em[t] ~ dnorm(0, 1/(sigma.em*sigma.em))
    } #t
    ################################
    # Likelihood for counts
    ################################
    for (t in 1:(n.yr-1)){
    # Number of wild born juvs
    N[1,t+1] ~ dpois(
    (NH[t]*eta.OSalpha[2,2,t]*eta.OBRalpha[2,t] +
    NW[t]*eta.OSalpha[2,1,t]*eta.OBRalpha[2,t] + # first year males to breeders
    NF[t]*eta.ASalpha[t]*eta.ABRalpha[2,2,t] + # nonbreeder male hacked adults to breeder
    NF[t]*eta.ASalpha[t]*eta.ABRalpha[2,1,t] + # nonbreeder male wild adults to breeder adults
    NB[t]*eta.BSalpha[t]*(1-eta.BARalpha[t])) * # breeders to breeders
    F[t+1]/2) 
    # first year hacked to nonbreeder adults
    N[2,t+1] ~ dpois(eta.OSalpha[2,2,t]*(1-eta.OBRalpha[2,t]) * NH[t] )
    # first year wild to nonbreeder adults
    N[3,t+1] ~ dpois(eta.OSalpha[2,1,t]*(1-eta.OBRalpha[2,t]) * NW[t] )
    # first year hacked to breeders
    N[4,t+1] ~ dpois(eta.OSalpha[2,2,t]*eta.OBRalpha[2,t] * NH[t] )
    # first year wild to breeders
    N[5,t+1] ~ dpois(eta.OSalpha[2,1,t]*eta.OBRalpha[2,t] * NW[t] )
    # Nonbreeder male wild adults to nonbreeder adults
    N[6,t+1] ~ dbin(eta.ASalpha[t]*(1-eta.ABRalpha[2,1,t]), NF[t] )
    # Nonbreeder male wild adults to breeders
    N[7,t+1] ~ dbin(eta.ASalpha[t]*eta.ABRalpha[2,1,t], NF[t] )
    # Breeders to nonbreeder adults
    N[8,t+1] ~ dbin(eta.BSalpha[t]*(eta.BARalpha[t]), NB[t] )
    # Breeders to breeders
    N[9,t+1] ~ dbin(eta.BSalpha[t]*(1-eta.BARalpha[t]), NB[t] )    
    # Nonbreeder hacked adults to nonbreeder adults
    N[10,t+1] ~ dbin(eta.ASalpha[t]*(1-eta.ABRalpha[2,2,t]), NF[t] )
    # Nonbreeder hacked adults to breeders
    N[11,t+1] ~ dbin(eta.ASalpha[t]*eta.ABRalpha[2,2,t], NF[t] )
    # Immigrants to breeders
    N[12,t+1] ~ dpois(NB[t]*omegaB[t])
    # Immigrants to nonbreeders
    N[13,t+1] ~ dpois(NF[t]*omegaA[t])
    # Breeders to emigrants
    N[14,t+1] ~ dpois(em[t] * NB[t] )
    # NonBreeders to emigrants
    N[15,t+1] ~ dpois(em[t] * NF[t] )
    # First-year wild to emigrants
    N[16,t+1] ~ dpois(em[t] * N[1,t] )
    # First-year hacked to emigrants
    N[17,t+1] ~ dpois(em[t] * aug[t] )
    } # t
    
    for (t in 1:n.yr){
    Ntot[t] <- sum(N[c(1,2,3,4,5,6,7,8,9,10,11,12,13),t]) + aug[t] - sum(N[c(14,15,16,17),t]) # total number
    NB[t] <-  sum(N[c(4,5,7,9,11,12),t]) - N[14,t] # number of breeders
    NF[t] <-  sum(N[c(2,3,6,8,10,13),t]) - N[15,t] # number of nonbreeders
    NO[t] <-  N[1,t] + aug[t] - N[16,t] - N[17,t]# number of total first years and translocated first years
    NW[t] <-  N[1,t] - N[16,t] # Number of wild born 
    NH[t] <- aug[t] - N[17,t] # Number of hacked 
    NI[t] <- sum(N[c(12,13),t]) # number of immigrants
    NE[t] <- sum(N[c(14,15,16,17),t])
    } # t
    
    # Observation process    
    for (t in 2:n.yr){
    countBM[t] ~ dnegbin(pNB[t],rNB) # breeding males negative binomial
    pNB[t] <- rNB/(rNB+NB[t])
    countFM[t] ~ dpois(NF[t]) # nonbreeding adult males
    countOM[t] ~ dpois(N[1,t]) # first year males 
    } # t
    
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 2:n.yr){ 
    c.expB[t] <- NB[t] + 0.001 # expected counts adult breeder
    c.expA[t] <- NF[t] + 0.001 # nonbreeder
    c.expO[t] <- N[1,t] + 0.001 # first year
    c.obsB[t] <- countBM[t] + 0.001
    c.obsA[t] <- countFM[t] + 0.001
    c.obsO[t] <- countOM[t] + 0.001
    dssm.obsB[t] <- abs( ( (c.obsB[t]) - (c.expB[t]) ) / (c.obsB[t]+0.001)  )
    dssm.obsA[t] <- abs( ( (c.obsA[t]) - (c.expA[t]) ) / (c.obsA[t]+0.001)  )
    dssm.obsO[t] <- abs( ( (c.obsO[t]) - (c.expO[t]) ) / (c.obsO[t]+0.001)  )
    } # t
    dmape.obs[1] <- sum(dssm.obsB[2:n.yr])
    dmape.obs[2] <- sum(dssm.obsA[2:n.yr])
    dmape.obs[3] <- sum(dssm.obsO[2:n.yr])
    # Compute fit statistic for replicate data
    # Mean absolute error
    for (t in 2:n.yr){ 
    c.repB[t] ~ dnegbin(pNB[t],rNB) # expected counts
    c.repA[t] ~ dpois(NF[t] ) 
    c.repO[t] ~ dpois(N[1,t] )
    } # t
    for (t in 2:n.yr){ 
    dssm.repB[t] <- abs( ( (c.repB[t]) - (c.expB[t]) ) / (c.repB[t]+0.001) )
    dssm.repA[t] <- abs( ( (c.repA[t]) - (c.expA[t]) ) / (c.repA[t]+0.001) )
    dssm.repO[t] <- abs( ( (c.repO[t]) - (c.expO[t]) ) / (c.repO[t]+0.001) )
    } # t
    dmape.rep[1] <- sum(dssm.repB[2:n.yr])
    dmape.rep[2] <- sum(dssm.repA[2:n.yr])
    dmape.rep[3] <- sum(dssm.repO[2:n.yr])
    
    # Test statistic for number of turns
    for (t in 2:(n.yr-2)){
    tt1.obsB[t] <- step(countBM[t+2] - countBM[t+1])
    tt2.obsB[t] <- step(countBM[t+1] - countBM[t])
    tt3.obsB[t] <- equals(tt1.obsB[t] + tt2.obsB[t], 1)
    tt1.obsA[t] <- step(countFM[t+2] - countFM[t+1])
    tt2.obsA[t] <- step(countFM[t+1] - countFM[t])
    tt3.obsA[t] <- equals(tt1.obsA[t] + tt2.obsA[t], 1)
    tt1.obsO[t] <- step(countOM[t+2] - countOM[t+1])
    tt2.obsO[t] <- step(countOM[t+1] - countOM[t])
    tt3.obsO[t] <- equals(tt1.obsO[t] + tt2.obsO[t], 1)
    } # t
    tturn.obs[1] <- sum(tt3.obsB[2:(n.yr-2)])
    tturn.obs[2] <- sum(tt3.obsA[2:(n.yr-2)])
    tturn.obs[3] <- sum(tt3.obsO[2:(n.yr-2)])
    
    for (t in 2:(n.yr-2)){
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
    tturn.rep[1] <- sum(tt3.repB[2:(n.yr-2)])
    tturn.rep[2] <- sum(tt3.repA[2:(n.yr-2)])
    tturn.rep[3] <- sum(tt3.repO[2:(n.yr-2)])
    
    
    ################################
    # Likelihood for survival
    ################################
    for (i in 1:nind){
    for (t in 1:(n.yr-1)){
    #Survival
    sO[i,t] <- eta.OSalpha[sex[i],hacked[i],t] # first year
    sA[i,t] <- eta.ASalpha[t] # nonbreeder
    sB[i,t] <- eta.BSalpha[t] # breeder
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

get.first <- function(x) min(which(x!=5))
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=5))
l <- apply(datl$y, 1, get.last)
TFmat <- is.na(z.inits) & is.na(datl$z)
for (i in 1:dim(TFmat)[1]){  TFmat[i,1:f[i]] <- FALSE }
z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 

inits <- function(){list(z = z.inits)} 

params <- c( "rNB", 
             "sigma.F", "eps.F", "mu.F",
             "r", "l.mu.em", "mu.em", "sigma.em",
             "omegaA1",  "omegaB1", "sigma.omegaA", "sigma.omegaB",
             "OSalpha", "ASalpha", "BSalpha", "OBRalpha", "ABRalpha", "BARalpha",  "mu.pA", "mu.pB",
             "OSalpha1", "ASalpha1", "BSalpha1","OBRalpha1", "ABRalpha1", "BARalpha1","mu.pA1", "mu.pB1",
             "sigma.OS.s", "sigma.AS.s", "sigma.BS.s", "sigma.OBR.psi", "sigma.ABR.psi", "sigma.BAR.psi", "sigma.pA", "sigma.pB",
             "N", "NB", "NF", "NO", "NW", "NH", "NI",  "NE", "Ntot",
             "dmape.obs", "dmape.rep", "tvm.obs", "tvm.rep", "dd.obs", "dd.rep",
             "tturn.obs", "tturn.rep",
             "eps.OS.s", "eps.AS.s", "eps.BS.s", "eps.OBR.psi", "eps.ABR.psi", "eps.BAR.psi", "eps.pA", "eps.pB", 
             "eps.omegaA", "eps.omegaB", "eps.em",
             "eta.OSalpha", "eta.ASalpha", "eta.BSalpha", "eta.OBRalpha", "eta.ABRalpha", "eta.BARalpha", "eta.pA",  "eta.pB", 
             "em", "omegaA", "omegaB", "F"
)

# Set initial values of N 
# to small numbers for immigrant and emigrants
# because preliminary analyses suggested low rates
# This is necessary to run model and avoid initial values error. 
Ni <- array(NA, dim=c(17,26))
Ni[12:17, 2:26] <- 0

inits <- function(){list(z = z.inits, N=Ni)} 

# MCMC settings
ni <- 200000; nt <- 50; nb <- 150000; nc <- 3; na <- 1000 
out <- jags(datl, inits, params, modfl,  
            n.chains = nc, n.thin = nt, n.burnin = nb, n.adapt=na, n.iter=ni, 
            parallel=T, module=c("glm", "bugs"))
#save(file=paste("./", m, ".Rdata", sep=""), list="out")
