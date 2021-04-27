## ---- ipm2 --------
#########################################
# Set the working drive to science drive
##############################
path.df<- data.frame(PCs=NA, Path=NA)
path.df[1,] <- c("rolek.brian", "S:\\")
path.df[2,] <- c("brianrolek", "/run/user/1002/gvfs/smb-share:server=pfsbs,share=science/")
path.df[3,] <- c("cmcclure", "S:\\")
path.df[4,] <- c("nutting", "C:\\Dropbox\\R\\")
path.df[5,] <- c("user", "C:\\Dropbox\\R\\")
ind<- which(path.df$PCs==Sys.info()["user"])
if (any(path.df$PCs==Sys.info()["user"])){ # check PC is in df path.df
  setwd(path.df$Path[ind]) # if TRUE, set path from df
}  else {}
rm(list=c("path.df", "ind"))
setwd("APFA_IPM")
getwd()

#################################
# The model
################################
sink("Analysis/JAGS models/IPM_reduced_no imm_fecfix_effort_reform.txt")
cat("
    model{
    ####################################################
    # DATA
    # Counts
    #  log counts of breeders = l.countBM
    #  log counts of floaters (nonbreeders) = l.countFM
    #  log counts of first year = l.countJM
    # Survival
    # Observations (O) = y  
    # 1 seen Juvenile
    # 2 seen Floater (nonbreeder)
    # 3 seen Breeder
    # 4 not seen
    # States (S)
    # 1 alive Juvenile
    # 2 alive Floater (nonbreeder)
    # 3 alive Breeder
    # 4 dead
    # Groups
    # 1 wild born
    # 2 hacked
    # Sex
    # 1 Female
    # 2 Male
    # Effort
    # 1 Low
    # 2 High
    # Productivity = prod
    ###################################################
    # PARAMETERS
    # F: fecundity
    # phiJ: survival probability First Year (juvenile)
    # phiF: survival probability floaters (nonbreeders)
    # phiB: survival probability breeders
    # psiJB: movement probability from First Year (juvenile) to Breeder
    # psiFB: movement probability from Floater (nonbreeders) to Breeder
    # psiBF: movement probability from Breeder to Floater (nonbreeders) 
    # pJ: resight probability First Year (juvenile)
    # pF: resight probability Floater (nonbreeders)
    # pB: resight probability Breeder
    ###################################################
    ##############
    # Priors and constraints
    ##############
    # Counts for first year are known with certainty
    N[1,1] <- 0 # Juvs
    N[2,1] <- 0 # Juvs to nonbreeder
    N[3,1] <- 0  # Juvs to breeders
    N[4,1] <- 0  # Nonbreeder wild to nonbreeder
    N[5,1] <- 0  # Nonbreeder wild to breeders
    N[6,1] <- 0  # Breeders to nonbreeder
    N[7,1] <- 0  # Breeders to breeders
    N[8,1] <- 0  # Wild juvs
    N[9,1] <- 0  # Nonbreeder hacked to nonbreeder
    N[10,1] <- 0  # Nonbreeder hacked to breeder
    F[1] <- 0 # No breeding first year, all hacked birds

    sigma.BM ~ dunif(0,10)
    sigma.FM ~ dunif(0,10) 
    sigma.JM ~ dunif(0,10) 
    sigma.prod ~dunif(0,10)

    tau.BM <- 1/(sigma.BM*sigma.BM)
    tau.FM <- 1/(sigma.FM*sigma.FM)
    tau.JM <- 1/(sigma.JM*sigma.JM)    
    # tau.omegaB <- 1/(sigma.omegaB*sigma.omegaB)
    # tau.omegaF <- 1/(sigma.omegaF*sigma.omegaF)
    for (j in 1:n.ter){ eps.prod[j] ~ dnorm(0, tau.prod)  } #j
    tau.prod <- 1/(sigma.prod*sigma.prod)

  for (m in 1:2){
    l.mu.F[m] ~ dnorm(0, 0.01)I(-7,2.5)  
    tau.F[m] <-  1/(sigma.F[m]*sigma.F[m])
    sigma.F[m] ~ dunif(0,10)
  } # m
    
    # Survival 
    tau.FS.phi <- 1/(sigma.FS.phi*sigma.FS.phi)
    sigma.FS.phi ~ dunif(0,10)
    FSalpha <- logit(FSalpha1)
    FSalpha1 ~ dunif(0, 1)
    
    tau.BS.phi <- 1/(sigma.BS.phi*sigma.BS.phi)
    sigma.BS.phi ~ dunif(0,10)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    tau.BFR.psi <- 1/(sigma.BFR.psi*sigma.BFR.psi)
    sigma.BFR.psi ~ dunif(0,10)
    BFRalpha <- logit(BFRalpha1)
    BFRalpha1 ~ dunif(0, 1)
    
    tau.pB <- 1/(sigma.pB*sigma.pB)
    sigma.pB ~ dunif(0,10)

for (k in 1:2){
    mu.pB[k]<- logit(mu.pB1[k])
    mu.pB1[k] ~ dunif(0, 1)
} # k

    for (t in 1:(n.yr-1)){
    logit(eta.FSalpha[t]) <- FSalpha + eps.FS.phi[t]
    eps.FS.phi[t] ~ dnorm(0, tau.FS.phi)
    
    logit(eta.BSalpha[t]) <- BSalpha + eps.BS.phi[t]
    eps.BS.phi[t] ~ dnorm(0, tau.BS.phi)    
    
    logit(eta.BFRalpha[t]) <- BFRalpha + eps.BFR.psi[t]
    eps.BFR.psi[t] ~ dnorm(0, tau.BFR.psi)

    logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
    eps.pB[t] ~ dnorm(0, tau.pB)
    } #t
    
    for(s in 1:2){
    tau.JS.phi[s] <- 1/(sigma.JS.phi[s]*sigma.JS.phi[s])
    sigma.JS.phi[s] ~ dunif(0,10)
    JSalpha[s] <- logit(JSalpha1[s]) 
    JSalpha1[s] ~ dunif(0, 1)
    
    tau.JBR.psi[s] <- 1/(sigma.JBR.psi[s]*sigma.JBR.psi[s])
    sigma.JBR.psi[s] ~ dunif(0,10)
    JBRalpha[s] <- logit(JBRalpha1[s]) 
    JBRalpha1[s] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.JSalpha[s,t]) <- JSalpha[s] + eps.JS.phi[s,t]
    eps.JS.phi[s,t] ~ dnorm(0, tau.JS.phi[s])
    
    logit(eta.JBRalpha[s,t]) <- JBRalpha[s] + eps.JBR.psi[s,t]
    eps.JBR.psi[s,t] ~ dnorm(0, tau.JBR.psi[s])
    } #t
    
    for (h in 1:2){
    tau.FBR.psi[s,h] <- 1/(sigma.FBR.psi[s,h]*sigma.FBR.psi[s,h])
    sigma.FBR.psi[s,h] ~ dunif(0,10)
    FBRalpha[s,h] <- logit(FBRalpha1[s,h])
    FBRalpha1[s,h] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.FBRalpha[s,h,t]) <- FBRalpha[s,h] + eps.FBR.psi[s,h,t]
    eps.FBR.psi[s,h,t] ~ dnorm(0, tau.FBR.psi[s,h])
    } #t
    } } #h #s
    
    for (h in 1:2){
    tau.pF[h] <- 1/(sigma.pF[h]*sigma.pF[h])
    sigma.pF[h] ~ dunif(0,10)
    for (k in 1:2){
    mu.pF[k,h] <- logit(mu.pF1[k,h])
    mu.pF1[k,h] ~ dunif(0,1)
} # k

    for (t in 1:(n.yr-1)){
    logit(eta.pF[h,t]) <- mu.pF[effort[t], h] + eps.pF[h,t]
    eps.pF[h,t] ~ dnorm(0, tau.pF[h])
    }} #t #h
    
    #######################
    # Derived params
    #######################
for (t in 2:(n.yr-1)){    
    lambda[t] <-  Ntot[t+1]/(Ntot[t])
    loglambda[t] <- log(lambda[t])
} #t

for (k in 1:K){ prod[k] ~ dpois( prod.mu[year.p[k] , ter[k]] ) }
for (t in 3:(n.yr-1)){ 
  for (j in 1:n.ter){
    log(prod.mu[t,j]) <- log(F[t]) + eps.prod[j] # link productivity to fecundity
    } # j
} # n.yr

    ###############################
    # Likelihood for productivity
    ###############################
for (t in 1:(n.yr-1)){ 
    log(F[t+1]) <- l.mu.F[manage[t+1]] + eps.F[t+1]
    eps.F[t+1] ~ dnorm (0, tau.F[manage[t+1]])

    ################################
    # Likelihood for counts
    ################################
    # Number of wild born juvs
    N[1,t+1] ~ dpois( (NJ[t]*eta.JSalpha[2,t]*eta.JBRalpha[2,t] + # juveniles to breeders
                      NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,2,t] + # nonbreeder male hacked adults to breeder
                      NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,1,t] + # nonbreeder male wild adults to breeder adults
                      NB[t]*eta.BSalpha[t]*(1-eta.BFRalpha[t]) )* # breeders to breeders
                      F[t+1]/2) 
    # Juvs to nonbreeder adults
    N[2,t+1] ~ dbin(eta.JSalpha[2,t]*(1-eta.JBRalpha[2,t]), NJ[t] )
    # Juvs to breeders
    N[3,t+1] ~ dbin(eta.JSalpha[2,t]*eta.JBRalpha[2,t], NJ[t] )
    # Nonbreeder male wild adults to nonbreeder adults
    N[4,t+1] ~ dbin(eta.FSalpha[t]*(1-eta.FBRalpha[2,1,t]), NF[t] )
    # Nonbreeder male wild adults to breeders
    N[5,t+1] ~ dbin(eta.FSalpha[t]*eta.FBRalpha[2,1,t], NF[t] )
    # Breeders to nonbreeder adults
    N[6,t+1] ~ dbin(eta.BSalpha[t]*(eta.BFRalpha[t]), NB[t] )
    # Breeders to breeders
    N[7,t+1] ~ dbin(eta.BSalpha[t]*(1-eta.BFRalpha[t]), NB[t] )    
    # Nonbreeder hacked adults to nonbreeder adults
    N[8,t+1] ~ dbin(eta.FSalpha[t]*(1-eta.FBRalpha[2,2,t]), NF[t] )
    # Nonbreeder hacked adults to breeders
    N[9,t+1] ~ dbin(eta.FSalpha[t]*eta.FBRalpha[2,2,t], NF[t] )

    # Immigrants to breeders
    # N[10,t+1] ~ dpois(NB[t])*omegaB[t])
    # Immigrants to nonbreeders
    # N[11,t+1] ~ dpois(NF[t])*omegaF[t])
    } # t
    
    for (t in 1:n.yr){
    Ntot[t] <- sum(N[c(1,2,3,4,5,6,7,8,9),t]) + aug[t]
    NB[t] <-  sum(N[c(3,5,7,9),t])
    NF[t] <-  sum(N[c(2,4,6,8),t])
    NJ[t] <-  N[1,t] + aug[t]
    } # t
    
# Observation process    
    for (t in 2:n.yr){
    l.countBM[t] ~ dnorm(log(NB[t]), tau.BM) # breeding males
    l.countFM[t] ~ dnorm(log(NF[t]), tau.FM) # nonbreeding adult males
    l.countJM[t] ~ dnorm(log(NJ[t]), tau.JM) # juvenile males
} # t
    
    ################################
    # Likelihood for survival
    ################################
    for (i in 1:nind){
    for (t in 1:(n.yr-1)){
    #Survival
    phiJ[i,t] <- eta.JSalpha[sex[i],t] 
    phiF[i,t] <- eta.FSalpha[t] 
    phiB[i,t] <- eta.BSalpha[t]
    #Recruitment
    psiJB[i,t] <- eta.JBRalpha[sex[i],t] 
    psiFB[i,t] <- eta.FBRalpha[sex[i], hacked[i],t]
    psiBF[i,t] <- eta.BFRalpha[t] 
    #Re-encounter
    pF[i,t] <- eta.pF[hacked[i],t] 
    pB[i,t] <- eta.pB[t] 
    }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in first[i]:(n.yr-1)){
    ps[1,i,t,1] <- 0
    ps[1,i,t,2] <- phiJ[i,t] * (1-psiJB[i,t])
    ps[1,i,t,3] <- phiJ[i,t] * psiJB[i,t]
    ps[1,i,t,4] <- 1-phiJ[i,t]
    
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- phiF[i,t] * (1-psiFB[i,t])
    ps[2,i,t,3] <- phiF[i,t] * psiFB[i,t]
    ps[2,i,t,4] <- 1-phiF[i,t]
    
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- phiB[i,t] * psiBF[i,t]
    ps[3,i,t,3] <- phiB[i,t] * (1-psiBF[i,t])
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
    po[2,i,t,2] <- pF[i,t]
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 1-pF[i,t]
    
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
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,1:4])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1,1:4])
    } #t
    } #i
    
    } #model
    ",fill = TRUE)
sink()