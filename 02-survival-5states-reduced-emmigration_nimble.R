load("/scratch/brolek/aplo_ipm/data/final-data.Rdata")
#load("C:\\Users\\rolek.brian\\Documents\\Projects\\APLO IPM\\Data\\final-data.Rdata")
m<- c("surv5-reduced-fid-nim")

library (nimble)
#######################
# modelcode
######################
code <- nimbleCode(
  {
    ####################################################
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
    # sJ: survival probability juveniles (first year age 0)
    # sF: survival probability floaters (nonbreeders)
    # sB: survival probability breeders
    # psiJB: movement probability from Juvenile (first year age zero) to Breeder
    # psiFB: movement probability from Floater (nonbreeders) to Breeder
    # psiBF: movement probability from Breeder to Floater (nonbreeders) 
    # pJ: resight probability Juvenile (First Year)
    # pF: resight probability Floater (nonbreeders)
    # pB: resight probability Breeder
    ###################################################
    
    ##############
    # Priors and constraints
    ##############
    # Counts for first year are known with certainty
    # N[1,1] <- 0 # Wild born juvs
    # N[2,1] <- 0 # Juvs to nonbreeder
    # N[3,1] <- 0  # Juvs to breeders
    # N[4,1] <- 0  # Nonbreeder wild to nonbreeder
    # N[5,1] <- 0  # Nonbreeder wild to breeders
    # N[6,1] <- 0  # Breeders to nonbreeder
    # N[7,1] <- 0  # Breeders to breeders
    # N[8,1] <- 0 # Nonbreeder hacked adults to nonbreeder adults
    # N[9,1] <- 0  # Nonbreeder hacked to breeder
    # N[10,1] <- 0  # Immigrants to breeder
    # N[11,1] <- 0  # Immigrants to nonbreeder
    # F[1] <- 0 # No breeding first year, all hacked birds
    # 
    # sigma.BM ~ dunif(0,10)
    # sigma.FM ~ dunif(0,10)
    # sigma.JM ~ dunif(0,10)
    # sigma.prod ~ dunif(0,10)
    # sigma.omegaB ~ dunif(0,10)
    # sigma.omegaF ~ dunif(0,10)
    # 
    # tau.BM <- 1/(sigma.BM*sigma.BM)
    # tau.FM <- 1/(sigma.FM*sigma.FM)
    # tau.JM <- 1/(sigma.JM*sigma.JM)
    # tau.omegaB <- 1/(sigma.omegaB*sigma.omegaB)
    # tau.omegaF <- 1/(sigma.omegaF*sigma.omegaF)
    # for (j in 1:n.ter){ eps.prod[j] ~ dnorm(0, tau.prod)  } #j
    # tau.prod <- 1/(sigma.prod*sigma.prod)
    # 
    # l.omegaB ~ dnorm(0, 0.01)I(-7 ,2.5) # upper limit of 20 imm per breeder to help run model
    # l.omegaF ~ dnorm(0, 0.01)I(-7 ,2.5) # upper limit of 20 imm per nonbreeder to help run model
    # 
    # for (m in 1:2){
    # l.mu.F[m] ~ dnorm(0, 0.01)I(-7, 2.5) # limits to help run model
    # tau.F[m] <-  1/(sigma.F[m]*sigma.F[m])
    # sigma.F[m] ~ dunif(0,10)
    # } # m
    
    r ~ dunif(0, 1)
    Fid ~ dunif(0, 1)
    
    # Survival loops for demographic categories by sex, hacked, effort 
    tau.FS.s <- 1/(sigma.FS.s*sigma.FS.s)
    sigma.FS.s ~ dunif(0,10)
    FSalpha <- logit(FSalpha1)
    FSalpha1 ~ dunif(0, 1)
    
    tau.BS.s <- 1/(sigma.BS.s*sigma.BS.s)
    sigma.BS.s ~ dunif(0,10)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    tau.BFR.psi <- 1/(sigma.BFR.psi*sigma.BFR.psi)
    sigma.BFR.psi ~ dunif(0,10)
    BFRalpha <- logit(BFRalpha1)
    BFRalpha1 ~ dunif(0, 1)
    
    tau.pB <- 1/(sigma.pB*sigma.pB)
    sigma.pB ~ dunif(0,10)
    tau.pF <- 1/(sigma.pF*sigma.pF)
    sigma.pF ~ dunif(0,10)
    
    for (k in 1:2){
      mu.pB[k]<- logit(mu.pB1[k])
      mu.pB1[k] ~ dunif(0, 1)
      mu.pF[k] <- logit(mu.pF1[k])
      mu.pF1[k] ~ dunif(0,1)
    } # k
    
    for (t in 1:(n.yr-1)){
      logit(eta.FSalpha[t]) <- FSalpha + eps.FS.s[t]
      eps.FS.s[t] ~ dnorm(0, tau.FS.s)
      
      logit(eta.BSalpha[t]) <- BSalpha + eps.BS.s[t]
      eps.BS.s[t] ~ dnorm(0, tau.BS.s)    
      
      logit(eta.BFRalpha[t]) <- BFRalpha + eps.BFR.psi[t]
      eps.BFR.psi[t] ~ dnorm(0, tau.BFR.psi)
      
      logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
      eps.pB[t] ~ dnorm(0, tau.pB)
    } #t
    
    for(s in 1:2){
      tau.JS.s[s] <- 1/(sigma.JS.s[s]*sigma.JS.s[s])
      sigma.JS.s[s] ~ dunif(0,10)
      JSalpha[s] <- logit(JSalpha1[s]) 
      JSalpha1[s] ~ dunif(0, 1)
      
      tau.JBR.psi[s] <- 1/(sigma.JBR.psi[s]*sigma.JBR.psi[s])
      sigma.JBR.psi[s] ~ dunif(0,10)
      JBRalpha[s] <- logit(JBRalpha1[s]) 
      JBRalpha1[s] ~ dunif(0, 1)
      
      for (t in 1:(n.yr-1)){
        logit(eta.JSalpha[s,t]) <- JSalpha[s] + eps.JS.s[s,t]
        eps.JS.s[s,t] ~ dnorm(0, tau.JS.s[s])
        
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
        } } }   #t #h #s
    
    for (t in 1:(n.yr-1)){
      logit(eta.pF[t]) <- mu.pF[effort[t]] + eps.pF[t]
      eps.pF[t] ~ dnorm(0, tau.pF) 
    } #t 
    
    # #######################
    # # Derived params
    # #######################
    # for (t in 2:(n.yr-1)){
    # lambda[t] <-  Ntot[t+1]/(Ntot[t])
    # loglambda[t] <- log(lambda[t])
    # } #t
    # # 
    # # ###############################
    # # # Likelihood for productivity
    # # ###############################
    # for (k in 1:K){ prod[k] ~ dpois( prod.mu[year.p[k] , ter[k]] ) }
    # for (t in 3:(n.yr-1)){
    # for (j in 1:n.ter){
    # log(prod.mu[t,j]) <- log(F[t]) + eps.prod[j] # link productivity to fecundity
    # } # j
    # } # n.yr
    # for (t in 1:(n.yr-1)){
    # log(F[t+1]) <- l.mu.F[manage[t+1]] + eps.F[t+1]
    # eps.F[t+1] ~ dnorm (0, tau.F[manage[t+1]])
    # # ################################
    # # # Likelihood for immigration
    # # ###############################
    # log(omegaB[t+1]) <-  l.omegaB + eps.omegaB[t+1] # breeder
    # log(omegaF[t+1]) <- l.omegaF + eps.omegaF[t+1] # nonbreeder
    # eps.omegaB[t+1] ~ dnorm(0, tau.omegaB)
    # eps.omegaF[t+1] ~ dnorm(0, tau.omegaF)
    # # ################################
    # # # Likelihood for counts
    # # ################################
    # # # Number of wild born juvs
    # N[1,t+1] ~ dpois( (NJ[t]*eta.JSalpha[2,t]*eta.JBRalpha[2,t] + # first year males to breeders
    # NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,2,t] + # nonbreeder male hacked adults to breeder
    # NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,1,t] + # nonbreeder male wild adults to breeder adults
    # NB[t]*eta.BSalpha[t]*(1-eta.BFRalpha[t]) )* # breeders to breeders
    # F[t+1]/2)
    # # first year to nonbreeder adults
    # N[2,t+1] ~ dbin(eta.JSalpha[2,t]*(1-eta.JBRalpha[2,t]), NJ[t] )
    # # first year to breeders
    # N[3,t+1] ~ dbin(eta.JSalpha[2,t]*eta.JBRalpha[2,t], NJ[t] )
    # # Nonbreeder male wild adults to nonbreeder adults
    # N[4,t+1] ~ dbin(eta.FSalpha[t]*(1-eta.FBRalpha[2,1,t]), NF[t] )
    # # Nonbreeder male wild adults to breeders
    # N[5,t+1] ~ dbin(eta.FSalpha[t]*eta.FBRalpha[2,1,t], NF[t] )
    # # Breeders to nonbreeder adults
    # N[6,t+1] ~ dbin(eta.BSalpha[t]*(eta.BFRalpha[t]), NB[t] )
    # # Breeders to breeders
    # N[7,t+1] ~ dbin(eta.BSalpha[t]*(1-eta.BFRalpha[t]), NB[t] )
    # # Nonbreeder hacked adults to nonbreeder adults
    # N[8,t+1] ~ dbin(eta.FSalpha[t]*(1-eta.FBRalpha[2,2,t]), NF[t] )
    # # Nonbreeder hacked adults to breeders
    # N[9,t+1] ~ dbin(eta.FSalpha[t]*eta.FBRalpha[2,2,t], NF[t] )
    # # Immigrants to breeders
    # N[10,t+1] ~ dpois( (NB[t]+NF[t])*omega[t+1])
    # # Immigrants to nonbreeders
    # N[11,t+1] ~ dpois( (NB[t]+NF[t])*Fid[t+1])
    # } # t
    # 
    # for (t in 1:n.yr){
    # Ntot[t] <-  sum(N[c(1,2,3,4,5,6,7,8,9,10),t]) + aug[t] - N[11,t] # dsum(N[1,t], N[2,t], N[3,t], N[4,t], N[5,t], N[6,t], N[7,t], N[8,t], N[9,t], N[10,t], N[11,t], aug[t]) ## total number
    # NB[t] <- sum(N[c(3,5,7,9,10),t]) #~ dsum(N[3,t], N[5,t], N[7,t], N[9,t], N[10,t])  # # number of breeders
    # NF[t] <- sum(N[c(2,4,6,8),t]) #~ dsum(N[2,t], N[4,t], N[6,t], N[8,t], N[11,t])  # # number of nonbreeders
    # NJ[t] <- N[1,t] + aug[t] #~ dsum(N[1,t] , aug[t])   # # number of first years
    # } # t
    # 
    # # Observation process
    # for (t in 2:n.yr){
    # l.countBM[t] ~ dnorm(log(NB[t]), tau.BM) # breeding males
    # l.countFM[t] ~ dnorm(log(NF[t]), tau.FM) # nonbreeding adult males
    # l.countJM[t] ~ dnorm(log(NJ[t]), tau.JM) # first year males
    # } # t
    # 
    ################################
    # Likelihood for survival
    ################################
    for (i in 1:nind){
      for (t in 1:(n.yr-1)){
        #Survival
        sJ[i,t] <- eta.JSalpha[sex[i],t] # first year
        sF[i,t] <- eta.FSalpha[t] # nonbreeder
        sB[i,t] <- eta.BSalpha[t] # breeder
        #Recruitment
        psiJB[i,t] <- eta.JBRalpha[sex[i],t] # first year to breeder
        psiFB[i,t] <- eta.FBRalpha[sex[i], hacked[i],t] # nonbreederto breeder
        psiBF[i,t] <- eta.BFRalpha[t] # breeder to nonbreeder
        #Re-encounter
        pF[i,t] <- eta.pF[t] # resight of nonbreeders
        pB[i,t] <- eta.pB[t]  # resight of breeders
      }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of state S(t+1) given S(t)
      for (t in first[i]:(n.yr-1)){
        ps[1,i,t,1] <- 0
        ps[1,i,t,2] <- sJ[i,t] * (1-psiJB[i,t]) * Fid
        ps[1,i,t,3] <- sJ[i,t] * psiJB[i,t] * Fid
        ps[1,i,t,4] <- (1-sJ[i,t]) * r * Fid
        ps[1,i,t,5] <- (1-sJ[i,t]) * (1-r) * Fid
        ps[1,i,t,6] <- sJ[i,t] * (1-Fid)
        
        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- sF[i,t] * (1-psiFB[i,t]) * Fid
        ps[2,i,t,3] <- sF[i,t] * psiFB[i,t] * Fid
        ps[2,i,t,4] <- (1-sF[i,t]) * r * Fid
        ps[2,i,t,5] <- (1-sF[i,t]) * (1-r) * Fid
        ps[2,i,t,6] <- sF[i,t] * (1-Fid)
        
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- sB[i,t] * psiBF[i,t] * Fid
        ps[3,i,t,3] <- sB[i,t] * (1-psiBF[i,t]) * Fid
        ps[3,i,t,4] <- (1-sB[i,t]) * r * Fid
        ps[3,i,t,5] <- (1-sB[i,t]) * (1-r) * Fid
        ps[3,i,t,6] <- sB[i,t] * (1-Fid)
        
        ps[4,i,t,1] <- 0
        ps[4,i,t,2] <- 0
        ps[4,i,t,3] <- 0
        ps[4,i,t,4] <- 0
        ps[4,i,t,5] <- 1
        ps[4,i,t,6] <- 0
        
        ps[5,i,t,1] <- 0
        ps[5,i,t,2] <- 0
        ps[5,i,t,3] <- 0
        ps[5,i,t,4] <- 0
        ps[5,i,t,5] <- 1
        ps[5,i,t,6] <- 0
        
        ps[6,i,t,1] <- 0
        ps[6,i,t,2] <- 0
        ps[6,i,t,3] <- 0
        ps[6,i,t,4] <- 0
        ps[6,i,t,5] <- 0
        ps[6,i,t,6] <- 1
        
        # Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- 1 
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 0
        po[1,i,t,4] <- 0
        po[1,i,t,5] <- 0
        
        po[2,i,t,1] <- 0
        po[2,i,t,2] <- pF[i,t]
        po[2,i,t,3] <- 0
        po[2,i,t,4] <- 0
        po[2,i,t,5] <- 1-pF[i,t]
        
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
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,first[i]] <- y[i,first[i]]
      for (t in (first[i]+1):n.yr){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:6])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:5])
      } #t
    } #i
    
  } #model
) 


# Initial values
get.first <- function(x) min(which(x!=5))
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=5))
l <- apply(datl$y, 1, get.last)
TFmat <- is.na(z.inits) & is.na(datl$z)
for (i in 1:dim(TFmat)[1]){  TFmat[i,1:f[i]] <- FALSE }
z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 

inits <- function(){list(z = z.inits,
                         JSalpha1= runif(2), 
                         FSalpha1= runif(1),  
                         BSalpha1= runif(1),
                         JBRalpha1= runif(2), 
                         FBRalpha1= array(runif(4), dim=c(2,2)), 
                         BFRalpha1= runif(1),
                         mu.pF1=runif(2),
                         mu.pB1=runif(2),
                         #l.mu.F = c(0.055, 0.535),
                         #l.omegaB= log(0.01), l.omegaF= log(0.01),
                         #sigma.F = c(0.01, 0.01),
                         sigma.JS.s = runif(2), 
                         sigma.FS.s = runif(1), 
                         sigma.BS.s = runif(1),
                         sigma.JBR.psi = runif(2), 
                         sigma.FBR.psi = array(runif(4), dim=c(2,2)), 
                         sigma.BFR.psi = runif(1),
                         sigma.pF = runif(1), 
                         sigma.pB = runif(1),
                         #sigma.BM = 0.01, sigma.FM = 0.01, sigma.JM = 0.01,
                         #sigma.omegaB=0.01, sigma.omegaF=0.01,
                         #NB= ifelse(exp(datl$l.countBM)<1,1, exp(datl$l.countBM)),
                         #NF= ifelse(exp(datl$l.countFM)<1,1, exp(datl$l.countFM)),
                         #NJ= ifelse(exp(datl$l.countJM)<1,1, exp(datl$l.countJM)),
                         r= runif(1, 0.01, 0.1),
                         Fid= runif(1, 0.5, 1)
)} 

params <- c(
  #"sigma.BM", "sigma.FM", "sigma.JM",
  #"l.mu.F", "F", "eps.F", "sigma.F", "eps.prod", "sigma.prod", "prod.mu",
  #"l.omegaB", "l.omegaF", "omegaB", "omegaF", "eps.omegaB", "eps.omegaF",
  # Survival params
  "r", "Fid",
  "JSalpha", "FSalpha", "BSalpha", "JBRalpha", "FBRalpha", "BFRalpha",  "mu.pF", "mu.pB",
  "JSalpha1", "FSalpha1", "BSalpha1","JBRalpha1", "FBRalpha1", "BFRalpha1","mu.pF1", "mu.pB1",
  "eps.JS.s", "eps.FS.s", "eps.BS.s", "eps.JBR.psi", "eps.FBR.psi", "eps.BFR.psi", "eps.pF", "eps.pB", 
  "eta.JSalpha", "eta.FSalpha", "eta.BSalpha", "eta.JBRalpha", "eta.FBRalpha", "eta.BFRalpha", "eta.pF",  "eta.pB", 
  "sigma.JS.s", "sigma.FS.s", "sigma.BS.s", "sigma.JBR.psi", "sigma.FBR.psi", "sigma.BFR.psi", "sigma.pF", "sigma.pB"  
  # Population params
  #"lambda", "loglambda",
  #"NB", "NF", "NJ",  "N", "Ntot"
)

constl <- datl[ c(2,3,4) ] 
datl2 <- datl[ c(1,5,6,26,27)  ]
# MCMC settings
ni <- 100000; nt <- 50; nb <- 50000; nc <- 3; na <- 1000
#ni <- 10000; nt <- 5; nb <- 5000; nc <- 3; na <- 1000
#ni <- 50; nt <- 1; nb <- 25; nc <- 1; na <- 10
mod<- nimbleModel(code, calculate=T, constants = constl, 
                  data = datl2, inits = inits())

out <- nimbleMCMC(  model=mod,
                    code = code,
                    monitors = params,
                    nchains=nc,
                    thin = nt,
                    niter = ni,
                    nburnin = nb,
                    progressBar=T,
                    summary=T,
                    WAIC=F,
                    samplesAsCodaMCMC = T,
                    samples=T )

save(file=paste("/scratch/brolek/aplo_ipm/modeloutputs/", m, ".Rdata", sep=""), list="out")