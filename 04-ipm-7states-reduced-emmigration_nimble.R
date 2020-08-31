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

load("Data\\final-data.Rdata")
#load("/scratch/brolek/aplo_ipm/data/final-data.Rdata")
#load("C:\\Users\\rolek.brian\\Documents\\Projects\\APLO IPM\\Data\\final-data.Rdata")
m<- c("ipm7-reduced-em-nim")

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
    N[1,1] <- 0 # Wild born juvs
    N[2,1] <- 0 # Juvs to nonbreeder
    N[3,1] <- 0  # Juvs to breeders
    N[4,1] <- 0  # Nonbreeder wild to nonbreeder
    N[5,1] <- 0  # Nonbreeder wild to breeders
    N[6,1] <- 0  # Breeders to nonbreeder
    N[7,1] <- 0  # Breeders to breeders
    N[8,1] <- 0 # Nonbreeder hacked adults to nonbreeder adults
    N[9,1] <- 0  # Nonbreeder hacked to breeder
    N[10,1] <- 0  # Immigrants to breeder
    N[11,1] <- 0  # Immigrants to nonbreeder
    N[12,1] <- 0  # Breeder to emigrant
    N[13,1] <- 0  # Nonbreeder to emigrant
    F[1] <- 0 # No breeding first year, all hacked birds

    sigma.BM ~ dunif(0,10)
    sigma.FM ~ dunif(0,10)
    sigma.JM ~ dunif(0,10)
    sigma.prod ~ dunif(0,10)
    sigma.omegaB ~ dunif(0,10)
    sigma.omegaF ~ dunif(0,10)
    sigma.em ~ dunif(0,10)
    for (j in 1:n.ter){ eps.prod[j] ~ dnorm(0, sd=sigma.prod)  } #j

    l.omegaB ~ T(dnorm(0, 0.01), , 3) # upper limit of 20 imm per breeder to help run model
    l.omegaF ~ T(dnorm(0, 0.01), , 3) # upper limit of 20 imm per nonbreeder to help run model
    l.mu.em <- logit(mu.em)
    mu.em ~ dunif(0,1)
    for (m in 1:2){
    l.mu.F[m] ~ T(dnorm(0, 0.01), , 3) # limits to help run model
    sigma.F[m] ~ dunif(0,10)
    } # m
    
    # Survival loops for demographic categories by sex, hacked, effort 
    sigma.FS.s ~ dunif(0,10)
    FSalpha <- logit(FSalpha1)
    FSalpha1 ~ dunif(0, 1)
    
    sigma.BS.s ~ dunif(0,10)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    sigma.BFR.psi ~ dunif(0,10)
    BFRalpha <- logit(BFRalpha1)
    BFRalpha1 ~ dunif(0, 1)
    r ~ dunif(0, 1)
    sigma.pB ~ dunif(0,10)
    sigma.pF ~ dunif(0,10)
    
    for (k in 1:2){
      mu.pB[k]<- logit(mu.pB1[k])
      mu.pB1[k] ~ dunif(0, 1)
      mu.pF[k] <- logit(mu.pF1[k])
      mu.pF1[k] ~ dunif(0,1)
    } # k
    
    for (t in 1:(n.yr-1)){
      logit(eta.FSalpha[t]) <- FSalpha + eps.FS.s[t]
      eps.FS.s[t] ~ dnorm(0, sd=sigma.FS.s)
      logit(eta.BSalpha[t]) <- BSalpha + eps.BS.s[t]
      eps.BS.s[t] ~ dnorm(0, sd=sigma.BS.s)    
      logit(eta.BFRalpha[t]) <- BFRalpha + eps.BFR.psi[t]
      eps.BFR.psi[t] ~ dnorm(0, sd=sigma.BFR.psi)
      logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
      eps.pB[t] ~ dnorm(0, sd=sigma.pB)
    } #t
    
    for(s in 1:2){
      sigma.JS.s[s] ~ dunif(0,10)
      JSalpha[s] <- logit(JSalpha1[s]) 
      JSalpha1[s] ~ dunif(0, 1)
      sigma.JBR.psi[s] ~ dunif(0,10)
      JBRalpha[s] <- logit(JBRalpha1[s]) 
      JBRalpha1[s] ~ dunif(0, 1)
      
      for (t in 1:(n.yr-1)){
        logit(eta.JSalpha[s,t]) <- JSalpha[s] + eps.JS.s[s,t]
        eps.JS.s[s,t] ~ dnorm(0, sd=sigma.JS.s[s])
        logit(eta.JBRalpha[s,t]) <- JBRalpha[s] + eps.JBR.psi[s,t]
        eps.JBR.psi[s,t] ~ dnorm(0, sd=sigma.JBR.psi[s])
      } #t
      
      for (h in 1:2){
        sigma.FBR.psi[s,h] ~ dunif(0,10)
        FBRalpha[s,h] <- logit(FBRalpha1[s,h])
        FBRalpha1[s,h] ~ dunif(0, 1)
        
        for (t in 1:(n.yr-1)){
          logit(eta.FBRalpha[s,h,t]) <- FBRalpha[s,h] + eps.FBR.psi[s,h,t]
          eps.FBR.psi[s,h,t] ~ dnorm(0, sd=sigma.FBR.psi[s,h])
        } } }   #t #h #s
    
    for (t in 1:(n.yr-1)){
      logit(eta.pF[t]) <- mu.pF[effort[t]] + eps.pF[t]
      eps.pF[t] ~ dnorm(0, sd=sigma.pF) 
    } #t 
    
    # #######################
    # # Derived params
    # #######################
    for (t in 3:(n.yr-1)){
    lambda[t] <- Ntot[t+1]/(Ntot[t])
    #loglambda[t] <- log(lambda[t])
    } #t
    # 
    # ###############################
    # # Likelihood for productivity
    # ###############################
    for (k in 1:K){ prod[k] ~ dpois( prod.mu[year.p[k] , ter[k]] ) } # k
    for (t in 3:(n.yr-1)){
    for (j in 1:n.ter){
    log(prod.mu[t,j]) <- log(F[t]) + eps.prod[j] # link productivity to fecundity
    } # j
    } # n.yr
    for (t in 1:(n.yr-1)){
    log(F[t+1]) <- l.mu.F[manage[t+1]] + eps.F[t+1]
    eps.F[t+1] ~ dnorm (0, sd=sigma.F[manage[t+1]])
    # ################################
    # # Likelihood for immigration and emigration
    # ###############################
    log(omegaB[t]) <- eps.omegaB[t] # breeder
    log(omegaF[t]) <- eps.omegaF[t] # nonbreeder
    logit(em[t]) <-  eps.em[t]
    eps.omegaB[t] ~ dnorm(l.omegaB, sd=sigma.omegaB)
    eps.omegaF[t] ~ dnorm(l.omegaF, sd=sigma.omegaF)
    eps.em[t] ~ dnorm(l.mu.em, sd=sigma.em)
    # ################################
    # # Likelihood for counts
    # ################################
    # # Number of wild born juvs
    N[1,t+1] ~ dpois( (NJ[t]*eta.JSalpha[2,t]*eta.JBRalpha[2,t] + # first year males to breeders
    NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,2,t] + # nonbreeder male hacked adults to breeder
    NF[t]*eta.FSalpha[t]*eta.FBRalpha[2,1,t] + # nonbreeder male wild adults to breeder adults
    NB[t]*eta.BSalpha[t]*(1-eta.BFRalpha[t])  ) * F[t+1]/2) # breeders to breeders
    # first year to nonbreeder adults
    N[2,t+1] ~ dbin( eta.JSalpha[2,t]*(1-eta.JBRalpha[2,t]), NJ[t] )
    # first year to breeders
    N[3,t+1] ~ dbin( eta.JSalpha[2,t]*eta.JBRalpha[2,t], NJ[t] )
    # Nonbreeder male wild adults to nonbreeder adults
    N[4,t+1] ~ dbin( eta.FSalpha[t]*(1-eta.FBRalpha[2,1,t]), NF[t] )
    # Nonbreeder male wild adults to breeders
    N[5,t+1] ~ dbin( eta.FSalpha[t]*eta.FBRalpha[2,1,t], NF[t] )
    # Breeders to nonbreeder adults
    N[6,t+1] ~ dbin( eta.BSalpha[t]*(eta.BFRalpha[t]), NB[t] )
    # Breeders to breeders
    N[7,t+1] ~ dbin( eta.BSalpha[t]*(1-eta.BFRalpha[t]), NB[t] )
    # Nonbreeder hacked adults to nonbreeder adults
    N[8,t+1] ~ dbin( eta.FSalpha[t]*(1-eta.FBRalpha[2,2,t]), NF[t] )
    # Nonbreeder hacked adults to breeders
    N[9,t+1] ~ dbin( eta.FSalpha[t]*eta.FBRalpha[2,2,t], NF[t] )
    # Immigrants to breeders
    N[10,t+1] ~ dpois( NB[t]*omegaB[t] )
    # Immigrants to nonbreeders
    N[11,t+1] ~ dpois( NF[t]*omegaF[t] )
    # Breeders to emigrants
    N[12,t+1] ~ dbin( em[t], NB[t] )
    # NonBreeders to emigrants
    N[13,t+1] ~ dbin( em[t], NF[t] )
    # Juveniles to emigrants
    N[14,t+1] ~ dbin( em[t], NJ[t] )
    } # t

    for (t in 1:n.yr){
    Ntot[t] <-  sum(N[1:11,t]) + aug[t] - N[12,t] - N[13,t] - N[14,t] # total number
    NB[t] <- N[3,t] + N[5,t] + N[7,t] + N[9,t] + N[10,t] - N[12,t] # number of breeders
    NF[t] <- N[2,t] + N[4,t] + N[6,t] + N[8,t] + N[11,t] - N[13,t] # number of nonbreeders
    NJ[t] <- N[1,t] + aug[t] - N[14,t] # number of first years
    NE[t] <- N[12,t] + N[13,t] + N[14,t] # number of emigrants
    NI[t] <- N[10,t] + N[11,t] # number of immigrants
    } # t

    # Observation process
    for (t in 2:n.yr){
    countBM[t] ~ dpois(NB[t])#dnorm(log(NB[t]), sd=sigma.BM) # breeding males
    countFM[t] ~ dpois(NF[t])#dnorm(log(NF[t]), sd=sigma.FM) # nonbreeding adult males
    countJM[t] ~ dpois(NJ[t])#dnorm(log(NJ[t]), sd=sigma.JM) # first year males
    } # t
    
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
        ps[1,i,t,2] <- sJ[i,t] * (1-psiJB[i,t]) * (1-em[t])
        ps[1,i,t,3] <- sJ[i,t] * psiJB[i,t] * (1-em[t])
        ps[1,i,t,4] <- (1-sJ[i,t]) * r * (1-em[t])
        ps[1,i,t,5] <- (1-sJ[i,t]) * (1-r) * (1-em[t])
        ps[1,i,t,6] <- sJ[i,t] * em[t]
        ps[1,i,t,7] <- (1-sJ[i,t]) * (1-r) * em[t]
        
        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- sF[i,t] * (1-psiFB[i,t]) * (1-em[t])
        ps[2,i,t,3] <- sF[i,t] * psiFB[i,t] * (1-em[t])
        ps[2,i,t,4] <- (1-sF[i,t]) * r * (1-em[t])
        ps[2,i,t,5] <- (1-sF[i,t]) * (1-r) * (1-em[t])
        ps[2,i,t,6] <- sF[i,t] * em[t]
        ps[2,i,t,7] <- (1-sF[i,t]) * (1-r) * em[t]
        
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- sB[i,t] * psiBF[i,t] * (1-em[t])
        ps[3,i,t,3] <- sB[i,t] * (1-psiBF[i,t]) * (1-em[t])
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
                         JSalpha1 = runif(2), 
                         FSalpha1 = runif(1),  
                         BSalpha1 = runif(1),
                         JBRalpha1 = runif(2), 
                         FBRalpha1 = array(runif(4), dim=c(2,2)), 
                         BFRalpha1 = runif(1),
                         mu.pF1 = runif(2),
                         mu.pB1 = runif(2),
                         l.mu.F = runif(2),
                         l.omegaB = runif(-5, 0), 
                         l.omegaF = runif(-5, 0),
                         sigma.prod = runif(1),
                         sigma.F = c(0.01, 0.01),
                         sigma.JS.s = runif(2), 
                         sigma.FS.s = runif(1), 
                         sigma.BS.s = runif(1),
                         sigma.JBR.psi = runif(2), 
                         sigma.FBR.psi = array(runif(4), dim=c(2,2)), 
                         sigma.BFR.psi = runif(1),
                         sigma.pF = runif(1), 
                         sigma.pB = runif(1),
                         sigma.BM = runif(1), 
                         sigma.FM = runif(1), 
                         sigma.JM = runif(1),
                         sigma.omegaB = runif(1), 
                         sigma.omegaF = runif(1), 
                         sigma.em = runif(1),
                         NB = datl$countBM,
                         NF = datl$countFM,
                         NJ = datl$countJM,
                         r= runif(1, 0.01, 0.1),
                         mu.em = runif(1, 0.01, 0.1),
                         eps.FS.s= runif(datl$n.yr, -1, 1),
                         eps.BS.s= runif(datl$n.yr, -1, 1),
                         eps.JS.s= array(runif(2*datl$n.yr, -1, 1), dim=c(2,datl$n.yr)), 
                         eps.BFR.psi= runif(datl$n.yr, -1, 1),
                         eps.FBR.psi= array(runif(2*2*datl$n.yr, -1, 1), dim=c(2,2,datl$n.yr)), 
                         eps.JBR.psi= array(runif(2*datl$n.yr, -1, 1), dim=c(2,datl$n.yr)),
                         eps.pB= runif(datl$n.yr, -1, 1),
                         eps.pF= runif(datl$n.yr, -1, 1), 
                         eps.F= runif(datl$n.yr, -1, 1), 
                         eps.prod= runif(datl$n.ter, -1, 1),
                         eps.omegaB= runif(datl$n.yr-1, -1, 1),
                         eps.omegaF= runif(datl$n.yr-1, -1, 1),
                         eps.em= runif(datl$n.yr-1, -1, 1),
                         sB= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         sF= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         sJ= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         psiJB= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         psiFB= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         psiBF= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         em= runif(datl$n.yr),
                         pB= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         pF= array(runif(datl$nind*datl$n.yr) , dim=c(datl$nind, datl$n.yr) ),
                         eta.BSalpha= runif(datl$n.yr, 0.7, 0.9),
                         eta.FSalpha= runif(datl$n.yr, 0.6, 0.9),
                         eta.JSalpha= array(runif(datl$n.yr*2, 0.1, 0.3), dim=c(2,datl$n.yr)),
                         eta.BFRalpha= runif(datl$n.yr),
                         eta.FBRalpha= array(runif(datl$n.yr*2*2), dim=c(2,2,datl$n.yr)),
                         eta.JBRalpha= array(runif(datl$n.yr*2), dim=c(2,datl$n.yr)),
                         eta.pB= runif(datl$n.yr),
                         eta.pF= runif(datl$n.yr),
                         N= array(rpois(13*datl$n.yr, lambda=5), dim=c(13, datl$n.yr) )
)} 

params <- c(
              "sigma.BM", "sigma.FM", "sigma.JM",
              "l.mu.F", "F", "eps.F", "sigma.F", "eps.prod", "sigma.prod", "prod.mu",
              "l.omegaB", "l.omegaF", "omegaB", "omegaF", "l.mu.em", "mu.em", "em",
              "eps.omegaB", "eps.omegaF", "eps.em", "sigma.omegaB", "sigma.omegaF", "sigma.em",
              # Survival params
              "r",
              "JSalpha", "FSalpha", "BSalpha", "JBRalpha", "FBRalpha", "BFRalpha",  "mu.pF", "mu.pB",
              "JSalpha1", "FSalpha1", "BSalpha1","JBRalpha1", "FBRalpha1", "BFRalpha1","mu.pF1", "mu.pB1",
              "eps.JS.s", "eps.FS.s", "eps.BS.s", "eps.JBR.psi", "eps.FBR.psi", "eps.BFR.psi", "eps.pF", "eps.pB", 
              "eta.JSalpha", "eta.FSalpha", "eta.BSalpha", "eta.JBRalpha", "eta.FBRalpha", "eta.BFRalpha", "eta.pF",  "eta.pB", 
              "sigma.JS.s", "sigma.FS.s", "sigma.BS.s", "sigma.JBR.psi", "sigma.FBR.psi", "sigma.BFR.psi", "sigma.pF", "sigma.pB",  
              # Population params
              "lambda", 
              "NB", "NF", "NJ", "NE", "NI", "N", "Ntot"
)

constl <- datl[ c(2,3,4,17,24,25) ] 
datl2 <- datl[ c(1,5,6,12,13,16,20,21,22,23,26,27)  ]
# MCMC settings
ni <- 200000; nt <- 50; nb <- 100000; nc <- 3; na <- 1000
ni <- 10000; nt <- 5; nb <- 5000; nc <- 3; na <- 1000
#ni <- 50; nt <- 1; nb <- 25; nc <- 1; na <- 10
mod<- nimbleModel(code, calculate=T, constants = constl, 
                  data = datl2, inits = inits())
mod$initializeInfo()
# cmod <- compileNimble(mod)
# 
# AFSS <- configureMCMC(cmod)
# AFSS_mcmc <- buildMCMC(AFSS)
# AFSSc <- compileNimble(AFSS_mcmc, project = cmod, resetFunctions = TRUE)
# out <- runMCMC(AFSSc,
#         niter=ni, 
#         nchains=nc,
#         thin = nt,
#         nburnin = nb)

out <- nimbleMCMC(  model= mod,
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

#save(file=paste("/scratch/brolek/aplo_ipm/modeloutputs/", m, ".Rdata", sep=""), list="out")
save(file=paste("Analysis\\NIMBLE results\\", m, ".Rdata", sep=""), list="out")


# change MCMC sampler to 
# Automated Factor Slice Sampler
# to increase speed and efficiency
#AFSS2 <- configureMCMC(mod , autoBlock=T, autoIt=10000)

#AFSS$removeSamplers(params)

# tar <- c()
# for (i in 1:length(AFSS$getSamplers())){
# tar[i]  <- AFSS$getSamplers()[[i]]$targetAsScalar
# }
# tar <- unique(tar)
# AFSS$addSampler(target = c(tar, params[c(1:4,6:9,11:12,16:23,32:47,56:63)]), # exclude derived and intermediate params
#                 type = 'AF_slice')