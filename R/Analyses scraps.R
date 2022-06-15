## 1.5 Implementation of the Multi-state CJS Model with Reduced Covariates  

### 1.5.1 JAGS and R Code for the CJS Model with Reduced Covariates
We reduced the number of covariates using probability of direction calculations from the global model by eliminating covariates that were not important (e.g., emigration) when 85% HDIs of differences between group covariates did not overlap zero. 
```{r, include=FALSE, cache=FALSE}
knitr::read_chunk('R/02-run-survival-reduced-JAGS.R')
```
```{r, reduced, eval=FALSE}
```


# DELETE THIS
```{r comp surv IPMs, fig.width=6, fig.height=4, fig.cap= "Fig. S6. Comparisons of survival estimates from an IPM that included immigration and an IPM that excluded immigration."}
## ---- comp surv --------
# load IPM results with immigration
load(".\\outputs\\ipm-ie.Rdata") 
outi <- out
# load IPM results without immigration
load(".\\outputs\\ipm-e.Rdata")
library (MCMCvis)
library (ggplot2)
library(tidybayes)
library (tidyr)
# extract posteriors of IPMs-immigration 
oji <- as.data.frame(outi$sims.list$OSalpha1) # first years, name slightly different due to older run, use OSalpha1
ojil <- oji %>% pivot_longer(cols=c(1,2))
ofi <- as.data.frame(outi$sims.list$ASalpha1) # first years, name slightly different due to older run, use ASalpha1
obi <- as.data.frame(outi$sims.list$BSalpha1)

# extract posteriors of IPMS-no immigration 
# thin posteriors so they have equal lengths
keep <- seq(1, 3000, by=1)
oj <- as.data.frame(out$sims.list$OSalpha1[keep, , ]) # name difference see above
ojl <- oj %>% pivot_longer(cols=c(1,2))
of <- as.data.frame(out$sims.list$ASalpha1[keep]) # name difference see above
ob <- as.data.frame(out$sims.list$BSalpha1[keep])

ojil <- oji %>% pivot_longer(cols=c(1,2))
ojil$cat <- ifelse(ojil$name=="V1", "female", "male")
ojl <- oj %>% pivot_longer(cols=c(1,2))
ojl$cat <- ifelse(ojl$name=="V1", "female", "male")
ofi$cat <- obi$cat <- of$cat <- ob$cat <- "all"
obi$model <- ofi$model <- ojil$model <- "Immigration"
ob$model <- of$model <- ojl$model <- "No immigration"

ojil$stage <- paste("First-year", ojil$cat)
ojl$stage <- paste("First-year", ojl$cat) 
ofi$stage <- of$stage <- "Nonbreeder"
obi$stage <- ob$stage <- "Breeder"

ojil <- ojil[,-1] 
ojl <- ojl[,-1]
colnames(ofi)[1] <- colnames(of)[1] <- 
  colnames(obi)[1] <- colnames(ob)[1] <- "value"

df <- rbind(ojil, ojl, ofi, of, obi, ob)

plt<- df %>%
  ggplot(aes(x = value, y = model, colors=model)) +
  stat_halfeye() +
  facet_wrap(facets="stage", scales="free") +
  theme_classic() +
  ylab("Density") + xlab("Survival probability")
plt
```