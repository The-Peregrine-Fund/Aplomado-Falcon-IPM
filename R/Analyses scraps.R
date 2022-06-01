## 1.5 Implementation of the Multi-state CJS Model with Reduced Covariates  

### 1.5.1 JAGS and R Code for the CJS Model with Reduced Covariates
We reduced the number of covariates using probability of direction calculations from the global model by eliminating covariates that were not important (e.g., emigration) when 85% HDIs of differences between group covariates did not overlap zero. 
```{r, include=FALSE, cache=FALSE}
knitr::read_chunk('R/02-run-survival-reduced-JAGS.R')
```
```{r, reduced, eval=FALSE}
```
