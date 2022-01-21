## ---- comp surv --------
load("G:\\APLO_IPM\\outputs\\IPM_reduced_imm_fecfix_effort_reform.Rdata")
outi <- out
load("G:\\APLO_IPM\\outputs\\IPM_reduced_no imm_fecfix_effort_reform_update.Rdata")
library (MCMCvis)
library (ggplot2)
library(tidybayes)
library (tidyr)
oji <- as.data.frame(outi$sims.list$JSalpha1)
ojil <- oji %>% pivot_longer(cols=c(1,2))
ofi <- as.data.frame(outi$sims.list$FSalpha1)
obi <- as.data.frame(outi$sims.list$BSalpha1)

# thin posteriors so they have equal lengths
keep <- seq(1, 27000, by=9)
oj <- as.data.frame(out$sims.list$JSalpha1[keep, ])
ojl <- oj %>% pivot_longer(cols=c(1,2))
of <- as.data.frame(out$sims.list$FSalpha1[keep])
ob <- as.data.frame(out$sims.list$BSalpha1[keep])

ojil <- oji %>% pivot_longer(cols=c(1,2))
ojil$cat <- ifelse(ojil$name=="V1", "female", "male")
ojl <- oj %>% pivot_longer(cols=c(1,2))
ojl$cat <- ifelse(ojl$name=="V1", "female", "male")
ofi$cat <- obi$cat <- of$cat <- ob$cat <- "all"
ojil$model <- "Immigration"
ofi$model <- "Immigration"
obi$model <- "Immigration"
ojl$model <- "No immigration"
of$model <- "No immigration"
ob$model <- "No immigration"

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

ggsave("C:\\Users\\rolek.brian\\Documents\\Projects\\APLO IPM\\docs\\figs\\IPM_comparison-imm.tiff",
       device="tiff",
       plot=plt, width=6, height=4, units="in",
       dpi=320)
