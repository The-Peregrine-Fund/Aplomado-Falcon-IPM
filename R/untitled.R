load(".\\data\\data-7states.Rdata")
rm(list=c("model", "params", "inits"))
names(datl)
datl <- datl[-c(9,10,14,15,18,27)]


library (ggplot2)
df <- data.frame(num=c(out$mean$N[1,], datl$aug,  out$mean$NF, out$mean$NB),
           group=c(rep("NO wild",26), rep("Translocated",26), rep("NF",26), rep("NB",26)),
           year=c(1993:2018,1993:2018,1993:2018,1993:2018))
ggplot(df, aes(x=year, y=num, fill=group)) + 
  geom_area()
