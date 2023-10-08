
load("./data_Fig3b.RData")

library(ggpubr)
library(dplyr)
library(ggplot2)
library(spatialreg)
library(spdep)
library(ncf)   
require(doParallel)

## data
summary(gymno_pe_clim)

##run lm model
gymno_pe_lm <- lm(range_1.1(log10(PE_WE_P)) ~ range_1.1(log10(alt_range +1)) + 
                    range_1.1(mio_prec_anom) + 
                    range_1.1(bio1_curr) +
                    range_1.1(log10(bio12_curr +5))+
                    range_1.1(mio_temp_anom) + 
                    range_1.1(lgm_prec_anom)+  
                    range_1.1(lgm_temp_anom) ,
                  data = gymno_pe_clim)
car::vif(gymno_pe_lm) # checked the vif again
summary(gymno_pe_lm)
AIC(gymno_pe_lm)  ##obtain the AIC
BIC(gymno_pe_lm)
par(mar=c(5,5,2,0.1), mfrow=c(2,2))
plot(gymno_pe_lm)

cor.gymno_pe_lm <-correlog(gymno_pe_clim$Axis_0, gymno_pe_clim$Axis_1, z=residuals(gymno_pe_lm),
                           na.rm=T, increment=1, resamp=1)
#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))
#Plot correlogram
plot(cor.gymno_pe_lm$correlation[1:100], type="b", pch=1, cex=1.2, lwd=1.5,
     ylim=c(-0.5, 1), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2)
abline(h=0)
title(main="PE OLS model residuals", cex=1.5)
# annotate
legend(x=11, y=1, legend=c("PE OLS model residuals"), pch=c(1), bty="n", cex=1.2)


lm.morantest(gymno_pe_lm, spatial_weights_d1, zero.policy = T)## need to run spatial_weights_d1 first below


##############################################################
require(spdep)
#### Make spatial matrices #### 
#### Make a matrix of spatial coordinates (X and Y coordinates) #### 
sp <- SpatialPoints(data.frame(x=gymno_pe_clim$Axis_0, y=gymno_pe_clim$Axis_1))
sp@proj4string <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
coords <- coordinates(sp)
scheme <- "W"
####  get distance to nearest neighbour #### 
k1 = knn2nb(knearneigh(coords,k=1))
dist <- unlist(nbdists(k1, coords))

####  create a series of neighbour matrices based on different distances (d2 is in km) #### 
d1 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)) # min distance to nearest neighbours
d2 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)*2)
d3 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)*4)  
d4 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)*6) 
d5 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)*8) 
d6 <- dnearneigh(sp, longlat = NULL, d1=0, d2=min(dist)*10)


spatial_weights_d1 <- nb2listw(d1, zero.policy=TRUE, style=paste(scheme))
spatial_weights_d2 <- nb2listw(d2, zero.policy=TRUE, style=paste(scheme))
spatial_weights_d3 <- nb2listw(d3, zero.policy=TRUE, style=paste(scheme))
spatial_weights_d4 <- nb2listw(d4, zero.policy=TRUE, style=paste(scheme))
spatial_weights_d5 <- nb2listw(d5, zero.policy=TRUE, style=paste(scheme))
spatial_weights_d6 <- nb2listw(d6, zero.policy=TRUE, style=paste(scheme))



memory.limit(size=10000000)
set.ZeroPolicyOption(TRUE)# as some regions may have no neighbours e.g. islands
scheme <- "W"

#### SAR error model to run in parallel ####
model<-6; # number of models 
cores=detectCores()
cl <- makeCluster(cores[1]-2) 
registerDoParallel(cl)

###########################
#### for now, the model below can not be extract the output due to the package spedep, 
###  I need to use spatialreg package, and then everything works.
detach("package:spdep", unload = TRUE)
library(spatialreg)
RunModel <- function(model.i)
{
  if(model.i==1)
  {
    PE_error_d1 <- errorsarlm(gymno_pe_lm, listw = spatial_weights_d1, tol=1e-12,zero.policy=T)
    return(PE_error_d1)
  }
  if(model.i==2)
  {
    PE_error_d2 <- errorsarlm(gymno_pe_lm,listw = spatial_weights_d2, tol=1e-12,zero.policy=T)
    return(PE_error_d2)
  }
  if(model.i==3)
  {
    PE_error_d3 <<- errorsarlm(gymno_pe_lm,listw = spatial_weights_d3, tol=1e-12,zero.policy=T)
    return(PE_error_d3)
  }
  if(model.i==4)
  { 
    PE_error_d4 <- errorsarlm(gymno_pe_lm,listw = spatial_weights_d4, tol=1e-12,zero.policy=T)
    return(PE_error_d4)
  }
  if(model.i==5)
  { 
    PE_error_d5 <- errorsarlm(gymno_pe_lm,listw = spatial_weights_d5, tol=1e-12,zero.policy=T)
    return(PE_error_d5)
  }
  if(model.i==6)
  { 
    PE_error_d6 <- errorsarlm(gymno_pe_lm,listw = spatial_weights_d6, tol=1e-12,zero.policy=T)
    return(PE_error_d6)
  }
}
clusterExport(cl, list("errorsarlm","gymno_pe_lm","spatial_weights_d1","spatial_weights_d2",
                       "spatial_weights_d3","spatial_weights_d4","spatial_weights_d5","spatial_weights_d6"))
f<-parLapply(cl,1:model, RunModel)
stopCluster(cl)

#############################################################
####save output #### 
require(qpcR)

a<-AIC(f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]])[2]
b<-BIC(f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]])[2]
write.csv(a,"PE_full_akakie_score.csv")
write.csv(akaike.weights(a),"PE_full_akakie_weights.csv")
write.csv(b,"PE_full_BIC_score.csv")
fit<-summary(f[[3]],Nagelkerke = T)  # smallest AIC and BIC model will be selected.
write.csv(cbind(fit$Coef,fit$NK),"PE_full_SAR_summary.csv")
ss<-summary(gymno_pe_lm)
write.csv(cbind(ss$coefficients,ss$r.squared),"lm_PE_full_summary.csv")
## run the best model again

summary(f[[3]], correlation=F, Nagelkerke=TRUE, Hausman=TRUE)
plot(resid(f[[3]]))
BIC(f[[3]])
#Correlograms
library(ncf) 
cor.sem.nb3.gymno_pe_lm <-correlog(gymno_pe_clim$Axis_0, y=gymno_pe_clim$Axis_1,  z=residuals(f[[3]]), na.rm=T, increment=1, resamp=1)

#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))

#Plot correlogram
plot(cor.sem.nb3.gymno_pe_lm$correlation[1:20], type="b", pch=4, cex=1.2, lwd=1.5,
     ylim=c(-0.5, 1), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2)
abline(h=0)
title(main="PD climate anomaly noM2 SAR model residuals", cex=1.5)
library(spdep) # only for below code running 
moran.test(residuals(f[[3]]), spatial_weights_d3, randomisation=T, zero.policy = T)

##############################################################  
### extract the model output and plot   ######

# lm models
library(broom)
library(forcats)
library(ggplot2)

tidy(gymno_pe_lm, conf.int = T)
glance(gymno_pe_lm)
lm_gymno_pe <- as.data.frame(tidy(gymno_pe_lm, conf.int = TRUE))

rownames(lm_gymno_pe) <- c("Intercept", "Elevation range", "Miocene AP anomaly", "MAT", 
                           "AP", "Miocene MAT anomaly", "LGM AP anomaly", "LGM MAT anomaly")


lm_gymno_pe <- dplyr::as_data_frame(lm_gymno_pe, rownames = "Variable")

lm_gymno_pe1 <- as.data.frame(lm_gymno_pe[-1, -2 ])
colnames(lm_gymno_pe1)[2] <- "Estimate"
colnames(lm_gymno_pe1)[3] <- "SE"

lm_gymno_pe1$Estimate <- as.numeric(as.character(lm_gymno_pe1$Estimate))
lm_gymno_pe1$SE <- as.numeric(as.character(lm_gymno_pe1$SE))

lm_gymno_pe1 %>%
  mutate(Variable = fct_relevel(Variable, 
                                "Elevation range", "MAT", "AP","Miocene MAT anomaly","Miocene AP anomaly", 
                                "LGM MAT anomaly", "LGM AP anomaly" )) %>%
  ggplot(aes(x= Variable, y=Estimate)) + 
  geom_point( size=3.5) +
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width=.1)  +                  # Width of the error bars
  geom_hline(yintercept=0, size = 1,linetype="dashed") +
  scale_color_manual(values=c("#0073C2FF", "#FC4E07")) +   
  scale_shape_manual(values = c(15,16)) +
  theme_bw() + theme(#panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(colour= NA, fill="grey"),
    panel.background = element_rect(fill = NA, color = "black"),
    text = element_text(size = 16),
    legend.position = c(0.9,0.8)) +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14)) + coord_flip()

## sar model

tidy(sem.nb3.gymno_pe_lm) # does not work
tidy(f[[3]], conf.int = TRUE)
glance(f[[3]])
augment(f[[3]])

sarlm_gymno_pe <- as.data.frame(tidy(f[[3]], conf.int = TRUE))

rownames(sarlm_gymno_pe) <- c("Intercept", "Elevation range","Miocene AP anomaly", "MAT", 
                              "AP", "Miocene MAT anomaly", "LGM AP anomaly", "LGM MAT anomaly", "lambda")

sarlm_gymno_pe <- dplyr::as_data_frame(sarlm_gymno_pe, rownames = "Variable")

sarlm_gymno_pe1 <- as.data.frame(sarlm_gymno_pe[-c(1,9), -2])
colnames(sarlm_gymno_pe1)[2] <- "Estimate"
colnames(sarlm_gymno_pe1)[3] <- "SE"

sarlm_gymno_pe1$Estimate <- as.numeric(as.character(sarlm_gymno_pe1$Estimate))
sarlm_gymno_pe1$SE <- as.numeric(as.character(sarlm_gymno1$SE))

str(sarlm_gymno_pe1)
### plot the above results

### plot estimates of models withour M2
pd1 <- position_dodge(0.5)

sarlm_gymno_pe1 %>%
  mutate(Variable = fct_relevel(Variable, 
                                "Elevation range", "MAT","AP","Miocene MAT anomaly","Miocene AP anomaly",  
                                "LGM MAT anomaly","LGM AP anomaly")) %>%
  ggplot(aes(x= Variable, y=Estimate)) + 
  geom_point( size=3.5) +
  geom_errorbar(aes(ymin=  conf.low, ymax=  conf.high),
                width=.1) +
  geom_hline(yintercept=0, size = 0.5,linetype="dashed") +
  scale_color_manual(values=c("#0073C2FF", "#FC4E07")) +   
  scale_shape_manual(values = c(15,16)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour= NA, fill="grey"),
                     panel.background = element_rect(fill = NA, color = "black"),
                     text = element_text(size = 16),
                     legend.position = c(0.9,0.9)) +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14)) + coord_flip()

pdf("Fig.SAR gymno_SR5_PE.pdf", useDingbats=FALSE, width=9, height=6)

dev.off()