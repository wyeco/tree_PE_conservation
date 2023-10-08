## codes for plotting Fig. 4 & 5

load("./data_Fig4&5.RData")

library(raster)
library(RColorBrewer)
library(dplyr)
library(rgdal)
library(ggpubr)
library(ggplot2)



gymno_canape <- raster("./gymno_canape_p_SR5_mask_4.tif")
plot(gymno_canape,col = terrain.colors(4))
angio_canape <- raster("./angio_canape_p_SR5_mask_level4.tif")
plot(angio_canape,col = terrain.colors(4)) # 1 neo, 2 paleo; 3 not sig, 4 Mixed


angio_canape[angio_canape==3 ] <- NA
angio_canape[angio_canape==1] <- 1
angio_canape[angio_canape==2] <- 1
angio_canape[angio_canape==4] <- 1

gymno_canape[gymno_canape==3 ] <- NA
gymno_canape[gymno_canape==1] <- 3
gymno_canape[gymno_canape==2] <- 3
gymno_canape[gymno_canape==4] <- 3
my.palette <- brewer.pal(n = 3, name = "Paired")
plot(angio_canape, col = my.palette)
plot(gymno_canape, col = my.palette)

rsum = sum(angio_canape,gymno_canape,na.rm=TRUE)
rsum[rsum==0]  <- NA
plot(rsum,col= my.palette)

canape_overlap_df <- rsum %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)

#### prepare the nonhotspot cells 
gymno_canape1 <- raster("./gymno_canape_p_SR5_mask_4.tif")
plot(gymno_canape1,col = terrain.colors(4))
angio_canape1 <- raster("./angio_canape_p_SR5_mask_level4.tif")
plot(angio_canape1,col = terrain.colors(4)) # 1 neo, 2 paleo; 3 not sig, 4 Mixed


angio_canape1[angio_canape1==3 ] <- 5
angio_canape1[angio_canape1==1] <- NA
angio_canape1[angio_canape1==2] <- NA
angio_canape1[angio_canape1==4] <- NA

gymno_canape1[gymno_canape1==3 ] <- 3
gymno_canape1[gymno_canape1==1] <- NA
gymno_canape1[gymno_canape1==2] <- NA
gymno_canape1[gymno_canape1==4] <- NA
my.palette <- brewer.pal(n = 3, name = "Paired")
plot(angio_canape1, col = my.palette)
plot(gymno_canape1, col = my.palette)


rsum1 = sum(angio_canape1,gymno_canape1,na.rm=TRUE)
rsum1[rsum1==0]  <- NA
plot(rsum1,col= my.palette)

rsum1[rsum1==3] <- 5
rsum1[rsum1==8] <- 5
plot(rsum1,col= my.palette)
#convert raster to data frame for ggplot2 plotting
canape_nonhotpot_df <- rsum1 %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)

canape_nonhotpot_df$code <- paste0(round(canape_nonhotpot_df$x,2),";",round(canape_nonhotpot_df$y, 2))
canape_overlap_df$code <- paste0(round(canape_overlap_df$x,2),";",round(canape_overlap_df$y, 2))

canape_overlap_nonhotspot_df <- rbind(canape_overlap_df, canape_nonhotpot_df)
####  combine the climate with the obtained overlap layer  #####


summary(angio_gymno_overlap_data_clim_future)
## change the layer name
angio_gymno_overlap_data_clim_future1 <- angio_gymno_overlap_data_clim_future %>% dplyr::mutate(layer =
                                                                                                  case_when(layer == 1 ~ "Angiosperm", 
                                                                                                            layer == 3 ~ "Gymnosperm",
                                                                                                            layer == 4 ~ "Both",
                                                                                                            layer == 5 ~ "Non-hotspot")
)

##plot 
p_both_hmc <- angio_gymno_overlap_data_clim_future1 %>%
  mutate(layer = factor(layer, levels=c("Angiosperm", "Gymnosperm", "Both", "Non-hotspot")))%>%
  ggplot(aes(x=layer, y=gHM, fill=layer)) +
  geom_violin(width=1, trim=T) +
  scale_fill_manual(values = c("Angiosperm" = "#1B9E77",
                               "Gymnosperm" = "#D95F02",
                               "Both"="#7570B3",
                               "Non-hotspot"="#E6AB02")) +
  #scale_fill_brewer(palette="Dark2") +
  geom_boxplot(width=0.05, color="black", alpha=0.8,  fill="white") +
  #theme_classic()+
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour= NA, fill="grey"),
                     panel.background = element_rect(fill = NA, color = "black"),
                     text = element_text(size = 16),
                     legend.position = "none") +
  scale_y_continuous(limits = c(0, 0.85)) +
    theme(axis.title.y=element_text(face="bold",size=14),
        axis.text.x =element_blank(),
        axis.text.y = element_text(face="bold", size=14)) + 
  labs(title="(A) Human modification index",x="", y = "Human modification index")

p_both_future_mat <- angio_gymno_overlap_data_clim_future1 %>%
  mutate(layer = factor(layer, levels=c("Angiosperm", "Gymnosperm", "Both", "Non-hotspot")))%>%
  ggplot(aes(x=layer, y=(bio1_ssp370_median - bio1_curr), fill=layer)) +
  geom_violin(width=1, trim=T) +
  scale_fill_manual(values = c("Angiosperm" = "#1B9E77",
                               "Gymnosperm" = "#D95F02",
                               "Both"="#7570B3",
                               "Non-hotspot"="#E6AB02")) +
  geom_boxplot(width=0.05, color="black", alpha=0.8,  fill="white") +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour= NA, fill="grey"),
                     panel.background = element_rect(fill = NA, color = "black"),
                     text = element_text(size = 16),
                     legend.position = "none") +
  scale_y_continuous(limits = c(-5, 10)) +
  theme(axis.title.y=element_text(face="bold",size=14),
        axis.text.x =element_blank(),
        axis.text.y = element_text(face="bold", size=14)) + 
  labs(title="(B) Future MAT anomaly",x="", y = "Mean annual temperature (oC)")

p_both_future_ap <- angio_gymno_overlap_data_clim_future1 %>%
  mutate(layer = factor(layer, levels=c("Angiosperm", "Gymnosperm", "Both", "Non-hotspot")))%>%
  ggplot(aes(x=layer, y=(bio12_ssp370_median - bio12_curr), fill=layer)) +
  geom_violin(width=1, trim=T) +
  scale_fill_manual(values = c("Angiosperm" = "#1B9E77",
                               "Gymnosperm" = "#D95F02",
                               "Both"="#7570B3",
                               "Non-hotspot"="#E6AB02")) +
  geom_boxplot(width=0.05, color="black", alpha=0.8,  fill="white") +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour= NA, fill="grey"),
                     panel.background = element_rect(fill = NA, color = "black"),
                     text = element_text(size = 16),
                     legend.position = "none") +
  scale_y_continuous(limits = c(-1200, 3200)) +
  theme(axis.title.y=element_text(face="bold",size=14),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) +
  labs(title="(C) Future AP anomaly",x="", y = "Annual precipitation (mm)")



library(rcompanion)
# Human MI
PT_gHM_both = pairwisePermutationTest(gHM ~ layer, data = angio_gymno_overlap_data_clim_future1,
                                      method="fdr")
PT_gHM_both
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gHM_both,
        threshold  = 0.05)

PM_ghm = pairwisePermutationMatrix(gHM ~ layer, data = angio_gymno_overlap_data_clim_future1,
                                   method="fdr")
multcompLetters(PM_ghm$Adjusted,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)

# future MAT
PT_future_mat_both = pairwisePermutationTest(bio1_ssp370_median ~ layer,
                                             data = angio_gymno_overlap_data_clim_future1,
                                             method="fdr")
PT_future_mat_both
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_future_mat_both,
        threshold  = 0.05)

PM = pairwisePermutationMatrix(bio1_ssp370_median ~ layer, data = angio_gymno_overlap_data_clim_future1,
                               method="fdr")
library(multcompView)

multcompLetters(PM$Adjusted,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)

# future AP
PT_future_ap_both = pairwisePermutationTest(bio12_ssp370_median ~ layer, data = angio_gymno_overlap_data_clim_future1,
                                            method="fdr")
PT_future_ap_both
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_future_ap_both,
        threshold  = 0.05)

PM_future_ap = pairwisePermutationMatrix(bio12_ssp370_median ~ layer, data = angio_gymno_overlap_data_clim_future1,
                                         method="fdr")

multcompLetters(PM_future_ap$Adjusted,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)
# 

### obtain the protected coverage

PA_level4 <-  readOGR(dsn = "./WDPA", layer = "WDPA_lev4")  ##removed due to the size
plot(rsum,col= my.palette)
angio_sig <- rsum
gymno_sig <- rsum
both_sig <- rsum
plot(angio_sig)
angio_sig[angio_sig== 3 ]  <- NA
angio_sig[angio_sig== 4 ]  <- NA
gymno_sig[gymno_sig== 1 ]  <- NA
gymno_sig[gymno_sig== 4 ]  <- NA
both_sig[both_sig== 1 ]  <- NA
both_sig[both_sig== 3 ]  <- NA

rsum1## nonhotsopt

par(mar=c(5,5,2,0.1), mfrow=c(1,4))
plot(angio_sig)
plot(gymno_sig)
plot(both_sig)
plot(rsum1)

angio_sig_PA_area <- rasterize(PA_level4, angio_sig, getCover = TRUE)

plot(angio_sig_PA_area)

# inside angio
angio_canape_PA_area_mask <- mask(angio_sig_PA_area, angio_sig)
plot(angio_canape_PA_area_mask)
angio_canape_PA_area_mask_data <- as.data.frame(getValues(angio_canape_PA_area_mask))
hist(angio_canape_PA_area_mask_data$`getValues(angio_canape_PA_area_mask)`)

gumno_canape_PA_area_mask <- mask(angio_sig_PA_area,gymno_sig)
plot(gumno_canape_PA_area_mask)
gumno_canape_PA_area_mask_data <- as.data.frame(getValues(gumno_canape_PA_area_mask))
hist(gumno_canape_PA_area_mask_data$`getValues(gumno_canape_PA_area_mask)`)

both_canape_PA_area_mask <- mask(angio_sig_PA_area,both_sig)
plot(both_canape_PA_area_mask)
both_canape_PA_area_mask_data <- as.data.frame(getValues(both_canape_PA_area_mask))
hist(both_canape_PA_area_mask_data$`getValues(both_canape_PA_area_mask)`)

nonhotspot_canape_PA_area_mask <- mask(angio_sig_PA_area,rsum1)
plot(nonhotspot_canape_PA_area_mask)
nonhotspot_canape_PA_area_mask_data <- as.data.frame(getValues(nonhotspot_canape_PA_area_mask))
hist(nonhotspot_canape_PA_area_mask_data$`getValues(nonhotspot_canape_PA_area_mask)`)

names(angio_canape_PA_area_mask_data)[1] <- "PA"
names(gumno_canape_PA_area_mask_data)[1] <- "PA"
names(both_canape_PA_area_mask_data)[1] <- "PA"
names(nonhotspot_canape_PA_area_mask_data)[1] <- "PA"
angio_canape_PA_area_mask_data$group <- "Angiosperm"  
gumno_canape_PA_area_mask_data$group <- "Gymnosperm"  
both_canape_PA_area_mask_data$group <- "Both"  
nonhotspot_canape_PA_area_mask_data$group <- "Non-hotspot"  
## combine the two as one
all3_combine_canape_PA_data <- rbind(angio_canape_PA_area_mask_data,gumno_canape_PA_area_mask_data,
                                     both_canape_PA_area_mask_data, nonhotspot_canape_PA_area_mask_data)
summary(all3_combine_canape_PA_data)
all3_combine_canape_PA_clean_data <- na.omit(all3_combine_canape_PA_data)
summary(all3_combine_canape_PA_clean_data)
all3_combine_canape_PA_clean_data$group <- as.factor(all3_combine_canape_PA_clean_data$group)

plot(all3_combine_canape_PA_clean_data$PA ~ all3_combine_canape_PA_clean_data$group)
summary(aov(log10(all3_combine_canape_PA_clean_data$PA+2) ~ all3_combine_canape_PA_clean_data$group))


library(Rmisc)  ## to use summarySE function
summary_all3_combine_canape_PA_group <- summarySE(all3_combine_canape_PA_clean_data, measurevar = "PA",
                                                  groupvars=c("group"))
median_all3_combine_canape_PA_group <- ddply(all3_combine_canape_PA_clean_data, "group",
                                             summarise, median = median(PA))

library(ggplot2)
ggplot(all3_combine_canape_PA_clean_data, aes(group, PA, fill=group))  +
  geom_violin(alpha=.5) +
  geom_boxplot(width=.1)


##  plot  ##
p_both_PA <- all3_combine_canape_PA_clean_data %>%
  mutate(group = factor(group, levels=c("Angiosperm", "Gymnosperm", "Both", "Non-hotspot")))%>%
  ggplot(aes(group, PA*100, fill=group)) +
  geom_violin(width=1, trim=T) +
  scale_fill_manual(values = c("Angiosperm" = "#1B9E77",
                               "Gymnosperm" = "#D95F02",
                               "Both"="#7570B3",
                               "Non-hotspot"="#E6AB02")) +
  geom_boxplot(width=0.05, color="black", alpha=0.8,  fill="white") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar",
               width = 0.05, colour = "red",
               position = position_dodge(width = .2)
  ) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = "none") +
  # scale_y_continuous(limits = c(0, 0.8)) +
  theme(axis.title.y=element_text(face="bold",size=14),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) + 
  labs(x=NULL,y="Protected coverage (%)") +
  scale_y_break(c(18, 92), scales = 0.1)

library(rcompanion)
# elevation range
PT_pa_both = pairwisePermutationTest(PA ~ group, data = all3_combine_canape_PA_clean_data,
                                     method="fdr")
PT_pa_both
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_pa_both,
        threshold  = 0.05)

ph_PA = pairwisePermutationMatrix(PA ~ group, data = all3_combine_canape_PA_clean_data,
                                  method="fdr")

multcompLetters(ph_PA$Adjusted,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)


pdf(file = "both angio and gymno hotspot nonhotspot_4gruop.pdf", useDingbats=FALSE, width=7, height=10, onefile = F)
ggarrange(p_both_hmc, p_both_future_mat, p_both_future_ap,
          nrow=3, ncol = 1,align = "v") 
dev.off()

pdf(file = "both angio and gymno hotspot nonhotspot protected coverage_4group.pdf", useDingbats=FALSE, width=8, height=6, onefile = F)
p_both_PA
dev.off()

## import the top 17% priority areas  #############

# richness ranking 
rank_3d_original <- raster("./priority_areas_tree.tif")
MAINL <- readOGR(dsn = "./COUNTRIES", layer = "GSHHS_i_L1_simple")

rank_3d_mask <- raster::mask(rank_3d_original, MAINL)
plot(rank_3d_original)
plot(rank_3d_mask)

### change the conservation maps into binary

rank_3d_top17 <- rank_3d_mask
rank_3d_top17[rank_3d_top17 < 0.83] <- NA
rank_3d_top17[rank_3d_top17 >= 0.83] <- 1
plot(rank_3d_top17)


# inside angio
angio_canape_top17_area_mask <- mask(rank_3d_top17, angio_sig)
plot(angio_canape_top17_area_mask)
angio_canape_top17_area_mask_data <- as.data.frame(getValues(angio_canape_top17_area_mask))
angio_canape_top17_area_mask_data <- na.omit(angio_canape_top17_area_mask_data)
#957/1410=0.6787
gumno_canape_top17_area_mask <- mask(rank_3d_top17,gymno_sig)
plot(gumno_canape_top17_area_mask)
gumno_canape_top17_area_mask_data <- as.data.frame(getValues(gumno_canape_top17_area_mask))
gumno_canape_top17_area_mask_data <- na.omit(gumno_canape_top17_area_mask_data)
#9/25=0.36

both_canape_top17_area_mask <- mask(rank_3d_top17,both_sig)
plot(both_canape_top17_area_mask)
both_canape_top17_area_mask_data <- as.data.frame(getValues(both_canape_top17_area_mask))
both_canape_top17_area_mask_data <- na.omit(both_canape_top17_area_mask_data)
#99/126=0.7857

nonhotspot_canape_top17_area_mask <- mask(rank_3d_top17,rsum1)
plot(nonhotspot_canape_top17_area_mask)
nonhotspot_canape_top17_area_mask_data <- as.data.frame(getValues(nonhotspot_canape_top17_area_mask))
nonhotspot_canape_top17_area_mask_data <- na.omit(nonhotspot_canape_top17_area_mask_data)
##903/5945=0.1519


## check the overlap between PA and the four groups of endemics, just check overlap or not, not percentage as above. thus just 0/1 data
angio_sig_PA <- rasterize(PA_level4, angio_sig, getCover = F)

plot(angio_sig_PA)

# inside angio
angio_canape_PA_mask <- mask(angio_sig_PA, angio_sig)
plot(angio_canape_PA_mask)
angio_canape_PA_mask_data <- as.data.frame(getValues(angio_canape_PA_mask))
angio_canape_PA_mask_data <- na.omit(angio_canape_PA_mask_data)
hist(angio_canape_PA_mask_data$`getValues(angio_canape_PA_mask)`)
#105/1410=0.07446
gumno_canape_PA_mask <- mask(angio_sig_PA,gymno_sig)
plot(gumno_canape_PA_mask)
gumno_canape_PA_mask_data <- as.data.frame(getValues(gumno_canape_PA_mask))
gumno_canape_PA_mask_data <- na.omit(gumno_canape_PA_mask_data)
hist(gumno_canape_PA_mask_data$`getValues(gumno_canape_PA_mask)`)
#0
both_canape_PA_mask <- mask(angio_sig_PA,both_sig)
plot(both_canape_PA_mask)
both_canape_PA_mask_data <- as.data.frame(getValues(both_canape_PA_mask))
both_canape_PA_mask_data <- na.omit(both_canape_PA_mask_data)
hist(both_canape_PA_mask_data$`getValues(both_canape_PA_mask)`)
##11/126=0.0873
nonhotspot_canape_PA_mask <- mask(angio_sig_PA,rsum1)
plot(nonhotspot_canape_PA_mask)
nonhotspot_canape_PA_mask_data <- as.data.frame(getValues(nonhotspot_canape_PA_mask))
nonhotspot_canape_PA_mask_data <- na.omit(nonhotspot_canape_PA_mask_data)
hist(nonhotspot_canape_PA_mask_data$`getValues(nonhotspot_canape_PA_mask)`)
##361/5945= 0.0607

###obtain the list of the four groups
plot(angio_sig)

angio_canape_sig_data <- as.data.frame(getValues(angio_sig))
angio_canape_sig_data <- na.omit(angio_canape_sig_data)


plot(gymno_sig)
gymno_canape_sig_data <- as.data.frame(getValues(gymno_sig))
gymno_canape_sig_data <- na.omit(gymno_canape_sig_data)

plot(both_sig)
both_canape_sig_data <- as.data.frame(getValues(both_sig))
both_canape_sig_data <- na.omit(both_canape_sig_data)


plot(rsum1)
nonhotspot_canape_data <- as.data.frame(getValues(rsum1))
nonhotspot_canape_data <- na.omit(nonhotspot_canape_data)



# richness ranking 
rank_3d_50_binary <- raster("./rank_3d_top50 wgs84.tif")

rank_3d_50_mask <- raster::mask(rank_3d_50_binary, MAINL)
plot(rank_3d_50_binary)
plot(rank_3d_50_mask)

# inside angio
angio_canape_top50_area_mask <- mask(rank_3d_50_mask, angio_sig)
plot(angio_canape_top50_area_mask)
angio_canape_top50_area_mask_data <- as.data.frame(getValues(angio_canape_top50_area_mask))
angio_canape_top50_area_mask_data <- na.omit(angio_canape_top50_area_mask_data)
#1399/1410 =0.9922
gumno_canape_top50_area_mask <- mask(rank_3d_50_mask,gymno_sig)
plot(gumno_canape_top50_area_mask)
gumno_canape_top50_area_mask_data <- as.data.frame(getValues(gumno_canape_top50_area_mask))
gumno_canape_top50_area_mask_data <- na.omit(gumno_canape_top50_area_mask_data)
#24/25=0.96

both_canape_top50_area_mask <- mask(rank_3d_50_mask,both_sig)
plot(both_canape_top50_area_mask)
both_canape_top50_area_mask_data <- as.data.frame(getValues(both_canape_top50_area_mask))
both_canape_top50_area_mask_data <- na.omit(both_canape_top50_area_mask_data)
#126/126=100

nonhotspot_canape_top50_area_mask <- mask(rank_3d_50_mask,rsum1)
plot(nonhotspot_canape_top50_area_mask)
nonhotspot_canape_top50_area_mask_data <- as.data.frame(getValues(nonhotspot_canape_top50_area_mask))
nonhotspot_canape_top50_area_mask_data <- na.omit(nonhotspot_canape_top50_area_mask_data)
##3618/5945=0.6086

################### added on 17May2023, as requested by reviewer to check the top 30% priority areas   ########

rank_3d_top30 <- rank_3d_mask
rank_3d_top30[rank_3d_top30 < 0.70] <- NA
rank_3d_top30[rank_3d_top30 >= 0.70] <- 1
plot(rank_3d_top30)

# inside angio
angio_canape_top30_area_mask <- mask(rank_3d_top30, angio_sig)
plot(angio_canape_top30_area_mask)
angio_canape_top30_area_mask_data <- as.data.frame(getValues(angio_canape_top30_area_mask))
angio_canape_top30_area_mask_data <- na.omit(angio_canape_top30_area_mask_data)
#1275/1410=0.90425
gumno_canape_top30_area_mask <- mask(rank_3d_top30,gymno_sig)
plot(gumno_canape_top30_area_mask)
gumno_canape_top30_area_mask_data <- as.data.frame(getValues(gumno_canape_top30_area_mask))
gumno_canape_top30_area_mask_data <- na.omit(gumno_canape_top30_area_mask_data)
#17/25=0.68

both_canape_top30_area_mask <- mask(rank_3d_top30,both_sig)
plot(both_canape_top30_area_mask)
both_canape_top30_area_mask_data <- as.data.frame(getValues(both_canape_top30_area_mask))
both_canape_top30_area_mask_data <- na.omit(both_canape_top30_area_mask_data)
#116/126=0.9206

nonhotspot_canape_top30_area_mask <- mask(rank_3d_top30,rsum1)
plot(nonhotspot_canape_top30_area_mask)
nonhotspot_canape_top30_area_mask_data <- as.data.frame(getValues(nonhotspot_canape_top30_area_mask))
nonhotspot_canape_top30_area_mask_data <- na.omit(nonhotspot_canape_top30_area_mask_data)
##1968/5945=0.3310


#preparing the data for Fig. 5  #####
df.protection_change <- read.csv("./PA_top17_protection_precentage3.csv", header = T, sep = ',')

#"Angiosperm", "Gymnosperm", "Both", "Non-hotspot"
#"PA","PA_non", "Top17", "Top17_non"


ggplot(df.protection_change, aes(x = group, y = percentage)) +
  geom_col(aes(color = level, fill = level), position = position_dodge(0.8), width = 0.7) +
  scale_color_manual(values = c("#0073C2FF","#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#0073C2FF", "#00AFBB", "#E7B800", "#FC4E07"))+
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = c(0.9,0.8)) +
  # scale_y_continuous(limits = c(0, 100)) +
  theme(axis.title.y=element_text(face="bold",size=14),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14)) + 
  labs(x=NULL,y="Protected coverage (%)") 
pdf(file = "Fig. 5 protection status.pdf", useDingbats=FALSE, width=8, height=6, onefile = F)

dev.off()



