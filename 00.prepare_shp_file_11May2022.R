########################################################
########## code for endemism project  #################
##########   updated on 11 May 2022   ##############
#####################################################


library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(letsR)
library(tmap)
library(sf)
library(viridis)
library(phyloregion)

setwd("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022")
save.image("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022/prepare_shp_files_11May2022_2.RData")
#load("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022/prepare_shp_files_12May2022.RData")
#######   update the tree species list   using GTS v.1.6  11 May 2022   ###########
######## prepare the species list
# import the old tree species list and GTC v.1.6
old_list <- read.csv(file = "Tree_species_2group.csv", header = T, sep = ',')
gtc16_list <- read.csv(file = "Species_wenyong_pep.csv", header = T, sep = ',')
library(dplyr)

check_list <- left_join(old_list, gtc16_list, by ="species")
tree_gtc16_list <- check_list %>%  filter(delete == "FALSE") %>% droplevels()

tree_gtc16_angio <- tree_gtc16_list %>%  filter(group == "Angiosperm") %>% droplevels()
tree_gtc16_gymno <- tree_gtc16_list %>%  filter(group == "Gymnosperm") %>% droplevels()

TC_sp <-  readOGR(dsn = "C:/Users/Wenyong Guo/OneDrive/TreeChange_data/universal_scripts_data",
                  layer = "TC_sum6_46752sp_07Apr2020")

gymno_list <- tree_gtc16_gymno$species

gymno_sp_shp = subset(TC_sp, SOURCE_SHP %in% gymno_list)
plot(gymno_sp_shp)

angio_sp_shp <- subset(TC_sp, !SOURCE_SHP %in% gymno_list)
#plot(angio_sp_shp)

writeOGR(gymno_sp_shp, dsn = "C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022",
         layer = "gymno_sp_shp_11May2022", driver = "ESRI Shapefile")


writeOGR(angio_sp_shp, dsn = "C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022",
         layer = "angio_sp_shp_11May2022",driver = "ESRI Shapefile")

# ### below codes to prepare species list for Biodiverse    ###########

## update TC_sp using the new GTC1.6 list

cp_list <- tree_gtc16_list$species
TC_sp_gtc16 = subset(TC_sp, SOURCE_SHP %in% cp_list)

TC_sp_gtc16
print(bbox(TC_sp_gtc16), digits=12)
set_ll_warn(FALSE)
set_ll_TOL(0.2)
st_crs(TC_sp_gtc16)

colnames(TC_sp_gtc16@data)<- "binomial"
behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
tree.behrmann <- spTransform(TC_sp_gtc16,behrmann)

#To keep running for spatial data with minor range exceed
print(bbox(tree.behrmann), digits=12)
set_ll_warn(FALSE)
set_ll_TOL(0.2)

##the world map2 to limit the range
data(wrld_simpl)
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world,behrmann)
ext <- extent(world)

##Convert species ranges into presence/absence matrix using different resolution
# PAMA6hull_46752_110 <- lets.presab(tree.behrmann,xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],resol=110,
#                                    crs=behrmann,count = T, crs.grid = behrmann, cover = 0.5) ## define the cover bigger than 50%

#save(PAMA6hull_46752_110, file ="PAMA6hull_46752_110.r")# did NOT USE THIS DATA.



##### run code on server, then import it here
PAMA6hull_41835_110_nocover <- lets.presab(tree.behrmann,xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],resol=110,
                                           crs=behrmann,count = T, crs.grid = behrmann) 
save(PAMA6hull_41835_110_nocover, file ="PAMA6hull_41835_110_nocover.Rdata")  # THIS IS THE DATA USED

## hihest richness 1494
plot(PAMA6hull_41835_110_nocover$Richness_Raster)

#### to correct the cell with highest species richness
####Set row names of presence/absence matrix as the cell number of raster

tree.pam_new_110 <- PAMA6hull_41835_110_nocover

rownames(tree.pam_new_110[[1]]) <- cellFromXY(object=tree.pam_new_110$Rich,xy=tree.pam_new_110$Pre[,1:2])

##Detect the cell of outlier in richness. Remove species in the outlier cell, but not in the eight neighbour cells
#The cell of outlier and its eight neighbour cells
id.out <- which(as.vector(tree.pam_new_110$Rich>3500))
row.nr <- rowFromCell(tree.pam_new_110$Rich,id.out)
col.nr <- colFromCell(tree.pam_new_110$Rich,id.out)
id.out.nb <- cellFromRowColCombine(tree.pam_new_110$Rich,(row.nr-1):(row.nr+1),(col.nr-1):(col.nr+1))
id.out.nb <- id.out.nb[!id.out.nb %in% id.out]

#Remove species in the outlier cell, but not in the neighbour 8 cells
id.pre.out <- which(rownames(tree.pam_new_110$Pre) %in% id.out)
id.pre.out.nb <- which(rownames(tree.pam_new_110$Pre) %in% id.out.nb)
occ.out.nb <- colSums(tree.pam_new_110$Pre[id.pre.out.nb,-c(1:2)])
occ.out <- tree.pam_new_110$Pre[id.pre.out,-c(1:2)]
id.spe.out <- which(occ.out.nb==0 & occ.out==1)
tree.pam_new_110$Presence_and_Absence_Matrix[id.pre.out,id.spe.out+2] <- 0

#Update the richnes of outlier cell
rich.out.corrected <- sum(tree.pam_new_110$Pre[id.pre.out,-c(1:2)])
tree.pam_new_110[[2]][id.out] <- rich.out.corrected

## generate the list of all tree species for Biodiverse
presab <- tree.pam_new_110$P
# Print only the first 5 rows and 3 columns
presab[1:15, 1:5]
library(reshape2)

presab <- as.data.frame(presab)

presab_long <- melt(presab, id.vars = c("Longitude(x)", "Latitude(y)"))
library(dplyr)
presab_long_final <- presab_long %>% filter (value > 0) %>% droplevels() 
str(presab_long_final)
summary(presab_long_final)
presab_long_final <- presab_long_final[, -4]

write.csv(presab_long_final, file = "./presab_46752_long_final.csv")
 
# # prepare a copy of the presenceAbsence file, in case I need to run the above script again (it takes time)
# PAMAhull_6mod_new <- PAMAhull_6mod_50
# 
# # # global mainlands (not divided by country boundaries)
# MAINL <- readOGR(dsn = "D:/Wenyong/environment_data/global_tree_S/Data/COUNTRIES", layer = "GSHHS_i_L1_simple")
# MAINL <- spTransform(MAINL, CRSobj = behrmann)
# # # now use the mask function to remove the areas outside of the coastal line
# PAMA6hull_46752_110_mask <- mask(PAMA6hull_46752_110_richness, MAINL)
# hist(PAMA6hull_46752_110_mask)
# 
# # 
# # 
# PAMA6hull_46752_110_mask[PAMA6hull_46752_110_mask < 1] <- NA
# plot(PAMA6hull_46752_110_mask)

### I decide yo use angiosperm and gymnosperm separately, 
###  thus subset the above PAM for each group  ###############
#load the angiosperm species list
## there are 45 wrong tips for angymnosperm, thus I will remove them

angiosperm_wrong_tips <- read.csv("C:/Users/Wenyong Guo/OneDrive/TreeChange_data/universal_scripts_data/TC_phylogeny_wrong_tips_11Mar2020.csv", header = T, sep = ';') # there are 44 species showing wrong location in the phylogeny, thus I decide to remove them at all.
angiosperm_tree_right <- tree_gtc16_angio[!(tree_gtc16_angio$species %in% angiosperm_wrong_tips$species),]

PAM_angiosperm_12May2022 <-lets.subsetPAM(tree.pam_new_110, angiosperm_tree_right$species,remove.cells = T)
plot(PAM_angiosperm_12May2022$Richness_Raster)  ## highest value of 1493
save(PAM_angiosperm_12May2022, file ="PAM_angiosperm_12May2022_removing_WTipS.RData")  # THIS IS THE DATA USED
PAM_gymnosperm_12May2022 <-lets.subsetPAM(tree.pam_new_110, tree_gtc16_gymno$species,remove.cells = T)
plot(PAM_gymnosperm_12May2022$Richness_Raster)  ## highest value of 28
save(PAM_gymnosperm_12May2022, file ="PAM_gymnosperm_12May2022.RData")  # THIS IS THE DATA USED

## export raster for all three layers
PAMA6hull_TC_richness <- tree.pam_new_110$Richness_Raster

PAMA6hull_TC_richness[PAMA6hull_TC_richness < 1] <- NA

PAMA6hull_angio_richness <- PAM_angiosperm_12May2022$Richness_Raster
PAMA6hull_angio_richness[PAMA6hull_angio_richness < 5] <- NA

PAMA6hull_gymno_richness <- PAM_gymnosperm_12May2022$Richness_Raster
PAMA6hull_gymno_richness[PAMA6hull_gymno_richness < 5] <- NA

pdf("Fig.tree SR 12May2022.pdf", useDingbats=FALSE, width=6, height=4)
plot(PAMA6hull_TC_richness,col = viridis(256), main =" Tree richness")
dev.off()

pdf("Fig.tree angio and gymno SR 12May2022.pdf", useDingbats=FALSE, width=9, height=6)
par(mfrow=c(2,2))

plot(PAM_angiosperm_12May2022$Richness_Raster,col = viridis(256), main =" Angiosperm richness")
plot(PAM_gymnosperm_12May2022$Richness_Raster,col = viridis(256), main =" Gymnosperm richness")

plot(PAMA6hull_angio_richness,col = viridis(256), main =" Angiosperm richness < 5")
plot(PAMA6hull_gymno_richness,col = viridis(256), main =" Gymnosperm richness < 5")
cellStats(PAMA6hull_gymno_richness, 'sum')
dev.off()

writeRaster(PAMA6hull_TC_richness, filename="./PAM_treeSpecies_12May2022.tif", overwrite=TRUE)
writeRaster(PAM_angiosperm_12May2022$Richness_Raster, filename="./PAM_angiosperm_12May2022.tif", overwrite=TRUE)
writeRaster(PAM_gymnosperm_12May2022$Richness_Raster, filename="./PAM_gymnosperm_12May2022.tif", overwrite=TRUE)

######prepare data for biodiverse gymno and angio  #####
##angiosperm species

presab_angio <- PAM_angiosperm_12May2022$P
presab_gymno <- PAM_gymnosperm_12May2022$P
# Print only the first 5 rows and 3 columns
presab_angio[1:15, 1:5]
presab_gymno[1:15, 1:5]
library(reshape2)
library(dplyr)


presab_angio_df <- as.data.frame(presab_angio)
presab_angio_df_long <- melt(presab_angio_df, id.vars = c("Longitude(x)", "Latitude(y)"))
library(dplyr)
#presab_angio_df_long_final <- presab_angio_df_long %>% filter (value > 0) %>% droplevels() 
memory.limit(size=40000)  ## allocate more memory here
presab_angio_df_long_final <-  subset(presab_angio_df_long, value == 1, 
                                      select=c("Longitude(x)", "Latitude(y)","variable", "value"))
str(presab_angio_df_long_final)
summary(presab_angio_df_long_final)
presab_angio_df_long_final <- presab_angio_df_long_final[, -4]
write.csv(presab_angio_df_long_final, file = "./presab_angio_df_long_final_12May2022.csv")

## remove the cells with less than 5 species

keep <- rowSums(presab_angio_df[,-(1:2)]) >= 5 
presab_angio_SR5 <- presab_angio_df[keep,]
# or
#presab_angio_SR5_1 <-  presab_angio_df %>% filter(rowSums(.[3:41255]) >= 5)

presab_angio_SR5_long <- melt(presab_angio_SR5, id.vars = c("Longitude(x)", "Latitude(y)"))
presab_angio_SR5_long_final <-  subset(presab_angio_SR5_long, value == 1, 
                                      select=c("Longitude(x)", "Latitude(y)","variable", "value"))
str(presab_angio_SR5_long_final)
summary(presab_angio_SR5_long_final)
presab_angio_SR5_long_final <- presab_angio_SR5_long_final[, -4]
write.csv(presab_angio_SR5_long_final, file = "./presab_angio_SR5_long_final_12May2022.csv")



##  gymnosperm species
presab_gymno_df <- as.data.frame(presab_gymno)
presab_gymno_df_long <- melt(presab_gymno_df, id.vars = c("Longitude(x)", "Latitude(y)"))

library(dplyr)
presab_gymno_df_long_final <- presab_gymno_df_long %>% filter (value > 0) %>% droplevels() 
str(presab_gymno_df_long_final)
summary(presab_gymno_df_long_final)
presab_gymno_df_long_final <- presab_gymno_df_long_final[, -4]
write.csv(presab_gymno_df_long_final, file = "./presab_gymno_df_long_final_12May2022.csv")


### remove the cells with less than 5 species

presab_gymno_SR5 <-  presab_gymno_df %>% filter(rowSums(.[3:562]) >= 5)
presab_gymno_SR5_long <- melt(presab_gymno_SR5, id.vars = c("Longitude(x)", "Latitude(y)"))
presab_gymno_SR5_long_final <- presab_gymno_SR5_long %>% filter (value > 0) %>% droplevels() 
str(presab_gymno_SR5_long_final)
summary(presab_gymno_SR5_long_final)
presab_gymno_SR5_long_final <- presab_gymno_SR5_long_final[, -4]
write.csv(presab_gymno_SR5_long_final, file = "./presab_gymno_SR5_long_final_12May2022.csv")


### prepare phylogenies for both angio and gymno   #####
library(ape)
tc_phylogeny <- read.tree("C:/Users/Wenyong Guo/OneDrive/TreeChange_data/universal_scripts_data/tree_final_11Mar_46708sp.phy.tre")

angio_phylogeny_12May2022 <- keep.tip(tc_phylogeny, angiosperm_tree_right$species)
gymno_phylogeny_12May2022 <- keep.tip(tc_phylogeny, tree_gtc16_gymno$species)

write.tree(angio_phylogeny_12May2022, "angio_phylogeny_12May2022.tre")
write.tree(gymno_phylogeny_12May2022, "gymno_phylogeny_12May2022.tre")

## do not load the above data as it is too big
##### added on 12Nov2021
load("PAM_angiosperm_removing_WTipS.RData")
load("PAM_gymnosperm.RData")

plot(PAM_angiosperm$Richness_Raster)
plot(PAM_gymnosperm$Richness_Raster)
## I need to change all the biodiverse outputs into raster and then mask them with MANIL

# #### gymnosperm first no use ####
# ## make a copy for assign values
# pam_gyno_WE <- PAM_gymnosperm$Richness_Raster
# ## assign resid to the raster
# id_gymno<- as.numeric(rownames(PAM_gymnosperm$Presence_and_Absence_Matrix))
# pam_gyno_WE[id_gymno] <-biodiverse_results_concatenated$
#   plot(pam_gyno_WE)
# 
# pam_gyno_WE_mask <- mask(pam_gyno_WE,MAINL)
# plot(pam_gyno_WE_mask)

