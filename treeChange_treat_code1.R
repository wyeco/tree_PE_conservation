########################################################
##########  tree phylogenetic endemism project #########
##########   1: data preparation          ##############
########################################################


library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(letsR)
library(tmap)
library(sf)
library(viridis)
library(phyloregion)
library(dplyr)
library(reshape2)
library(ape)

setwd("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022")
save.image("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022/prepare_shp_files_final.RData")
# load("C:/file_EcoInfo_group_04June2021/PE_RPE_CANPE/new_11May2022/prepare_shp_files_12May2022.RData")
#######   Update the tree species list using GTS v.1.6  11 May 2022   ###########
######## Prepare the species list
# Import the old tree species list and GTC v.1.6
old_list <- read.csv(file = "Tree_species.v1.csv", header = T, sep = ',')
gtc16_list <- read.csv(file = "Tree_species_gtc.v1.6.csv", header = T, sep = ',')


check_list <- left_join(old_list, gtc16_list, by ="species")
tree_gtc16_list <- check_list %>%  filter(delete == "FALSE") %>% droplevels()

tree_gtc16_angio <- tree_gtc16_list %>%  filter(group == "Angiosperm") %>% droplevels()
tree_gtc16_gymno <- tree_gtc16_list %>%  filter(group == "Gymnosperm") %>% droplevels()

## import the alphahull shp file
TC_sp <-  readOGR(dsn = "./Rdata",
                  layer = "TC_sum6_46752sp_07Apr2020")


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


##### run code on server, then import it here
PAMA6hull_41835_110_nocover <- lets.presab(tree.behrmann,xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],resol=110,
                                           crs=behrmann,count = T, crs.grid = behrmann) 
save(PAMA6hull_41835_110_nocover, file ="PAMA6hull_41835_110_nocover.Rdata")  # THIS IS THE DATA USED

plot(PAMA6hull_41835_110_nocover$Richness_Raster)
tree.pam_new_110 <- PAMA6hull_41835_110_nocover

### separate angiosperm and gymnosperm species and subset the above PAM for each group  ###############
angiosperm_wrong_tips <- read.csv("./Rdata/TC_phylogeny_wrong_tips.csv", header = T, sep = ';') 
angiosperm_tree_right <- tree_gtc16_angio[!(tree_gtc16_angio$species %in% angiosperm_wrong_tips$species),]

PAM_angiosperm <-lets.subsetPAM(tree.pam_new_110, angiosperm_tree_right$species,remove.cells = T)
PAM_gymnosperm <-lets.subsetPAM(tree.pam_new_110, tree_gtc16_gymno$species,remove.cells = T)


######prepare data for biodiverse for gymno and angio species #####
##angiosperm species

presab_angio <- PAM_angiosperm$P
presab_gymno <- PAM_gymnosperm$P
# Print only the first 5 rows and 3 columns
presab_angio[1:15, 1:5]
presab_gymno[1:15, 1:5]

presab_angio_df <- as.data.frame(presab_angio)
presab_angio_df_long <- melt(presab_angio_df, id.vars = c("Longitude(x)", "Latitude(y)"))
memory.limit(size=40000) 
presab_angio_df_long_final <-  subset(presab_angio_df_long, value == 1, 
                                      select=c("Longitude(x)", "Latitude(y)","variable", "value"))
str(presab_angio_df_long_final)
summary(presab_angio_df_long_final)
presab_angio_df_long_final <- presab_angio_df_long_final[, -4]
write.csv(presab_angio_df_long_final, file = "./presab_angio_df_long_final.csv")

## remove the cells with less than 5 species
presab_angio_SR5 <-  presab_angio_df %>% filter(rowSums(.[3:41255]) >= 5)
presab_angio_SR5_long <- melt(presab_angio_SR5, id.vars = c("Longitude(x)", "Latitude(y)"))
presab_angio_SR5_long_final <-  subset(presab_angio_SR5_long, value == 1, 
                                      select=c("Longitude(x)", "Latitude(y)","variable", "value"))
str(presab_angio_SR5_long_final)
summary(presab_angio_SR5_long_final)
presab_angio_SR5_long_final <- presab_angio_SR5_long_final[, -4]
write.csv(presab_angio_SR5_long_final, file = "./presab_angio_SR5_long_final.csv")

##  gymnosperm species
presab_gymno_df <- as.data.frame(presab_gymno)
presab_gymno_df_long <- melt(presab_gymno_df, id.vars = c("Longitude(x)", "Latitude(y)"))
presab_gymno_df_long_final <- presab_gymno_df_long %>% filter (value > 0) %>% droplevels() 
str(presab_gymno_df_long_final)
summary(presab_gymno_df_long_final)
presab_gymno_df_long_final <- presab_gymno_df_long_final[, -4]
write.csv(presab_gymno_df_long_final, file = "./presab_gymno_df_long_final.csv")

### remove the cells with less than 5 species
presab_gymno_SR5 <-  presab_gymno_df %>% filter(rowSums(.[3:562]) >= 5)
presab_gymno_SR5_long <- melt(presab_gymno_SR5, id.vars = c("Longitude(x)", "Latitude(y)"))
presab_gymno_SR5_long_final <- presab_gymno_SR5_long %>% filter (value > 0) %>% droplevels() 
str(presab_gymno_SR5_long_final)
summary(presab_gymno_SR5_long_final)
presab_gymno_SR5_long_final <- presab_gymno_SR5_long_final[, -4]
write.csv(presab_gymno_SR5_long_final, file = "./presab_gymno_SR5_long_final.csv")

### prepare phylogeny for each of angio and gymno   #####
tc_phylogeny <- read.tree("./Rdata/tree_final_46708sp.phy.tre")
angio_phylogeny <- keep.tip(tc_phylogeny, angiosperm_tree_right$species)
gymno_phylogeny <- keep.tip(tc_phylogeny, tree_gtc16_gymno$species)

write.tree(angio_phylogeny, "angio_phylogeny.tre")
write.tree(gymno_phylogeny, "gymno_phylogeny.tre")


