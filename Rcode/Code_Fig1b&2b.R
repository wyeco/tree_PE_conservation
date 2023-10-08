## Using Biodiverse R pipeline to run Biodiverse and plot the results
## this pipeline is from https://github.com/NunzioKnerr/biodiverse_pipeline ###
####   biodiverse_path_reference.R   ######
Sys.setenv(PATH="C:/Strawberry/perl/bin;C:/Strawberry/c/bin;%PATH%") #set path to srawberryperl for 64 bit version
Sys.setenv(PERL5LIB="./biodiverse/lib")
biodiverse_install_folder <- "./biodiverse/"
biodiverse_pipeline_install_folder <- "./biodiverse_pipeline/"

setwd("./data_biodiverse")
save.image("./data_Fig1b&2b.RData")
load("./data_Fig1b&2b.RData")


#This script checks and installs a number of packages if they are not already installed on the system
#to use it just select all of the text and run it in R. It basically speeds up the process of getting an R install up and running so you don't need to manually install packages.
#
#Nunzio.Knerr@csiro.au
#updated Date:19/06/2015

#below is a list of packages to check and install
stdpkgs <- c("sp", "maptools", "XML", "gridExtra", "Cairo", "rgdal", "extrafont", "grid",  "rgeos", "raster", "plyr", "dplyr", "tidyr", "ggplot2", "RColorBrewer", "colorspace", "colorRamps", "spacetime", "aqp", "spatstat", "scales", "stringr", "gWidgets", "phytools", "ape", "apTreeshape") 
otherpkgs <- c("plotKML", "maps", "mapdata") 
#
#Do not edit below

#check and install libraries
# script to install needed packages for R Windows
#  -- set CRAN mirror locally
#setInternet2() # this sets an alternative internet. good for proxy if stupid intenet blocking is done
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.csiro.au/"
options(repos=r)
})

allpkgs <- union(stdpkgs, otherpkgs)

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

lapply(allpkgs, pkgTest)


### I used standalone software "Biodiverse" to prepare the data, thus, here I will not use the 2ed-5th steps in the pipeline  ####

######open_biodiverse_run_analyses.R  RUN THIS SCRIPT 6th#####
#This script uses a perl script to open biodiverse file and run the analyses you want.
#Inputs:
#input_bds_file: The full path to the biodiverse file to use (.bds)
#input_bds_file: The full path to the tree file to use (.bts)
#calcs: A comma delimited list of the calculations/analyses you want to compute.
#       List of calc names can be found on the biodiverse indicies page
#       https://code.google.com/p/biodiverse/wiki/Indices
#       listed as Subroutine: 
#       For example "Subroutine: calc_endemism_central" 
#Nunzio.Knerr@csiro.au
#Date:30/05/2015
#

# source("./R_release/biodiverse_path_reference.R")
input_bds_file <- paste0("./gymno_SR5.bds")
input_bts_file <- paste0("./gymno_phylogeny.bts")
calcs <- paste("calc_endemism_whole,calc_pd,calc_pe,calc_phylo_rpd1,calc_phylo_rpd2,calc_phylo_rpe1,calc_phylo_rpe2")#,calc_phylo_rpe2_branch_stats
#calcs = paste("calc_numeric_label_stats")
cmd <- paste ("perl ", biodiverse_pipeline_install_folder, "perl/run_analyses.pl", sep="")
#
###### do not edit below 
cmd <- paste(cmd, "--input_bds_file", shQuote(input_bds_file), "--input_bts_file", shQuote(input_bts_file), "--calcs", shQuote(calcs), sep =" ")

print(cmd)

system(cmd, intern=T, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE) 

######   open_biodiverse_run_randomisations.R RUN THIS SCRIPT 7th #####
#This script uses a perl script to open biodiverse and run randomisations.
#Note: the perl script is part of the biodiverse distribution and resides there.
#Inputs:
#csv_file: the file conating the spatial data to load in
#out_file: the full path to the biodiverse basedata file (.bds) that should be generated.
#label_column_number: the column number for the labels, typically the taxon name
#group_column_number_x: the x co-ordinate column, could be latitude or metres
#group_column_number_y: the y co-ordinate column, could be longitude or metres
#cell_size_x: the x size of the group cell (in the same unit as the x co-ordinates)
#cell_size_y: the y size of the group cell (in the same unit as the y co-ordinates)
#
#
#Nunzio.Knerr@csiro.au
#Date:30/05/2015
#
# source("./R_release/biodiverse_path_reference.R")


basedata <- paste0("./gymno_SR5_analysed.bds")
rand_name <- paste0("rand")
iterations <- 999
args <- paste0("function=rand_structured max_iters=999")
#
###### do not edit below
cmd1 <- paste("perl ", biodiverse_install_folder, "bin/run_randomisation.pl", sep="")

cmd1 <- paste(cmd1, "--basedata", shQuote(basedata), "--rand_name", shQuote(rand_name), "--iterations", shQuote(iterations), "--args", args)

print(cmd1)

system(cmd1, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE) 


########open_biodiverse_export_results_from_bds # RUN THIS SCRIPT 8th########################

#This script uses a perl script to open biodiverse file and output results as csv.
#Inputs:
#input_bds_file: The biodiverse file (.bds) to load and export results from
#output_csv_prefix: the text to add to the begingin of the csv file name
#Nunzio.Knerr@csiro.au
#Date:1/07/2015
#source("./R_release/biodiverse_path_reference.R")

cmd2 <- paste ("perl ", biodiverse_pipeline_install_folder, "perl/load_bds_and_export_results.pl", sep="")
input_bds_file <-  paste0("./gymno_SR5_analysed.bds")
output_csv_prefix <-  paste0("./gymno_SR5_analysed_output")

###### do not edit below 
cmd2 <- paste(cmd2, "--input_bds_file", shQuote(input_bds_file), "--output_csv_prefix", shQuote(output_csv_prefix), sep =" ")

print(cmd2)
system(cmd2, intern=T, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE) 

########load_biodiverse_results_and_plot_maps_on_world_map # RUN THIS SCRIPT 9th########################

#
#This script loads the biodiverse results output files (csv's) calculates the CANAPE significances and plots maps.
#Note: This script downlaods a world shape file and then allows you to choose the region that you want to plot. 
#Inputs:
#input_csv_files: The biodiverse results .csv's, 2 files, 1 of spatial results, the other randomisations
#several input parameters: 
#   base_dir <- base direcory of the pipeline eg. 'C:/biodiverse_pipeline/'
#   data_dir <- the directory where your data is stored eg. 'C:/my_documents/my_data/'
#   output_folder  <- the directory you want the results save in eg. 'C:/my_documents/my_data/output/'
#   print_seperate_images <- "TRUE/FALSE"
#   output_PNG <- "TRUE/FALSE"
#   output_PDF <- "TRUE/FALSE"
#   NOTE: data and map projection need to be the same if you want it to line up and map properly
#   newproj <- paste0(" +init=epsg:3310") the epsg code for the map to be projected into, should match your data and be in equal area coordinates
#   region_to_map <-  paste0("USA-3521") # specify the region you want to plot has to match the field in the shape file used
#   use: summary(worldLowres) to see fields
#   edit : single.region.map <- worldLowres[worldLowres$adm1_cod_1 == region_to_map,] to match field you want
#outputs: several maps in the specified output directory either png,pdf or both
#
#
#Nunzio.Knerr@csiro.au
#Date:30/03/2016
#
##########################################################################################################

library(sp)    
library(maptools) 
library(raster)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(gridExtra)
library(Cairo)
library(extrafont)
library(rgdal)

#loadfonts()

myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available

base_dir <- './biodiverse_pipeline/'  
#  CHANGE THIS PATH TO POINT TO YOUR INSTALLATION
data_dir <- paste0("./gymno/")
output_folder  <- paste0("./gymno/output/")
setwd(data_dir)

fname_spatial_results <- "gymno_SR5_analysed_output_SPATIAL_RESULTS.csv"
fname_rand_results    <- "gymno_SR5_analysed_output_rand--SPATIAL_RESULTS.csv"

observed_data_file  <- paste0(data_dir, fname_spatial_results)
rand_results_file <- paste0(data_dir, fname_rand_results)

observed_data <- read.csv(observed_data_file)
rand_results  <- read.csv(rand_results_file)
biodiverse_results_concatenated <- cbind(observed_data, rand_results)
#View(biodiverse_results_concatenated)

print_seperate_images <- TRUE
output_PNG <- TRUE
output_PDF <- TRUE

#############################################################################
# Functions for populating new columns with text of significance
#############################################################################
#Standard 2 tailed test for RPD
significance_fun <- function(x){
  if (x >= 0.99) {
    return("Very Highly Sig")
  } else if (x >= 0.975){
    return ("Highly Sig")
  } else if (x <= 0.01){
    return ("Very Sig Low")
  } else if (x <= 0.025){
    return ("Sig Low")
  } else {
    return("Not Sig")
  }
}


## I modified the code to remove the super group
significance_super_fun <- function(x, y, z){
  if (x < 0.95 & y < 0.95) {
    return("Not Sig")
  } else if (z <= 0.025){
    return ("Neo")
  } else if (z >= 0.975){
    return ("Paleo")
    # } else if (x >= 0.99 & y >= 0.99){
    #   return ("Super")
  } else {
    return("Mixed")
  }
}


# # global mainlands (not divided by country boundaries)
library(tidyverse)
library(sf)

behrmann <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
MAINL <- readOGR(dsn = "./Data/COUNTRIES", layer = "GSHHS_i_L1_simple")  ## coastline
MAINL <- spTransform(MAINL, CRSobj = behrmann)

map_data <- MAINL
plot(map_data)


###################################################################################
#Create new columns in dataframe and populate them using the functions above
###################################################################################

targets <- c("PHYLO_RPD1", "PHYLO_RPD2", "PD_P", "PE_WE_P", "PD_P_per_taxon", "PHYLO_RPE2")

for (name in targets) {
  colname <- paste0("P_", name)  #  prepend the P_ since we want the proportions, saves some typing above
  new_colname = paste0(colname, "_SIG")
  trait_index <- match (colname, colnames(biodiverse_results_concatenated))
  # Apply the function to every row of column with index "trait_index" 
  #  and generate a column in the dataframe showing significant cells
  if (!is.na(trait_index)) {
    biodiverse_results_concatenated[[new_colname]] <- apply (biodiverse_results_concatenated[trait_index],  MARGIN=c(1), significance_fun) 
  } else {
    print (paste("Cannot find index", colname, "in data frame"))
  }
}

#This uses the 2 pass test to pull out palaeo, neo and super for RPE
#  SWL - could be refactored into a function
biodiverse_results_concatenated$P_PHYLO_RPE1_CANAPE_SIG <- sapply(
  1:nrow(biodiverse_results_concatenated), 
  function(x) significance_super_fun(
    biodiverse_results_concatenated$P_PE_WE_P[x], 
    biodiverse_results_concatenated$P_PHYLO_RPE_NULL1[x], 
    biodiverse_results_concatenated$P_PHYLO_RPE1[x]
  )
)

biodiverse_results_concatenated$P_PHYLO_RPE2_CANAPE_SIG <- sapply(
  1:nrow(biodiverse_results_concatenated), 
  function(x) significance_super_fun(
    biodiverse_results_concatenated$P_PE_WE_P[x], 
    biodiverse_results_concatenated$P_PHYLO_RPE_NULL2[x], 
    biodiverse_results_concatenated$P_PHYLO_RPE2[x]
  )
)

summary(as.factor(biodiverse_results_concatenated$P_PHYLO_RPE2_CANAPE_SIG))

###########################

theme_pipeline <- function() {
  theme(text = element_text(family=myFont),
        strip.background=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key.height=unit(1.1,"cm"),
        legend.margin=unit(2,"cm"),
        legend.key.width=unit(6.2,"cm"),
        legend.position=c(.5,0.93),
        legend.direction='horizontal',
        legend.title=element_text(colour='black',angle=0,size=rel(4),family=myFont),
        legend.text=element_text(colour='black',angle=0,size=rel(4),family=myFont),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.background= element_blank(),
        plot.margin=unit(c(0,0,-0.61,-0.61),"line")
  )
}


##the world map2 to limit the range
data(wrld_simpl)  
world <- wrld_simpl[wrld_simpl@data$LAT>= -70,]
world <- spTransform(world,behrmann)
ext <- extent(world)

max_x <- ext[2]# extent of map + space for legend 
min_x <- ext[1] # other extent of map 
max_y <- ext[4] # extent of map + space for legend
min_y <- ext[3] # other extent of map


## generate a module for all result layers
r_obj <- raster(xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y, resolution=c(110,110))   

########## SR spcies richness #####
# use rasterize to create desired raster
r_data <- rasterize(x=biodiverse_results_concatenated[, 2:3], # lon-lat data
                    y=r_obj, # raster object
                    field=biodiverse_results_concatenated$ENDW_RICHNESS) # vals to fill raster with
plot(r_data)
crs(r_data) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
pam_gyno_SR_mask <- mask(r_data,MAINL)
plot(pam_gyno_SR_mask)
writeRaster(pam_gyno_SR_mask, filename="./output/pam_gyno_SR5_SR_mask.tif", overwrite=TRUE)

# convert raster to dataframe
# for ggplot2 plotting
pam_gyno_SR_mask_df <- pam_gyno_SR_mask %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)
### plot 

df1 <- pam_gyno_SR_mask_df # dataframe to use
Axis_0 <- "x" #x axis column
Axis_1 <- "y" #y axis column
sigplot <- "layer" # grid to plot column

yellowOrangeBrown_colours <- brewer.pal(9, "YlOrBr") # use colourbrewer to make colours
colours <- yellowOrangeBrown_colours
legend_text <- paste("Taxon Richness", sep="") # text for the legend
#legend_text <- paste0(sigplot)
legend_position <- paste("bottom", sep="") # position of the legend
rounding_digits <- 0 # rounding digits to use

#Create vectors for legend text using the max and min vaules in the raster and a number in between
legend_sequence.a <- seq(0,max(df1[,sigplot]),length.out=5)
legend_sequence.a.round <- round(legend_sequence.a, digits=rounding_digits)
legend_sequence.a.max_val <- round(max(df1[,sigplot]), digits=0)
legend_sequence.a.limits_set <- c(0,legend_sequence.a.max_val+0.1)

p1_new <- ggplot(data=df1) +  #xlim(min_x, max_x) +  ylim(min_y, max_y) +  
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x= Axis_0, y= Axis_1, fill=sigplot)) + 
  scale_fill_gradientn(name = NULL,#name = legend_text, 
                       limits = legend_sequence.a.limits_set, colours = colours, breaks= legend_sequence.a.round,
                       guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust=0.5,
                                               title.vjust=0.9, label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, raster=FALSE)) + 
  coord_fixed() +
  theme_pipeline() +
  #theme(legend.position=c(.6,0.01)) original code without plot the title, if I want to change back later
  theme(legend.position=c(.6,0.01), plot.title = element_text(size = 40, face = "bold", hjust = 0.04, vjust = 0)) +
  labs(title = "(A) Species richness (Gymnosperm)")

print(p1_new)
# print(p1)

if (print_seperate_images == TRUE){
  if (output_PNG == TRUE){
    CairoPNG(width = 2325, height = 2246, file = paste(output_folder, "figure_1_a.png", sep=""), canvas="white", bg = "white", units="px", dpi=72, title = "Figure 1 a") #
    print(p1_new)
    dev.off()
  }
  
  if (output_PDF == TRUE){
    CairoPDF(width = 36, height = 34, file = paste(output_folder, "figure_1_a.pdf",sep=""), pointsize=40, bg = "white", title = "Figure 1 a", version = "1.7", paper = "special", pagecentre=TRUE) #
    print(p1_new)
    dev.off()
  }
}

if (output_PDF == TRUE){
  pdf(file = paste(output_folder, "figure_1_a_new.pdf", sep=""), useDingbats=FALSE, width = 32, height = 16 )
  grid.arrange(p1_new) 
  dev.off()
}
########## WE  wighted species endemism  #####

# use rasterize to create desired raster
r_data_we <- rasterize(x=biodiverse_results_concatenated[, 2:3], # lon-lat data
                       y=r_obj, # raster object
                       field=biodiverse_results_concatenated$ENDW_WE) # vals to fill raster with
plot(r_data_we)
crs(r_data_we) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
pam_gyno_WE_mask <- mask(r_data_we,MAINL)
plot(pam_gyno_WE_mask)
writeRaster(pam_gyno_WE_mask, filename="./output/pam_gyno_SR5_WE_mask.tif", overwrite=TRUE)

# convert raster to dataframe
# for ggplot2 plotting
pam_gyno_WE_mask_df <- pam_gyno_WE_mask %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)
### plot 

df2 <- pam_gyno_WE_mask_df # dataframe to use
Axis_0 <- "x" #x axis column
Axis_1 <- "y" #y axis column
sigplot <- "layer" # grid to plot column

yellowOrangeBrown_colours <- brewer.pal(9, "YlOrBr") # use colourbrewer to make colours
colours <- yellowOrangeBrown_colours
legend_text <- paste("Weighted Endemism", sep="")# text for the legend
#legend_text <- paste0(sigplot)
legend_position <- paste("bottom", sep="") # position of the legend
rounding_digits <- 2 # rounding digits to use

#Create vectors for legend text using the max and min vaules in the raster and a number in between
legend_sequence.b <- seq(0,max(df2[,sigplot]),length.out=5)
legend_sequence.b.round <- round(legend_sequence.b, digits=rounding_digits)
legend_sequence.b.max_val <- round(max(df2[,sigplot]), digits=rounding_digits)
legend_sequence.b.limits_set <- c(0,legend_sequence.b.max_val)

p2_new <- ggplot(data=df2) +  #xlim(min_x, max_x) +  ylim(min_y, max_y) +  
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x= Axis_0, y= Axis_1, fill=sigplot)) + 
  scale_fill_gradientn(name = NULL, #name = legend_text, 
                       limits = legend_sequence.b.limits_set, colours = colours, breaks= legend_sequence.b.round,
                       guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust=0.5,
                                               title.vjust=0.9, label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, raster=FALSE)) + 
  coord_fixed() +
  theme_pipeline() +
  theme(legend.position=c(.6,0.01), 
        plot.title = element_text(size = 40, face = "bold", hjust = 0.04, vjust = 0)) +
  labs(title = "(B) Weighted endemism (Gymnosperm)")

print(p2_new)
# print
if (print_seperate_images == TRUE){
  if (output_PNG == TRUE){
    CairoPNG(width = 2325, height = 2246, file = paste(output_folder, "figure_1_b.png", sep=""), canvas="white", bg = "white", units="px", dpi=72, title = "Figure 1 b") #
    print(p2_new)
    dev.off()
  }
  
  if (output_PDF == TRUE){
    CairoPDF(width = 36.74, height = 39.19, file = paste(output_folder, "figure_1_b.pdf",sep=""), pointsize=40, bg = "white", title = "Figure 1 b", version = "1.7", paper = "special", pagecentre=TRUE) #
    print(p2_new)
    dev.off()
  }
}
if (output_PDF == TRUE){
  pdf(file = paste(output_folder, "figure_1_b_new_WE.pdf", sep=""), useDingbats=FALSE, width = 32, height = 16 )
  grid.arrange(p2_new) 
  dev.off()
}

########### PE,  phylogenetic endemism  ####
# use rasterize to create desired raster
r_data_pe <- rasterize(x=biodiverse_results_concatenated[, 2:3], # lon-lat data
                       y=r_obj, # raster object
                       field=biodiverse_results_concatenated$PE_WE_P) # vals to fill raster with
plot(r_data_pe)
crs(r_data_pe) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
pam_gyno_pe_mask <- mask(r_data_pe,MAINL)
plot(pam_gyno_pe_mask)
writeRaster(pam_gyno_pe_mask, filename="./output/pam_gyno_SR5_pe_mask.tif", overwrite=TRUE)
# convert raster to dataframe
# for ggplot2 plotting
pam_gyno_pe_mask_df <- pam_gyno_pe_mask %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)
### plot 

df4 <- pam_gyno_pe_mask_df # dataframe to use
Axis_0 <- "x" #x axis column
Axis_1 <- "y" #y axis column
sigplot <- "layer" # grid to plot column

legend_position <- paste("bottom", sep="") # position of the legend
rounding_digits <- 3 # rounding digits to use
yellowOrangeBrown_colours <- brewer.pal(9, "YlOrBr")
colours <- yellowOrangeBrown_colours

legend_sequence.d <- seq(0,max(df4[,sigplot]),length.out=5)
legend_sequence.d.round <- round(legend_sequence.d, digits=3)
legend_sequence.d.max_val <- round(max(df4[,sigplot]), digits=rounding_digits)
legend_sequence.d.limits_set <- c(0,legend_sequence.d.max_val)
#legend_text <- paste("Phylogenetic endemism", sep="")

p4_new <- ggplot(data=df4)+ 
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) + 
  scale_fill_gradientn(name = NULL, #name = legend_text,  not show the legend, as I add a title on top of the figure
                       limits = legend_sequence.d.limits_set, colours = colours, breaks= legend_sequence.d.round, 
                       guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust=0.5, title.vjust=0.9, 
                                               label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, raster=FALSE)) + 
  coord_fixed() +
  theme_pipeline()+
  theme(legend.position=c(.6,0.15), plot.title = element_text(size = 40, face = "bold", hjust = 0.04, vjust = 0)) + labs(title = "(D) Phylogenetic endemism")

print(p4_new)
# print
if (print_seperate_images == TRUE){
  if (output_PNG == TRUE){
    CairoPNG(width = 2000, height = 1000, file = paste(output_folder, "figure_1_d.png", sep=""), canvas="white", bg = "white", units="px", dpi=72, title = "Figure 1 b") #
    print(p4_new)
    dev.off()
  }
  
  ## pdf produced using below code can be edited
  if (output_PDF == TRUE){
    pdf(file = paste(output_folder, "figure_1_d.pdf", sep=""), useDingbats=FALSE, width = 32, height = 16 )
    grid.arrange(p4_new) 
    dev.off()
  }
}


## final Fig. 1b: gymnosperm PE zoomed in  ####
library(viridis)
library(dplyr)
p4_new1 <- ggplot(data=df4)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80,na.value='#f5f5f2')+
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(direction = 'vertical',  ## transform legend
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 30))

## Cory suggested zoom in certain regions

tmp_NA <- df4 %>% 
  filter(x< -4710.9375, x> -16716.7969, y<7141.3177, y>1079.2254) ## NA   1279.2254


gymno_pe_NA <-ggplot(tmp_NA)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(-16716.7969, -4710.9375),  ylim = c(1079.2254, 7141.3177), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0.2,0.2)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6)) + labs(fill = "PE") 

tmp_SA <- df4 %>% 
  filter(x< -3304.6875, x> -11830.0781, y<2736.7759, y>-6084.8972) ## SA 3036.7759


gymno_pe_SA <-ggplot(tmp_SA)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(-11830.0781, -3304.6875),  ylim = c(-6084.8972, 2536.7759), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0.1,0.5)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6))+ labs(fill = "PE")


tmp1 <- df4 %>%  
  filter(x< 6486.3281, x> -2003.9063, y<4027.4769, y>-4603.1332) ## Africa3977.4769


gymno_pe_AF <- ggplot(tmp1)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(-2003.9063, 6486.3281),  ylim = c(-4603.1332, 4027.4769), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0.1,0.2)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6))+ labs(fill = "PE")

tmp2 <- df4 %>% 
  filter(x< 6006.3281, x> -2003.9063, y< 7130.0691, y>3077.4769) ## Europe


gymno_pe_EU <-ggplot(tmp2)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74",size = 0.1, fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(-2003.9063, 6006.3281),  ylim = c(3077.4769, 7130.0691), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0,0.2)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6))+ labs(fill = "PE")

tmp_AS <- df4 %>% 
  filter(x< 16284.3750, x> 5986.3281, y<7360.4940, y> -1506.8851) ## Asia


gymno_pe_AS <-ggplot(tmp_AS)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(5986.3281, 16284.3750),  ylim = c(-1506.8851, 7360.4940), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0.8,0.4)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6))+ labs(fill = "PE")

tmp_AU <- df4 %>% 
  filter(x< 17610.9375, x> 10603.9063, y< -1006.8851, y> -5606.8851) ## Australia -5706.8851


gymno_pe_AU <-ggplot(tmp_AU)+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) +
  scale_fill_viridis(alpha=0.80)+
  coord_fixed(xlim = c(10603.9063, 17610.9375),  ylim = c(-5606.8851, -1006.8851), expand = c(0,0)) +
  theme(panel.background = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#f5f5f2", fill=NA, size=1),
        legend.justification=c(0,0), legend.position=c(0.7,0.2)) +
  guides(fill = guide_colourbar(direction = 'vertical',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='#f5f5f2',
                                ticks.linewidth=2,
                                barwidth = 0.5,
                                barheight = 6))+ labs(fill = "PE")


library(cowplot)

# Arranging multiple plots in a figure
gymno_pe_full <- plot_grid(gymno_pe_NA,gymno_pe_EU,gymno_pe_AS,
                           gymno_pe_SA,gymno_pe_AF,gymno_pe_AU,
                           rel_widths = c(1, 1,1),rel_heights = c(1, 1,1),
                           nrow = 2)

pdf("fig_gymno_pe.pdf", useDingbats=FALSE, width=20, height=16)
gymno_pe_full
dev.off()



################### Fig. 2b CANEPE 4 groups#############
# use rasterize to create desired raster  P_PHYLO_RPE2_CANAPE_SIG
##prepare a df for coordinates and biodiverse_results_concatenated$P_PHYLO_RPE2_CANAPE_SIG
canpe_gymno <- biodiverse_results_concatenated[, c("Axis_0","Axis_1", "P_PHYLO_RPE2_CANAPE_SIG")]
### add a new column for the 
library(dplyr)
canpe_gymno1 <- canpe_gymno %>% dplyr::mutate(Group =
                                                case_when(P_PHYLO_RPE2_CANAPE_SIG == "Neo" ~ 1, 
                                                          P_PHYLO_RPE2_CANAPE_SIG == "Paleo" ~ 2,
                                                          P_PHYLO_RPE2_CANAPE_SIG == "Not Sig"  ~ 3,
                                                          P_PHYLO_RPE2_CANAPE_SIG == "Mixed" ~ 4)
)

r_data_canpe_p <- rasterize(x=canpe_gymno1[, 1:2], # lon-lat data
                            y=r_obj, # raster object
                            field=canpe_gymno1$Group) # vals to fill raster with
plot(r_data_canpe_p)
crs(r_data_canpe_p) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
r_data_canpe_p_mask <- mask(r_data_canpe_p,MAINL)
plot(r_data_canpe_p_mask)
writeRaster(r_data_canpe_p_mask, filename="./output/gymno_SR5_canpe_p_mask_4_CANAPE_groups.tif", overwrite=TRUE)

# convert raster to data frame
# for ggplot2 plotting
r_data_canpe_p_mask_df <- r_data_canpe_p_mask %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)
## add a new variable
r_data_canpe_p_mask_df_df1 <- r_data_canpe_p_mask_df %>% 
  dplyr::mutate(P_PHYLO_RPE2_CANAPE_SIG = layer)

## pick the variable
legend_order_canpo <- r_data_canpe_p_mask_df_df1 %>% select(P_PHYLO_RPE2_CANAPE_SIG) #pick the variable 
library(dplyr)
## tags
labelcapno <-c("Neo","Paleo", "Not Sig", "Mixed")

legend_order_labels_canpo <- as_tibble(legend_order_canpo) %>% 
  dplyr::mutate(tag = case_when(
    P_PHYLO_RPE2_CANAPE_SIG == 1 ~ labelcapno[1],
    P_PHYLO_RPE2_CANAPE_SIG == 2  ~ labelcapno[2],
    P_PHYLO_RPE2_CANAPE_SIG == 3  ~ labelcapno[3],
    P_PHYLO_RPE2_CANAPE_SIG == 4 ~ labelcapno[4]
  ))
summary(as.factor(legend_order_labels_canpo$tag))
dat_canpo <- as.data.frame(legend_order_labels_canpo)
dat_canpo$tag <-as.factor(dat_canpo$tag)
r_data_canpe_p_mask_df_final <- cbind(r_data_canpe_p_mask_df, dat_canpo)

map_text <- "Categorical Analysis of Neo- And Paleo- Endemism"
sigplot <- "tag"

col_scheme1 <- c("Not Sig" = "lightgoldenrodyellow","Neo" = "red","Mixed" = "#33A02C", "Paleo" = "purple") #,"Super" = "#9D00FF"

legend_order <-c("Neo","Paleo", "Not Sig", "Mixed") #, "Super
legend_labels <- c("Neo"="Neo","Paleo"="Paleo", "Not Sig"="Not significant", "Mixed"="Mixed") #, "Super"="Super"

r_data_canpe_p_mask_df_final[, "sigplot"] <- factor(r_data_canpe_p_mask_df_final[, sigplot], levels=legend_order)
Axis_0 <- "x"
Axis_1 <- "y"

map_plot_5 <- ggplot(data=r_data_canpe_p_mask_df_final) +
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot))+ # using aes_string allows variables to be passed to ggplot
  #scale_fill_manual(values = col_scheme) +  
  scale_fill_manual(values = col_scheme1, labels=legend_labels, name=map_text, 
                    guide = guide_legend(direction = "horizontal", title.position = "top", title.hjust=0.5, 
                                         label.position="bottom", label.hjust = 0.5, label.vjust = 0.5))+ #label.theme = element_text(angle = 90), label.hjust = 0.5, label.vjust = 0.5
  theme_minimal() + 
  coord_fixed() +
  theme_pipeline()+
  theme(legend.title = element_text(colour = 'black', angle = 0, size=rel(3.8), family = myFont),
        legend.text = element_text(colour = 'black', angle = 0, size=rel(3.8), family = myFont))+
  #theme(legend.position=c(.6,0.03))
  theme(legend.position=c(.6,0.15), plot.title = element_text(size = 40, face = "bold", hjust = 0.04, vjust = 0)) #+ labs(title = "Categorical Analysis of Neo- And Paleo- Endemism")

print(map_plot_5)


if (print_seperate_images == TRUE){
  if (output_PNG == TRUE){
    CairoPNG(width = 2000, height = 1000, file = paste(output_folder, "figure_3a_4groups_16Feb2023.png", sep=""), canvas="white", bg = "white", units="px", dpi=72, title = "Figure 1 b") #
    print(map_plot_5)
    dev.off()
  }
  
  ## pdf produced using below code can be edited
  if (output_PDF == TRUE){
    pdf(file = paste(output_folder, "figure_3a_4gruops_16Feb2023.pdf", sep=""), useDingbats=FALSE,  width = 32, height = 16  )
    grid.arrange(map_plot_5) 
    dev.off()
  }
}

########### PE without coastline masked, Fig. S3b ####
# use rasterize to create desired raster
r_data_pe <- rasterize(x=biodiverse_results_concatenated[, 2:3], # lon-lat data
                       y=r_obj, # raster object
                       field=biodiverse_results_concatenated$PE_WE_P) # vals to fill raster with
plot(r_data_pe)
crs(r_data_pe) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
pam_gyno_pe_mask <- mask(r_data_pe,MAINL)
plot(pam_gyno_pe_mask)
writeRaster(r_data_pe, filename="./output/pam_gyno_SR5_pe_noCoastline_mask.tif", overwrite=TRUE)
# convert raster to dataframe
# for ggplot2 plotting
pam_gyno_pe_df <- r_data_pe %>%
  rasterToPoints %>%
  as.data.frame() %>%
  filter(layer!=0)
### plot 

df4_nc <- pam_gyno_pe_df # dataframe to use
Axis_0 <- "x" #x axis column
Axis_1 <- "y" #y axis column
sigplot <- "layer" # grid to plot column

legend_position <- paste("bottom", sep="") # position of the legend
rounding_digits <- 3 # rounding digits to use
yellowOrangeBrown_colours <- brewer.pal(9, "YlOrBr")
colours <- yellowOrangeBrown_colours

legend_sequence.d_nc <- seq(0,max(df4_nc[,sigplot]),length.out=5)
legend_sequence.d.round_nc <- round(legend_sequence.d_nc, digits=3)
legend_sequence.d.max_val_nc <- round(max(df4_nc[,sigplot]), digits=rounding_digits)
legend_sequence.d.limits_set_nc <- c(0,legend_sequence.d.max_val_nc)
#legend_text <- paste("Phylogenetic endemism", sep="")

p4_new_nc <- ggplot(data=df4_nc)+ 
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="gray74", fill="white") +
  geom_tile(aes_string(x=Axis_0, y=Axis_1, fill=sigplot)) + 
  scale_fill_gradientn(name = NULL, #name = legend_text,  not show the legend, as I add a title on top of the figure
                       limits = legend_sequence.d.limits_set_nc, colours = colours, breaks= legend_sequence.d.round_nc, 
                       guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust=0.5, title.vjust=0.9, 
                                               label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, raster=FALSE)) + 
  coord_fixed() +
  theme_pipeline()+
  theme(legend.position=c(.6,0.15), plot.title = element_text(size = 40, face = "bold", hjust = 0.04, vjust = 0)) + labs(title = "(D) Phylogenetic endemism")

print(p4_new_nc)
# print
if (print_seperate_images == TRUE){
  if (output_PNG == TRUE){
    CairoPNG(width = 2000, height = 1000, file = paste(output_folder, "figure_pe_noCoastline_masked.png", sep=""), canvas="white", bg = "white", units="px", dpi=72, title = "Figure 1 b") #
    print(p4_new_nc)
    dev.off()
  }
  
  ## pdf produced using below code can be edited
  if (output_PDF == TRUE){
    pdf(file = paste(output_folder, "figure_pe_noCoastline_masked.pdf", sep=""), useDingbats=FALSE, width = 32, height = 16 )
    grid.arrange(p4_new_nc) 
    dev.off()
  }
}
