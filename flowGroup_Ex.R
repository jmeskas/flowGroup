# Written by Justin Meskas
# Last updated on September 2019
# alpha version 0.99.5

library("flowCore")
library("flowDensity")
library("foreach")

library("Cairo")
library("MASS") # kde2d
library("doMC")
library("stringr") # str_pad

library("OpenImageR")
library("flowType")
library("Rtsne")

# library("raster")
library("oce")
library("plot3D")

source("~/code/flowGroup/flowGroup.R")
source("~/code/flowGroup/flowGroup_functions.R")

# need to input
path.output <- ""

suppressWarnings ( dir.create ( paste0(path.output,"/flowGroup")))

# get names of files that have CD4 data saved.
all.cd4.files <- list.files(paste0(path.output, "/flowGroup_data/CD4"),full.names = T)
# load list of flowFrames
cd4.files <- lapply(all.cd4.files, function(x){load(file = x); list(ff.cd4) })
# morph into flowSet
cd4.files <- as(unlist(cd4.files), "flowSet")

# need to load the gating thresholds to incorporate in plot
gates1 <- 
gates2 <- 

labels <- sapply(all.cd4.files, function(x){
    temp <- strsplit(x, split = "/CD4/")[[1]][2]
    temp <- strsplit(temp, split = "__")[[1]][1]
})
names(labels) <- NULL


res.flowGroup.cd4 <- flowGroup(files=cd4.files, xMark=14, yMark=9, experiment_name=paste0(path.output, "flowGroup/CD4"),
                      HOG_cex=1, partitions_flowType = 5, vuc=0, vuc_d = 0, hog=0, flwT = 1, DoOverlay=T,
	              plot_groups = T, verbose=T, xlim=c(-0.5,4.5), ylim=c(-0.5,4.5), vert_line=gates1, horiz_line=gates2, k2run = 8, Labels = labels)




