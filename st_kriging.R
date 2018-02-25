library(gstat)
library(spacetime)
library(sp)
library(raster)
library(rgdal)
library(rgeos) 
library(maptools)
library(dplyr)
library(reshape2)
library(alphahull)
library(lubridate)
library(stringr)
library(parallel)
library(microbenchmark)


setwd("~/Documents/Projects/efn_2017/rstats/krige")
source("st_kriging_functions.R") # Read functions held in a seperate script

# Read & prepare the data ------------------------------------------------------

# Read the filenames in a folder for a given study area:
filenames <- list.files("/Users/sab/Dropbox/EFN\ Phase\ 2/Channel\ centric\ projection/SA2\ dimensionless/", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv) # Creat a list from the data files with each element representing a survey year
res <- lapply(ldf, summary) # Summerize each year of data
names(res) <- substr(filenames, 6, 30)

survey_ldf <- lapply(ldf,filter,hub != "RBH" & hub != "LBH") # Filter the survey points
hubs_ldf <- lapply(ldf,filter,hub == "RBH" | hub == "LBH") # Filter the hubs

# There may be duplicate shots in a given survey year of data (especially where
# two CS's cross). These need to be removed using the function
# "remove_duplicates" before a variogram can be built or the data can be kriged:
remove_duplicates <- function(survey_ldf) {
        coordinates(survey_ldf) <- ~ Easting + Northing # Generate a SpatialPoints data frame within each list element
        survey_ldf <- remove.duplicates(survey_ldf, zero = 0.0, remove.second = TRUE, memcmp = TRUE) # Remove duplicates in a given element
        survey_ldf <- as.data.frame(survey_ldf) # Convert the SpatialPoints data frame back to a regular df
        survey_ldf # Return the list with duplicates removed
}

survey_ldf <- lapply(survey_ldf,remove_duplicates) # Return the list with duplicates removed
survey_data <- do.call(rbind.data.frame, survey_ldf) # Create a df from the list

# Subset all the active bed points (include BR, BL, and EoV):
feature_logic1 <- str_detect(survey_data$attribs,"BL")
feature_logic2 <- str_detect(survey_data$attribs,"BR")
feature_logic3 <- str_detect(survey_data$attribs,"EOV")
survey_data$feature_logic <-
        ifelse(
                feature_logic1 == TRUE |
                        feature_logic2 == TRUE |
                        feature_logic3 == TRUE |
                        survey_data$attribs == "" |
                        survey_data$attribs == "t",
                "keep",
                "delete"
        )

bed_pts <- filter(survey_data,feature_logic == "keep")
bed_pts$feature_logic <- NULL

# ST variogram does not work with "year" time-steps, so the next three lines
# create a fake date/time using "seconds" time-step (e.g. 1971 becomes
# 1971-01-01 00:00:01, while 1981 becomes 1971-01-01 00:00:10, etc.)
bed_pts$Year <- bed_pts$Year - 1970
bed_pts$time <- ifelse(bed_pts$Year < 10, paste("1971-01-01 0:0:0", bed_pts$Year, sep=""), paste("1971-01-01 0:0:", bed_pts$Year, sep=""))
bed_pts$time <- ymd_hms(bed_pts$time)

bed_pts_all <- bed_pts
bed_pts <- sample_n(bed_pts,1500)
time <- bed_pts$time
elev <- data.frame(bed_pts$Elevation) # Elevation needs to be in a single column df

# Create a spatial points df
data_coords <- bed_pts[,5:6]
coordinates(data_coords) <- ~ Easting + Northing
#data_coords <- elide(data_coords, shift=c(-353000, -5419000)) # Translate the shapefile since we need more numbers than R will allow in a floating point
proj4string(data_coords)=CRS("+init=epsg:26910") # Set the EPSG code for NAD83 / UTM zone 10N
timeDF <- STIDF(data_coords,time,elev) # Spatial-temporal data frame
stplot(timeDF)



time <- bed_pts_all$time
elev <- data.frame(bed_pts_all$Elevation) # Elevation needs to be in a single column df

# Create a spatial points df
data_coords <- bed_pts_all[,5:6]
coordinates(data_coords) <- ~ Easting + Northing
#data_coords <- elide(data_coords, shift=c(-353000, -5419000)) # Translate the shapefile since we need more numbers than R will allow in a floating point
proj4string(data_coords)=CRS("+init=epsg:26910") # Set the EPSG code for NAD83 / UTM zone 10N
timeDF_all <- STIDF(data_coords,time,elev) # Spatial-temporal data frame
stplot(timeDF_all)








# Fit a variogram --------------------------------------------------------------

varg <- variogramST(bed_pts.Elevation~1,data=timeDF,tunit="secs",assumeRegular=F,na.omit=T,tlags=0:15,progress=T)
# res <- microbenchmark(variogramST(bed_pts.Elevation~1,data=timeDF,tunit="secs",assumeRegular=F,na.omit=T,tlags=0:15,progress=T),times=1)
saveRDS(varg, "~/Documents/Projects/efn_2017/rstats/krige/sa2_st_variogram.rds")
varg <- readRDS("~/Documents/Projects/efn_2017/rstats/krige/sa2_st_variogram.rds")

# Evaluate the distribution:
plot(varg,map=F)
plot(varg,map=T)
plot(varg, wireframe=T, scales=list(arrows=F))

# Estimate the low-end of the parameters (for the optim function)
pars.l <- c(sill.s = .06, range.s = .75, nugget.s = 0.0,
            sill.t = 0.02, range.t = 6, nugget.t = 0,
            sill.st = 0.1, range.st = 7.5, nugget.st = 0, anis = 170)
# Estimate the upper-end of the parameters (for the optim function)
pars.u <- c(sill.s = .11, range.s = 1.2, nugget.s = 0.01,
            sill.t = 0.04, range.t = 12, nugget.t = .001,
            sill.st = .25, range.st = 15, nugget.st = 0.1,anis = 190)


# Fit different models to the variogram and choose the one with the best fit (i.e. the lowest MSE):

# Some models first require seperate space and time vgm's"
space_vgm <- vgm(psill = 0.2, "Exp", range = 15, anis = c(180,1), Err = 0.02^2)
time_vgm <- vgm(psill = 0.1, "Exp", range = 10, Err = 0.25^2)

# Separable model
separable <- vgmST("separable", space = space_vgm,time = time_vgm, sill=0.2)
plot(varg,separable,map=F)
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=0)
attr(separable_Vgm,"MSE")
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
attr(separable_Vgm, "MSE")
plot(varg,separable_Vgm,map=F) 
extractPar(separable_Vgm)

# Product-sum model
prodSumModel <- vgmST("productSum",space = space_vgm,time = time_vgm,k = 5)
plot(varg,prodSumModel,map=F)
prodSumModel_Vgm <- fit.StVariogram(varg, prodSumModel,method = "L-BFGS-B",lower=pars.l)
attr(prodSumModel_Vgm, "MSE")
plot(varg,prodSumModel_Vgm,map=F) 
extractPar(prodSumModel_Vgm)

# Metric model
metric <- vgmST("metric", joint = vgm(0.15,"Sph",12,0), stAni=.1)
plot(varg,metric,map=F)
metric_Vgm <- fit.StVariogram(varg, metric, method="L-BFGS-B",lower=pars.l,upper=pars.u)
attr(metric_Vgm, "MSE")
plot(varg,metric_Vgm,map=F) 
extractPar(metric_Vgm)

# Sum metric
sumMetric <- vgmST("sumMetric", space = space_vgm,time = time_vgm, joint = vgm(psill=0.2,"Sph", range=12, nugget=0), stAni=0.1) 
plot(varg,sumMetric,map=F)
sumMetric_Vgm <- fit.StVariogram(varg, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="secs")
attr(sumMetric_Vgm, "MSE")
plot(varg,sumMetric_Vgm,map=F) 
extractPar(sumMetric_Vgm)

# Simple sum metric
SimplesumMetric <- vgmST("simpleSumMetric",space = space_vgm,time = time_vgm, joint = vgm(0.5,"Sph", 12, 0), nugget=0, stAni=0.1) 
plot(varg,SimplesumMetric,map=F)
SimplesumMetric_Vgm <- fit.StVariogram(varg, SimplesumMetric,method = "L-BFGS-B",lower=pars.l)
attr(SimplesumMetric_Vgm, "MSE")
plot(varg,SimplesumMetric_Vgm,map=F) 
extractPar(SimplesumMetric_Vgm)

# Plot all the models
plot(varg,list(separable_Vgm, prodSumModel_Vgm, metric_Vgm, sumMetric_Vgm, SimplesumMetric_Vgm),all=T,wireframe=T) 

vgm_mse <- c(attr(separable_Vgm, "MSE"),attr(prodSumModel_Vgm, "MSE"),attr(metric_Vgm, "MSE"),attr(sumMetric_Vgm, "MSE"),attr(SimplesumMetric_Vgm, "MSE"))
models <- c("seperable","product-sum","metric","sum-metric","simple sum-metric")
vgm_results <- data.frame(models,vgm_mse)

# Choose the model with the lowest MSE:
print(vgm_results)
m <- metric_Vgm # User needs to set the correct model

# Create the temporal and spatial prediction grid ------------------------------

points <- survey_data[,5:7] # Subset the easting, northing, and elevation
coordinates(points) <- ~ Easting + Northing
points <- remove.duplicates(points, zero = 0.0, remove.second = TRUE, memcmp = TRUE)
points <- as.data.frame(points)

# Define a function that creates the smallest rectangle that includes all data
# points otherwise known as the minimum bounding region (MBR):

MBR <- function(points) {
        # Analyze the convex hull edges                       
        a <- ashape(points, alpha=1000) # One way to get a convex hull...
        e <- a$edges[, 5:6] - a$edges[, 3:4] # Edge directions
        norms <- apply(e, 1, function(x) sqrt(x %*% x)) # Edge lengths
        v <- diag(1/norms) %*% e # Unit edge directions
        w <- cbind(-v[,2], v[,1]) # Normal directions to the edges
        
        # Find the MBR
        vertices <- (points) [a$alpha.extremes, 1:2] # Convex hull vertices
        minmax <- function(x) c(min(x), max(x)) # Computes min and max
        x <- apply(vertices %*% t(v), 2, minmax) # Extremes along edges
        y <- apply(vertices %*% t(w), 2, minmax) # Extremes normal to edges
        areas <- (y[1,]-y[2,])*(x[1,]-x[2,]) # Areas
        k <- which.min(areas) # Index of the best edge (smallest area)
        
        # Form a rectangle from the extremes of the best edge
        cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])
}

points <- as.matrix(points[,1:2])
mbr <- MBR(points) # The minimum bounding region

# Plot the hull, MBR, and points. Note that the hull gives the boundary
# where no extrapolation occurs.
limits <- apply(mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
plot(ashape(points, alpha=1000), col="Gray", pch=20, 
     xlim=limits[,1], ylim=limits[,2],asp=1) # The hull
lines(mbr, col="Blue", lwd=3) # The MBR
points(points, pch=19) # The points
print(mbr) # Point listing

# Plot the MBR with the survey data
x1 <- mbr[1:4,1]
y1 <- mbr[1:4,2]
c1 <- cbind(x1, y1)
r1 <- rbind(c1, c1[1, ])  # join
P1 <- Polygon(r1)
Ps1 <- Polygons(list(P1), ID = "a")
# Create a spatial polygons df:
border <- SpatialPolygons(list(Ps1))
#border <- elide(border, shift=c(-353000, -5419000)) # Translate the shapefile since we need more numbers than R will allow in a floating point 
proj4string(border)=CRS("+init=epsg:26910") # Adjust the epsg code as necassary
plot(border)
plot(data_coords,add=T)

# Create a spatial prediction grid using the "border" shape file:
vals <- border@bbox
delta.x <- as.integer((vals[1,2] - vals[1,1]) + 1.5)
delta.y <- as.integer((vals[2,2] - vals[2,1]) + 1.5)
grid.res <- .7/15   # Change this value to change the grid size (in metres)
grid.size.x <- delta.x / grid.res
grid.size.y <- delta.y / grid.res
grd <- GridTopology(vals[,1],c(grid.res,grid.res),c(grid.size.x,grid.size.y))
pts <- SpatialPoints(coordinates(grd))
proj4string(pts)=CRS("+init=epsg:26910")  # Remember to set the coordinate system also for the prediction grid

# Create a temporal prediction grid:
tm.grid <- seq(ymd_hms('1971-01-01 00:00:01'),ymd_hms('1971-01-01 00:00:03'),length.out=2)

# Create a spatail-temporal prediction grid: 
grid.ST <- STF(pts,tm.grid)


# Start Cluster -----------------------------------------------------------

# This code block was modified from two sources:
# https://gis.stackexchange.com/questions/237672/how-to-achieve-parallel-kriging-in-r-to-speed-up-the-process/237940
# and,
# https://github.com/guzmanlopez/parallelizingAndClusteringInR/blob/master/parallel-kriging-function.R
# The former is straight-forward and adaopted here for simplicity. The latter is
# probably better (faster) and allows for clusters (i.e. running on multiple computers). The only modification
# I made was in the subsetting and rebuilding of the STFDF. The orginal code was
# written for spatial kriging only.

no_cores <- detectCores() - 1 # Calculate the number of cores
cl <- makeCluster(no_cores) # Initiate cluster
parts <- split(x = 1:9600, f = 1:no_cores)
clusterExport(cl = cl, varlist = c("timeDF", "grid.ST", "m", "parts"), envir = .GlobalEnv)
clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
# Set-up the kriging function:
parallelX <- parLapply(cl = cl, 
                       X = 1:no_cores, 
                       fun = function(x) krigeST(formula = bed_pts.Elevation~1,data=timeDF, modelList=m, newdata=grid.ST[c(parts[[x]]),],computeVar=TRUE)
)

# "parallelX" is a list of "n" STFDF's (where "n" = number of cores). These
# need to be merged into a single STFDF. I can't find a function to merge two
# STFDF's. Instead, the following lines of code convert all the SDFDF clusters
# to a single data frame, and then a final STFDF is built from scratch:

clusters_df <- lapply(parallelX,as.data.frame) # Convert the STFDF's created for each cluster to a list of data frames
clusters_df <- do.call("rbind", clusters_df) # Convert the list of data frames into one data frame
clusters_df <- coordinates(clusters_df[,1:2])

length(clusters_df)

ti1 <- filter(clusters_df, timeIndex == 1)
ti2 <- filter(clusters_df, timeIndex == 2)

ti12 <- cbind(ti1,ti2)
ti12 <- ti12[,c(1,2,4,12,7)]




clusters_df <- arrange(clusters_df,timeIndex) # Sort by time
kriged_coord <- clusters_df[,1:2]
coordinates(kriged_coord) <- ~ x + y
proj4string(kriged_coord)=CRS("+init=epsg:26910")
kriged_time <- clusters_df[1,4]
kriged_time2 <- clusters_df[6400,4]
kriged_time3 <- c(kriged_time,kriged_time2)
kriged_elev <- data.frame(clusters_df[,7])
kriged_data <- STFDF(kriged_coord,kriged_time,kriged_elev)
stplot(kriged_data)





filter(ti1, sp.ID==2)

time2 <- 




xcx = stConstruct(clusters_df, c("x", "y"), "time", interval = FALSE)
class(x)
stplot(x[,,"unemp"])

stConstruct(Produc, c("x", "y"), "time", interval = TRUE)

stConstruct(kriged_elev,)

data_coords <- bed_pts_all[,5:6]
coordinates(data_coords) <- ~ Easting + Northing
#data_coords <- elide(data_coords, shift=c(-353000, -5419000)) # Translate the shapefile since we need more numbers than R will allow in a floating point
proj4string(data_coords)=CRS("+init=epsg:26910") # Set the EPSG code for NAD83 / UTM zone 10N
timeDF_all <- STIDF(data_coords,time,elev) # Spatial-temporal data frame
stplot(timeDF_all)

a <- remove.duplicates(data_coords, zero = 0.0, remove.second = TRUE, memcmp = TRUE)














































################################################################################



# lower and upper bounds
#pars.l <- c(sill.s = .2, range.s = 5, nugget.s = 0.0,sill.t = 0.05, range.t = 6, nugget.t = 0,sill.st = 0.21, range.st = 7.5, nugget.st = 0, anis = 290)
#pars.u <- c(sill.s = .3, range.s = 20, nugget.s = 0.1,sill.t = 0.15, range.t = 12, nugget.t = .001,sill.st = .35, range.st = 24, nugget.st = 0.1,anis = 330) 
pars.l <- c(sill.s = .2, range.s = 5, nugget.s = 0.0,sill.t = 0.05, range.t = 6, nugget.t = 0,sill.st = 0.21, range.st = 7.5, nugget.st = 0, anis = 170)
pars.u <- c(sill.s = .3, range.s = 20, nugget.s = 0.1,sill.t = 0.15, range.t = 12, nugget.t = .001,sill.st = .35, range.st = 24, nugget.st = 0.1,anis = 190) 
space_vgm <- vgm(psill = 0.5, "Exp", range = 15, anis = c(180,1), Err = 0.02^2)
time_vgm <- vgm(psill = 0.1, "Exp", range = 10, Err = 0.25^2)

# Separable model
separable <- vgmST("separable", space = space_vgm,time = time_vgm, sill=0.5)
plot(varg,separable,map=F)
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=0)
attr(separable_Vgm,"MSE")
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
attr(separable_Vgm, "MSE")
plot(varg,separable_Vgm,map=F) 
extractPar(separable_Vgm)

# Product-sum model
prodSumModel <- vgmST("productSum",space = space_vgm,time = time_vgm,k = 5)
plot(varg,prodSumModel,map=F)
prodSumModel_Vgm <- fit.StVariogram(varg, prodSumModel,method = "L-BFGS-B",lower=pars.l)
attr(prodSumModel_Vgm, "MSE")
plot(varg,prodSumModel_Vgm,map=F) 
extractPar(prodSumModel_Vgm)

# Metric model
metric <- vgmST("metric", joint = vgm(0.15,"Exp",12,0), stAni=.1)
plot(varg,metric,map=F)
metric_Vgm <- fit.StVariogram(varg, metric, method="L-BFGS-B",lower=pars.l,upper=pars.u)
attr(metric_Vgm, "MSE")
plot(varg,metric_Vgm,map=F) 
extractPar(metric_Vgm)

# Sum metric
sumMetric <- vgmST("sumMetric", space = space_vgm,time = time_vgm, joint = vgm(psill=0.15,"Sph", range=12, nugget=0), stAni=1) 
plot(varg,sumMetric,map=F)
sumMetric_Vgm <- fit.StVariogram(varg, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="secs")
attr(sumMetric_Vgm, "MSE")
plot(varg,sumMetric_Vgm,map=F) 
extractPar(sumMetric_Vgm)

# Simple sum metric
SimplesumMetric <- vgmST("simpleSumMetric",space = space_vgm,time = time_vgm, joint = vgm(0.5,"Sph", 15, 0.1), nugget=0.1, stAni=1) 
plot(varg,SimplesumMetric,map=F)
SimplesumMetric_Vgm <- fit.StVariogram(varg, SimplesumMetric,method = "L-BFGS-B",lower=pars.l)
attr(SimplesumMetric_Vgm, "MSE")
plot(varg,SimplesumMetric_Vgm,map=F) 
extractPar(SimplesumMetric_Vgm)

plot(varg,list(separable_Vgm, prodSumModel_Vgm, metric_Vgm, SimplesumMetric_Vgm),all=T,wireframe=T) 














# Calculate temporal autocorrelation and cross correlation ---------------------

# Find coincident points with at least one spatial overlap. Then produce a STI
# irregular layout with pts that repeat at least once:
active_channel <- filter(survey_data, attribs == "" | attribs == "t")
active_channel$Easting <- round(active_channel$Easting,3) 
active_channel$Northing <- round(active_channel$Northing,3) 
coordinates(active_channel) <- c("Easting", "Northing") #Set the Easting and Northing to coordinates
active_channel$unique_id <- zerodist(active_channel, zero = 0.0, unique.ID = TRUE, memcmp = TRUE) # Assign unique ID
active_channel <- as.data.frame(active_channel) # Back to a df
active_channel <- arrange(active_channel,unique_id) # Arrange the df by unique ID
active_channel$dup <- duplicated(active_channel$unique_id) # "TRUE" flags a repeated ID no. (first pt in seq. flagged as "FALSE")
repeat_pts <- data.frame(active_channel$unique_id,active_channel$dup)
names(repeat_pts) <- c("unique_id","dup")
repeat_pts <- filter(repeat_pts,dup == "TRUE") # Filter the repeated points
repeat_pts <- unique(repeat_pts) # All the unique ID's that repeat at least once
repeat_pts$dup <- NULL
repeat_pts <- left_join(repeat_pts,active_channel) # df with all repeating pts (even the first pt in a repeating sequence) 
repeat_pts <- arrange(repeat_pts,time)
no_repeats <- length(repeat_pts$unique_id) # Total number of repeat observations in the study area
# For testing, it helps to take a small subset of the data
set.seed(42)
repeatsub_pts <- sample_n(repeat_pts, 1500)

# Create a spatial points df
repeat_coords <- repeat_pts[,6:7]
coordinates(repeat_coords) <- ~ Easting + Northing
repeat_coords <- elide(repeat_coords, shift=c(-353000, -5419000)) #
proj4string(repeat_coords)=CRS("+init=epsg:26910") #Set the EPSG code for NAD83 / UTM zone 10N
repeat_elev <- data.frame(repeat_pts$Elevation)
repeat_time <- repeat_pts$time
mydata <- repeat_elev
repeat_timeDF <- STIDF(repeat_coords,repeat_time,repeat_elev)
stplot(repeat_timeDF) 



repeat_df <- as.data.frame(repeat_timeDF)
repeat_locations <- unique(repeat_df[c("x","y")]) # Create a factor for each unique x,y point
repeat_locations$unique_pt <- as.factor(seq_along(repeat_locations$x))
repeat_df <- left_join(repeat_df,repeat_locations)
repeat_df <- repeat_df %>% add_count(unique_pt)
arrange(repeat_df,desc(n))
hist(repeat_df$n)

big_repeats <- as.data.frame(filter(repeat_df, n >= 7))

repeat_df <- as.data.frame(repeat_df)
#coordinates(repeat_df) <- ~ x + y
#proj4string(repeat_df)=CRS("+init=epsg:26910") #Set the EPSG 
#spplot(repeat_df,c("n"))



tot <- tabulate(repeat_df$n)
tot_seq <- seq_along(tot)
tot/tot_seq

repeat_df7 <- filter(repeat_df,n == 7)
repeat_df7 <- arrange(repeat_df7,unique_pt)
repeat_df7_469 <- filter(repeat_df7,unique_pt == 469)


date2 <- repeat_df7_469$time #extract the date and time (in this example, the date and time are in a single column from Excel)
date3 <- strptime(date2, format = "%Y-%m-%d %H:%M:%S", tz = "EST") #convert the date into time format
q1 <- zoo(repeat_df7_469$repeat_pts.Elevation, repeat_df7_469$time) #use "Zoo" package to turn the data into an time object

acf(na.omit(as(r5to10[4:10,], "xts")))

acf(na.omit(q1))

# Fit a variogram --------------------------------------------------------------

varg <- variogramST(repeat_pts.Elevation~1,data=timeDF,tunit="secs",assumeRegular=F,na.omit=T)
var3 <- variogramST(survey_data.Elevation~1,data=timeDF_all,tunit="secs",assumeRegular=T,na.omit=T)
saveRDS(varg, "~/Documents/Projects/efn_2017/rstats/krige/varsa2.rds")
saveRDS(var2, "~/Documents/Projects/efn_2017/rstats/krige/var2sa2.rds")
varg <- readRDS("~/Documents/Projects/efn_2017/rstats/krige/varsa2.rds")


varg <- var3

plot(var3,map=F)
plot(var3,map=T)
plot(var3,wireframe=T)

# lower and upper bounds
pars.l <- c(sill.s = 0.4, range.s = 5, nugget.s = 0.1,sill.t = 0, range.t = 1, nugget.t = 0,sill.st = 0, range.st = 5, nugget.st = 0, anis = 0)
pars.u <- c(sill.s = 0.6, range.s = 25, nugget.s = 0.3,sill.t = 30, range.t = 30, nugget.t = 1,sill.st = 1, range.st = 30, nugget.st = 1,anis = 0) 
space_vgm <- vgm(psill = 0.5, "Lin", range = 15, 0.15)
time_vgm <- vgm(psill = 0.5, "Lin", range = 15, 0.15)

# Separable model
separable <- vgmST("separable", space = space_vgm,time = time_vgm, sill=0.5)
plot(varg,separable,map=F)
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=0)
attr(separable_Vgm,"MSE")
separable_Vgm <- fit.StVariogram(varg, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
attr(separable_Vgm, "MSE")
plot(varg,separable_Vgm,map=F) 
extractPar(separable_Vgm)

# Product-sum model
prodSumModel <- vgmST("productSum",space = space_vgm,time = time_vgm,k = 5)
plot(varg,prodSumModel,map=F)
prodSumModel_Vgm <- fit.StVariogram(varg, prodSumModel,method = "L-BFGS-B",lower=pars.l)
attr(prodSumModel_Vgm, "MSE")
plot(varg,prodSumModel_Vgm,map=F) 
extractPar(prodSumModel_Vgm)

# Metric model
metric <- vgmST("metric", joint = vgm(0.25,"Sph",15,0.1), stAni=0.5)
plot(varg,metric,map=F)
metric_Vgm <- fit.StVariogram(varg, metric, method="L-BFGS-B",lower=pars.l,upper=pars.u)
attr(metric_Vgm, "MSE")
plot(varg,metric_Vgm,map=F) 
extractPar(metric_Vgm)

# Sum metric
sumMetric <- vgmST("sumMetric", space = space_vgm,time = time_vgm, joint = vgm(psill=0.5,"Sph", range=15, nugget=0.1), stAni=1) 
plot(varg,sumMetric,map=F)
sumMetric_Vgm <- fit.StVariogram(varg, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="secs")
attr(sumMetric_Vgm, "MSE")
plot(varg,sumMetric_Vgm,map=F) 
extractPar(sumMetric_Vgm)

# Simple sum metric
SimplesumMetric <- vgmST("simpleSumMetric",space = space_vgm,time = time_vgm, joint = vgm(0.5,"Sph", 15, 0.1), nugget=0.1, stAni=1) 
plot(varg,SimplesumMetric,map=F)
SimplesumMetric_Vgm <- fit.StVariogram(varg, SimplesumMetric,method = "L-BFGS-B",lower=pars.l)
attr(SimplesumMetric_Vgm, "MSE")
plot(varg,SimplesumMetric_Vgm,map=F) 
extractPar(SimplesumMetric_Vgm)

plot(varg,list(separable_Vgm, prodSumModel_Vgm, metric_Vgm, SimplesumMetric_Vgm),all=T,wireframe=T) 






df_pts <- c(1500,1500,2500,2500,2500,2500,3500)
time_lag <- c(2,3,15,2,10,6,15)
t_units <- c("sec")
cutoff <- c(0)
res_time <- c(85.35831,137.9748,592.4871,86.06239,426.7184,269.9019,)
mem_var <- c(1.14,1.15,1.1,1.1,1.1,1.1,1.22)
mean_d <- c(15.353,15.353,15.353,15.333,15.333,15.333,15.333)
max_d <- c(29.704,29.704,29.704,29.704,29.704,29.704,29.704)
