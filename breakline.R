# Breakline Interpolation ------------------------------------------------------

# This script is for interpolating 3D lines from a sparse set of survey data. The
# interpolated lines may serve as useful breaklines in a seperate geostatitiscal 
# analysis. Input data is held in a data frame consisting of eastings, northings,
# elelvations, and feature codes. The data are first interpolated in x,y space
# using Akima splines. The z coordinate is then interpolated in a subsequent
# step. The user must choose the point density of the interpolated line. There 
# are sereval plots to review at each step to ensure the optmial performance of
# the algorythm and for blunder detection.

# by Steve Bird
# https://github.com/stephen-bird/
# R version 3.4.3 (2017-11-30) -- "Kite-Eating Tree"
# 2018-01-22


# Load packages ----------------------------------------------------------------

library(plyr) # Check to see if plyr can be eliminated with dplyr
library(dplyr)
library(akima)
library(rgeos)
library(alphahull)
library(MASS)
library(sp)
library(scatterplot3d)
library(GEOmap)
library(zoo)
library(ggplot2)


# Set-up user-specific preliminariers ------------------------------------------

setwd("~/Documents/Projects/efn_2017/rstats") # Set working directory
source("break_line/breakline_functions_2017-10-17.R") # Read functions held in a seperate script
survey_data <- read.csv("data/sa2_r_py_2.csv", header = TRUE) # Read data
survey_data <- survey_data[,1:4]
# Data must be restricted to four columns with the following header: x, y, z, code
# Required feature codes are given in the accompanying text file: feature_codes_2018-01-17.txt
names(survey_data) <- c("x","y","z","code")




