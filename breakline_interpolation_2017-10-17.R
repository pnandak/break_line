################################################################################
#breakline_interp.R

#This script is for interpolating 3D lines from a sparse set of survey data. The
#interpolated lines may serve as useful breaklines in a seperate geostatitiscal 
#analysis. Input data is held in a data frame consisting of eastings, northings,
#elelvations, and feature codes. The data are first interpolated in x,y space
#using Akima splines. The z coordinate is then interpolated in a subsequent
#step. The user must choose the point density of the interpolated line. There are
#sereval plots to review at each step to ensure the optmial performance of the
#algorythm and for blunder detection.

#by Steve Bird
#2017-10-02

################################################################################
source("rstats/sa9_2016_2017-11-02/total_station_dems-2017-10-19/breakline_functions.R")

#load required packages
library(dplyr)
library(akima)
library(rgeos)
library(alphahull)
library(MASS)
library(sp)
library(scatterplot3d)
library(GEOmap)
library(plyr)
library(zoo)
library(ggplot2)

################################################################################
#Import data file after processing the SS data frame in the Python script to
#parse the feature codes:
ss.parsed <- read.csv("data/sa9_r_py_2.csv",header = TRUE)
head(ss.parsed)
#Top Right######################################################################
#TR
#Manually select the feature to process (e.g. TR, BR, EOV, etc.). Some feature
#codes can be combined if needed (e.g. TR and EOVR)
tr <- ss.parsed[,c(1,2,3,6)]
tr$TOP_R <- ifelse(tr$TOP_R=="",NA,1)
good <- complete.cases(tr)
tr <- tr[good, ]
tr <- tr[,1:3]
tr$no <- seq_along(tr$X)
tr$desc <- "TR"

# 
eov.r <- ss.parsed[,c(1,2,3,10)]
eov.r$EOV <- ifelse(eov.r$EOV=="",NA,1)
good <- complete.cases(eov.r)
eov.r <- eov.r[good, ]
eov.r <- eov.r[,1:3]
eov.r$no <- seq_along(eov.r$X)
eov.r$desc <- "eov"
# 

tl <- ss.parsed[,c(1,2,3,7)]
tl$TOP_L <- ifelse(tl$TOP_L=="",NA,1)
good <- complete.cases(tl)
tl <- tl[good, ]
tl <- tl[,1:3]
tl$no <- seq_along(tl$X)
tl$desc <- "TL"

breakline <- rbind(tr,eov.r,tl)
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
breakline$X - lead(breakline$X)




breakline.master <- breakline


#plot the data and determine if the points need to be re-ordered
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(eov.r$X,eov.r$Y,col="red")
points(tl$X,tl$Y,col="blue")
text(breakline$X, breakline$Y,labels=breakline$no)
text(tl$X,tl$Y,labels=tl$no,col="blue")

#Top#####################################################################
#TR1a
#manually re-order the points
breakline <- breakline.master[c(1:8,240),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TR1a <- xyz.int
TR1a.name <- deparse(substitute(TR1a))
len.bl <- length(TR1a.name)
bl.name <- rep(TR1a.name,len.bl)
TR1a[[4]] <- bl.name
TR1a <- as.data.frame(TR1a)
names(TR1a) <- c("x","y","z","line")
################################################################################
#TR1b
#manually re-order the points
breakline <- breakline.master[9:16,]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TR1b <- xyz.int
TR1b.name <- deparse(substitute(TR1b))
len.bl <- length(TR1b.name)
bl.name <- rep(TR1b.name,len.bl)
TR1b[[4]] <- bl.name
TR1b <- as.data.frame(TR1b)
names(TR1b) <- c("x","y","z","line")
################################################################################
#TR1b2
#manually re-order the points
breakline <- breakline.master[c(16,241:242),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TR1b2 <- xyz.int
TR1b2.name <- deparse(substitute(TR1b2))
len.bl <- length(TR1b2.name)
bl.name <- rep(TR1b2.name,len.bl)
TR1b2[[4]] <- bl.name
TR1b2 <- as.data.frame(TR1b2)
names(TR1b2) <- c("x","y","z","line")

TR1b <- rbind(TR1b2,TR1b)
plot(TR1b$x,TR1b$y,type="l")
#################################################################################

#TR1c
#manually re-order the points

breakline <- breakline.master[c(17:33,54:91,107:140,211:213,141:156),]
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TR1c <- xyz.int
TR1c.name <- deparse(substitute(TR1c))
len.bl <- length(TR1c.name)
bl.name <- rep(TR1cd.name,len.bl)
TR1c[[4]] <- bl.name
TR1c <- as.data.frame(TR1c)
names(TR1c) <- c("x","y","z","line")
################################################################################
#TR1d
#manually re-order the points
breakline <- breakline.master[c(214:216,157:173),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TR1d <- xyz.int
TR1d.name <- deparse(substitute(TR1d))
len.bl <- length(TR1d.name)
bl.name <- rep(TR1cd.name,len.bl)
TR1d[[4]] <- bl.name
TR1d <- as.data.frame(TR1d)
names(TR1d) <- c("x","y","z","line")
################################################################################
#TR1cd
counter <- seq_along(TR1c[[1]])
TR1c <- cbind(TR1c,counter)
TR1c <- TR1c[order(TR1c[,5], decreasing = TRUE),]
TR1c <- TR1c[,1:4]

counter <- seq_along(TR1d[[1]])
TR1d <- cbind(TR1d,counter)
TR1d <- TR1d[order(TR1d[,5], decreasing = FALSE),]
TR1d <- TR1d[,1:4]

TR1cd <- rbind(TR1c,TR1d)
plot(TR1cd$x,TR1cd$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1,type="l")
#################################################################################
#TLa
#manually re-order the points
breakline <- breakline.master[c(219:224),]
#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(eov.r$X,eov.r$Y,col="red")
points(tl$X,tl$Y,col="blue")
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
text(tl$X,tl$Y,labels=tl$no,col="blue")
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TLa <- xyz.int
TLa.name <- deparse(substitute(TLa))
len.bl <- length(TLa.name)
bl.name <- rep(TLa.name,len.bl)
TLa[[4]] <- bl.name
TLa <- as.data.frame(TLa)
names(TLa) <- c("x","y","z","line")
#################################################################################
#TLb
#manually re-order the points
breakline <- breakline.master[c(224:231,174:177,232,178:180,233:239,243:258,181:182,261:260,262:269,183:188,296:308,315:323,330:358),]
#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(eov.r$X,eov.r$Y,col="red")
points(tl$X,tl$Y,col="blue")
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.5)
text(tl$X,tl$Y,labels=tl$no,col="blue")
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TLb <- xyz.int
TLb.name <- deparse(substitute(TLb))
len.bl <- length(TLb.name)
bl.name <- rep(TLb.name,len.bl)
TLb[[4]] <- bl.name
TLb <- as.data.frame(TLb)
names(TLb) <- c("x","y","z","line")
#################################################################################
#TLc
#manually re-order the points
breakline <- breakline.master[c(358:370,389:385,217:218,383:374),]
#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(eov.r$X,eov.r$Y,col="red")
points(tl$X,tl$Y,col="blue")
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.5)
text(tl$X,tl$Y,labels=tl$no,col="blue")
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TLc <- xyz.int
TLc.name <- deparse(substitute(TLc))
len.bl <- length(TLc.name)
bl.name <- rep(TLc.name,len.bl)
TLc[[4]] <- bl.name
TLc <- as.data.frame(TLc)
names(TLc) <- c("x","y","z","line")
#################################################################################
#TLd
#manually re-order the points
breakline <- breakline.master[c(371:374),]
#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(eov.r$X,eov.r$Y,col="red")
points(tl$X,tl$Y,col="blue")
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
text(tl$X,tl$Y,labels=tl$no,col="blue")
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TLd <- xyz.int
TLd.name <- deparse(substitute(TLd))
len.bl <- length(TLd.name)
bl.name <- rep(TLd.name,len.bl)
TLd[[4]] <- bl.name
TLd <- as.data.frame(TLd)
names(TLd) <- c("x","y","z","line")
################################################################################
#TL
counter <- seq_along(TLa[[1]])
TLa <- cbind(TLa,counter)
TLa.s <- TLa[order(TLa[,5], decreasing = TRUE),]
TLa.s <- TLa.s[,1:4]

counter <- seq_along(TLb[[1]])
TLb <- cbind(TLb,counter)
TLb.s <- TLb[order(TLb[,5], decreasing = FALSE),]
TLb.s <- TLb.s[,1:4]

counter <- seq_along(TLc[[1]])
TLc <- cbind(TLc,counter)
TLc.s <- TLc[order(TLc[,5], decreasing = TRUE),]
TLc.s <- TLc.s[,1:4]

counter <- seq_along(TLd[[1]])
TLd <- cbind(TLd,counter)
TLd.s <- TLd[order(TLc[,5], decreasing = FALSE),]
TLd.s <- TLd.s[,1:4]

TL <- rbind(TLa.s,TLb.s,TLc.s,TLd.s)
plot(TL$x, TL$y,asp=1,type="l")

################################################################################
#Island 1a

breakline <- breakline.master[c(189:192,50:53),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island1 <- xyz.int
island1.name <- deparse(substitute(island1))
len.bl <- length(island1.name)
bl.name <- rep(island1.name,len.bl)
island1[[4]] <- bl.name
island1a <- as.data.frame(island1)
names(island1a) <- c("x","y","z","line")
################################################################################
#Island 1b

breakline <- breakline.master[c(53,270:279),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island1 <- xyz.int
island1.name <- deparse(substitute(island1))
len.bl <- length(island1.name)
bl.name <- rep(island1.name,len.bl)
island1[[4]] <- bl.name
island1b <- as.data.frame(island1)
names(island1b) <- c("x","y","z","line")
################################################################################
#Island 1

counter <- seq_along(island1a[[1]])
island1a <- cbind(island1a,counter)
island1a <- island1a[order(island1a[,5], decreasing = TRUE),]
island1a <- island1a[,1:4]

counter <- seq_along(island1b[[1]])
island1b <- cbind(island1b,counter)
island1b <- island1b[order(island1b[,5], decreasing = FALSE),]
island1b <- island1b[,1:4]

island1 <- rbind(island1a,island1b)
plot(island1$x, island1$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1,type="l")
################################################################################
#Island 2a

breakline <- breakline.master[c(35:49,99:95),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island2a <- xyz.int
island2a.name <- deparse(substitute(island2a))
len.bl <- length(island2a.name)
bl.name <- rep(island2a.name,len.bl)
island2a[[4]] <- bl.name
island2a <- as.data.frame(island2a)
names(island2a) <- c("x","y","z","line")
################################################################################
#Island 2b

breakline <- breakline.master[c(95:92,196:194,295:287,280:285,35),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island2b <- xyz.int
island2b.name <- deparse(substitute(island2b))
len.bl <- length(island2b.name)
bl.name <- rep(island2b.name,len.bl)
island2b[[4]] <- bl.name
island2b <- as.data.frame(island2b)
names(island2b) <- c("x","y","z","line")
################################################################################
#Island 2

counter <- seq_along(island2a[[1]])
island2a <- cbind(island2a,counter)
island2a <- island2a[order(island2a[,5], decreasing = FALSE),]
island2a <- island2a[,1:4]

counter <- seq_along(island2b[[1]])
island2b <- cbind(island2b,counter)
island2b <- island2b[order(island2b[,5], decreasing = TRUE),]
island2b <- island2b[,1:4]

island2 <- rbind(island2b,island2a)
plot(island2$x, island2$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1,type="l")
################################################################################
#Island 3a

breakline <- breakline.master[c(312:310,200,309,101:105,204:201),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island3a <- xyz.int
island3a.name <- deparse(substitute(island3a))
len.bl <- length(island3a.name)
bl.name <- rep(island3a.name,len.bl)
island3a[[4]] <- bl.name
island3a <- as.data.frame(island3a)
names(island3a) <- c("x","y","z","line")
################################################################################
#Island 3b

breakline <- breakline.master[c(312:314,210:205,329:324),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
island3b <- xyz.int
island3b.name <- deparse(substitute(island3b))
len.bl <- length(island3b.name)
bl.name <- rep(island3b.name,len.bl)
island3b[[4]] <- bl.name
island3b <- as.data.frame(island3b)
names(island3b) <- c("x","y","z","line")
################################################################################
#Island 3c

breakline <- breakline.master[c(324,201),]

#check plot to make sure there're no ordering problems
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no,pos=4,offset = 0.1,cex=.3)
xyz.int <- interp2pts(breakline, 0.1)
island3c <- xyz.int
island3c.name <- deparse(substitute(island3c))
len.bl <- length(island3c.name)
bl.name <- rep(island3c.name,len.bl)
island3c[[4]] <- bl.name
island3c <- as.data.frame(island3c)
names(island3c) <- c("x","y","z","line")
#Island 3

counter <- seq_along(island3a[[1]])
island3a <- cbind(island3a,counter)
island3a <- island3a[order(island3a[,5], decreasing = FALSE),]
island3a <- island3a[,1:4]

counter <- seq_along(island3b[[1]])
island3b <- cbind(island3b,counter)
island3b <- island3b[order(island3b[,5], decreasing = TRUE),]
island3b <- island3b[,1:4]

island3 <- rbind(island3a,island3b,island3c)
plot(island3$x, island3$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1,type="l")

################################################################################
#TH
#Manually select the feature to process (e.g. TR, BR, EOV, etc.). Some feature
#codes can be combined if needed (e.g. TR and EOVR)
thal <- ss.parsed[,c(1,2,3,16)]
thal$THAL <- ifelse(thal$THAL=="",NA,1)
good <- complete.cases(thal)
thal <- thal[good, ]
thal <- thal[,1:3]

breakline <- thal
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
breakline$X - lead(breakline$X)

breakline <- breakline[c(1:25,27,26,28:39,63:45),]
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$new.index)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TH <- xyz.int
TH.name <- deparse(substitute(TH))
len.bl <- length(TH[[3]])
bl.name <- rep(TH.name,len.bl)
TH[[4]] <- bl.name
TH.1 <- as.data.frame(TH)
names(TH.1) <- c("x","y","z","line")
################################################################################
#TH
#Manually select the feature to process (e.g. TR, BR, EOV, etc.). Some feature
#codes can be combined if needed (e.g. TR and EOVR)
thal <- ss.parsed[,c(1,2,3,16)]
thal$THAL <- ifelse(thal$THAL=="",NA,1)
good <- complete.cases(thal)
thal <- thal[good, ]
thal <- thal[,1:3]

breakline <- thal
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
breakline$X - lead(breakline$X)

breakline <- breakline[c(45:40,64:82,84:100,115:101),]
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$new.index)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
TH <- xyz.int
TH.name <- deparse(substitute(TH))
len.bl <- length(TH[[3]])
bl.name <- rep(TH.name,len.bl)
TH[[4]] <- bl.name
TH.2 <- as.data.frame(TH)
names(TH.2) <- c("x","y","z","line")
################################################################################
#TH merge
counter <- seq_along(TH.2[[1]])
TH.2 <- cbind(TH.2,counter)
TH.2 <- TH.2[order(TH.2[,5], decreasing = TRUE),]
TH.2 <- TH.2[,1:4]
TH <- rbind(TH.1,TH.2)
plot(TH$x, TH$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1,type="l")

################################################################################
#Misc breaklines:
unique(ss.parsed$B_LN)
#BLA  BLB  BLC  BLD  BLE  BLF  BLG  BLH  BLI  BLJ  BLK  BLL  BLM  BLN  BLO  BLP  BLQ  BLR  BLS  BLT  BLU  BLV  BLW  BLX  BLY  BLZ  BLAA BLAB BLAC BLAD BLAE BLAF
#BLAG BLAH BLAI BLAJ BLAK BLAL BLAM BLAN BLAO BLAP BLAQ BLAR BLAS BLAT BLAU BLAV BLBA BLBB BLBC BLBD BLBE BLBF BLBG BLBH BLBI BLBJ BLBK BLBL BLBM BLBN BLBO BLBP BLBQ
#BLBR BLBS BLBT BLBU BLBV BLBW BLBX BLBY BLBZ BLCA BLCB BLCC BLCD BLCE BLCF BLCG BLCH BLCI BLCJ BLCK BLCL BLCM BLCN BLCO BLCP BLCQ

bln <- ss.parsed[,c(1,2,3,37)]
breakline <- filter(bln, B_LN == "BLA")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)

#breakline <- breakline[1:10,]
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLA <- xyz.int
BLA.name <- deparse(substitute(BLA))
len.bl <- length(BLA[[3]])
bl.name <- rep(BLA.name,len.bl)
BLA[[4]] <- bl.name
BLA <- as.data.frame(BLA)
names(BLA) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLB")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLB <- xyz.int
BLB.name <- deparse(substitute(BLB))
len.bl <- length(BLB[[3]])
bl.name <- rep(BLB.name,len.bl)
BLB[[4]] <- bl.name
BLB <- as.data.frame(BLB)
names(BLB) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLC")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLC <- xyz.int
BLC.name <- deparse(substitute(BLC))
len.bl <- length(BLC[[3]])
bl.name <- rep(BLC.name,len.bl)
BLC[[4]] <- bl.name
BLC <- as.data.frame(BLC)
names(BLC) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLD")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)

xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLD <- xyz.int
BLD.name <- deparse(substitute(BLD))
len.bl <- length(BLD[[3]])
bl.name <- rep(BLD.name,len.bl)
BLD[[4]] <- bl.name
BLD <- as.data.frame(BLD)
names(BLD) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLE")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLE <- xyz.int
BLE.name <- deparse(substitute(BLE))
len.bl <- length(BLE[[3]])
bl.name <- rep(BLE.name,len.bl)
BLE[[4]] <- bl.name
BLE <- as.data.frame(BLE)
names(BLE) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLF")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLF <- xyz.int
BLF.name <- deparse(substitute(BLF))
len.bl <- length(BLF[[3]])
bl.name <- rep(BLF.name,len.bl)
BLF[[4]] <- bl.name
BLF <- as.data.frame(BLF)
names(BLF) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLG")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
breakline <- breakline[1:3,]
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLG <- xyz.int
BLG.name <- deparse(substitute(BLG))
len.bl <- length(BLG[[3]])
bl.name <- rep(BLG.name,len.bl)
BLG[[4]] <- bl.name
BLG <- as.data.frame(BLG)
names(BLG) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLH")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLH <- xyz.int
BLH.name <- deparse(substitute(BLH))
len.bl <- length(BLH[[3]])
bl.name <- rep(BLH.name,len.bl)
BLH[[4]] <- bl.name
BLH <- as.data.frame(BLH)
names(BLH) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLI")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLI <- xyz.int
BLI.name <- deparse(substitute(BLI))
len.bl <- length(BLI[[3]])
bl.name <- rep(BLI.name,len.bl)
BLI[[4]] <- bl.name
BLI <- as.data.frame(BLI)
names(BLI) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLJ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLJ <- xyz.int
BLJ.name <- deparse(substitute(BLJ))
len.bl <- length(BLJ[[3]])
bl.name <- rep(BLJ.name,len.bl)
BLJ[[4]] <- bl.name
BLJ <- as.data.frame(BLJ)
names(BLJ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLK")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLK <- xyz.int
BLK.name <- deparse(substitute(BLK))
len.bl <- length(BLK[[3]])
bl.name <- rep(BLK.name,len.bl)
BLK[[4]] <- bl.name
BLK <- as.data.frame(BLK)
names(BLK) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLL")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLL <- xyz.int
BLL.name <- deparse(substitute(BLL))
len.bl <- length(BLL[[3]])
bl.name <- rep(BLL.name,len.bl)
BLL[[4]] <- bl.name
BLL <- as.data.frame(BLL)
names(BLL) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLM")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLM <- xyz.int
BLM.name <- deparse(substitute(BLM))
len.bl <- length(BLM[[3]])
bl.name <- rep(BLM.name,len.bl)
BLM[[4]] <- bl.name
BLM <- as.data.frame(BLM)
names(BLM) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLN")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLN <- xyz.int
BLN.name <- deparse(substitute(BLN))
len.bl <- length(BLN[[3]])
bl.name <- rep(BLN.name,len.bl)
BLN[[4]] <- bl.name
BLN <- as.data.frame(BLN)
names(BLN) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLO")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLO <- xyz.int
BLO.name <- deparse(substitute(BLO))
len.bl <- length(BLO[[3]])
bl.name <- rep(BLO.name,len.bl)
BLO[[4]] <- bl.name
BLO <- as.data.frame(BLO)
names(BLO) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLP")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLP <- xyz.int
BLP.name <- deparse(substitute(BLP))
len.bl <- length(BLP[[3]])
bl.name <- rep(BLP.name,len.bl)
BLP[[4]] <- bl.name
BLP <- as.data.frame(BLP)
names(BLP) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLQ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLQ <- xyz.int
BLQ.name <- deparse(substitute(BLQ))
len.bl <- length(BLQ[[3]])
bl.name <- rep(BLQ.name,len.bl)
BLQ[[4]] <- bl.name
BLQ <- as.data.frame(BLQ)
names(BLQ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLR")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLR <- xyz.int
BLR.name <- deparse(substitute(BLR))
len.bl <- length(BLR[[3]])
bl.name <- rep(BLR.name,len.bl)
BLR[[4]] <- bl.name
BLR <- as.data.frame(BLR)
names(BLR) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLS")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLS <- xyz.int
BLS.name <- deparse(substitute(BLS))
len.bl <- length(BLS[[3]])
bl.name <- rep(BLS.name,len.bl)
BLS[[4]] <- bl.name
BLS <- as.data.frame(BLS)
names(BLS) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLT")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLT <- xyz.int
BLT.name <- deparse(substitute(BLT))
len.bl <- length(BLT[[3]])
bl.name <- rep(BLT.name,len.bl)
BLT[[4]] <- bl.name
BLT <- as.data.frame(BLT)
names(BLT) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLU")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLU <- xyz.int
BLU.name <- deparse(substitute(BLU))
len.bl <- length(BLU[[3]])
bl.name <- rep(BLU.name,len.bl)
BLU[[4]] <- bl.name
BLU <- as.data.frame(BLU)
names(BLU) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLV")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLV <- xyz.int
BLV.name <- deparse(substitute(BLV))
len.bl <- length(BLV[[3]])
bl.name <- rep(BLV.name,len.bl)
BLV[[4]] <- bl.name
BLV <- as.data.frame(BLV)
names(BLV) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLW")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLW <- xyz.int
BLW.name <- deparse(substitute(BLW))
len.bl <- length(BLW[[3]])
bl.name <- rep(BLW.name,len.bl)
BLW[[4]] <- bl.name
BLW <- as.data.frame(BLW)
names(BLW) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLX")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLX <- xyz.int
BLX.name <- deparse(substitute(BLX))
len.bl <- length(BLX[[3]])
bl.name <- rep(BLX.name,len.bl)
BLX[[4]] <- bl.name
BLX <- as.data.frame(BLX)
names(BLX) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLY")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLY <- xyz.int
BLY.name <- deparse(substitute(BLY))
len.bl <- length(BLY[[3]])
bl.name <- rep(BLY.name,len.bl)
BLY[[4]] <- bl.name
BLY <- as.data.frame(BLY)
names(BLY) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLZ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLZ <- xyz.int
BLZ.name <- deparse(substitute(BLZ))
len.bl <- length(BLZ[[3]])
bl.name <- rep(BLZ.name,len.bl)
BLZ[[4]] <- bl.name
BLZ <- as.data.frame(BLZ)
names(BLZ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAA")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAA <- xyz.int
BLAA.name <- deparse(substitute(BLAA))
len.bl <- length(BLAA[[3]])
bl.name <- rep(BLAA.name,len.bl)
BLAA[[4]] <- bl.name
BLAA <- as.data.frame(BLAA)
names(BLAA) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAB")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAB <- xyz.int
BLAB.name <- deparse(substitute(BLAB))
len.bl <- length(BLAB[[3]])
bl.name <- rep(BLAB.name,len.bl)
BLAB[[4]] <- bl.name
BLAB <- as.data.frame(BLAB)
names(BLAB) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAC")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLAC <- xyz.int
BLAC.name <- deparse(substitute(BLAC))
len.bl <- length(BLAC[[3]])
bl.name <- rep(BLAC.name,len.bl)
BLAC[[4]] <- bl.name
BLAC <- as.data.frame(BLAC)
names(BLAC) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAD")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAD <- xyz.int
BLAD.name <- deparse(substitute(BLAD))
len.bl <- length(BLAD[[3]])
bl.name <- rep(BLAD.name,len.bl)
BLAD[[4]] <- bl.name
BLAD <- as.data.frame(BLAD)
names(BLAD) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAE")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAE <- xyz.int
BLAE.name <- deparse(substitute(BLAE))
len.bl <- length(BLAE[[3]])
bl.name <- rep(BLAE.name,len.bl)
BLAE[[4]] <- bl.name
BLAE <- as.data.frame(BLAE)
names(BLAE) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAF")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAF <- xyz.int
BLAF.name <- deparse(substitute(BLAF))
len.bl <- length(BLAF[[3]])
bl.name <- rep(BLAF.name,len.bl)
BLAF[[4]] <- bl.name
BLAF <- as.data.frame(BLAF)
names(BLAF) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAG")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAG <- xyz.int
BLAG.name <- deparse(substitute(BLAG))
len.bl <- length(BLAG[[3]])
bl.name <- rep(BLAG.name,len.bl)
BLAG[[4]] <- bl.name
BLAG <- as.data.frame(BLAG)
names(BLAG) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAH")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAH <- xyz.int
BLAH.name <- deparse(substitute(BLAH))
len.bl <- length(BLAH[[3]])
bl.name <- rep(BLAH.name,len.bl)
BLAH[[4]] <- bl.name
BLAH <- as.data.frame(BLAH)
names(BLAH) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAI")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAI <- xyz.int
BLAI.name <- deparse(substitute(BLAI))
len.bl <- length(BLAI[[3]])
bl.name <- rep(BLAI.name,len.bl)
BLAI[[4]] <- bl.name
BLAI <- as.data.frame(BLAI)
names(BLAI) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAJ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAJ <- xyz.int
BLAJ.name <- deparse(substitute(BLAJ))
len.bl <- length(BLAJ[[3]])
bl.name <- rep(BLAJ.name,len.bl)
BLAJ[[4]] <- bl.name
BLAJ <- as.data.frame(BLAJ)
names(BLAJ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAK")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAK <- xyz.int
BLAK.name <- deparse(substitute(BLAK))
len.bl <- length(BLAK[[3]])
bl.name <- rep(BLAK.name,len.bl)
BLAK[[4]] <- bl.name
BLAK <- as.data.frame(BLAK)
names(BLAK) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAL")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAL <- xyz.int
BLAL.name <- deparse(substitute(BLAL))
len.bl <- length(BLAL[[3]])
bl.name <- rep(BLAL.name,len.bl)
BLAL[[4]] <- bl.name
BLAL <- as.data.frame(BLAL)
names(BLAL) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAM")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLAM <- xyz.int
BLAM.name <- deparse(substitute(BLAM))
len.bl <- length(BLAM[[3]])
bl.name <- rep(BLAM.name,len.bl)
BLAM[[4]] <- bl.name
BLAM <- as.data.frame(BLAM)
names(BLAM) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAN")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAN <- xyz.int
BLAN.name <- deparse(substitute(BLAN))
len.bl <- length(BLAN[[3]])
bl.name <- rep(BLAN.name,len.bl)
BLAN[[4]] <- bl.name
BLAN <- as.data.frame(BLAN)
names(BLAN) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAO")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLAO <- xyz.int
BLAO.name <- deparse(substitute(BLAO))
len.bl <- length(BLAO[[3]])
bl.name <- rep(BLAO.name,len.bl)
BLAO[[4]] <- bl.name
BLAO <- as.data.frame(BLAO)
names(BLAO) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAP")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAP <- xyz.int
BLAP.name <- deparse(substitute(BLAP))
len.bl <- length(BLAP[[3]])
bl.name <- rep(BLAP.name,len.bl)
BLAP[[4]] <- bl.name
BLAP <- as.data.frame(BLAP)
names(BLAP) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAQ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAQ <- xyz.int
BLAQ.name <- deparse(substitute(BLAQ))
len.bl <- length(BLAQ[[3]])
bl.name <- rep(BLAQ.name,len.bl)
BLAQ[[4]] <- bl.name
BLAQ <- as.data.frame(BLAQ)
names(BLAQ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAR")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAR <- xyz.int
BLAR.name <- deparse(substitute(BLAR))
len.bl <- length(BLAR[[3]])
bl.name <- rep(BLAR.name,len.bl)
BLAR[[4]] <- bl.name
BLAR <- as.data.frame(BLAR)
names(BLAR) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAS")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAS <- xyz.int
BLAS.name <- deparse(substitute(BLAS))
len.bl <- length(BLAS[[3]])
bl.name <- rep(BLAS.name,len.bl)
BLAS[[4]] <- bl.name
BLAS <- as.data.frame(BLAS)
names(BLAS) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAT")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAT <- xyz.int
BLAT.name <- deparse(substitute(BLAT))
len.bl <- length(BLAT[[3]])
bl.name <- rep(BLAT.name,len.bl)
BLAT[[4]] <- bl.name
BLAT <- as.data.frame(BLAT)
names(BLAT) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAU")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLAU <- xyz.int
BLAU.name <- deparse(substitute(BLAU))
len.bl <- length(BLAU[[3]])
bl.name <- rep(BLAU.name,len.bl)
BLAU[[4]] <- bl.name
BLAU <- as.data.frame(BLAU)
names(BLAU) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLAV")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLAV <- xyz.int
BLAV.name <- deparse(substitute(BLAV))
len.bl <- length(BLAV[[3]])
bl.name <- rep(BLAV.name,len.bl)
BLAV[[4]] <- bl.name
BLAV <- as.data.frame(BLAV)
names(BLAV) <- c("x","y","z","line")
# ################################################################################
# breakline <- filter(bln, B_LN == "BLAW")
# no <- seq_along(breakline$X)
# breakline <- cbind(no,breakline)
# plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
# text(breakline$X, breakline$Y,labels=breakline$no)
# xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
# BLAW <- xyz.int
# BLAW.name <- deparse(substitute(BLAW))
# len.bl <- length(BLAW[[3]])
# bl.name <- rep(BLAW.name,len.bl)
# BLAW[[4]] <- bl.name
# BLAW <- as.data.frame(BLAW)
# names(BLAW) <- c("x","y","z","line")
################################################################################
# breakline <- filter(bln, B_LN == "BLAX")
# no <- seq_along(breakline$X)
# breakline <- cbind(no,breakline)
# plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
# text(breakline$X, breakline$Y,labels=breakline$no)
# xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
# BLAX <- xyz.int
# BLAX.name <- deparse(substitute(BLAX))
# len.bl <- length(BLAX[[3]])
# bl.name <- rep(BLAX.name,len.bl)
# BLAX[[4]] <- bl.name
# BLAX <- as.data.frame(BLAX)
# names(BLAX) <- c("x","y","z","line")
################################################################################
# breakline <- filter(bln, B_LN == "BLAZ")
# no <- seq_along(breakline$X)
# breakline <- cbind(no,breakline)
# plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
# text(breakline$X, breakline$Y,labels=breakline$no)
# breakline <- breakline[1:3,]
# xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
# BLAZ <- xyz.int
# BLAZ.name <- deparse(substitute(BLAZ))
# len.bl <- length(BLAZ[[3]])
# bl.name <- rep(BLAZ.name,len.bl)
# BLAZ[[4]] <- bl.name
# BLAZ <- as.data.frame(BLAZ)
# names(BLAZ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBA")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBA <- xyz.int
BLBA.name <- deparse(substitute(BLBA))
len.bl <- length(BLBA[[3]])
bl.name <- rep(BLBA.name,len.bl)
BLBA[[4]] <- bl.name
BLBA <- as.data.frame(BLBA)
names(BLBA) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBB")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBB <- xyz.int
BLBB.name <- deparse(substitute(BLBB))
len.bl <- length(BLBB[[3]])
bl.name <- rep(BLBB.name,len.bl)
BLBB[[4]] <- bl.name
BLBB <- as.data.frame(BLBB)
names(BLBB) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBC")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- interp2pts(breakline, 0.1)
BLBC <- xyz.int
BLBC.name <- deparse(substitute(BLBC))
len.bl <- length(BLBC[[3]])
bl.name <- rep(BLBC.name,len.bl)
BLBC[[4]] <- bl.name
BLBC <- as.data.frame(BLBC)
names(BLBC) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBD")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBD <- xyz.int
BLBD.name <- deparse(substitute(BLBD))
len.bl <- length(BLBD[[3]])
bl.name <- rep(BLBD.name,len.bl)
BLBD[[4]] <- bl.name
BLBD <- as.data.frame(BLBD)
names(BLBD) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBE")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBE <- xyz.int
BLBE.name <- deparse(substitute(BLBE))
len.bl <- length(BLBE[[3]])
bl.name <- rep(BLBE.name,len.bl)
BLBE[[4]] <- bl.name
BLBE <- as.data.frame(BLBE)
names(BLBE) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBF")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBF <- xyz.int
BLBF.name <- deparse(substitute(BLBF))
len.bl <- length(BLBF[[3]])
bl.name <- rep(BLBF.name,len.bl)
BLBF[[4]] <- bl.name
BLBF <- as.data.frame(BLBF)
names(BLBF) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBG")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBG <- xyz.int
BLBG.name <- deparse(substitute(BLBG))
len.bl <- length(BLBG[[3]])
bl.name <- rep(BLBG.name,len.bl)
BLBG[[4]] <- bl.name
BLBG <- as.data.frame(BLBG)
names(BLBG) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBH")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBH <- xyz.int
BLBH.name <- deparse(substitute(BLBH))
len.bl <- length(BLBH[[3]])
bl.name <- rep(BLBH.name,len.bl)
BLBH[[4]] <- bl.name
BLBH <- as.data.frame(BLBH)
names(BLBH) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBI")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBI <- xyz.int
BLBI.name <- deparse(substitute(BLBI))
len.bl <- length(BLBI[[3]])
bl.name <- rep(BLBI.name,len.bl)
BLBI[[4]] <- bl.name
BLBI <- as.data.frame(BLBI)
names(BLBI) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBJ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBJ <- xyz.int
BLBJ.name <- deparse(substitute(BLBJ))
len.bl <- length(BLBJ[[3]])
bl.name <- rep(BLBJ.name,len.bl)
BLBJ[[4]] <- bl.name
BLBJ <- as.data.frame(BLBJ)
names(BLBJ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBK")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBK <- xyz.int
BLBK.name <- deparse(substitute(BLBK))
len.bl <- length(BLBK[[3]])
bl.name <- rep(BLBK.name,len.bl)
BLBK[[4]] <- bl.name
BLBK <- as.data.frame(BLBK)
names(BLBK) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBL")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBL <- xyz.int
BLBL.name <- deparse(substitute(BLBL))
len.bl <- length(BLBL[[3]])
bl.name <- rep(BLBL.name,len.bl)
BLBL[[4]] <- bl.name
BLBL <- as.data.frame(BLBL)
names(BLBL) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBM")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBM <- xyz.int
BLBM.name <- deparse(substitute(BLBM))
len.bl <- length(BLBM[[3]])
bl.name <- rep(BLBM.name,len.bl)
BLBM[[4]] <- bl.name
BLBM <- as.data.frame(BLBM)
names(BLBM) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBN")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBN <- xyz.int
BLBN.name <- deparse(substitute(BLBN))
len.bl <- length(BLBN[[3]])
bl.name <- rep(BLBN.name,len.bl)
BLBN[[4]] <- bl.name
BLBN <- as.data.frame(BLBN)
names(BLBN) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBO")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBO <- xyz.int
BLBO.name <- deparse(substitute(BLBO))
len.bl <- length(BLBO[[3]])
bl.name <- rep(BLBO.name,len.bl)
BLBO[[4]] <- bl.name
BLBO <- as.data.frame(BLBO)
names(BLBO) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBP")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBP <- xyz.int
BLBP.name <- deparse(substitute(BLBP))
len.bl <- length(BLBP[[3]])
bl.name <- rep(BLBP.name,len.bl)
BLBP[[4]] <- bl.name
BLBP <- as.data.frame(BLBP)
names(BLBP) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBQ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBQ <- xyz.int
BLBQ.name <- deparse(substitute(BLBQ))
len.bl <- length(BLBQ[[3]])
bl.name <- rep(BLBQ.name,len.bl)
BLBQ[[4]] <- bl.name
BLBQ <- as.data.frame(BLBQ)
names(BLBQ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBR")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBR <- xyz.int
BLBR.name <- deparse(substitute(BLBR))
len.bl <- length(BLBR[[3]])
bl.name <- rep(BLBR.name,len.bl)
BLBR[[4]] <- bl.name
BLBR <- as.data.frame(BLBR)
names(BLBR) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBS")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBS <- xyz.int
BLBS.name <- deparse(substitute(BLBS))
len.bl <- length(BLBS[[3]])
bl.name <- rep(BLBS.name,len.bl)
BLBS[[4]] <- bl.name
BLBS <- as.data.frame(BLBS)
names(BLBS) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBT")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBT <- xyz.int
BLBT.name <- deparse(substitute(BLBT))
len.bl <- length(BLBT[[3]])
bl.name <- rep(BLBT.name,len.bl)
BLBT[[4]] <- bl.name
BLBT <- as.data.frame(BLBT)
names(BLBT) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBU")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBU <- xyz.int
BLBU.name <- deparse(substitute(BLBU))
len.bl <- length(BLBU[[3]])
bl.name <- rep(BLBU.name,len.bl)
BLBU[[4]] <- bl.name
BLBU <- as.data.frame(BLBU)
names(BLBU) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBV")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBV <- xyz.int
BLBV.name <- deparse(substitute(BLBV))
len.bl <- length(BLBV[[3]])
bl.name <- rep(BLBV.name,len.bl)
BLBV[[4]] <- bl.name
BLBV <- as.data.frame(BLBV)
names(BLBV) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBW")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBW <- xyz.int
BLBW.name <- deparse(substitute(BLBW))
len.bl <- length(BLBW[[3]])
bl.name <- rep(BLBW.name,len.bl)
BLBW[[4]] <- bl.name
BLBW <- as.data.frame(BLBW)
names(BLBW) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBX")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBX <- xyz.int
BLBX.name <- deparse(substitute(BLBX))
len.bl <- length(BLBX[[3]])
bl.name <- rep(BLBX.name,len.bl)
BLBX[[4]] <- bl.name
BLBX <- as.data.frame(BLBX)
names(BLBX) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBY")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBY <- xyz.int
BLBY.name <- deparse(substitute(BLBY))
len.bl <- length(BLBY[[3]])
bl.name <- rep(BLBY.name,len.bl)
BLBY[[4]] <- bl.name
BLBY <- as.data.frame(BLBY)
names(BLBY) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLBZ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLBZ <- xyz.int
BLBZ.name <- deparse(substitute(BLBZ))
len.bl <- length(BLBZ[[3]])
bl.name <- rep(BLBZ.name,len.bl)
BLBZ[[4]] <- bl.name
BLBZ <- as.data.frame(BLBZ)
names(BLBZ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCA")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCA <- xyz.int
BLCA.name <- deparse(substitute(BLCA))
len.bl <- length(BLCA[[3]])
bl.name <- rep(BLCA.name,len.bl)
BLCA[[4]] <- bl.name
BLCA <- as.data.frame(BLCA)
names(BLCA) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCB")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCB <- xyz.int
BLCB.name <- deparse(substitute(BLCB))
len.bl <- length(BLCB[[3]])
bl.name <- rep(BLCB.name,len.bl)
BLCB[[4]] <- bl.name
BLCB <- as.data.frame(BLCB)
names(BLCB) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCC")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCC <- xyz.int
BLCC.name <- deparse(substitute(BLCC))
len.bl <- length(BLCC[[3]])
bl.name <- rep(BLCC.name,len.bl)
BLCC[[4]] <- bl.name
BLCC <- as.data.frame(BLCC)
names(BLCC) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCD")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCD <- xyz.int
BLCD.name <- deparse(substitute(BLCD))
len.bl <- length(BLCD[[3]])
bl.name <- rep(BLCD.name,len.bl)
BLCD[[4]] <- bl.name
BLCD <- as.data.frame(BLCD)
names(BLCD) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCE")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCE <- xyz.int
BLCE.name <- deparse(substitute(BLCE))
len.bl <- length(BLCE[[3]])
bl.name <- rep(BLCE.name,len.bl)
BLCE[[4]] <- bl.name
BLCE <- as.data.frame(BLCE)
names(BLCE) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCF")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCF <- xyz.int
BLCF.name <- deparse(substitute(BLCF))
len.bl <- length(BLCF[[3]])
bl.name <- rep(BLCF.name,len.bl)
BLCF[[4]] <- bl.name
BLCF <- as.data.frame(BLCF)
names(BLCF) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCG")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCG <- xyz.int
BLCG.name <- deparse(substitute(BLCG))
len.bl <- length(BLCG[[3]])
bl.name <- rep(BLCG.name,len.bl)
BLCG[[4]] <- bl.name
BLCG <- as.data.frame(BLCG)
names(BLCG) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCH")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCH <- xyz.int
BLCH.name <- deparse(substitute(BLCH))
len.bl <- length(BLCH[[3]])
bl.name <- rep(BLCH.name,len.bl)
BLCH[[4]] <- bl.name
BLCH <- as.data.frame(BLCH)
names(BLCH) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCI")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCI <- xyz.int
BLCI.name <- deparse(substitute(BLCI))
len.bl <- length(BLCI[[3]])
bl.name <- rep(BLCI.name,len.bl)
BLCI[[4]] <- bl.name
BLCI <- as.data.frame(BLCI)
names(BLCI) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCJ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCJ <- xyz.int
BLCJ.name <- deparse(substitute(BLCJ))
len.bl <- length(BLCJ[[3]])
bl.name <- rep(BLCJ.name,len.bl)
BLCJ[[4]] <- bl.name
BLCJ <- as.data.frame(BLCJ)
names(BLCJ) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCK")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCK <- xyz.int
BLCK.name <- deparse(substitute(BLCK))
len.bl <- length(BLCK[[3]])
bl.name <- rep(BLCK.name,len.bl)
BLCK[[4]] <- bl.name
BLCK <- as.data.frame(BLCK)
names(BLCK) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCL")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCL <- xyz.int
BLCL.name <- deparse(substitute(BLCL))
len.bl <- length(BLCL[[3]])
bl.name <- rep(BLCL.name,len.bl)
BLCL[[4]] <- bl.name
BLCL <- as.data.frame(BLCL)
names(BLCL) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCM")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCM <- xyz.int
BLCM.name <- deparse(substitute(BLCM))
len.bl <- length(BLCM[[3]])
bl.name <- rep(BLCM.name,len.bl)
BLCM[[4]] <- bl.name
BLCM <- as.data.frame(BLCM)
names(BLCM) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCN")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCN <- xyz.int
BLCN.name <- deparse(substitute(BLCN))
len.bl <- length(BLCN[[3]])
bl.name <- rep(BLCN.name,len.bl)
BLCN[[4]] <- bl.name
BLCN <- as.data.frame(BLCN)
names(BLCN) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCO")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCO <- xyz.int
BLCO.name <- deparse(substitute(BLCO))
len.bl <- length(BLCO[[3]])
bl.name <- rep(BLCO.name,len.bl)
BLCO[[4]] <- bl.name
BLCO <- as.data.frame(BLCO)
names(BLCO) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCP")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCP <- xyz.int
BLCP.name <- deparse(substitute(BLCP))
len.bl <- length(BLCP[[3]])
bl.name <- rep(BLCP.name,len.bl)
BLCP[[4]] <- bl.name
BLCP <- as.data.frame(BLCP)
names(BLCP) <- c("x","y","z","line")
################################################################################
breakline <- filter(bln, B_LN == "BLCQ")
no <- seq_along(breakline$X)
breakline <- cbind(no,breakline)
plot(breakline$X, breakline$Y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
text(breakline$X, breakline$Y,labels=breakline$no)
xyz.int <- auto.breakline(breakline, 0.5, 0.1, 0.001)
BLCQ <- xyz.int
BLCQ.name <- deparse(substitute(BLCQ))
len.bl <- length(BLCQ[[3]])
bl.name <- rep(BLCQ.name,len.bl)
BLCQ[[4]] <- bl.name
BLCQ <- as.data.frame(BLCQ)
names(BLCQ) <- c("x","y","z","line")









all.breaklines <- rbind(TR1a,TR1b,TR1cd,TL,island1,island2,island3,TH,BLA,BLB,BLC,BLD,BLE,BLF,BLG,BLH,BLI,BLJ,BLK,BLL,BLM,BLN,BLO,BLP,BLQ,BLR,BLS,BLT,BLU,BLV,BLW,BLX,BLY,BLZ,BLAA,BLAB,BLAC,BLAD,BLAE,BLAF,BLAG,BLAH,BLAI,BLAJ,BLAK,BLAL,BLAM,BLAN,BLAO,BLAP,BLAQ,BLAR,BLAS,BLAT,BLAU,BLAV,BLBA,BLBB,BLBC,BLBD,BLBE,BLBF,BLBG,BLBH,BLBI,BLBJ,BLBK,BLBL,BLBM,BLBN,BLBO,BLBP,BLBQ,BLBR,BLBS,BLBT,BLBU,BLBV,BLBW,BLBX,BLBY,BLBZ,BLCA,BLCB,BLCC,BLCD,BLCE,BLCF,BLCG,BLCH,BLCI,BLCJ,BLCK,BLCL,BLCM,BLCN,BLCO,BLCP,BLCQ)


no.breakline.pts <- length(all.breaklines$x)
survey.metrics.sa9 <- list(rms.sa9,linear.clos.error,precision,sa9.pts,no.breakline.pts)
names(survey.metrics.sa9) <- c("rms","linear.error","precision","no.pts","no.bl.pts")

bed.spots <- filter(ss.parsed, CTRL == "" & TOP_R == "" & TOP_L == ""  & EOV == "" & is.na(THAL)  )
plot(bed.spots$X,bed.spots$Y)
bed.spots <- bed.spots[,1:3]
names(bed.spots) <- c("x","y","z")
sa9.16 <- rbind(bed.spots,all.breaklines[,1:3])

plot(all.breaklines$x,all.breaklines$y,pch = ".",col="red",cex=3,xlab="Easting",ylab="Northing",asp=1)
points(bed.spots$X,bed.spots$Y,col="black",pch=20)
spots <- ss.parsed[,1:3]

all.bed.spots <- filter(ss.parsed, CTRL == "")
all.bed.spots <- all.bed.spots[,1:3]
good <- complete.cases(all.bed.spots)
all.bed.spots <- all.bed.spots[good,]
all.bed.spots
all.bed.spots <- all.bed.spots[,1:3]
names(all.bed.spots) <- c("x","y","z")
all.sa9.16 <- all.bed.spots
plot(all.sa9.16$x,all.sa9.16$y)


survey.pts.plot.sa9 <- ggplot() +
        geom_path(data = all.breaklines, aes(x, y, group = line),color = "red")  +
        geom_point(data = all.bed.spots,aes(x,y,stroke=0.1,colour = z)) + 
        scale_colour_gradient(low = "black", high = "white",name="Elevation (m asl)") +
        labs(x = "Easting (m)", y = "Northing (m)") + 
        #scale_x_continuous(breaks = c(seq(from=353860, to=353905,by=20)))  +
        #scale_y_continuous(breaks = c(seq(from=5420120, to=5420200,by=20)))  +
        coord_equal() +
        theme_light()

print(survey.pts.plot.sa9)
saveRDS(survey.pts.plot.sa9, "output/survey.pts.plot.sa9.rds")
