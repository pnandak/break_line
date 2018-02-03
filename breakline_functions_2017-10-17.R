# This function creates a 3D interpolation of X,Y,Z points using Akima interpolation

auto.breakline <- function(breakline,pt.spacing,user.max.spacing,min.pt.spacing) {

        
# PART A: A series of functions that plug into the main function ---------------

# A.1 Roation function ---------------------------------------------------------

# This function takes a raw breakline stored as a df and rotates the X and Y
# columns with the objective of eliminating any switchbacks along the X axis. The
# function returns a new df with only the rotated X and Y columns included for
# further processesing.
        
        rot.df <- function(x) {
                rot <- function(x) {
                        rotated.east <- origin.e + (temp.east - origin.e) * cos(x) - (temp.north - origin.n) * sin(x) # Rotate easting by x radians
                        rotated.north <- origin.n + (temp.east - origin.e) * sin(x) + (temp.north - origin.n) * cos(x) # Rotate northing by x radians
                        # Calaculate the number of switchbacks:
                        # Use "lead" in dplyr to subtract the score from one line to the next
                        easting.dir.logical <- ifelse(rotated.east - lag(rotated.east) > 0, 1, 0) # 1=positve, 0=negative
                        easting.dir.logical <- easting.dir.logical[2:length(easting.dir.logical)]
                        east2west <- length(easting.dir.logical[easting.dir.logical==0]) # Count the number of westerly moves
                        west2east <- length(easting.dir.logical[easting.dir.logical==1]) # Count the number of westerly moves
                        min.moves <- ifelse(east2west > west2east, west2east, east2west) # If the direction is E to W, then minimize the number of west2east moves
                        return(min.moves)
                }
                
                origin.e <- breakline$X[1] # Set the rotation point as the first point in the series
                origin.n <- breakline$Y[1] # Set the rotation point as the first point in the series
                temp.east <- breakline$X
                temp.north <- breakline$Y
                rot.angle <- optimize(rot, c(-pi, pi))
                optim.angle <- rot.angle$minimum # This is the optimal rotation angle that minimizes the number of sign changes
                remaining.switchbacks <- rot.angle$objective # Number of switchbacks remaining
                # Rotate the breakline for interpolation
                rotated.east <- origin.e + (temp.east - origin.e) * cos(optim.angle) - (temp.north - origin.n) * sin(optim.angle) # Rotate easting by x radians
                rotated.north <- origin.n + (temp.east - origin.e) * sin(optim.angle) + (temp.north - origin.n) * cos(optim.angle) # Rotate northing by x radians
                rotated.df <- data.frame(rotated.east,rotated.north)
                rotated.list <- list(rotated.df,optim.angle,remaining.switchbacks)
                # The data are still rotated at this point and will be rotated back to the original orientation at the end of the script
                rotated.list
        }

        
# A.2 Stretch function ---------------------------------------------------------

# If there are switchbacks remainging in the data, this next block of code will
# simply stretch the data out before running the interpolation.
        
        stretch.x <- function(rotated.east) {
                
                easting.dir.logical <- ifelse(rotated.east - lag(rotated.east) >0,1,0) # 1=positve, 0=negative
                easting.dir.logical[1] <- easting.dir.logical[2] # Get rid of the NA
                # Depending if the line is surveyed east to west OR west to east, most numbers should be 1 OR 0
                east2west <- length(easting.dir.logical[easting.dir.logical==1]) # Count the numer of easterly moves
                west2east <- length(easting.dir.logical[easting.dir.logical==0]) # Count the numer of westerly moves
                direction <- ifelse(east2west>west2east,1,0) # Binary output: the most moves in one direction gives the overall direction
                # In order to simplify the code, if the direction of the line is negative then
                # the line is reflected here to a positve direction and then reflected back to
                # negative at the end of the operation.
                x.length <- length(easting.dir.logical)
                max.X <- max(rotated.east) + 1 # Reflection factor
                max.X.rep <- rep(max.X,x.length)
                direction.rep <- rep(direction,x.length)
                rotated.east <- ifelse(direction.rep==0,max.X.rep-rotated.east,rotated.east)
                
                dist2shot <- rotated.east - lag(rotated.east) # Distance between sucessive shots along the x axis
                # Shift all switchback points by the distance of the switchback
                stretched.dist <- abs(dist2shot)
                shift <- stretched.dist - dist2shot # All unadjusted points == 0; switchbacks will be adjusted
                shift[1] <- 0
                uniqe.shifts <- unique(shift) # Find all the unique shifts in the data
                uniqe.shifts.len <- length(uniqe.shifts)
                cum.uniqe.shifts <-cumsum(uniqe.shifts) # These are the cumulative shifts in the data that get applied to the line
                shift.index <- match(shift,shift) # Find the vector indicies that match the location of a shift
                unique.shift.bins <- unique(shift.index) # Find the unique shifts
                index.rows <- rep(1:x.length) # Create a vector with the unique shifts
                unique.shift.bins <- c(unique.shift.bins,x.length)
                index.factor <- .bincode(index.rows, unique.shift.bins, FALSE,FALSE)
                index.factor[x.length] <- max(uniqe.shifts.len) # Replace the NA at the end
                index.factor <- data.frame(index.factor)
                levels <- as.integer(seq_along(1:uniqe.shifts.len))
                levels <- data.frame(levels,cum.uniqe.shifts)
                names(levels) <- c("index.factor","cum.uniqe.shifts")
                final.shift <- left_join(index.factor, levels)
                adjusted.dist <- final.shift$cum.uniqe.shifts + rotated.east
                adjustment.df <- data.frame(final.shift,adjusted.dist)
                min.adjust.class <- ddply(adjustment.df,"index.factor", function(df)min(df$adjusted.dist))
                max.adjust.class <- ddply(adjustment.df,"index.factor", function(df)max(df$adjusted.dist))
                adjust.class <- data.frame(min.adjust.class,max.adjust.class,cum.uniqe.shifts)
                adjust.class$index.factor.1 <- NULL
                names(adjust.class) <- c("index","min.X","max.X","cum.adjust")
                all.adjusts.bins <- c(adjust.class$min.X,adjust.class$max.X)
                all.adjusts.bins <- sort(all.adjusts.bins)
                stretch.list <- list(final.shift,adjusted.dist,all.adjusts.bins,adjust.class,uniqe.shifts.len,max.adjust.class,max.X)
                stretch.list
        }

        
# A.3 Maximum point spacing function -------------------------------------------
        
# Function to measure the maximum point spacing along an interplated line based
# on the number of points produced by the aspline function. The User sets the
# average point spacing (pt.spacing). Note that the aspline function simplifies
# the line (remove vertex's) along straight sections of the line. This is good for
# plotting a line in 2D but not for interpolating a surface in 3D.
        
        max.spacing <- function(pt.spacing) {
                breakline.rough <- aspline(rotated.east, rotated.north, n=200)
                # Convert the list holding the spline to a matrix
                breakline.matrix <- cbind(breakline.rough$x,breakline.rough$y)
                # Calculate the distance between each vertex in the spline
                breakline.vertex<- (spDists(breakline.matrix,segments = TRUE))
                # Sum all the vertex distances to get the breakline length (m)
                breakline.length <- sum(breakline.vertex)
                pts <- as.integer(breakline.length/pt.spacing) # Total number of points for aspline function
                breakline <- aspline(rotated.east,rotated.north, n=pts) # Akima interpolation
                breakline <- data.frame(breakline)
                breakline.matrix <- cbind(breakline$x,breakline$y)
                # Calculate the distance between each vertex in the spline
                breakline.vertex <- (spDists(breakline.matrix,segments = TRUE))
                target.spacing <- max(breakline.vertex)
                target.spacing
        }

        
# A.4 Final point spacing function ---------------------------------------------
        
# The max.spacing function calculates the maximum spacing for a given average
# pt.spacing set by the User. This function adjusts the User's estimate and
# derives a max.spacing that meets a value set by the User.
        
        final.pt.spacing <- function(min.pt.spacing) {
                while(min.pt.spacing < .09) {
                        # Max.spacing(pt.spacing)
                        # Print(target.spacing)
                        refined.pt.spacing <- min.pt.spacing
                        min.pt.spacing <- max.spacing(min.pt.spacing)
                        min.pt.spacing <- min.pt.spacing + 0.0001 #increase the density for next interation
                        #print(refined.pt.spacing)
                }
                refined.pt.spacing
        }
        
        
# A.5 Shrink function ----------------------------------------------------------
        
# This function shrinks the interpolated line IF the line was stretched using
# the Stretch function.
        
        shrink.x <- function(input) {
                interp.bins <- .bincode(stretch.x.interpolated, all.adjusts.bins, TRUE)
                interp.bins[1] <- interp.bins[2] # Get rid of the first NA
                breakline.interp <- data.frame(stretch.x.interpolated,breakline.interp$y,interp.bins)
                names(breakline.interp) <- c("x","y","index")
                adjust.iso <- adjust.class$cum.adjust
                adjust.iso <- adjust.iso[2:uniqe.shifts.len]
                adjust.iso <- rep(adjust.iso,each = 2)
                adjust.iso <- c(0,adjust.iso)
                index <- seq_along(1:length(adjust.iso))
                adjustment.index <- data.frame(adjust.iso,index)
                breakline.interp <- left_join(breakline.interp, adjustment.index)
                
                special.adjust <- max.adjust.class[1:(uniqe.shifts.len-1),]
                special.adjust$index.factor <- special.adjust$index.factor + 1
                special.adjust$index <- seq(from=2,to = (uniqe.shifts.len-1)*2,by=2)
                #max.special.adjust <- all.adjusts.bins - lag(all.adjusts.bins)
                #max.special.adjust <- max.special.adjust[2:(length(max.special.adjust))]
                #special.adjust$max.special.adjust <- max.special.adjust[seq(2,length(max.special.adjust),2)]
                #special.adjust$cum.adjust <- adjust.class$cum.adjust[2:length(adjust.class$cum.adjust)]
                breakline.interp <- left_join(breakline.interp, special.adjust)
                breakline.interp$x.diff <- lead(breakline.interp$x)-breakline.interp$x
                #fn <- function(x) {cumsum(x)}
                breakline.interp <- ddply(breakline.interp,"index", transform,cum = cumsum(x.diff))
                breakline.interp$final.x <- ifelse(is.na(breakline.interp$V1),
                                                   breakline.interp$x-breakline.interp$adjust.iso,
                                                   breakline.interp$V1 - breakline.interp$cum)
                reflect.back <- rep(max.X,length(breakline.interp$x))
                breakline.interp$final.x.flip <- reflect.back - breakline.interp$final.x
                #plot(breakline.interp$final.x,breakline.interp$y,asp=1,type="line")
                breakline.interp
        }
        
# A.6 Rotate back to original function -----------------------------------------
        
# This function reverses the function that rotated the data to minimize the
# number of switchbacks.
        
        rot.back <- function(x) {
                origin.e <- breakline$X[1] # Set the rotation point as the first point in the series
                origin.n <- breakline$Y[1] # Set the rotation point as the first point in the series
                original.east <- origin.e + (breakline.interp$x - origin.e) * cos(x) - (breakline.interp$y - origin.n) * sin(x) # Rotate easting by x radians
                original.north <- origin.n + (breakline.interp$x - origin.e) * sin(x) + (breakline.interp$y - origin.n) * cos(x) # Rotate northing by x radians
                final.interp <- data.frame(original.east,original.north)
                final.interp
        }

        
# A.7 Reduce point spacing function --------------------------------------------
        
#Function to thin a line to a regular spacing after the line had been densified
        
        line.thinner <- function(user.max.spacing) {
                breakline.final.matrix <- as.matrix(breakline.final)
                # Calculate the distance between each vertex in the spline
                breakline.final.vertex <- (spDists(breakline.final.matrix,segments = TRUE))
                # Sum all the vertex distances to get the breakline length (m)
                breakline.final.length <- sum(breakline.final.vertex)
                # Calculate the number of points and convert to an interger
                
                # The preceeding code ensured the maximum spacing between points
                # matched the User's desired spacing. However, now we have too
                # many points spaced to closely together. The next step is to
                # thin the line to an even spacing that matches the user's
                # desired spacing.
                cumsum.line.len <- cumsum(breakline.final.vertex) # Cumsum of line length
                bins <- seq(from = 0, to = max(cumsum.line.len) + user.max.spacing, by = user.max.spacing) #create bins along the line
                breakline.bins <- .bincode(cumsum.line.len, bins, TRUE)
                bin.index <- lead(breakline.bins)-breakline.bins
                bin.index <- c(bin.index,1) # Lead produces an NA at the end of the series - replace with a 1
                bin.index[1] <- 1 # Also, make sure we select the first shot in the series
                breakline.final$index <- bin.index
                breakline.final$index <- ifelse(breakline.final$index==1,1,NA)
                good <- complete.cases(breakline.final)
                breakline <- breakline.final[good,1:2]
                breakline
        }
        

# A.8 Function to calculate the minimum bounding region ------------------------

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

        
# A.9 Function to pad the minimum bounding region ------------------------------
        
        breakline.pad <- function() {
                plane.model <- lm(Z ~ X + Y, data=breakline)
                points <-  as.matrix(breakline[,c(1,2)]) # Data need to be in matrix form
                mbr <- MBR(points) # Apply the mbr function
                
                # Plot the hull, the MBR, and the points:
                limits <- apply(mbr, 2, function(x) c(min(x),max(x))) # Plotting limits
                plot(ashape(points, alpha=1000), col="Gray", pch=20, 
                     xlim=limits[,1], ylim=limits[,2]) #T he hull
                lines(mbr, col="Blue", lwd=1) # The MBR
                points(points, pch=1) # The points

                # Dave's code to convert the mbr to a spatial polygon
                outline.points = as.matrix(points)
                # Do some magic 
                outline.polygon = Polygon(mbr)
                outline.polys = Polygons(list(outline.polygon), 1)
                # This creates a "spatial" object
                mbr.outline = SpatialPolygons(list(outline.polys), proj4string = CRS("+proj=utm +zone=10"))
                plot(mbr.outline)

                # Buffer the mbr
                mbr.buffer <- gBuffer(mbr.outline, width=.1,capStyle='SQUARE',joinStyle='MITRE',mitreLimit=10)
                # I'm not sure if there is a default setting to get this function to give a
                # rectangle. I've set the mitre limit high to get a rectangular buffer (adjust as
                # needed). (It turns out it a mitre edge doesn't affect the outcome of the
                # script)

                plot(mbr.buffer)
                # Dave's trick to extract x,y points from spatial polygon
                buffer.xy <- mbr.buffer@polygons[[1]]@Polygons[[1]]@coords
                buffer.xy <- as.data.frame(buffer.xy)
                buffer.xy <- buffer.xy[1:4,] # The subset works if the buffer is a rectangle
                names(buffer.xy) <- c('X','Y')
                buffer.z <- predict(plane.model,new=buffer.xy)
                buffer.xyz <- cbind(buffer.xy,buffer.z)
                names(buffer.xyz) <- c('X','Y','Z')
                breakline.padded <- rbind(buffer.xyz,breakline[,1:3]) # Original data + 4 interpolated z points
                breakline.padded

        }
        
# A.10 Function to interpolate an elevation along a 2D line --------------------
        
        interp.z.NA <- function() {
                easting.dir.logical <- ifelse(rotated.east - lag(rotated.east) >0,1,0) # 1=positve, 0=negative
                easting.dir.logical[1] <- easting.dir.logical[2] # Get rid of the NA
                # Depending if the line is surveyed east to west OR west to east, most numbers should be 1 OR 0
                east2west <- length(easting.dir.logical[easting.dir.logical==1]) # Count the numer of easterly moves
                west2east <- length(easting.dir.logical[easting.dir.logical==0]) # Count the numer of westerly moves
                direction <- ifelse(east2west>west2east,1,0) # Binary output: the most moves in one direction gives the overall direction
                # If the opening and/or closing elevation contains an NA, then paste in the original surveyed Z before linear interpolation
                one.end.z <- head(breakline,1)
                other.end.z <- tail(breakline,1)
                z.original <- ifelse(direction==0,c(other.end.z$Z,one.end.z$Z),c(one.end.z$Z,other.end.z$Z))
                xyz.int$z[1] <- ifelse(is.na(head(xyz.int$z,1))==TRUE, z.original[1], head(xyz.int$z,1))
                len.z <- length(xyz.int$z)
                xyz.int$z[len.z] <- ifelse(is.na(tail(xyz.int$z,1))==TRUE, z.original[2],tail(xyz.int$z,1))
                xyz.int$z <- na.approx(xyz.int$z) # Use the zoo function to do a linear interpolation of NA's
                xyz.int
        }


#PART B: Main function that runs all the preceeding funcitons ------------------

# B.1 Rotate the line to minimize the number of switchbacks --------------------
        rotated.list <- rot.df(breakline)
        # Unpack the data.frame so I don't have to rewrite a bunch of code:
        rotated.df <- as.data.frame(rotated.list[1])
        rotated.east <- rotated.df$rotated.east
        rotated.north <- rotated.df$rotated.north
        rot.angle <- as.numeric(rotated.list[2]) # Optimized rotation angle in rads
        no.switchbacks <- as.numeric(rotated.list[3]) # Number of switchbacks remaining after rotation
        # If there are still switchbacks in the data, resort to a simple linear stretch to fix the rest:
        if(no.switchbacks > 0) {
                stretch.list <- stretch.x(rotated.east)
                # If the data was stretched, then assign each item in the resulting list to a variable:
                final.shift <- stretch.list[[1]]
                rotated.east <- stretch.list[[2]] # This is the "adjusted distance" in the function
                all.adjusts.bins <- stretch.list[[3]]
                adjust.class <- stretch.list[[4]]
                uniqe.shifts.len <- stretch.list[[5]]
                max.adjust.class <- stretch.list[[6]]
                max.X <- stretch.list[[7]]
        }
        
# B.2 Interpolate breakline ----------------------------------------------------
        
        pt.spacing <- final.pt.spacing(min.pt.spacing) # This is the pt spacing used in the final interpolation
        breakline.rough <- aspline(rotated.east, rotated.north, n=200)
        # Convert the list holding the spline to a matrix
        breakline.matrix <- cbind(breakline.rough$x,breakline.rough$y)
        # Calculate the distance between each vertex in the spline
        breakline.vertex<- (spDists(breakline.matrix,segments = TRUE))
        # Sum all the vertex distances to get the breakline length (m)
        breakline.length <- sum(breakline.vertex)
        pts <- as.integer(breakline.length/pt.spacing) # Total number of points for aspline function
        breakline.interp <- aspline(rotated.east,rotated.north, n=pts) #
        plot(breakline.interp,asp=1)
        
        # Remove the stretch if one was applied
        if(no.switchbacks > 0) {
                stretch.x.interpolated <- breakline.interp$x
                shrinked.temp <- shrink.x(stretch.x.interpolated)
                breakline.interp$x <- shrinked.temp$final.x.flip
        }
        
        # Rotate back to the original orientation
        breakline.final <- rot.back(-1*rotated.list[[2]])
        breakline.xy <- line.thinner(user.max.spacing) #thins the line to the desired spacing
        
        # Plot the breakline in 2D
        breakline.xy <- as.data.frame(breakline.xy)
        plot(breakline.xy,type="l",asp=1)
        points(breakline$X,breakline$Y,col="red")

        # Interpolate z for each point along the new line:
        breakline.padded <- breakline.pad()
        x <- breakline.padded$X
        y <- breakline.padded$Y
        z <- breakline.padded$Z
        x.int <- breakline.xy$original.east # Extract the x coordinate of the interpolated plan line to a vector
        y.int <- breakline.xy$original.north # Extract the y coordinate of the interpolated plan line to a vector
        # Interpolate the z coordinate for the plan line
        options( warn = -1 )
                xyz.int <-
                interpp(x, y, z, x.int, y.int, linear = TRUE, extrap = FALSE,
                        duplicate = "error", dupfun = NULL,
                        jitter = 10 ^ -12, jitter.iter = 6, jitter.random = FALSE
                )
        
                
                
                
                
                
                
                
                
                
                
        ################################################################################
        
        xyz.int <- interp.z.NA()
        s3d <- scatterplot3d(
                x = xyz.int$x, y = xyz.int$y, z = xyz.int$z,xlab = "Easting",ylab = "Northing",zlab =
                        "Elevation", main = "Breakline", col.grid = "lightblue", pch = 20,angle = 120)

        s3d$points3d(x,y,z,col="red",pch=NULL,type="p")
        xyz.int
        
}



interp.z.NA <- function() {
        easting.dir.logical <- ifelse(rotated.east - lag(rotated.east) >0,1,0) #1=positve, 0=negative
        easting.dir.logical[1] <- easting.dir.logical[2] #get rid of the NA
        #depending if the line is surveyed east to west OR west to east, most numbers should be 1 OR 0
        east2west <- length(easting.dir.logical[easting.dir.logical==1]) #count the numer of easterly moves
        west2east <- length(easting.dir.logical[easting.dir.logical==0]) #count the numer of westerly moves
        direction <- ifelse(east2west>west2east,1,0) #binary output: the most moves in one direction gives the overall direction
        #if the opening and/or closing elevation contains an NA, then paste in the original surveyed Z before linear interpolation
        one.end.z <- head(breakline,1)
        other.end.z <- tail(breakline,1)
        z.original <- ifelse(direction==0,c(other.end.z$Z,one.end.z$Z),c(one.end.z$Z,other.end.z$Z))
        xyz.int$z[1] <- ifelse(is.na(head(xyz.int$z,1))==TRUE, z.original[1], head(xyz.int$z,1))
        len.z <- length(xyz.int$z)
        xyz.int$z[len.z] <- ifelse(is.na(tail(xyz.int$z,1))==TRUE, z.original[2],tail(xyz.int$z,1))
        xyz.int$z <- na.approx(xyz.int$z) #use the zoo function to do a linear interpolation of NA's
        xyz.int
}



################################################################################
#If the breakline has only two points, then the only way to interpolate is by using a linear model:
interp2pts <- function(breakline,pt.spacing) {
        breakline.matrix <- cbind(breakline$X,breakline$Y)
        #calculate the distance between each vertex in the spline
        breakline.vertex<- (spDists(breakline.matrix,segments = TRUE))
        #sum all the vertex distances to get the breakline length (m)
        breakline.length <- sum(breakline.vertex)
        n.pts <- as.integer(breakline.length/pt.spacing)
        xy <- approx(breakline$X, breakline$Y, n=n.pts)
        xz <- approx(breakline$X, breakline$Z, n=n.pts)
        xy <- as.data.frame(xy)
        xz <- as.data.frame(xz)
        names(xz) <- c("x","z")
        xyz <- cbind(xy,xz$z)
        names(xyz) <- c("x","y","z")
        xyz <- as.list(xyz)
        return(xyz)
}






