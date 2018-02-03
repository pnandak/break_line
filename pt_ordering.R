library("dplyr")
source("break_line/breakline_functions_2017-10-17.R") # Read functions held in a seperate script

xs_data <- read.csv("/Users/sab/Dropbox/EFN\ Phase\ 2/Final\ recoded\ data/SA2/SA21971recode.csv", header = TRUE)
hubs <- filter(xs_data, hub == "LBH" | hub == "RBH")
hubs <- arrange(hubs, CS,hub,desc(Elevation)) # Sort the hubs

duplicate_hub <- duplicated(hubs$Easting,hubs$Northing) # Look for ground shot at hub (comes w/ lower Z but same X,Y)
hubs[duplicated(!duplicated(hubs$Easting,hubs$Northing)), ]
hubs <- hubs[-c(40,42),]

# Calculate the equation of a line for each cross-section ----------------------

hubs <- arrange(hubs, hub, CS) # Sort the hubs by bank and then cross section
# Rearrange so that LB and RB hubs are in seperate columns matched by cross section number:
lb_hubs <- filter(hubs, hub == "LBH")
lb_hubs <- lb_hubs[,-c(1,8:10)] # Delete uneeded columns
names(lb_hubs) <- c("site","year","cs","e_lb","n_lb","elv_lb") 
rb_hubs <- filter(hubs, hub == "RBH")
rb_hubs <- rb_hubs[,-c(1:3,8:10)]  # Delete uneeded columns
names(rb_hubs) <- c("cs_rb","e_rb","n_rb","elv_rb") 
hubs <- cbind(lb_hubs,rb_hubs)
stopifnot((hubs$cs == hubs$cs_rb)) # Sanity check: hubs are not ordered and/or paired correctly
hubs <- hubs[,-c(7)] # Delete uneeded columns
hubs$m <- (hubs$n_rb - hubs$n_lb) / (hubs$e_rb - hubs$e_lb) # Find the slopes
hubs$b <- hubs$n_lb - hubs$m * hubs$e_lb  # Find the intercepts

# Create a df holding all combinations of cross sections w/ required data ------

combos <- combn(hubs$cs,2,simplify = TRUE)  # Generate all combinations of cross sections (1 vs 2, 1 vs 3..., etc.)
upper_cs <- combos[1,] # Upstream cross section
lower_cs <- combos[2,] # Downstream cross setion
cs_grid <- data.frame(upper_cs,lower_cs)
cs_grid$upper_cs <- as.numeric(cs_grid$upper_cs)
cs_grid$lower_cs <- as.numeric(cs_grid$lower_cs)
# Match df containing all combinations of cross sections with data from the "hubs" df:
hubs$upper_cs <- as.numeric(hubs$cs)
intersect_upper_cs <- left_join(cs_grid, hubs)
hubs$lower_cs <- as.numeric(hubs$cs)
hubs$upper_cs <- NULL
intersect_lower_cs <- left_join(cs_grid, hubs)
# Create a new df with only the colums required to complete the analysis
cs_int <- data.frame(intersect_upper_cs$e_lb,intersect_upper_cs$n_lb,intersect_upper_cs$e_rb,intersect_upper_cs$n_rb,intersect_upper_cs$upper_cs,intersect_upper_cs$lower_cs,intersect_upper_cs$m,intersect_upper_cs$b,intersect_lower_cs$m,intersect_lower_cs$b)
names(cs_int) <- c("e_lb","n_lb","e_rb","n_rb","upper_cs","lower_cs","upper_m","upper_b","lower_m","lower_b")

# Determine if the lines intersect within the study area -----------------------

# Find the point where the lines cross:
cs_int$x <- (cs_int$lower_b - cs_int$upper_b) / (cs_int$upper_m - cs_int$lower_m) # Common x coordinate
cs_int$y <- cs_int$upper_m * cs_int$x + cs_int$upper_b # Common y coordinate
# Data needs to be ordered from smallest to largest since the subsequent function doesn't understand coordinates:
cs_int$e1_order <- ifelse(cs_int$e_lb > cs_int$e_rb,cs_int$e_rb,cs_int$e_lb)
cs_int$e2_order <- ifelse(cs_int$e_lb > cs_int$e_rb,cs_int$e_lb,cs_int$e_rb)
cs_int$n1_order <- ifelse(cs_int$n_lb > cs_int$n_rb,cs_int$n_rb,cs_int$n_lb)
cs_int$n2_order <- ifelse(cs_int$n_lb > cs_int$n_rb,cs_int$n_lb,cs_int$n_rb)

# Determine if the lines cross within the study area or at some point beyond the hubs:
findInt <- function(value, start, end) {
        start < value & end > value
}

cs_int$x_intersect <- findInt(cs_int$x, cs_int$e1_order, cs_int$e2_order)
cs_int$y_intersect <- findInt(cs_int$y, cs_int$n1_order, cs_int$n2_order)
cs_int$intersect <- ifelse(cs_int$y_intersect == TRUE & cs_int$x_intersect == TRUE,TRUE,FALSE)
# Clean up the df
cs_int$e1_order <- NULL
cs_int$e2_order <- NULL
cs_int$n1_order <- NULL
cs_int$n2_order <- NULL
cs_int$x_intersect <- NULL
cs_int$y_intersect <- NULL
# Generate a df populated w/ cross section data from cross sections that intersect within the study area
intersecting_cs <- filter(cs_int,intersect == TRUE)
intersecting_cs$intersect <- NULL
colnames(intersecting_cs)[colnames(intersecting_cs) == 'x'] <- 'x_intersect'
colnames(intersecting_cs)[colnames(intersecting_cs) == 'y'] <- 'y_intersect'

# Order the cross sections for plotting ----------------------------------------

# Calculate distance from the LB hub to the point of intersection:
intersecting_cs$lb_hub2inter <- sqrt((intersecting_cs$y_intersect - intersecting_cs$n_lb)^2 +  (intersecting_cs$x_intersect - intersecting_cs$e_lb)^2 )
# Subset the feature of interest
feat_raw <- filter(xs_data,attribs == "t")
colnames(feat_raw)[colnames(feat_raw) == 'CS'] <- 'upper_cs'
feat <- left_join(feat_raw, intersecting_cs) # Join the feature to the order data
feat$lb_hub2pt <- sqrt((feat$Northing - feat$n_lb)^2 + (feat$Easting - feat$e_lb)^2) # Distance of the feature from LB hub
feat$in_order <- ifelse(feat$lb_hub2inter>feat$lb_hub2pt,TRUE,FALSE) # Is the feature on the left of the intersection?
feat_sub <- data.frame(feat$upper_cs,feat$lower_cs,feat$in_order)
feat_sub$feat.in_order <- ifelse(is.na(feat_sub$feat.in_order),TRUE,feat_sub$feat.in_order) # Convert "NA" to "TRUE"
# If "FALSE", then a re-ordering is needed:
feat_sub$reorder <- ifelse(feat_sub$feat.in_order == "TRUE",feat_sub$feat.upper_cs,feat_sub$feat.lower_cs)
re_order <- unique(feat_sub$reorder) # Remove the duplicate cross sections from all the generated combinations
re_order <- data.frame(re_order)
names(re_order) <- c("upper_cs")
t_BL <- left_join(re_order,feat_raw)
t_BL <- t_BL[,-c(1,2,3,4,8,10)]
names(t_BL) <- c("X","Y","Z","desc")

breakline <- t_BL
interp_bl <- auto.breakline(breakline,0.5, 0.1, 0.001)
interp_bl <- as.data.frame(interp_bl)


