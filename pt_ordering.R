library("dplyr")
source("break_line/breakline_functions_2017-10-17.R") # Read functions held in a seperate script

xs_data <- read.csv("/Users/sab/Dropbox/EFN\ Phase\ 2/Final\ recoded\ data/SA2/SA21971recode.csv", header = TRUE)
hubs <- filter(xs_data, hub == "LBH" | hub == "RBH")
hubs <- arrange(hubs, CS,hub,desc(Elevation)) # Sort the hubs

duplicate_hub <- duplicated(hubs$Easting,hubs$Northing) # Look for ground shot at hub (comes w/ lower Z but same X,Y)
hubs[duplicated(!duplicated(hubs$Easting,hubs$Northing)), ]
#hubs <- hubs[-c(40,42),]

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
inter_upper_cs <- left_join(cs_grid, hubs)

hubs$lower_cs <- as.numeric(hubs$cs)
hubs$upper_cs <- NULL
inter_lower_cs <- left_join(cs_grid, hubs)
# Create a new df with only the columns required to complete the analysis
cs_int <- data.frame(inter_upper_cs$upper_cs,inter_upper_cs$e_lb,inter_upper_cs$n_lb,inter_upper_cs$e_rb,inter_upper_cs$n_rb,inter_upper_cs$m,inter_upper_cs$b,inter_upper_cs$lower_cs,inter_lower_cs$e_lb,inter_lower_cs$n_lb,inter_lower_cs$e_rb,inter_lower_cs$n_rb,inter_lower_cs$m,inter_lower_cs$b)
names(cs_int) <- c("upper_cs","e_lb_up","n_lb_up","e_rb_up","n_rb_up","upper_m","upper_b","lower_cs","e_lb_low","n_lb_low","e_rb_low","n_rb_low","lower_m","lower_b")

# Determine if the lines intersect within the study area -----------------------

# Find the point where the lines cross:
cs_int$x <- (cs_int$lower_b - cs_int$upper_b) / (cs_int$upper_m - cs_int$lower_m) # Common x coordinate
cs_int$y <- cs_int$upper_m * cs_int$x + cs_int$upper_b # Common y coordinate
# Data needs to be ordered from smallest to largest since the subsequent function doesn't understand coordinates:
cs_int$e1_up_order <- ifelse(cs_int$e_lb_up > cs_int$e_rb_up,cs_int$e_rb_up,cs_int$e_lb_up)
cs_int$e2_up_order <- ifelse(cs_int$e_lb_up > cs_int$e_rb_up,cs_int$e_lb_up,cs_int$e_rb_up)
cs_int$n1_up_order <- ifelse(cs_int$n_lb_up > cs_int$n_rb_up,cs_int$n_rb_up,cs_int$n_lb_up)
cs_int$n2_up_order <- ifelse(cs_int$n_lb_up > cs_int$n_rb_up,cs_int$n_lb_up,cs_int$n_rb_up)

cs_int$e1_low_order <- ifelse(cs_int$e_lb_low > cs_int$e_rb_low,cs_int$e_rb_low,cs_int$e_lb_low)
cs_int$e2_low_order <- ifelse(cs_int$e_lb_low > cs_int$e_rb_low,cs_int$e_lb_low,cs_int$e_rb_low)
cs_int$n1_low_order <- ifelse(cs_int$n_lb_low > cs_int$n_rb_low,cs_int$n_rb_low,cs_int$n_lb_low)
cs_int$n2_low_order <- ifelse(cs_int$n_lb_low > cs_int$n_rb_low,cs_int$n_lb_low,cs_int$n_rb_low)

# Determine if the lines cross within the study area or at some point beyond the hubs:
findInt <- function(value, start, end) {
        start < value & end > value
}

cs_int$x_up_intersect <- findInt(cs_int$x, cs_int$e1_up_order, cs_int$e2_up_order)
cs_int$y_up_intersect <- findInt(cs_int$y, cs_int$n1_up_order, cs_int$n2_up_order)

cs_int$x_low_intersect <- findInt(cs_int$x, cs_int$e1_low_order, cs_int$e2_low_order)
cs_int$y_low_intersect <- findInt(cs_int$y, cs_int$n1_low_order, cs_int$n2_low_order)

# Only one condition actually needs to be true, but if only one is true then there's an error somewhere

cs_int$intersect <- ifelse(cs_int$y_up_intersect == TRUE & cs_int$x_up_intersect == TRUE & cs_int$y_low_intersect == TRUE & cs_int$x_low_intersect == TRUE,TRUE,FALSE)


# Clean up the df
cs_int$e1_up_order <- NULL
cs_int$e2_up_order <- NULL
cs_int$n1_up_order <- NULL
cs_int$n2_up_order <- NULL
cs_int$e1_low_order <- NULL
cs_int$e2_low_order <- NULL
cs_int$n1_low_order <- NULL
cs_int$n2_low_order <- NULL
cs_int$x_intersect <- NULL
cs_int$y_intersect <- NULL
# Generate a df populated w/ cross section data from cross sections that intersect within the study area
intersecting_cs <- filter(cs_int,intersect == TRUE)
intersecting_cs$intersect <- NULL
colnames(intersecting_cs)[colnames(intersecting_cs) == 'x'] <- 'x_intersect'
colnames(intersecting_cs)[colnames(intersecting_cs) == 'y'] <- 'y_intersect'

# Order the cross sections for plotting ----------------------------------------


# Calculate distance from the LB hub to the point of intersection:
intersecting_cs$lb_hub2inter <- sqrt((intersecting_cs$y_intersect - intersecting_cs$n_lb_up)^2 +  (intersecting_cs$x_intersect - intersecting_cs$e_lb_up)^2 )
# Subset the feature of interest
feat_raw <- filter(xs_data,attribs == "t")
# Remove any duplicate shots with the same x,y coordinates:
feat_raw <- arrange(feat_raw,CS,desc(hub))
bad <- duplicated(feat_raw$Easting,feat_raw$Northing,feat_raw$Elevation)
feat_raw <- feat_raw[!bad, ]
colnames(feat_raw)[colnames(feat_raw) == 'CS'] <- 'upper_cs'
feat <- left_join(feat_raw, intersecting_cs) # Join the feature to the order data
feat$lb_hub2pt <- sqrt((feat$Northing - feat$n_lb_up)^2 + (feat$Easting - feat$e_lb_up)^2) # Distance of the feature from LB hub
feat$in_order <- ifelse(feat$lb_hub2inter>feat$lb_hub2pt,TRUE,FALSE) # Is the feature on the left of the intersection?
feat_sub <- data.frame(feat$upper_cs,feat$lower_cs,feat$in_order)
feat_sub$feat.in_order <- ifelse(is.na(feat_sub$feat.in_order),TRUE,feat_sub$feat.in_order) # Convert "NA" to "TRUE"
out_order <- filter(feat_sub,feat_sub$feat.in_order == FALSE)
no_out_order <- length(out_order$feat.upper_cs) # number of cross sections out-of-order

# If features are out of order, then re-order them here:
if(no_out_order > 0) {
        out_order_cs <- c(out_order$feat.upper_cs,out_order$feat.lower_cs)
        out_order_cs <- unique(out_order_cs)
        out_order_cs <- sort(out_order_cs)
        # Count the number of times a cross section intersects another cross section in a given direction (upstream or downstream)
        no_inters_up_cs <- aggregate(data.frame(count = out_order$feat.upper_cs), list(value = out_order$feat.upper_cs), length)
        no_inters_do_cs <- aggregate(data.frame(count = out_order$feat.lower_cs), list(value = out_order$feat.lower_cs), length)
        # Join the two counts together
        net_no_intersections <- full_join(no_inters_up_cs,no_inters_do_cs, by = "value")
        net_no_intersections$count.x <- ifelse(is.na(net_no_intersections$count.x),0,net_no_intersections$count.x)
        net_no_intersections$count.y <- ifelse(is.na(net_no_intersections$count.y),0,net_no_intersections$count.y)
        names(net_no_intersections) <- c("cs_no","positive","negative")
        net_no_intersections$net_shift <- net_no_intersections$positive - net_no_intersections$negative
        net_no_intersections$new_order <- net_no_intersections$net_shift + net_no_intersections$cs_no
        colnames(feat_raw)[colnames(feat_raw) == 'upper_cs'] <- 'cs_no'
        feat_raw <- full_join(net_no_intersections,feat_raw)
        feat_raw$new_order <- ifelse(is.na(feat_raw$new_order),feat_raw$cs_no,feat_raw$new_order)
        feat_raw <- arrange(feat_raw,new_order)
        breakline <- feat_raw[,-c(1:8,12,14)]
}

if(no_out_order == 0) {
        breakline <- feat_raw[,-c(1:4,8,10)]
}

# Interpolate the breakline ----------------------------------------------------

names(breakline) <- c("X","Y","Z","desc")
interp <- auto.breakline(breakline,0.5, 0.1, 0.001)
bl_t <- as.data.frame(interp)






