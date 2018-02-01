library("dplyr")

xs_data <- read.csv("/Users/sab/Dropbox/EFN\ Phase\ 2/Final\ recoded\ data/SA2/SA21971recode.csv", header = TRUE)
hubs <- filter(xs_data, hub == "LBH" | hub == "RBH")
hubs <- arrange(hubs, CS,hub,desc(Elevation)) # Sort the hubs

duplicate_hub <- duplicated(hubs$Easting,hubs$Northing) # Look for ground shot at hub (comes w/ lower Z but same X,Y)
hubs[duplicated(!duplicated(hubs$Easting,hubs$Northing)), ]
hubs <- hubs[-c(40,42),]

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

# Find the point where cross section lines intersect ---------------------------

combos <- combn(hubs$cs,2,simplify = TRUE)  # Generate an array of all combinations of cross sections (1 vs 2, 1 vs 3..., etc.)
upper_cs <- combos[1,]
lower_cs <- combos[2,]
cs_grid <- data.frame(upper_cs,lower_cs)
cs_grid$upper_cs <- as.numeric(cs_grid$upper_cs)
cs_grid$lower_cs <- as.numeric(cs_grid$lower_cs)
hubs$upper_cs <- as.numeric(hubs$cs)
intersect_upper_cs <- left_join(cs_grid, hubs)
hubs$lower_cs <- as.numeric(hubs$cs)
hubs$upper_cs <- NULL
intersect_lower_cs <- left_join(cs_grid, hubs)
cs_int <- data.frame(intersect_upper_cs$e_lb,intersect_upper_cs$n_lb,intersect_upper_cs$e_rb,intersect_upper_cs$n_rb,intersect_upper_cs$upper_cs,intersect_upper_cs$lower_cs,intersect_upper_cs$m,intersect_upper_cs$b,intersect_lower_cs$m,intersect_lower_cs$b)
names(cs_int) <- c("e_lb","n_lb","e_rb","n_rb","upper_cs","lower_cs","upper_m","upper_b","lower_m","lower_b")
# Find the point where the lines cross
cs_int$x <- (cs_int$lower_b - cs_int$upper_b) / (cs_int$upper_m - cs_int$lower_m)
cs_int$y <- cs_int$upper_m * cs_int$x + cs_int$upper_b
# Determine if the lines cross within the study area or at some point beyond the hubs:
cs_int$e1_order <- ifelse(cs_int$e_lb > cs_int$e_rb,cs_int$e_rb,cs_int$e_lb)
cs_int$e2_order <- ifelse(cs_int$e_lb > cs_int$e_rb,cs_int$e_lb,cs_int$e_rb)
cs_int$n1_order <- ifelse(cs_int$n_lb > cs_int$n_rb,cs_int$n_rb,cs_int$n_lb)
cs_int$n2_order <- ifelse(cs_int$n_lb > cs_int$n_rb,cs_int$n_lb,cs_int$n_rb)


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



