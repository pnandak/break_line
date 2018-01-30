
# Code hubs as c(x_lb, y_lb,x_rb,y_rb)
xs_us <- c(2,2,4,4)
xs_ds <- c(4,2,2,4)

x_us_lb <- c(xs_us[1],xs_us[3])
y_us_lb <- c(xs_us[2],xs_us[4])
x_ds_lb <- c(xs_ds[1],xs_ds[3])
y_ds_lb <- c(xs_ds[2],xs_ds[4])


# Find the slopes
m_1 <- (xs_us[4]-xs_us[2])/(xs_us[3]-xs_us[1])
m_2 <- (xs_ds[4]-xs_ds[2])/(xs_ds[3]-xs_ds[1])

# Find the intercepts
b_1 <- xs_us[2] - m_1 * xs_us[1]
b_2 <- xs_ds[2] - m_2 * xs_ds[1]

# Find the point where the lines cross
x <- (b_2 - b_1) / (m_1 - m_2)
y <- m_1 * x + b_1

# Find the distace from the LB hub to the point where the lines cross
dist2pt <- sqrt((x - xs_us[1])^2 + (y - xs_us[2])^2)

xs_us_allx <- c(2.1,2.2,2.3,2.4,2.5,3.5,3.6,3.8)
xs_us_ally <- c(2.1,2.2,2.3,2.4,2.5,3.5,3.6,3.8)
xs_ds_allx <- c(2.1,2.2,2.3,2.4,2.5,3.5,3.6,3.8)
xs_ds_ally <- c(3.8,3.6,3.5,2.5,2.4,2.3,2.2,2.1)
xs_us_all <- cbind(xs_us_allx,xs_us_ally)
xs_ds_all <- cbind(xs_ds_allx,xs_ds_ally)

dist_hd <- sqrt((xs_us[1] - xs_us_all[,1])^2 + (xs_us[2]  - xs_us_all[,2])^2)
ifelse(dist_hd <= dist2pt,1,0)








# Plot the line
x_1p <- c(xs_us[1],xs_us[3])
y_1p <- c(xs_us[2],xs_us[4])
x_2p <- c(xs_ds[1],xs_ds[3])
y_2p <- c(xs_ds[2],xs_ds[4])

plot(x_1p,y_1p,type="l")
lines(x_2p,y_2p,type="l")
