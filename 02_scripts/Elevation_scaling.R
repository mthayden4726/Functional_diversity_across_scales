# Also make figures of scaling output
dat <- fread('03_output/TEAK_functional_richness_output.csv')
scale_plot <- ggplot(data = subset(dat, fric != "in progress"), aes(x = window, y = as.numeric(fric),
                                                                    color = file)) +
  geom_point(size = 4) +
  geom_smooth(method = lm, formula = y ~ log(x), se = FALSE) +
  theme_bw() +
  ylim(0, 400000) +
  labs(x = "Window Size)", y = "Average Functional Richness Across Scene",
       color = "Scene") +
  #annotate("text", x=2, y=380000, label="R2 = 0.64", size=12, color="black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26))
scale_plot

summary(lm(as.numeric(fric) ~ log(window), data = subset(dat, fric != "in progress")))

# Get scaling of elevation for TEAK
library(terra)
library(raster)


#first import all files in a single folder as a list 
rastlist <- list.files(path = "NEON_lidar-elev/NEON.D17.TEAK.DP3.30024.001.2019-06/", pattern='.tif', 
                       all.files=TRUE, full.names=FALSE)

#import all raster files in folder using lapply
setwd("NEON_lidar-elev/NEON.D17.TEAK.DP3.30024.001.2019-06/")
allrasters <- lapply(rastlist, raster)

#to check the index numbers of all imported raster list elements
allrasters

#to run a function on an individual raster e.g., plot 
plot(allrasters[[100]])

windows <- c(9, 27, 73, 153, 275, 393, 465)
means_all <- data.frame(rast = as.numeric(), window = as.numeric(), mean = as.numeric())

for (i in 1:length(allrasters)){
  rast <- allrasters[[i]]
  for (j in 1:length(windows)){
    r <- terra::rast(rast) # needs to be SpatRaster
    r_window <- focal(r, w = windows[j], fun = "mean", na.rm=FALSE) #the default is to not consider NAs
    means_loc <- data.frame(rast = rastlist[i], 
                            window = windows[j], 
                            mean = as.numeric(global(r_window, mean, na.rm = TRUE)))
    means_all <- rbind(means_all, means_loc)
  }
}