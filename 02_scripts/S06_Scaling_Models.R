# Also make figures of scaling output
dat <- fread('03_output/TEAK_functional_richness_output1.csv')
scale_plot <- ggplot(data = subset(dat, fric != "in progress"), aes(x = window, y = as.numeric(fric),
                                                                    color = file)) +
  geom_point(size = 4) +
  geom_smooth(method = loess, formula = y ~ log(x), se = FALSE) +
  theme_bw() +
  facet_wrap(~file, scales = "free") +
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

dat2 <- fread('03_output/SRER_functional_richness_output1.csv')
scale_plot <- ggplot(data = subset(dat2, fric != "in progress"), aes(x = window, y = as.numeric(fric),
                                                                    color = file)) +
  geom_point(size = 4) +
  #geom_smooth(method = lm, formula = y ~ log(x), se = FALSE) +
  theme_bw() +
  facet_wrap(~file, scales = "free") +
  labs(x = "Window Size", y = "Average Functional Richness Across Scene",
       color = "Scene") +
  #annotate("text", x=2, y=380000, label="R2 = 0.64", size=12, color="black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26))
scale_plot

dat3 <- fread('03_output/SERC_functional_richness_output.csv')
scale_plot <- ggplot(data = subset(dat3, fric != "in progress"), aes(x = window, y = as.numeric(fric),
                                                                     color = file)) +
  geom_point(size = 4) +
  #geom_smooth(method = lm, formula = y ~ log(x), se = FALSE) +
  theme_bw() +
  facet_wrap(~file, scales = "free") +
  labs(x = "Window Size", y = "Average Functional Richness Across Scene",
       color = "Scene") +
  #annotate("text", x=2, y=380000, label="R2 = 0.64", size=12, color="black") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26))
scale_plot

summary(lm(as.numeric(fric) ~ log(window), data = subset(dat2, fric != "in progress")))

dat$site <- "TEAK"
dat$fric <- as.numeric(dat$fric)
dat2$site <- "SRER"
dat3$site <- "SERC"
dat_all <- rbind(dat, dat2, dat3, fill = TRUE)
dat_all <- dat_all[, c(1:4)]


ggplot(data = dat_all, aes(x = window, y = fric, color = site)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE) +
  facet_wrap(~site, scales = "free") +
  theme_bw(base_size = 20) +
  labs(x = "Window Size (x by x m)", y = "Functional Richness", 
       color = "Site") +
  theme(legend.position = "none")

means <- dat_all %>%
  group_by(site, window) %>%
  summarise(mean = mean(fric, na.rm = TRUE),
            sd = sd(fric, na.rm = TRUE))

means_dat <- dat %>%
  filter(is.numeric(fric)) %>%
  group_by(window) %>%
  summarise(mean = mean(fric, na.rm = TRUE),
            sd = sd(fric, na.rm = TRUE))

means_dat$site <- "TEAK"

means_dat2 <- dat2 %>%
  group_by(window) %>%
  summarise(mean = mean(fric, na.rm = TRUE),
            sd = sd(fric, na.rm = TRUE))

means_dat2$site <- "SRER"

means_dat3 <- dat3 %>%
  group_by(window) %>%
  summarise(mean = mean(fric, na.rm = TRUE),
            sd = sd(fric, na.rm = TRUE))

means_dat3$site <- "SERC"

means_all <- rbind(means_dat, means_dat2, means_dat3)

ggplot(data = means_all, aes(x = window, y = mean,
                             ymin = mean - sd,
                             ymax = mean + sd, 
                             color = site)) +
  geom_pointrange(size = 1, linewidth = 1) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE) +
  facet_wrap(~site, scales = "free") +
  theme_bw(base_size = 20) +
  labs(x = "Window Size (x by x m)", y = "Functional Richness", 
       color = "Site") +
  theme(legend.position = "none")

dat_sar <- dat3[, c(2,3)]
#dat_sar <- dat_sar[-10, ]
fit <- sar_loga(data = dat_sar, grid_start = "partial") 
summary(fit) 
plot(fit)
fitC <- sar_multi(data = dat_sar, obj = c("power", "loga", "monod"))
plot(fitC)

fit <- sar_average(data= dat_sar, obj =c("power","loga","koba","logistic","monod",
                                         "negexpo","chapman","weibull3","asymp"),
                   grid_start = "none", normaTest = "none", homoTest = "none", neg_check = FALSE, 
                   confInt = TRUE, ciN = 50) #a message is provided indicating that one model

par(mfrow = c(3,1)) #plot all model fits and the multimodel SAR curve as a separate curve on top
plot(fit, ModTitle = "a) Multimodel SAR", mmSep = TRUE)
summary(fit)
fit <- lin_pow(dat_sar, con = 1)
summary(fit)
plot(fit)

#plot the multimodel SAR curve (with confidence intervals; see explanation
#in the main text, above) on its own 
plot(fit, allCurves = FALSE, ModTitle =
       "c) Multimodel SAR with confidence intervals", confInt = TRUE)

#Barplot of the information criterion weights of each model 
plot(fit, type = "bar", ModTitle = "b) Model weights", cex.lab = 1.3)







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
