library(data.table)
library(ggplot2)
data <- fread("for_R.csv")

ggplot(data = data, aes(x = scale, y = fric, color = as.factor(scene))) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 18) +
  labs(x = "Window size (x by x meters)", y = "Average functional richness across scene",
       color = "Scene") +
  theme(legend.position = "top")

data_11 <- data[scene == "1.1", ]
summary(lm(fric ~ scale, data = data_11))

data_12 <- data[scene == "1.2", ]
summary(lm(fric ~ scale, data = data_12))

data_21 <- data[scene == "2.1", ]
summary(lm(fric ~ scale, data = data_21))

data_31 <- data[scene == "3.1", ]
summary(lm(fric ~ scale, data = data_31))

data_32 <- data[scene == "3.2", ]
summary(lm(fric ~ scale, data = data_32))

data_41 <- data[scene == "4.1", ]
summary(lm(fric ~ scale, data = data_41))

data_42 <- data[scene == "4.2", ]
summary(lm(fric ~ scale, data = data_42))

data_51 <- data[scene == "5.1", ]
summary(lm(fric ~ scale, data = data_51))

data_52 <- data[scene == "5.2", ]
summary(lm(fric ~ scale, data = data_52))
