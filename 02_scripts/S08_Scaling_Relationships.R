library(data.table)
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(sars)
setwd("/Users/meha3816/Downloads")
## Processing sites -----
## BART ----

# FRIC
fric_bart_1 <- fread("_BART_fric_veg_015.csv")
fric_bart_2 <- fread("_BART_fric_veg_013.csv")
fric_bart_3 <- fread("_BART_fric_veg_026.csv")
fric_bart_4 <- fread("_BART_fric_veg_029.csv")
fric_bart_5 <- fread("_BART_fric_veg_012.csv")
fric_bart_6 <- fread("_BART_fric_veg_027.csv")

fric_bart_1$plot <- '015'
fric_bart_2$plot <- '013'
fric_bart_3$plot <- '026'
fric_bart_4$plot <- '029'
fric_bart_5$plot <- '012'
fric_bart_6$plot <- '027'

fric_bart_all <- rbind(fric_bart_1,
                       fric_bart_2,
                       fric_bart_3,
                       fric_bart_4,
                       fric_bart_5,
                       fric_bart_6)

fric_bart_all$site <- "bart"

fric_bart_sum <- fric_bart_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_bart_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_bart_1 <- fread("_BART_fdiv_veg_015.csv")
fdiv_bart_2 <- fread("_BART_fdiv_veg_013.csv")
fdiv_bart_3 <- fread("_BART_fdiv_veg_026.csv")
fdiv_bart_4 <- fread("_BART_fdiv_veg_029.csv")
fdiv_bart_5 <- fread("_BART_fdiv_veg_012.csv")
fdiv_bart_6 <- fread("_BART_fdiv_veg_027.csv")

fdiv_bart_1$plot <- '015'
fdiv_bart_2$plot <- '013'
fdiv_bart_3$plot <- '026'
fdiv_bart_4$plot <- '029'
fdiv_bart_5$plot <- '012'
fdiv_bart_6$plot <- '027'

fdiv_bart_all <- rbind(fdiv_bart_1,
                       fdiv_bart_2,
                       fdiv_bart_3,
                       fdiv_bart_4,
                       fdiv_bart_5,
                       fdiv_bart_6)

fdiv_bart_all$site <- "bart"

fdiv_bart_sum <- fdiv_bart_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_bart_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_bart_all, "FRic_BART.csv")
#write.csv(fdiv_bart_all, "FDiv_BART.csv")

## CLBJ ----

# FRIC
fric_clbj_1 <- fread("_CLBJ_fric_veg_003.csv")
fric_clbj_2 <- fread("_CLBJ_fric_veg_035.csv")
fric_clbj_3 <- fread("_CLBJ_fric_veg_041.csv")
fric_clbj_4 <- fread("_CLBJ_fric_veg_047.csv")
fric_clbj_5 <- fread("_CLBJ_fric_veg_048.csv")
fric_clbj_6 <- fread("_CLBJ_fric_veg_050.csv")
fric_clbj_7 <- fread("_CLBJ_fric_veg_051.csv")
fric_clbj_8 <- fread("_CLBJ_fric_veg_053.csv")
fric_clbj_9 <- fread("_CLBJ_fric_veg_054.csv")
fric_clbj_10 <- fread("_CLBJ_fric_veg_056.csv")

fric_clbj_1$plot <- '003'
fric_clbj_2$plot <- '035'
fric_clbj_3$plot <- '041'
fric_clbj_4$plot <- '047'
fric_clbj_5$plot <- '048'
fric_clbj_6$plot <- '050'
fric_clbj_7$plot <- '051'
fric_clbj_8$plot <- '053'
fric_clbj_9$plot <- '054'
fric_clbj_10$plot <- '056'

fric_clbj_all <- rbind(fric_clbj_1,
                       fric_clbj_2,
                       fric_clbj_3,
                       fric_clbj_4,
                       fric_clbj_5,
                       fric_clbj_6,
                       fric_clbj_7,
                       fric_clbj_8,
                       fric_clbj_9,
                       fric_clbj_10)

fric_clbj_all$site <- "clbj"

fric_clbj_sum <- fric_clbj_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_clbj_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_clbj_1 <- fread("_CLBJ_fdiv_veg_003.csv")
fdiv_clbj_2 <- fread("_CLBJ_fdiv_veg_035.csv")
fdiv_clbj_3 <- fread("_CLBJ_fdiv_veg_041.csv")
fdiv_clbj_4 <- fread("_CLBJ_fdiv_veg_047.csv")
fdiv_clbj_5 <- fread("_CLBJ_fdiv_veg_048.csv")
fdiv_clbj_6 <- fread("_CLBJ_fdiv_veg_050.csv")
fdiv_clbj_7 <- fread("_CLBJ_fdiv_veg_051.csv")
fdiv_clbj_8 <- fread("_CLBJ_fdiv_veg_053.csv")
fdiv_clbj_9 <- fread("_CLBJ_fdiv_veg_054.csv")
fdiv_clbj_10 <- fread("_CLBJ_fdiv_veg_056.csv")

fdiv_clbj_1$plot <- '003'
fdiv_clbj_2$plot <- '035'
fdiv_clbj_3$plot <- '041'
fdiv_clbj_4$plot <- '047'
fdiv_clbj_5$plot <- '048'
fdiv_clbj_6$plot <- '050'
fdiv_clbj_7$plot <- '051'
fdiv_clbj_8$plot <- '053'
fdiv_clbj_9$plot <- '054'
fdiv_clbj_10$plot <- '056'

fdiv_clbj_all <- rbind(fdiv_clbj_1,
                       fdiv_clbj_2,
                       fdiv_clbj_3,
                       fdiv_clbj_4,
                       fdiv_clbj_5,
                       fdiv_clbj_6,
                       fdiv_clbj_7,
                       fdiv_clbj_8,
                       fdiv_clbj_9,
                       fdiv_clbj_10)

fdiv_clbj_all$site <- "clbj"

fdiv_clbj_sum <- fdiv_clbj_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_clbj_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_clbj_all, "FRic_CLBJ.csv")
#write.csv(fdiv_clbj_all, "FDiv_CLBJ.csv")

## HEAL -------
# FRIC
fric_heal_1 <- fread("_HEAL_fric_veg_002_v2.csv")
fric_heal_2 <- fread("_HEAL_fric_veg_004_v2.csv")
fric_heal_3 <- fread("_HEAL_fric_veg_005_v2.csv")
fric_heal_4 <- fread("_HEAL_fric_veg_013_v2.csv")
fric_heal_5 <- fread("_HEAL_fric_veg_015_v2.csv")
fric_heal_6 <- fread("_HEAL_fric_veg_018_v2.csv")
fric_heal_7 <- fread("_HEAL_fric_veg_024_v2.csv")
fric_heal_8 <- fread("_HEAL_fric_veg_026_v2.csv")

fric_heal_1$plot <- '002'
fric_heal_2$plot <- '004'
fric_heal_3$plot <- '005'
fric_heal_4$plot <- '013'
fric_heal_5$plot <- '015'
fric_heal_6$plot <- '018'
fric_heal_7$plot <- '024'
fric_heal_8$plot <- '026'

fric_heal_all <- rbind(fric_heal_1,
                       fric_heal_2,
                       fric_heal_3,
                       fric_heal_4,
                       fric_heal_5,
                       fric_heal_6,
                       fric_heal_7,
                       fric_heal_8)

fric_heal_all$site <- "HEAL"

fric_heal_sum <- fric_heal_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_heal_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_heal_1 <- fread("_HEAL_fdiv_veg_002_v2.csv")
fdiv_heal_2 <- fread("_HEAL_fdiv_veg_004_v2.csv")
fdiv_heal_3 <- fread("_HEAL_fdiv_veg_005_v2.csv")
fdiv_heal_4 <- fread("_HEAL_fdiv_veg_013_v2.csv")
fdiv_heal_5 <- fread("_HEAL_fdiv_veg_015_v2.csv")
fdiv_heal_6 <- fread("_HEAL_fdiv_veg_018_v2.csv")
fdiv_heal_7 <- fread("_HEAL_fdiv_veg_024_v2.csv")
fdiv_heal_8 <- fread("_HEAL_fdiv_veg_026_v2.csv")

fdiv_heal_1$plot <- '002'
fdiv_heal_2$plot <- '004'
fdiv_heal_3$plot <- '005'
fdiv_heal_4$plot <- '013'
fdiv_heal_5$plot <- '015'
fdiv_heal_6$plot <- '018'
fdiv_heal_7$plot <- '024'
fdiv_heal_8$plot <- '026'

fdiv_heal_all <- rbind(fdiv_heal_1,
                       fdiv_heal_2,
                       fdiv_heal_3,
                       fdiv_heal_4,
                       fdiv_heal_5,
                       fdiv_heal_6,
                       fdiv_heal_7,
                       fdiv_heal_8)

fdiv_heal_all$site <- "HEAL"

fdiv_heal_sum <- fdiv_heal_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_heal_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_heal_all, "FRic_HEAL.csv")
#write.csv(fdiv_heal_all, "FDiv_HEAL.csv")

## KONZ -----
# FRIC
fric_KONZ_1 <- fread("_KONZ_fric_veg_005.csv")
fric_KONZ_2 <- fread("_KONZ_fric_veg_007.csv")
fric_KONZ_3 <- fread("_KONZ_fric_veg_009.csv")
fric_KONZ_4 <- fread("_KONZ_fric_veg_019.csv")

fric_KONZ_1$plot <- '005'
fric_KONZ_2$plot <- '007'
fric_KONZ_3$plot <- '009'
fric_KONZ_4$plot <- '019'


fric_KONZ_all <- rbind(fric_KONZ_1,
                       fric_KONZ_2,
                       fric_KONZ_3,
                       fric_KONZ_4)

fric_KONZ_all$site <- "KONZ"

fric_KONZ_sum <- fric_KONZ_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_KONZ_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_KONZ_1 <- fread("_KONZ_fdiv_veg_005.csv")
fdiv_KONZ_2 <- fread("_KONZ_fdiv_veg_007.csv")
fdiv_KONZ_3 <- fread("_KONZ_fdiv_veg_009.csv")
fdiv_KONZ_4 <- fread("_KONZ_fdiv_veg_019.csv")

fdiv_KONZ_1$plot <- '005'
fdiv_KONZ_2$plot <- '007'
fdiv_KONZ_3$plot <- '009'
fdiv_KONZ_4$plot <- '019'


fdiv_KONZ_all <- rbind(fdiv_KONZ_1,
                       fdiv_KONZ_2,
                       fdiv_KONZ_3,
                       fdiv_KONZ_4)

fdiv_KONZ_all$site <- "KONZ"

fdiv_KONZ_sum <- fdiv_KONZ_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_KONZ_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_KONZ_all, "FRic_KONZ.csv")
#write.csv(fdiv_KONZ_all, "FDiv_KONZ.csv")

## NIWO ----
# FRIC
fric_niwo_1 <- fread("_NIWO_fric_veg_004.csv")
fric_niwo_2 <- fread("_NIWO_fric_veg_006.csv")
fric_niwo_3 <- fread("_NIWO_fric_veg_007.csv")
fric_niwo_4 <- fread("_NIWO_fric_veg_016.csv")
fric_niwo_5 <- fread("_NIWO_fric_veg_021.csv")
fric_niwo_6 <- fread("_NIWO_fric_veg_030.csv")

fric_niwo_1$plot <- '004'
fric_niwo_2$plot <- '006'
fric_niwo_3$plot <- '007'
fric_niwo_4$plot <- '016'
fric_niwo_5$plot <- '021'
fric_niwo_6$plot <- '030'

fric_niwo_all <- rbind(fric_niwo_1,
                       fric_niwo_2,
                       fric_niwo_3,
                       fric_niwo_4,
                       fric_niwo_5,
                       fric_niwo_6)

fric_niwo_all$site <- "NIWO"

fric_niwo_sum <- fric_niwo_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_niwo_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_niwo_1 <- fread("_NIWO_fdiv_veg_004.csv")
fdiv_niwo_2 <- fread("_NIWO_fdiv_veg_006.csv")
fdiv_niwo_3 <- fread("_NIWO_fdiv_veg_007.csv")
fdiv_niwo_4 <- fread("_NIWO_fdiv_veg_016.csv")
fdiv_niwo_5 <- fread("_NIWO_fdiv_veg_021.csv")
fdiv_niwo_6 <- fread("_NIWO_fdiv_veg_030.csv")

fdiv_niwo_1$plot <- '004'
fdiv_niwo_2$plot <- '006'
fdiv_niwo_3$plot <- '007'
fdiv_niwo_4$plot <- '016'
fdiv_niwo_5$plot <- '021'
fdiv_niwo_6$plot <- '030'

fdiv_niwo_all <- rbind(fdiv_niwo_1,
                       fdiv_niwo_2,
                       fdiv_niwo_3,
                       fdiv_niwo_4,
                       fdiv_niwo_5,
                       fdiv_niwo_6)

fdiv_niwo_all$site <- "NIWO"

fdiv_niwo_sum <- fdiv_niwo_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_niwo_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_niwo_all, "FRic_NIWO.csv")
#write.csv(fdiv_niwo_all, "FDiv_NIWO.csv")

## ONAQ --------

fric_onaq_2 <- fread("_ONAQ_fric_veg_005_v2.csv")
fric_onaq_3 <- fread("_ONAQ_fric_veg_008_v2.csv")
fric_onaq_4 <- fread("_ONAQ_fric_veg_010_v2.csv")
fric_onaq_5 <- fread("_ONAQ_fric_veg_011_v2.csv")
fric_onaq_6 <- fread("_ONAQ_fric_veg_018_v2.csv")
fric_onaq_7 <- fread("_ONAQ_fric_veg_019_v2.csv")
fric_onaq_8 <- fread("_ONAQ_fric_veg_021_v2.csv")
fric_onaq_9 <- fread("_ONAQ_fric_veg_024_v2.csv")
fric_onaq_10 <- fread("_ONAQ_fric_veg_030_v2.csv")
fric_onaq_11 <- fread("_ONAQ_fric_veg_043_v2.csv")
fric_onaq_12 <- fread("_ONAQ_fric_veg_073_v2.csv")

fric_onaq_2$plot <- "005"
fric_onaq_3$plot <- "008"
fric_onaq_4$plot <- "010"
fric_onaq_5$plot <- "011"
fric_onaq_6$plot <- "018"
fric_onaq_7$plot <- "019"
fric_onaq_8$plot <- "021"
fric_onaq_9$plot <- "024"
fric_onaq_10$plot <- "030"
fric_onaq_11$plot <- "043"
fric_onaq_12$plot <- "073"

fric_onaq_all <- rbind(
                       fric_onaq_2,
                       fric_onaq_3,
                       fric_onaq_4,
                       fric_onaq_5,
                       fric_onaq_6,
                       fric_onaq_7,
                       fric_onaq_8,
                       fric_onaq_9,
                       fric_onaq_10,
                       fric_onaq_11,
                       fric_onaq_12)

fric_onaq_all$site <- "ONAQ"

fric_onaq_sum <- fric_onaq_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_onaq_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

fdiv_onaq_2 <- fread("_ONAQ_fdiv_veg_005_v2.csv")
fdiv_onaq_3 <- fread("_ONAQ_fdiv_veg_008_v2.csv")
fdiv_onaq_4 <- fread("_ONAQ_fdiv_veg_010_v2.csv")
fdiv_onaq_5 <- fread("_ONAQ_fdiv_veg_011_v2.csv")
fdiv_onaq_6 <- fread("_ONAQ_fdiv_veg_018_v2.csv")
fdiv_onaq_7 <- fread("_ONAQ_fdiv_veg_019_v2.csv")
fdiv_onaq_8 <- fread("_ONAQ_fdiv_veg_021_v2.csv")
fdiv_onaq_9 <- fread("_ONAQ_fdiv_veg_024_v2.csv")
fdiv_onaq_10 <- fread("_ONAQ_fdiv_veg_030_v2.csv")
fdiv_onaq_11 <- fread("_ONAQ_fdiv_veg_043_v2.csv")
fdiv_onaq_12 <- fread("_ONAQ_fdiv_veg_073_v2.csv")

fdiv_onaq_2$plot <- "005"
fdiv_onaq_3$plot <- "008"
fdiv_onaq_4$plot <- "010"
fdiv_onaq_5$plot <- "011"
fdiv_onaq_6$plot <- "018"
fdiv_onaq_7$plot <- "019"
fdiv_onaq_8$plot <- "021"
fdiv_onaq_9$plot <- "024"
fdiv_onaq_10$plot <- "030"
fdiv_onaq_11$plot <- "043"
fdiv_onaq_12$plot <- "073"

fdiv_onaq_all <- rbind(
                       fdiv_onaq_2,
                       fdiv_onaq_3,
                       fdiv_onaq_4,
                       fdiv_onaq_5,
                       fdiv_onaq_6,
                       fdiv_onaq_7,
                       fdiv_onaq_8,
                       fdiv_onaq_9,
                       fdiv_onaq_10,
                       fdiv_onaq_11,
                       fdiv_onaq_12)

fdiv_onaq_all$site <- "ONAQ"

fdiv_onaq_sum <- fdiv_onaq_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_onaq_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median fdiv +/- sd")


#write.csv(fric_onaq_all, "FRic_ONAQ.csv")
#write.csv(fdiv_onaq_all, "FDiv_ONAQ.csv")

## OSBS ----
# FRIC
fric_osbs_1 <- fread("_OSBS_fric_veg_002.csv")
fric_osbs_2 <- fread("_OSBS_fric_veg_005.csv")
fric_osbs_3 <- fread("_OSBS_fric_veg_007.csv")
fric_osbs_4 <- fread("_OSBS_fric_veg_010.csv")
fric_osbs_5 <- fread("_OSBS_fric_veg_011.csv")
fric_osbs_6 <- fread("_OSBS_fric_veg_027.csv")
fric_osbs_7 <- fread("_OSBS_fric_veg_048.csv")
fric_osbs_8 <- fread("_OSBS_fric_veg_051.csv")

fric_osbs_1$plot <- '002'
fric_osbs_2$plot <- '005'
fric_osbs_3$plot <- '007'
fric_osbs_4$plot <- '010'
fric_osbs_5$plot <- '011'
fric_osbs_6$plot <- '027'
fric_osbs_7$plot <- '048'
fric_osbs_8$plot <- '051'

fric_osbs_all <- rbind(fric_osbs_1,
                       fric_osbs_2,
                       fric_osbs_3,
                       fric_osbs_4,
                       fric_osbs_5,
                       fric_osbs_6,
                       fric_osbs_7,
                       fric_osbs_8)

fric_osbs_all$site <- "osbs"

fric_osbs_sum <- fric_osbs_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_osbs_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_osbs_1 <- fread("_OSBS_fdiv_veg_002.csv")
fdiv_osbs_2 <- fread("_OSBS_fdiv_veg_005.csv")
fdiv_osbs_3 <- fread("_OSBS_fdiv_veg_007.csv")
fdiv_osbs_4 <- fread("_OSBS_fdiv_veg_010.csv")
fdiv_osbs_5 <- fread("_OSBS_fdiv_veg_011.csv")
fdiv_osbs_6 <- fread("_OSBS_fdiv_veg_027.csv")
fdiv_osbs_7 <- fread("_OSBS_fdiv_veg_048.csv")
fdiv_osbs_8 <- fread("_OSBS_fdiv_veg_051.csv")

fdiv_osbs_1$plot <- '002'
fdiv_osbs_2$plot <- '005'
fdiv_osbs_3$plot <- '007'
fdiv_osbs_4$plot <- '010'
fdiv_osbs_5$plot <- '011'
fdiv_osbs_6$plot <- '027'
fdiv_osbs_7$plot <- '048'
fdiv_osbs_8$plot <- '051'

fdiv_osbs_all <- rbind(fdiv_osbs_1,
                       fdiv_osbs_2,
                       fdiv_osbs_3,
                       fdiv_osbs_4,
                       fdiv_osbs_5,
                       fdiv_osbs_6,
                       fdiv_osbs_7,
                       fdiv_osbs_8)

fdiv_osbs_all$site <- "osbs"

fdiv_osbs_sum <- fdiv_osbs_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_osbs_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_osbs_all, "FRic_OSBS.csv")
#write.csv(fdiv_osbs_all, "FDiv_OSBS.csv")
## PUUM ----
# FRIC
fric_PUUM_1 <- fread("_PUUM_fric_veg_005.csv")
fric_PUUM_2 <- fread("_PUUM_fric_veg_020.csv")
fric_PUUM_3 <- fread("_PUUM_fric_veg_017.csv")
fric_PUUM_4 <- fread("_PUUM_fric_veg_004.csv")
fric_PUUM_5 <- fread("_PUUM_fric_veg_013.csv")
fric_PUUM_6 <- fread("_PUUM_fric_veg_032.csv")

fric_PUUM_1$plot <- '005'
fric_PUUM_2$plot <- '020'
fric_PUUM_3$plot <- '017'
fric_PUUM_4$plot <- '004'
fric_PUUM_5$plot <- '013'
fric_PUUM_6$plot <- '032'

fric_PUUM_all <- rbind(fric_PUUM_1,
                       fric_PUUM_2,
                       fric_PUUM_3,
                       fric_PUUM_4,
                       fric_PUUM_5,
                       fric_PUUM_6)

fric_PUUM_all$site <- "PUUM"

fric_PUUM_sum <- fric_PUUM_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_PUUM_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_PUUM_1 <- fread("_PUUM_fdiv_veg_005.csv")
fdiv_PUUM_2 <- fread("_PUUM_fdiv_veg_020.csv")
fdiv_PUUM_3 <- fread("_PUUM_fdiv_veg_017.csv")
fdiv_PUUM_4 <- fread("_PUUM_fdiv_veg_004.csv")
fdiv_PUUM_5 <- fread("_PUUM_fdiv_veg_013.csv")
fdiv_PUUM_6 <- fread("_PUUM_fdiv_veg_032.csv")

fdiv_PUUM_1$plot <- '005'
fdiv_PUUM_2$plot <- '020'
fdiv_PUUM_3$plot <- '017'
fdiv_PUUM_4$plot <- '004'
fdiv_PUUM_5$plot <- '013'
fdiv_PUUM_6$plot <- '032'

fdiv_PUUM_all <- rbind(fdiv_PUUM_1,
                       fdiv_PUUM_2,
                       fdiv_PUUM_3,
                       fdiv_PUUM_4,
                       fdiv_PUUM_5,
                       fdiv_PUUM_6)

fdiv_PUUM_all$site <- "PUUM"

fdiv_PUUM_sum <- fdiv_PUUM_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_PUUM_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_PUUM_all, "FRic_PUUM.csv")
#write.csv(fdiv_PUUM_all, "FDiv_PUUM.csv")

## SERC ----
# FRIC
fric_SERC_1 <- fread("_SERC_fric_veg_044.csv")
fric_SERC_2 <- fread("_SERC_fric_veg_012.csv")
fric_SERC_3 <- fread("_SERC_fric_veg_010.csv")
fric_SERC_4 <- fread("_SERC_fric_veg_009.csv")
fric_SERC_5 <- fread("_SERC_fric_veg_005.csv")
fric_SERC_6 <- fread("_SERC_fric_veg_004.csv")
fric_SERC_7 <- fread("_SERC_fric_veg_001.csv")

fric_SERC_1$plot <- '044'
fric_SERC_2$plot <- '012'
fric_SERC_3$plot <- '010'
fric_SERC_4$plot <- '009'
fric_SERC_5$plot <- '005'
fric_SERC_6$plot <- '004'
fric_SERC_7$plot <- '001'

fric_SERC_all <- rbind(fric_SERC_1,
                       fric_SERC_2,
                       fric_SERC_3,
                       fric_SERC_4,
                       fric_SERC_5,
                       fric_SERC_6,
                       fric_SERC_7)

fric_SERC_all$site <- "SERC"

fric_SERC_sum <- fric_SERC_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_SERC_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_SERC_1 <- fread("_SERC_fdiv_veg_044.csv")
fdiv_SERC_2 <- fread("_SERC_fdiv_veg_012.csv")
fdiv_SERC_3 <- fread("_SERC_fdiv_veg_010.csv")
fdiv_SERC_4 <- fread("_SERC_fdiv_veg_009.csv")
fdiv_SERC_5 <- fread("_SERC_fdiv_veg_005.csv")
fdiv_SERC_6 <- fread("_SERC_fdiv_veg_004.csv")
fdiv_SERC_7 <- fread("_SERC_fdiv_veg_001.csv")

fdiv_SERC_1$plot <- '044'
fdiv_SERC_2$plot <- '012'
fdiv_SERC_3$plot <- '010'
fdiv_SERC_4$plot <- '009'
fdiv_SERC_5$plot <- '005'
fdiv_SERC_6$plot <- '004'
fdiv_SERC_7$plot <- '001'

fdiv_SERC_all <- rbind(fdiv_SERC_1,
                       fdiv_SERC_2,
                       fdiv_SERC_3,
                       fdiv_SERC_4,
                       fdiv_SERC_5,
                       fdiv_SERC_6,
                       fdiv_SERC_7)

fdiv_SERC_all$site <- "SERC"

fdiv_SERC_sum <- fdiv_SERC_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_SERC_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd,
                                 color = plot)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

write.csv(fric_SERC_all, "FRic_SERC.csv")
write.csv(fdiv_SERC_all, "FDiv_SERC.csv")

## SRER ----
# FRIC
fric_SRER_1 <- fread("_SRER_fric_veg_002_EPSG.csv")
fric_SRER_2 <- fread("_SRER_fric_veg_003_EPSG.csv")
fric_SRER_3 <- fread("_SRER_fric_veg_006_EPSG.csv")
fric_SRER_4 <- fread("_SRER_fric_veg_014_EPSG.csv")
fric_SRER_5 <- fread("_SRER_fric_veg_021_EPSG.csv")
fric_SRER_6 <- fread("_SRER_fric_veg_023_EPSG.csv")
fric_SRER_7 <- fread("_SRER_fric_veg_026_EPSG.csv")
fric_SRER_8 <- fread("_SRER_fric_veg_027_EPSG.csv")
fric_SRER_9 <- fread("_SRER_fric_veg_028_EPSG.csv")

fric_SRER_1$plot <- '002'
fric_SRER_2$plot <- '003'
fric_SRER_3$plot <- '006'
fric_SRER_4$plot <- '014'
fric_SRER_5$plot <- '021'
fric_SRER_6$plot <- '023'
fric_SRER_7$plot <- '026'
fric_SRER_8$plot <- '027'
fric_SRER_9$plot <- '028'


fric_SRER_all <- rbind(fric_SRER_1,
                       fric_SRER_2,
                       fric_SRER_3,
                       fric_SRER_4,
                       fric_SRER_5,
                       fric_SRER_6,
                       fric_SRER_7,
                       fric_SRER_8,
                       fric_SRER_9)

fric_SRER_all$site <- "SRER"

fric_SRER_sum <- fric_SRER_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_SRER_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_SRER_1 <- fread("_SRER_fdiv_veg_002_EPSG.csv")
fdiv_SRER_2 <- fread("_SRER_fdiv_veg_003_EPSG.csv")
fdiv_SRER_3 <- fread("_SRER_fdiv_veg_006_EPSG.csv")
fdiv_SRER_4 <- fread("_SRER_fdiv_veg_014_EPSG.csv")
fdiv_SRER_5 <- fread("_SRER_fdiv_veg_021_EPSG.csv")
fdiv_SRER_6 <- fread("_SRER_fdiv_veg_023_EPSG.csv")
fdiv_SRER_7 <- fread("_SRER_fdiv_veg_026_EPSG.csv")
fdiv_SRER_8 <- fread("_SRER_fdiv_veg_027_EPSG.csv")
fdiv_SRER_9 <- fread("_SRER_fdiv_veg_028_EPSG.csv")

fdiv_SRER_1$plot <- '002'
fdiv_SRER_2$plot <- '003'
fdiv_SRER_3$plot <- '006'
fdiv_SRER_4$plot <- '014'
fdiv_SRER_5$plot <- '021'
fdiv_SRER_6$plot <- '023'
fdiv_SRER_7$plot <- '026'
fdiv_SRER_8$plot <- '027'
fdiv_SRER_9$plot <- '028'

fdiv_SRER_all <- rbind(fdiv_SRER_1,
                       fdiv_SRER_2,
                       fdiv_SRER_3,
                       fdiv_SRER_4,
                       fdiv_SRER_5,
                       fdiv_SRER_6,
                       fdiv_SRER_7,
                       fdiv_SRER_8,
                       fdiv_SRER_9)

fdiv_SRER_all$site <- "SRER"

fdiv_SRER_sum <- fdiv_SRER_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_SRER_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_SRER_all, "FRic_SRER.csv")
#write.csv(fdiv_SRER_all, "FDiv_SRER.csv")

## TALL -----
# FRIC
fric_tall_1 <- fread("_TALL_fric_veg_001.csv")
fric_tall_2 <- fread("_TALL_fric_veg_002.csv")
fric_tall_3 <- fread("_TALL_fric_veg_009.csv")
fric_tall_4 <- fread("_TALL_fric_veg_012.csv")
fric_tall_5 <- fread("_TALL_fric_veg_013.csv")
fric_tall_6 <- fread("_TALL_fric_veg_020.csv")
fric_tall_7 <- fread("_TALL_fric_veg_026.csv")
fric_tall_8 <- fread("_TALL_fric_veg_027.csv")
fric_tall_9 <- fread("_TALL_fric_veg_032.csv")
fric_tall_10 <- fread("_TALL_fric_veg_044.csv")

fric_tall_1$plot <- "001"
fric_tall_2$plot <- "002"
fric_tall_3$plot <- "009"
fric_tall_4$plot <- "012"
fric_tall_5$plot <- "013"
fric_tall_6$plot <- "020"
fric_tall_7$plot <- "026"
fric_tall_8$plot <- "027"
fric_tall_9$plot <- "032"
fric_tall_10$plot <- "044"

fric_tall_all <- rbind(fric_tall_1,
                       fric_tall_2,
                       fric_tall_3,
                       fric_tall_4,
                       fric_tall_5,
                       fric_tall_6,
                       fric_tall_7,
                       fric_tall_8,
                       fric_tall_9,
                       fric_tall_10
                       )

fric_tall_all$site <- "TALL"

fric_tall_sum <- fric_tall_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_tall_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_tall_1 <- fread("_TALL_fdiv_veg_001.csv")
fdiv_tall_2 <- fread("_TALL_fdiv_veg_002.csv")
fdiv_tall_3 <- fread("_TALL_fdiv_veg_009.csv")
fdiv_tall_4 <- fread("_TALL_fdiv_veg_012.csv")
fdiv_tall_5 <- fread("_TALL_fdiv_veg_013.csv")
fdiv_tall_6 <- fread("_TALL_fdiv_veg_020.csv")
fdiv_tall_7 <- fread("_TALL_fdiv_veg_026.csv")
fdiv_tall_8 <- fread("_TALL_fdiv_veg_027.csv")
fdiv_tall_9 <- fread("_TALL_fdiv_veg_032.csv")
fdiv_tall_10 <- fread("_TALL_fdiv_veg_044.csv")

fdiv_tall_1$plot <- "001"
fdiv_tall_2$plot <- "002"
fdiv_tall_3$plot <- "009"
fdiv_tall_4$plot <- "012"
fdiv_tall_5$plot <- "013"
fdiv_tall_6$plot <- "020"
fdiv_tall_7$plot <- "026"
fdiv_tall_8$plot <- "027"
fdiv_tall_9$plot <- "032"
fdiv_tall_10$plot <- "044"

fdiv_tall_all <- rbind(fdiv_tall_1,
                       fdiv_tall_2,
                       fdiv_tall_3,
                       fdiv_tall_4,
                       fdiv_tall_5,
                       fdiv_tall_6,
                       fdiv_tall_7,
                       fdiv_tall_8,
                       fdiv_tall_9,
                       fdiv_tall_10
)

fdiv_tall_all$site <- "TALL"

fdiv_tall_sum <- fdiv_tall_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_tall_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_tall_all, "FRic_TALL.csv")
#write.csv(fdiv_tall_all, "FDiv_TALL.csv")
## TALL v3 -----
# FRIC
fric_tall_1 <- fread("_TALL_fric_veg_001_v3.csv")
fric_tall_2 <- fread("_TALL_fric_veg_002_v3.csv")
fric_tall_3 <- fread("_TALL_fric_veg_009_v3.csv")
fric_tall_4 <- fread("_TALL_fric_veg_012_v3.csv")
fric_tall_5 <- fread("_TALL_fric_veg_013_v3.csv")
fric_tall_6 <- fread("_TALL_fric_veg_020_v3.csv")
fric_tall_7 <- fread("_TALL_fric_veg_026_v3.csv")
fric_tall_8 <- fread("_TALL_fric_veg_027_v3.csv")
fric_tall_9 <- fread("_TALL_fric_veg_032_v3.csv")
fric_tall_10 <- fread("_TALL_fric_veg_044_v3.csv")

fric_tall_1$plot <- "001"
fric_tall_2$plot <- "002"
fric_tall_3$plot <- "009"
fric_tall_4$plot <- "012"
fric_tall_5$plot <- "013"
fric_tall_6$plot <- "020"
fric_tall_7$plot <- "026"
fric_tall_8$plot <- "027"
fric_tall_9$plot <- "032"
fric_tall_10$plot <- "044"

fric_tall_all <- rbind(fric_tall_1,
                       fric_tall_2,
                       fric_tall_3,
                       fric_tall_4,
                       fric_tall_5,
                       fric_tall_6,
                       fric_tall_7,
                       fric_tall_8,
                       fric_tall_9,
                       fric_tall_10
)

fric_tall_all$site <- "TALL"

fric_tall_sum <- fric_tall_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_tall_sum1, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_tall_1 <- fread("_TALL_fdiv_veg_001_v3.csv")
fdiv_tall_2 <- fread("_TALL_fdiv_veg_002_v3.csv")
fdiv_tall_3 <- fread("_TALL_fdiv_veg_009_v3.csv")
fdiv_tall_4 <- fread("_TALL_fdiv_veg_012_v3.csv")
fdiv_tall_5 <- fread("_TALL_fdiv_veg_013_v3.csv")
fdiv_tall_6 <- fread("_TALL_fdiv_veg_020_v3.csv")
fdiv_tall_7 <- fread("_TALL_fdiv_veg_026_v3.csv")
fdiv_tall_8 <- fread("_TALL_fdiv_veg_027_v3.csv")
fdiv_tall_9 <- fread("_TALL_fdiv_veg_032_v3.csv")
fdiv_tall_10 <- fread("_TALL_fdiv_veg_044_v3.csv")

fdiv_tall_1$plot <- "001"
fdiv_tall_2$plot <- "002"
fdiv_tall_3$plot <- "009"
fdiv_tall_4$plot <- "012"
fdiv_tall_5$plot <- "013"
fdiv_tall_6$plot <- "020"
fdiv_tall_7$plot <- "026"
fdiv_tall_8$plot <- "027"
fdiv_tall_9$plot <- "032"
fdiv_tall_10$plot <- "044"

fdiv_tall_all <- rbind(fdiv_tall_1,
                       fdiv_tall_2,
                       fdiv_tall_3,
                       fdiv_tall_4,
                       fdiv_tall_5,
                       fdiv_tall_6,
                       fdiv_tall_7,
                       fdiv_tall_8,
                       fdiv_tall_9,
                       fdiv_tall_10
)

fdiv_tall_all$site <- "TALL"

fdiv_tall_sum <- fdiv_tall_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_tall_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_tall_all, "FRic_TALL.csv")
#write.csv(fdiv_tall_all, "FDiv_TALL.csv")
## TEAK ----
# FRIC
fric_TEAK_1 <- fread("_TEAK_fric_veg_site_1.csv")
fric_TEAK_2 <- fread("_TEAK_fric_veg_site_2.csv")
fric_TEAK_3 <- fread("_TEAK_fric_veg_site_3.csv")
fric_TEAK_4 <- fread("_TEAK_fric_veg_site_4.csv")
fric_TEAK_5 <- fread("_TEAK_fric_veg_site_10.csv")
fric_TEAK_6 <- fread("_TEAK_fric_veg_site_12.csv")
fric_TEAK_7 <- fread("_TEAK_fric_veg_site_5.csv")
fric_TEAK_8 <- fread("_TEAK_fric_veg_site_7.csv")
fric_TEAK_9 <- fread("_TEAK_fric_veg_site_9.csv")
fric_TEAK_10 <- fread("_TEAK_fric_veg_site_0.csv")

fric_TEAK_1$plot <- '001'
fric_TEAK_2$plot <- '002'
fric_TEAK_3$plot <- '003'
fric_TEAK_4$plot <- '004'
fric_TEAK_5$plot <- '010'
fric_TEAK_6$plot <- '012'
fric_TEAK_7$plot <- '005'
fric_TEAK_8$plot <- '007'
fric_TEAK_9$plot <- '009'
fric_TEAK_10$plot <- '000'

fric_TEAK_all <- rbind(fric_TEAK_1,
                       fric_TEAK_2,
                       fric_TEAK_3,
                       fric_TEAK_4,
                       fric_TEAK_5,
                       fric_TEAK_6,
                       fric_TEAK_7,
                       fric_TEAK_8,
                       fric_TEAK_9,
                       fric_TEAK_10)

fric_TEAK_all$site <- "TEAK"

fric_TEAK_sum <- fric_TEAK_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_TEAK_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_TEAK_1 <- fread("_TEAK_fdiv_veg_site_1.csv")
fdiv_TEAK_2 <- fread("_TEAK_fdiv_veg_site_2.csv")
fdiv_TEAK_3 <- fread("_TEAK_fdiv_veg_site_3.csv")
fdiv_TEAK_4 <- fread("_TEAK_fdiv_veg_site_4.csv")
fdiv_TEAK_5 <- fread("_TEAK_fdiv_veg_site_10.csv")
fdiv_TEAK_6 <- fread("_TEAK_fdiv_veg_site_12.csv")
fdiv_TEAK_7 <- fread("_TEAK_fdiv_veg_site_5.csv")
fdiv_TEAK_8 <- fread("_TEAK_fdiv_veg_site_7.csv")
fdiv_TEAK_9 <- fread("_TEAK_fdiv_veg_site_9.csv")
fdiv_TEAK_10 <- fread("_TEAK_fdiv_veg_site_0.csv")

fdiv_TEAK_1$plot <- '001'
fdiv_TEAK_2$plot <- '002'
fdiv_TEAK_3$plot <- '003'
fdiv_TEAK_4$plot <- '004'
fdiv_TEAK_5$plot <- '010'
fdiv_TEAK_6$plot <- '012'
fdiv_TEAK_7$plot <- '005'
fdiv_TEAK_8$plot <- '007'
fdiv_TEAK_9$plot <- '009'
fdiv_TEAK_10$plot <- '000'

fdiv_TEAK_all <- rbind(fdiv_TEAK_1,
                       fdiv_TEAK_2,
                       fdiv_TEAK_3,
                       fdiv_TEAK_4,
                       fdiv_TEAK_5,
                       fdiv_TEAK_6,
                       fdiv_TEAK_7,
                       fdiv_TEAK_8,
                       fdiv_TEAK_9,
                       fdiv_TEAK_10)

fdiv_TEAK_all$site <- "TEAK"

fdiv_TEAK_sum <- fdiv_TEAK_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_TEAK_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd,
                                 color = plot)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

# write.csv(fric_TEAK_all, "FRic_TEAK.csv")
# write.csv(fdiv_TEAK_all, "FDiv_TEAK.csv")

## TEAK (larger boundaries, v2)----
# FRIC
fric_TEAK_1 <- fread("_TEAK_fric_veg_002.csv")
fric_TEAK_2 <- fread("_TEAK_fric_veg_004.csv")
fric_TEAK_3 <- fread("_TEAK_fric_veg_006.csv")
fric_TEAK_4 <- fread("_TEAK_fric_veg_007.csv")
fric_TEAK_5 <- fread("_TEAK_fric_veg_011.csv")
fric_TEAK_6 <- fread("_TEAK_fric_veg_014.csv")

fric_TEAK_1$plot <- '002'
fric_TEAK_2$plot <- '004'
fric_TEAK_3$plot <- '006'
fric_TEAK_4$plot <- '007'
fric_TEAK_5$plot <- '011'
fric_TEAK_6$plot <- '014'

fric_TEAK_all <- rbind(fric_TEAK_1,
                       fric_TEAK_2,
                       fric_TEAK_3,
                       fric_TEAK_4,
                       fric_TEAK_5,
                       fric_TEAK_6)

fric_TEAK_all$site <- "TEAK"

fric_TEAK_sum <- fric_TEAK_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_TEAK_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_TEAK_1 <- fread("_TEAK_fdiv_veg_002.csv")
fdiv_TEAK_2 <- fread("_TEAK_fdiv_veg_004.csv")
fdiv_TEAK_3 <- fread("_TEAK_fdiv_veg_006.csv")
fdiv_TEAK_4 <- fread("_TEAK_fdiv_veg_007.csv")
fdiv_TEAK_5 <- fread("_TEAK_fdiv_veg_011.csv")
fdiv_TEAK_6 <- fread("_TEAK_fdiv_veg_014.csv")

fdiv_TEAK_1$plot <- '002'
fdiv_TEAK_2$plot <- '004'
fdiv_TEAK_3$plot <- '006'
fdiv_TEAK_4$plot <- '007'
fdiv_TEAK_5$plot <- '011'
fdiv_TEAK_6$plot <- '014'

fdiv_TEAK_all <- rbind(fdiv_TEAK_1,
                       fdiv_TEAK_2,
                       fdiv_TEAK_3,
                       fdiv_TEAK_4,
                       fdiv_TEAK_5,
                       fdiv_TEAK_6)

fdiv_TEAK_all$site <- "TEAK"

fdiv_TEAK_sum <- fdiv_TEAK_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_TEAK_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd,
                                 color = plot)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_TEAK_all, "FRic_TEAK.csv")
#write.csv(fdiv_TEAK_all, "FDiv_TEAK.csv")

## TOOL ----
# FRIC
fric_TOOL_1 <- fread("_TOOL_fric_veg_023.csv")
fric_TOOL_2 <- fread("_TOOL_fric_veg_022.csv")
fric_TOOL_3 <- fread("_TOOL_fric_veg_043.csv")
fric_TOOL_4 <- fread("_TOOL_fric_veg_028.csv")
fric_TOOL_5 <- fread("_TOOL_fric_veg_071.csv")
fric_TOOL_6 <- fread("_TOOL_fric_veg_024.csv")
fric_TOOL_7 <- fread("_TOOL_fric_veg_010.csv")
fric_TOOL_8 <- fread("_TOOL_fric_veg_014.csv")
fric_TOOL_9 <- fread("_TOOL_fric_veg_026.csv")
fric_TOOL_10 <- fread("_TOOL_fric_veg_018.csv")
fric_TOOL_11 <- fread("_TOOL_fric_veg_003.csv")
fric_TOOL_12 <- fread("_TOOL_fric_veg_020.csv")

fric_TOOL_1$plot <- '023'
fric_TOOL_2$plot <- '022'
fric_TOOL_3$plot <- '043'
fric_TOOL_4$plot <- '028'
fric_TOOL_5$plot <- '071'
fric_TOOL_6$plot <- '024'
fric_TOOL_7$plot <- '010'
fric_TOOL_8$plot <- '014'
fric_TOOL_9$plot <- '026'
fric_TOOL_10$plot <- '018'
fric_TOOL_11$plot <- '003'
fric_TOOL_12$plot <- '020'

fric_TOOL_all <- rbind(fric_TOOL_1,
                       fric_TOOL_2,
                       fric_TOOL_3,
                       fric_TOOL_4,
                       fric_TOOL_5,
                       fric_TOOL_6,
                       fric_TOOL_7,
                       fric_TOOL_8,
                       fric_TOOL_9,
                       fric_TOOL_10,
                       fric_TOOL_11,
                       fric_TOOL_12)

fric_TOOL_all$site <- "TOOL"

fric_TOOL_sum <- fric_TOOL_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_TOOL_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_TOOL_1 <- fread("_TOOL_fdiv_veg_023.csv")
fdiv_TOOL_2 <- fread("_TOOL_fdiv_veg_022.csv")
fdiv_TOOL_3 <- fread("_TOOL_fdiv_veg_043.csv")
fdiv_TOOL_4 <- fread("_TOOL_fdiv_veg_028.csv")
fdiv_TOOL_5 <- fread("_TOOL_fdiv_veg_071.csv")
fdiv_TOOL_6 <- fread("_TOOL_fdiv_veg_024.csv")
fdiv_TOOL_7 <- fread("_TOOL_fdiv_veg_010.csv")
fdiv_TOOL_8 <- fread("_TOOL_fdiv_veg_014.csv")
fdiv_TOOL_9 <- fread("_TOOL_fdiv_veg_026.csv")
fdiv_TOOL_10 <- fread("_TOOL_fdiv_veg_018.csv")
fdiv_TOOL_11 <- fread("_TOOL_fdiv_veg_003.csv")
fdiv_TOOL_12 <- fread("_TOOL_fdiv_veg_020.csv")

fdiv_TOOL_1$plot <- '023'
fdiv_TOOL_2$plot <- '022'
fdiv_TOOL_3$plot <- '043'
fdiv_TOOL_4$plot <- '028'
fdiv_TOOL_5$plot <- '071'
fdiv_TOOL_6$plot <- '024'
fdiv_TOOL_7$plot <- '010'
fdiv_TOOL_8$plot <- '014'
fdiv_TOOL_9$plot <- '026'
fdiv_TOOL_10$plot <- '018'
fdiv_TOOL_11$plot <- '003'
fdiv_TOOL_12$plot <- '020'

fdiv_TOOL_all <- rbind(fdiv_TOOL_1,
                       fdiv_TOOL_2,
                       fdiv_TOOL_3,
                       fdiv_TOOL_4,
                       fdiv_TOOL_5,
                       fdiv_TOOL_6,
                       fdiv_TOOL_7,
                       fdiv_TOOL_8,
                       fdiv_TOOL_9,
                       fdiv_TOOL_10,
                       fdiv_TOOL_11,
                       fdiv_TOOL_12)

fdiv_TOOL_all$site <- "TOOL"

fdiv_TOOL_sum <- fdiv_TOOL_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_TOOL_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_TOOL_all, "FRic_TOOL.csv")
#write.csv(fdiv_TOOL_all, "FDiv_TOOL.csv")

## UNDE ----
# FRIC
fric_unde_1 <- fread("_UNDE_fric_veg_001 (1).csv")
fric_unde_2 <- fread("_UNDE_fric_veg_003.csv")
fric_unde_3 <- fread("_UNDE_fric_veg_006.csv")
fric_unde_4 <- fread("_UNDE_fric_veg_011.csv")
fric_unde_5 <- fread("_UNDE_fric_veg_020.csv")
fric_unde_6 <- fread("_UNDE_fric_veg_023.csv")
fric_unde_7 <- fread("_UNDE_fric_veg_038.csv")

fric_unde_1$plot <- '001'
fric_unde_2$plot <- '003'
fric_unde_3$plot <- '006'
fric_unde_4$plot <- '011'
fric_unde_5$plot <- '020'
fric_unde_6$plot <- '023'
fric_unde_7$plot <- '038'

fric_unde_all <- rbind(fric_unde_1,
                       fric_unde_2,
                       fric_unde_3,
                       fric_unde_4,
                       fric_unde_5,
                       fric_unde_6,
                       fric_unde_7)

fric_unde_all$site <- "unde"

fric_unde_sum <- fric_unde_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_unde_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_unde_1 <- fread("_UNDE_fdiv_veg_001.csv")
fdiv_unde_2 <- fread("_UNDE_fdiv_veg_003.csv")
fdiv_unde_3 <- fread("_UNDE_fdiv_veg_006.csv")
fdiv_unde_4 <- fread("_UNDE_fdiv_veg_011.csv")
fdiv_unde_5 <- fread("_UNDE_fdiv_veg_020.csv")
fdiv_unde_6 <- fread("_UNDE_fdiv_veg_023.csv")
fdiv_unde_7 <- fread("_UNDE_fdiv_veg_038.csv")

fdiv_unde_1$plot <- '001'
fdiv_unde_2$plot <- '003'
fdiv_unde_3$plot <- '006'
fdiv_unde_4$plot <- '011'
fdiv_unde_5$plot <- '020'
fdiv_unde_6$plot <- '023'
fdiv_unde_7$plot <- '038'

fdiv_unde_all <- rbind(fdiv_unde_1,
                       fdiv_unde_2,
                       fdiv_unde_3,
                       fdiv_unde_4,
                       fdiv_unde_5,
                       fdiv_unde_6,
                       fdiv_unde_7)

fdiv_unde_all$site <- "unde"

fdiv_unde_sum <- fdiv_unde_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_unde_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_unde_all, "FRic_UNDE.csv")
#write.csv(fdiv_unde_all, "FDiv_UNDE.csv")
## WOOD ----
# FRIC
fric_wood_1 <- fread("_WOOD_fric_veg_008.csv")
fric_wood_2 <- fread("_WOOD_fric_veg_012.csv")
fric_wood_3 <- fread("_WOOD_fric_veg_019.csv")
fric_wood_4 <- fread("_WOOD_fric_veg_023.csv")
fric_wood_5 <- fread("_WOOD_fric_veg_024.csv")

fric_wood_1$plot <- '008'
fric_wood_2$plot <- '012'
fric_wood_3$plot <- '019'
fric_wood_4$plot <- '023'
fric_wood_5$plot <- '024'

fric_wood_all <- rbind(fric_wood_1,
                       fric_wood_2,
                       fric_wood_3,
                       fric_wood_4,
                       fric_wood_5)

fric_wood_all$site <- "WOOD"

fric_wood_sum <- fric_wood_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_wood_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_wood_1 <- fread("_WOOD_fdiv_veg_008.csv")
fdiv_wood_2 <- fread("_WOOD_fdiv_veg_012.csv")
fdiv_wood_3 <- fread("_WOOD_fdiv_veg_019.csv")
fdiv_wood_4 <- fread("_WOOD_fdiv_veg_023.csv")
fdiv_wood_5 <- fread("_WOOD_fdiv_veg_024.csv")

fdiv_wood_1$plot <- '008'
fdiv_wood_2$plot <- '012'
fdiv_wood_3$plot <- '019'
fdiv_wood_4$plot <- '023'
fdiv_wood_5$plot <- '024'

fdiv_wood_all <- rbind(fdiv_wood_1,
                       fdiv_wood_2,
                       fdiv_wood_3,
                       fdiv_wood_4,
                       fdiv_wood_5)

fdiv_wood_all$site <- "WOOD"

fdiv_wood_sum <- fdiv_wood_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_wood_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  geom_pointrange() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fdiv_wood_all, "FDiv_WOOD.csv")
#write.csv(fric_wood_all, "FRic_WOOD.csv")

## WREF ----
# FRIC
fric_WREF_1 <- fread("_WREF_fric_veg_070.csv")
fric_WREF_2 <- fread("_WREF_fric_veg_029.csv")
fric_WREF_3 <- fread("_WREF_fric_veg_026.csv")
fric_WREF_4 <- fread("_WREF_fric_veg_023.csv")
fric_WREF_5 <- fread("_WREF_fric_veg_021.csv")
fric_WREF_6 <- fread("_WREF_fric_veg_020.csv")
fric_WREF_7 <- fread("_WREF_fric_veg_015.csv")
fric_WREF_8 <- fread("_WREF_fric_veg_012.csv")
fric_WREF_9 <- fread("_WREF_fric_veg_006.csv")
fric_WREF_10 <- fread("_WREF_fric_veg_002.csv")

fric_WREF_1$plot <- '070'
fric_WREF_2$plot <- '029'
fric_WREF_3$plot <- '026'
fric_WREF_4$plot <- '023'
fric_WREF_5$plot <- '021'
fric_WREF_6$plot <- '020'
fric_WREF_7$plot <- '015'
fric_WREF_8$plot <- '012'
fric_WREF_9$plot <- '006'
fric_WREF_10$plot <- '002'

fric_WREF_all <- rbind(fric_WREF_1,
                       fric_WREF_2,
                       fric_WREF_3,
                       fric_WREF_4,
                       fric_WREF_5,
                       fric_WREF_6,
                       fric_WREF_7,
                       fric_WREF_8,
                       fric_WREF_9,
                       fric_WREF_10)

fric_WREF_all$site <- "WREF"

fric_WREF_sum <- fric_WREF_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_WREF_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_WREF_1 <- fread("_WREF_fdiv_veg_070.csv")
fdiv_WREF_2 <- fread("_WREF_fdiv_veg_029.csv")
fdiv_WREF_3 <- fread("_WREF_fdiv_veg_026.csv")
fdiv_WREF_4 <- fread("_WREF_fdiv_veg_023.csv")
fdiv_WREF_5 <- fread("_WREF_fdiv_veg_021.csv")
fdiv_WREF_6 <- fread("_WREF_fdiv_veg_020.csv")
fdiv_WREF_7 <- fread("_WREF_fdiv_veg_015.csv")
fdiv_WREF_8 <- fread("_WREF_fdiv_veg_012.csv")
fdiv_WREF_9 <- fread("_WREF_fdiv_veg_006.csv")
fdiv_WREF_10 <- fread("_WREF_fdiv_veg_002.csv")

fdiv_WREF_1$plot <- '070'
fdiv_WREF_2$plot <- '029'
fdiv_WREF_3$plot <- '026'
fdiv_WREF_4$plot <- '023'
fdiv_WREF_5$plot <- '021'
fdiv_WREF_6$plot <- '020'
fdiv_WREF_7$plot <- '015'
fdiv_WREF_8$plot <- '012'
fdiv_WREF_9$plot <- '006'
fdiv_WREF_10$plot <- '002'

fdiv_WREF_all <- rbind(fdiv_WREF_1,
                       fdiv_WREF_2,
                       fdiv_WREF_3,
                       fdiv_WREF_4,
                       fdiv_WREF_5,
                       fdiv_WREF_6,
                       fdiv_WREF_7,
                       fdiv_WREF_8,
                       fdiv_WREF_9,
                       fdiv_WREF_10)

fdiv_WREF_all$site <- "WREF"

fdiv_WREF_sum <- fdiv_WREF_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_WREF_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_WREF_all, "FRic_WREF.csv")
#write.csv(fdiv_WREF_all, "FDiv_WREF.csv")

## YELL ----
# FRIC
fric_YELL_1 <- fread("_YELL_fric_veg_003.csv")
fric_YELL_2 <- fread("_YELL_fric_veg_005.csv")
fric_YELL_3 <- fread("_YELL_fric_veg_011.csv")
fric_YELL_4 <- fread("_YELL_fric_veg_012.csv")
fric_YELL_5 <- fread("_YELL_fric_veg_013.csv")
fric_YELL_6 <- fread("_YELL_fric_veg_014.csv")
fric_YELL_7 <- fread("_YELL_fric_veg_015.csv")
fric_YELL_8 <- fread("_YELL_fric_veg_019.csv")
fric_YELL_9 <- fread("_YELL_fric_veg_021.csv")
fric_YELL_10 <- fread("_YELL_fric_veg_022.csv")
fric_YELL_11 <- fread("_YELL_fric_veg_024.csv")
fric_YELL_12 <- fread("_YELL_fric_veg_026.csv")
fric_YELL_13 <- fread("_YELL_fric_veg_046.csv")

fric_YELL_1$plot <- '001'
fric_YELL_2$plot <- '003'
fric_YELL_3$plot <- '004'
fric_YELL_4$plot <- '005'
fric_YELL_5$plot <- '007'
fric_YELL_6$plot <- '014'
fric_YELL_7$plot <- '015'
fric_YELL_8$plot <- '019'
fric_YELL_9$plot <- '021'
fric_YELL_10$plot <- '022'
fric_YELL_11$plot <- '024'
fric_YELL_12$plot <- '026'
fric_YELL_13$plot <- '046'

fric_YELL_all <- rbind(fric_YELL_1,
                       fric_YELL_2,
                       fric_YELL_3,
                       fric_YELL_4,
                       fric_YELL_5,
                       fric_YELL_6,
                       fric_YELL_7,
                       fric_YELL_8,
                       fric_YELL_9,
                       fric_YELL_10,
                       fric_YELL_11,
                       fric_YELL_12,
                       fric_YELL_13)

fric_YELL_all$site <- "YELL"

fric_YELL_sum <- fric_YELL_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(Hull_Volume), sd = sd(Hull_Volume))

ggplot(data = fric_YELL_sum, aes(x = Window_Size^2, y = median)) +
  geom_point() +
  facet_wrap(~plot, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FRic +/- sd")

# FDIV
fdiv_YELL_1 <- fread("_YELL_fdiv_veg_003.csv")
fdiv_YELL_2 <- fread("_YELL_fdiv_veg_005.csv")
fdiv_YELL_3 <- fread("_YELL_fdiv_veg_011.csv")
fdiv_YELL_4 <- fread("_YELL_fdiv_veg_012.csv")
fdiv_YELL_5 <- fread("_YELL_fdiv_veg_013.csv")
fdiv_YELL_6 <- fread("_YELL_fdiv_veg_014.csv")
fdiv_YELL_7 <- fread("_YELL_fdiv_veg_015.csv")
fdiv_YELL_8 <- fread("_YELL_fdiv_veg_019.csv")
fdiv_YELL_9 <- fread("_YELL_fdiv_veg_021.csv")
fdiv_YELL_10 <- fread("_YELL_fdiv_veg_022.csv")
fdiv_YELL_11 <- fread("_YELL_fdiv_veg_024.csv")
fdiv_YELL_12 <- fread("_YELL_fdiv_veg_026.csv")
fdiv_YELL_13 <- fread("_YELL_fdiv_veg_046.csv")

fdiv_YELL_1$plot <- '001'
fdiv_YELL_2$plot <- '003'
fdiv_YELL_3$plot <- '004'
fdiv_YELL_4$plot <- '005'
fdiv_YELL_5$plot <- '007'
fdiv_YELL_6$plot <- '014'
fdiv_YELL_7$plot <- '015'
fdiv_YELL_8$plot <- '019'
fdiv_YELL_9$plot <- '021'
fdiv_YELL_10$plot <- '022'
fdiv_YELL_11$plot <- '024'
fdiv_YELL_12$plot <- '026'
fdiv_YELL_13$plot <- '046'

fdiv_YELL_all <- rbind(fdiv_YELL_1,
                       fdiv_YELL_2,
                       fdiv_YELL_3,
                       fdiv_YELL_4,
                       fdiv_YELL_5,
                       fdiv_YELL_6,
                       fdiv_YELL_7,
                       fdiv_YELL_8,
                       fdiv_YELL_9,
                       fdiv_YELL_10,
                       fdiv_YELL_11,
                       fdiv_YELL_12,
                       fdiv_YELL_13)

fdiv_YELL_all$site <- "YELL"

fdiv_YELL_sum <- fdiv_YELL_all %>%
  group_by(plot, Window_Size) %>%
  summarize(median = median(FDiv), sd = sd(FDiv))

ggplot(data = fdiv_YELL_sum, aes(x = Window_Size^2, y = median,
                                 ymin = median - sd,
                                 ymax = median + sd)) +
  facet_wrap(~plot, scales = "free") +
  geom_pointrange() +
  theme_bw(base_size = 12) +
  labs(x = "Extent (m2)", y = "Median FDiv +/- sd")

#write.csv(fric_YELL_all, "FRic_YELL.csv")
#write.csv(fdiv_YELL_all, "FDiv_YELL.csv")

## Visualizing sites -----
fric_bart <- fread("FRic_BART.csv")
fdiv_bart <- fread("FDiv_BART.csv")
fric_clbj <- fread("FRic_CLBJ.csv")
fdiv_clbj <- fread("FDiv_CLBJ.csv")
fric_heal <- fread("FRic_HEAL.csv")
fdiv_heal <- fread("FDiv_HEAL.csv")
fric_konz <- fread("FRic_KONZ.csv")
fdiv_konz <- fread("FDiv_KONZ.csv")
fric_niwo <- fread("FRic_NIWO.csv")
fdiv_niwo <- fread("FDiv_NIWO.csv")
fric_onaq <- fread("FRic_ONAQ.csv")
fdiv_onaq <- fread("FDiv_ONAQ.csv")
fric_osbs <- fread("FRic_OSBS.csv")
fdiv_osbs <- fread("FDiv_OSBS.csv")
fric_puum <- fread("FRic_PUUM.csv")
fdiv_puum <- fread("FDiv_PUUM.csv")
fric_serc <- fread("FRic_SERC.csv")
fdiv_serc <- fread("FDiv_SERC.csv")
fric_srer <- fread("FRic_SRER.csv")
fdiv_srer <- fread("FDiv_SRER.csv")
fric_teak <- fread("FRic_TEAK.csv")
fdiv_teak <- fread("FDiv_TEAK.csv")
fric_tall <- fread("FRic_TALL.csv")
fdiv_tall <- fread("FDiv_TALL.csv")
fric_tool <- fread("FRic_TOOL.csv")
fdiv_tool <- fread("FDiv_TOOL.csv")
fric_unde <- fread("FRic_UNDE.csv")
fdiv_unde <- fread("FDiv_UNDE.csv")
fric_wood <- fread("FRic_WOOD.csv")
fdiv_wood <- fread("FDiv_WOOD.csv")
fric_wref <- fread("FRic_WREF.csv")
fdiv_wref <- fread("FDiv_WREF.csv")
fric_yell <- fread("FRic_YELL.csv")
fdiv_yell <- fread("FDiv_YELL.csv")

# Bind into one frame
fric_results <- rbind(fric_heal,
                      fric_konz,
                      fric_niwo,
                      fric_onaq,
                      fric_puum,
                      fric_srer,
                      fric_serc,
                      fric_teak,
                      fric_tall,
                      fric_tool,
                      fric_unde,
                      fric_clbj,
                      fric_osbs,
                      fric_wood,
                      fric_wref,
                      fric_yell, fric_bart)

fdiv_results <- rbind(fdiv_clbj,
                      fdiv_heal,
                      fdiv_konz,
                      fdiv_niwo,
                      fdiv_onaq,
                      fdiv_osbs,
                      fdiv_puum,
                      fdiv_srer,
                      fdiv_serc,
                      fdiv_teak,
                      fdiv_tall,
                      fdiv_tool,
                      fdiv_unde,
                      fdiv_wood,
                      fdiv_wref,
                      fdiv_yell, fdiv_bart)

# Summarize
fric_plot <- fric_results %>%
  mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
  group_by(site, plot, Window_Size) %>%
  summarise(median = median(Hull_Volume, na.rm = TRUE),
            sd = sd(Hull_Volume, na.rm = TRUE))

summary(lm(data = fric_plot, log(median) ~ log(I(Window_Size^2))))

fdiv_plot <- fdiv_results %>%
  group_by(site, plot, Window_Size) %>%
  summarise(median = median(FDiv, na.rm = TRUE),
            sd = sd(FDiv, na.rm = TRUE))

# FRic and FDiv by Window Size across sites
ggplot(data = fric_plot, aes(x = Window_Size^2, y = median, color = site)) +
  geom_point(size = 3) +
  #facet_wrap(~site, scales = "free") +
  theme_bw(base_size = 12) +
  xlim(0, 5e+06) +
  labs(x = "Window Area (m2)", y = "Median FRic")

ggplot(data = fdiv_plot, aes(x = Window_Size^2, y = median, color = site)) +
  geom_point(size = 3) +
  #facet_wrap(~site, scales = "free") +
  theme_bw(base_size = 12) +
  xlim(0, 5e+06) +
  labs(x = "Window Area (m2)", y = "Median FDiv")

## Model fitting - power and log -----
# Loop through sites
sites <- unique(fric_results$site)

# At plot level
results <- data.table(site = as.numeric(),
                      plot = as.numeric(),
                      fit = as.numeric(),
                      a = as.numeric(),
                      b = as.numeric(),
                      r2 = as.numeric())
sites <- sites[c(1, 3:13)] # keep ONAQ out while re-running
for (i in 1:length(sites)){
  fric_site <- fric_results %>%
    mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
    subset(site == sites[i])
  plots <- unique(fric_site$plot)
  for (j in 1:length(plots)){
    fric_plot <- fric_site %>%
      subset(plot == plots[j])
    lm <- lm(Hull_Volume ~ log(Window_Size^2), data = fric_plot)
    results_plot <- data.table(site = sites[i],
                               plot = plots[j],
                               fit = "ln",
                               a = coef(lm)[1], # exponential of intercept
                               b = coef(lm)[2], # slope is the exponent
                               r2 = summary(lm)$r.squared)
    results <- rbind(results, results_plot)
  }
}

# At site level
results_site <- data.table(site = as.numeric(),
                      fit = as.numeric(),
                      a = as.numeric(),
                      b = as.numeric(),
                      r2 = as.numeric())
for (i in 1:length(sites)){
  fric_site <- fric_results %>%
    mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
    subset(site == sites[i])
  lm <- lm(Hull_Volume ~ log(Window_Size^2), data = fric_site)
  results_sub <- data.table(site = sites[i],
                               fit = "power",
                               a = exp(coef(lm)[1]),
                               b = coef(lm)[2],
                               r2 = summary(lm)$r.squared)
  results_site <- rbind(results_site, results_sub)
}
x = fric_tall$Window_Size^2
y = fric_tall$Hull_Volume
plot(x,y, col = alpha("black", 0.4))
curve(723.244939*x^0.4575935, 0, 5e06, add = TRUE, col = "blue")

sar_fric <- fric_results %>%
  mutate(x = Window_Size^2) %>%
  group_by(site, plot, x) %>%
  summarize(mean = mean(Hull_Volume, na.rm = TRUE),
            median = median(Hull_Volume, na.rm = TRUE)) %>%
  rename(y = median)
# Define the power function
power_function <- function(x, a, b) {
  a * x^b
}

# Fit the model
fit <- nls(y ~ power_function(x, a, b), data = data.frame(sar_fric), start = list(a = exp(5), b = 0.3))
fit_glm <- glm(log(y) ~ log(x), data = sar_fric)
summary(fit_glm)
# Fit a GAM with a power function
gam_model <- gam(y ~ s(x, bs = "cr", k = 3), data = sar_fric)
summary(gam_model)
# Extract coefficients
coefficients <- coef(fit)
deviance(fit)

# Print coefficients
print(coefficients)

# Assess the fit (e.g., plot the data and the fitted curve)
plot(sar_fric$x, sar_fric$y, col = "blue", pch = 16, xlab = "x", ylab = "y")
curve(power_function(x, a = coefficients["a"], b = coefficients["b"]), add = TRUE, col = "red")
legend("topleft", legend = "Fitted Curve", col = "red", lty = 1)


sum_tall <- fric_tall %>%
  group_by(Window_Size, plot) %>%
  summarise(mean = mean(Hull_Volume))
ggplot(data = fric_tall, aes(x = Window_Size^2, y = Hull_Volume)) +
  geom_point(data = sum_tall, aes(x = Window_Size^2, y = mean),
             alpha = 0.5, size = 3) +
  geom_function(fun = function(x) 723.244939*x^0.4575935, color = "black") +
  theme_bw(base_size = 12) +
  labs(x = "Area", y = "FRic")


# At plot level
results_power_window <- data.table(site = as.numeric(),
                      plot = as.numeric(),
                      fit = as.numeric(),
                      c = as.numeric(),
                      z = as.numeric(),
                      resid_dev = as.numeric())
for (i in 1:length(sites)){
  fric_site <- fric_results %>%
    mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
    mutate(Area = Window_Size^2) %>%
    subset(site == sites[i])
  plots <- unique(fric_site$plot)
  for (j in 1:length(plots)){
    fric_plot <- fric_site %>%
      subset(plot == plots[j])
    lm <- glm(log(Hull_Volume) ~ log(Area), data = fric_plot)
    results_plot <- data.table(site = sites[i],
                               plot = plots[j],
                               fit = "power",
                               c = exp(coef(lm)[1]),
                               z = coef(lm)[2],
                               resid_dev = 1 -(summary(lm)$deviance/summary(lm)$null.deviance))
    results_power_window <- rbind(results_power_window, results_plot)
  }
}

# At site level
results_power_site_window <- data.table(site = as.numeric(),
                           fit = as.numeric(),
                           c = as.numeric(),
                           z = as.numeric(),
                           resid_dev = as.numeric())
for (i in 1:length(sites)){
  fric_site <- fric_results %>%
    mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
    mutate(Area = Window_Size^2) %>%
    subset(site == sites[i])
  lm <- glm(log(Hull_Volume) ~ log(Area), family = gaussian, data = fric_site)
  results_sub <- data.table(site = sites[i],
                            fit = "power",
                            c = exp(coef(lm)[1]),
                            z = coef(lm)[2],
                            resid_dev = summary(lm)$deviance)
  results_power_site <- rbind(results_power_site, results_sub)
}

lm <- lm(Hull_Volume ~ Window_Size, data = fric_site)
summary(lm)

ggplot(data = results_power, aes(x = exp(c), y = site, color = resid_dev)) +
  geom_point(alpha = 0.5, size = 2, aes(color = resid_dev)) +
  #geom_point(data = results_power_site, aes(x = c, y = site, color = resid_dev),
  #           size = 3) +
  #facet_wrap(~site, scales = "free") +
  scale_color_gradient(low = "red", high = "green") +
  theme_bw(base_size = 12) +
  xlim(0,2500) +
  labs(x = "Scaling parameter, c", y = "Site", title = "Area",
       color = "Residual \nDeviance")

# Investigate outliers
ggplot(data = subset(fric, site == "ONAQ"), aes(x = Window_Size^2, y = median)) +
  #facet_wrap(~plot, scales = "free") +
  geom_point(alpha = 0.3) +
  geom_function(fun = function(x) 832.021267*x^0.3683893, color = "black") +
  theme_bw(base_size = 12) +
  labs(x = "Area", y = "FRic")
ggplot(data = subset(fric, site == "ONAQ"), aes(x = Window_Size^2, y = median)) +
  #facet_wrap(~plot, scales = "free") +
  geom_point(alpha = 0.3) +
  geom_function(fun = function(x) 115*x^0.5, color = "black") +
  theme_bw(base_size = 12) +
  labs(x = "Area", y = "FRic")

## Model Fitting with standardized z ----------
get_c <- function(df, z){
  df$Area <- df$Window_Size^2
  c_table <- data.frame(c = seq(0,2000,1), mse = NA)
  for (i in 1:nrow(c_table)){
    df_predict <- df %>%
      mutate(pred = c_table$c[i]*Area^z)
    c_table$mse[i] <- mean((df_predict$pred - df_predict$Hull_Volume)^2)
  }
  return(c_table$c[which(c_table$mse == min(c_table$mse))])
}

get_c_lin <- function(df, z){
  df$Area <- df$Window_Size^2
  c_table <- data.frame(c = seq(0,3000,3), mse = NA)
  for (i in 1:nrow(c_table)){
    df_predict <- df %>%
      mutate(pred = log(c_table$c[i]) + z*log(Area))
    c_table$mse[i] <- mean((df_predict$pred - log(df_predict$Hull_Volume))^2)
  }
  return(c_table$c[which(c_table$mse == min(c_table$mse))])
}

get_mse_lin <- function(df, z){
  df$Area <- df$Window_Size^2
  c_table <- data.frame(c = seq(0,3000,3), mse = NA)
  for (i in 1:nrow(c_table)){
    df_predict <- df %>%
      mutate(pred = log(c_table$c[i]) + z*log(Area))
    c_table$mse[i] <- mean((df_predict$pred - log(df_predict$Hull_Volume))^2)
  }
  return(c_table$mse[which(c_table$mse == min(c_table$mse))])
}

get_deviance_lin <- function(df, z){
  df$Area <- df$Window_Size^2
  c_table <- data.frame(c = seq(0,2000,1), mse = NA)
  for (i in 1:nrow(c_table)){
    df_predict <- df %>%
      mutate(pred = log(c_table$c[i]) + z*log(Area))
    c_table$mse[i] <- mean((df_predict$pred - log(df_predict$Hull_Volume))^2)
  }
  return(c_table$mse[which(c_table$mse == min(c_table$mse))])
}

z = 0.4873393 # averaged across plot power functions
sites <- unique(fric_results$site) # if not already done elsewhere

c_results <- data.table(site = as.numeric(),
                        plot = as.numeric(),
                        z = as.numeric(),
                        c = as.numeric(),
                        mse = as.numeric()
)

for (i in 1:length(sites)){
    fric_site <- fric_results %>%
      mutate(Hull_Volume = as.numeric(Hull_Volume)) %>%
      subset(site == sites[i])
    plots <- unique(fric_site$plot)
    for (j in 1:length(plots)){
      fric_plot <- fric_site %>%
        subset(plot == plots[j])
      opt_c <- get_c_lin(fric_plot, z)
      opt_mse <- get_mse_lin(fric_plot, z)
      c_plot <- data.table(site = sites[i],
                           plot = plots[j],
                           z = z,
                           c = opt_c,
                           mse = opt_mse)
      c_results <- rbind(c_results, c_plot)
    }
}

# Played around with different z values --> rank order of sites by c value does NOT change
# Decided to stick with average z from power functions
ggplot(data = c_results, aes(x = c, y = site, color = as.factor(z))) +
  geom_point() +
  scale_color_manual(values = c("red3", "darkgreen", "lightblue",
                                "orange", "purple")) +
  theme_bw(base_size = 12) +
  labs(x = "Scaling parameter, c", y = "Site")

rank_0.5 <- c_results %>%
  mutate(site = as.factor(site)) %>%
  filter(z == 0.5) %>%
  group_by(site) %>%
  summarise(mean = mean(c)) %>%
  arrange(mean)

rank_0.48 <- c_results %>%
  mutate(site = as.factor(site)) %>%
  filter(z == 0.48) %>%
  group_by(site) %>%
  summarise(mean = mean(c)) %>%
  arrange(mean)

rank_0.52 <- c_results %>%
  mutate(site = as.factor(site)) %>%
  filter(z == 0.52) %>%
  group_by(site) %>%
  summarise(mean = mean(c)) %>%
  arrange(mean)

merged_ranks <- merge(rank, rank_0.48, by = "site")
merged_ranks <- merge(merged_ranks, rank_0.52, by = "site")

ggplot(data = merged_ranks, aes(x = mean.x, y = mean.y, color = mean)) +
  geom_point()

## Compare to NEON environmental covariates ------
# Read in files
dtm_bart <- fread("BART_DTM_Summary.csv")
chm_bart <- fread("BART_CHM_Summary.csv")
dtm_clbj <- fread("CLBJ_DTM_Summary.csv")
chm_clbj <- fread("CLBJ_CHM_Summary.csv")
dtm_heal <- fread("HEAL_DTM_Summary.csv")
chm_heal <- fread("HEAL_CHM_Summary.csv")
#dtm_konz <- fread("KONZ_DTM_Summary.csv")
chm_konz <- fread("KONZ_CHM_Summary.csv")
dtm_niwo <- fread("NIWO_DTM_Summary.csv")
chm_niwo <- fread("NIWO_CHM_Summary.csv")
#dtm_onaq <- fread("ONAQ_DTM_Summary.csv")
#chm_onaq <- fread("ONAQ_CHM_Summary.csv")
dtm_osbs <- fread("OSBS_DTM_Summary.csv")
chm_osbs <- fread("OSBS_CHM_Summary.csv")
dtm_puum <- fread("PUUM_DTM_Summary.csv")
chm_puum <- fread("PUUM_CHM_Summary.csv")
dtm_serc <- fread("SERC_DTM_Summary.csv")
chm_serc <- fread("SERC_CHM_Summary.csv")
dtm_srer <- fread("SRER_DTM_Summary.csv")
chm_srer <- fread("SRER_CHM_Summary.csv")
dtm_teak <- fread("TEAK_DTM_Summary (1).csv")
chm_teak <- fread("TEAK_CHM_Summary (1).csv")
dtm_tall <- fread("TALL_DTM_Summary.csv")
chm_tall <- fread("TALL_CHM_Summary.csv")
dtm_tool <- fread("TOOL_DTM_Summary.csv")
chm_tool <- fread("TOOL_CHM_Summary.csv")
dtm_unde <- fread("UNDE_DTM_Summary.csv")
chm_unde <- fread("UNDE_CHM_Summary.csv")
dtm_wood <- fread("WOOD_DTM_Summary.csv")
chm_wood <- fread("WOOD_CHM_Summary.csv")
dtm_wref <- fread("WREF_DTM_Summary.csv")
chm_wref <- fread("WREF_CHM_Summary.csv")
dtm_yell <- fread("YELL_DTM_Summary.csv")
chm_yell <- fread("YELL_CHM_Summary.csv")

# Bind into one frame
dtm_results <- rbind(dtm_heal,
                      #dtm_konz,
                      dtm_niwo,
                      #dtm_onaq,
                      dtm_puum,
                      dtm_srer,
                      dtm_serc,
                      dtm_teak,
                      dtm_tall,
                      dtm_tool,
                      dtm_unde,
                      dtm_clbj,
                      dtm_osbs,
                      dtm_wood,
                      dtm_wref,
                      dtm_yell)

chm_results <- rbind(chm_clbj,
                      chm_heal,
                      chm_konz,
                      chm_niwo,
                      #chm_onaq,
                      chm_osbs,
                      chm_puum,
                      chm_srer,
                      chm_serc,
                      chm_teak,
                      chm_tall,
                      chm_tool,
                      chm_unde,
                      chm_wood,
                      chm_wref,
                      chm_yell)

env_results <- rbind(chm_results, dtm_results)
env_results$Std <- as.numeric(env_results$Std)


env_results$plot <- str_pad(env_results$Plot, 3, pad =
                              "0")

env_results <- env_results %>%
  rename(site = Site)

env_results <- merge(env_results, coords_all, by = c("site", "plot"))

# Visualize
# CHM
ggplot(data = subset(env_results, Env == "CHM" & Site != "KONZ"),
       aes(x = Std, y = Site)) +
  geom_point(alpha = 0.3, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Std Canopy Height (m)")
# DTM
ggplot(data = subset(env_results, Env == "DTM" & Site != "KONZ"),
       aes(x = Std, y = Site)) +
  geom_point(alpha = 0.3, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Std Elevation (m)")

# To wide for looking at both covariates at once
env_results <- env_results %>%
  mutate(Site_plot = paste(Site, "_", Plot))
env_wide <- env_results %>%
  pivot_wider(
    names_from = Env,
    values_from = c(Mean, Median, Max, Min, Var, Std)
  )

# Visualize in 2D space
ggplot(data = subset(env_wide, Site != "KONZ"),
       aes(x = Mean_CHM, y = Mean_DTM, color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Canopy Height (m)", y = "Mean Elevation (m)")
ggplot(data = subset(env_wide,site != "KONZ"),
       aes(x = ppt_annual, y = temp_annual, color = Mean_CHM)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Annual Precip (mm)", y = "Mean Annual Temperature (C)")
ggplot(data = subset(env_wide, Site != "KONZ"),
       aes(x = Std_CHM, y = Std_DTM, color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Std Canopy Height (m)", y = "Std Elevation (m)")

# Save env results
#write.csv(env_wide, "NEON_env.csv")

env_wide <- fread("NEON_env.csv")

# Compare to power functions
results_power_window <- results_power_window %>%
  #rename(Site = site, Plot = plot) %>%
  mutate(Plot = as.character(Plot))
env_power <- merge(results_power_window, env_wide, by = c("Site",
                                                         "Plot"))
# CHM and power function fits
ggplot(data = subset(env_power, Site != "KONZ"),
       aes(x = Mean_CHM, y = log(resid_dev), color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Canopy Height (m)", y = "log(Deviance)")

# DTM and power function fits
ggplot(data = subset(env_power, Site != "KONZ"),
       aes(x = Mean_DTM, y = z, color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Elevation (m)", y = "Scaling exponent, z")

# Compare to standardized z power functions
c_results <- c_results %>%
  rename(Site = site, Plot = plot)

env_standard_power <- merge(env_wide, c_results, by = c("Site",
                                                        "Plot"))
# DTM and standard power function fits
ggplot(data = subset(env_standard_power, Site != "KONZ"),
       aes(x = Mean_DTM, y = mse, color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Elevation (m)", y = "MSE")

# CHM and standard power function fits
ggplot(data = subset(env_standard_power, Site != "KONZ"),
       aes(x = Mean_CHM, y = c, color = Site)) +
  geom_point(alpha = 0.6, size = 3) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Canopy Height (m)", y = "Scaling parameter, c")

## Compare to PRISM climate data -----
climate <- fread("NEON_PRISM_data.csv")
c_clim <- merge(sr_c_sum, climate, by = c("site", "plot"))
c_clim$index <- c_clim$temp_annual*c_clim$ppt_annual

temp <- ggplot(data = c_clim, aes(x = temp_annual, y = c_stand, color = veg_cat_new2)) +
  geom_point(size = 3, alpha = 0.5) +
  #geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 12) +
  scale_color_viridis_d() +
  ylim(0, 1500) +
  labs(x = "Mean Annual Temperature (C)", y = expression("Standardized scaling coefficient, c"[0.48]),
       color = "Vegetation Type")
ppt <- ggplot(data = c_clim, aes(x = ppt_annual, y = c_stand,  color = veg_cat_new2)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  ylim(0, 1500) +
  scale_color_viridis_d() +
  theme_bw(base_size = 12) +
  labs(x = "Mean Annual Precipitation (mm)", y = expression("Standardized scaling coefficient, c"[0.48]),
       color = "Vegetation Type")
index <- ggplot(data = c_clim, aes(x = index, y = c_stand,  color = veg_cat_new2)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  ylim(0, 1500) +
  scale_color_viridis_d() +
  theme_bw(base_size = 12) +
  labs(x = "Climate Index (mm*C)", y = expression("Standardized scaling coefficient, c"[0.48]),
       color = "Vegetation Type")
ggpubr::ggarrange(ppt, temp + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),
                  index + theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank()), ncol = 3, common.legend = TRUE)
summary(lm(c_stand ~ index, data = c_clim))

# Look at site level
c_clim_sum <- c_clim %>%
  group_by(site) %>%
  summarize(mean_c = mean(c, na.rm = TRUE),
            mean_ppt = mean(ppt_annual, na.rm = TRUE),
            mean_temp = mean(temp_annual, na.rm = TRUE))

temp_site <- ggplot(data = c_clim_sum, aes(x = mean_temp, y = mean_c, color = site)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 12) +
  ylim(0, 1500) +
  labs(x = "Mean Annual Temperature (C)", y = "Scaling parameter, c (z = 0.48)",
       color = "Site")
ppt_site <- ggplot(data = c_clim_sum, aes(x = mean_ppt, y = mean_c, color = site)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  ylim(0, 1500) +
  theme_bw(base_size = 12) +
  labs(x = "Mean Annual Precipitation (in)", y = "Scaling parameter, c (z = 0.48)",
       color = "Site")
ggpubr::ggarrange(ppt_site, temp_site, common.legend = TRUE)

c_clim_sum <- c_clim %>%
  group_by(site) %>%
  summarize(mean_c = mean(c, na.rm = TRUE),
            mean_ppt = mean(ppt_annual, na.rm = TRUE),
            mean_temp = mean(temp_annual, na.rm = TRUE))

## Compare with species richness ------
fric$Area <- fric$Window_Size^2
sr <- fread("NEON_Species_Richness_in_situ.csv")
sr_veg <- sr %>%
  select(site, plot, veg_cat) %>%
  distinct()
sr_fric <- merge(sr, fric_plot, by = c("site", "plot"))
sr_fdiv <- merge(sr, fdiv_plot, by = c("site", "plot"))
sr_fric$Area <- sr_fric$Window_Size^2
sr_fdiv$Area <- sr_fdiv$Window_Size^2
lm <- glm(log(median) ~ log(Area), data = sr_fric)
summary(lm)
predicted.values <- predict(lm, type = "response")
preds <- exp(predicted.values)

veg_count <- sr_fric %>%
  group_by(veg_cat_new) %>%
  distinct(site, plot, .keep_all = TRUE) %>%
  summarise(n = n())

sr_fric <- merge(sar_fric, sr, by = c("site", "plot"))
sr_fric$Area <- sr_fric$Window_Size^2

sar_fric <- fric_results %>%
  group_by(site, plot, Window_Size) %>%
  summarize(mean = mean(Hull_Volume, na.rm = TRUE),
            median = median(Hull_Volume, na.rm = TRUE))
sr_fric <- merge(sar_fric, sr, by = c("site", "plot"))
sr_fric$Area <- sr_fric$Window_Size^2

scale_parameters_sr <- scale_parameters_sr %>%
  #mutate(Area = Window_Size^2) %>%
  #mutate(site_plot = paste0(site, "_", plot)) %>%
  mutate(veg_cat_new = ifelse(veg_cat == "sedgeHerbaceous" |veg_cat == "grasslandHerbaceous" | veg_cat == "emergentHerbaceousWetlands" | veg_cat == "woodyWetlands",
                              "Herbaceous", veg_cat)) %>%
  mutate(veg_cat_new = ifelse(veg_cat == "dwarfScrub" | veg_cat == "shrubScrub", "Scrub", veg_cat_new)) %>%
  mutate(veg_cat_new = ifelse(veg_cat == "deciduousForest" | veg_cat == "evergreenForest" | veg_cat == "mixedForest", "Forest",veg_cat_new)) %>%
  mutate(veg_cat_new = ifelse(site == "SERC" | site_plot == "WOOD_12", "Cropland",veg_cat_new))

fdiv_plot <- ggplot(data = sr_fdiv, aes(x = Area, y = median, color = veg_cat_new)) +
  geom_point(alpha = 0.2, size = 3) +
  #facet_wrap(~veg_cat_new, scales = "fixed") +
  stat_smooth(method = 'nls', formula = 'y~a*x^b',
             method.args = list(start=c(a=25, b=0.5)), se=FALSE, color = "black") +
  stat_smooth(method = 'nls', formula = 'y~a*x^b',
              method.args = list(start=c(a=25, b=0.5)),se=FALSE, alpha = 0.2,
              linetype = 2) +
  theme_bw(base_size = 12) +
  coord_cartesian(clip = "off") +
  labs(x = expression('Area (m'^'2' *')'), y = "Functional Divergence",
       color = "Vegetation Type")

fric_plot <- ggplot(data = sr_fric, aes(x = x, y = y, color = veg_cat_new)) +
  geom_point(alpha = 0.2, size = 3) +
  #facet_wrap(~veg_cat_new, scales = "fixed") +
  stat_smooth(method = 'nls', formula = 'y~a*x^b',
              method.args = list(start=c(a=25, b=0.5)), se=FALSE, color = "black") +
  stat_smooth(method = 'nls', formula = 'y~a*x^b',
              method.args = list(start=c(a=25, b=0.5)),se=FALSE, alpha = 0.2,
              linetype = 2) +
  theme_bw(base_size = 12) +
  coord_cartesian(clip = "off") +
  labs(x = expression('Area (m'^'2' *')'), y = "Functional Richness",
       color = "Vegetation Type")

ggpubr::ggarrange(fric_plot, fdiv_plot, common.legend = TRUE)

summary(glm(data = sr_fdiv, log(median) ~ log(Area)))

sr_unique <- sr_unique %>%
  mutate(veg_cat_new = ifelse(veg_cat == "sedgeHerbaceous" |veg_cat == "grasslandHerbaceous" | veg_cat == "emergentHerbaceousWetlands",
                              "Herbaceous (n = 15)", veg_cat)) %>%
  mutate(veg_cat_new = ifelse(veg_cat == "dwarfScrub" | veg_cat == "shrubScrub", "Scrub (n = 30)", veg_cat_new)) %>%
  mutate(veg_cat_new = ifelse(veg_cat == "deciduousForest" | veg_cat == "evergreenForest" | veg_cat == "mixedForest", "Forest (n = 29)",veg_cat_new))

power_veg <- merge(sr_unique, results_power_window, by = c("site", "plot"))
stand_veg <- merge(sr_unique, c_results, by = c("site", "plot"))

power_veg_sum <- power_veg %>%
  group_by(veg_cat_new) %>%
  filter(c < 2500) %>%
  summarise(mean_z = mean(z),
         mean_c = mean(c))

stand_veg_sum <- stand_veg %>%
  group_by(veg_cat_new) %>%
  filter(c < 2500) %>%
  summarise(mean_z = mean(z),
            mean_c = mean(c))

power_z <- ggplot(data = power_veg, aes(x = z, y = site, color = veg_cat_new)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_vline(xintercept = 0.5422, color = "red3", linetype = 2) +
  geom_vline(xintercept = 0.5944, color = "green3", linetype = 2) +
  geom_vline(xintercept = 0.4242, color = "blue3", linetype = 2) +
  theme_bw(base_size = 12) +
  labs(x = "Scaling exponent, z", y = "Site", color = "Vegetation Type")
power_c <-ggplot(data = power_veg, aes(x = c, y = site, color = veg_cat_new)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_vline(xintercept = 505, color = "red3", linetype = 2) +
  geom_vline(xintercept = 185, color = "green3", linetype = 2) +
  geom_vline(xintercept = 451, color = "blue3", linetype = 2) +
  theme_bw(base_size = 12) +
  xlim(0,2500) +
  labs(x = "Scaling parameter, c", y = "Site",color = "Vegetation Type")
power_cn <- ggplot(data = stand_veg, aes(x = c, y = site, color = veg_cat_new)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_vline(xintercept = 453, color = "red3", linetype = 2) +
  geom_vline(xintercept = 376, color = "green3", linetype = 2) +
  geom_vline(xintercept = 327, color = "blue3", linetype = 2) +
  theme_bw(base_size = 12) +
  xlim(0,1500) +
  labs(x = "Standardized scaling parameter, c[n]", y = "Site",
       color = "Vegetation Type")
ggpubr::ggarrange(power_z, power_c + theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank()),
                  power_cn + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank()),
                  nrow = 1,
                  labels = "auto",
                  align = "v",
                  common.legend = TRUE)

sr_z_alpha <- ggplot(data = power_veg, aes(x = z, y = ave_alpha, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  labs(x = "Scaling exponent, z", y = "Alpha diversity", color = "Vegetation Type")

sr_z_beta <- ggplot(data = power_veg, aes(x = z, y = beta, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  labs(x = "Scaling exponent, z", y = "Beta diversity", color = "Vegetation Type")

sr_z_gamma <- ggplot(data = power_veg, aes(x = z, y = gamma, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  labs(x = "Scaling exponent, z", y = "Gamma diversity", color = "Vegetation Type")

ggpubr::ggarrange(sr_z_alpha  + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank())
                  , sr_z_beta + theme(axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.title.x = element_blank()),
                  sr_z_gamma,
                  ncol = 1,
                  labels = "auto",
                  common.legend = TRUE)

sr_c_alpha <- ggplot(data = stand_veg, aes(x = c, y = ave_alpha, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  xlim(0,1500) +
  labs(x = "Standardized scaling parameter, c", y = "Alpha diversity", color = "Vegetation Type")

sr_c_beta <- ggplot(data = stand_veg, aes(x = c, y = beta, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  xlim(0,1500) +
  labs(x = "Standardized scaling parameter, c", y = "Beta diversity", color = "Vegetation Type")

sr_c_gamma <- ggplot(data = stand_veg, aes(x = c, y = gamma, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw(base_size = 12) +
  xlim(0,1500) +
  labs(x = "Standardized scaling parameter, c", y = "Gamma diversity", color = "Vegetation Type")

ggpubr::ggarrange(sr_c_alpha  + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank())
                  , sr_c_beta + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title.x = element_blank()),
                  sr_c_gamma,
                  ncol = 1, align = "v",
                  labels = "auto",
                  common.legend = TRUE)


ggplot(data = subset(sr_fric_merge, veg_cat_new != "Wetland"),
       aes(x = z, y = log(c), color = veg_cat_new)) +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~veg_cat_new, scales = "fixed", nrow = 2) +
  theme_bw(base_size = 12) +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("red3", "green4")) +
  labs(x = "Scaling exponent, z", y = "Scaling parameter, log(c)")

z <- ggplot(data = sr_fric_merge,
       aes(x = z, y = site, color = veg_cat_new)) +
  geom_point(alpha = 0.8, size = 2) +
  #facet_wrap(~site, scales = "fixed") +
  theme_bw(base_size = 12) +
  theme(legend.position = "NA") +
  labs(x = "Scaling exponent, z", y = "Site")

c <- ggplot(data = sr_fric_merge,
            aes(x = log(c), y = site, color = veg_cat_new)) +
  geom_point(alpha = 0.8, size = 2) +
  #facet_wrap(~site, scales = "fixed") +
  theme_bw(base_size = 12) +
  labs(x = "Scaling exponent, log(c)", y = "Site",
       color = "Vegetation \nType")

ggpubr::ggarrange(z, c, common.legend = TRUE)
cowplot::plot_grid(z,
                   c + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),
                   nrow = 1,
                   labels = "auto", align = "v")

env <- fread("NEON_env.csv")
env <- env[,-1]
env <- env %>%
  mutate(CV_CHM = Std_CHM/Mean_CHM,
         CV_DTM = Std_DTM/Mean_DTM) %>%
 # rename(site = Site, plot = Plot) %>%
  mutate(plot = as.numeric(plot))
c_env <- merge(c_results, env, by = c("site", "plot"))
c_env_sr <- merge(c_env, sr, by = c("site", "plot"))

ggplot(data = c_env, aes(x = CV_DTM, y = c)) +
  geom_point(alpha = 0.2, size = 3) +
  #facet_wrap(~veg_cat_new, scales = "fixed") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 12) +
  xlim(0,0.25) +
  labs(x = "CV of Elevation", y = "Scaling parameter, c (z = 0.48)")

ggplot(data = c_env_sr, aes(x = CV_DTM, y = c)) +
  geom_point(alpha = 0.2, size = 3) +
  #facet_wrap(~veg_cat_new, scales = "fixed") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw(base_size = 12) +
  xlim(0,0.25) +
  labs(x = "Scaling parameter, c (z = 0.48)",
       y = "MSE")

summary(lm(data = subset(sr_fric, veg_cat_new == "Forest"), log(median) ~ log(Area)))

stand_env <- merge(stand_veg, env, by = c("site", "plot"))
chm <- ggplot(data = stand_env, aes(x = log(CV_CHM), y = c, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw(base_size = 12) +
  labs(y = "Standardized scaling parameter, c",
       x = "log(Canopy Height CV)", color = "Vegetation Type")
dtm <- ggplot(data = stand_env, aes(x = log(CV_DTM), y = c, color = veg_cat_new)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw(base_size = 12) +
  labs(y = "Standardized scaling parameter, c",
       x = "log(Elevation CV)", color = "Vegetation Type")

ggpubr::ggarrange(dtm,
                   chm + theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank()),
                   nrow = 1,
                   labels = "auto", align = "v",
                   common.legend = TRUE)

## Compare NLS and LM model fits ------
# From Xiao et al. 2011 (Ecology)
## Loading required packages. May need to be installed from CRAN.
library(nlrwr)
library(boot)
## Function to carry out the analytical process described in General Guidelines
## Input: x - vector, explanatory variable
##        y - vector, response variable
##        CI_boot - whether confidence intervals should be computed using bootstrapping
##                  for model averaging. Default is TRUE.
##        diagno - residual plots as visual diagnostics to check the assumptions of
##                 normality and homoscedasticity for LR and NLR. Default is TRUE.
##        output_plot - whether a scatter plot with fitted relationship is desired.
##                      Default is FALSE.
## Output: method - method used in analysis
##         a, b - estimated parameters
##         a_confint, b_confint - 95% confidence intervals for a & b
##                                (optional if method is "Model Averaging")
x = sar_fric$x
y = sar_fric$y
power_analysis(x, y, CI_boot = TRUE, diagno = TRUE, output_plot = TRUE)
power_analysis = function(x, y, CI_boot = TRUE, diagno = TRUE, output_plot = FALSE){
  ## Step 1: Likelihood analysis
  model_lr = lm(log(y) ~ log(x))
  a_lr = exp(coef(summary(model_lr))[1, 1])
  b_lr = coef(summary(model_lr))[2, 1]
  sd_lr = sd(log(y) - (log(a_lr) + b_lr * log(x)))

  model_nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a_lr, a2 = b_lr),
                  control = nls.control(maxiter = 2000, warnOnly = TRUE))
  a_nlr = coef(summary(model_nlr))[1, 1]
  b_nlr = coef(summary(model_nlr))[2, 1]
  sd_nlr = sd(y - a_nlr * x ^ b_nlr)

  l_logn = sum(log(dlnorm(y, log(a_lr * x ^ b_lr), sd_lr)))
  l_norm = sum(log(dnorm(y, a_nlr * x ^ b_nlr, sd_nlr)))

  n = length(x)
  k = 3
  AICc_logn = 2 * k - 2 * l_logn + 2 * k * (k + 1) / (n - k - 1)
  AICc_norm = 2 * k - 2 * l_norm + 2 * k * (k + 1) / (n - k - 1)
  delta_AICc = AICc_norm - AICc_logn
  writeLines(paste("AICc_logn: ", AICc_logn, "\nAICc_norm: ", AICc_norm))
  w_logn = exp(-(AICc_logn - min(AICc_logn, AICc_norm)) / 2)
  w_norm = exp(-(AICc_norm - min(AICc_logn, AICc_norm)) / 2)
  weight_logn = w_logn / (w_logn + w_norm)
  weight_norm = w_norm / (w_logn + w_norm)

  ## Step 2a: Analysis with NLR
  if (delta_AICc < - 2){
    writeLines("The assumption of additive normal error is better supported. \nProceed with NLR.")
    method = "NLR"
    a = a_nlr
    b = b_nlr
    a_confint = confint2(model_nlr)[1, ]
    b_confint = confint2(model_nlr)[2, ]
  }
  ## Step 2b: Analysis with LR
  else if (delta_AICc > 2){
    writeLines("The assumption of multiplicative log-normal error is better supported. \nProceed with LR.")
    method = "LR"
    a = a_lr
    b = b_lr
    a_confint = confint(model_lr)[1, ]
    b_confint = confint(model_lr)[2, ]
  }
  ## Step 2c: Analysis with model averaging
  else {
    writeLines("The two error distributions have similar support. \nProceed with model averaging.")
    method = "Model Averaging"
    a = a_lr * weight_logn + a_nlr * weight_norm
    b = b_lr * weight_logn + b_nlr * weight_norm
    if (!CI_boot){
      a_confint = NA
      b_confint = NA
    }
    else {
      boot.est=function(dat, indices) {
        dat.sub=dat[indices, ]
        names(dat.sub) = c("x", "y")
        model.lr = lm(log(y) ~ log(x), dat = dat.sub)
        a.lr = exp(coef(summary(model.lr))[1, 1])
        b.lr = coef(summary(model.lr))[2, 1]
        sd.lr = sd(log(dat.sub$y) - (log(a.lr) + b.lr * log(dat.sub$x)))
        a.lr.CI = confint(model.lr)[1, ]
        b.lr.CI = confint(model.lr)[2, ]
        model.nlr = nls(y ~ a1 * x ^ a2, start = list(a1 = a.lr, a2 = b.lr), dat = dat.sub,
                        control = nls.control(maxiter = 2000, warnOnly = TRUE))
        a.nlr = coef(summary(model.nlr))[1, 1]
        b.nlr = coef(summary(model.nlr))[2, 1]
        sd.nlr = sd(dat.sub$y - a.nlr * dat.sub$x ^ b.nlr)
        a.nlr.CI = confint2(model.nlr)[1, ]
        b.nlr.CI = confint2(model.nlr)[2, ]

        l.logn = sum(log(dlnorm(dat.sub$y, log(a.lr * dat.sub$x ^ b.lr), sd.lr)))
        l.norm = sum(log(dnorm(dat.sub$y, a.nlr * dat.sub$x ^ b.nlr, sd.nlr)))
        AICc.logn = 2 * k - 2 * l.logn + 2 * k * (k + 1) / (n - k - 1)
        AICc.norm = 2 * k - 2 * l.norm + 2 * k * (k + 1) / (n - k - 1)
        AICc.min = min(AICc.logn, AICc.norm)
        weight.logn = exp(-(AICc.logn - AICc.min)/2)
        weight.norm = exp(-(AICc.norm - AICc.min)/2)
        logn.w = weight.logn / (weight.logn + weight.norm)
        norm.w = weight.norm / (weight.logn + weight.norm)

        a.boot = a.lr * logn.w + a.nlr * norm.w
        b.boot = b.lr * logn.w + b.nlr * norm.w
        return(c(a.boot, b.boot))
      }
      dat.boot=boot(data = as.data.frame(cbind(x, y)), statistic = boot.est, R = 1000)
      a_confint = boot.ci(dat.boot, index = 1, type = "perc")$perc[4:5]
      b_confint = boot.ci(dat.boot, index = 2, type = "perc")$perc[4:5]
    }
  }
  writeLines(paste("a: ", a, "\nb: ", b))

  ## Step 3: residual plots as visual diagnostics
  if (diagno){
    lr_res = resid(model_lr)
    lr_pred = predict(model_lr)
    nlr_res = resid(model_nlr)
    nlr_pred = predict(model_nlr)
    par(mfrow = c(2, 2), oma = c(0, 4, 2, 0))
    hist(lr_res, freq = FALSE, main = "", xlab = "Residuals")
    curve(dnorm(x, mean(lr_res), sd(lr_res)), add = TRUE)
    plot(lr_pred, lr_res, xlab = "Predicted y", ylab = "Residuals")
    abline(h = 0)
    hist(nlr_res, freq = FALSE, main = "", xlab = "Residuals")
    curve(dnorm(x, mean(nlr_res), sd(nlr_res)), add = TRUE)
    plot(nlr_pred, nlr_res, xlab = "Predicted y", ylab = "Residuals")
    abline(h = 0)
    mtext("      Normality                                            Homoscedasticity",
          cex = 1.5, side = 3, outer = TRUE)
    mtext("NLR                                        LR", cex = 1.5, side = 2, outer = TRUE)
  }
  if (output_plot){
    par(mfrow = c(1, 2), oma = c(0, 4, 2, 0))
    plot(x, y, log = "xy", xlab = "x", ylab = "y", pch = 20, main = "Logarithmic Scale")
    curve(a * x ^ b, add = TRUE, lty = "dashed")
    curve(an * x ^ bn, add = TRUE, lty = "dotted")
    plot(x, y, xlab = "x", ylab = "y", pch = 20, main = "Arithmetic Scale")
    curve(a * x ^ b, add = TRUE, lty = "dashed")
    curve(an * x ^ bn, add = TRUE, lty = "dotted")
    title(main = paste("Fitting Power-Law Data with ", method), outer = TRUE)
  }

  return (list(method = method, a = a, b = b, a_confint = a_confint, b_confint = b_confint))
}

sar_fric_site <- sar_fric %>%
  group_by(site, x) %>%
  mutate(median_fric = median(y)) %>%
  distinct(median_fric, .keep_all= TRUE) %>%
  group_by(x) %>%
  mutate(window_rank = order(order(median_fric, decreasing = T)))

c_site <- c_results %>%
  group_by(site) %>%
  summarise(mean_c = mean(c)) %>%
  mutate(window_rank = order(order(mean_c, decreasing = T)),
         fake_area = 6e06)

ggplot(data = sar_fric_site, aes(x = as.factor(x), y = window_rank, color = site,
                                 group = site)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_text(aes(label = window_rank), vjust = 0.6, size = 3.5, color = 'black') +
  geom_point(data = c_site, aes(x = as.factor(fake_area), y = window_rank, color = site), size = 6, alpha = 0.6) +
  geom_text(data = c_site, aes(x = as.factor(fake_area), y = window_rank,label = window_rank), vjust = 0.6, size = 3.5, color = 'black') +
  geom_line(size = 0.3) +
  theme_bw(base_size = 20) +
  labs(x = expression('Area (m'^'2' *')'), y = "Rank Order (FRic)", color = "Site")

ggplot(data = sar_fric_site, aes(x = x, y = window_rank, color = site)) +
  geom_point() +
  geom_line(size = 0.3) +
  theme_bw(base_size = 12) +
  labs(x = expression('Area (m'^'2' *')'), y = "Rank Order (FRic)", color = "Site")

ggplot(data = sar_fric_site, aes(x = x, y = window_rank, color = site)) +
  geom_point() +
  geom_line(size = 0.3) +
  theme_bw(base_size = 12) +
  labs(x = "Area", y = "Rank Order (FRic)", color = "Site")

  geom_line(size = .3) +
  geom_point(size = 6) +
  geom_text(aes(label = rank), vjust = 0.6, size = 3.5, color = 'white') +
  scale_color_manual(values = my_colors) +
  scale_y_reverse(breaks = 1:nrow(df)) +
  labs(x = "", y = "") +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                legend.title = element_blank(), legend.position = 'bottom',    panel.border = element_blank()
  )  + scale_y_reverse(
    breaks = 1:n_distinct(df$car),
    labels = df %>% filter(year == 'year1') %>% arrange(rank) %>% pull(car),
    sec.axis = sec_axis(
      trans = I,
      breaks = 1:n_distinct(df$car),
      labels = df %>% filter(year == 'year4') %>%
        arrange(rank) %>% pull(car))) +
  scale_x_discrete(expand = c(0,.1))
  
## Rank order within sites
sar_fric_plot <- fric_results %>%
    group_by(site, plot, Window_Size) %>%
    mutate(median_fric = median(Hull_Volume)) %>%
    distinct(median_fric, .keep_all= TRUE) %>%
    ungroup() %>%
    group_by(site, Window_Size) %>%
    mutate(window_rank = order(order(median_fric, decreasing = T))) %>%
    mutate(Area = Window_Size^2)
  
c_rank_plot <- c_results %>%
    group_by(site, plot) %>%
    summarise(median_c = median(c)) %>%
    ungroup() %>%
    group_by(site) %>%
    mutate(window_rank = order(order(median_c, decreasing = T)),
           fake_area = 6e06)

sr_rank_plot <- sr %>%
  filter(sample_unit == 1) %>%
  group_by(site, plot) %>%
  summarise(median_g = median(gamma)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(window_rank = order(order(median_g, decreasing = T)))

rank_plot <- merge(sar_fric_plot, c_rank_plot, by = c("site", "plot"))

sar_fric_heal <- sar_fric_plot %>%
  filter(site =="HEAL")

c_rank_heal <- c_rank_plot %>%
  filter(site =="HEAL")

ggplot(data = sar_fric_heal, aes(x = as.factor(Area), y = window_rank, color = as.factor(plot),
                                 group = plot)) +
  geom_point(size = 6, alpha = 0.6) +
  geom_text(aes(label = window_rank), vjust = 0.6, size = 3.5, color = 'black') +
  geom_point(data = c_rank_heal, aes(x = as.factor(fake_area), y = window_rank, color = as.factor(plot)), size = 6, alpha = 0.6) +
  geom_text(data = c_rank_heal, aes(x = as.factor(fake_area), y = window_rank,label = window_rank), vjust = 0.6, size = 3.5, color = 'black') +
  geom_line(size = 0.3) +
  theme_bw(base_size = 20) +
  labs(x = expression('Area (m'^'2' *')'), y = "Rank Order (FRic)", color = "Plot")
