### Analyze data from CC mice (SARS1 exposure study) - sent from Sarah in Prism format
library(pzfx)
library(tidyverse)
library(readxl)
#library(reshape2)
#library(patchwork)
library(ggpubr)
library(cowplot)
source('code-dependencies/qtl_functions.R')

load_themes()

strain_ids <- c("16211" = "CC005",
                "16750" = "CC006",
                "3252" = "CC011",
                "8024" = "CC016",
                "16513" = "CC019",
                "8054" = "CC020",
                "8043" = "CC023",
                "8004" = "CC024",
                "4410" = "CC044",
                "16768" = "CC068")

custom_colors <- scale_color_manual(values = c("CC005" = "#F781BF",
                                               "CC006" = "#E41A1C",
                                               "CC011" = "#984EA3",
                                               "CC016" = "#FF7F00",
                                               "CC019" = "#000000",
                                               "CC020" = "#A65628",
                                               "CC023" = "#4DAF4A",
                                               "CC024" = "#999999",
                                               "CC044" = "#377EB8",
                                               "CC068" = "#98F5FF"))

#-----------------------------------Weight-------------------------------------#
CC_weight <- read_pzfx("source_data/screen/CC_mice.pzfx", "CC_weight")

for (i in c(1:length(strain_ids))){
  cc_id <- names(strain_ids)[i]
  cc_strain <- strain_ids[cc_id]
  df <- CC_weight[,which(str_detect(colnames(CC_weight), cc_id))]
  
  # data formatted for analysis 
  df2 <- t(df)
  colnames(df2) <- c(0,1,2,3,4)
  df2 <- as_tibble(df2, rownames = 'mouse_ID') 
  df2$Strain = cc_strain
  # combine data
  if (i==1){
    analysisdf <- df2
  } else {
    analysisdf <- rbind(analysisdf, df2)
  }
  
  # data formatted for plotting (just means and sds)
  df$mean <- rowMeans(df, na.rm = T)
  df$sd <- apply(df[1:4], MARGIN = 1, FUN = sd, na.rm = T)
  df$dpi = seq(0,4)
  df$Strain = cc_strain
  # combine data 
  if (i==1){
    plotdf <- df[,c(5:8)]
  } else {
    plotdf <- rbind(plotdf, df[,c(5:8)])
  }
}

### Plotting 
plotdf$Strain <- as.factor(plotdf$Strain)

# plot avg weight trajectories of all 10 CC strains 
wt_plot <- ggplot(data = plotdf, mapping = aes(x = dpi, y = mean, color = Strain)) + 
  geom_point() + 
  geom_line(size = 1) + 
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) + 
  labs(x = "Days post infection", y = "% of starting weight", title = "") + 
  ylim(77,102) + 
  custom_colors +
  bw_big_theme 
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_wt.png", width = 7, height = 5)
# 
# ggplot(plotdf, aes(x = dpi, y = mean, color = Strain)) + 
#   geom_point(aes(alpha = Strain %in% c("CC006", "CC044"))) + 
#   geom_line(aes(alpha = Strain %in% c("CC006", "CC044")), size = 1) + 
#   labs(x = "Days post infection", y = "% of starting weight", title = "") + 
#   ylim(77, 102) + 
#   custom_colors +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
#   bw_big_theme
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_wt2.png", width = 7, height = 5)

ensure_directory('figures/screen/')
png("figures/screen/parent_weight.png", width = 600)
print(wt_plot)
dev.off()

# plot avg weight trajectories of CC006 and CC044 mice 
f2_plotdf <- plotdf[which(plotdf$Strain %in% c('CC006', 'CC044')),]

wt_plot_f2 <- ggplot(data = f2_plotdf, mapping = aes(x = dpi, y = mean, color = Strain)) + 
  geom_point() + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), size=1, width=0.3, alpha=0.4) + 
  labs(x = "Days post infection", y = "% of starting weight", title = "") + 
  ylim(77,102) +
  bw_big_theme

png("figures/screen/f2_parent_weight.png", width = 600)
print(wt_plot_f2)
dev.off()


### Analysis 

# How much weight do CC006 mice lose on average? 
cc006id <- "16750"
cc006df <- CC_weight[,which(str_detect(colnames(CC_weight), cc006id))]
cc006pct <- 1-(cc006df[5,]/cc006df[1,])
# 2 mice died by 3dpi 
mean(as.numeric(cc006pct)[1:2]) # mean weight loss in 2 mice surviving to d4

# How much weight do CC044 mice lose on average? 
cc044id <- "4410"
cc044df <- CC_weight[,which(str_detect(colnames(CC_weight), cc044id))]
cc044pct <- 1-(cc044df[5,]/cc044df[1,])
mean(as.numeric(cc044pct))

# calculate auc statistic for Tukey's HSD
analysisdf$Strain <- as.factor(analysisdf$Strain)
adf <- na.omit(analysisdf)

for (i in 1:nrow(adf)){
  line <- adf[i,c('0', '1', '2', '3', '4')]
  auc <- auc(x = c(0:4), y = line)
  aac <- line[1]*4 - auc # 400 bc line[1] is always =100
  adf[i,'weight_aac'] <- aac 
}

fit <- aov(weight_aac ~ Strain, adf)
hsd <- TukeyHSD(fit)
print(hsd)



#-------------------------------------HS---------------------------------------#
CC_HS <- read_pzfx("source_data/screen/CC_mice.pzfx", "CC_HS")

for (i in c(1:length(strain_ids))){
  cc_id <- names(strain_ids)[i]
  df <- CC_HS[,which(str_detect(colnames(CC_HS), cc_id))]
  strainID <- str_split(colnames(df)[1], "_")[[1]][1]
  df$Strain = strain_ids[strainID]
  
  dflong <- pivot_longer(df, cols = c(1:4), names_to = "rep", values_to = "HS")
  
  # combine mean, sd, strain for plotting 
  if (i==1){
    hsdf <- dflong
  } else {
    hsdf <- rbind(hsdf, dflong)
  }
}

hsdf$Strain <- as.factor(hsdf$Strain)

# create summary table 
hs_sum <- hsdf %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(HS, na.rm=TRUE), sd = sd(HS, na.rm=TRUE))

# plot HS for all 10 strains 
hs_plot <- ggplot(data = hsdf, mapping = aes(x = Strain, y = HS, color = Strain)) + 
  geom_jitter(width=0.1, size=2) + 
  geom_errorbar(data = hs_sum, aes(x = Strain, y = mean, ymin = mean, ymax = mean), size=1, width=0.3, alpha=0.4) + 
  geom_errorbar(data = hs_sum, aes(x = Strain, y = mean, ymin = mean-sd, ymax = mean+sd), size=1, width=0.5, alpha=0.4) +
  labs(x = "Strain", y = "Congestion score", title = "") + 
  theme(legend.position = "none") +
  custom_colors +
  ylim(-0.2,4.5) +
  bw_big_theme2 + 
  guides(x = guide_axis(angle = 30))

png("figures/screen/parent_hs.png")
print(hs_plot)
dev.off()

# p_jit <- position_jitter(width = 0.1, height = 0, seed = 4)
# ggplot(hsdf, aes(x = Strain, y = HS, color = Strain)) + 
#   geom_point(position = p_jit, size = 2) +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.5, size = 1) +
#   stat_summary(fun = mean, geom = "point", size = 3) +
#   theme(legend.position = "none") +
#   custom_colors +
#   coord_cartesian(ylim = c(-0.2, 4.5)) +
#   bw_big_theme2 +
#   theme(axis.title = element_text(size = 20), 
#         axis.text = element_text(size = 18)) +
#   guides(x = guide_axis(angle = 30)) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
#   labs(x = "Strain", y = "Congestion score", title = "")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_hs.png", width = 5, height = 5)
# 
# ggplot(hsdf, aes(x = Strain, y = HS, color = Strain)) + 
#   geom_point(aes(alpha = Strain %in% c("CC006", "CC044")), position = p_jit, size = 2) +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
#                aes(alpha = Strain %in% c("CC006", "CC044")), width = 0.5, size = 1) +
#   stat_summary(fun = mean, geom = "point", aes(alpha = Strain %in% c("CC006", "CC044")), size = 3) +
#   theme(legend.position = "none") +
#   custom_colors +
#   coord_cartesian(ylim = c(-0.2, 4.5)) +
#   bw_big_theme2 +
#   theme(axis.title = element_text(size = 20), 
#         axis.text = element_text(size = 18)) +
#   guides(x = guide_axis(angle = 30)) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
#   labs(x = "Strain", y = "Congestion score", title = "")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_hs2.png", width = 5, height = 5)

# plot HS for F2 parent strains
hs_f2 <- hsdf[which(hsdf$Strain %in% c('CC006', 'CC044')),]

hs_sum_f2 <- hs_f2 %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(HS, na.rm=TRUE), sd = sd(HS, na.rm=TRUE))

hs_plot_f2 <- ggplot(data = hs_f2, mapping = aes(x = Strain, y = HS, color = Strain)) + 
  geom_jitter(width=0.1, size=2) + 
  geom_errorbar(data=hs_sum_f2, aes(x = Strain, y = mean, ymin = mean, ymax = mean), size=1, width=0.2, alpha=0.4) + 
  geom_errorbar(data=hs_sum_f2, aes(x = Strain, y = mean, ymin = mean-sd, ymax = mean+sd), size=1, width=0.3, alpha=0.4) +
  labs(x = "Strain", y = "Congestion score", title = "") + 
  theme(legend.position = "none") +
  custom_colors +
  ylim(-0.2,4.5) +
  bw_big_theme

png("figures/screen/f2_parent_hs.png", width = 280)
print(hs_plot_f2)
dev.off()

# statistical test 
fit <- aov(HS ~ Strain, data = hsdf)
hsd <- TukeyHSD(fit)
print(hsd)



#-----------------------------------Titer--------------------------------------#
CC_titer <- read_pzfx("source_data/screen/CC_mice.pzfx", "CC_titer")

for (i in c(1:length(strain_ids))){
  cc_id <- names(strain_ids)[i]
  df <- CC_titer[,which(str_detect(colnames(CC_titer), cc_id))]
  strainID <- str_split(colnames(df)[1], "_")[[1]][1]
  df$Strain = strain_ids[strainID]
  
  dflong <- pivot_longer(df, cols = c(1:4), names_to = "rep", values_to = "Titer")
  dflong$Titer <- log10(dflong$Titer+1)
  
  # combine mean, sd, strain for plotting 
  if (i==1){
    titerdf <- dflong
  } else {
    titerdf <- rbind(titerdf, dflong)
  }
}

titerdf$Strain <- as.factor(titerdf$Strain)

# create summary table 
titer_sum <- titerdf %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(Titer, na.rm=TRUE), sd = sd(Titer, na.rm=TRUE))

# Plot titer for all 10 CC strains 
titer_plot <- ggplot(data=titerdf, mapping=aes(x=Strain, y=Titer, color=Strain)) +
  geom_jitter(width=0.1, size=2) + 
  geom_errorbar(data=titer_sum, aes(x=Strain, y=mean, ymin=mean, ymax=mean), size=1, width=0.3, alpha=0.4) + 
  geom_errorbar(data=titer_sum, aes(x=Strain, y=mean, ymin=mean-sd, ymax=mean+sd), size=1, width=0.5, alpha=0.4) +
  labs(y = "Log titer (PFU/lobe)", title = "") + 
  theme(legend.position = "none") +
  ylim(-1,7.1) + 
  custom_colors +
  bw_big_theme2 + 
  guides(x = guide_axis(angle = 30))

png("figures/screen/parent_titer.png")
print(titer_plot)
dev.off()


# p_jit <- position_jitter(width = 0.1, height = 0, seed = 6)
# ggplot(titerdf, aes(x = Strain, y = Titer, color = Strain)) + 
#   geom_point(position = p_jit, size = 2) +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.5, size = 1) +
#   stat_summary(fun = mean, geom = "point", size = 3) +
#   theme(legend.position = "none") +
#   custom_colors +
#   coord_cartesian(ylim = c(-1, 7.1)) +
#   bw_big_theme2 +
#   theme(axis.title = element_text(size = 20), 
#         axis.text = element_text(size = 18)) +
#   guides(x = guide_axis(angle = 30)) +
#   labs(x = "Strain", y = "Log titer (PFU/lobe)", title = "")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_titer.png", width = 5, height = 5)
# 
# ggplot(titerdf, aes(x = Strain, y = Titer, color = Strain)) + 
#   geom_point(aes(alpha = Strain %in% c("CC006", "CC044")), position = p_jit, size = 2) +
#   stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
#                aes(alpha = Strain %in% c("CC006", "CC044")), width = 0.5, size = 1) +
#   stat_summary(fun = mean, geom = "point", aes(alpha = Strain %in% c("CC006", "CC044")), size = 3) +
#   theme(legend.position = "none") +
#   custom_colors +
#   coord_cartesian(ylim = c(-1, 7.1)) +
#   bw_big_theme2 +
#   theme(axis.title = element_text(size = 20), 
#         axis.text = element_text(size = 18)) +
#   guides(x = guide_axis(angle = 30)) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
#   labs(x = "Strain", y = "Log titer (PFU/lobe)", title = "")
# ggsave("~/Documents/ValdarFerris/defense/figures/chapter3/screen_titer2.png", width = 5, height = 5)



# plot titer for F2 parent CC strains 
titer_f2 <- titerdf[which(hsdf$Strain %in% c('CC006', 'CC044')),]

titer_sum_f2 <- titer_f2 %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(Titer, na.rm=TRUE), sd = sd(Titer, na.rm=TRUE))

titer_plot_f2 <- ggplot(data=titer_f2, mapping=aes(x=Strain, y=Titer, color=Strain)) +
  geom_jitter(width=0.1, size=2) + 
  geom_errorbar(data=titer_sum_f2, aes(x=Strain, y=mean, ymin=mean, ymax=mean), size=1, width=0.2, alpha=0.4) + 
  geom_errorbar(data=titer_sum_f2, aes(x=Strain, y=mean, ymin=mean-sd, ymax=mean+sd), size=1, width=0.3, alpha=0.4) +
  labs(y = "Log titer (PFU/lobe)", title = "") + 
  theme(legend.position = "none") +
  ylim(0,7) +
  custom_colors +
  bw_big_theme

png("figures/screen/f2_parent_titer.png", width = 280)
print(titer_plot_f2)
dev.off()

# statistical test 
cc006titer <- as.numeric(CC_titer[,which(str_detect(colnames(CC_titer), cc006id))])
cc044titer <- as.numeric(CC_titer[,which(str_detect(colnames(CC_titer), cc044id))])
wilcox.test(x = cc006titer, y = cc044titer)
t.test(cc006titer, cc044titer) # assumes normality, data are not normal

# Tukey's HSD 
fit <- aov(Titer ~ Strain, data = titerdf)
hsd <- TukeyHSD(fit)
print(hsd)

#----------------------------------Figure--------------------------------------#
png("figures/screen/screen-figure.png", height = 400, width = 1300)
plot_grid(wt_plot, hs_plot, titer_plot, align = "h", ncol = 3, 
          rel_widths = c(0.42, 0.3, 0.3), labels = c('A', 'B', 'C'), label_size = 24)
dev.off()

png("figures/screen/screen-figure_f2.png", height = 400, width = 900)
plot_grid(wt_plot_f2, hs_plot_f2, titer_plot_f2, align = "h", ncol = 3, 
          rel_widths = c(0.5, 0.25, 0.25), labels = c('A', 'B', 'C'), label_size = 24)
dev.off()










#-----------------------------Repeat experiment--------------------------------#
repdat <- read_xlsx("source_data/screen/CC006_CC044_repeatinfection.xlsx", 
                    sheet = "forR")

repwt <- repdat[,c('Cage', 'Strain', 'Animal', 'd0', 'd1', 'd2', 'd3', 'd4')]
reptiter <- repdat[,c('Cage', 'Strain', 'Animal', 'Titer')]
  
rephs <- read_pzfx("source_data/screen/CC006_CC044_repeatinfection.pzfx", table = "HS")
# use bw_big_theme2

#-----------------------------------Weight-------------------------------------#
repwt[repwt == 'DIC'] <- NA # get rid of DIC
repwt$d4 <- as.numeric(repwt$d4)
# Calculate percentage of starting weight 
repwt$pd0 <- (repwt$d0/repwt$d0)*100
repwt$pd1 <- (repwt$d1/repwt$d0)*100
repwt$pd2 <- (repwt$d2/repwt$d0)*100
repwt$pd3 <- (repwt$d3/repwt$d0)*100
repwt$pd4 <- (repwt$d4/repwt$d0)*100

# pivot data long
colnames(repwt)[9:13] <- c(0,1,2,3,4)
repwtlong <- pivot_longer(repwt, cols = c('0', '1', '2', '3', '4'), 
                          names_to = 'dpi', values_to = 'pct_weight')
repwtlong <- repwtlong %>% dplyr::select(-c('d0', 'd1', 'd2', 'd3', 'd4')) # remove raw weight columns

# calculate averages and sds
repwtsum <- repwtlong %>% 
  group_by(Strain, dpi) %>%
  dplyr::summarise(mean = mean(pct_weight, na.rm = TRUE), 
                   sd = sd(pct_weight, na.rm = TRUE), 
                   se = sd(pct_weight, na.rm = TRUE)/sqrt(n()))

# plot avg weight trajectories of CC006, CC023, CC044 
repwtplot <- ggplot(data = repwtsum, mapping = aes(x = dpi, y = mean, group = Strain, color = Strain)) + 
  geom_point() + 
  geom_line(size = 1) + 
  # geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) + 
  labs(x = "Days post-infection", y = "% of starting weight", title = "") + 
  ylim(77,102) +
  custom_colors +
  bw_big_theme 

png("figures/screen/rep_parent_weight.png", width = 600)
print(repwtplot)
dev.off()


# calculate auc statistic for Tukey's HSD
repwt$Strain <- as.factor(repwt$Strain)
arepwt <- na.omit(repwt)
arepwt$weight_aac <- rep(NA, nrow(arepwt))

for (i in 1:nrow(arepwt)){
  line <- arepwt[i,c('0', '1', '2', '3', '4')]
  auc <- auc(x = c(0:4), y = line)
  aac <- line[1]*4 - auc # 400 bc line[1] is always =100
  arepwt[i,'weight_aac'] <- aac 
}

fit <- aov(weight_aac ~ Strain, arepwt)
hsd <- TukeyHSD(fit)
print(hsd)

#-------------------------------------HS---------------------------------------#
colnames(rephs) <- c('CC006', 'CC023', 'CC044')
rephslong <- pivot_longer(rephs, cols = c(1:3), names_to = "Strain", values_to = "HS") %>%
  arrange(Strain) %>% 
  filter((Strain == 'CC006') | (Strain %in% c('CC023', 'CC044') & !is.na(HS))) %>% # remove NAs except for CC006 bc those were mice that died  
  add_row(Strain = 'CC044', HS = NA) # one CC044 NA was real 

# create summary table 
rephssum <- rephslong %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(HS, na.rm=TRUE), sd = sd(HS, na.rm=TRUE))

# plot
rephsplot <- ggplot(data = rephslong, mapping = aes(x = Strain, y = HS, color = Strain)) + 
  geom_jitter(width=0.1, size=2) + 
  geom_errorbar(data = rephssum, aes(x = Strain, y = mean, ymin = mean, ymax = mean), size=1, width=0.3, alpha=0.4) + 
  geom_errorbar(data = rephssum, aes(x = Strain, y = mean, ymin = mean-sd, ymax = mean+sd), size=1, width=0.5, alpha=0.4) +
  labs(x = "Strain", y = "Congestion score", title = "") + 
  theme(legend.position = "none") +
  ylim(-0.2,4.5) +
  custom_colors +
  bw_big_theme2 
# guides(x = guide_axis(angle = 30))
# ylim(-0.2,4.5) 

png("figures/screen/rep_parent_hs.png")
print(rephsplot)
dev.off()

# statistical test 
fit <- aov(HS ~ Strain, data = rephslong)
hsd <- TukeyHSD(fit)
print(hsd)

#-----------------------------------Titer--------------------------------------#
#Log-transform titer
reptiter$Titer <- log(reptiter$Titer+1, base=10)

# create summary table 
reptitersum <- reptiter %>% 
  group_by(Strain) %>% 
  dplyr::summarize(mean = mean(Titer), sd = sd(Titer))

# plot HS for all 10 strains 
reptiterplot <- ggplot(data = reptiter, mapping = aes(x = Strain, y = Titer, color = Strain)) + 
  geom_point(size=2) + 
  geom_errorbar(data = reptitersum, aes(x = Strain, y = mean, ymin = mean, ymax = mean, color = Strain), 
                size=1, width=0.3, alpha=0.4) +
  geom_errorbar(data = reptitersum, aes(x = Strain, y = mean, ymin = mean-sd, ymax = mean+sd), 
                size=1, width=0.5, alpha=0.4) +
  labs(x = "Strain", y = "Log titer (PFU/lobe)", title = "") +
  ylim(-1,7.1) + 
  custom_colors +
  theme(legend.position = "none") + 
  bw_big_theme2 

png("figures/screen/rep_parent_titer.png")
print(reptiterplot)
dev.off()

# statistical test 
fit <- aov(Titer ~ Strain, data = reptiter)
hsd <- TukeyHSD(fit)
print(hsd)


#----------------------------------Figure--------------------------------------#
png("figures/screen/rep-screen-figure.png", height = 400, width = 1300)
plot_grid(repwtplot, rephsplot, reptiterplot, align = "h", ncol = 3, 
          rel_widths = c(0.42, 0.3, 0.3), labels = c('A', 'B', 'C'), label_size = 24)
dev.off()



ogscreen <- plot_grid(wt_plot, hs_plot, titer_plot, align = "h", ncol = 3,
                      rel_widths = c(0.42, 0.3, 0.3), label_size = 24,
                      labels = c('A', 'B', 'C'))

ogtitle <- ggdraw() + 
  draw_label("Initial screen", fontface = 'bold', x = 0, hjust = 0, size = 24) + 
  theme(plot.margin = margin(0,0,0,7))

rescreen <- plot_grid(repwtplot, rephsplot, reptiterplot, align = "h", ncol = 3, 
                      rel_widths = c(0.42, 0.3, 0.3), label_size = 24,
                      labels = c('D', 'E', 'F'))

retitle <- ggdraw() + 
  draw_label("Repeat experiment", fontface = 'bold', x = 0, hjust = 0, size = 24) + 
  theme(plot.margin = margin(0,0,0,7))



png("figures/screen/both-screens-v2.png", heigh = 850, width = 1300)
plot_grid(ogtitle, ogscreen, retitle, rescreen, ncol = 1, rel_heights = c(0.1,1,0.1,1))
dev.off()




