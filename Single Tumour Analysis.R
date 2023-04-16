rm(list = ls())

install.packages("devtools")
install.packages("ggpubr")
install.packages("patchwork")
install.packages("wesanderson")

library(dplyr)
library(usethis)
library(devtools)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(plotly)
library(cowplot)
library(boot)
library(ggthemr)
library(ggpubr)
library(pastecs)
library(gdata)
library(tidyr)
library(reshape2)
library(forcats)
library(scales)
library(magrittr)
library(graphics)
library(viridis)
library(ggsci)
library(patchwork)
library(wesanderson)

dapi_gt <- read.csv("/Users/peterqiang/Documents/ex vivo/EVD12/EVD12 Full Set GT Results.csv")
dapi_dl <- read.csv("/Users/peterqiang/Documents/ex vivo/EVD12/EVD12 Full Set I Results.csv")
well_template_UT <- read.csv("/Users/peterqiang/Documents/ex vivo/well_template_UT.csv")

pal <- wes_palette("Zissou1", 100, type = "continuous")
pal

str(dapi_gt)
str(dapi_dl)

evd_dapi <- lapply(list(dapi_gt, dapi_dl),
                   function(rn) rn %>%
                     rename(well = ImageSceneContainerName..Image.Scene.Container.Name,
                            ImageSceneName = ImageSceneName..Image.Scene.Name,
                            ParentID = ParentID..ID.of.the.parent..I,
                            ID = ID..ID..I,
                            Channel = ImageChannelName..Image.Channel.Name,
                            IMean_DAPI = IntensityMean_DAPI..Intensity.Mean.Value.of.channel..DAPI...R,
                            IMax_DAPI = IntensityMaximum_DAPI..Intensity.Maximum.of.channel..DAPI...R,
                            IMin_DAPI = IntensityMinimum_DAPI..Intensity.Minimum.of.channel..DAPI...R,
                            Isd_DAPI = IntensityStd_DAPI..Intensity.Standard.Deviation.of.channel..DAPI...R,
                            Area = Area..Area..R,
                            Circularity = Circularity..Circularity..R
                     ))
str(evd_dapi)

dapicount <- function(evd) {
  evd <- lapply(evd, function(rnan) na.omit(rnan))
  no.na <- lapply(evd, function(nnan) sum(is.na(nnan))) == 0 
  if (sum(no.na) == 2) {
    dapicount_evd <- lapply(evd, function(g) g %>%
                              group_by(well) %>%
                              dplyr::summarise(freq = n()))
    dapicount_evd[[1]] <- dapicount_evd[[1]] %>% rename(freq_gt = freq)
    dapicount_evd[[2]] <- dapicount_evd[[2]] %>% rename(freq_dl = freq)
    dapicount_evd[[1]]$freq_gt <- as.numeric(dapicount_evd[[1]]$freq_gt)
    dapicount_evd[[2]]$freq_dl <- as.numeric(dapicount_evd[[2]]$freq_dl)
  }
  dapicount_evd
}

dapi_evd <- dapicount(evd_dapi)
str(dapi_evd)
dapi.vs <- reduce(dapi_evd, inner_join, by = "well")
str(dapi.vs)

#plot(factor(dapicount_evd[[1]]$well), dapicount_evd[[1]]$freq_gt, 
#     pch = ".", col = "red", xlab = "well", ylab = "dapicount")
#points(factor(dapicount_evd[[1]]$well), dapicount_EVD[[2]]$freq_dl, 
#       pch = "*", col = "blue")

master.dapi <- merge(well_template_UT, dapi.vs, by = "well", all = FALSE) 
master.dapi[,c("row", 
               "col", 
               "order")] <- sapply(master.dapi[,c("row", "col", "order")], as.numeric)
str(master.dapi)

gtdl.heatmap <- function(master) {
  ht.gt <- ggplot(master, aes(let, row, fill = freq_gt))
  ht.dl <- ggplot(master, aes(let, row, fill = freq_dl))
  ht <- lapply(list(ht.gt, ht.dl), function(ht) ht + 
                 geom_tile(aes(x = let, y = row)) +
                 labs(x = "Column", y = "Row", fill = "DAPI count") +
                 theme_classic() +
                 scale_fill_gradientn(colours = pal, limits = c(0, 2500)) +
                 scale_x_discrete(position = "top") +
                 scale_y_reverse(breaks = c(1:24)) +
                 coord_fixed())
  ht.comp <- ggarrange(ht[[1]], ht[[2]])
  return(ht.comp)
}
heatmap.evd <- gtdl.heatmap(master.dapi)
heatmap.evd
ggsave2("Heatmap EVD DAPI counts Deep Learning vs Global Thresholding.png", dpi = 700)

#lin reg
cor.test(master.dapi$freq_gt, master.dapi$freq_dl)
lm_all <- lm(freq_dl~freq_gt, master.dapi)
summary(lm_all)
plot(master.dapi$freq_gt, 
     master.dapi$freq_dl, 
     xlab = "GT counts", ylab = "I counts",
     pch = "*", main = "Correlation between freq_GT and freq_I")
abline(lm_all, col = "red")

df_paired <- data.frame(
  Method = c(rep("Global Thresholding", times = 308), 
             rep("CNN Deep Learning", times = 308)),
  Frequency = c(master.dapi$freq_gt, master.dapi$freq_dl))
df_paired

box.jitter <- ggboxplot(df_paired, x = "Method", y = "Frequency",
                   color = "Method", 
                   add = "jitter", ylab = "DAPI counts") +
  scale_colour_manual(breaks = unique(df_paired$Method), 
                      values = c(pal[100], pal[1])) +
  scale_y_continuous(limits = c(0, 2500)) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, 
                     label.x = 1.35, label.y = 2500)
box.jitter

linebox <- ggpaired(df_paired, x = "Method", y = "Frequency", 
                      color = "Method", line.size = 0.1, line.color = "Grey",
                      xlab = "Method", ylab = "DAPI counts") +
  scale_colour_manual(breaks = unique(df_paired$Method), 
                      values = c(pal[100], pal[1])) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, 
                     label.x = 1.35, label.y = 2400)
linebox

#wilcox
wilcox.test(master.dapi$freq_gt, 
            master.dapi$freq_dl, 
            paired = TRUE, alternative = "two.sided")
#shapiro
shapiro.test(master.dapi$freq_gt)
shapiro.test(master.dapi$freq_dl)

#well template info
compound <- as.vector(unique(master.dapi$compound))
compound_info <- sapply(compound, function(d) sum(master.dapi$compound == d))
compound_info
compound

control <- c("DMSO",
             "Staurosporine",
             "Aphidicolin",
             "Media")
dapi.ctrl <- master.dapi %>% filter(compound %in% control)
str(dapi.ctrl)

drugs <- compound[!compound %in% control]
drugs

dapi.dmso <- dapi.ctrl %>% filter(compound == "DMSO")
stats_dmso_gt <- stat.desc(dapi.dmso$freq_gt)
stats_dmso_dl <- stat.desc(dapi.dmso$freq_dl)
stats_dmso_gt
stats_dmso_dl

dapi_drugs0 <- bind_rows(replicate(22, dapi.dmso, simplify = FALSE)) 
dapi_drugs0$compound <- rep(drugs, each = 28)
dapi_drugs0

drug_profile <- rbind(master.dapi %>% filter(! compound %in% control), 
                      dapi_drugs0)
dapi_drugs <- master.dapi %>% filter(!compound %in% control)
str(dapi_drugs)

################################################################################
#with 0 uM, no log trans
b <- ggplot(drug_profile, aes(x = conc, y = freq_gt)) + 
  facet_wrap(~ compound, ncol = 4, scales = "free_x") +
  stat_summary(aes(x = conc, y = freq_gt),
               fun = mean, geom ="line", size = 0.7) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2),
               geom = "errorbar", width = 0.05) +
  labs(x = "[uM]", y = "DAPI counts") +
  ggtitle("Global Thresholding Drug Response")
b
c <- ggplot(drug_profile, aes(x = conc, y = freq_I)) + 
  facet_wrap(~ compound, ncol = 4, scales = "free_x") +
  stat_summary(aes(x = conc, y = freq_dl),
               fun = mean, geom ="line", size = 0.7) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2),
               geom = "errorbar", width = 0.05) + 
  labs(x = "[uM]", y = "DAPI counts") +
  ggtitle("Intellesis Drug Response")
c
ggarrange(b, c)

gtdl.drp <- function(dp) {
  drp.gt <- ggplot(dp, aes(x = conc, y = freq_gt),
                   environment = environment()) + 
    stat_summary(mapping = aes(x = conc, y = freq_gt), 
                 fun = mean, geom ="line", size = 0.7) 
  drp.dl <- ggplot(dp, aes(x = conc, y = freq_dl), 
                   environment = environment()) + 
    stat_summary(mapping = aes(x = conc, y = freq_dl),
                 fun = mean, geom ="line", size = 0.7) 
  drp <- lapply(list(drp.gt, drp.dl), function (d) d +
                  facet_wrap(~ compound, ncol = 4, scales = "free_x") + 
                  stat_summary(fun.data = "mean_sdl", 
                               fun.args = list(mult = 2),
                               labs(x = "[uM]", y = "DAPI counts")) +
                  ggtitle("Intellesis Drug Response"))
  return(ggarrange(drp[[1]], drp[[2]]))
}
gtdl.drp(drug_profile)                              

#in one
drp_gtdl <- ggplot() +
  facet_wrap(~ compound, ncol = 4, scales = "free_x") +
  stat_summary(data = drug_profile, 
               aes(x = conc, y = freq_GT, colour = "Global Thresholding"),
               fun = mean, geom ="line", size = 0.7) +
  stat_summary(data = drug_profile, 
               aes(x = conc, y = freq_GT, colour = "Global Thresholding"),
               fun.data = "mean_sdl", fun.args = list(mult = 2),
               geom = "errorbar", width = 0.1) +
  stat_summary(data = drug_profile, 
               aes(x = conc, y = freq_I, colour = "CNN Deep Learning"),
               fun = mean, geom ="line", size = 0.7) +
  stat_summary(data = drug_profile, 
               aes(x = conc, y = freq_I, colour = "CNN Deep Learning"),
               fun.data = "mean_sdl", fun.args = list(mult = 2),
               geom = "errorbar", width = 0.1) + 
  theme(legend.position = c(0.8, 0.05)) +
  scale_colour_manual(values = c("Global Thresholding" = pal[100],
                                 "CNN Deep Learning" = pal[1])) +
  labs(x = "[uM]", y = "DAPI count", colour = "Model") +
  ggtitle("EVD12 Drug Response Profile")
drp_gtdl
ggsave2("EVD12 Drug Response Profile.png", dpi = 700)

#without 0 uM log transformed
#ggplot(dapi_drugs, aes(x = conc, y = freq_GT)) + 
#  facet_wrap(~ compound, ncol = 4, scales = "free_x") +
#  stat_summary(aes(x = conc, y = freq_GT),
#               fun = mean, geom ="line", size = 0.7) + 
#  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2),
#               geom = "errorbar", width = 0.02) +
#  scale_x_continuous(trans = "log10", breaks = dapi_drugs$conc) +
#  labs(x = "[log10 uM]", y = "DAPI counts") +
#  ggtitle("GT Drug Response")

#by drug, without 0 uM, log transformed
##side by side
comparison.drug <- function(count, cpd) {
  gt.drug <- ggplot(subset(count, compound == cpd), 
                    aes(x = conc, y = freq_gt)) + 
    stat_summary(aes(x = conc, y = freq_gt),
                 fun = mean, geom ="line", size = 0.7, colour = pal[100]) + 
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2),
                 geom = "errorbar", width = 0.02, colour = pal[100]) +
    ggtitle(paste0("GT Response Profile ", cpd))
  dl.drug <- ggplot(subset(count, compound == cpd),
                    aes(x = conc, y = freq_dl)) +
    stat_summary(aes(x = conc, y = freq_dl),
                 fun = mean, geom ="line", size = 0.7,  colour = pal[1]) + 
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 2),
                 geom = "errorbar", width = 0.02,  colour = pal[1]) +
    ggtitle(paste0("DL Response Profile ", cpd))
  drug <- lapply(list(gt.drug, dl.drug), function(d) d + 
                   labs(x = "[uM]", y = "DAPI count") +
                   scale_x_continuous(trans = "log10", breaks = c(0, 2.5, 5, 10, 20)) +
                   ylim(0, 2000))
  return(ggarrange(drug[[1]], drug[[2]]))
}
comparison.drug(master.dapi, "Cisplatin")

##in one
graph.drugcomp <- function(count, cpd) {
  ggplot() +
    stat_summary(data = subset(count, compound == cpd), 
               aes(x = conc, y = freq_gt, colour = "Global Thresholding"),
               fun = mean, geom ="line", size = 0.7) +
    stat_summary(data = subset(count, compound == cpd), 
                 aes(x = conc, y = freq_gt, colour = "Global Thresholding"),
                 fun.data = "mean_sdl", fun.args = list(mult = 2),
                 geom = "errorbar", width = 0.1) +
    stat_summary(data = subset(count, compound == cpd), 
                 aes(x = conc, y = freq_dl, colour = "Deep Learning"),
                 fun = mean, geom ="line", size = 0.7) +
    stat_summary(data = subset(count, compound == cpd), 
                 aes(x = conc, y = freq_dl, colour = "Deep Learning"),
                 fun.data = "mean_sdl", fun.args = list(mult = 2),
                 geom = "errorbar", width = 0.1) +
    scale_x_continuous(trans = "log10",
                       breaks = subset(master.dapi, 
                                       compound == cpd)$conc) +  
    scale_colour_manual(values = c("Global Thresholding" = pal[100],
                                   "Deep Learning" = pal[1])) +
    labs(x = "[uM]", y = "DAPI count") +
    ggtitle(paste0(cpd, " Response Profile GT vs DL"))
}
graph.drugcomp(master.dapi, "Cisplatin")

################################################################################
#GR metrics
install.packages("BiocManager")
BiocManager::install("GRmetrics")
browseVignettes("GRmetrics")
library(GRmetrics)
## GT Results
mDMSO_GT <- mean(dapi.dmso$freq_gt)
mAphi_GT <- mean(filter(dapi.ctrl, compound == "Aphidicolin")$freq_gt)
gr_table_GT <- dapi_drugs
names(gr_table_GT)[names(gr_table_GT) == "conc"] <- "concentration"
names(gr_table_GT)[names(gr_table_GT) == "freq_gt"] <- "cell_count"

gr_table_GT$time <- rep(96, each = 264)
gr_table_GT$cell_count__ctrl <- rep(mDMSO_GT, each = 264)
gr_table_GT$cell_count__time0 <- rep(mAphi_GT, each = 264) 
str(gr_table_GT)

drc_output_GT = GRfit(gr_table_GT, groupingVariables = c('compound'))
View(GRgetMetrics(drc_output_GT))
View(GRgetDefs(drc_output_GT))
View(GRgetValues(drc_output_GT))

write.csv(GRgetMetrics(drc_output_GT), file = "EVD GRMetrics GT DAPI.csv")

GRdrawDRC(drc_output_GT)
## I Results
mDMSO_DL <- mean(dapi.dmso$freq_dl)
mAphi_DL <- mean(filter(dapi.ctrl, compound == "Aphidicolin")$freq_dl)
gr_table_DL <- dapi_drugs
names(gr_table_DL)[names(gr_table_DL) == "conc"] <- "concentration"
names(gr_table_DL)[names(gr_table_DL) == "freq_dl"] <- "cell_count"

gr_table_DL$time <- rep(96, each = 264)
gr_table_DL$cell_count__ctrl <- rep(mDMSO_DL, each = 264)
gr_table_DL$cell_count__time0 <- rep(mAphi_DL, each = 264)
str(gr_table_DL)

drc_output_DL = GRfit(gr_table_DL, groupingVariables = c('compound'))
View(GRgetMetrics(drc_output_DL))
View(GRgetDefs(drc_output_DL))
View(GRgetValues(drc_output_DL))

write.csv(GRgetMetrics(drc_output_DL), file = "EVD GRMetrics I DAPI.csv")

GRdrawDRC(drc_output_DL)

#IC50
GT_IC50 <- GRgetMetrics(drc_output_GT)
GT_IC50$minus_log10IC50 <- -(log10(GT_IC50$IC50))
GT_IC50$minus_GR50 <- -(log10(GT_IC50$GR50)) 
GT_IC50$compound <- as.factor(GT_IC50$compound)
str(GT_IC50)

DL_IC50 <- GRgetMetrics(drc_output_DL)
DL_IC50$minus_log10IC50 <- -(log10(DL_IC50$IC50))
DL_IC50$minus_GR50 <- -(log10(DL_IC50$GR50))
DL_IC50$compound <- as.factor(DL_IC50$compound)
str(DL_IC50)

#AUC comparison
GT_IC50 %>%
  mutate(compound = fct_reorder(compound, AUC)) %>%
  ggplot(aes(compound, 1/AUC)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("EVD waterfall plot by 1/AUC, Global Thresholding") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DL_IC50 %>%
  mutate(compound = fct_reorder(compound, AUC)) %>%
  ggplot(aes(compound, 1/AUC)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("EVD waterfall plot by 1/AUC, Deep Learning") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


auc_gt <- GT_IC50 %>% mutate(compound = fct_reorder(compound, AUC))
auc_dl <- DL_IC50 %>% mutate(compound = fct_reorder(compound, AUC))
auc_gtvi <- data.frame(compound = fct_reorder(DL_IC50$compound, DL_IC50$AUC),
                       auc_dl = DL_IC50$AUC,
                       auc_gt = GT_IC50$AUC) %>%
  melt(id = "compound", value.name = "AUC", variable.name = "Model")
ggplot(auc_gtvi, 
       aes(x = compound, AUC, y = 1/AUC, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("auc_dl" = pal[1],
                               "auc_gt" = pal[100]),
                    labels = c("CNN Deep Learning", "Global Thresholding")) +
  labs(fill = "Model") +
  ggtitle("EVD12 waterfall plot by 1/AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave2("EVD12 waterfall plot by 1/AUC I vs GT.png", dpi = 700)

cor <- ggplot(data = data.frame(x = GT_IC50$AUC, y = DL_IC50$AUC), 
                aes(x = x, y = y)) +
  geom_point(colour = pal[50]) +
  ylim(0.4, 1) +
  xlim(0.4, 1) +
  xlab("AUC by Global Thresholding") +
  ylab("AUC by CNN Deep Learning") +
  stat_cor(method = "pearson") +
  ggtitle("EVD")
cor
cor.test(GT_IC50$AUC, DL_IC50$AUC)

#not run
GT_IC50 %>%
  mutate(compound = fct_reorder(compound, dplyr::desc(minus_log10IC50))) %>%
  ggplot(aes(compound, minus_log10IC50)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("waterfall plot by -log10IC50, GT") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DL_IC50 %>%
  mutate(compound = fct_reorder(compound, dplyr::desc(minus_log10IC50))) %>%
  ggplot(aes(compound, minus_log10IC50)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("waterfall plot by -log10IC50, DL") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

GT_IC50 %>%
  mutate(compound = fct_reorder(compound, dplyr::desc(minus_GR50))) %>%
  ggplot(aes(compound, minus_GR50)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("waterfall plot by -log10GR50, GT") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DL_IC50 %>%
  mutate(compound = fct_reorder(compound, dplyr::desc(minus_GR50))) %>%
  ggplot(aes(compound, minus_GR50)) +
  geom_bar(stat="identity", position = "dodge") +
  ggtitle("waterfall plot by -log10GR50, DL") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

inv_AUC12 <- data.frame(invAUC12_GT = 1/GT_IC50$AUC,
                        invAUC12_DL = 1/DL_IC50$AUC)

################################################################################
#exploratory
ggplot(data = master.dapi) +
  geom_point(aes(x = fct(well), y = freq_gt, color = "GT"), size = 1) +
  geom_point(aes(x = fct(well), y = freq_dl, color = "DL"), size = 1) +
  scale_colour_manual(values = c("GT" = pal[100],
                                 "DL" = pal[1])) +
  labs(x = "well", y = "DAPI count")

master.dapi$difference <- master.dapi$freq_I - master.dapi$freq_GT
write.csv(stat.desc(master.dapi$difference), file = "differences stats.csv")

sum_diff <- sum(master.dapi$difference)
mean_diff <- mean(master.dapi$difference)
sd_diff <- sd(master.dapi$difference)
med_diff <- median(master.dapi$difference)
iqr_diff <- quantile(master.dapi$difference)

ggplot(data = master.dapi, aes(x = difference)) +
  geom_density(fill = "lightgray") +
  stat_overlay_normal_density(color = "red")

#mean +/- 2sd
ggplot(data = master.dapi, aes(x = fct(well), y = difference, colour = difference)) +
  geom_point() +
  scale_colour_gradient(low = pal[100], high = pal[50]) +
  labs(x = "well", y = "difference (DL - GT)")

#median and iqr
ggplot(data = master.dapi, aes(x = fct(well), y = difference, 
                             color = difference)) +
  geom_point() +
  labs(x = "well", y = "mod_trans(difference (DL - GT))") +
  scale_colour_gradient(low = pal[100], high = pal[50])

master.dapi$diff_sq <- master.dapi$difference^2

stat.desc(master.dapi$diff_sq)

sum_diff_sq <- sum(master.dapi$diff_sq)
mean_diff_sq <- mean(master.dapi$diff_sq)
sd_diff_sq <- sd(master.dapi$diff_sq)
med_diff_sq <- median(master.dapi$diff_sq)
iqr_diff_sq <- quantile(master.dapi$diff_sq)

#mean +/- 2sd
ggplot(data = master.dapi, aes(x = fct(well), y = diff_sq)) +
  geom_point() +
  labs(x = "well", y = "mod_trans(squared difference (DL - GT))") +
  geom_hline(yintercept = mean_diff_sq, colour = "green") +
  scale_y_continuous(trans = modulus_trans(0))

#median and iqr
ggplot(data = master.dapi, aes(x = fct(well), y = diff_sq)) +
  geom_point() +
  labs(x = "well", y = "mod_trans(squared difference (DL - GT))") +
  geom_hline(yintercept = med_diff_sq, colour = "blue") +
  geom_hline(yintercept = iqr_diff_sq[[2]], colour = "yellow") +
  geom_hline(yintercept = iqr_diff_sq[[4]], colour = "yellow") +
  scale_y_continuous(trans = modulus_trans(0))

ggplot(data = master.dapi, aes(x = fct(well), y = diff_sq, colour = diff_sq)) +
  geom_point() +
  scale_colour_gradientn(colours = pal)

widediff <- DAPIcount %>% filter(diff_sq > mean_diff_sq)
widediff[order(widediff$diff_sq, decreasing = TRUE),]




































