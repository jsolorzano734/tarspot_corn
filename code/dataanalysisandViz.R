#Tar spot of corn inoculation method
#========================packages===============================================
library("tidyverse") #includes ggplot and for data management
library("rcompanion") #to check normality & transformation
library("olsrr") #to check normality
library("multcompView") #to use multcompLetters4
library("agricolae")
library("rstatix") #rstatix provides pipe-friendly R functions for easy statistical analyses.
library("MetBrewer") #colors

#===========================data================================================
#load data
setwd("~/inoculationTarSpot")
signs <- read.csv("data.csv") #will be provided
str(signs)
summary(signs)

signs$hybrid <- factor(signs$hybrid,levels = c("goldCountry", 
                                               "nk9610", 
                                               "nk9653", 
                                               "sweetCorn"), 
                       labels = c("H3", "H2", "H1", "H4"))

signs$place <- factor(signs$place , levels = c("greenHouse", "growthChamber"), 
                      labels = c("Greenhouse", "Growth Chamber"))

signs$essay <- factor(signs$essay , levels = c("trial1", "trial2", "trial3", "trial4"), 
                      labels = c("Experiment 1", "Experiment 2", "Experiment 3", "Experiment 4"))

#=========================stats for trial 1=====================================
#stats trial1_stats
trial1_stats <- signs[which(signs$essay %in% c('Experiment 1')),]
plotNormalHistogram(trial1_stats$stromata)
shapiro.test(trial1_stats$stromata)

#transformation
trial1_stats.signs.t <- transformTukey(trial1_stats$stromata)
trial1_stats.signs.df.t_pre <- cbind(trial1_stats, trial1_stats.signs.t)
str(trial1_stats.signs.df.t_pre)
#View(trial1_stats.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
trial1_stats.signs.df.t.lm <- lm(trial1_stats.signs.t ~ place + leaves, data = trial1_stats.signs.df.t_pre) #Just one hybrid trial 1
summary(trial1_stats.signs.df.t.lm)

#anova
trial1_stats.signs.df.t.lm.av <- aov(trial1_stats.signs.df.t.lm)
summary(trial1_stats.signs.df.t.lm.av)

#residuals
trial1_stats.signs.df.t.lm.av.residuals <- residuals(object = trial1_stats.signs.df.t.lm.av)
hist( x = trial1_stats.signs.df.t.lm.av.residuals)
ks.test(trial1_stats.signs.df.t.lm.av.residuals, "pnorm", mean(trial1_stats.signs.df.t.lm.av.residuals), sd(trial1_stats.signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#No groups since there are no statistical differences
#=========================plot trial 1==========================================
#total infected leaves
#  count(hybrid)
trial1_stats %>%
  group_by(hybrid, place) %>%
  summarise(top = max(stromata), n=n()) #top is to know where to put the N= in the graph

subset_tria1 <- data.frame(trial1_stats$hybrid, trial1_stats$place, trial1_stats$stromata)

#calculate sample size
ndata_t1 <- subset_tria1 %>%
  group_by(trial1_stats.hybrid, trial1_stats.place)%>%
  summarise(stroma = sum(trial1_stats.stromata)); ndata_t1

#function for 
give.n_t1 <- function(x){
  return(c(y = 40, label = length(x)))  #y is to change where to put the number in the graph
}

pos1 <- position_dodge(0.85) #for annotations N= 
#level_order <- c('H1', 'H2', 'H3', 'H4')
pdf.options(family="Helvetica")
pdf("experiment1.pdf", pointsize=30, width =5)
trial1 <- ggplot(trial1_stats,
                 aes(x = factor(hybrid, level =level_order), y = stromata, fill=place)) +
  geom_boxplot(lwd=1,
               position=position_dodge(width=.85), outlier.size=3) +
  scale_fill_manual(values=c('#f7f7f7', '#969696'),
                    name="Condition",
                    labels=c("Greenhouse", "Growth Chamber")) + 
  scale_y_continuous(breaks = seq(0, 40, 10), limits = c(0, 40)) +
  stat_summary(fun.data = give.n_t1, geom = "text", label = paste("N=", ndata_t1$stroma, sep = ""), size = 8,
               position = pos1) +
  labs(x = "", y = "Number of Stromata/Leaf") + 
  #geom_violin(alpha=0.3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=1),
    axis.ticks.y = element_line(linewidth=1),
    axis.ticks.x = element_line(linewidth=1),
    legend.title = element_text(size=30),
    legend.text = element_text(size=30),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  );trial1
dev.off()
#=========================stats trial 2=========================================
#stats trial2_stats
trial2_stats <- signs[which(signs$essay %in% c('Experiment 2')),]
plotNormalHistogram(trial2_stats$stromata)
shapiro.test(trial2_stats$stromata)

#combine hybrid with condition
trial2_stats$hybridCondition <- paste(trial2_stats$hybrid, 
                                      trial2_stats$place); head(trial2_stats)

#regular stats
trial2_stats.signs.t <- transformTukey(trial2_stats$stromata)
trial2_stats.signs.df.t <- cbind(trial2_stats, trial2_stats.signs.t)
head(trial2_stats.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
trial2_stats.signs.df.t.lm <- lm(trial2_stats.signs.t ~ leaves + hybridCondition, data = trial2_stats.signs.df.t) #GC and GH
summary(trial2_stats.signs.df.t.lm)

#anova
trial2_stats.signs.df.t.lm.av <- aov(trial2_stats.signs.df.t.lm)
summary(trial2_stats.signs.df.t.lm.av)

#residuals
trial2_stats.signs.df.t.lm.av.residuals <- residuals(object = trial2_stats.signs.df.t.lm.av)
hist( x = trial2_stats.signs.df.t.lm.av.residuals)
ks.test(trial2_stats.signs.df.t.lm.av.residuals, "pnorm", mean(trial2_stats.signs.df.t.lm.av.residuals), sd(trial2_stats.signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.
#the warning means that, that different data should not have equal values. 

#groups
#Tukey to check the differences
trial2_stats.signs.posthoc <- TukeyHSD(trial2_stats.signs.df.t.lm.av)
trial2_stats.signs.posthoc

#getting significance letters hybrid
trial2_stats.signs.significance.letters <- multcompLetters4(trial2_stats.signs.df.t.lm.av, trial2_stats.signs.posthoc)
trial2_stats.signs.significance.letters <- as.data.frame.list(trial2_stats.signs.significance.letters$hybrid) #add significance.letters to column | use this to make the plot
trial2_stats.signs.significance.letters
#=========================plot trial 2==========================================
#trial2_stats %>% 
#  count(hybrid)
trial2_stats %>%
  group_by(hybrid, place) %>%
  summarise(top = max(stromata), n=n()) #top is to know where to put the N= in the graph

subset_tria2 <- data.frame(trial2_stats$hybrid, trial2_stats$place, trial2_stats$stromata)

#calculate sample size
nData_t2 <- subset_tria2 %>%
  group_by(trial2_stats.hybrid, trial2_stats.place)%>%
  summarise(stroma = sum(trial2_stats.stromata)); nData_t2

#organize
#is bad with R, so
org_t2 <- c(3,3, 2,2, 1,1, 4)
n_t2 <- cbind(nData_t2, org_t2); n_t2
ndata_t2 <- with(n_t2, n_t2[order(org_t2), ]); ndata_t2

#function for 
give.n_t2 <- function(x){
  return(c(y = 98, label = length(x)))  #y is to change where to put the number in the graph
}

pos <- position_dodge(0.86) #for annotations N= 

pdf.options(family="Helvetica")
pdf("experiment2.pdf", pointsize=30, width=10.5) #width=10
trial2 <- ggplot(trial2_stats,
                 aes(x = factor(hybrid, level=level_order), y = stromata, fill=place)) +
  geom_boxplot(lwd=1,
               position=position_dodge(width=.85), outlier.size=3) +
  scale_fill_manual(values=c('#f7f7f7', '#969696'),
                    name="Condition",
                    labels=c("Greenhouse", "Growth Chamber")) + 
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  stat_summary(fun.data = give.n_t2, geom = "text", label = paste("N=", ndata_t2$stroma, sep = ""), angle = 0, size = 8,
               position = pos) +
  labs(x = "Hybrid", y = "Number of Stromata/Leaf") + 
  #geom_violin(alpha=0.3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=1),
    axis.ticks.y = element_line(linewidth=1),
    axis.ticks.x = element_line(linewidth=1),
    legend.title = element_text(size=30),
    legend.text = element_text(size=30),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  ); trial2 #+ coord_flip() 
dev.off()
#===================================hybrids trial 2=============================
level_order <- c('H1', 'H2', 'H3', 'H4')
trial2_hibrids <- ggplot(trial2_stats,
                 aes(x = factor(hybrid, level =level_order), y = stromata, fill=hybrid)) +
  scale_fill_manual(labels=c("n = 12", "n = 25", "n = 61", "n = 19"), 
                    values=c('#f7f7f7', '#969696', '#cccccc', '#525252'),
                    name="Leaves",
                    breaks=level_order) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  geom_boxplot(outlier.shape=NA) + 
  labs(x = "Hybrid", y = "Number of Stromata/Leaf") +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Roboto Mono", size = 12,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  )

T2 <- trial2_hibrids +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(size = 2, alpha = 0.5, width = 0.1) +
  #facet_wrap(~ essay) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 

#==============================stats for trial 3================================
#stats trial3_stats
trial3_stats <- signs[which(signs$essay %in% c('Experiment 3')),]
plotNormalHistogram(trial3_stats$stromata)
shapiro.test(trial3_stats$stromata)

#regular stats
trial3_stats.signs.t <- transformTukey(trial3_stats$stromata)
trial3_stats.signs.df.t <- cbind(trial3_stats, trial3_stats.signs.t)
head(trial3_stats.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
trial3_stats.signs.df.t.lm <- lm(trial3_stats.signs.t ~ hybrid + leaves, data = trial3_stats.signs.df.t) #no GC
summary(trial3_stats.signs.df.t.lm)

#anova
trial3_stats.signs.df.t.lm.av <- aov(trial3_stats.signs.df.t.lm)
summary(trial3_stats.signs.df.t.lm.av)

#residuals
trial3_stats.signs.df.t.lm.av.residuals <- residuals(object = trial3_stats.signs.df.t.lm.av)
hist( x = trial3_stats.signs.df.t.lm.av.residuals)
ks.test(trial3_stats.signs.df.t.lm.av.residuals, "pnorm", mean(trial3_stats.signs.df.t.lm.av.residuals), sd(trial3_stats.signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#groups
#Tukey to check the differences
trial3_stats.signs.posthoc <- TukeyHSD(trial3_stats.signs.df.t.lm.av)
trial3_stats.signs.posthoc

#getting significance letters hybrid
trial3_stats.signs.significance.letters <- multcompLetters4(trial3_stats.signs.df.t.lm.av, trial3_stats.signs.posthoc)
trial3_stats.signs.significance.letters <- as.data.frame.list(trial3_stats.signs.significance.letters$hybrid) #add significance.letters to column | use this to make the plot
trial3_stats.signs.significance.letters
#=============================plot trial 3======================================
#hybrids
#calculate sample size
trial3_stats %>%
  group_by(hybrid) %>%
  summarise(top = max(stromata), n=n()) #top is to know where to put the N= in the graph

#count number of stromata
#subset_tria3 <-trial3_stats %>% 
#  select(hybrid, stromata); subset_tria3 

subset_tria3 <- data.frame(trial3_stats$hybrid, trial3_stats$stromata)

#calculate sample size
nData_t3 <- subset_tria3 %>%
  group_by(trial3_stats.hybrid)%>%
  summarise(stroma = sum(trial3_stats.stromata)); nData_t3

#organize
#is bad with R, so
org_t3 <- c(3, 2, 1, 4)
n_t3 <- cbind(nData_t3, org_t3); n_t3
ndata_t3 <- with(n_t3, n_t3[order(org_t3), ]); ndata_t3

#function for 
give.n_t3 <- function(x){
  return(c(y = 60, label = length(x)))  #y is to change where to put the number in the graph
}

#level_order <- c('H1', 'H2', 'H3', 'H4')
pdf.options(family="Helvetica")
pdf("experiment3.pdf", pointsize=30)
trial3_hibrids <- ggplot(trial3_stats,
                         aes(x = factor(hybrid, level =level_order), y = stromata, fill=place)) +
  geom_boxplot(lwd=1, outlier.size=3) +
  scale_fill_manual(values=c('#f7f7f7'),
                    name="Condition",
                    labels=c("Greenhouse", "Growth Chamber")) + 
  scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0, 60)) +
  stat_summary(fun.data = give.n_t3, geom = "text", label = paste("N=", ndata_t3$stroma, sep = ""), size = 8) +
  labs(x = "", y = "Number of Stromata/Leaf") + 
  #geom_violin(alpha=0.3)  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=1),
    axis.ticks.y = element_line(linewidth=1),
    axis.ticks.x = element_line(linewidth=1),
    legend.title = element_text(size=30),
    legend.text = element_text(size=30),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  ); trial3_hibrids
dev.off()
#==============================stats for trial 4================================
#stats trial4_stats
trial4_stats <- signs[which(signs$essay %in% c('Experiment 4')),]
str(trial4_stats)
plotNormalHistogram(trial4_stats$stromata)
shapiro.test(trial4_stats$stromata)

#regular stats
trial4_stats.signs.t <- transformTukey(trial4_stats$stromata)
trial4_stats.signs.df.t <- cbind(trial4_stats, trial4_stats.signs.t)
head(trial4_stats.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
trial4_stats.signs.df.t.lm <- lm(trial4_stats.signs.t ~ hybrid, data = trial4_stats.signs.df.t) #no GC
summary(trial4_stats.signs.df.t.lm)

#anova
trial4_stats.signs.df.t.lm.av <- aov(trial4_stats.signs.df.t.lm)
summary(trial4_stats.signs.df.t.lm.av)
#no statistical differences were observed

#=============================plot trial 4======================================
#hybrids
trial4_stats %>%
  group_by(hybrid) %>%
  summarise(top = max(stromata), n=n()) #top is to know where to put the N= in the graph

#count number of stromata
#subset_tria4 <-trial4_stats %>% 
#  select(hybrid, stromata); subset_tria4 

subset_tria4 <- data.frame(trial4_stats$hybrid, trial4_stats$stromata)

#calculate sample size
nData_t4 <- subset_tria4 %>%
  group_by(trial4_stats.hybrid)%>%
  summarise(stroma = sum(trial4_stats.stromata)); nData_t4

#organize
#is bad with R, so
org_t4 <- c(3, 2, 1, 4)
n_t4 <- cbind(nData_t4, org_t4); n_t4
ndata_t4 <- with(n_t4, n_t4[order(org_t4), ]); ndata_t4

#function for 
give.n <- function(x){
  return(c(y = 10, label = length(x)))  #y is to change where to put the number in the graph
}

#level_order <- c('H1', 'H2', 'H3', 'H4') #no need to do it if this was done before
pdf.options(family="Helvetica")
pdf("experiment4.pdf", pointsize=30)
trial4_hibrids <- ggplot(trial4_stats,
                         aes(x = factor(hybrid, level =level_order), y = stromata, fill=place)) +
  geom_boxplot(lwd=1, outlier.size=3) +
  scale_fill_manual(values=c('#f7f7f7'),
                    name="Condition",
                    labels=c("Greenhouse", "Growth Chamber")) + 
  scale_y_continuous(breaks = seq(0, 10, 5), limits = c(0, 10)) +
  labs(x = "", y = "Number of Stromata/Leaf") + 
  stat_summary(fun.data = give.n, geom = "text", label = paste("N=", ndata_t4$strom, sep = ""), size = 8) +
  #geom_violin(alpha=0.3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=1),
    axis.ticks.y = element_line(linewidth=1),
    axis.ticks.x = element_line(linewidth=1),
    legend.title = element_text(size=30),
    legend.text = element_text(size=30),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  ); trial4_hibrids
dev.off()

#==============================stats general====================================
#stats general
#GH
GH <- signs[which(signs$place %in% c('Greenhouse')),]

plotNormalHistogram(GH$stromata)
shapiro.test(GH$stromata)

GH.t <- transformTukey(GH$stromata)
GH.df.t <- cbind(GH, GH.t)
head(GH.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
GH.df.t.lm <- lm(GH.t ~ hybrid + essay + leaves, data = GH.df.t)
summary(GH.df.t.lm)

#anova
GH.df.t.lm.av <- aov(GH.df.t.lm)
summary(GH.df.t.lm.av)

#residuals
GH.df.t.lm.av.residuals <- residuals(object = GH.df.t.lm.av)
hist( x = GH.df.t.lm.av.residuals)
ks.test(GH.df.t.lm.av.residuals, "pnorm", mean(GH.df.t.lm.av.residuals), sd(GH.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#groups
#Tukey to check the differences
GH.posthoc <- TukeyHSD(GH.df.t.lm.av)
GH.posthoc

#trial
GH.significance.letters_trial <- multcompLetters4(GH.df.t.lm.av, GH.posthoc)
GH.significance.letters_trial <- as.data.frame.list(GH.significance.letters_trial$essay) #add significance.letters to column | use this to make the plot
GH.significance.letters_trial

#hybrid
GH.significance.letters_hybrid <- multcompLetters4(GH.df.t.lm.av, GH.posthoc)
GH.significance.letters_hybrid <- as.data.frame.list(GH.significance.letters_hybrid$hybrid) #add significance.letters to column | use this to make the plot
GH.significance.letters_hybrid


#growth chamber
plotNormalHistogram(Gchamber$stromata)
shapiro.test(Gchamber$stromata)

#GH
Gchamber <- signs[which(signs$place %in% c('Growth Chamber')),]

plotNormalHistogram(Gchamber$stromata)
shapiro.test(Gchamber$stromata)

Gchamber.t <- transformTukey(Gchamber$stromata)
Gchamber.df.t <- cbind(Gchamber, Gchamber.t)
head(Gchamber.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
Gchamber.df.t.lm <- lm(Gchamber.t ~ hybrid + essay + leaves, data = Gchamber.df.t)
summary(Gchamber.df.t.lm)

#anova
Gchamber.df.t.lm.av <- aov(Gchamber.df.t.lm)
summary(Gchamber.df.t.lm.av)

#residuals
Gchamber.df.t.lm.av.residuals <- residuals(object = Gchamber.df.t.lm.av)
hist( x = Gchamber.df.t.lm.av.residuals)
ks.test(Gchamber.df.t.lm.av.residuals, "pnorm", mean(Gchamber.df.t.lm.av.residuals), sd(Gchamber.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#general To assess GH vs GC
plotNormalHistogram(signs$stromata)
shapiro.test(signs$stromata)

#combine assay with condition
signs$placeCondition <- paste(signs$essay, 
                              signs$place); head(signs)

#transformation
signs.t <- transformTukey(signs$stromata)
signs.df.t <- cbind(signs, signs.t)
head(signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
signs.df.t.lm <- lm(signs.t ~ placeCondition + leaves, data = signs.df.t)
summary(signs.df.t.lm)

#anova
signs.df.t.lm.av <- aov(signs.df.t.lm)
summary(signs.df.t.lm.av)

#residuals
signs.df.t.lm.av.residuals <- residuals(object = signs.df.t.lm.av)
hist( x = signs.df.t.lm.av.residuals)
ks.test(signs.df.t.lm.av.residuals, "pnorm", mean(signs.df.t.lm.av.residuals), sd(signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#groups
#Tukey to check the differences
signs.posthoc <- TukeyHSD(signs.df.t.lm.av)
signs.posthoc

#trial
signs.significance.letters_trial <- multcompLetters4(signs.df.t.lm.av, signs.posthoc)
signs.significance.letters_trial <- as.data.frame.list(signs.significance.letters_trial$placeCondition) #add significance.letters to column | use this to make the plot
signs.significance.letters_trial

#leaves
signs.significance.letters_leaves <- multcompLetters4(signs.df.t.lm.av, signs.posthoc)
signs.significance.letters_leaves <- as.data.frame.list(signs.significance.letters_leaves$leaves) #add significance.letters to column | use this to make the plot
signs.significance.letters_leaves

#==========================plot General=========================================
#calculate sample size
str(signs)
#just to catch the top # of stromata
signs %>%
  group_by(essay) %>%
  summarise(top = max(stromata), n=n()) #top is to know where to put the N= in the graph


subset_signs <- data.frame(signs$essay, signs$stromata)

#calculate sample size
nData_signs <- subset_signs %>%
  group_by(signs.essay)%>%
  summarise(stroma = sum(signs.stromata)); nData_signs

#function for 
give.signs <- function(x){
  return(c(y = 100, label = length(x)))  #y is to change where to put the number in the graph
}

#met.brewer("Kandinsky")
allplot_general <- ggplot(signs, aes(x=essay, y=stromata, fill=essay)) +
  geom_boxplot(lwd=1, outlier.size=3) +
  scale_fill_manual(values=met.brewer('Manet')) + 
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  stat_summary(fun.data = give.signs, geom = "text", label = paste("N =", nData_signs$stroma), size = 8) +
  labs(x = "Hybrid", y = "Number of Stromata/Leaf") + 
  #geom_violin(alpha=0.3) + 
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=1),
    axis.ticks.y = element_line(linewidth=1),
    axis.ticks.x = element_line(linewidth=1),
    legend.title = element_text(size=30),
    legend.text = element_text(size=30),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.text.x = element_text(size = 30,
                               angle = 0, hjust=0.5, color = "black"), #y axis since we are flipping the graph
    panel.grid = element_blank()
  ) + coord_flip() ; allplot_general

#=========================General stats GC============================================
#stats GC
GC <- signs[which(signs$place %in% c('Growth Chamber')),]
plotNormalHistogram(GC$stromata)
shapiro.test(GC$stromata)

GC.signs.t <- transformTukey(GC$stromata)
GC.signs.df.t <- cbind(GC, GC.signs.t)
head(GC.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
GC.signs.df.t.lm <- lm(GC.signs.t ~ hybrid + essay, data = GC.signs.df.t)
summary(GC.signs.df.t.lm)

#anova
GC.signs.df.t.lm.av <- aov(GC.signs.df.t.lm)
summary(GC.signs.df.t.lm.av)

#residuals
GC.signs.df.t.lm.av.residuals <- residuals(object = GC.signs.df.t.lm.av)
hist( x = GC.signs.df.t.lm.av.residuals)
ks.test(GC.signs.df.t.lm.av.residuals, "pnorm", mean(GC.signs.df.t.lm.av.residuals), sd(GC.signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#groups
#no post hoc since there is no anova differences
#=========================General stats GH======================================
#stats GH
GH <- signs[which(signs$place %in% c('Greenhouse')),]

plotNormalHistogram(GH$stromata)
shapiro.test(GH$stromata)

GH.signs.t <- transformTukey(GH$stromata)
GH.signs.df.t <- cbind(GH, GH.signs.t)
head(GH.signs.df.t)

#stats
#to conduct tukey HSD we need to conduct a lm first
GH.signs.df.t.lm <- lm(GH.signs.t ~ hybrid + essay + leaves, data = GH.signs.df.t)
summary(GH.signs.df.t.lm)

#anova
GH.signs.df.t.lm.av <- aov(GH.signs.df.t.lm)
summary(GH.signs.df.t.lm.av)

#residuals
GH.signs.df.t.lm.av.residuals <- residuals(object = GH.signs.df.t.lm.av)
hist( x = GH.signs.df.t.lm.av.residuals)
ks.test(GH.signs.df.t.lm.av.residuals, "pnorm", mean(GH.signs.df.t.lm.av.residuals), sd(GH.signs.df.t.lm.av.residuals))
#if p-value is greater than 0.05, it suggests that the error terms are normally distributed
#Hence each group was sampled from normally distributed population.

#groups
#Tukey to check the differences
GH.signs.posthoc <- TukeyHSD(GH.signs.df.t.lm.av)
GH.signs.posthoc

#getting significance letters hybrid
GH.signs.significance.letters <- multcompLetters4(GH.signs.df.t.lm.av, GH.signs.posthoc)
GH.signs.significance.letters <- as.data.frame.list(GH.signs.significance.letters$hybrid) #add significance.letters to column | use this to make the plot
GH.signs.significance.letters

#trial
GH.signs.significance.letters_trial <- multcompLetters4(GH.signs.df.t.lm.av, GH.signs.posthoc)
GH.signs.significance.letters_trial <- as.data.frame.list(GH.signs.significance.letters_trial$essay) #add significance.letters to column | use this to make the plot
GH.signs.significance.letters_trial
#=======================end
