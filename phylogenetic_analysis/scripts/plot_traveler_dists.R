# Libs
library(ggplot2)
library(cowplot)

colorset = c('Georgia, USA (traveler)' ='#D070B9' , 'Georgia, USA (others)'= '#8FBCBB', 
  'Louisiana, USA' = 'lightgoldenrod2', 'Colorado, USA' = 'lightgoldenrod2', 'North Carolina, USA' = 'lightgoldenrod2',
  'Mississippi, USA' = 'lightgoldenrod2', 'Costa Rica' = 'lightgoldenrod2', 'Switzerland' = 'lightgoldenrod2', 'Italy' = 'darkolivegreen4')

## -------------##
#### FIGURE 4 ####
## -------------##
#064L
dat_064L_Italy = read.table('data/travel_analysis/GA-EHC-064L_Italy_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_064L_Switzerland = read.table('data/travel_analysis/GA-EHC-064L_Switzerland_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_064L_Georgia= read.table('data/travel_analysis/GA-EHC-064L_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_064L <- unique(rbind(dat_064L_Italy, dat_064L_Switzerland, dat_064L_Georgia))
dat_064L$Date <- as.Date(dat_064L$Date)
dat_064L$Location <- gsub("USA", ", USA", dat_064L$Location)
dat_064L[dat_064L$accession=='GA-EHC-064L', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_064L[dat_064L$accession=='GA-EHC-064L', "Location"])
dat_064L[(dat_064L$accession!='GA-EHC-064L') & 
  (dat_064L$Location=="Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_064L[(dat_064L$accession!='GA-EHC-064L') & 
    (dat_064L$Location=="Georgia, USA"), "Location"])

p1 <- ggplot(dat_064L, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_064L[dat_064L$Distance == 7, colnames(dat_064L)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_064L[dat_064L$Distance == 8, colnames(dat_064L)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_064L[dat_064L$Distance == 9, colnames(dat_064L)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(7,9,1), limits = c(6.5,9.5)) +
  labs( title="Patient 22 (Sample 064L)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
  theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")

#016P
dat_016P_Louisiana <- 
  read.table('data/travel_analysis/GA-EHC-016P_LouisianaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_016P_Georgia <- 
  read.table('data/travel_analysis/GA-EHC-016P_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_016P = unique(rbind(dat_016P_Louisiana, dat_016P_Georgia))
dat_016P$Date <- as.Date(dat_016P$Date)
dat_016P$Location <- gsub("USA", ", USA", dat_016P$Location)


dat_016P[dat_016P$accession=='GA-EHC-016P', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_016P[dat_016P$accession=='GA-EHC-016P', "Location"])
dat_016P[(dat_016P$accession!='GA-EHC-016P') & (dat_016P$Location == "Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_016P[(dat_016P$accession!='GA-EHC-016P') & (dat_016P$Location == "Georgia, USA"), "Location"])


p2 <- ggplot(dat_016P, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_016P[dat_016P$Distance == 7, colnames(dat_016P)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_016P[dat_016P$Distance == 8, colnames(dat_016P)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_016P[dat_016P$Distance == 9, colnames(dat_016P)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(7,9,1), limits = c(6.5,9.5)) +
  labs( title="Patient 2 (Sample 016P)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
 theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")

# combine
figure_4 = plot_grid(p1, p2, labels = c('A', 'B'), label_size = 14)
ggsave("figures/figure_4.pdf", figure_4, width=14, height=4)



## -------------##
#### FIGURE S6 ####
## -------------##
# 083E
dat_083E_CostaRica <- 
  read.table('data/travel_analysis/GA-EHC-083E_Costa Rica_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_083E_Georgia <- 
  read.table('data/travel_analysis/GA-EHC-083E_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_083E = unique(rbind(dat_083E_CostaRica, dat_083E_Georgia))
dat_083E$Date <- as.Date(dat_083E$Date)
dat_083E$Location <- gsub("USA", ", USA", dat_083E$Location)
dat_083E[dat_083E$accession=='GA-EHC-083E', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_083E[dat_083E$accession=='GA-EHC-083E', "Location"])
dat_083E[(dat_083E$accession!='GA-EHC-083E') & (dat_083E$Location == "Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_083E[(dat_083E$accession!='GA-EHC-083E') & (dat_083E$Location == "Georgia, USA"), "Location"])

#083E
p1 <- ggplot(dat_083E, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_083E[dat_083E$Distance == 4, colnames(dat_083E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_083E[dat_083E$Distance == 5, colnames(dat_083E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_083E[dat_083E$Distance == 6, colnames(dat_083E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(4,6,1), limits = c(3.5,6.5)) +
  labs( title="Patient 23 (Sample 083E)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
  theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")

# 043Q
dat_043Q_Mississippi <- 
  read.table('data/travel_analysis/GA-EHC-043Q_MississippiUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_043Q_Georgia <- 
  read.table('data/travel_analysis/GA-EHC-043Q_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_043Q = unique(rbind(dat_043Q_Mississippi, dat_043Q_Georgia))
dat_043Q$Date <- as.Date(dat_043Q$Date)
dat_043Q$Location <- gsub("USA", ", USA", dat_043Q$Location)
dat_043Q[dat_043Q$accession=='GA-EHC-043Q', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_043Q[dat_043Q$accession=='GA-EHC-043Q', "Location"])
dat_043Q[(dat_043Q$accession!='GA-EHC-043Q') & (dat_043Q$Location == "Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_043Q[(dat_043Q$accession!='GA-EHC-043Q') & (dat_043Q$Location == "Georgia, USA"), "Location"])


p2 <- ggplot(dat_043Q, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_043Q[dat_043Q$Distance == 11, colnames(dat_043Q)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_043Q[dat_043Q$Distance == 12, colnames(dat_043Q)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_043Q[dat_043Q$Distance == 13, colnames(dat_043Q)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(11,13,1), limits = c(10.5,13.5)) +
  labs( title="Patient 14 (Sample 043Q)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
  theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")

# 031E
dat_031E_Colorado <- 
  read.table('data/travel_analysis/GA-EHC-031E_ColoradoUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_031E_Georgia <- 
  read.table('data/travel_analysis/GA-EHC-031E_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_031E= unique(rbind(dat_031E_Colorado, dat_031E_Georgia))
dat_031E$Date <- as.Date(dat_031E$Date)
dat_031E$Location <- gsub("USA", ", USA", dat_031E$Location)
dat_031E[dat_031E$accession=='GA-EHC-031E', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_031E[dat_031E$accession=='GA-EHC-031E', "Location"])
dat_031E[(dat_031E$accession!='GA-EHC-031E') & (dat_031E$Location == "Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_031E[(dat_031E$accession!='GA-EHC-031E') & (dat_031E$Location == "Georgia, USA"), "Location"])


p3 <- ggplot(dat_031E, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_031E[dat_031E$Distance == 5, colnames(dat_031E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_031E[dat_031E$Distance == 6, colnames(dat_031E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_031E[dat_031E$Distance == 7, colnames(dat_031E)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(5,7,1), limits = c(4.5,7.5)) +
  labs( title="Patient 5 (Sample 031E)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
  theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")

# finally 110f
dat_110F_NorthCarolina <- 
  read.table('data/travel_analysis/GA-EHC-110F_North CarolinaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_110F_Georgia <- 
  read.table('data/travel_analysis/GA-EHC-110F_GeorgiaUSA_dists.tsv', 
    sep='\t', col.names=c('accession', 'Location', 'Date', 'Distance', 'substitutions'))
dat_110F = unique(rbind(dat_110F_NorthCarolina, dat_110F_Georgia))
dat_110F$Date <- as.Date(dat_110F$Date)
dat_110F$Location <- gsub("USA", ", USA", dat_110F$Location)
dat_110F[dat_110F$accession=='GA-EHC-110F', "Location"] <- 
  gsub(", USA", ", USA (traveler)", dat_110F[dat_110F$accession=='GA-EHC-110F', "Location"])
dat_110F[(dat_110F$accession!='GA-EHC-110F') & (dat_110F$Location == "Georgia, USA"), "Location"] <- 
  gsub(", USA", ", USA (others)", 
    dat_110F[(dat_110F$accession!='GA-EHC-110F') & (dat_110F$Location == "Georgia, USA"), "Location"])

p4 <- ggplot(dat_110F, aes(y=Date, x=Distance, fill=Location)) + 
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_110F[dat_110F$Distance == 7, colnames(dat_110F)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_110F[dat_110F$Distance == 8, colnames(dat_110F)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    data = dat_110F[dat_110F$Distance == 9, colnames(dat_110F)],
    stackratio=1,
    dotsize = 1,
    method = 'histodot',
    color = "#575c66",
    stackgroups = T,
    binwidth = 1) +
  scale_y_date(breaks="5 days", date_labels="%m/%d", limits = as.Date(c('2020-02-24','2020-03-31'))) +
  scale_x_continuous(breaks = seq(7,9,1), limits = c(6.5,9.5)) +
  labs( title="Patient 27 (Sample 110F)",
        x="SNP Distance From Wuhan/Hu-1", y="Date (2020)") +
  theme(legend.position = "right", panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_blank(), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12)) + 
  theme(axis.text=element_text(size=12), 
        axis.title =element_text(size=12, face="bold")) +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_manual(values=colorset) +
  coord_flip(clip="off")


figure_s6 <- plot_grid(p1, p2, p3, p4, 
  labels = c('A', 'B', 'C', 'D'), label_size = 12)

ggsave("figures/figure_s6.pdf", figure_s6, width=14, height=8)

