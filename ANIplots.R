```R
library(ggthemes)
library(ggplot2)
library(RColorBrewer)
library(yarrr)
library(dplyr)
setwd("~/Documents/phages/spades/NEW/08_clusters")
# First load all the information from vContACT quality
cytoclus <- read.csv(file="RAW_cytoscape_clusters.csv",sep=",") 
colnames(cytoclus)
cytonominal <- unique(subset(cytoclus,select = c(VCluster,HostNEW,HostSimple)))
cy2plot <- cytoclus %>% 
  select(VCluster,covMedian,covQCD,AdjPvalue,Cohesiveness,Quality,TopologyConfidenceScore) %>%
  group_by(VCluster) %>%
  summarize(covMedian=median(covMedian),covQCD=median(covQCD),AdjPvalue=mean(AdjPvalue),mean(Cohesiveness),mean(Quality),mean(TopologyConfidenceScore))
cytoplot <- merge(cy2plot,cytonominal,by = "VCluster")
write.csv(cytoplot,file="cytoplot.csv")
colnames(cytoplot)
allcytoFULL <- read.csv(file="cytoplot.csv",sep=",",header = TRUE)
colnames(allcyto)
allcyto$HostNEW
allcyto <- subset(allcytoFULL, MeanANI > 60 & avgAAI > 60)
subset(allcytoFULL, MeanANI == 100 & avgAAI == 100)
colourCount = length(unique(allcyto$HostSimple))
getPalette = colorRampPalette(piratepal(palette = "basel", trans = .2)) 
hist_top <- ggplot(allcyto, aes(MeanANI,avgAAI)) +
  theme_gray() +
  geom_point(aes(size=totCtg,color=HostSimple),alpha=0.55) +
  scale_colour_manual(values = getPalette(n = colourCount )) +
  scale_size(range = c(2,8))
ggExtra::ggMarginal(hist_top, type = "histogram", groupFill=FALSE, groupColour = FALSE, size = 4, xparams = list(bins=15), yparams = list(bins=15))

# Plot histograms
library(ggpubr)
summary(allcyto)
ani <- ggplot(allcyto, aes(x=MeanANI)) +
  geom_histogram(position="identity",alpha=0.2,color="black", fill = "steelblue",binwidth = 2)+
  geom_vline(xintercept=median(allcyto$MeanANI),color="red")
aai <- ggplot(allcyto, aes(x=avgAAI)) +
  geom_histogram(position="identity",alpha=0.2,color="black", fill = "steelblue",binwidth = 5)+
  geom_vline(xintercept=median(allcyto$avgAAI),color="red")
cvM <-ggplot(allcyto, aes(x=covMedian)) +
  geom_histogram(position="identity",alpha=0.2,color="black", fill = "steelblue",binwidth = 500)+
  geom_vline(xintercept=median(allcyto$covMedian),color="red")
tcs <- ggplot(allcyto, aes(x=TopologyConfidenceScore)) +
  geom_histogram(position="identity",alpha=0.2,color="black", fill = "steelblue",binwidth = 0.1)+
  geom_vline(xintercept=median(allcyto$TopologyConfidenceScore),color="red")
qcd <- ggplot(allcyto, aes(x=covQCD)) +
  geom_histogram(position="identity",alpha=0.2,color="black", fill = "steelblue",binwidth = 0.05)+
  geom_vline(xintercept=mean(allcyto$covQCD),color="red")
ggarrange(ani,aai,tcs,qcd, ncol = 2,  nrow = 2) 
```
