library(reshape2)
library(ggplot2)
library(ggrepel)
library(lme4)
library(emmeans)
library(sommer)
library(AlphaSimR)
library(qvalue)

#load all the input files.
geno <- readRDS("source/SWB_geno_012.RDS")
line.info <- read.csv("source/SWB_line_info.csv", header=T, as.is=T)
trait <- read.csv("source/DUS_traits_20200430.csv", header=T, as.is=T)
trait.SASA <- read.csv("source/DUS_traits_20200430_SASA.csv", header=T, as.is=T)
trait.names <- read.csv("source/DUS_trait_names.csv", header=F, as.is=T)
trait.dmyld <- read.csv("source/YLD_20200604.csv", header=T, as.is=T)

#remove lines without AFP number.
line.info <- line.info[!is.na(line.info$AFP),]
trait <- trait[!is.na(trait$AFP),]
trait.SASA <- trait.SASA[!is.na(trait.SASA$AFP),]

#merge the line information and trait data by AFP.
trait <- merge(line.info, trait, by="AFP", all=F)
trait <- trait[order(trait$AFP),]

#check and fix any discrepancies.
trait[which(!(trait$Year==trait$year)),c(1,3,5,17)]
#    AFP     Line Year year
#13  721 Corniche 1984 1983
#34 1011  Derkado 1989 1988
#based on AFP, the correct ones are 1984 for Corniche and 1989 for Derkado.

table(trait[,c(7,30)])
#       t12
#RowType   1   2
#     2R 667   0
#     6R   0  39

table(trait[,c(6,46)])
#      t28
#Season   1   2   3
#    AB   0   5   0
#    SB   0   0 370
#    WB 335   0   0

trait <- trait[,c(2,3,1,5:13,19:46)] #keep only neccessary columns.

#set variables as factor in the DMYLD dataset.
trait.dmyld$Line <- as.factor(trait.dmyld$Line)
trait.dmyld$Year <- as.factor(trait.dmyld$Year)
trait.dmyld$Management <- as.factor(trait.dmyld$Management)
trait.dmyld$LocationID <- as.factor(trait.dmyld$LocationID)

#run LMM and calculate the BLUE for dmyld.
dmyld.lmm <- lmer(data=trait.dmyld,
                  DMYLD ~ Line + (1|Management) + (1|Management:Year) + (1|Management:Year:Line) + (1|Management:Year:LocationID), na.action=na.omit)
dmyld.blue <- data.frame(emmeans(dmyld.lmm, ~Line))[,1:2]
colnames(dmyld.blue)[2] <- "BLUE_DMYLD"
temp <- unique(trait.dmyld[,1:2])
dmyld.blue <- merge(dmyld.blue, temp, by="Line", sort=F)
dmyld.blue <- dmyld.blue[,c(1,3,2)]
dmyld.blue$Line <- as.character(dmyld.blue$Line)
write.csv(dmyld.blue, "source/BLUE_DMYLD.csv", quote=F, row.names=F)

#keep only DMYLD for lines with AFP.
trait.dmyld <- dmyld.blue[!is.na(dmyld.blue$AFP),]

#reduce the marker genotype to 710 lines with DUS trait data, and convert from 0/1/2 to -1/0/1.
geno <- geno[trait$genoID,] - 1

#remove markers with missing data.
geno <- geno[,!(colSums(is.na(geno))>0)]

#remove monomorphic markers.
geno <- geno[,!(colSums(geno==-1)==nrow(geno) | colSums(geno==0)==nrow(geno) | colSums(geno==1)==nrow(geno))]

dim(geno)
#710 40065

#FigS01-DUS trait histograms#####
#plot the distribution for all 28 DUS traits.
temp <- trait[,c(5,13:40)]
colnames(temp)[-1] <- paste(trait.names$V1, " - ", trait.names$V2, sep="")
temp <- melt(temp, id.vars="Season")
temp$value <- as.factor(temp$value)
temp$Season <- as.factor(temp$Season)
temp$Season <- factor(temp$Season, levels=c("SB","AB","WB"))  

ggplot() +
  geom_bar(data=temp, aes(value, fill=Season, color=Season), position="stack") +
  facet_wrap(vars(variable), scales="free", drop=F, ncol=4, labeller=label_wrap_gen()) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(strip.background=element_blank(), strip.text=element_text(size=7)) +
  theme(axis.text=element_text(size=7), axis.title=element_text(size=8)) +
  theme(legend.title=element_text(size=7), legend.text=element_text(size=7)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="#505050") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="#505050") +
  scale_fill_manual(values=c("#8080FF","#FF8080","#80FF80")) +
  scale_color_manual(values=c("#8080FF","#FF8080","#80FF80")) +
  xlab("trait score")

ggsave("output/sfig01.png",
       height=8,
       width=6.5,
       units="in",
       dpi=600)


#Tab01/TabS02-DUS trait heritabilities#####
#calculate the GRM for 710 lines.
grm <- A.mat(geno)
grm.SB <- A.mat(geno[trait$Season=="SB",])
grm.WB <- A.mat(geno[trait$Season=="WB",])

#run mixed models for each trait (combined, spring only, winter only).
out.univar.C <- out.univar.SB <- out.univar.WB <- replicate(28,list())
for(i in 1:28){
  temp <- trait[,c(1,4,5,12+i)]
  colnames(temp)[4] <- "t"
  out.univar.C[[i]] <- mmer(t~Year + Season,
                            random=~vs(genoID, Gu=grm),
                            rcov=~units,
                            data=temp)
  
  temp <- trait[trait$Season=="SB",c(1,4,5,12+i)]
  colnames(temp)[4] <- "t"
  if(length(c(table(temp$t))) > 1) {
    out.univar.SB[[i]] <- mmer(t~Year,
                               random=~vs(genoID, Gu=grm.SB),
                               rcov=~units,
                               data=temp)
  }
  
  temp <- trait[trait$Season=="WB",c(1,4,5,12+i)]
  colnames(temp)[4] <- "t"
  if(length(c(table(temp$t))) > 1) {
    out.univar.WB[[i]] <- mmer(t~Year,
                               random=~vs(genoID, Gu=grm.WB),
                               rcov=~units,
                               data=temp)
  }
}

#summarize the variance components and heritabilities for all traits.
out.varcomp <- vector()
for(i in 1:28){
  if(length(out.univar.C[[i]])>0) {
    temp.C <- cbind(pin(out.univar.C[[i]], ~V1),
                    pin(out.univar.C[[i]], ~V1+V2),
                    pin(out.univar.C[[i]], ~V1/(V1+V2)))
    temp.C <- c(t(temp.C))
  } else {
    temp.C <- rep(NA,6)
  }
  
  if(length(out.univar.SB[[i]])>0) {
    temp.SB <- cbind(pin(out.univar.SB[[i]], ~V1),
                     pin(out.univar.SB[[i]], ~V1+V2),
                     pin(out.univar.SB[[i]], ~V1/(V1+V2)))
    temp.SB <- c(t(temp.SB))
  } else {
    temp.SB <- rep(NA,6)
  }
  
  if(length(out.univar.WB[[i]])>0) {
    temp.WB <- cbind(pin(out.univar.WB[[i]], ~V1),
                     pin(out.univar.WB[[i]], ~V1+V2),
                     pin(out.univar.WB[[i]], ~V1/(V1+V2)))
    temp.WB <- c(t(temp.WB))
  } else {
    temp.WB <- rep(NA,6)
  }
  
  out.varcomp <- rbind(out.varcomp, c(temp.C, temp.SB, temp.WB))
}

#clean up the univariate outputs.
rownames(out.varcomp) <- colnames(trait)[13:40]
colnames(out.varcomp) <- c("VgC","VgC_SE","VpC","VpC_SE","H2C","H2C_SE",
                           "VgS","VgS_SE","VpS","VpS_SE","H2S","H2S_SE",
                           "VgW","VgW_SE","VpW","VpW_SE","H2W","H2W_SE")

#export the univariate output for making Table 1.
write.csv(out.varcomp, "temp/univariate_varcomp.csv", quote=F, row.names=T)



#Fig01-DUS trait score difference#####
#subset the trait data for which they are in common between NIAB and SASA.
temp.NIAB <- trait[trait$AFP%in%trait.SASA$AFP,]
rownames(temp.NIAB) <- temp.NIAB$AFP
temp.NIAB <- as.matrix(temp.NIAB[,13:40])
temp.SASA <- as.matrix(trait.SASA[,6:33])
rownames(temp.SASA) <- trait.SASA$AFP
temp.SASA <- temp.SASA[rownames(temp.NIAB),]

#since t03, t23 and t26 are binary 1 or 9, we need to change it to binary 1 or 2.
temp.NIAB[which(temp.NIAB[,3]==9), 3] <- 2
temp.NIAB[which(temp.NIAB[,23]==9), 23] <- 2
temp.NIAB[which(temp.NIAB[,26]==9), 26] <- 2
temp.SASA[which(temp.SASA[,3]==9), 3] <- 2
temp.SASA[which(temp.SASA[,23]==9), 23] <- 2
temp.SASA[which(temp.SASA[,26]==9), 26] <- 2

#calculate the trait score differences between NIAB and SASA datasets.
trait.diff <- abs(temp.NIAB - temp.SASA)

#summarize the trait score differences.
diff.byLine <- sapply(0:8, FUN=function(x) rowSums(trait.diff==x, na.rm=T))
diff.byLine <- diff.byLine/rowSums(diff.byLine)
colnames(diff.byLine) <- 0:8

#Fig01A - boxplot of trait score differences.
dat.plot <- data.frame(AFP=rownames(diff.byLine), diff.byLine)
colnames(dat.plot)[-1] <- 0:8
dat.plot <- melt(dat.plot, id.vars="AFP")
dat.plot$AFP <- factor(dat.plot$AFP, levels=unique(dat.plot$AFP))

ggplot() +
  geom_boxplot(data=dat.plot, aes(x=variable, y=value), outlier.size=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  xlab("DUS trait score difference") +
  ylab("Proportion")

ggsave("output/Fig01A.svg",
       height=2.5,
       width=3.5,
       units="in",
       scale=4/3)
ggsave("output/Fig01A.png",
       height=2.5,
       width=3.5,
       units="in",
       dpi=600)

#Fig01B - barplot of trait score differences over time.
ggplot() +
  geom_bar(data=dat.plot, aes(x=AFP, y=value, fill=variable), position="stack", stat="identity", width=1, color=NA) +
  scale_fill_manual(values=c("#CCCCCC", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"), name="DUS trait score difference") +
  ylab("Proportion") +
  xlab("Line") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  theme(legend.position="bottom", legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  guides(fill=guide_legend(title.position="top", label.position="bottom", nrow=1)) +
  theme(legend.title=element_text(size=9), legend.text=element_text(size=8)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,1.5))

ggsave("output/Fig01B.svg",
       height=3,
       width=3.5,
       units="in",
       scale=4/3)
ggsave("output/Fig01B.png",
       height=3,
       width=3.5,
       units="in",
       dpi=600)

#Fig01C/D - Rolling mean of Manhattan distances of DUS traits over time.
trait.SB <- trait[trait$Season=="SB",]
rolldist.SB <- sapply(10:360, FUN=function(x) mean(c(dist(x=trait.SB[(x-9):(x+10), 13:40], method="manhattan"))/28))
plot(1:length(rolldist.SB), rolldist.SB)

trait.WB <- trait[trait$Season=="WB",]
rolldist.WB <- sapply(10:325, FUN=function(x) mean(c(dist(x=trait.WB[(x-9):(x+10), 13:40], method="manhattan"))/28))
plot(1:length(rolldist.WB), rolldist.WB)

trait.SB <- trait[trait$Season=="SB",]
trait.SB <- trait.SB[!(rowSums(is.na(trait.SB[,13:39])) > 4), ]
rolldist.SB <- sapply(10:323, FUN=function(x) mean(c(dist(x=trait.SB[(x-9):(x+10), 13:40], method="manhattan"))/27))
rolldist.SB <- data.frame(dist=rolldist.SB, Time=1:length(rolldist.SB))

ggplot() +
  geom_point(data=rolldist.SB, aes(x=Time, y=dist), size=1) +
  scale_y_continuous(name="Manhattan distance") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999")
ggsave("output/Fig01C.svg",
       height=2,
       width=3.5,
       units="in",
       scale=4/3)
ggsave("output/Fig01C.png",
       height=2,
       width=3.5,
       units="in",
       dpi=600)

trait.WB <- trait[trait$Season=="WB",]
trait.WB <- trait.WB[!(rowSums(is.na(trait.WB[,13:39])) > 4), ]
rolldist.WB <- sapply(10:271, FUN=function(x) mean(c(dist(x=trait.WB[(x-9):(x+10), 13:40], method="manhattan"))/27))
rolldist.WB <- data.frame(dist=rolldist.WB, Time=1:length(rolldist.WB))
ggplot() +
  geom_point(data=rolldist.WB, aes(x=Time, y=dist), size=1) +
  scale_y_continuous(name="Manhattan distance") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999")
ggsave("output/Fig01D.svg",
       height=2,
       width=3.5,
       units="in",
       scale=4/3)
ggsave("output/Fig01D.png",
       height=2,
       width=3.5,
       units="in",
       dpi=600)

#get some information about the lines.
trait.info <- trait[,c(1,3,5)]
rownames(trait.info) <- trait$AFP
trait.info <- trait.info[rownames(trait.diff),]

#import the previously calculated heritabilities.
out.varcomp <- read.csv("temp/univariate_varcomp.csv", header=T, row.names=1, as.is=T)

#Fig01E - scatter plot of proportion trait score difference vs heritability of each trait.
dat.plot <- data.frame(diff=c(colSums(trait.diff>0, na.rm=T)/colSums(trait.diff>=0, na.rm=T),
                              colSums(trait.diff[trait.info$Season=="SB",]>0, na.rm=T)/colSums(trait.diff[trait.info$Season=="SB",]>=0, na.rm=T),
                              colSums(trait.diff[trait.info$Season=="WB",]>0, na.rm=T)/colSums(trait.diff[trait.info$Season=="WB",]>=0, na.rm=T)),
                       h2=c(out.varcomp$H2C, out.varcomp$H2S, out.varcomp$H2W),
                       pop=c(rep("Combined",28), rep("Spring",28), rep("Winter",28)),
                       trait=rep(1:28,3)) #colnames(trait.diff)
dat.plot <- dat.plot[!(is.na(dat.plot$diff) | is.na(dat.plot$h2)),]

cor(dat.plot$diff[dat.plot$pop=="Combined"], dat.plot$h2[dat.plot$pop=="Combined"]) #-0.67
cor(dat.plot$diff[dat.plot$pop=="Spring"], dat.plot$h2[dat.plot$pop=="Spring"]) #-0.61
cor(dat.plot$diff[dat.plot$pop=="Winter"], dat.plot$h2[dat.plot$pop=="Winter"]) #-0.59

temp <- data.frame(pop=c("Combined","Spring","Winter"),
                   cor=c("r = -0.67", "r = -0.61", "r = -0.59"))

ggplot() +
  geom_text_repel(data=dat.plot, aes(x=h2, y=diff, label=trait), box.padding=unit(0.1,"lines"), point.padding=unit(0,"lines"), size=3, min.segment.length=unit(0,"lines"), segment.size=0.2) +
  geom_text(data=temp, aes(x=1, y=1, label=cor), size=3) +
  facet_wrap(vars(pop), ncol=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_y_continuous(limits=c(0,1), expand=c(0,0.1)) +
  scale_x_continuous(expand=c(0,0.1)) +
  ylab("Proportion") +
  xlab(bquote("h"^2))

ggsave("output/Fig01C.svg",
       height=6,
       width=3.5,
       units="in",
       scale=4/3)
ggsave("output/Fig01C.png",
       height=6,
       width=3.5,
       units="in",
       dpi=600)
       
#SFig02-Histograms of trait difference.
trait.diff <- temp.NIAB - temp.SASA
temp <- which(colSums(!is.na(trait.diff)) > 0)
temp.mean <- t(sapply(temp, FUN=function(x) c(x, round(mean(trait.diff[,x], na.rm=T), 2))))
temp.mean <- data.frame(temp.mean)
colnames(temp.mean) <- c("trait", "mean")

dat.plot <- vector()
for(i in temp){
  dat.plot <- rbind(dat.plot,
                    data.frame(x=names(c(table(trait.diff[,i]))),
                               y=c(table(trait.diff[,i])),
                               trait=i))
}
dat.plot$x <- as.integer(dat.plot$x)

ggplot() +
  geom_segment(data=temp.mean, aes(x=mean, xend=mean, y=-Inf, yend=Inf), color="#FF0000") +
  geom_point(data=dat.plot, aes(x=x, y=y)) +
  facet_wrap(vars(trait), ncol=5, scales="free") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
  xlab("trait difference") +
  ylab("count")

ggsave("output/SFig02.svg",
       height=6,
       width=7,
       units="in",
       scale=4/3)



#Fig02-Correlations between DUS traits and yield#####
#merge the trait information between DUS and DMYLD by AFP.
trait2 <- merge(trait, trait.dmyld[,c(2,3)], by="AFP", all=F, sort=F) #192 spring barley lines.

#get the marker data for these 192 lines and remove monomorphic markers.
geno2 <- geno[trait2$genoID,]
geno2 <- geno2[,!(colSums(geno2==-1)==nrow(geno2) | colSums(geno2==0)==nrow(geno2) | colSums(geno2==1)==nrow(geno2))]

#calculate the grm for these 192 lines.
grm2 <- A.mat(geno2)

#run bivariate DUS-DMYLD mmer with year as fixed effect.
out.multivar <- replicate(28, list())
for(i in 1:28){
  temp <- trait2[,c(2,4,12+i,41)]
  colnames(temp)[3:4] <- c("x","y")
  if(length(c(table(temp$x))) > 1){
    out.multivar[[i]] <- mmer(cbind(x,y)~Year,
                              random=~vs(genoID, Gu=grm2),
                              rcov=~units,
                              data=temp)
  }
}

#save the multivariate runs.
saveRDS(out.multivar, "temp/out_multivar.RDS")

#identify which mmer runs went OK.
temp <- which(sapply(1:28, FUN=function(x) if(length(out.multivar[[x]]) > 0) nrow(summary(out.multivar[[x]])$varcomp) else 0)==6)

#collect the bivariate mmer results.
out.covar <- matrix(NA, nrow=28, ncol=22)
for(i in temp){
  out.covar[i,] <- c(c(t(pin(out.multivar[[i]], ~V1))),
                     c(t(pin(out.multivar[[i]], ~V2))),
                     c(t(pin(out.multivar[[i]], ~V3))),
                     c(t(pin(out.multivar[[i]], ~V1+V4))),
                     c(t(pin(out.multivar[[i]], ~V2+V5))),
                     c(t(pin(out.multivar[[i]], ~V3+V6))),
                     c(t(pin(out.multivar[[i]], ~V1/(V1+V4)))),
                     c(t(pin(out.multivar[[i]], ~V3/(V3+V6)))),
                     c(t(pin(out.multivar[[i]], ~V2/sqrt(V1*V3)))),
                     c(t(pin(out.multivar[[i]], ~(V2+V5)/sqrt((V1+V4)*(V3+V6))))),
                     c(t(pin(out.multivar[[i]], ~V5/sqrt(V4*V6)))))
}

colnames(out.covar) <- c("VgX", "VgX_SE", "Covg", "Covg_SE", "VgY", "VgY_SE",
                         "VpX", "VpX_SE", "Covp", "Covp_SE", "VpY", "VpY_SE",
                         "H2X", "H2X_SE", "H2Y", "H2Y_SE",
                         "Corg", "Corg_SE", "Corp", "Corp_SE", "CorE", "CorE_SE")
rownames(out.covar) <- colnames(trait2)[13:40]

write.csv(out.covar, "temp/multivariate_varcomp.csv", quote=F, row.names=T)

#plot the correlations.
dat.plot <- data.frame(trait=rep(1:28,2), rbind(out.covar[,17:18], out.covar[,19:20]), Type=c(rep("Genetic",28),rep("Phenotypic",28)))
colnames(dat.plot)[2:3] <- c("value", "SE")

ggplot(data=dat.plot, aes(x=trait, y=value, shape=Type, color=Type)) +
  annotate("rect", xmin=seq(0.5,26.5,2), xmax=seq(1.5,28.5,2), ymin=-Inf, ymax=Inf, color=NA, fill="#EEEEEE", linetype=2) +
  annotate("segment", x=-Inf, xend=Inf, y=0, yend=0, color="#808080", linetype=2, size=0.3) +
  annotate("segment", x=-Inf, xend=Inf, y=-1, yend=-1, color="#808080", linetype=2, size=0.3) +
  annotate("segment", x=-Inf, xend=Inf, y=1, yend=1, color="#808080", linetype=2, size=0.3) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), position=position_dodge(0.8), width=0.3, color="#505050", size=0.3) +
  geom_point(position=position_dodge(width=0.8), na.rm=T) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_x_continuous(breaks=1:28, limits=c(0.5,28.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1.6,1.6)) +
  scale_shape_manual(values=c(16,15)) +
  theme(legend.position="bottom", legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  theme(legend.title=element_text(size=9), legend.text=element_text(size=8)) +
  ylab("Correlation") +
  xlab("Trait")

ggsave("output/Fig02.png",
       height=3,
       width=7,
       units="in",
       dpi=600)
ggsave("output/Fig02.svg",
       height=3,
       width=7,
       units="in",
       scale=4/3)


#Fig03/TabS06-DUS markers simulation#####
#load the original marker data for 805 lines and DUS marker information.
geno.DUS <- readRDS("source/SWB_geno_012.RDS")
marker.DUS <- read.csv("source/DUS_markers.csv", header=T, as.is=T)

#reduce the marker data to only 39 DUS markers.
geno.DUS <- geno.DUS[,marker.DUS$NewMarkerName]

#identify lines that have marker data for both parents.
lineF <- line.info[which(line.info$CrossType==1),c(1,12:13)]
lineF <- lineF[!(rowSums(is.na(lineF))>0),]
lineF <- as.matrix(lineF)

#loop to simulate progeny from each parent pair.
out.U <- vector()
out.NS <- vector()
out.S <- vector()
out.B <- vector()
out.sim.DH <- list()
out.real.DH <- vector()
out.sim.F6 <- list()
out.real.F6 <- vector()
for(i in 1:nrow(lineF)){
  #extract the marker data.
  temp.geno <- geno.DUS[lineF[i,],]
  temp.marker <- marker.DUS
  rownames(temp.marker) <- temp.marker$NewMarkerName
  
  #count the number of unknown alleles in the progeny.
  out.U <- c(out.U,
             sum((temp.geno[1,]==2 & temp.geno[2,]==0 & temp.geno[3,]==0) | (temp.geno[1,]==0 & temp.geno[2,]==2 & temp.geno[3,]==2), na.rm=T))
  
  #count the number of non-segregating markers.
  out.NS <- c(out.NS,
              sum((temp.geno[1,]==0 & temp.geno[2,]==0 & temp.geno[3,]==0) | (temp.geno[1,]==2 & temp.geno[2,]==2 & temp.geno[3,]==2), na.rm=T))
  
  #count the number of segregating markers.
  out.S <- c(out.S,
             sum((temp.geno[2,]==0 & temp.geno[3,]==2) | (temp.geno[2,]==2 & temp.geno[3,]==0), na.rm=T))
  
  #count the number of bad markers.
  out.B <- c(out.B,
             ncol(temp.geno) - out.U[i] - out.NS[i] - out.S[i])

  #keep only segregating markers for simulation.
  temp.geno <- temp.geno[,which((temp.geno[2,]==0 & temp.geno[3,]==2) | (temp.geno[2,]==2 & temp.geno[3,]==0))]
  temp.marker <- temp.marker[colnames(temp.geno),]

  #separate the marker data by chromosome.
  parentHap <- lapply(unique(temp.marker$Chr), FUN=function(x) temp.geno[-1, temp.marker$Chr==x, drop=F])
  genMap <- lapply(unique(temp.marker$Chr), FUN=function(x) temp.marker$GenPos[temp.marker$Chr==x]/100)
  genMap <- lapply(1:length(genMap), FUN=function(x) genMap[[x]] - genMap[[x]][1])

  #transpose the marker data for easier subsequent uses.
  temp.geno <- t(temp.geno)
    
  #initiate the founder haplotype based on the two parents.
  founder <- newMapPop(genMap=genMap,
                       haplotypes=parentHap,
                       inbred=T,
                       ploidy=2L)
  SP <- SimParam$new(founder)

  #create the founder population.
  founder <- newPop(founder, simParam=SP)
  
  #make a cross between the two parents.
  F1 <- randCross(pop=founder,
                  nCrosses=1,
                  nProgeny=1,
                  simParam=SP)
  
  #simulate 10,000 double haploids from the cross.
  DH <- makeDH(pop=F1,
               nDH=10000,
               simParam=SP)
  
  #extract the marker genotypes for the simulated DHs.
  geno.DH <- pullSegSiteGeno(pop=DH,
                             chr=NULL,
                             simParam=SP)
  
  
  #tabulate the occurences of unique haplotypes in the simulated DHs.
  unique.DH <- unique(geno.DH)
  geno.DH <- t(geno.DH)
  match.sim.DH <- sapply(1:nrow(unique.DH), FUN=function(x) sum(colSums(geno.DH==unique.DH[x,])==nrow(geno.DH)))*100/ncol(geno.DH)
  out.sim.DH <- c(out.sim.DH, list(match.sim.DH))
  
  #compare the simulated DHs to the progeny and two parents.
  match.real.DH <- cbind((colSums(geno.DH==temp.geno[,1], na.rm=T)==nrow(geno.DH)-sum(is.na(temp.geno[,1]))),
                         (colSums(geno.DH==temp.geno[,2], na.rm=T)==nrow(geno.DH)),
                         (colSums(geno.DH==temp.geno[,3], na.rm=T)==nrow(geno.DH)))
  
  #bootstrap for 1,000 cycles to get the mean and variance for match.real.DH.
  bs.DH <- vector()
  for(j in 1:1000){
    bs.idx <- sample(1:nrow(match.real.DH), nrow(match.real.DH), replace=T)
    bs.DH <- rbind(bs.DH, colSums(match.real.DH[bs.idx,]))
  }
  bs.DH <- bs.DH*100/nrow(match.real.DH)
  out.real.DH <- rbind(out.real.DH,
                       c(sapply(1:3, FUN=function(x) mean(bs.DH[,x])),
                         sapply(1:3, FUN=function(x) var(bs.DH[,x]))))
  
  #simulate 10,000 F6 progeny from the cross as an alternative to DH.
  #self the F1 and create 10,000 progeny.
  F2 <- self(pop=F1,
             nProgeny=10000,
             simParam=SP)
  
  #self the F2 and keep 1 progeny each.
  F3 <- self(pop=F2,
             nProgeny=1,
             simParam=SP)
  
  #self the F3 and keep 1 progeny each.
  F4 <- self(pop=F3,
             nProgeny=1,
             simParam=SP)
  
  #self the F4 and keep 1 progeny each.
  F5 <- self(pop=F4,
             nProgeny=1,
             simParam=SP)
  
  #self the F5 and keep 1 progeny each.
  F6 <- self(pop=F5,
             nProgeny=1,
             simParam=SP)
  
  #extract the marker genotypes for the simulated progeny.
  geno.F6 <- pullSegSiteGeno(pop=F6,
                             chr=NULL,
                             simParam=SP)
  
  #tabulate the occurences of unique haplotypes in the simulated F6s.
  unique.F6 <- unique(geno.F6)
  geno.F6 <- t(geno.F6)
  match.sim.F6 <- sapply(1:nrow(unique.F6), FUN=function(x) sum(colSums(geno.F6==unique.F6[x,])==nrow(geno.F6)))*100/ncol(geno.F6)
  out.sim.F6 <- c(out.sim.F6, list(match.sim.F6))
  
  #compare the simulated DHs to the progeny and two parents.
  match.real.F6 <- cbind((colSums(geno.F6==temp.geno[,1], na.rm=T)==nrow(geno.F6)-sum(is.na(temp.geno[,1]))),
                         (colSums(geno.F6==temp.geno[,2], na.rm=T)==nrow(geno.F6)),
                         (colSums(geno.F6==temp.geno[,3], na.rm=T)==nrow(geno.F6)))
  
  #bootstrap for 1,000 cycles to get the mean and variance for match.real.F6.
  bs.F6 <- vector()
  for(j in 1:1000){
    bs.idx <- sample(1:nrow(match.real.F6), nrow(match.real.F6), replace=T)
    bs.F6 <- rbind(bs.F6, colSums(match.real.F6[bs.idx,]))
  }
  bs.F6 <- bs.F6*100/nrow(match.real.F6)
  out.real.F6 <- rbind(out.real.F6,
                       c(sapply(1:3, FUN=function(x) mean(bs.F6[,x])),
                         sapply(1:3, FUN=function(x) var(bs.F6[,x]))))
  
  message(i)
  
}

save.image("temp/temp_20200714.RData")


#prepare the simulated data marker comparison for plotting.
temp.AFP <- line.info$AFP[line.info$genoID%in%lineF[,1]]
temp.Season <- line.info$Season[line.info$genoID%in%lineF[,1]]
dat.plot.sim.DH <- data.frame(AFP=rep(temp.AFP, 2),
                              Season=rep(temp.Season, 2),
                              Match=c(sapply(1:212, FUN=function(x) min(out.sim.DH[[x]])),
                                      sapply(1:212, FUN=function(x) max(out.sim.DH[[x]]))))
dat.plot.sim.F6 <- data.frame(AFP=rep(temp.AFP, 2),
                              Season=rep(temp.Season, 2),
                              Match=c(sapply(1:212, FUN=function(x) min(out.sim.F6[[x]])),
                                      sapply(1:212, FUN=function(x) max(out.sim.F6[[x]]))))
dat.plot.real.DH <- data.frame(AFP=temp.AFP,
                               Season=temp.Season,
                               out.real.DH[,1:3])
colnames(dat.plot.real.DH)[3:5] <- c("Progeny", "Parent1", "Parent2")
dat.plot.real.F6 <- data.frame(AFP=temp.AFP,
                               Season=temp.Season,
                               out.real.F6[,1:3])
colnames(dat.plot.real.F6)[3:5] <- c("Progeny", "Parent1", "Parent2")

#separate out the SB and WB.
dat.plot.sim.DH.SB <- dat.plot.sim.DH[dat.plot.sim.DH$Season=="SB",c(1,3)]
dat.plot.sim.DH.WB <- dat.plot.sim.DH[dat.plot.sim.DH$Season=="WB",c(1,3)]
dat.plot.sim.F6.SB <- dat.plot.sim.F6[dat.plot.sim.F6$Season=="SB",c(1,3)]
dat.plot.sim.F6.WB <- dat.plot.sim.F6[dat.plot.sim.F6$Season=="WB",c(1,3)]
dat.plot.real.DH.SB <- dat.plot.real.DH[dat.plot.real.DH$Season=="SB",c(1,3:5)]
dat.plot.real.DH.WB <- dat.plot.real.DH[dat.plot.real.DH$Season=="WB",c(1,3:5)]
dat.plot.real.F6.SB <- dat.plot.real.F6[dat.plot.real.F6$Season=="SB",c(1,3:5)]
dat.plot.real.F6.WB <- dat.plot.real.F6[dat.plot.real.F6$Season=="WB",c(1,3:5)]

#add dummy levels to the WB so the plots will come out with the same spacing between SB and WB.
dat.plot.real.DH.WB <- rbind(dat.plot.real.DH.WB,
                             data.frame(AFP=9001:9017, Progeny=NA, Parent1=NA, Parent2=NA))
dat.plot.real.F6.WB <- rbind(dat.plot.real.F6.WB,
                             data.frame(AFP=9001:9017, Progeny=NA, Parent1=NA, Parent2=NA))

#set AFP as factors.
dat.plot.sim.DH.SB$AFP <- as.factor(dat.plot.sim.DH.SB$AFP)
dat.plot.sim.DH.WB$AFP <- as.factor(dat.plot.sim.DH.WB$AFP)
dat.plot.sim.F6.SB$AFP <- as.factor(dat.plot.sim.F6.SB$AFP)
dat.plot.sim.F6.WB$AFP <- as.factor(dat.plot.sim.F6.WB$AFP)
dat.plot.real.DH.SB$AFP <- as.factor(dat.plot.real.DH.SB$AFP)
dat.plot.real.DH.WB$AFP <- as.factor(dat.plot.real.DH.WB$AFP)
dat.plot.real.F6.SB$AFP <- as.factor(dat.plot.real.F6.SB$AFP)
dat.plot.real.F6.WB$AFP <- as.factor(dat.plot.real.F6.WB$AFP)

#change the column names from "Parent1" and "Parent2" to "Parent".
colnames(dat.plot.real.DH.SB)[3:4] <- "Parent"
colnames(dat.plot.real.DH.WB)[3:4] <- "Parent"
colnames(dat.plot.real.F6.SB)[3:4] <- "Parent"
colnames(dat.plot.real.F6.WB)[3:4] <- "Parent"

#melt the data.frame for real data marker comparison.
dat.plot.real.DH.SB <- melt(dat.plot.real.DH.SB, id.vars="AFP")
dat.plot.real.DH.WB <- melt(dat.plot.real.DH.WB, id.vars="AFP")
dat.plot.real.F6.SB <- melt(dat.plot.real.F6.SB, id.vars="AFP")
dat.plot.real.F6.WB <- melt(dat.plot.real.F6.WB, id.vars="AFP")

#rename the columns.
colnames(dat.plot.real.DH.SB)[2:3] <- c("Type", "Match")
colnames(dat.plot.real.DH.WB)[2:3] <- c("Type", "Match")
colnames(dat.plot.real.F6.SB)[2:3] <- c("Type", "Match")
colnames(dat.plot.real.F6.WB)[2:3] <- c("Type", "Match")

#F6 + SB
ggplot() +
  geom_line(data=dat.plot.sim.F6.SB, aes(x=AFP, y=Match), color="#505050") +
  geom_point(data=dat.plot.real.F6.SB, aes(x=AFP, y=Match, color=Type), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_discrete(expand=c(0,0.5)) +
  scale_y_continuous(expand=c(0,0.2)) +
  xlab("Line") +
  ylab("Match Percentage") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_color_manual(values=c("#FF5050", "#5050FF")) +
  theme(legend.position=c(0.07, 0.85), legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  theme(legend.key=element_rect(fill=NA), legend.title=element_blank(), legend.text=element_text(size=8))

ggsave("output/Fig03A.png",
       width=7,
       height=2,
       units="in",
       dpi=600)
ggsave("output/Fig03A.svg",
       width=7,
       height=2,
       units="in",
       scale=4/3)

#F6 + WB
ggplot() +
  geom_line(data=dat.plot.sim.F6.WB, aes(x=AFP, y=Match), color="#505050") +
  geom_point(data=dat.plot.real.F6.WB, aes(x=AFP, y=Match, color=Type), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_discrete(expand=c(0,0.5)) +
  scale_y_continuous(expand=c(0,0.2), limits=c(0,8)) +
  xlab("Line") +
  ylab("Match Percentage") +
  annotate("rect", xmin=-Inf, xmax=97.5, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_color_manual(values=c("#FF5050", "#5050FF")) +
  theme(legend.position=c(0.07, 0.85), legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  theme(legend.key=element_rect(fill=NA), legend.title=element_blank(), legend.text=element_text(size=8))

ggsave("output/Fig03B.png",
       width=7,
       height=2,
       units="in",
       dpi=600)
ggsave("output/Fig03B.svg",
       width=7,
       height=2,
       units="in",
       scale=4/3)

#DH simulation results are highly similar to F6.
#DH simulation is not included in the manuscript.
#DH + SB
ggplot() +
  geom_line(data=dat.plot.sim.DH.SB, aes(x=AFP, y=Match), color="#505050") +
  geom_point(data=dat.plot.real.DH.SB, aes(x=AFP, y=Match, color=Type), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_discrete(expand=c(0,0.5)) +
  scale_y_continuous(expand=c(0,0.2)) +
  xlab("Line") +
  ylab("Match Percentage") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_color_manual(values=c("#FF5050", "#5050FF")) +
  theme(legend.position=c(0.07, 0.85), legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  theme(legend.key=element_rect(fill=NA), legend.title=element_blank(), legend.text=element_text(size=8))

ggsave("output/SFig00A.png",
       width=7,
       height=2,
       units="in",
       dpi=600)
ggsave("output/SFig00A.svg",
       width=7,
       height=2,
       units="in")

#DH + WB
ggplot() +
  geom_line(data=dat.plot.sim.DH.WB, aes(x=AFP, y=Match), color="#505050") +
  geom_point(data=dat.plot.real.DH.WB, aes(x=AFP, y=Match, color=Type), size=1) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_discrete(expand=c(0,0.5)) +
  scale_y_continuous(expand=c(0,0.2), limits=c(0,10)) +
  xlab("Line") +
  ylab("Match Percentage") +
  annotate("rect", xmin=-Inf, xmax=97.5, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_color_manual(values=c("#FF5050", "#5050FF")) +
  theme(legend.position=c(0.07, 0.85), legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  theme(legend.key=element_rect(fill=NA), legend.title=element_blank(), legend.text=element_text(size=8))

ggsave("output/SFig00B.png",
       width=7,
       height=2,
       units="in",
       dpi=600)
ggsave("output/SFig00B.svg",
       width=7,
       height=2,
       units="in")

#TabS06 - small marker set simulation in table.
check <- data.frame(Line=temp.Line,
                    AFP=temp.AFP,
                    Season=temp.Season,
                    out.real.F6[,1:3],
                    nseg=out.S,
                    nunique=sapply(1:length(out.sim.F6), FUN=function(x) length(out.sim.F6[[x]])),
                    matchpct=sapply(1:length(out.sim.F6), FUN=function(x) sum(out.sim.F6[[x]][out.sim.F6[[x]] > 1])))
colnames(check)[4:6] <- c("Progeny", "Parent1", "Parent2")
check <- check[order(check$Season, check$AFP),]
write.csv(check, "output/TableS06.csv", row.names=F, quote=F)

check[check$Season=="SB" & check$Parent2 > 5,]
#          Line  AFP Season Progeny Parent1 Parent2 nseg nunique matchpct
#209 LG Goddess 2999     SB 7.50758 7.97118 7.83114    4      72     88.4


#Tab02/FigS02,03,04/TabS03,04,05-GWAS#####
#read in marker position.
marker <- read.csv("source/SWB_marker_info.csv", header=T, as.is=T)
rownames(marker) <- marker$NewMarkerName
temp <- unname(c(0, sapply(unique(marker$Chr)[-7], FUN=function(x) max(marker$PhyPos[marker$Chr==x]) + 3e7)))
temp <- sapply(1:length(temp), FUN=function(x) sum(temp[1:x]))
marker$Pos <- unlist(sapply(1:7, FUN=function(x) marker$PhyPos[marker$Chr==unique(marker$Chr)[x]] + temp[x]))

#calculate the allele frequencies for the minor allele in each dataset.
x0 <- geno; f0 <- (colSums(x0==1) + 0.5*colSums(x0==0))/nrow(x0)
x1 <- geno[trait$genoID[trait$Season=="SB"],]; f1 <- (colSums(x1==1) + 0.5*colSums(x1==0))/nrow(x1)
x2 <- geno[trait$genoID[trait$Season=="WB"],]; f2 <- (colSums(x2==1) + 0.5*colSums(x2==0))/nrow(x2)
marker$FreqW <- marker$FreqS <- marker$FreqC <- NA
marker[names(f0),"FreqC"] <- f0
marker[names(f1),"FreqS"] <- f1
marker[names(f2),"FreqW"] <- f2

#function to generate Manhattan plot from sommer GWAS output.
plotM <- function(g, marker, h=1, w=1.75, filename){
  
  #get the -log10p scores.
  out <- t(g$scores)[,2]
  
  #match the marker position information with the -log10p score.
  out <- data.frame(Score=out,
                    q=qvalue(10^-out)$qvalues,
                    Chr=marker[names(out), "Chr"],
                    Pos=marker[names(out), "Pos"],
                    Type="Combined")
  
  #set the alternating colors for each chromosome.
  out$Color <- "A"
  out$Color[out$Chr%in%c("2H","4H","6H")] <- "B"
  
  #get the FDR = 0.05 thresholds, ignore if there is 0 or 1 hit above the threshold.
  if(sum(out$q<0.05)>1){
    t <- out[out$q<0.05,]
    t <- t[order(t$q),]
    t <- t$Score[nrow(t)]
    t <- c(-Inf, Inf, t)
  } else {
    t <- rep(Inf, 3)
  }
  
  #identify the positions on x-axis to display chromosome names.
  temp <- sapply(unique(marker$Chr), FUN=function(x) (max(marker$Pos[marker$Chr==x])+min(marker$Pos[marker$Chr==x]))/2)
  
  p <- ggplot() +
    geom_point(data=out, aes(x=Pos, y=Score, color=Color), size=0.5) +
    annotate("segment", x=t[1], xend=t[2], y=t[3], yend=t[3], color="#FF5050") +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text=element_text(size=7, margin=margin(0,0,0,0))) +
    theme(axis.title.y=element_text(size=7, margin=margin(0,0,0,0))) +
    scale_x_continuous(expand=c(0.01,0), breaks=temp, labels=names(temp)) +
    scale_color_manual(values=c("#202020", "#AAAAAA")) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
    ylab(bquote("-log"[10]~"p"))
  
  ggsave(plot=p,
         filename=filename,
         height=h,
         width=w,
         units="in",
         dpi=600)
}

#GWAS on t01-t27 in the Combined dataset.
for(k in 1:27){
  y <- trait[!is.na(trait[,12+k]), c(1,4,5,12+k)]
  colnames(y)[4] <- "t"
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year + Season,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Combined_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Combined_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Combined_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#GWAS on t28 in the Combined dataset (without Season covariate).
for(k in 28){
  y <- trait[!is.na(trait[,12+k]), c(1,4,5,12+k)]
  colnames(y)[4] <- "t"
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Combined_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Combined_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Combined_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#GWAS on t01,t02,t04-t11,t13-t27 in the Spring dataset.
for(k in c(1:2,4:11,13:27)){
  y <- trait[!is.na(trait[,12+k]) & trait$Season=="SB", c(1,4,12+k)]
  colnames(y)[3] <- "t"
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Spring_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Spring_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Spring_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#GWAS on t19 in the Spring dataset (sommer bug, recoded trait from 1/2 to 0/1).
for(k in c(19)){
  y <- trait[!is.na(trait[,12+k]) & trait$Season=="SB", c(1,4,12+k)]
  colnames(y)[3] <- "t"
  y$t <- y$t - 1
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Spring_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Spring_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Spring_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#GWAS on t01-t11, t13-t19, t21-t26 in the Winter dataset.
for(k in c(1:11,13:19,21:26)){
  y <- trait[!is.na(trait[,12+k]) & trait$Season=="WB", c(1,4,6,12+k)]
  colnames(y)[4] <- "t"
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year + RowType,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Winter_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Winter_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Winter_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#GWAS on t12 and t20 in the Winter dataset (without RowType covariate).
for(k in c(12,20)){
  y <- trait[!is.na(trait[,12+k]) & trait$Season=="WB", c(1,4,6,12+k)]
  colnames(y)[4] <- "t"
  x <- geno[y$genoID,]
  x <- x[,!(colSums(x==-1)==nrow(x) | colSums(x==0)==nrow(x) | colSums(x==1)==nrow(x))]
  z <- A.mat(x)
  g <- GWAS(t~Year,
            random=~vs(genoID, Gu=z),
            rcov=~units,
            data=y,
            M=x,
            gTerm="u:genoID")
  plotM(g, marker, h=1, w=1.75, paste(c("output/GWAS/GWAS_Winter_Small_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  plotM(g, marker, h=2, w=3.50, paste(c("output/GWAS/GWAS_Winter_Large_t",rep(0,2-nchar(k)),k,".png"), collapse=""))
  saveRDS(g, file=paste(c("temp/GWAS/GWAS_Winter_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
}

#compile all GWAS results for the Combined dataset.
for(k in 1:28){
  g <- readRDS(paste(c("temp/GWAS/GWAS_Combined_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
  g <- data.frame(t(g$scores)[,1:2])
  colnames(g) <- c("Effect", "Score")
  g <- data.frame(marker[rownames(g), c("Chr","PhyPos","major","minor","FreqC","FreqS","FreqW")], g)
  g$q <- qvalue(10^-g$Score)$qvalues
  g <- g[g$q<0.10,]
  write.csv(g, paste(c("temp/GWAS/Hits_Combined_t",rep(0,2-nchar(k)),k,".csv"), collapse=""), quote=F, row.names=T)
}

#compile all GWAS results for the Spring dataset.
for(k in c(1:2,4:11,13:27)){
  g <- readRDS(paste(c("temp/GWAS/GWAS_Spring_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
  g <- data.frame(t(g$scores)[,1:2])
  colnames(g) <- c("Effect", "Score")
  g <- data.frame(marker[rownames(g), c("Chr","PhyPos","major","minor","FreqC","FreqS","FreqW")], g)
  g$q <- qvalue(10^-g$Score)$qvalues
  g <- g[g$q<0.10,]
  write.csv(g, paste(c("temp/GWAS/Hits_Spring_t",rep(0,2-nchar(k)),k,".csv"), collapse=""), quote=F, row.names=T)
}

#compile all GWAS results for the Winter dataset.
for(k in 1:26){
  g <- readRDS(paste(c("temp/GWAS/GWAS_Winter_t",rep(0,2-nchar(k)),k,".RDS"), collapse=""))
  g <- data.frame(t(g$scores)[,1:2])
  colnames(g) <- c("Effect", "Score")
  g <- data.frame(marker[rownames(g), c("Chr","PhyPos","major","minor","FreqC","FreqS","FreqW")], g)
  g$q <- qvalue(10^-g$Score)$qvalues
  g <- g[g$q<0.10,]
  write.csv(g, paste(c("temp/GWAS/Hits_Winter_t",rep(0,2-nchar(k)),k,".csv"), collapse=""), quote=F, row.names=T)
}


#Fig04-Manhattan distances from trait and marker data#####
#calculate Manhattan distance from all marker data.
temp <- unique(c(round(10^seq(0,4.5,0.1)), ncol(geno)))
dist.geno <- sapply(temp, FUN=function(y) c(dist(x=geno[,sample(1:ncol(geno), y, replace=F)], method="manhattan"))/y)
saveRDS(dist.geno, "temp/dist_geno.RDS")

#calculate Manhattan distance from 28 DUS traits.
dist.trait <- c(dist(x=as.matrix(trait[,13:40]), method="manhattan"))/28
saveRDS(dist.trait, "temp/dist_trait.RDS")

#identify traits with low and high h2.
out.varcomp <- read.csv("temp/univariate_varcomp.csv", header=T, as.is=T)
temp.lo <- which(out.varcomp$H2C < 0.5)
temp.hi <- which(out.varcomp$H2C >= 0.5)

#calculate Manhattan distance from 15 low h2 DUS traits.
dist.trait.lo <- c(dist(x=as.matrix(trait[,12+temp.lo]), method="manhattan"))/15
saveRDS(dist.trait.lo, "temp/dist_trait_lo.RDS")

#calculate Manhattan distance from 13 high h2 DUS traits.
dist.trait.hi <- c(dist(x=as.matrix(trait[,12+temp.hi]), method="manhattan"))/13
saveRDS(dist.trait.hi, "temp/dist_trait_hi.RDS")

#calculate correlations between the Manhattan distances from marker and trait data.
cor.all <- sapply(1:ncol(dist.geno), FUN=function(x) cor(dist.geno[,x], dist.trait))
cor.lo <- sapply(1:ncol(dist.geno), FUN=function(x) cor(dist.geno[,x], dist.trait.lo, use="complete.obs"))
cor.hi <- sapply(1:ncol(dist.geno), FUN=function(x) cor(dist.geno[,x], dist.trait.hi))

#prepare the data for plotting.
dat.plot <- data.frame(marker=temp,
                       cor.all=cor.all,
                       cor.lo=cor.lo,
                       cor.hi=cor.hi)
dat.plot <- melt(dat.plot, id.vars="marker")

ggplot() +
  geom_point(data=dat.plot, aes(x=marker, y=value, color=variable), size=1) +
  scale_x_continuous(trans="log10") +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  theme(legend.key=element_rect(fill=NA), legend.title=element_text(size=9), legend.text=element_text(size=8)) +
  theme(legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1)) +
  scale_color_manual(name="Traits", labels=c("All", "Low h2", "High h2"), values=c("#505050", "#FF5050", "#5050FF")) +
  xlab("Marker Number") +
  ylab("Correlation")

ggsave(filename="output/Fig04A.png",
       height=1.5,
       width=3.4,
       units="in",
       dpi=600)
ggsave(filename="output/Fig04A.svg",
       height=1.5,
       width=3.4,
       units="in",
       scale=4/3)

#retrieve marker data and information.
geno.full <- readRDS("source/SWB_geno_012.RDS")
geno.full <- geno.full[line.info$genoID, colSums(is.na(geno.full))==0]
marker <- read.csv("source/SWB_marker_info.csv", header=T, as.is=T)
rownames(marker) <- marker$NewMarkerName
marker <- marker[colnames(geno.full),]

#identify which parent-pair has the lowest and highest Manhattan distances.
dist.geno <- dist(x=geno.full, method="manhattan")/ncol(geno.full)
dist.geno <- as.matrix(dist.geno)
dist.geno[upper.tri(dist.geno, diag=T)] <- NA
dist.geno <- melt(dist.geno, na.rm=T)

lineF <- line.info[which(line.info$CrossType==1),c(1:2,4:8,10:13)]
lineF <- lineF[!(rowSums(is.na(lineF[,10:11]))>0),]
lineF$temp <- paste(lineF$ParentID1, lineF$ParentID2, sep="_")

temp.dist <- rep(dist.geno$value,2)
names(temp.dist) <- c(paste(dist.geno$Var1, dist.geno$Var2, sep="_"),
                      paste(dist.geno$Var2, dist.geno$Var1, sep="_"))
temp.dist <- temp.dist[lineF$temp]
lineF$dist <- unname(temp.dist)
lineF <- lineF[,-12]

saveRDS(lineF, "temp/dist_parents.RDS")

#let's use SB as example.
lineF <- lineF[lineF$Season=="SB",]

lineF[which.min(lineF$dist),]
#       genoID   Line  AFP Year Season RowType   Breeder Parent1 Parent2 ParentID1 ParentID2      dist
#675 ID_32A487 Acumen 2623 2011     SB      2R Ackermann Propino  Quench ID_22A288  ID_1A132 0.2041519

lineF[which.max(lineF$dist),]
#     genoID    Line  AFP Year Season RowType Breeder Parent1 Parent2 ParentID1 ParentID2      dist
#245 ID_6A80 Berwick 1474 1997     SB      2R    RAGT Riviera  Cooper  ID_1A136  ID_6A108 0.5878537


###simulate 1,000 progeny from two parents with lowest Manhattan distance (Propino x Quench).
#get the marker data for Propino and Quench.
temp.geno <- geno.full[c("ID_22A288","ID_1A132"),]

#identify heterozygous markers and temporarily assign them as 0, and set to NA at the end of the simulation.
het.idx <- which(colSums(temp.geno==1) > 0)
temp.geno[,het.idx] <- 0

#separate the marker data by chromosome.
parentHap <- lapply(unique(marker$Chr), FUN=function(x) temp.geno[, marker$Chr==x, drop=F])
genMap <- lapply(unique(marker$Chr), FUN=function(x) marker$GenPos[marker$Chr==x]/100)
genMap <- lapply(1:length(genMap), FUN=function(x) genMap[[x]] - genMap[[x]][1])
  
#initiate the founder haplotype based on the two parents.
founder <- newMapPop(genMap=genMap,
                     haplotypes=parentHap,
                     inbred=T,
                     ploidy=2L)
SP <- SimParam$new(founder)

#create the founder population.
founder <- newPop(founder, simParam=SP)
  
#make a cross between the two parents.
F1 <- randCross(pop=founder,
                nCrosses=1,
                nProgeny=1,
                simParam=SP)
  
#simulate 1,000 F6 progeny from the cross.
#self the F1 and create 1,000 progeny.
F2 <- self(pop=F1,
           nProgeny=1000,
           simParam=SP)
  
#self the F2 and keep 1 progeny each.
F3 <- self(pop=F2,
           nProgeny=1,
           simParam=SP)
  
#self the F3 and keep 1 progeny each.
F4 <- self(pop=F3,
           nProgeny=1,
           simParam=SP)
  
#self the F4 and keep 1 progeny each.
F5 <- self(pop=F4,
           nProgeny=1,
           simParam=SP)
  
#self the F5 and keep 1 progeny each.
F6 <- self(pop=F5,
           nProgeny=1,
           simParam=SP)
  
#extract the marker genotypes for the simulated progeny.
geno.F6 <- pullSegSiteGeno(pop=F6,
                           chr=NULL,
                           simParam=SP)
colnames(geno.F6) <- colnames(temp.geno)
geno.F6[,het.idx] <- NA

#combine the marker data from the full 805 "reference" panel and 1,000 simulated progeny.
geno.F6 <- rbind(geno.full, geno.F6)
dist.F6.low <- as.matrix(dist(geno.F6, method="manhattan"))/ncol(geno.F6)
saveRDS(dist.F6.low, "temp/dist_Propino_Quench_F6.RDS")

#simulate 1,000 BC4S2 progeny from the cross.
#backcross F1 to parent 1 and create 1000 progeny.
BC1 <- makeCross2(females=founder,
                  males=F1,
                  crossPlan=matrix(c(1,1),1,2),
                  nProgeny=1000,
                  simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S1 <- self(pop=BC1,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S2 <- self(pop=BC1S1,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S3 <- self(pop=BC1S2,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S4 <- self(pop=BC1S3,
              nProgeny=1,
              simParam=SP)


#extract the marker genotypes for the simulated progeny.
geno.BC1S4 <- pullSegSiteGeno(pop=BC1S4,
                              chr=NULL,
                              simParam=SP)
colnames(geno.BC1S4) <- colnames(temp.geno)
geno.BC1S4[,het.idx] <- NA

#combine the marker data from the full 805 "reference" panel and 1,000 simulated progeny.
geno.BC1S4 <- rbind(geno.full, geno.BC1S4)
dist.BC1S4.low <- as.matrix(dist(geno.BC1S4, method="manhattan"))/ncol(geno.BC1S4)
saveRDS(dist.BC1S4.low, "temp/dist_Propino_Quench_BC1S4.RDS")





###simulate 1,000 progeny from two parents with highest Manhattan distance (Riviera x Cooper).
#get the marker data for Riviera and Cooper.
temp.geno <- geno.full[c("ID_1A136", "ID_6A108"),]

#identify heterozygous markers and temporarily assign them as 0, and set to NA at the end of the simulation.
het.idx <- which(colSums(temp.geno==1) > 0)
temp.geno[,het.idx] <- 0

#separate the marker data by chromosome.
parentHap <- lapply(unique(marker$Chr), FUN=function(x) temp.geno[, marker$Chr==x, drop=F])
genMap <- lapply(unique(marker$Chr), FUN=function(x) marker$GenPos[marker$Chr==x]/100)
genMap <- lapply(1:length(genMap), FUN=function(x) genMap[[x]] - genMap[[x]][1])

#initiate the founder haplotype based on the two parents.
founder <- newMapPop(genMap=genMap,
                     haplotypes=parentHap,
                     inbred=T,
                     ploidy=2L)
SP <- SimParam$new(founder)

#create the founder population.
founder <- newPop(founder, simParam=SP)

#make a cross between the two parents.
F1 <- randCross(pop=founder,
                nCrosses=1,
                nProgeny=1,
                simParam=SP)

#simulate 1,000 F6 progeny from the cross.
#self the F1 and create 1,000 progeny.
F2 <- self(pop=F1,
           nProgeny=1000,
           simParam=SP)

#self the F2 and keep 1 progeny each.
F3 <- self(pop=F2,
           nProgeny=1,
           simParam=SP)

#self the F3 and keep 1 progeny each.
F4 <- self(pop=F3,
           nProgeny=1,
           simParam=SP)

#self the F4 and keep 1 progeny each.
F5 <- self(pop=F4,
           nProgeny=1,
           simParam=SP)

#self the F5 and keep 1 progeny each.
F6 <- self(pop=F5,
           nProgeny=1,
           simParam=SP)

#extract the marker genotypes for the simulated progeny.
geno.F6 <- pullSegSiteGeno(pop=F6,
                           chr=NULL,
                           simParam=SP)
colnames(geno.F6) <- colnames(temp.geno)
geno.F6[,het.idx] <- NA

#combine the marker data from the full 805 "reference" panel and 1,000 simulated progeny.
geno.F6 <- rbind(geno.full, geno.F6)
dist.F6.high <- as.matrix(dist(geno.F6, method="manhattan"))/ncol(geno.F6)
saveRDS(dist.F6.high, "temp/dist_Riviera_Cooper_F6.RDS")

#simulate 1,000 BC4S2 progeny from the cross.
#backcross F1 to parent 1 and create 1000 progeny.
BC1 <- makeCross2(females=founder,
                  males=F1,
                  crossPlan=matrix(c(1,1),1,2),
                  nProgeny=1000,
                  simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S1 <- self(pop=BC1,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S2 <- self(pop=BC1S1,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S3 <- self(pop=BC1S2,
              nProgeny=1,
              simParam=SP)

#self the BC1 and keep 1 progeny each.
BC1S4 <- self(pop=BC1S3,
              nProgeny=1,
              simParam=SP)


#extract the marker genotypes for the simulated progeny.
geno.BC1S4 <- pullSegSiteGeno(pop=BC1S4,
                              chr=NULL,
                              simParam=SP)
colnames(geno.BC1S4) <- colnames(temp.geno)
geno.BC1S4[,het.idx] <- NA

#combine the marker data from the full 805 "reference" panel and 1,000 simulated progeny.
geno.BC1S4 <- rbind(geno.full, geno.BC1S4)
dist.BC1S4.high <- as.matrix(dist(geno.BC1S4, method="manhattan"))/ncol(geno.BC1S4)
saveRDS(dist.BC1S4.high, "temp/dist_Riviera_Cooper_BC1S4.RDS")


#set the upper triangle and diagonal to NA.
dist.F6.low[upper.tri(dist.F6.low, diag=T)] <- NA
dist.BC1S4.low[upper.tri(dist.BC1S4.low, diag=T)] <- NA
dist.F6.high[upper.tri(dist.F6.high, diag=T)] <- NA
dist.BC1S4.high[upper.tri(dist.BC1S4.high, diag=T)] <- NA

#separate out the distances.
#1 is distances among the 805 "reference" panel.
#2 is distances between 1,000 simulated progeny and 805 "reference" panel.
#3 is distances among the 1,000 simulated progeny.
dist.F6.low.1 <- dist.F6.low[1:805, 1:805]
dist.F6.low.2 <- dist.F6.low[806:1805, 1:805]
dist.F6.low.3 <- dist.F6.low[806:1805, 806:1805]

dist.BC1S4.low.1 <- dist.BC1S4.low[1:805, 1:805]
dist.BC1S4.low.2 <- dist.BC1S4.low[806:1805, 1:805]
dist.BC1S4.low.3 <- dist.BC1S4.low[806:1805, 806:1805]

dist.F6.high.1 <- dist.F6.high[1:805, 1:805]
dist.F6.high.2 <- dist.F6.high[806:1805, 1:805]
dist.F6.high.3 <- dist.F6.high[806:1805, 806:1805]

dist.BC1S4.high.1 <- dist.BC1S4.high[1:805, 1:805]
dist.BC1S4.high.2 <- dist.BC1S4.high[806:1805, 1:805]
dist.BC1S4.high.3 <- dist.BC1S4.high[806:1805, 806:1805]

#Fig04B - distribution of distances among the 805 "reference" panel.
ggplot() +
  geom_density(aes(dist.F6.low.1[lower.tri(dist.F6.low.1, diag=F)]), fill="#FFFF55", alpha=0.5, color=NA) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_x_continuous(breaks=seq(0,1,0.2), name="Manhattan distance") +
  scale_y_continuous(expand=c(0.02,0)) +
  annotate("segment", x=0.05, xend=0.05, y=-Inf, yend=Inf, color="#FF5555") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999")

ggsave(filename="output/Fig04B.svg",
       width=3.4,
       height=1.5,
       units="in",
       scale=4/3)
ggsave(filename="output/Fig04B.png",
       width=3.4,
       height=1.5,
       units="in",
       dpi=600)

#Fig04C - distribution of minimum distances between 1,000 simulated progeny and 805 "reference" panel.
dat.plot <- data.frame(sim1=sapply(1:nrow(dist.F6.low.2), FUN=function(x) min(dist.F6.low.2[x,])),
                       sim2=sapply(1:nrow(dist.BC1S4.low.2), FUN=function(x) min(dist.BC1S4.low.2[x,])),
                       sim3=sapply(1:nrow(dist.F6.high.2), FUN=function(x) min(dist.F6.high.2[x,])),
                       sim4=sapply(1:nrow(dist.BC1S4.high.2), FUN=function(x) min(dist.BC1S4.high.2[x,])))
dat.plot <- melt(dat.plot)

ggplot() +
  geom_density(data=dat.plot, aes(x=value, fill=variable), color=NA, alpha=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_x_continuous(breaks=seq(0,1,0.1), name="Manhattan distance") +
  scale_y_continuous(expand=c(0.02,0)) +
  annotate("segment", x=0.05, xend=0.05, y=-Inf, yend=Inf, color="#FF5555") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_fill_manual(values=c("#55FFFF", "#55FF55", "#5555FF", "#FF55FF"), labels=c("F6 (low)", "BC1S4 (low)", "F6 (high)", "BC1S4 (high)"), name="Simulation") +
  theme(legend.title=element_text(size=9), legend.text=element_text(size=8)) +
  theme(legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1))  

ggsave(filename="output/Fig04C.svg",
       width=3.4,
       height=1.5,
       units="in",
       scale=4/3)
ggsave(filename="output/Fig04C.png",
       width=3.4,
       height=1.5,
       units="in",
       dpi=600)

#Fig04D - distribution of distances among the 805 "reference" panel.
dat.plot <- data.frame(sim1=dist.F6.low.3[lower.tri(dist.F6.low.3, diag=F)],
                       sim2=dist.BC1S4.low.3[lower.tri(dist.BC1S4.low.3, diag=F)],
                       sim3=dist.F6.high.3[lower.tri(dist.F6.high.3, diag=F)],
                       sim4=dist.BC1S4.high.3[lower.tri(dist.BC1S4.high.3, diag=F)])
dat.plot <- melt(dat.plot)

ggplot() +
  geom_density(data=dat.plot, aes(x=value, fill=variable), color=NA, alpha=0.5) +
  theme(panel.background=element_blank(), panel.grid=element_blank()) +
  scale_x_continuous(breaks=seq(0,1,0.1), name="Manhattan distance") +
  scale_y_continuous(expand=c(0.02,0)) +
  annotate("segment", x=0.05, xend=0.05, y=-Inf, yend=Inf, color="#FF5555") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#999999") +
  scale_fill_manual(values=c("#55FFFF", "#55FF55", "#5555FF", "#FF55FF"), labels=c("F6 (low)", "BC1S4 (low)", "F6 (high)", "BC1S4 (high)"), name="Simulation") +
  theme(legend.title=element_text(size=9), legend.text=element_text(size=8)) +
  theme(legend.key.size=unit(0.8,"lines"), legend.margin=margin(1,1,1,1))  

ggsave(filename="output/Fig04D.svg",
       width=3.4,
       height=1.5,
       units="in",
       scale=4/3)
ggsave(filename="output/Fig04D.png",
       width=3.4,
       height=1.5,
       units="in",
       dpi=600)


#Check for current varieties that have Manhattan distances below 0.05.
check <- dist.F6.low.1
rownames(check) <- colnames(check) <- paste(line.info$Season, line.info$Line, line.info$AFP, sep="_")
check <- melt(check, na.rm=T)
check <- check[order(check$value),]
check[check$value<0.05,]
#                  Var1              Var2      value
#427288 WB_KWS Joy_2521  WB_Wintmalt_2135 0.03632916
#227385   SB_Class_1777  SB_Prestige_1566 0.03825041
#58036   WB_Angora_1165   WB_Melanie_1160 0.03952293
#572315  WB_Mackie_2955 WB_KWS Tower_2725 0.04393932

#Check for percentage of simulated progeny vs "reference" panel with Manhattan distances below 0.05.
check <- data.frame(sim1=sapply(1:nrow(dist.F6.low.2), FUN=function(x) min(dist.F6.low.2[x,])),
                    sim2=sapply(1:nrow(dist.BC1S4.low.2), FUN=function(x) min(dist.BC1S4.low.2[x,])),
                    sim3=sapply(1:nrow(dist.F6.high.2), FUN=function(x) min(dist.F6.high.2[x,])),
                    sim4=sapply(1:nrow(dist.BC1S4.high.2), FUN=function(x) min(dist.BC1S4.high.2[x,])))
colSums(check<0.05)/nrow(check)*100
#sim1 sim2 sim3 sim4 
#13.0 59.6  0.0  4.9


#pull out the distribution of Manhattan distances for SB-SB, SB-WB, WB-WB.
check <- dist.F6.low.1
rownames(check) <- colnames(check) <- line.info$Season
check <- melt(check, na.rm=T)
check$pair <- 9
check$pair[rowSums(check[,1:2]=="SB")==2] <- 0
check$pair[rowSums(check[,1:2]=="WB")==2] <- 2
check$pair[rowSums(check[,1:2]=="SB")==1 & rowSums(check[,1:2]=="WB")==1] <- 1

summary(check$value[check$pair==0]) #SB-SB
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03825 0.40748 0.45464 0.44864 0.49673 0.68926

summary(check$value[check$pair==2]) #WB-WB
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03633 0.45399 0.51068 0.50979 0.57166 0.86853

summary(check$value[check$pair==1]) #SB-WB
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4391  0.6910  0.7240  0.7317  0.7597  0.9657


temp.Line <- line.info$Line[line.info$genoID%in%lineF[,1]]
temp.AFP <- line.info$AFP[line.info$genoID%in%lineF[,1]]
temp.Season <- line.info$Season[line.info$genoID%in%lineF[,1]]

#TabS07 - Manhattan distances for different markers.
dist.geno <- readRDS("temp/dist_geno.RDS")
check <- t(sapply(1:ncol(dist.geno), FUN=function(x) c(mean(dist.geno[,x]), var(dist.geno[,x]))))
check <- data.frame(nmarker=unique(c(round(10^seq(0,4.5,0.1)), ncol(geno))),
                    mean=check[,1],
                    var=check[,2])
write.csv(check, "output/TableS07.csv", row.names=F, quote=F)