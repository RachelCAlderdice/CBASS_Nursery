library(emmeans)
library(drc)
library(plyr)
library(gplots)
library(lattice)

#Sort input file
pamdata<-read.delim("Inputfiles/Metadata_DRC_modelling_ahyc.csv", header=T, sep = ",")
pamdata$Geno<-as.factor(pamdata$Geno)
pamdata$Temp<-as.numeric(pamdata$Temp)
pamdata$Geno_Source <- paste(pamdata$Geno, pamdata$Source, sep = "_")
pamdata$Geno_Source<-as.factor(pamdata$Geno_Source)
aggregate(PAM ~ Species_Source + Temp, data=pamdata, summary)

# Checking sample sizes
aggregate(PAM ~ Species_Source + Temp, data=pamdata, length)
aggregate(PAM ~ Species_Source + Temp + Geno, data=pamdata, length)

# DRC Curve Fitting considering all colonies
DRCpam = drm(PAM ~ Temp, data = pamdata, curveid = Geno_Source,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpam)
compParm(DRCpam, 'ed50')
compParm(DRCpam, 'ed50', "-")
#pdf("Results/CBASS_nursery_ALL_colonies.pdf",10,7)
#plot(DRCpam)
#title(main="CBASS_nursery_ALL_colonies")
#dev.off()

plot(DRCpam)
points(pamdata$Temp, pamdata$PAM)
ED(DRCpam, c(50))

# fit to each Species individually - Donor 
DRCpamAhyc = drm(PAM ~ Temp, data = pamdata[pamdata$Source=="D",],
                 fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAhyc)
DRCpamAhyc$coefficients[3]
ED(DRCpamAhyc, c(50))

# genotype-specific curve fits - Donor
DRCpamAhycgeno = drm(PAM ~ Temp, data = pamdata[pamdata$Source=="D",], curveid=Geno,
                     fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAhycgeno)
DRCpamAhycgeno$coefficients[11:15]
ED(DRCpamAhycgeno, c(50))
#pdf("Results/CBASS_nursery_D_colonies.pdf",10,7)
#plot(DRCpamAhycgeno)
#title(main="CBASS_nursery_D_colonies")
#dev.off()

# fit to each Species individually - Nursery
DRCpamAhyc2 = drm(PAM ~ Temp, data = pamdata[pamdata$Source=="N",],
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAhyc2)
DRCpamAhyc2$coefficients[3]
ED(DRCpamAhyc2, c(50))

# genotype-specific curve fits - Nursery
DRCpamAhycgeno2 = drm(PAM ~ Temp, data = pamdata[pamdata$Source=="N",], curveid=Geno,
                      fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamAhycgeno2)
DRCpamAhycgeno2$coefficients[11:15]
ED(DRCpamAhycgeno2, c(50))
#pdf("Results/CBASS_nursery_N_colonies.pdf",10,7)
#plot(DRCpamAhycgeno2)
#title(main="CBASS_nursery_N_colonies")
#dev.off()

#collate coefficents
Coeffs<-c(DRCpamAhyc$coefficients[3], DRCpamAhyc2$coefficients[3])
GenoCoeffs<-data.frame("ED50"=c(DRCpamAhycgeno$coefficients[11:15],DRCpamAhycgeno2$coefficients[11:15]))
GenoCoeffs$Source=c(rep("D",5), rep("N",5))
aggregate(ED50 ~ Source, data=GenoCoeffs, mean)

#statistics
aggregate(ED50 ~ Source, data=GenoCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Source, data=GenoCoeffs)
t.test(ED50 ~ Source, GenoCoeffs)

#tabulate
aggregate(ED50 ~ Source, data=GenoCoeffs, summary)
GenoCoeffs
SumaryStats<-data.frame("Source"=names(tapply(GenoCoeffs$ED50,GenoCoeffs$Source, mean)), "MeanED50"=tapply(GenoCoeffs$ED50,GenoCoeffs$Source, mean), "ED50StdDev"=tapply(GenoCoeffs$ED50,GenoCoeffs$Source, sd), "ED50StdErr"=tapply(GenoCoeffs$ED50,GenoCoeffs$Source, sd)/sqrt(tapply(GenoCoeffs$ED50,GenoCoeffs$Source, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
write.table(data.frame("Genet"=gsub("ed50:", "", row.names(GenoCoeffs)), GenoCoeffs), file="Results/CBASS_nursery_ahya_ED50_genets.txt", quote=F, sep="\t", row.names=F)
write.table(SumaryStats, file="Results/CBASS_nursery_ahya_ED50_summarystats.txt", quote=F, sep="\t", row.names=F)

# Summary curve plot 
temp_x<- seq(28, 40, length = 100)
pdf("Results/CBASS_nursery_ahya_curve_plot.pdf",7,7)
line_width=3
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)

i<-1 #Donors
matplot(temp_x, predict(DRCpamAhyc, data.frame(Temp = temp_x), interval="confidence"),
        type="l",col="steelblue1",lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(29.5,39.5),ylim=c(-0.15 ,0.9), cex.axis=1.5, cex.lab=1.5)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamAhyc, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamAhyc, data.frame(Temp = temp_x), interval="confidence")[,3])), col= adjustcolor("steelblue1", alpha.f=0.1), density = 30)
with(pamdata[pamdata$Source=="D",],matpoints(Temp-offsets[i],PAM,pch=20, col="steelblue1", cex=2))

i<-1.5 #Nursery
matpoints(temp_x, predict(DRCpamAhyc2, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col="blue",lty=c(1,3,3),lwd=line_width)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamAhyc2, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamAhyc2, data.frame(Temp = temp_x), interval="confidence")[,3])), col= adjustcolor("blue", alpha.f=0.1), density = 30)
with(pamdata[pamdata$Source=="N",],matpoints(Temp-offsets[i],PAM,pch=20, col="blue", cex=2))

legend("bottomleft",c("D ED50: 37.26°C","N ED50: 36.44°C"),pch=c(20,18), col=c("steelblue1","blue"),pt.cex=2, bty="n",cex=1.5)
dev.off()

############################# boxplot of ED50s ###########################################

input <-read.delim("Inputfiles/Metadata_ED50_boxplot.txt", header= T)
col <- c("steelblue1","blue")

input$Source <- as.character(input$Source)
input$Source <- factor(input$Source, levels= unique(input$Source))

pdf("Results/Boxplot_ed50s_genotype_ahya.pdf", 5,5)
ggplot(data=input, aes(x=Source, y=ED50, fill= Source))+
  geom_boxplot()+
  scale_fill_manual(values= col)+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=22))+
  xlab(paste0("Source"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())+
  theme(text = element_text(size=20), axis.text.x = element_text(vjust = 1, hjust=1))

dev.off()

