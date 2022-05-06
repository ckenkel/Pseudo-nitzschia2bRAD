#devtools::install_github("baptiste/egg")
devtools::install_github("road2stat/ggsci")

library(lme4)
library(adegenet)
library(ggplot2)
library(car)
library(multcomp)
library(gtools)
library(plyr)
library(quantreg)
library(calibrate)
library(MASS)
library(AICcmodavg)
library(e1071)
library(nlme)
library(MCMCglmm)
library(labdsv)
library(vegan)
library(plotrix)
library(pgirmess)
library(gridExtra)
library(egg)
library(scales)
library(vegan)
library(MCMC.OTU)
library(ggfortify)
library(cluster)
library(labdsv)
library(wesanderson)
library(tidyverse)
library(IRanges)
library(RColorBrewer)
library(ggsci)

#############################

setwd("") #location of input files - fill in here

#######################
#First, analysis of contaminants from blast match
########################

dat=read.table("Paus_bcgIN1sam.fasta.br.all.stripped.tab",sep="\t",header=FALSE, stringsAsFactors=TRUE)

summary(dat$V4) #use this to look at contaminant counts then to add to summary spreadsheet; note that 'no match' is difference between identified contaminants and number of original tags in the fasta file 

write.csv(summary(dat$V4),file="Ppuncontams.csv",quote=FALSE)

#### sort and tabulate contaminant categories offline and reload summed csv

tab<-read.csv("AllContaminantCounts.csv",stringsAsFactors=TRUE)

summary(tab)

###############pie charts for contaminant breakdown

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

spp<-tab[tab$Species=="Psub",2:3] #change to appropriate spp or subset

pie<-ggplot(spp,aes(x="",y=Ntags,fill=Category))+
geom_bar(width=1,stat="identity")+
coord_polar("y",start=0)

pie+scale_fill_brewer(palette="Dark2")+blank_theme+ #change color theme
theme(axis.text.x=element_blank()) +
  geom_text(aes(label = paste(Ntags),x=1.3),
            position = position_stack(vjust = 0.5))

#to plot as a percent instead use:
geom_text(aes(label = paste(round(Ntags / sum(Ntags) * 100, 1), "%"), x=1.3),
            position = position_stack(vjust = 0.5))


spp<-tab[tab$Species=="PdelOther",2:3] #change to appropriate spp or subset

pie<-ggplot(spp,aes(x="",y=Ntags,fill=Category))+
geom_bar(width=1,stat="identity")+
coord_polar("y",start=0)

pie+scale_fill_brewer(palette="Purples")+blank_theme+ #change color theme
theme(axis.text.x=element_blank()) +
  geom_text(aes(label = paste(Ntags)),
            position = position_stack(vjust = 0.5))
            
############## 
#######################
#Next, analysis of Accuracy, precision, FP rates, etc
########################
#FP rate

fp<-read.csv("FPrate.csv",stringsAsFactors=TRUE)
head(fp)
str(fp)
summary(fp)

fpH<-fp[fp$nHQreads>1000,] #remove samples with less than 1000 HQ reads remaining

str(fpH)
summary(fpH)

pd <- position_dodge(.2)
p2<-ggplot(fpH,aes(factor(MixBinary),FPR_relToHQMapReads))+
	geom_boxplot(outlier.shape=NA, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="Mix Category")+
	ylab(NULL)+
	theme_classic()+scale_color_discrete(name="Mix Ratio")

p1<-ggplot(fpH,aes(nHQreads,FPR_relToHQMapReads))+
	geom_point(aes(colour=Mix.type),show.legend=FALSE)+
	labs(y="FPR relative to mapped reads",x="Read Depth")+
	theme_classic()

ggarrange(p1,p2,ncol=2,labels=c('a','b'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))

#is FPrate related to read depth?
summary(lm(FPR_relToHQMapReads~nHQreads,fpH)) #no

#is FPrate significantly higher in pools vs pures?

pool<-subset(fpH, MixBinary=="pool")
summary(pool)
pure<-subset(fpH, MixBinary=="single")
summary(pure)

t.test(FPR_relToHQMapReads~MixBinary,fpH) #yes, significantly higher in pool

p1<-ggplot(fpH,aes(factor(MixBinary),FPR_Paus))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. australis")+
#	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(legend.position="none")+
	scale_y_continuous(limits=c(0,0.04))

p2<-ggplot(fpH,aes(factor(MixBinary),FPR_Pdel))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. delicatissima")+
	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(legend.position="none")+
	scale_y_continuous(limits=c(0,0.04))

p3<-ggplot(fpH,aes(factor(MixBinary),FPR_Pmul))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. multiseries")+
	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(legend.position="none")+
	scale_y_continuous(limits=c(0,0.001))

p4<-ggplot(fpH,aes(factor(MixBinary),FPR_Pmstri))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. multistriata")+
	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(legend.position="none")+
	scale_y_continuous(limits=c(0,0.001))

p5<-ggplot(fpH,aes(factor(MixBinary),FPR_Ppun))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. pungens")+
	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(legend.position="none")+
	scale_y_continuous(limits=c(0,0.04))

p6<-ggplot(fpH,aes(factor(MixBinary),FPR_Psub))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="FPR relative to mapped reads",x="P. subpacifica")+
	ylab(NULL)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))+
	scale_y_continuous(limits=c(0,0.04))+scale_color_discrete(name="Mix Ratio")


ggarrange(p1,p2,p3,p4,p5,p6,ncol=6,labels=c('a','b','c','d','e','f'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))

#is FPrate significantly higher in pools vs pures BY FOCAL SPP?

t.test(FPR_Paus~MixBinary,alternative="greater",fpH) #trend, p=0.09
t.test(FPR_Pdel~MixBinary,alternative="greater",fpH) #yes, p=3.959e-11
t.test(FPR_Pmul~MixBinary,alternative="greater",fpH) #yes, p=0.008 
t.test(FPR_Pmstri~MixBinary,alternative="greater",fpH) #yes, p=0.004
t.test(FPR_Ppun~MixBinary,alternative="greater",fpH) #yes, p=0.005
t.test(FPR_Psub~MixBinary,alternative="greater",fpH) #trend, p=0.09


################################################
acc<-read.csv("Accuracy.csv",stringsAsFactors=TRUE)
head(acc)
str(acc)
summary(acc)

accH<-acc[acc$nHQMreads>1000,] #remove samples with less than 1000 HQ reads remaining
summary(accH)

accH2<-subset(accH,FocalSpp!="Pdel") #killing Pdel for plot bc there are no intermdiate mixes

p1<-ggplot(accH2,aes(x=ExpPercOfPNCommunity,y=ObsPercOfPNcommunity,color=FocalSpp))+
	geom_point(aes(colour=FocalSpp))+
	#geom_jitter(position = position_jitter(width=.01)) +
	labs(y="Observed",x="Expected")+
	geom_abline(intercept = 0, slope = 1, linetype=2,size=0.75,colour="grey60") +
	geom_smooth(method=lm,   # Add linear regression lines
                se=FALSE,
                fullrange=FALSE)+
	theme_classic() +scale_color_manual(
    values=c("#FC8D62","#66C2A5","#8DA0CB"),
    name="Focal Sp.",
    labels=c("Paus","Ppun","Psub"))+
	theme(legend.position="none")

#what is overall accuracy?

summary(lm(ObsPercOfPNcommunity~ExpPercOfPNCommunity,accH)) #yes

#Is accuracy a function of focal spp?
	
summary(lm(ObsPercOfPNcommunity~ExpPercOfPNCommunity+FocalSpp,accH)) #yes

# Call:
# lm(formula = ObsPercOfPNcommunity ~ ExpPercOfPNCommunity + FocalSpp, 
    # data = accH)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.188263 -0.051334  0.000184  0.052269  0.196211 

# Coefficients:
                     # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.08147    0.01844   4.419 3.24e-05 ***
# ExpPercOfPNCommunity  0.94512    0.02314  40.837  < 2e-16 ***
# FocalSppPdel         -0.02683    0.04792  -0.560   0.5772    
# FocalSppPpun         -0.11182    0.02095  -5.338 9.43e-07 ***
# FocalSppPsub         -0.06057    0.02119  -2.858   0.0055 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.07624 on 76 degrees of freedom
# Multiple R-squared:   0.96,	Adjusted R-squared:  0.9579 
# F-statistic: 456.1 on 4 and 76 DF,  p-value: < 2.2e-16

spp<-subset(accH,FocalSpp=="Ppun")
summary(lm(ObsPercOfPNcommunity~ExpPercOfPNCommunity,spp))

p2<-ggplot(accH,aes(FocalSpp,Accuracy))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="Accuracy",x="Focal Species")+
	#ylab(NULL)+
	theme_classic()+scale_color_discrete(name="Mix Ratio")

p3<-ggplot(accH,aes(Mix.type,Accuracy))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=FocalSpp),position=pd)+
	labs(y="Accuracy",x="Mix Ratio")+
	#ylab(NULL)+
	theme_classic()+scale_color_manual(
    values=c("#FC8D62","#FFD92F","#66C2A5","#8DA0CB"),
    name="Focal Sp.",
    labels=c("Paus","Pdel","Ppun","Psub"))

#Is accuracy a function of focal species, strain (culture), or mix type?
summary(aov(Accuracy~FocalSpp,accH)) #no

summary(aov(Accuracy~FocalStrain,accH)) #yes
TukeyHSD(aov(Accuracy~FocalStrain,accH))

summary(aov(Accuracy~Mix.type,accH)) #yes
TukeyHSD(aov(Accuracy~Mix.type,accH))

ggplot(accH,aes(FocalStrain,Accuracy))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="Accuracy",x="Focal Species")+
	#ylab(NULL)+
	theme_classic()

p4<-ggplot(accH,aes(FocalStrain,Accuracy))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=Mix.type),position=pd)+
	labs(y="Accuracy",x="Focal Culture")+
	#ylab(NULL)+
#	scale_colour_hue(l=50) +
	theme_classic()+scale_color_discrete(name="Mix Ratio")
	
ggarrange(p1,p2,ncol=2,labels=c('a','b'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))
ggarrange(p3,ncol=1,labels=c('c'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))
ggarrange(p4,ncol=1,labels=c('d'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))




###########################################################################
##################### Now for Natural Sample Libraries via 2bRAD

nsWQ<-read.csv("NatSamLibs.csv",stringsAsFactors=TRUE)
head(nsWQ)
str(nsWQ)
summary(nsWQ)
nsWQ$Date<-as.Date(nsWQ$Date, format = "%m/%d/%Y")
nsWQ$Year<-as.factor(nsWQ$Year)


nsWQyear<-subset(nsWQ,Year=="2019")
summary(nsWQyear)

pd <- position_dodge(.2)
#adjust limits in plots below to match N monitoring weeks by year
p1<-ggplot(nsWQyear,aes(Week,Domoic.Acid..ng.mL.))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
	geom_line(aes(),position=pd)+ #colour=factor(Year)
#	scale_fill_manual(values=c("grey60"))+
#    ylim=(0,0.3)+
	labs(title="2019",y="Domoic Acid ng/mL")+
	xlab(NULL)+
	theme_classic()+scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+scale_y_continuous(limits=c(0,0.35))+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())


p2<-ggplot(nsWQyear,aes(Week,Pseudo.nitzschia.delicatissima.group..cells.L.))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
	geom_line(aes(),position=pd)+ #colour=factor(Year)
	labs(y="P. del-like cells/L")+
	xlab(NULL)+
	theme_classic()+scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+scale_y_continuous(limits=c(0,3.6e+05))+
	theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())


p3<-ggplot(nsWQyear,aes(Week,Pseudo.nitzschia.seriata.group..cells.L.))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
	geom_line(aes(),position=pd)+ #colour=factor(Year)
	labs(y="P. ser-like cells/L")+
	xlab(NULL)+
	theme_classic()+scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+scale_y_continuous(limits=c(0,6.8e+05))+
	theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())


#ggarrange(p1,p2,p3,ncol=1,labels=c('a','b','c'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))

p5<-ggplot(nsWQyear,aes(Week,PercMappedReads))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
	geom_line(aes(),position=pd)+ #colour=factor(Year)
	labs(y="% Mapped Reads")+
	xlab(NULL)+
	theme_classic()+scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+scale_y_continuous(limits=c(0,1))+
	theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())



####################Looking at spp abundance
nsWQH<-nsWQ[nsWQ$nHQreads>1000,] #remove samples with less than 1000 HQ reads remaining
str(nsWQH)
summary(nsWQH) #tossed 10 samples for insufficient read depth
names(nsWQH)

nsWQHstack=otuStack(nsWQH,count.columns=c(49:54),condition.columns=c(1:5))[1:276,]
nsWQHstack$Week<-as.integer(nsWQHstack$Week)

year<-subset(nsWQHstack,Year=="2019")

#pal<-wes_palette("Zissou1",5,type="")

brewer.pal(5,"Set2")

p4<-ggplot(year,aes(x=Week,y=count,fill = factor(otu))) + 
    geom_bar(position = "stack",stat="identity")+
    labs(y="Relative Abundance",x="Monitoring Week")+
    scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+
    theme_bw() + scale_fill_manual(
    values=c("#FC8D62","#FFD92F","#E78AC3","#A6D854","#66C2A5","#8DA0CB"),
    name="Species",
    labels=c("Paus","Pdel","Pmser","Pmstri","Ppun","Psub"))+
    theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())
   
#    scale_fill_manual(values=c(wes_palette(name = "Darjeeling1")))+#FFD92F
        

ggarrange(p1,p2,p3,p5,p4,ncol=1,labels=c('d','h','l','p','t'),label.args=list(gp=grid::gpar(fint=4,cex=1)))

#################################
#IS DA correlated with the relative abundance of Paustralis?
head(nsWQH)
ggplot(nsWQH,aes(FracPaus,Domoic.Acid..ng.mL.))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
    labs(y="Domoic Acid ng/mL", x="Relative % Abundance of Paus")+
    theme_bw()
#	geom_line(aes(),position=pd)
	
summary(lm(Domoic.Acid..ng.mL.~FracPaus,nsWQH)) #not a strong relationship

#what about just by year?
summary(aov(FracPaus~Year,nsWQH)) #no

#IS % mapping correlated with the relative abundance of Paustralis?
summary(lm(PercMappedReads~FracPaus,nsWQH)) #yes! suggests Paus is dominating community at these times?
ggplot(nsWQH,aes(FracPaus,PercMappedReads))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(Year)
    labs(y="% Mapped Reads", x="Relative % Abundance of Paus")+
    theme_bw()

#######################################################################
################################ now compare with ARISA results

arisa<-read.csv("ARISA.csv",stringsAsFactors=TRUE)
head(arisa)
str(arisa)
summary(arisa)
arisa$Date<-as.Date(arisa$Date, format = "%m/%d/%Y")
arisa$Year<-as.factor(arisa$Year)

arisaStack=otuStack(arisa,count.columns=c(2:18),condition.columns=c(1,19:22))[1:850,]
arisaStack$Week<-as.integer(arisaStack$Week)

year<-subset(arisaStack,Year=="2017")

# brewer.pal(5,"Set2")
# brewer.pal(12,"Dark2")
#pal_simpsons("springfield")(16)
#pal_futurama("planetexpress")(4)

p6<-ggplot(year,aes(x=Week,y=count,fill = factor(otu))) + 
    geom_bar(position = "stack",stat="identity")+
    labs(y="Relative Abundance",x="Monitoring Week")+
    scale_x_continuous(limits=c(0.5,12.5),breaks=c(1:12))+
    theme_bw()+ scale_fill_manual(values=c("#FC8D62", "#71D0F5FF","#FFD92F", "#8A9197FF" ,"#D2AF81FF" ,"#D5E4A2FF", "#197EC0FF", "#E78AC3","#66C2A5","#8DA0CB","#46732EFF","#71D0F5FF", "#370335FF" ,"#075149FF" ,"#C80813FF", "#91331FFF","#8A4198FF"),
    name="Species")+
	theme(legend.position="none")

p7<-ggplot(year,aes(x=Week,y=count,fill = factor(otu))) + 
    geom_bar(position = "stack",stat="identity")+
    labs(y="Relative Abundance",x="Monitoring Week")+
    scale_x_continuous(limits=c(0.5,17.5),breaks=c(1:17))+
    theme_bw()+ scale_fill_manual(values=c("#FC8D62", "#71D0F5FF","#FFD92F", "#8A9197FF" ,"#D2AF81FF" ,"#D5E4A2FF", "#197EC0FF", "#E78AC3","#66C2A5","#8DA0CB","#46732EFF","#71D0F5FF", "#370335FF" ,"#075149FF" ,"#C80813FF", "#91331FFF","#8A4198FF"),
    name="Species")+
    theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())+
	theme(legend.position="none")

p8<-ggplot(year,aes(x=Week,y=count,fill = factor(otu))) + 
    geom_bar(position = "stack",stat="identity")+
    labs(y="Relative Abundance",x="Monitoring Week")+
    scale_x_continuous(limits=c(0.5,9.5),breaks=c(1:9))+
    theme_bw()+ scale_fill_manual(values=c("#FC8D62", "#71D0F5FF","#FFD92F", "#8A9197FF" ,"#D2AF81FF" ,"#D5E4A2FF", "#197EC0FF", "#E78AC3","#66C2A5","#8DA0CB","#46732EFF","#71D0F5FF", "#370335FF" ,"#075149FF" ,"#C80813FF", "#91331FFF","#8A4198FF"),
    name="Species")+
	theme(legend.position="none")
    	  
p9<-ggplot(year,aes(x=Week,y=count,fill = factor(otu))) + 
    geom_bar(position = "stack",stat="identity")+
    labs(y="Relative Abundance",x="Monitoring Week")+
    scale_x_continuous(limits=c(0.5,18.5),breaks=c(1:18))+
    theme_bw()+ scale_fill_manual(values=c("#FC8D62", "#71D0F5FF","#FFD92F", "#8A9197FF" ,"#D2AF81FF" ,"#D5E4A2FF", "#197EC0FF", "#E78AC3","#66C2A5","#8DA0CB","#46732EFF","#71D0F5FF", "#370335FF" ,"#075149FF" ,"#C80813FF", "#91331FFF","#8A4198FF"),
    name="Species")+
    theme(axis.title.y=element_blank(),
    	  axis.text.y=element_blank(),
    	  axis.ticks.y=element_blank())

ggarrange(p6,p7,p8,p9,ncol=2,labels=c('a','b','c','d'),label.args=list(gp=grid::gpar(fint=4,cex=1)))

### now for correlations with 2bRAD data

arisaH<-arisa[arisa$nHQreads>1000,] # toss samples with insufficient 2bRAD read depth

summary(lm(FracPaus~P.australis_P.seriata_150,arisaH))# P<0.05, R2=0.096

summary(lm(FracPdel~P.delicatissima_168,arisaH)) # NS

summary(lm(FracPpun~P.pungens_142,arisaH)) # NS

summary(lm(FracPsub~P.subpacifica_196,arisaH)) # NS

summary(lm(FracPmul~P.multiseries_144,arisaH)) # NS

ggplot(arisaH,aes(x=P.australis_P.seriata_150,y=FracPaus))+
	geom_point()+
	#geom_jitter(position = position_jitter(width=.01)) +
	labs(y="P.australis 2bRAD",x="P.australis/P.seriata ARISA")+
	geom_abline(intercept=0,slope=1,linetype=2,size=0.75,colour="grey60") + 
	geom_smooth(method=lm,   # Add linear regression lines
                se=FALSE,
                fullrange=FALSE)+theme_classic()

#what if we model as simple presence/absence data? For 2bRAD, require minimim of 1% abundance to call a spp present

RA1<-arisaH[,c(2,4,9:11)] #for the ARISA data, anything above 0 becomes 1
RA1[RA1>0]<-1
RA1

RA2<-arisaH[,c(23:28)]
RA2[RA2>=0.015]<-1
RA2
RA2[RA2<0.015]<-0
RA2

RA<-cbind(RA1,RA2)

#now fishers exact test accuracy of simple presence/absence calls

test<-RA[,c(1,6)]
test$BothPres = as.factor(apply(test,1,function(x){sum(x>=1)})) 
summary(test)

Paus <-
matrix(c(17, 0, 24, 1),
       nrow = 2,
       dimnames =
       list(c("2bRAD_pres", "2bRAD_abs"),
            c("ARISA_pres", "ARISA_abs")))

Paus
fisher.test(Paus) #NS

test<-RA[,c(2,7)]
test$BothPres = as.factor(apply(test,1,function(x){sum(x>=1)})) 
summary(test)

Pdel <-
matrix(c(0, 2, 2, 38),
       nrow = 2,
       dimnames =
       list(c("2bRAD_pres", "2bRAD_abs"),
            c("ARISA_pres", "ARISA_abs")))

Pdel
fisher.test(Pdel) #NS

test<-RA[,c(3,9)]
test$BothPres = as.factor(apply(test,1,function(x){sum(x>=1)})) 
summary(test)

Pmul <-
matrix(c(0, 5, 2, 35),
       nrow = 2,
       dimnames =
       list(c("2bRAD_pres", "2bRAD_abs"),
            c("ARISA_pres", "ARISA_abs")))

Pmul
fisher.test(Pmul) #NS

test<-RA[,c(4,10)]
test$BothPres = as.factor(apply(test,1,function(x){sum(x>=1)})) 
summary(test)

Ppun <-
matrix(c(32, 2, 8, 0),
       nrow = 2,
       dimnames =
       list(c("2bRAD_pres", "2bRAD_abs"),
            c("ARISA_pres", "ARISA_abs")))

Ppun
fisher.test(Ppun) #NS

test<-RA[,c(5,11)]
test$BothPres = as.factor(apply(test,1,function(x){sum(x>=1)})) 
summary(test)

Psub <-
matrix(c(26, 3, 12, 1),
       nrow = 2,
       dimnames =
       list(c("2bRAD_pres", "2bRAD_abs"),
            c("ARISA_pres", "ARISA_abs")))

Psub
fisher.test(Psub) #NS


#################################
# Now for pop gen on allele freqs
#################################

#plot histogram of coverage 

plot_multi_histogram <- function(df, feature, label_column) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
#    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column))
}

cov<-read.csv("PpunCoverage.csv")
plot_multi_histogram(cov,'Ntags','coverage')

set <- subset(cov,coverage=="80x")
summary(set)

################ read in CMH test results

cmh<- read.csv("CMH_analysis_ByYear.csv")

head(cmh)

cmh_good<-subset(cmh,good=="yes")

cmh_good
nrow(cmh_good)

id=c(1:2) #column with species & tag designation
maf=c(19:22) #columns with MAF data to stack
cmh_good_st=stack(cmh_good[,maf]) 
cmh_good_st[,3:4]=cmh_good[,id] #dat
head(cmh_good_st)
names(cmh_good_st)=c("MAF","year","species","tag")

cmh_good_st$DA<-c(rep("low",14),rep("high",14),rep("low",14),rep("high",14))
head(cmh_good_st)


brewer.pal(5,"Set2")


cmh_good_st$DA<-factor(cmh_good_st$DA,levels=c("low","high"))
ggplot(cmh_good_st,aes(x=DA, y=MAF, color = year))+
#	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(),position=pd)+ #colour=factor(year)+
    labs(y="Minor Allele Frequency", x="Relative Domoic Acid")+
    facet_wrap(~tag,scales="fixed",ncol=5)+theme_bw()+ scale_color_manual(
    values=c("#66C2A5","#FC8D62","#A6D854","#E78AC3"),
    name="Year",
    labels=c("2015","2017","2018","2019"))
