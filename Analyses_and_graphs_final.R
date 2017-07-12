#### Clear Workspace, Load packages ####
rm(list = ls())

library(lme4)
library(lmerTest) # will enhance the lme4 summary by df's and p-values
library(pbkrtest) # required by lmerTest for Kenward-Roger's
library(ggplot2)
library(GGally)
library(plyr)
library(lsmeans)
library(gridExtra)
library(Rmisc)
library(psychometric)
library(boot)
library(ICC)
library(car)
library(grid)
library(lsr)

#### USER INPUT REQUIRED: PATH FOR df
# Enter folder where Data_Carry_Public.dat is located HERE:
dfFolder="/Users/matthiaszunhammer/Dropbox/CarryFollow/1_Paper/SciTranslMed/Revision/"
setwd(dfFolder)

#### Define helper-functions for bootstrapping ####
# R's bootstrapping (requires functions with two arguments: data(x) and cases (d))
resamples=1000
# Bootstrapped mean and confidence intervals for plots
bmean <- function(x, d) {
  return(mean(x[d],na.rm=TRUE))
}
booty_lo <- function(y,d) {
  if(length(unique(y)) > 1) {
    b <- boot(data= y[d], statistic= bmean, R= resamples)
    bci <- boot.ci(b,conf = 0.95,type="bca")
    ci<-bci$bca[4]
    }
  else 
  {ci<- bmean(y)
  }
  return(ci)
}
booty_hi <- function(y,d) {
  if(length(unique(y)) > 1) {
    b <- boot(data= y[d], statistic= bmean, R= resamples)
    bci <- boot.ci(b,conf = 0.95,type="bca")
    ci<-bci$bca[5]
    }
  else 
  {ci<- bmean(y)
  }
  return(ci)
}

#### Define helper-function for residual analysis####
residualAnalysis <- function(lmod,df,varname) {
  df$resids=resid(lmod,na.action =na.exclude) # GET Residuals
  df$zresids=scale(df$resids,center=TRUE,scale=TRUE) # Scale Residuals
  df$predicts=predict(lmod,na.action =na.exclude) # GET Predicted values
  df$zpredicts=scale(df$predicts,center=TRUE,scale=TRUE) # Scale Residuals
  plot(df[,varname],df$predicts,
       xlab=varname,ylab=paste("Predicted", varname, sep = " ")) # Compare ratings and predicted
  abline(lm(df$predicts~df[,varname]))
  plot(df$zpredicts,df$zresids,
       xlab=paste("Predicted", varname, "(z-Score)", sep = " "),ylab="Standardized residuals") # Compare standardised predicted values and standardised residuals
  abline(0,0)
  hist(df$zresids,prob=TRUE) # Histogram of standardized residuals
  curve(dnorm(x, mean=0, sd=1), add=TRUE)
  
  #qqplot: plots actual qq-plot of residuals (in black) 
  # + a line running through quartiles (.25th & .75th quantile) in black
  # + the ideal diagonal (in red)
  # + 100 random draws from the normal distribution (in grey) which further allow to eyeball how extreme the deviation from normality is
  par(pty="s")
  for (i in 1:100){
    rsample=rnorm(length(df$zresids), mean = 0, sd = 1)
    qqnorm(rsample,yaxt='n',xaxt='n',ylab="",xlab="",xlim=range(-4,4),ylim=range(-4,4),main=NULL,col='grey')
    #qqline(rsample,yaxt='n',xaxt='n',ylab=NULL,xlab=NULL,main=NULL,col='grey')
    par(new=TRUE)
  }
  qqnorm(df$zresids,xlim=range(-4,4),ylim=range(-4,4))
  qqline(df$zresids)
  abline(0,1,col='red')
}

#### Basic settings for all graphs (colours, fonts) ####
groupcolor=c(rgb(1,.580,.902),rgb(.627,.706,1),rgb(.780,.275,.910),rgb(.455,.392,.910))
grouplinecolor=c(rgb(1,.580,.902),rgb(.627,.706,1),rgb(.780,.275,.910),rgb(.455,.392,.910))
grouplabels=c("Positive experience\nNo change in RoA",
              "Negative experience\nNo change in RoA",
              "Positive experience\nChange in RoA",
              "Negative experience\nChange in RoA")

groupcolor_bw=c(rgb(.8,.8,.8),rgb(.2,.2,.2),rgb(.8,.8,.8),rgb(.2,.2,.2))
grouplinecolor_bw=c(rgb(.0,.0,.0),rgb(.0,.0,.0),rgb(.0,.0,.0),rgb(.0,.0,.0))

titleformat=element_text(size = 9,face="bold",family='Myriad Pro')
legendtitleformat = element_text(size = 7,face="bold",family='Myriad Pro')
legendformat = element_text(size = 7,family='Myriad Pro')
axistitleformat=element_text(size=7,family='Myriad Pro')
axistextformat=element_text(size=7,family='Myriad Pro')
facetlabelformat = element_text(size = 7,face="bold",family='Myriad Pro')
  
#### Load original data, define variables, outliers#####
# Load data
df <- read.table("Data_Carry_Publication.dat",header=TRUE,sep="\t")
df <- within(df, {
                id       <- factor(id)
                positive  <- factor(positive, levels = 1:0, labels = c("Exp+","Exp-"))
                RoAchange       <- factor(RoAchange, levels = 0:1, labels = c("RoA-", "RoA+"))
                day       <- factor(day, levels = 1:3, labels = c("Day 1", "Day 2","Day 3"),ordered=TRUE)
                treatsite      <- factor(treatsite, levels = 0:1, labels = c("control", "treat"))
                male      <- factor(male, levels = 0:1, labels = c("female", "male"))
                })
df$group<-factor(df$positive:df$RoAchange)
# IMPORTANT! Replace NaN's by NA... otherwise problems.
df[is.na(df)]=NA 
# Exclude Outliers.
df=df[!df$id=="O_102",] #Excluded: Temperature data for day 3 were entered ambiguously/illegibly.
df=df[!df$id=="O_103",] #Excluded: Temperature data for day 3 were entered ambiguously/illegibly..
df=df[!df$id=="S_007",] #Excluded: Temperature difference at day 3 was extremely high (+3C)
# z-Transform continuous data
df$zrating=as.numeric(scale(df$rating, center = TRUE, scale = TRUE))
df$ztemp=as.numeric(scale(df$temp, center = TRUE, scale = TRUE))
df$zexpect=as.numeric(scale(df$expect, center = TRUE, scale = TRUE))
df$zage=as.numeric(scale(df$age, center = TRUE, scale = TRUE))
#### Create data-frame contrasting contol-treat for Day 3 #####
dfDiffD3<-df[df$day=="Day 3"&df$treatsite=="treat",]
#Main outcome
dfDiffD3$rating_diff<-df$rating[df$day=="Day 3"&df$treatsite=="control"]-df$rating[df$day=="Day 3"&df$treatsite=="treat"]
#Temperature differences (potential confound to correct for)
dfDiffD3$temp_diff<-df$temp[df$day=="Day 3"&df$treatsite=="control"]-df$temp[df$day=="Day 3"&df$treatsite=="treat"]
#Conditioning strength (potential confound to correct for)
#Represents actual strength of conditioning achieved.
#Calculated as mean rating difference Control-treat across Day 1&2.
#Is HIGHLY collinear with factor "positive" since conditioning of treatment history was quite on target overal.
#However, has more power to detect effects of treatment history since it accounts for participants deviating from target despite temperature calibration.
dfDiffD3$condi_strengthd2=df$rating[df$day=="Day 2"&df$treatsite=="control"]-df$rating[df$day=="Day 2"&df$treatsite=="treat"]
dfDiffD3$condi_strengthd1=df$rating[df$day=="Day 1"&df$treatsite=="control"]-df$rating[df$day=="Day 1"&df$treatsite=="treat"]
dfDiffD3$condi_strength=rowMeans(cbind(dfDiffD3$condi_strengthd2,dfDiffD3$condi_strengthd1),na.rm = TRUE)
#z-Transform
dfDiffD3$zrating_diff<-scale(dfDiffD3$rating_diff, center = TRUE, scale = TRUE)
dfDiffD3$ztemp_diff<-scale(dfDiffD3$temp_diff, center = TRUE, scale = TRUE)
dfDiffD3$zcondi_strength<-scale(dfDiffD3$condi_strength, center = TRUE, scale = TRUE)
# Auxiliaries: Change from Day2 to Day 3 in rating_diff (separate for positive and negative groups needed for descriptives)
dfDiffD3$abs_r_change_d2d3<-dfDiffD3$rating_diff-dfDiff$rating_diff[dfDiff$day=="Day 2"]
dfDiffD3$abs_r_change_d2d3[dfDiffD3$positive=='Exp+']<-(dfDiffD3$rating_diff[dfDiffD3$positive=='Exp+']-dfDiff$rating_diff[dfDiff$day=="Day 2"&dfDiff$positive=='Exp+'])*-1

#### TABLE S1: n, age, sex for sub-studies #####
# RECRUITED N's
Nrec_simon_pre= 10 # Study 1
Nrec_simon   = 40 # Study 2
Nrec_mareile = 30 # Study 3
Nrec_oliver  = 37 # Study 4
Nrec_johanna = 47 # Study 5
Nrec_lotta   = 72 # Study 6
N_recruited=Nrec_lotta+Nrec_mareile+Nrec_johanna+Nrec_oliver+Nrec_simon+Nrec_simon_pre
# Included vs tested N's by gender
n_male=sum(df$male=="male"&df$day=="Day 3"&df$treatsite=="treat")
n_female=sum(!df$male=="male"&df$day=="Day 3"&df$treatsite=="treat")
n_male/(n_female+n_male)
# Included vs tested N's by study
N_simon_pre =length(unique(df$id[df$investigator=="Z"])) # Study 1
N_simon   =length(unique(df$id[df$investigator=="S"])) # Study 2
N_mareile =length(unique(df$id[df$investigator=="M"])) # Study 3
N_oliver  =length(unique(df$id[df$investigator=="O"])) # Study 4
N_johanna =length(unique(df$id[df$investigator=="J"])) # Study 5
N_lotta   =length(unique(df$id[df$investigator=="L"])) # Study 6
N=N_lotta+N_mareile+N_johanna+N_oliver+N_simon+N_simon_pre
# Tested N's by group/investigator
#by(df,df$investigator,length(unique))
# Tested N's by group
n_pos_noRoAchange=by(df$group[df$group=="Exp+:RoA-"&df$day=="Day 3"&df$treatsite=="treat"],
               df$investigator[df$group=="Exp+:RoA-"&df$day=="Day 3"&df$treatsite=="treat"],
               length)
n_neg_noRoAchange=by(df$group[df$group=="Exp-:RoA-"&df$day=="Day 3"&df$treatsite=="treat"],
               df$investigator[df$group=="Exp-:RoA-"&df$day=="Day 3"&df$treatsite=="treat"],
               length)
n_pos_RoAchange=by(df$group[df$group=="Exp+:RoA+"&df$day=="Day 3"&df$treatsite=="treat"],
             df$investigator[df$group=="Exp+:RoA+"&df$day=="Day 3"&df$treatsite=="treat"],
             length)
n_neg_RoAchange=by(df$group[df$group=="Exp-:RoA+"&df$day=="Day 3"&df$treatsite=="treat"],
             df$investigator[df$group=="Exp-:RoA+"&df$day=="Day 3"&df$treatsite=="treat"],
             length)
#Column sums for groups
studyN=rowSums(cbind(n_pos_noRoAchange,n_neg_noRoAchange,n_pos_RoAchange,n_neg_RoAchange),na.rm = TRUE)
studyN/sum(studyN)
groupN=colSums(cbind(n_pos_noRoAchange,n_neg_noRoAchange,n_pos_RoAchange,n_neg_RoAchange),na.rm = TRUE)
groupN/sum(groupN)
#% Male
histogram(df$male,na.rm = TRUE)
table(df$male)/length(df)
a=by(df$male,df$group,table)
lapply(a,prop.table)
table(df$male)/nrow(df)
#Mean Age (Range)
histogram(df$age,na.rm = TRUE)
mean(df$age)
range(df$age)
agemean=by(df$age,df$group,mean)
agemean
agerange=by(df$age,df$group,range)
agerange

#### MAIN RESULTS TEXT, TABLE S2 (left): EXPECTATION RATINGS TESTING DAY 3: ####
#MEAN TREATMENT EXPERIENCE Day 1
mean(df$expect[df$treatsite=='treat'&df$day=='Day 1'],na.rm = TRUE)
booty_lo(df$expect[df$treatsite=='treat'&df$day=='Day 1'])
booty_hi(df$expect[df$treatsite=='treat'&df$day=='Day 1'])
#MEAN TREATMENT EXPERIENCE Day 2 negative
mean(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp-'],na.rm = TRUE)
booty_lo(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp-'])
booty_hi(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp-'])
#MEAN TREATMENT EXPERIENCE Day 2 positive
mean(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp+'],na.rm = TRUE)
booty_lo(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp+'])
booty_hi(df$expect[df$treatsite=='treat'&df$day=='Day 2'&df$positive=='Exp+'])
#Testing changes from Day 2 to Day 3 in EXPECTATIONS
#Negative group
t.test(df$expect[df$day=="Day 2"&df$treatsite=="treat"&df$positive=='Exp-'],df$expect[df$day=="Day 3"&df$treatsite=="treat"&df$positive=='Exp-'],paired=TRUE)
#Positive group
t.test(df$expect[df$day=="Day 2"&df$treatsite=="treat"&df$positive=='Exp+'],df$expect[df$day=="Day 3"&df$treatsite=="treat"&df$positive=='Exp+'],paired=TRUE)

# Here (and in all other analyses) conditioning strength is used in lieu of factor positive to account for deviations of intended vs actual conditioning.
# Each model is calculated as a version with actual conditioning strenght (a) and per-protocol conditioning (b) strenght
# (Overall, it makes no big difference, but just to make sure we miss nothing)
lm1a<-lm(zexpect~investigator+(zcondi_strength+RoAchange)^2,
        data=dfDiffD3,na.action =na.exclude)
AIC(lm1a)
summary(lm1a)
confint(lm1a)
alias(lm1a)
Anova(lm1a,type='II') # Anova package is used to allow for Type-II ANOVA
etaSquared(lm1a,type=2,anova=TRUE)

lm1b<-lm(expect~investigator+(positive+RoAchange)^2,
        data=dfDiffD3,na.action =na.exclude)
m=lsmeans(lm1b,~positive)
postHoc_positive=show(pairs(m))
postHoc_positive$CIlo<-postHoc_positive$estimate-(postHoc_positive$SE*qt(.975, df=lm1b$df))
postHoc_positive$CIhi<-postHoc_positive$estimate+(postHoc_positive$SE*qt(.975, df=lm1b$df))
ms=show(m)
postHoc_positive
# Perform post hoc-tests on group-analysis
m=lsmeans(lm1b,~RoAchange:positive)
postHoc_RoAchange=show(pairs(m))
postHoc_RoAchange$CIlo<-postHoc_RoAchange$estimate-(postHoc_RoAchange$SE*qt(.975, df=lm1b$df))
postHoc_RoAchange$CIhi<-postHoc_RoAchange$estimate+(postHoc_RoAchange$SE*qt(.975, df=lm1b$df))
ms=show(m)
postHoc_RoAchange
# Analyze residuals
residualAnalysis(lm1a,dfDiffD3,"expect")


#### MAIN RESULTS TEXT: Descriptive rating results, t-test for day-by-day changes ####
# Descriptive: Rating differences Day 1 and Day 2 (i.e. conditioning strength)
mean(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp-'])
booty_lo(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp-'])
booty_hi(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp-'])

mean(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp+'])
booty_lo(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp+'])
booty_hi(dfDiffD3$condi_strength[dfDiffD3$positive=='Exp+'])

#ALL DAY 3
bmean(dfDiffD3$rating_diff)
booty_lo(dfDiffD3$rating_diff)
booty_hi(dfDiffD3$rating_diff)

#### MAIN RESULTS, TABLE S2 (right): PAIN RATINGS TESTING DAY 3: ####
lm2a<-lm(zrating_diff~investigator+ztemp_diff+(zcondi_strength+RoAchange)^2 # Everything up to two-way interaction effects
        ,data=dfDiffD3,na.action =na.exclude)
AIC(lm2a)
summary(lm2a)
#anova(lm4) # use for getting Type I sums of squares
confint(lm2a)
alias(lm2a)
Anova(lm2a,type='II')
etaSquared(lm2a,type=2,anova=TRUE)

lm2b<-lm(rating_diff~investigator+ztemp_diff+(positive+RoAchange)^2 # Everything up to two-way interaction effects
         ,data=dfDiffD3,na.action =na.exclude)

m=lsmeans(lm2b,~positive)
postHoc_positive=show(pairs(m))
postHoc_positive$CIlo<-postHoc_positive$estimate-(postHoc_positive$SE*qt(.975, df=lm2b$df))
postHoc_positive$CIhi<-postHoc_positive$estimate+(postHoc_positive$SE*qt(.975, df=lm2b$df))
ms=show(m)
postHoc_positive

m=lsmeans(lm2b,~RoAchange)
postHoc_RoAchange=show(pairs(m))
postHoc_RoAchange$CIlo<-postHoc_RoAchange$estimate-(postHoc_RoAchange$SE*qt(.975, df=lm2b$df))
postHoc_RoAchange$CIhi<-postHoc_RoAchange$estimate+(postHoc_RoAchange$SE*qt(.975, df=lm2b$df))
ms=show(m)
postHoc_RoAchange
residualAnalysis(lm2a,dfDiffD3,"rating_diff")


#### PLOT MAIN RESULTS, FIGURE 2 ####
# Estimate corrected single-subject values for plot 2: Simple models accounting for effects of study and temperature #
lm0.1<-lm(rating_diff~investigator+ztemp,
          data=dfDiffD3,na.action =na.exclude)
dfDiffD3$prating_diff=resid(lm0.1,na.action =na.exclude)+lm0.1$coefficients[1] # GET Residuals
lm0.2<-lm(expect~investigator,
          data=dfDiffD3,na.action =na.exclude)
dfDiffD3$pexpect=resid(lm0.2,na.action =na.exclude)+lm0.2$coefficients[1] # GET Residuals
#Expectations and pain rating differences 
#Note that prating_diff, pactual, and pexpect are used for graphing (to correct for study and temperature-differences)
# + TEMPERATURE DIFFERENCES (Expectations: STUDY only)
# OTHERWISE THEY DON'T MATCH THE STATS
plotFigure2 <- function(dfDiffD3) {
  dodge <- position_dodge(width=0.5)
  jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
  #EXPECTATIONS DAY 3
  e<- ggplot(dfDiffD3, aes(x=RoAchange:positive,y=pexpect,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
    #geom_dotplot(position=dodge, binaxis='y', stackdir='center',binwidth = 1)+ #, binwidth=1
    geom_point(stat = "summary", fun.y=bmean,position=dodge)+
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    # Scale ect.
    #scale_x_discrete(labels= c("Exp+\nRoA-", "Exp+\nRoA+", "Exp-\nRoA-", "Exp-\nRoA+"))+
    scale_colour_manual("Group:\n",values =  grouplinecolor_bw,
                        labels= grouplabels)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    scale_fill_manual("Group:\n",values =  groupcolor_bw,
                      labels= grouplabels)+ 
    scale_y_continuous(breaks=c(0, 25, 50, 75, 100),limits = c(-10, 100),expand=c(0.005, 0))+
    labs(y ="Expectation Rating")+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat,
          axis.text.x=element_blank())+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(#axis.ticks =element_blank(),
      axis.line = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x= element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'transparent'),
      panel.background = element_blank(),
      plot.title = titleformat,
      plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"),
      legend.position='none')+
    ggtitle("a) Day 3 Expectation Ratings")+
    coord_fixed(ratio=0.20)
  #PAIN RATING DIFFERENCES DAY 3
  rd<- ggplot(dfDiffD3, aes(x=RoAchange:positive,y=prating_diff,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
    #geom_dotplot(position=dodge, binaxis='y', stackdir='center',binwidth = 1)+
    geom_point(stat = "summary", fun.y=bmean,position=dodge)+
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    # Scale ect.
    #scale_x_discrete(labels= c("Exp+\nRoA-", "Exp+\nRoA+", "Exp-\nRoA-", "Exp-\nRoA+"))+
    scale_colour_manual("Group:\n",values =  grouplinecolor_bw,
                        labels= grouplabels)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    scale_fill_manual("Group:\n",values =  groupcolor_bw,
                      labels= grouplabels)+ 
    scale_y_continuous(breaks=c(0, 30, 60),limits = c(-5, 90),expand=c(0.005, 0))+
    labs(y ="Pain Rating Difference")+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat,
          axis.text.x=element_blank())+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(#axis.ticks =element_blank(),
      axis.line = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks.x = element_blank(),
      #panel.grid.major.y = element_line(colour = "black"), #element_line(colour = "#555555"), #fix grid display
      #panel.grid.minor.y = element_blank(), #fix grid display
      panel.grid.major.x = element_blank(),
      axis.title.x= element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'transparent'),
      panel.background = element_blank(),
      plot.title = titleformat,
      plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
    theme(legend.position='none',                 #place legend
          legend.direction="vertical",                    #stack format
          legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
          legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
          legend.title = legendtitleformat,
          legend.text = legendformat)+
    ggtitle("b) Day 3 Pain Ratings (Control-treat)")+
    coord_fixed(ratio=0.20)
  #lay=rbind(c(1,2),
  #          c(3,NA))
  #g=arrangeGrob(rd,a,e,rde,ae, ncol=2,nrow=3)#tc,tp,c,p,
  #g<-grid.arrange(e,rd,rde, layout_matrix=lay, widths=c(1,1.2))
  ggsave(file="Figure_2a.svg",e,width = 3.5,height = 3.5,units = c("in"),dpi = 1200,scale=1)
  ggsave(file="Figure_2b.svg",rd,width = 3.5,height = 3.5,units = c("in"),dpi = 1200,scale=1)
}
plotFigure2(dfDiffD3) # Combine and prepare figures for publication format in inkscape

#### PLOT FIGURE S3: EXPECTATIONS VS RATING DIFFERENCES DAY 3 ####
plotFigureS3 <- function(dfDiffD3) {
emin=min(dfDiffD3$pexpect, na.rm ='T')
emax=max(dfDiffD3$pexpect, na.rm ='T')
rmin=min(dfDiffD3$prating_diff, na.rm ='T')
rmax=max(dfDiffD3$prating_diff, na.rm ='T')
rde<- ggplot(dfDiffD3, aes(x=pexpect,y=prating_diff,color=RoAchange:positive,fill=RoAchange:positive))+
  geom_point()+
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x, alpha=0.5)+
  # SCALE
  scale_colour_manual("Group:\n",values = grouplinecolor,labels= grouplabels)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
  scale_fill_manual("Group:\n",values = groupcolor,labels= grouplabels)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
  scale_x_continuous(breaks=c(0, 25, 50, 75, 100),limits = c(emin, emax))+
  scale_y_continuous(breaks=c(0, 20, 40, 60, 80),limits = c(rmin, rmax))+
  labs(y ="Treatment outcome [VAS]", x = "Treatment expectation [VAS]")+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=axistextformat, #fix axis font size
        axis.title=axistitleformat)+ #fix axis title font size
  theme(strip.text.x = facetlabelformat)+
  theme(#axis.ticks =element_blank(),
    axis.line = element_line(colour = "black"),
    #axis.line.y = element_line(colour = "black"),
    #panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
    #panel.grid.minor.y = element_blank(), #fix grid display
    panel.grid.major.x = element_blank(),
    axis.title.x= axistextformat,
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'transparent'),
    panel.background = element_blank(),
#    plot.title = titleformat,
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(#legend.position=c(1.75,0.5),                 #place legend
    legend.direction="vertical",                    #stack format
    legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
    legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
    legend.title = legendtitleformat,
    legend.text = legendformat)+
#  ggtitle("Day 3 Expectation & Pain Difference")+
  coord_fixed(ratio=0.3333)
ggsave(file="Figure_S3.svg",rde,width = 12,height = 12,units = c("cm"),dpi = 1200,scale=2)
}
plotFigureS3(dfDiffD3)

#### FIGURE S2: Plot expectations, applied temperatures and ratings for all days ####
plotTempNRatings <- function(df) {
  # Select data
  dftreatsite=df[df$treatsite=='treat',] #,
  dftreatsite$rating_control<-df$rating[df$treatsite=='control']
  dftreatsite$rating_treat<-df$rating[df$treatsite=='treat']
  dftreatsite$temp_control<-df$temp[df$treatsite=='control']
  dftreatsite$temp_treat<-df$temp[df$treatsite=='treat']
  dftreatsite$diff_rating<-df$rating[df$treatsite=='control']-df$rating[df$treatsite=='treat']
  dftreatsite$diff_temp<-df$temp[df$treatsite=='control']-df$temp[df$treatsite=='treat']
  dftreatsite$mean_temp<-(df$temp[df$treatsite=='control']+df$temp[df$treatsite=='treat'])/2
  dftreatsite$expect<-df$expect[df$treatsite=='treat']  
  
  dodge <- position_dodge(width=0.5)
  jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
  
  dotalpha=0.4
  linealpha=1
  
  e<- ggplot(dftreatsite, aes(x=day,y=expect,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_line(aes(group=RoAchange:positive),stat='summary', fun.y=bmean,position=dodge,alpha=linealpha)+
    #geom_violin(adjust = 0.5,position=dodge,alpha=0.2)+ #
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    geom_dotplot(position=dodge, binwidth=1, binaxis='y', stackdir='center',alpha=dotalpha)+
    #geom_jitter(position=jdodge,alpha=0.2)+ #
    scale_colour_manual("Group:\n",values =  groupcolor,
                        labels= grouplabels)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    scale_fill_manual("Group:\n",values =  groupcolor,
                      labels= grouplabels)+
    scale_y_continuous(breaks=c(0, 20, 50, 80, 100),limits = c(0, 100))+
    labs(y ="Treatment expectations [VAS]")+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat)+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "#111111",size=0.2), #element_line(colour = "#555555"), #fix grid display
          panel.grid.minor.y = element_blank(), #fix grid display
          panel.grid.major.x = element_blank(),
          axis.title.x= element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'transparent'),
          panel.background = element_blank(),
          plot.title = titleformat,
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+ # set margins around plots to control multiplot-spacing
    theme(legend.position=c(0.5,0.5),                 #place legend
          legend.direction="vertical",                    #stack format
          legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
          legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
          legend.title = legendtitleformat,
          legend.text = legendformat)+
    ggtitle("a) Expectation Ratings")+
    coord_fixed(ratio=0.33)

  tc<- ggplot(dftreatsite, aes(x=day,y=temp_control,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_line(aes(group=RoAchange:positive),stat='summary', fun.y=bmean,position=dodge,alpha=linealpha)+
    #geom_violin(adjust = 0.5,position=dodge,alpha=0.2)+ #
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    geom_dotplot(position=dodge, binwidth=.05, binaxis='y', stackdir='center',alpha=dotalpha)+
    #geom_jitter(position=jdodge,plotFigure2=0.2)+ #
    scale_colour_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    scale_fill_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    ylim(43, 49)+
    labs(y = "Temperature [??C]")+
    #scale_y_continuous(breaks=c(44,46,48))+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat)+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "#111111",size=0.2), #element_line(colour = "#555555"), #fix grid display
          panel.grid.minor.y = element_blank(), #fix grid display
          panel.grid.major.x = element_blank(),
          axis.title.x= element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'transparent'),
          panel.background = element_blank(),
          legend.position='none',
          plot.title = titleformat,
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+#top right bottom left margins
    ggtitle("b) Temperatures Control Site")+
    coord_fixed(ratio=0.33)
  tc
  tp<- ggplot(dftreatsite, aes(x=day,y=temp_treat,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_line(aes(group=RoAchange:positive),stat='summary', fun.y=bmean,position=dodge,alpha=linealpha)+
    #geom_violin(adjust = 0.5,position=dodge,alpha=0.2)+ #
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    geom_dotplot(position=dodge, binwidth=.05, binaxis='y', stackdir='center',alpha=dotalpha)+
    #geom_jitter(position=jdodge,alpha=0.2)+ #
    scale_colour_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    scale_fill_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c("Exp-+ROA_topical_only", "Exp-+ROA_change_to_oral", "Exp++ROA_topical_only", "Exp++ROA_change_to_oral")
    ylim(43, 49)+
    labs(y = "Temperature [??C]")+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat)+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "#111111",size=0.2), #element_line(colour = "#555555"), #fix grid display
          panel.grid.minor.y = element_blank(), #fix grid display
          panel.grid.major.x = element_blank(),
          axis.title.x= element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'transparent'),
          panel.background = element_blank(),
          legend.position='none',
          plot.title = titleformat,
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+ # set margins around plots to control multiplot-spacing
    ggtitle("c) Temperatures treat Site")+ 
    coord_fixed(ratio=0.33)
  c<- ggplot(dftreatsite, aes(x=day,y=rating_control,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_line(aes(group=RoAchange:positive),stat='summary', fun.y=bmean,position=dodge,alpha=linealpha)+
    #geom_violin(adjust = 0.5,position=dodge,alpha=0.2)+ #
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    geom_dotplot(position=dodge, binwidth=1, binaxis='y', stackdir='center',alpha=dotalpha)+
    #geom_jitter(position=jdodge,alpha=0.2)+ #
    scale_colour_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c("Exp+ and RoA-", "Exp+ and RoA+", "Exp- and RoA-", "Exp- and RoA+")
    scale_fill_manual("Group:",values =  groupcolor)+   #Change Legend Nameing and Position ,labels= c( "Exp+ and RoA-", "Exp+ and RoA+","Exp- and RoA-", "Exp- and RoA+")
    labs(y ="Pain Rating [VAS]")+
    scale_y_continuous(breaks=c(0, 20, 50, 80, 100),limits = c(0, 100))+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat)+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "#111111",size=0.2), #element_line(colour = "#555555"), #fix grid display
          panel.grid.minor.y = element_blank(), #fix grid display
          panel.grid.major.x = element_blank(),
          axis.title.x= element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'transparent'),
          panel.background = element_blank(),
          legend.position='none',
          plot.title = titleformat,#theme(axis.title.y=element_blank()+
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+ # set margins around plots to control multiplot-spacing
          ggtitle("d) Pain Ratings Control Site")+
          coord_fixed(ratio=0.33)
  #scale_x_discrete(breaks=NULL)
  #ylim(0, 100)+
  p<- ggplot(dftreatsite, aes(x=day,y=rating_treat,color=RoAchange:positive,fill=RoAchange:positive))+
    geom_line(aes(group=RoAchange:positive),stat='summary', fun.y=bmean,position=dodge,alpha=linealpha)+
    #geom_violin(adjust = 0.5,position=dodge,alpha=0.2)+ #
    geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
    geom_dotplot(position=dodge, binwidth=1, binaxis='y', stackdir='center',alpha=dotalpha)+
    #geom_jitter(position=jdodge,alpha=0.2)+ #
    scale_colour_manual("Group:",values =  groupcolor)+
    scale_fill_manual("Group:",values =  groupcolor)+
    scale_y_continuous(breaks=c(0, 20, 50, 80, 100),limits = c(0, 100))+
    labs(y ="Pain Rating [VAS]")+
    theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
    theme(axis.text=axistextformat, #fix axis font size
          axis.title=axistitleformat)+ #fix axis title font size
    theme(strip.text.x = facetlabelformat)+
    theme(axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(colour = "#111111",size=0.2), #element_line(colour = "#555555"), #fix grid display
          panel.grid.minor.y = element_blank(), #fix grid display
          panel.grid.major.x = element_blank(),
          axis.title.x= element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = 'transparent'),
          panel.background = element_blank(),
          legend.position='none',
          plot.title = titleformat,
          plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+ # set margins around plots to control multiplot-spacing
    ggtitle("e) Pain Ratings treat Site")+
    coord_fixed(ratio=0.33)
#  lay=rbind(c(1,NA),
#            c(2,3),
#            c(4,5))
#  g<-grid.arrange(e,tc,tp,c,p, layout_matrix=lay, widths=c(1,1.2))
 
  ggsave(file="Figure_S2A.svg",e,width = 8,height = 8,units = c("cm"),dpi = 1200,scale=1.5)
  ggsave(file="Figure_S2B.svg",tc,width = 8,height = 8,units = c("cm"),dpi = 1200,scale=1.5)
  ggsave(file="Figure_S2C.svg",tp,width = 8,height = 8,units = c("cm"),dpi = 1200,scale=1.5)
  ggsave(file="Figure_S2D.svg",c,width = 8,height = 8,units = c("cm"),dpi = 1200,scale=1.5)
  ggsave(file="Figure_S2E.svg",p,width = 8,height = 8,units = c("cm"),dpi = 1200,scale=1.5)
}
plotTempNRatings(df)
warnings()

#### MAIN RESULTS, TABLE S3: PAIN RATINGS PREDICTED BY EXPECTATIONS TESTING DAY 3:####
lm3a<-lm(zrating_diff~investigator+ztemp_diff+zcondi_strength+RoAchange+zcondi_strength:RoAchange+zexpect # Everything up to two-way interaction effects
        ,data=dfDiffD3,na.action =na.exclude)
AIC(lm3a)
summary(lm3a)
confint(lm3a)
alias(lm3a)
#anova(lm3a) # Yields fully orthagonalized results
Anova(lm3a,type='II')
#Anova(lm3a,type='III')
etaSquared(lm3a,type=2,anova=TRUE)
residualAnalysis(lm3a,dfDiffD3,"rating_diff")

# Just to exclude that expectations have an effect after orthagonalizing with all groups:
# lm3b<-lm(zrating_diff~investigator+ztemp_diff+group+zexpect # Everything up to two-way interaction effects
#          ,data=dfDiffD3,na.action =na.exclude)
# anova(lm3b)

#### MAIN RESULTS, TABLE S4: Age and Gender:####
lm4a<-lm(zrating_diff~investigator+ztemp_diff+zcondi_strength+RoAchange+zcondi_strength:RoAchange+male+zage # Everything up to two-way interaction effects
         ,data=dfDiffD3,na.action =na.exclude)
AIC(lm4a)
summary(lm4a)
confint(lm4a)
alias(lm4a)
#anova(lm4a) # Yields fully orthagonalized results
Anova(lm4a,type='II')
#Anova(lm4a,type='III')
etaSquared(lm4a,type=2,anova=TRUE)
residualAnalysis(lm4a,dfDiffD3,"rating_diff")

#### AUXILIARY ANALYSIS: PAIN RATINGS PREDICTED BY EXPECTATIONS TESTING DAY 3, DIFFERENT EXPECTATION EFFECTS BY GROUP:####
lm3b<-lm(zrating_diff~investigator+ztemp_diff+(zcondi_strength+RoAchange+zexpect)^2 # Everything up to two-way interaction effects
         ,data=dfDiffD3,na.action =na.exclude)
AIC(lm3b)
summary(lm3b)
#anova(lm5)
confint(lm3b)
alias(lm3b)
Anova(lm3b,type='II')
etaSquared(lm3b,type=2,anova=TRUE)
residualAnalysis(lm3b,dfDiffD3,"rating_diff")


#### AUXILIARY ANALYSES: TEMPERATURE DIFFERENCES BETWEEN GROUPS####
#ALL GROUPS CONTROL SITE DAY 1 and 2
mean(df$temp[df$treatsite=='control'&df$day!='Day 3'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='control'&df$day!='Day 3'])
booty_hi(df$temp[df$treatsite=='control'&df$day!='Day 3'])

#NEGATIVE GROUPS CONTROL SITE DAY 1 and 2
mean(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp-'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp-'])
booty_hi(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp-'])
#NEGATIVE GROUPS TREATMENT SITE DAY 1 and 2
mean(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp-'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp-'])
booty_hi(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp-'])

#POSITIVE GROUPS CONTROL SITE DAY 1 and 2
mean(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp+'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp+'])
booty_hi(df$temp[df$treatsite=='control'&df$day!='Day 3'&df$positive=='Exp+'])
#POSITIVE GROUPS TREATMENT SITE DAY 1 and 2
mean(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp+'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp+'])
booty_hi(df$temp[df$treatsite=='treat'&df$day!='Day 3'&df$positive=='Exp+'])

#All GROUPs Control SITE DAY 3
mean(df$temp[df$treatsite=='control'&df$day=='Day 3'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='control'&df$day=='Day 3'])
booty_hi(df$temp[df$treatsite=='control'&df$day=='Day 3'])
#All GROUPs TREATMENT SITE DAY 3
mean(df$temp[df$treatsite=='treat'&df$day=='Day 3'],na.rm = TRUE)
booty_lo(df$temp[df$treatsite=='treat'&df$day=='Day 3'])
booty_hi(df$temp[df$treatsite=='treat'&df$day=='Day 3'])

#TEST DAY 3 TEMPERATURE DIFFERENCES
lm_temp1<-lm(temp_diff~(positive+RoAchange)^2,
             data=dfDiffD3,na.action =na.exclude)
AIC(lm_temp1)
summary(lm_temp1)
confint(lm_temp1)
anova(lm_temp1)
confint(lm_temp1)
lsmeans(lm_temp1,~positive)
lsmeans(lm_temp1,~positive:RoAchange)
etaSquared(lm_temp1,type=2,anova=TRUE)





