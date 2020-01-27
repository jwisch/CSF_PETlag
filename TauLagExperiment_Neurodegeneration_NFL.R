library(ggplot2)
library(dplyr)
library(ggseg)
library(gridExtra)
library(grid)
library(psych)
library(tableone)

FILEPATH_CODE<-"C:/Users/julie.wisch/Documents/Transition/InProgress_CSFPETlag/Code/"
FILEPATH_DATA<-"C:/Users/julie.wisch/Documents/Transition/InProgress_CSFPETlag/Data/"

source(paste(FILEPATH_CODE, "CommonFuncs.R", sep = ""))
source(paste(FILEPATH_CODE, "TauLagExperiment_Functions.R", sep = ""))


CSF<-read.csv(paste(FILEPATH_DATA,"OI_Schindler_FullNFL_2019_06_10.csv", sep = ""))
demogs<-read.csv(paste(FILEPATH_DATA,"DR_demographics_20190122.csv", sep = ""))
CDR<-read.csv(paste(FILEPATH_DATA, "DR_clinical_20190122.csv", sep = ""))
MRI<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR15/HASD_ACS_DR15_3TMR.csv")

###################################################################################################
#Data Cleaning#####################################################################################
###################################################################################################

CSF<-CSF[CSF$ELECSYS_BATCH == 1500,]
CSF<-CSF[,c("MAP_ID", "LP_date", "NFL", "E_ptau", "E_ab42")]
CSF$LP_date<-as.Date(CSF$LP_date, format = "%m/%d/%Y")
colnames(CSF)[1]<-"ID"
CSF<-CSF[complete.cases(CSF$NFL),]
CSF$E_tau<-log(CSF$NFL)

MRI$MR_Date<-as.Date(MRI$MR_Date, format= "%m/%d/%Y")
colnames(MRI)[5]<-"ID"
MRI<-MRI[complete.cases(MRI$ID),]
MRI<-MRI[complete.cases(MRI[,"MR_TOTV_INTRACRANIAL"]),]
MRI<-MRI[with(MRI, order(ID, MR_Date)),]
MRI<-MRI[!(MRI$Scanner == "CCIR Vida"),]
MRI<-aggregate(MRI, list(MRI$ID), FUN = tail, 1)
MRI_regions<-MRI
MRI<-MRI[,c("ID", "MR_Date", "MR_LT_INFRTMP", "MR_LT_MIDTMP", "MR_LT_SUPERTMP", "MR_RT_INFRTMP", "MR_RT_MIDTMP",
            "MR_RT_SUPERTMP", "MR_LT_INFRPRTL", "MR_LT_SUPERPRTL", "MR_RT_INFRPRTL", "MR_RT_SUPERPRTL",
            "MR_LT_ENTORHINAL", "MR_RT_ENTORHINAL", "MR_LT_PRECUNEUS", "MR_RT_PRECUNEUS",
            "MR_LV_HIPPOCAMPUS", "MR_RV_HIPPOCAMPUS", "MR_TOTV_INTRACRANIAL")]
MRI<-MRI[MRI$ID %in% CSF$ID,]
ScaleandNormalize<-function(COLUMN, COLNAME){
  MEAN<-mean(MRI$MR_TOTV_INTRACRANIAL)
  model<-lm(COLUMN ~ MR_TOTV_INTRACRANIAL, data = MRI)
  MRI[,paste("Normalized", COLNAME, sep = "")]<-scale(COLUMN - (coef(model)[2] * (MRI$MR_TOTV_INTRACRANIAL - MEAN)))
  return(MRI)}

for(i in 17:18){
  MRI<-ScaleandNormalize(MRI[,i], names(MRI[i]))
}


MRI$ADsig<-(rowSums(scale(MRI[,3:16])))/14
MRI$HippoVol<-(MRI$NormalizedMR_LV_HIPPOCAMPUS + MRI$NormalizedMR_RV_HIPPOCAMPUS )/2
MRI<-MRI[,c("ID", "MR_Date", "ADsig", "HippoVol")]


#Since we're doing everything by amyloid status, have to pull this data, too
Amyloid<-CleanPET(paste(FILEPATH_DATA, "HASD_ACS_DR15_PIB.csv", sep = ""), "PIB_fSUVR_rsf_TOT_CORTMEAN", "PIBpos", 1.42)
Amyloid<-Amyloid[,c("ID", "PET_Date", "PIBpos")]
colnames(Amyloid)[2]<-"PIB_Date"
MRI<-MatchbyNearestDate(MRI, Amyloid, "ID", "MR_Date", "PIB_Date")


df<-merge(MRI, CSF, by = "ID", all = FALSE)
df$TimeBetween<-as.numeric( df$LP_date - df$MR_Date)/365

df$TimeCut<-cut(df$TimeBetween, breaks = c( -13, -10, -8, -6, -4, -2, 0, 2))
#Now making sure we only have 1 visit per individual in a given timecut
df<-df[with(df, order(ID, TimeCut)),]
df<-df [!duplicated(df[c("ID", "TimeCut")]),]


###################################################################################################
###################################################################################################



###################################################################################################
###################################################################################################
#Getting Demographics info for this group
###################################################################################################
###################################################################################################

demogs<-merge(demogs, df, by = "ID")
demogs$BIRTH<-format(as.Date(demogs$BIRTH, "%d-%b-%y"), "19%y-%m-%d")
demogs$BIRTH<-as.Date(demogs$BIRTH, format = "%Y-%m-%d")
demogs$AgeatLP<-as.numeric(demogs$LP_date - demogs$BIRTH)/365
demogs$AgeatPET<-as.numeric(demogs$MR_Date - demogs$BIRTH)/365

CDR<-CDR[,c("ID", "CDR", "TESTDATE")]
CDR$TESTDATE<-as.Date(CDR$TESTDATE, format = "%Y-%m-%d")

demogs<-MatchbyNearestDate(demogs, CDR, "ID", "MR_Date", "TESTDATE")
demogs$apoe4<-ifelse(demogs$apoe == "34" | demogs$apoe == "44" | demogs$apoe == "24", 1, 0)


myVars <- c("GENDER", "EDUC", "apoe4", "race2", "AgeatLP", "AgeatPET", "TimeBetween", "CDR", "TauPos", "PIBpos")
catVars <- c("apoe4",  "race2", "GENDER", "CDR", "TauPos", "PIBpos")
CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))
write.csv(print(CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))),
          paste(FILEPATH_DATA, "DemogsTable_Neurodegen.csv", sep = ""), row.names = TRUE)
###################################################################################################
###################################################################################################




###################################################################################################
#Making Figures
#######################################################################################################

p1_all<-gplot2(df[df$TimeCut == "(-4,-2]",], "E_tau", "ADsig", "Log Transformed NFL",
          "Normalized AD Signature Cortical Thickness", c(-2, 2),c(0.58, .85), 1.4, 0.85, 0.80)

p1_pos<-gplot2(df[df$TimeCut == "(0,2]"& df$PIBpos == 1,], "E_tau", "ADsig", "Log Transformed NFL",
              "Normalized AD Signature Cortical Thickness",c(-2, 2),c(0.58, .85), 1.4, 0.85, 0.80)

p1_neg<-gplot2(df[df$TimeCut == "(-10,-8]"& df$PIBpos == 0,], "E_tau", "ADsig", "Log Transformed NFL",
               "Normalized AD Signature Cortical Thickness",c(-2, 2),c(0.58, .85), 1.4, 0.85, 0.80)

#Doing Fisher R to Z to compare the correlations
#https://onlinelibrary.wiley.com/doi/full/10.1002/0471667196.ess5050.pub2


ForPlot<-data.frame(rbind(CreatePlotDF(df, "(-13,-10]", "ADsig", "E_tau"), 
                          CreatePlotDF(df, "(-10,-8]", "ADsig", "E_tau"), 
                          CreatePlotDF(df, "(-8,-6]", "ADsig", "E_tau"),
                          CreatePlotDF(df, "(-6,-4]", "ADsig", "E_tau"),
                          CreatePlotDF(df, "(-4,-2]", "ADsig", "E_tau"),
                          CreatePlotDF(df, "(-2,0]", "ADsig", "E_tau"),
                          CreatePlotDF(df, "(0,2]", "ADsig", "E_tau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(2, 1, 6, 5, 4, 3, 7)])
p_overall<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("AD Cortical Thickness - CSF NFL Correlation")+ ggtitle("All Participants")+ ylim(c(-0.65, 0.65))+ theme_classic()



ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 1,], "(-10,-8]", "ADsig", "E_tau"), 
                          CreatePlotDF(df[df$PIBpos == 1,], "(-8,-6]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-6,-4]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-4,-2]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-2,0]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(0,2]", "ADsig", "E_tau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(1, 5, 4, 3, 2, 6)])
p_pos<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("AD Cortical Thickness - CSF NFL Correlation")+ ggtitle("Amyloid Positive Individuals")+ ylim(c(-0.65, 0.65))+ theme_classic()


ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 0,], "(-10,-8]", "ADsig", "E_tau"), 
                          CreatePlotDF(df[df$PIBpos == 0,], "(-8,-6]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-6,-4]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-4,-2]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-2,0]", "ADsig", "E_tau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(0,2]", "ADsig", "E_tau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(1, 5, 4, 3, 2, 6)])
p_neg<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("AD Cortical Thickness - CSF NFL Correlation")+ ggtitle("Amyloid Negative Individuals")+ ylim(c(-0.65, 0.65))+ theme_classic()



lay = rbind(c(1, 2, 2, 3), c(1, 2, 2, 3), c(4, 5, 5, 6))
grid.arrange(p_neg, p_overall, p_pos, p1_neg, p1_all, p1_pos, layout_matrix = lay)
grid.arrange(p_overall, p1_all, nrow = 1)
grid.arrange(p_neg, p_pos, p1_neg, p1_pos, nrow = 2)

rm(p_pos, p_overall, p1_all, p1_pos, lay)
###################################################################################################
#######################################################################################################################
#Doing the cutpoint analysis
###################################################################################################
###################################################################################################
library(cutpointr)
df$MRIPos<-ifelse(df$ADsig < 0, 1, 0)
GetCutpoint<-function(TIMECUT){
cp <- cutpointr(df[df$TimeCut == TIMECUT,], E_tau, MRIPos, 
                 method = maximize_boot_metric, metric = youden, boot_runs = 1000)
result<-c("TimeCut" = TIMECUT, "Mean" = as.data.frame(summary(cp)[[8]])[2,6],
          "Bottom5" = as.data.frame(summary(cp)[[8]])[2,3],
          "Top5" = as.data.frame(summary(cp)[[8]])[2,8])
return(result)}


Cutpoints<-data.frame(rbind(GetCutpoint("(-13,-10]"),
                            GetCutpoint("(-10,-8]"),
                            GetCutpoint("(-8,-6]"),
                            GetCutpoint("(-6,-4]"),
                            GetCutpoint("(-4,-2]"),
                            GetCutpoint("(-2,0]"),
                            GetCutpoint("(0,2]")))
Cutpoints$TimeCut = factor(Cutpoints$TimeCut,levels(Cutpoints$TimeCut)[c(2, 1, 6, 5, 4, 3, 7)])
for(i in 2:4){Cutpoints[,i]<-as.numeric(as.character(Cutpoints[,i]))}

library(ggplot2)
ggplot(Cutpoints, aes(x = TimeCut, y = Mean)) +  geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(Bottom5), ymax=(Top5), width=.1))+
  xlab("Time Window") + ylab("AUC")+ ggtitle("Predictive Power for Neurodegeneration Identification \nConsidered by Time Lag")+ 
  ylim(c(0.1, 1.0)) + theme_classic()

cp0 <- cutpointr(df[df$TimeCut == "(0,2]",], E_tau, MRIPos, 
                method = maximize_boot_metric, metric = youden, boot_runs = 1000)
plot(cp0, theme_classic())

###################################################################################################
###################################################################################################
  

###################################################################################################
###################################################################################################
#################Looking at the regional progression of correlations 

ScaleandNormalize2<-function(COLUMN, COLNAME){
  MEAN<-mean(MRI_regions$MR_TOTV_INTRACRANIAL)
  model<-lm(COLUMN ~ MR_TOTV_INTRACRANIAL, data = MRI_regions)
  MRI_regions[,paste("Normalized", COLNAME, sep = "")]<-scale(COLUMN - (coef(model)[2] * (MRI_regions$MR_TOTV_INTRACRANIAL - MEAN)))
  return(MRI_regions)}
for(i in 50:135){
  MRI_regions<-ScaleandNormalize2(MRI_regions[,i], names(MRI_regions[i]))
}



MRI_regions<-MRI_regions[,c(6:7, 346:431)]
  NAMES<-names(MRI_regions)[rep(c(TRUE,FALSE),44)]
  NAMES<-NAMES[-1]

#averaging every pair of columns together so i have a single volume for each
MRI_regions <- data.frame(cbind(MRI_regions[,1:2], byapply(MRI_regions[,3:88], 2, rowMeans)))
colnames(MRI_regions)[3:45]<-NAMES


MRI_regions<-merge(MRI_regions, df, by = c("ID", "MR_Date"), all = FALSE)

RegionalResult<-data.frame("Region" = NA,
                           "Correlation" = NA,
                           "P" = NA, "TimeCut" = NA)


TimeCutList<-levels(df$TimeCut)
k<-0
for(i in 3:45){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(MRI_regions, i, TimeCutList[j], "E_tau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult[k, 1]<-names(MRI_regions[i])
    RegionalResult[k, 2]<-fisherz(hold$estimate)
    RegionalResult[k, 3]<-hold$p.value
    RegionalResult[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("NormalizedMR_LV")
Regions<-substr(Regions, start = 35, stop = nchar(Regions))
RegionalResult$Region<-substr(RegionalResult$Region, start = 17, stop = nchar(RegionalResult$Region))

RegionalResult_plot<-RegionalResult[RegionalResult$Region %in% Regions,]


b <- c(-0.6, 0, 0.9)
MakeBrainPlot("(-13,-10]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(-10,-8]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(-8,-6]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(-6,-4]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(-4,-2]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(-2,0]", b, -0.6, 0.9, RegionalResult_plot)
MakeBrainPlot("(0,2]", b, -0.6, 0.9, RegionalResult_plot)




#Now repeating for amyloid positives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                           "Correlation" = NA,
                           "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1)]
k<-0

for(i in 3:45){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(MRI_regions[MRI_regions$PIBpos == 1,], i, TimeCutList[j], "E_tau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(MRI_regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


RegionalResult_PIBpos$Region<-substr(RegionalResult_PIBpos$Region, start = 17, stop = nchar(RegionalResult_PIBpos$Region))

RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]



MakeBrainPlot("(-10,-8]", b, -0.6, 0.6, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-8,-6]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)



#Amyloid Negative
#Now repeating for amyloid positives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                                  "Correlation" = NA,
                                  "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1)]
k<-0
R_COR<-function(DF, index, Window, CSF){cor.test(DF[DF[,"TimeCut"] == Window,index], DF[DF[,"TimeCut"] == Window,CSF], method = "spearman")}

for(i in 3:24){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(MRI_regions[MRI_regions$PIBpos == 0,], i, TimeCutList[j], "E_tau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(MRI_regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


RegionalResult_PIBpos$Region<-substr(RegionalResult_PIBpos$Region, start = 17, stop = nchar(RegionalResult_PIBpos$Region))

RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]



MakeBrainPlot("(-10,-8]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-8,-6]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -0.6, 0.9, RegionalResult_plot_PIBpos)


