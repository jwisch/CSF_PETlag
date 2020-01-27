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

###################################################################################################
#Data Cleaning#####################################################################################
###################################################################################################

CSF<-CSF[CSF$ELECSYS_BATCH == 1500,]
CSF<-CSF[,c("MAP_ID", "LP_date", "E_tau", "E_ptau", "E_ab42")]
CSF$LP_date<-as.Date(CSF$LP_date, format = "%m/%d/%Y")
colnames(CSF)[1]<-"ID"
CSF<-CSF[complete.cases(CSF$E_ab42),]

TAU<-CleanPET(paste(FILEPATH_DATA, "HASD_ACS_DR15_TAU.csv", sep = ""), "Tauopathy", "TauPos", 1.22)

#Since we're doing everything by amyloid status, have to pull this data, too
Amyloid<-CleanPET(paste(FILEPATH_DATA, "HASD_ACS_DR15_PIB.csv", sep = ""), "PIB_fSUVR_rsf_TOT_CORTMEAN", "PIBpos", 1.42)
Amyloid<-Amyloid[,c("ID", "PET_Date", "PIBpos")]
colnames(Amyloid)[2]<-"PIB_Date"
TAU<-MatchbyNearestDate(TAU, Amyloid, "ID", "PET_Date", "PIB_Date")


df<-merge(TAU, CSF, by = "ID", all = FALSE)
df$TimeBetween<-as.numeric( df$LP_date - df$PET_Date)/365

df$TimeCut<-cut(df$TimeBetween, breaks = c( -13, -10, -8, -6, -4, -2, 0, 2))
#Now making sure we only have 1 visit per individual in a given timecut
df<-df[with(df, order(ID, TimeCut)),]
df<-df [!duplicated(df[c("ID", "TimeCut")]),]

Tau_Regions<-df[,c(1, 7, 807:831, 868:880, 916:925, 26, 1066, 1068:1074)]
df<-df[,c("ID", "Group.1", "PET_Date", "Tauopathy", "TauPos", "LP_date", "E_tau", "E_ptau",
          "E_ab42", "TimeBetween", "TimeCut", "PIBpos")]
TAU<-TAU[,c("ID", "PET_Date", "Tauopathy", "TauPos", "PIBpos")]

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
demogs$AgeatPET<-as.numeric(demogs$PET_Date - demogs$BIRTH)/365

CDR<-CDR[,c("ID", "CDR", "TESTDATE")]
CDR$TESTDATE<-as.Date(CDR$TESTDATE, format = "%Y-%m-%d")

demogs<-MatchbyNearestDate(demogs, CDR, "ID", "PET_Date", "TESTDATE")
demogs$apoe4<-ifelse(demogs$apoe == "34" | demogs$apoe == "44" | demogs$apoe == "24", 1, 0)


myVars <- c("GENDER", "EDUC", "apoe4", "race2", "AgeatLP", "AgeatPET", "TimeBetween", "CDR", "TauPos", "PIBpos")
catVars <- c("apoe4",  "race2", "GENDER", "CDR", "TauPos", "PIBpos")
CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))
write.csv(print(CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))),
          paste(FILEPATH_DATA, "DemogsTable_Tau.csv", sep = ""), row.names = TRUE)
###################################################################################################
###################################################################################################




###################################################################################################
#Making Figures
#######################################################################################################

p1_all<-gplot(df[df$TimeCut == "(-6,-4]",], "E_ptau", "Tauopathy", "Log Transformed pTau",
          "Log Transformed Tauopathy",c(-0.2, 0.8), c(2.0, 5.0), 0.6, 4.8, 4.5) + ggtitle("B. Correlations associated with the (-6, -4] Bin")

p1_pos<-gplot(df[df$TimeCut == "(-8,-6]"& df$PIBpos == 1,], "E_ptau", "Tauopathy", "Log Transformed pTau",
              "Log Transformed Tauopathy",c(-0.2, 0.8), c(2.0, 5.0), 0.6, 4.8, 4.5)+ ggtitle("F. Correlations associated with the (-8, -6] Bin")

p1_neg<-gplot(df[df$TimeCut == "(-6,-4]"& df$PIBpos == 0,], "E_ptau", "Tauopathy", "Log Transformed pTau",
              "Log Transformed Tauopathy",c(-0.2, 0.8), c(2.0, 5.0), 0.6, 4.8, 4.5)+ ggtitle("D. Correlations associated with the (-6, -4] Bin")

#Doing Fisher R to Z to compare the correlations
#https://onlinelibrary.wiley.com/doi/full/10.1002/0471667196.ess5050.pub2


ForPlot<-data.frame(rbind(CreatePlotDF(df, "(-10,-8]", "Tauopathy", "E_ptau"), 
                          CreatePlotDF(df, "(-8,-6]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df, "(-6,-4]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df, "(-4,-2]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df, "(-2,0]", "Tauopathy", "E_ptau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(1, 5, 4, 3, 2)])
p_overall<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Tau - CSF pTau \nCorrelation")+ ggtitle("A. All Participants")+ ylim(c(-0.5, 1))+ theme_classic()



ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 1,], "(-10,-8]", "Tauopathy", "E_ptau"), 
                          CreatePlotDF(df[df$PIBpos == 1,], "(-8,-6]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-6,-4]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-4,-2]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-2,0]", "Tauopathy", "E_ptau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(1, 5, 4, 3, 2)])
p_pos<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Tau - CSF pTau \nCorrelation")+ ggtitle("E. Amyloid Positive Individuals")+ ylim(c(-0.5, 1))+ theme_classic()



ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 0,], "(-10,-8]", "Tauopathy", "E_ptau"), 
                          CreatePlotDF(df[df$PIBpos == 0,], "(-8,-6]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-6,-4]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-4,-2]", "Tauopathy", "E_ptau"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-2,0]", "Tauopathy", "E_ptau")))

colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(1, 5, 4, 3, 2)])
p_neg<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Tau - CSF pTau \nCorrelation")+ ggtitle("C. Amyloid Negative Individuals")+ ylim(c(-0.5, 1))+ theme_classic()



lay = rbind(c(1, 2, 2, 3), c(1, 2, 2, 3), c(4, 5, 5, 6))
grid.arrange(p_neg, p_overall, p_pos, p1_neg, p1_all, p1_pos, layout_matrix = lay)
grid.arrange(p_overall, p1_all, nrow = 1)
grid.arrange(p_neg, p_pos, p1_neg, p1_pos, nrow = 2)

lay = rbind(c(1, 1, 4), c(2, 2, 5), c(3, 3, 6))
grid.arrange(p_overall,p_neg,  p_pos, p1_all, p1_neg,  p1_pos, layout_matrix = lay)

lay = rbind(c(1, 4), c(2, 5), c(3, 6))

rm(p_pos, p_overall, p1_all, p1_pos, lay)
###################################################################################################
#######################################################################################################################
#Doing the cutpoint analysis
###################################################################################################
###################################################################################################
library(cutpointr)

GetCutpoint<-function(TIMECUT){
cp <- cutpointr(df[df$TimeCut == TIMECUT,], E_ptau, TauPos, 
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
Cutpoints$TimeCut = factor(Cutpoints$TimeCut,levels(Cutpoints$TimeCut)[c(4, 3, 2, 1, 5, 6)])
for(i in 2:4){Cutpoints[,i]<-as.numeric(as.character(Cutpoints[,i]))}
Cutpoints$TimeCut = factor(Cutpoints$TimeCut,levels(Cutpoints$TimeCut)[c(3, 4, 6, 5, 1, 2)])
Cutpoints$TimeCut<-as.character(Cutpoints$TimeCut)
Cutpoints[7,1]<-"(0,2]"
Cutpoints$TimeCut<-as.factor(Cutpoints$TimeCut)
Cutpoints$TimeCut = factor(Cutpoints$TimeCut,levels(Cutpoints$TimeCut)[c(2, 1, 6, 5, 4, 3, 7)])

library(ggplot2)
ggplot(Cutpoints, aes(x = TimeCut, y = Mean)) +  geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(Bottom5), ymax=(Top5), width=.1))+
  xlab("Time Window") + ylab("AUC")+ ggtitle("Predictive Power for Tau Positivity \nConsidered by Time Lag")+ 
  ylim(c(0.2, 1.0)) + theme_classic()

cp0 <- cutpointr(df[df$TimeCut == "(-6,-4]",], E_ptau, TauPos, 
                method = maximize_boot_metric, metric = youden, boot_runs = 1000)
plot(cp0, theme_classic())

###################################################################################################
###################################################################################################
  

###################################################################################################
###################################################################################################
#################Looking at the regional progression of correlations 


RegionalResult<-data.frame("Region" = NA,
                           "Correlation" = NA,
                           "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
k<-0
for(i in 3:50){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(Tau_Regions, i, TimeCutList[j], "E_ptau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult[k, 1]<-names(Tau_Regions[i])
    RegionalResult[k, 2]<-fisherz(hold$estimate)
    RegionalResult[k, 3]<-hold$p.value
    RegionalResult[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("TAU")

RegionalResult_plot<-RegionalResult[RegionalResult$Region %in% Regions,]


b <- c(-0.7, 0, 1.3)
MakeBrainPlot("(-13,-10]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(-10,-8]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(-8,-6]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(-6,-4]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(-4,-2]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(-2,0]", b, -0.7, 1.3, RegionalResult_plot)
MakeBrainPlot("(0,2]", b, -0.7, 1.3, RegionalResult_plot)




#Now repeating for amyloid positives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                           "Correlation" = NA,
                           "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1,7)]
k<-0
for(i in 3:50){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(Tau_Regions[Tau_Regions$PIBpos == 1,], i, TimeCutList[j], "E_ptau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(Tau_Regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("TAU")
RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]

MakeBrainPlot("(-8,-6]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(2,4]", b, -0.7, 1.3,  RegionalResult_plot_PIBpos)


#Now repeating for amyloid negatives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                                  "Correlation" = NA,
                                  "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1,7)]
k<-0
for(i in 3:50){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(Tau_Regions[Tau_Regions$PIBpos == 0,], i, TimeCutList[j], "E_ptau")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(Tau_Regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("TAU")
RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]
MakeBrainPlot("(-10,-8]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-8,-6]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -0.7, 1.3, RegionalResult_plot_PIBpos)

