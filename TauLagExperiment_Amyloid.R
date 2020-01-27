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
Amyloid<-read.csv(paste(FILEPATH_DATA, "HASD_ACS_DR15_PIB.csv", sep = ""))
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


colnames(Amyloid)[5]<-"ID"
Amyloid$PET_Date<-as.Date(Amyloid$PET_Date, format= "%m/%d/%Y")
Amyloid<-Amyloid[complete.cases(Amyloid$ID),]
Amyloid<-Amyloid[complete.cases(Amyloid$PIB_fSUVR_rsf_TOT_CORTMEAN),]
Amyloid$PIBpos<-ifelse(Amyloid$PIB_fSUVR_rsf_TOT_CORTMEAN > 1.42, 1, 0)
colnames(Amyloid)[6]<-"PIB_Date"

#Keeping most recent amyloid scan only
Amyloid<-Amyloid[with(Amyloid, order(ID, PIB_Date)),]
Amyloid<-aggregate(Amyloid, list(Amyloid$ID), FUN = tail, 1)



df<-merge(Amyloid, CSF, by = "ID", all = FALSE)
df$TimeBetween<-as.numeric( df$LP_date - df$PIB_Date)/365

df$TimeCut<-cut(df$TimeBetween, breaks = c( -10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8))
#Now making sure we only have 1 visit per individual in a given timecut
df<-df[with(df, order(ID, TimeCut)),]
df<-df [!duplicated(df[c("ID", "TimeCut")]),]

Amyloid_Regions<-df[,c(1, 7, 802:826, 862:874, 910:919, 1060, 1065:1071)]
df<-df[,c("ID", "Group.1", "PIB_Date", "PIB_fSUVR_rsf_TOT_CORTMEAN", "PIBpos", "LP_date", "E_tau", "E_ptau",
          "E_ab42", "TimeBetween", "TimeCut")]
Amyloid<-Amyloid[,c("ID", "PIB_Date", "PIB_fSUVR_rsf_TOT_CORTMEAN", "PIBpos")]

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
demogs$AgeatPET<-as.numeric(demogs$PIB_Date - demogs$BIRTH)/365

CDR<-CDR[,c("ID", "CDR", "TESTDATE")]
CDR$TESTDATE<-as.Date(CDR$TESTDATE, format = "%Y-%m-%d")

demogs<-MatchbyNearestDate(demogs, CDR, "ID", "PIB_Date", "TESTDATE")
demogs$apoe4<-ifelse(demogs$apoe == "34" | demogs$apoe == "44" | demogs$apoe == "24", 1, 0)


myVars <- c("GENDER", "EDUC", "apoe4", "race2", "AgeatLP", "AgeatPET", "TimeBetween", "CDR", "PIBpos")
catVars <- c("apoe4",  "race2", "GENDER", "CDR", "PIBpos")
CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))
write.csv(print(CreateTableOne(vars = myVars, data = demogs, factorVars = catVars, strata = c("TimeCut"))),
paste(FILEPATH_DATA,"DemogsTable_Amyloid.csv", sep = ""), row.names = TRUE)
###################################################################################################
###################################################################################################



###################################################################################################
#Making Figures
#######################################################################################################
p1<-gplot(df[df$TimeCut == "(-10,-8]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p2<-gplot(df[df$TimeCut == "(-8,-6]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p3<-gplot(df[df$TimeCut == "(-6,-4]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p4<-gplot(df[df$TimeCut == "(-4,-2]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p5_neg<-gplot(df[df$TimeCut == "(-2,0]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p6<-gplot(df[df$TimeCut == "(0,2]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p7<-gplot(df[df$TimeCut == "(2,4]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p8<-gplot(df[df$TimeCut == "(4,6]",], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2,
             top = textGrob("All Participants",gp=gpar(fontsize=14)))


p1_pos<-gplot(df[df$TimeCut == "(-8,-6]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
              "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)
p2<-gplot(df[df$TimeCut == "(-6,-4]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.7)
p3<-gplot(df[df$TimeCut == "(-4,-2]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.7)
p4<-gplot(df[df$TimeCut == "(-2,0]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.7)
p5<-gplot(df[df$TimeCut == "(0,2]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.6)
p6<-gplot(df[df$TimeCut == "(2,4]"& df$PIBpos == 1,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
          "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.7)
grid.arrange(p1_pos, p2, p3, p4, p5, p6, nrow = 2,
             top = textGrob("Amyloid Positive Participants",gp=gpar(fontsize=14)))

p1_neg<-gplot(df[df$TimeCut == "(-2,0]"& df$PIBpos == 0,], "E_ab42", "PIB_fSUVR_rsf_TOT_CORTMEAN", "Log Transformed ab42",
              "Log Transformed Cortical Amyloid",c(-0.16, 1.66), c(5.0, 8.5), 1.5, 8, 7.5)

#Doing Fisher R to Z to compare the correlations
#https://onlinelibrary.wiley.com/doi/full/10.1002/0471667196.ess5050.pub2

ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 0,], "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-2,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(0,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))

ForPlot<-data.frame(rbind(CreatePlotDF(df, "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(-2,-1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(-1,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(0,1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(1,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df, "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))



colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(5, 4, 3, 2, 1, 6, 7, 8)])
p_overall<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Amyloid - CSF AB42 Correlation")+ ggtitle("All Participants")+ ylim(c(-1.35, 0.85))+ theme_classic()



ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 1,], "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-2,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(0,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))

ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 1,], "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-2,-1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(-1,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(0,1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(1,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 1,], "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))



colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(5, 4, 3, 2, 1, 6, 7, 8)])
p_pos<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Amyloid - CSF AB42 Correlation")+ ggtitle("Amyloid Positive Individuals")+ ylim(c(-0.9, 0.35))+ theme_classic()


ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 0,], "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-2,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(0,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))

ForPlot<-data.frame(rbind(CreatePlotDF(df[df$PIBpos == 0,], "(-8,-6]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-6,-4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-4,-2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-2,-1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(-1,0]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(0,1]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(1,2]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42"),
                          CreatePlotDF(df[df$PIBpos == 0,], "(2,4]", "PIB_fSUVR_rsf_TOT_CORTMEAN", "E_ab42")))


colnames(ForPlot)<-c("TimeCut", "Estimate", "CI_lower", "CI_upper")
for(i in 2:4){ForPlot[,i]<-as.numeric(as.character(ForPlot[,i]))}

ForPlot$TimeCut<-factor(ForPlot$TimeCut, levels(ForPlot$TimeCut)[c(5, 4, 3, 2, 1, 6, 7, 8)])
p_neg<-ggplot(ForPlot, aes(x = TimeCut, y = Estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(CI_lower + Estimate)/1.96, ymax=(CI_upper+Estimate)/1.96), width=.1)+
  xlab("Time Window") + ylab("PET Amyloid - CSF AB42 Correlation")+ ggtitle("Amyloid Negative Individuals")+ ylim(c(-0.9, 0.35))+ theme_classic()


lay = rbind(c(1, 2, 2, 3), c(1, 2, 2, 3), c(4, 5, 5, 6))
grid.arrange(p_neg, p_overall, p_pos, p1_neg, p5_neg, p1_pos, layout_matrix = lay)

grid.arrange(p_overall, p5_neg, nrow = 1)

grid.arrange(p_neg, p_pos, p1_neg, p1_pos, nrow = 2)
rm(p1, p2, p3, p4, p5, p6, p7, p8, p_overall, p_pos, p5_neg, p1_pos, p1_neg, p_neg, lay)
###################################################################################################
#######################################################################################################################
#Doing the cutpoint analysis
###################################################################################################
###################################################################################################
library(cutpointr)

GetCutpoint<-function(TIMECUT){
cp <- cutpointr(df[df$TimeCut == TIMECUT,], E_ab42, PIBpos, 
                 method = maximize_boot_metric, metric = youden, boot_runs = 1000)
result<-c("TimeCut" = TIMECUT, "Mean" = as.data.frame(summary(cp)[[8]])[2,6],
          "Bottom5" = as.data.frame(summary(cp)[[8]])[2,3],
          "Top5" = as.data.frame(summary(cp)[[8]])[2,8])
return(result)}


Cutpoints<-data.frame(rbind(GetCutpoint("(-8,-6]"),
                            GetCutpoint("(-6,-4]"),
                            GetCutpoint("(-4,-2]"),
                            GetCutpoint("(-2,0]"),
                            GetCutpoint("(0,2]"),
                            GetCutpoint("(2,4]")))
Cutpoints$TimeCut = factor(Cutpoints$TimeCut,levels(Cutpoints$TimeCut)[c(4, 3, 2, 1, 5, 6)])
for(i in 2:4){Cutpoints[,i]<-as.numeric(as.character(Cutpoints[,i]))}

library(ggplot2)
ggplot(Cutpoints, aes(x = TimeCut, y = Mean)) +  geom_point(size = 3) + 
  geom_errorbar(aes(ymin=(Bottom5), ymax=(Top5), width=.1))+
  xlab("Time Window") + ylab("AUC")+ ggtitle("Predictive Power for Amyloid Positivity \nConsidered by Time Lag")+ 
  ylim(c(0.50, 1.0)) + theme_classic()

cp0 <- cutpointr(df[df$TimeCut == "(-2,0]",], E_ab42, PIBpos, 
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
    hold<-R_COR(Amyloid_Regions, i, TimeCutList[j], "E_ab42")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult[k, 1]<-names(Amyloid_Regions[i])
    RegionalResult[k, 2]<-fisherz(hold$estimate)
    RegionalResult[k, 3]<-hold$p.value
    RegionalResult[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("PIB")

RegionalResult_plot<-RegionalResult[RegionalResult$Region %in% Regions,]


b <- c(-1, 0, 0.5)
MakeBrainPlot("(-8,-6]", b, -1, 0.5, RegionalResult_plot)
MakeBrainPlot("(-6,-4]", b, -1, 0.5, RegionalResult_plot)
MakeBrainPlot("(-4,-2]", b, -1, 0.5, RegionalResult_plot)
MakeBrainPlot("(-2,0]", b, -1, 0.5, RegionalResult_plot)
MakeBrainPlot("(0,2]", b, -1, 0.5, RegionalResult_plot)
MakeBrainPlot("(2,4]", b, -1, 0.5, RegionalResult_plot)




#Now repeating for amyloid positives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                           "Correlation" = NA,
                           "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1,8,9)]
k<-0
for(i in 3:50){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(Amyloid_Regions[Amyloid_Regions$PIBpos == 1,], i, TimeCutList[j], "E_ab42")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(Amyloid_Regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("PIB")
RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]
b <- c(-1, 0, 0.5)

MakeBrainPlot("(-8,-6]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(2,4]", b, -1, 0.5, RegionalResult_plot_PIBpos)

AmyloidAcc<-data.frame("names" = c("medial orbito frontal", "anterior cingulate", "posterior cingulate",
                                   "isthmus cingulate", "precuneus", "superior temporal", "inferior temporal"),
                       "em" = c(-1, -1, -1, -1, -1, -0.5, -0.5))
AmyloidAcc %>% 
  ggseg(mapping=aes(fill=as.numeric(em)), position="stacked", hemisphere = "left") +
  scale_fill_gradientn(limits = c(GradientMin, GradientMax),
                       colours=c("darkred", "white", "darkgreen"),
                       breaks=b, labels=format(b))

AmyloidAcc <- data.frame(cbind(area=c("medial orbito frontal", "caudal anterior cingulate", "rostral anterior cingulate", 
                                      "posterior cingulate",
                                      "isthmus cingulate", "precuneus", "superior temporal", "inferior temporal"),
                              em=c(-1, -1, -1, -1, -1, -1, -0.5, -0.5)),
                     stringsAsFactors=F)

AmyloidAcc %>% 
  ggseg(mapping=aes(fill=as.numeric(em)), position="stacked", hemisphere = "left")+
  scale_fill_gradientn(limits = c(-1, 0.5),
                       colours=c("darkred", "white", "darkgreen"),
                       breaks=b, labels=format(b))


#Now repeating for amyloid negatives only
RegionalResult_PIBpos<-data.frame("Region" = NA,
                                  "Correlation" = NA,
                                  "P" = NA, "TimeCut" = NA)

TimeCutList<-levels(df$TimeCut)
TimeCutList<-TimeCutList[-c(1,8,9)]
k<-0
for(i in 3:50){  
  for(j in 1:length(TimeCutList)){
    k<-k+1
    hold<-R_COR(Amyloid_Regions[Amyloid_Regions$PIBpos == 0,], i, TimeCutList[j], "E_ab42")
    #hold<-R_COR(PET[PET$PIBpos == 1,], i, TimeCutList[j])
    
    RegionalResult_PIBpos[k, 1]<-names(Amyloid_Regions[i])
    RegionalResult_PIBpos[k, 2]<-fisherz(hold$estimate)
    RegionalResult_PIBpos[k, 3]<-hold$p.value
    RegionalResult_PIBpos[k, 4]<-TimeCutList[j]
  }}


Regions<-MakeRegionList("PIB")
RegionalResult_plot_PIBpos<-RegionalResult_PIBpos[RegionalResult_PIBpos$Region %in% Regions,]
b <- c(-1, 0, 0.5)

MakeBrainPlot("(-8,-6]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-6,-4]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-4,-2]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(-2,0]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(0,2]", b, -1, 0.5, RegionalResult_plot_PIBpos)
MakeBrainPlot("(2,4]", b, -1, 0.5, RegionalResult_plot_PIBpos)
