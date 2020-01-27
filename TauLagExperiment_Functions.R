#Cleaning up the PET scan dataframe
CleanPET<-function(FILEPATH, SUMMARY, POS, CUTOFF){
  DF<-read.csv(FILEPATH)
  colnames(DF)[5]<-"ID"
  DF$PET_Date<-as.Date(DF$PET_Date, format= "%m/%d/%Y")
  DF<-DF[complete.cases(DF$ID),]
  DF<-DF[complete.cases(DF[,SUMMARY]),]
  DF[,POS]<-ifelse(DF[,SUMMARY] > CUTOFF, 1, 0)
  colnames(DF)[6]<-"PET_Date"
  
  #Keeping most recent amyloid scan only
  DF<-DF[with(DF, order(ID, PET_Date)),]
  DF<-aggregate(DF, list(DF$ID), FUN = tail, 1)
  return(DF)}


#Function to make the linear regression plots and label with Spearman COrrelation info
gplot<-function(DF, CSF, PET, CSF_TITLE, PET_TITLE, XLIM, YLIM, XPOS, YPOS1, YPOS2){
  DF[,CSF]<-log(DF[,CSF])
  DF[, PET]<-log(DF[, PET])
  hold<-cor.test(DF[,PET], DF[,CSF], method = "spearman")
  ggplot(DF, aes_string(x = PET, y = CSF)) + geom_point() + geom_smooth(method = "lm", colour = "grey37", fullrange= TRUE) +
    ylab(CSF_TITLE) + xlab(PET_TITLE) + ggtitle(DF[1, "TimeCut"]) +
    xlim(XLIM) + ylim(YLIM) +
    annotate("text", x = XPOS, y = YPOS1, fontface = 'italic', label = paste("R = ", round(as.numeric(hold$estimate), 3), sep = "")) +
    annotate("text", x = XPOS, y = YPOS2, fontface = 'italic',
             label = paste("p ", ifelse(round(as.numeric(hold$p.value), 3) == 0, "< 0.001", 
                                        paste(" = ", round(as.numeric(hold$p.value), 3), sep = "")))) + theme_classic()}


gplot2<-function(DF, CSF, PET, CSF_TITLE, PET_TITLE, XLIM, YLIM, XPOS, YPOS1, YPOS2){
  DF[,CSF]<-log(DF[,CSF])
  hold<-cor.test(DF[,PET], DF[,CSF], method = "spearman")
  ggplot(DF, aes_string(x = PET, y = CSF)) + geom_point() + geom_smooth(method = "lm", colour = "grey37", fullrange= TRUE) +
    ylab(CSF_TITLE) + xlab(PET_TITLE) + ggtitle(DF[1, "TimeCut"]) +
    xlim(XLIM) + ylim(YLIM) +
    annotate("text", x = XPOS, y = YPOS1, fontface = 'italic', label = paste("R = ", round(as.numeric(hold$estimate), 3), sep = "")) +
    annotate("text", x = XPOS, y = YPOS2, fontface = 'italic',
             label = paste("p ", ifelse(round(as.numeric(hold$p.value), 3) == 0, "< 0.001", 
                                        paste(" = ", round(as.numeric(hold$p.value), 3), sep = "")))) + theme_classic()}


#Function to make dataframe for correlations over time
CreatePlotDF<-function(DF, TIMECUT, PET, CSF){
  Result<-c("TimeCut" = TIMECUT, "Estimate" = cor.test(DF[DF[,"TimeCut"] == TIMECUT, PET],
                                                       DF[DF[,"TimeCut"] == TIMECUT, CSF], method = "spearman")$estimate,
            "CI_lower" = r.con(cor.test(DF[DF[,"TimeCut"] == TIMECUT, PET],
                                        DF[DF[,"TimeCut"] == TIMECUT, CSF], method = "spearman")$estimate,
                               length(DF[DF[,"TimeCut"]  == TIMECUT, "ID"]))[1],
            "CI_upper" = r.con(cor.test(DF[DF[,"TimeCut"] == TIMECUT, PET],
                                        DF[DF[,"TimeCut"] == TIMECUT, CSF], method = "spearman")$estimate,
                               length(DF[DF[,"TimeCut"]  == TIMECUT, "ID"]))[2])
  return(Result)
}

R_COR<-function(DF, index, Window, CSF){cor.test(DF[DF[,"TimeCut"] == Window,index], DF[DF[,"TimeCut"] == Window,CSF], method = "spearman")}


#Tedious list making for ggseg
MakeRegionList<-function(PET){
  Regions<-c(paste(PET, "_fSUVR_rsf_TOT_CTX_SSTSBANK", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_CAUDANTCNG", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_CAUDMIDFRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_CUNEUS", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_ENTORHINAL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_FRNPOLE", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_FUSIFORM", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_INFERPRTL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_INFERTMP", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_INSULA", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_ISTHMUSCNG", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_LATOCC", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_LATORBFRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_LINGUAL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_MEDORBFRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_MIDTMP", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PARACNTRL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PARAHPCMPL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PARSOPCLRS", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PARSORBLS", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PARSTRNGLS", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PERICLCRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_POSTCNTRL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_POSTCNG", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PRECNTRL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_PRECUNEUS", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_ROSANTCNG", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_ROSMIDFRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_SUPERFRN", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_SUPERPRTL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_SUPERTMP", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_SUPRAMRGNL", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_TMPPOLE", sep = ""),
             paste(PET, "_fSUVR_rsf_TOT_CTX_TRANSTMP", sep = ""))
  return(Regions)}


#Getting ready for ggseg
PreptoPlot<-function(TIMECUT, DF){
  hold<-DF[DF[,"TimeCut"] == TIMECUT,]
  results = data.frame(cbind(area=c("banks superior temporal", "caudal anterior cingulate", "caudal middle frontal",
                                    "cuneus", "entorhinal", "frontal pole", "fusiform",
                                    "inferior parietal", "inferior temporal", "insula",
                                    "isthmus cingulate", "lateral occipital", #"lateral orbitofrontal",
                                    "lingual", "medial orbito frontal", "middle temporal",
                                    "para central", "parahippocampal", "pars opercularis",
                                    "pars orbitalis", "pars triangularis",
                                    "pericalcarine", "post central", "posterior cingulate",
                                    "pre central", "precuneus", "rostral anterior cingulate",
                                    "rostral middle frontal", "superior frontal", "superior parietal", "superior temporal",
                                    "supramarginal", "temporal pole", "transverse temporal"),em=hold$Correlation),
                       stringsAsFactors=F)
  results$hemi = "left"
  return(results)}

#Actually making the plot
MakeBrainPlot<-function(TIMECUT, b, GradientMin, GradientMax, DF){
  results<-PreptoPlot(TIMECUT, DF)
  results %>% 
    ggseg(mapping=aes(fill=as.numeric(em)), position="stacked", hemisphere = "left") +
    scale_fill_gradientn(limits = c(GradientMin, GradientMax),
                         colours=c("darkred", "white", "darkgreen"),
                         breaks=b, labels=format(b))}

byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}
