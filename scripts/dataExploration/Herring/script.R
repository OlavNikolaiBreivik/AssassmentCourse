library(stockassessment)

#For windows user
if(Sys.info()['sysname'] == 'Windows'){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

source('../../utils.R')

#Run standard SAM
cn<-read.ices("../../../data/herring/cn.dat")
cw<-read.ices("../../../data/herring/cw.dat")
dw<-read.ices("../../../data/herring/dw.dat")
lf<-read.ices("../../../data/herring/lf.dat")
lw<-read.ices("../../../data/herring/lw.dat")
mo<-read.ices("../../../data/herring/mo.dat")
nm<-read.ices("../../../data/herring/nm.dat")
pf<-read.ices("../../../data/herring/pf.dat")
pm<-read.ices("../../../data/herring/pm.dat")
sw<-read.ices("../../../data/herring/sw.dat")
surveys<-read.ices("../../../data/herring/survey.dat")



#Some diagnostic
CatchPlot(cn*cw)


#THis needs some more work
PlotCohort(tmp = cn,title = 'CAA')
PlotCohort(tmp = surveys$SpawninggroundsalongtheNorwegiancoast,title = 'Spawning')
PlotCohort(tmp=surveys$FeedingareasintheNorwegianSeainMay,title = 'May')


#internal consistancy plot
internalConsistancy(dat<-surveys$FeedingareasintheNorwegianSeainMay)
internalConsistancy(dat<-surveys$SpawninggroundsalongtheNorwegiancoast)




#read stox replicate
getVariance(replicate = read.csv('../../../data/Herring/Variance/gyt_replicate.txt'),survey = surveys[[1]],name='S1')
getVariance(replicate = read.csv('../../../data/Herring/Variance/bar_replicate.txt'),survey = surveys[[2]],name='S2')
getVariance(replicate = read.csv('../../../data/Herring/Variance/mai_replicate.txt'),survey = surveys[[3]],name='S3')


#Do the same for catch
getVarianceCn(replicate = read.csv('../../../data/Herring/Variance/caa_replicate.txt'),cn)

