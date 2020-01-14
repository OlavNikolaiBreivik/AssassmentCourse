setwd(dirname(rstudioapi::getSourceEditorContext()$path))


library(stockassessment)
library(ggplot2)
source('script/utils.r')


cn<-read.ices("Herring/cn.dat")
cw<-read.ices("Herring/cw.dat")
dw<-read.ices("Herring/dw.dat")
lf<-read.ices("Herring/lf.dat")
lw<-read.ices("Herring/lw.dat")
mo<-read.ices("Herring/mo.dat")
nm<-read.ices("Herring/nm.dat")
pf<-read.ices("Herring/pf.dat")
pm<-read.ices("Herring/pm.dat")
sw<-read.ices("Herring/sw.dat")
surveys<-read.ices("Herring/survey.dat")

#read stox replicate
getVariance(replicate = read.csv('Herring/Variance/gyt_replicate.txt'),survey = surveys[[1]],name='S1')
getVariance(replicate = read.csv('Herring/Variance/bar_replicate.txt'),survey = surveys[[2]],name='S2')
getVariance(replicate = read.csv('Herring/Variance/mai_replicate.txt'),survey = surveys[[3]],name='S3')
getVarianceCn(replicate = read.csv('Herring/Variance/caa_replicate.txt'),cn)