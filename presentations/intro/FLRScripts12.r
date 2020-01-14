#R-version used:  2.3.1 (2006-06-01)
#First load these libraries from within R before proceeding
library(FLCore)             #version 1.3-3
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(FLAssess)           #version 1.2-1
#library(FLSURBA)            #version 1.2-2 (not used below)
#library(FLXSA)              #version 1.2-2 (not used below)
#library(FLSTF)              #version 1.3-1 (not used below)
#library(FLEDA)              #version 1.3-2 (not used below)
library(boot)
library(MASS)
library(lattice)
library(grid)
old.par <- par(no.readonly = TRUE)
source("cdd Survey scripts.r")
source("Corplots_CMillar.r")
#setwd("C:\\Users\\jd02\\Documents\\current\\ICES WGs\\wgnssk\\2015\\prescreening\\FLRplots\\")
path<-""
path2<-""
path3<-""

cod <- readFLStock(paste(path2,"Cod347.idx",sep=""))
cod.tun <- read.FLIndices(paste(path3,"Cod347_2016.dat",sep=""))     # Q3 age set to 5 for extra plots

units(harvest(cod)) <- "f"
for(i in 1:length(cod.tun)) cod.tun[[i]]@type <- "numbers"
for(i in 1:length(cod.tun)) cod.tun[[i]]@range["plusgroup"] <- NA

win.graph()
par(mfrow=c(2,2))

single_cpue_plot_year(FLIndices(cod.tun[[1]]))
single_cpue_plot_cohort(FLIndices(cod.tun[[1]]))
single_catchcurve_index(FLIndices(cod.tun[[1]]))
single_catchcurvegrad_index(FLIndices(cod.tun[[1]]),list(c(2,4)))
savePlot(file=paste(path,"Fig 14.4a ibtsq1_gam_catchcurve",sep=""),type="png",restoreConsole=T)

c2<-cod.tun[[2]]
single_cpue_plot_year(FLIndices(c2))
single_cpue_plot_cohort(FLIndices(c2))
single_catchcurve_index(FLIndices(c2))
single_catchcurvegrad_index(FLIndices(cod.tun[[2]]),list(c(2,4)))
savePlot(file=paste(path,"Fig 14.4b ibtsq3_gam_catchcurve",sep=""),type="png",restoreConsole=T)

corplotswithin_index(cod.tun,1)
savePlot(file=paste(path,"Fig 14.5a ibtsq1_gam_within",sep=""),type="png",restoreConsole=T)
corplotswithin_index(cod.tun,2)
savePlot(file=paste(path,"Fig 14.5b ibtsq3_gam_within",sep=""),device=dev.cur(),type="png",restoreConsole=T)
#savePlot(file=paste(path,"ibtsq3_gam_within_a",sep=""),device=dev.cur()-1,type="png",restoreConsole=T)
#savePlot(file=paste(path,"ibtsq3_gam_within_b",sep=""),device=dev.cur(),type="png",restoreConsole=T)

corplotsbetween_index(cod.tun,1,2)
savePlot(file=paste(path,"Fig 14.5c ibts_gam_between_a",sep=""),device=dev.cur()-1,type="png",restoreConsole=T)
savePlot(file=paste(path,"Fig 14.5c ibts_gam_between_b",sep=""),device=dev.cur(),type="png",restoreConsole=T)

windows(width=6, height=9)
par(mfrow=c(2,1))
catchcurve(no.discards(cod),c(0,12))
catchcurvegrad(no.discards(cod),c(2,4))
savePlot(file=paste(path,"Fig 14.8 catch_catchcurve",sep=""),type="png",restoreConsole=T)

plot.index.corr(list(FLIndex(catch.n=no.discards(cod)@catch.n[1:12,], name="commercial catch")))
savePlot(file=paste(path3,"catch_withincor",sep=""),type="png",restoreConsole=T)
corplotswithin(no.discards(cod))
savePlot(file=paste(path3,"catch_within_a",sep=""),device=dev.cur()-1,type="png",restoreConsole=T)
savePlot(file=paste(path3,"catch_within_b",sep=""),device=dev.cur(),type="png",restoreConsole=T)

# compute standardized catch proportion at age
cod.spay <- spay(log(trim(no.discards(cod)@catch.n,age=1:10)))
# fine tune
ttl <- list(label="Standardized log-catch proportion at age for cod in IV, IIIa & VIId", cex=1)
yttl <- list(label="age", cex=0.8)
xttl <- list(cex=0.8)
ax <- list(cex=0.7)
# plot
bubbles(age~year|unit, cod.spay,  main=ttl, ylab=yttl, xlab=xttl, scales=ax, bub.scale=5)
#apply(spay(no.discards(cod)@catch.n),c(1,3:5),sum)
#cod.spay <- spay(trim(no.discards(cod)@catch.n),year=1971:2006,age=1:10)

