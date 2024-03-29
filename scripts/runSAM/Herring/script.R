devtools::install_github("fishfollower/SAM/stockassessment",ref = "initReset")
library(stockassessment)

#Run standard SAM
cn<-read.ices("data/herring/cn.dat")
cw<-read.ices("data/herring/cw.dat")
dw<-read.ices("data/herring/dw.dat")
lf<-read.ices("data/herring/lf.dat")
lw<-read.ices("data/herring/lw.dat")
mo<-read.ices("data/herring/mo.dat")
nm<-read.ices("data/herring/nm.dat")
pf<-read.ices("data/herring/pf.dat")
pm<-read.ices("data/herring/pm.dat")
sw<-read.ices("data/herring/sw.dat")
surveys<-read.ices("data/herring/survey.dat")


dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)

conf = loadConf(dat,"scripts/runSAM/Herring/confStandard.cfg")
par<-defpar(dat,conf)
fitStandard<-sam.fit(dat,conf,par)

#Run SAM with XSAM-options
varC = as.matrix(read.table("data/herring/varC.txt", sep = " "))
attributes(cn)$weight = 1/(varC)
varS1 = as.matrix(read.table("data/herring/varS1.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/(varS1)
varS2 = as.matrix(read.table("data/herring/varS2.txt", sep = " "))
attributes(surveys[[2]])$weight = 1/(varS2)
varS3 = as.matrix(read.table("data/herring/varS3.txt", sep = " "))
attributes(surveys[[3]])$weight = 1/(varS3)

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)


conf = loadConf(dat,"scripts/runSAM/Herring/model.cfg")
par<-defpar(dat,conf)
par$logSdLogN = c(-0.35, -5)
map = list(logSdLogN = as.factor(c(0,NA)))
fitCurrent<-sam.fit(dat,conf,par,map =map)



#Validation of assessment and settings configurations
AIC(fitStandard,fitCurrent)  #AIC

resStandard = residuals(fitStandard) #OSA residuals
resCurrent = residuals(fitCurrent)
plot(resStandard)
mtext("OSA residuals standard settings", line = 1,at = 0, cex = 2.5)
plot(resCurrent)
mtext("OSA residuals current settings", line = 1,at = 0, cex = 2.5)

retroStandard = retro(fitStandard,year =7) #Retro
retroCurrent = retro(fitCurrent,year = 7)
plot(retroStandard)
mtext("Retro standard settings", line = 2,at = 2003, cex = 2.5)
plot(retroCurrent)
mtext("Retro current settings",  line = 2,at = 2003, cex = 2.5)


#NB, leavout currently not working when fixing model parameters
leaveoutStandard= leaveout(fitStandard)
plot(leaveoutStandard)
mtext("Leaveout standard settings",  line = 2,at = 2003, cex = 2.5)
#leaveoutCurrent= leaveout(fitCurrent)
#plot(leaveoutCurrent)
#mtext("Leaveout current settings",  line = 2,at = 2003, cex = 2.5)




#Forecast #Note needs developmentversion of SAM
#devtools::install_github("fishfollower/SAM/stockassessment",ref = "forecastWeightedF")
set.seed(12345)
forecast(fitCurrent,catchval.exact = c(773.750,NA,NA),fBarW = c(NA,0.14,0.14),nosim = 1000,ave.years = c(2016,2017,2018))



#New configurations
cn<-read.ices("data/herring/cn.dat")
cw<-read.ices("data/herring/cw.dat")
dw<-read.ices("data/herring/dw.dat")
lf<-read.ices("data/herring/lf.dat")
lw<-read.ices("data/herring/lw.dat")
mo<-read.ices("data/herring/mo.dat")
nm<-read.ices("data/herring/nm.dat")
pf<-read.ices("data/herring/pf.dat")
pm<-read.ices("data/herring/pm.dat")
sw<-read.ices("data/herring/sw.dat")
surveys<-read.ices("data/herring/survey.dat")

varC = as.matrix(read.table("data/herring/varC.txt", sep = " "))
attributes(cn)$weight = 1/(varC)
varS1 = as.matrix(read.table("data/herring/varS1.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/(varS1)
varS2 = as.matrix(read.table("data/herring/varS2.txt", sep = " "))
attributes(surveys[[2]])$weight = 1/(varS2)
varS3 = as.matrix(read.table("data/herring/varS3.txt", sep = " "))
attributes(surveys[[3]])$weight = 1/(varS3)

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)


conf = loadConf(dat,"scripts/runSAM/Herring/modelNew.cfg") #Task: each make their own suggestion
par<-defpar(dat,conf)
fitNew<-sam.fit(dat,conf,par)

AIC(fitStandard,fitCurrent,fitNew)
corplot(fitNew)#Plot the internally estimated correlation structure
