library(stockassessment)

#Standard SAM
cn<-read.ices("data/NEAcod/cn.dat")
cw<-read.ices("data/NEAcod/cw.dat")
dw<-read.ices("data/NEAcod/dw.dat")
lf<-read.ices("data/NEAcod/lf.dat")
lw<-read.ices("data/NEAcod/lw.dat")
mo<-read.ices("data/NEAcod/mo.dat")
nm<-read.ices("data/NEAcod/nm.dat")
pf<-read.ices("data/NEAcod/pf.dat")
pm<-read.ices("data/NEAcod/pm.dat")
sw<-read.ices("data/NEAcod/sw.dat")
surveys<-read.ices("data/NEAcod/survey.dat")

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
conf = defcon(dat)
par<-defpar(dat,conf)
fitStandard<-sam.fit(dat,conf,par)


#Use external variance estimates
varC = as.matrix(read.table("data/NEAcod/variances/varCatch.txt", sep = " "))
attributes(cn)$weight = 1/varC
varS1 = as.matrix(read.table("data/NEAcod/variances/varFLT15Cod.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/varS1
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
conf = defcon(dat)
par<-defpar(dat,conf)
fitWithVar<-sam.fit(dat,conf,par)

#Validate assessment
AIC(fitStandard, fitWithVar)

resStandard = residuals(fitStandard)
resWithVar = residuals(fitWithVar)
retroStandard = retro(fitStandard,year = 7)
retroWithVar = retro(fitWithVar,year = 7)

plot(retroStandard)
mtext("Retro standard settings", line = 2,at = 1980, cex = 2.5)
plot(retroWithVar)
mtext("Retro external variances",  line = 2,at = 1980, cex = 2.5)

plot(resStandard)
mtext("OSA residuals standard settings", line = 2,at = 0, cex = 2.5)
plot(resWithVar)
mtext("OSA residuals external variances", line = 2,at = 0, cex = 2.5)





#Fit SAM with currently used configurations
cn<-read.ices("data/NEAcod/cn.dat")
cw<-read.ices("data/NEAcod/cw.dat")
dw<-read.ices("data/NEAcod/dw.dat")
lf<-read.ices("data/NEAcod/lf.dat")
lw<-read.ices("data/NEAcod/lw.dat")
mo<-read.ices("data/NEAcod/mo.dat")
nm<-read.ices("data/NEAcod/nm.dat")
pf<-read.ices("data/NEAcod/pf.dat")
pm<-read.ices("data/NEAcod/pm.dat")
sw<-read.ices("data/NEAcod/sw.dat")
surveys<-read.ices("data/NEAcod/survey.dat")


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
conf<-loadConf(dat,"scripts/runSAM/NEAcod/model.cfg", patch=TRUE)
par<-defpar(dat,conf)
fitCurrent<-sam.fit(dat,conf,par)
AIC(fitStandard, fitWithVar,fitCurrent)


#Find advised TAC
F_sq = 697412
F_tr = 0.503
set.seed(12345)
forecast(fitCurrent,catchval = c(F_sq,NA,NA,NA),fval = c(NA,F_tr,F_tr,F_tr))





#Fit SAM with new configurations
cn<-read.ices("data/NEAcod/cn.dat")
cw<-read.ices("data/NEAcod/cw.dat")
dw<-read.ices("data/NEAcod/dw.dat")
lf<-read.ices("data/NEAcod/lf.dat")
lw<-read.ices("data/NEAcod/lw.dat")
mo<-read.ices("data/NEAcod/mo.dat")
nm<-read.ices("data/NEAcod/nm.dat")
pf<-read.ices("data/NEAcod/pf.dat")
pm<-read.ices("data/NEAcod/pm.dat")
sw<-read.ices("data/NEAcod/sw.dat")
surveys<-read.ices("data/NEAcod/survey.dat")

varC = as.matrix(read.table("data/NEAcod/variances/varCatch.txt", sep = " "))
attributes(cn)$weight = 1/varC
varS1 = as.matrix(read.table("data/NEAcod/variances/varFLT15Cod.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/varS1

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
conf<-loadConf(dat,"scripts/runSAM/NEAcod/modelNew.cfg", patch=TRUE) #Task: suggest configurations
par<-defpar(dat,conf)
fitNew<-sam.fit(dat,conf,par)
AIC(fitStandard, fitWithVar,fitCurrent,fitNew)

