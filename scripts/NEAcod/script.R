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
#save(fitStandard,file = "markdown/Herring/fitStandard.Rda")


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
conf = defcon(dat)
par<-defpar(dat,conf)
fitWithVar<-sam.fit(dat,conf,par)



#Validate
AIC(fitStandard, fitWithVar)

resStandard = residuals(fitStandard) 
resfitWithVar = residuals(fitWithVar)  
plot(resStandard)
plot(resWithVar)

retroStandard = retro(fitStandard,year = 7)
retroWithVar = retro(fitWithVar,year = 7)
plot(retroStandard)
plot(retroWithVar)

















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
conf<-loadConf(dat,"scripts/NEAcod/model.cfg", patch=TRUE)
par<-defpar(dat,conf)
fitCurrent<-sam.fit(dat,conf,par)


AIC(fitStandard, fitWithVar,fitCurrent)

resCurrent = residuals(fitCurrent) 
retroCurrent = retro(fitCurrent,year = 7)

plot(retroWithVar)






F_sq = 697000
F_tr = 0.6
forecast(fit,catchval = c(F_sq,NA,NA,NA),fval = c(NA,F_tr,F_tr,F_tr))

fitOnlyCatch = leaveout(fit, fleet = list(c(2,3,4,5)))


