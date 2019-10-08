library(stockassessment)
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

conf<-loadConf(dat,"scripts/NEAcod/model.cfg", patch=TRUE)
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)

F_sq = 697000
F_tr = 0.6
forecast(fit,catchval = c(F_sq,NA,NA,NA),fval = c(NA,F_tr,F_tr,F_tr))

fitOnlyCatch = leaveout(fit, fleet = list(c(2,3,4,5)))


