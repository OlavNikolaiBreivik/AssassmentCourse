library(stockassessment)


#For windows user
if(Sys.info()['sysname'] == 'Windows'){
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  }


source('../utils.R')

#Run standard SAM
cn<-read.ices("../../data/herring/cn.dat")
cw<-read.ices("../../data/herring/cw.dat")
dw<-read.ices("../../data/herring/dw.dat")
lf<-read.ices("../../data/herring/lf.dat")
lw<-read.ices("../../data/herring/lw.dat")
mo<-read.ices("../../data/herring/mo.dat")
nm<-read.ices("../../data/herring/nm.dat")
pf<-read.ices("../../data/herring/pf.dat")
pm<-read.ices("../../data/herring/pm.dat")
sw<-read.ices("../../data/herring/sw.dat")
surveys<-read.ices("../../data/herring/survey.dat")


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

conf = loadConf(dat,"../../scripts/Herring/confStandard.cfg")
par<-defpar(dat,conf)
fitStandard<-sam.fit(dat,conf,par)

#Run SAM with XSAM-options
cn<-read.ices("../../data/herring/cn.dat")
cw<-read.ices("../../data/herring/cw.dat")
dw<-read.ices("../../data/herring/dw.dat")
lf<-read.ices("../../data/herring/lf.dat")
lw<-read.ices("../../data/herring/lw.dat")
mo<-read.ices("../../data/herring/mo.dat")
nm<-read.ices("../../data/herring/nm.dat")
pf<-read.ices("../../data/herring/pf.dat")
pm<-read.ices("../../data/herring/pm.dat")
sw<-read.ices("../../data/herring/sw.dat")
surveys<-read.ices("../../data/herring/survey.dat")

varC = as.matrix(read.table("../../data/herring/varC.txt", sep = " "))
attributes(cn)$weight = 1/(varC)
varS1 = as.matrix(read.table("../../data/herring/varS1.txt", sep = " "))
attributes(surveys[[1]])$weight = 1/(varS1)
varS2 = as.matrix(read.table("../../data/herring/varS2.txt", sep = " "))
attributes(surveys[[2]])$weight = 1/(varS2)
varS3 = as.matrix(read.table("../../data/herring/varS3.txt", sep = " "))
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


conf = loadConf(dat,"../../scripts/Herring/model.cfg")
par<-defpar(dat,conf)
par$logSdLogN = c(-0.35, -5)
map = list(logSdLogN = as.factor(c(0,NA)))
fitCurrent<-sam.fit(dat,conf,par,map =map)


#Fit SAM with settings given in modelModified.cfg
conf = loadConf(dat,"../../scripts/Herring/modelModified.cfg")
par<-defpar(dat,conf)
fitNew<-sam.fit(dat,conf,par)



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

leaveoutStandard= leaveout(fitStandard)
plot(leaveoutStandard)
mtext("Leaveout standard settings",  line = 2,at = 2003, cex = 2.5)




#Forecast
set.seed(12345)
forecast(fitCurrent, catchval = c(773.750,0)) #TODO: forecast fails with master-version when sd_N is fixed
set.seed(12345)
forecast(fitCurrent,catchval.exact = c(773.750,NA),fval = c(NA,0.14),nosim = 1000,ave.years = c(2016,2017,2018))






#Some diagnostic
CatchPlot(cn*cw)


#THis needs some more work
PlotCohort(cn,title = 'CAA')
PlotCohort(surveys$SpawninggroundsalongtheNorwegiancoast,title = 'Spawning')
PlotCohort(surveys$FeedingareasintheNorwegianSeainMay,title = 'May')


#internal consistancy plot
internalConsistancy(dat<-surveys$FeedingareasintheNorwegianSeainMay)
internalConsistancy(dat<-surveys$SpawninggroundsalongtheNorwegiancoast)




#read stox replicate
getVariance(replicate = read.csv('../../data/Herring/Variance/gyt_replicate.txt'),survey = surveys[[1]],name='S1')
getVariance(replicate = read.csv('../../data/Herring/Variance/bar_replicate.txt'),survey = surveys[[2]],name='S2')
getVariance(replicate = read.csv('../../data/Herring/Variance/mai_replicate.txt'),survey = surveys[[3]],name='S3')






replicate = read.csv('../../data/Herring/Variance/caa_replicate.txt')



tmp<-cn
tmprep<-replicate
#The ECA point estimates are in the column 'mean'. Add a new column for the 'caa' point estimates used in the modelfitting
#norcaa
tmprep$caa<-NA
maxage<-max(as.numeric(attributes(tmp)$dimnames[[2]]))
minage<-min(as.numeric(attributes(tmp)$dimnames[[2]]))

for(i in 1:nrow(tmprep)){
  if(tmprep$alder[i]<maxage & tmprep$alder[i]>=minage){
    tmprep$caa[i]<-tmp[paste(tmprep$aar[i]),paste(tmprep$alder[i])]
    }
}





taylorvar<-function(alfa,beta,n,k,mu){
  #Function that gives the variance of mu=k*mu' where mu' is Norwegian catches or indices from StoX and alpha and beta is estimated for mu' based on the Taylor variance function
  k^(2-beta)*(alfa/n)*(mu)^beta
}


replicate<-replicate[!is.na(replicate$alder),]
replicate$i <- NA
#------------------------------------------------------------------------------#
#This first function scales the variance to the correct scale
for(i in 1:nrow(replicate)){
  if(replicate$age[i]>=min(as.integer(attributes(survey)$dimnames[[2]])) & 
     replicate$age[i]<=max(as.integer(attributes(survey)$dimnames[[2]]))){
    k<-survey[paste(replicate$year[i]),paste(replicate$age[i])]
    replicate$i[i]<-k
  }
}


#---------------------------------------------------------------------#
#Establish scaling factor; i.e. to match scaling deviations in the two sets
#------------------------------------------------------------------------------#
k<-lm(i~-1+mean,data=replicate)
k<-k$coef[1]


#------------------------------------------------------------------------------#
#Establish alpha and beta in Taylor variance relationship
#------------------------------------------------------------------------------#
#Taylor model for variance
fi<-lm(log(v)~log(mean),data=replicate[replicate$mean>0,])




#Establish structure of data
sd_I<-survey*NA
years <- attributes(sd_I)$dimnames[[1]]
ages <- attributes(sd_I)$dimnames[[2]]
for(y in years){
  for(a in ages){
    mu<-survey[paste(y),paste(a)]
    v<-taylorvar(alfa=exp(fi$coef[1]),beta=fi$coef[2],n=1,k=k,mu=mu)
    logv<-log(v/mu^2+1)
    sd_I[which(years==y),which(ages==a)]<-(logv)
  }
}

sd_I[is.na(sd_I)]<-1
write.table(sd_I,paste0('../../data/herring/variance/var',name,'_taylor.dat'),sep=' ')


#An example of replacing with empirical if available
sd_I<-survey*NA
years <- attributes(sd_I)$dimnames[[1]]
ages <- attributes(sd_I)$dimnames[[2]]
for(y in years){
  for(a in ages){
    tmpl<-replicate[replicate$year==y & replicate$age==a,]
    if(nrow(tmpl)>0){
      sd_I[which(years==y),which(ages==a)]<-(log(tmpl$RSE^2+1))
    }
  }
}

sd_I[is.na(sd_I)]<-1
