
library(ggplot2)
library(FLCore)




CatchPlot <- function(data){
  #Plotting function for the catch curves
  temp<- c()
  temp$caton <-rowSums(data,na.rm = TRUE)
  temp$years <- as.integer(names(rowSums(data,na.rm = TRUE)))
  temp<- as.data.frame(temp)
  p<-ggplot(data= temp,aes(x=years,y=caton/1000))+geom_bar(stat="identity")+ylab('Catch')+xlab('Year')+ggtitle('Catch')+
    theme(text = element_text(size=15))
  show(p)
}




arr2coh <- function(object){ 
  # dimensions and co
  
  ages <- as.numeric(attributes(object)$dimnames[[1]])
  years <- as.numeric(attributes(object)$dimname[[2]])
  temp <- c()
  temp$age <- c()
  temp$value <- c()
  temp$cohort <- c()
  temp$survey_year <- c()
  dobj <- dim(object)
  dflq <- dobj
  dflq[2] <- dobj[2] + dobj[1] - 1
  flq <- array(NA, dim = dflq)
  for(i in 1:dflq[1]){ 
    flq[i,(dobj[1] - i + 1):(dflq[2] + 1 - i)] <- as.numeric(object[i,])
    temp$value <- c(temp$value,as.numeric((object[i,])))
    temp$age <- c(temp$age,replicate(length(object[i,]),ages[i]))
    temp$cohort <- c(temp$cohort,years-ages[i])
    temp$survey_year <- c(temp$survey_year,years)
  }
  
  return(temp)
} 




PlotCohort <- function(tmp,cohort_min=1990, ignore0 = TRUE,title=''){
  #plotting functin for the cohort and catch curves  
  tmp <-arr2coh(t(tmp))
  
  tmp <- as.data.frame(tmp)
  
  test <- tmp[tmp$cohort>=cohort_min,]
  
  if(ignore0){
    test$value[test$value==0] = NA
  }
  ggplot(data=test, aes(x=survey_year,log(value),label=age,group=cohort))+
    geom_line()+xlab('Year')+ylab('')+
    xlim(min(test$survey_year),max(test$survey_year))+
    geom_point()+
    ylim(0,1+log(max(test$value,na.rm = TRUE)))+
    facet_wrap(~cohort, scales = "free_y")+ggtitle(title)
}



internalConsistancy <- function(dat){
  
  
  
  ### ======================================================================================================
  ### Custom plotting functions
  ### ======================================================================================================
  pairwiseConsistency <- function(x,show.scales=FALSE,log.scales=TRUE,...) {
    require(grid)
    #Convert to Cohorts
    flc		<-		FLCohort(x@index)
    
    #Convert to log scales if necessary
    if(log.scales) {flc <- log10(flc)}
    
    #Convert to data frame, filter NAs, non-positive values
    cohort.df   	<-  as.data.frame(flc)
    cohort.df       <-  cohort.df[is.finite(cohort.df$data),]
    cohort.df       <-  cohort.df[(cohort.df$data>0),]
    #Subset by age
    age.list        <-  as.list( dimnames(x@index)$age)
    n.ages          <-  length(age.list)
    index.by.age.l	<-	lapply(age.list,function(d)
    {subset(cohort.df,age==d)   })
    #Matched sequential ages together
    paired.ages.l   <-  lapply(as.list(1:(n.ages-1)),function(i) {
      merge(index.by.age.l[[i]],index.by.age.l[[i+1]],by=c("cohort","unit","season","area","iter"),all=FALSE)
    })
    paired.ages     <-  do.call(rbind,paired.ages.l)
    paired.ages$id  <-  paste("Age",paired.ages$age.x,"vs",paired.ages$age.y)
    #Do plot
    p	<-	xyplot((data.y) ~ (data.x) | id, type = "p",
                data = paired.ages,
                scales = list(relation = "free",draw=show.scales),
                as.table=TRUE,
                sub=list(label="Dotted lines are 95% confidence interval for the mean.",cex=0.7),
                panel=function(x,y,...) {
                  #Plot data points
                  #panel.xyplot(x,y,...)
                  cols<-rep("blue",length(x))
                  cols[length(x)]<-"red"
                  cols[length(x)-1]<-"orange"
                  cols[length(x)-2]<-"green"
                  panel.xyplot(x,y,col=cols,...)
                  #Fit linear model and plot, along with confidence limits
                  #but only if there are sufficient data points
                  g   <-  lm(y~x)
                  x.rng   <-  data.frame(x=seq(min(pretty(x)),max(pretty(x)),length.out=10))
                  ci  <-  data.frame(x.rng,predict(g,x.rng,interval="confidence"))
                  panel.xyplot(ci$x,ci$fit,lwd=2,type="l")
                  panel.xyplot(ci$x,ci$upr,lwd=1,lty=2,type="l")
                  panel.xyplot(ci$x,ci$lwr,lwd=1,lty=2,type="l")
                  #Put Rsq on plot in top left corner
                  rsq <-  sprintf("%4.3f",round(summary(g)$r.squared,3))
                  grid::grid.text(label = bquote(r^2 == .(rsq)),x = unit(1, "npc") - unit(0.25, "lines"), y = unit(1,"lines"),just="right",gp=gpar(cex=0.8))
                },
                xlab = if(log.scales) {expression(paste(Log[10]," (Younger Age)"))} else {"Younger Age"},
                ylab = if(log.scales) {expression(paste(Log[10]," (Older Age)"))} else { "Older Age"},
                ...)
  }
  
  
  
  
  # internal consistency  {{{
  plotInternalConsistency <-  function(idx,log.scales=TRUE,
                                       cols=c("white", "yellow", "red"),use.rsq=TRUE,...)
  {
    require(grid)   #Addition of text to panels requires the grid package
    
    # Define colour function
    colFn <- colorRamp(colors=cols)
    
    #Number of ages
    ages <- dimnames(idx@index)[[1]]
    
    #Convert to Cohorts, reshape into appropriate format for splom
    flc <-  if(log.scales) {log10(idx@index)}
    else {idx@index}
    flc  <- as.data.frame(FLCohort(flc))
    flc.wide <-  reshape(flc,direction="wide",timevar=names(flc)[1],idvar=names(flc)[2:6])
    names(flc.wide) <-  gsub("data.","",names(flc.wide))
    
    #Default plot settings
    plot.args <- list(~flc.wide[ages],data=flc.wide, pscales=0,varname.font=2,varname.cex=1.5,
                      xlab = if(log.scales) {expression(paste(Log[10]," (Index Value)"))}
                      else {"Index Value"},
                      ylab = if(log.scales) {expression(paste(Log[10]," (Index Value)"))}
                      else { "Index Value"},
                      sub=list(if(use.rsq) {expression(paste("Lower right panels show the Coefficient of Determination (",italic(r^2),")"))}
                               else { expression(paste("Lower right panels show the Coefficient of Correlation (",italic(r),")"))},cex=0.7),
                      upper.panel=function(x,y,...)
                      {
                        # Filter out NAs
                        both.points  <-  is.finite(x) & is.finite(y)
                        x.filtered <-  x[both.points]
                        y.filtered <-  y[both.points]
                        # Only plot lmline if there is more than one point - colour panel according to rsq.
                        if(length(x.filtered)>2)
                        {
                          r <-  cor(y.filtered,x.filtered)
                          if(use.rsq) {
                            panel.fill(col = rgb(colFn(r^2),max=255))   #Colour panel based on the coefficient of determination (r^2)
                          } else {
                            panel.fill(col = rgb(colFn(0.5*r+0.5),max=255))   #Colour panel based on the correlation coefficient (r)
                          }
                          panel.lmline(x.filtered,y.filtered,lwd=2)
                        }
                        panel.splom(x.filtered,y.filtered,col="black",...)
                      },
                      lower.panel=function(x, y, ...)
                      {
                        #Filter out NAs
                        both.points  <-  is.finite(x) & is.finite(y)
                        x.filtered <-  x[both.points]
                        y.filtered <-  y[both.points]
                        
                        #Calculate r squared - but only if there is enough data to do so
                        if(length(x.filtered)>2)
                        {
                          r <-  cor(y.filtered,x.filtered)
                          if(use.rsq) {
                            panel.colour <- r^2           #Colour & number panel based on the coefficient of determination (r^2)
                            panel.number <- round(r^2,3)
                          } else {
                            panel.colour <- 0.5*r+0.5  #Colour & number panel based on the correlation coefficient (r)
                            panel.number <- round(r,3)
                          }
                          panel.fill(col = rgb(colFn(panel.colour),max=255))   #Colour panel according to rsq
                          grid::grid.text(label =sprintf("%4.3f",panel.number),x = unit(0.5, "npc"),
                                          y = unit(0.5,"npc"),just="center")}})
    
    #Passed settings
    passed.args   <- list(...)
    plot.args[names(passed.args)] <- passed.args
    
    #Do plot
    p <- do.call(splom,plot.args)
    print(p)
    return(p)
  }   # }}}
  
  
  
  #Append source path to list of files
  dat          <-  t(dat)         #Transpose dat so as to line up with the standard dimensions of the quant
  idx			 <-  FLIndex(name='')
  idx@index    <-  FLQuant(as.vector(as.matrix(dat)),dim=dim(dat),dimnames=list(age=gsub("X","",rownames(dat)),year=colnames(dat)))
  
  
  p	<-	pairwiseConsistency(idx,show.scales=TRUE,pch=19,main=idx@name)
  print(p)
#Plot internal consistency for individual index
#This is a candidate function to submit to FLR
  p	<-	plotInternalConsistency(idx,use.rsq=FALSE,cols=c("white","white","white","yellow","red"))
  print(p)
}





getVariance <- function(replicate,survey = surveys[[1]],name='S1'){
  
  
  taylorvar<-function(alfa,beta,n,k,mu){
    #Function that gives the variance of mu=k*mu' where mu' is Norwegian catches or indices from StoX and alpha and beta is estimated for mu' based on the Taylor variance function
    k^(2-beta)*(alfa/n)*(mu)^beta
  }
  
  
  replicate<-replicate[!is.na(replicate$age),]
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
  write.table(sd_I,paste0('Herring/Variance/var',name,'_taylor.dat'),sep=' ')
  
  
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
  write.table(sd_I,paste0('Herring/Variance/var',name,'_emp.dat'),sep=' ')
}



getVarianceCn <- function(replicate= read.csv('Herring/Variance/caa_replicate.txt'),cn = cn){
  
  
  taylorvar<-function(alfa,beta,n,k,mu){
    #Function that gives the variance of mu=k*mu' where mu' is Norwegian catches or indices from StoX and alpha and beta is estimated for mu' based on the Taylor variance function
    k^(2-beta)*(alfa/n)*(mu)^beta
  }
  
  
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
  
  
  #------------------------------------------------------------------------------#
  #Establish scaling factor; i.e. from Norwegian to total catches...
  #------------------------------------------------------------------------------#
  k<-lm(caa~-1+cn,data=tmprep)
  k<-k$coef[1]
  
  #------------------------------------------------------------------------------#
  #Establish alpha and beta in Taylor variance relationship
  #------------------------------------------------------------------------------#
  #Taylor model for variance
  fcaa<-lm(log(v)~log(mean),data=tmprep)
  
  
  
  #Establish structure of data
  sd_C<-cn*NA
  years <- attributes(sd_C)$dimnames[[1]]
  ages <- attributes(sd_C)$dimnames[[2]]
  for(y in years){
    for(a in ages){
      mu<-cn[paste(y),paste(a)]
      v<-taylorvar(alfa=exp(fcaa$coef[1]),beta=fcaa$coef[2],n=1,k=k,mu=mu)
      logv<-log(v/mu^2+1)
      sd_C[which(years==y),which(ages==a)]<-(logv)
    }
  }
  
  sd_C[is.na(sd_C)]<-1
  write.table(sd_C,paste0('Herring/Variance/var_Catch_taylor.dat'),sep=' ')
  
}
