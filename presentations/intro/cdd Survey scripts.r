#Based on Jose's version of Robs scripts
##------------cpue_plot function--------------------------------
single_cpue_plot_cohort <- function(cpue) {
##	program written by rds, but slightly adapted by jdo
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("FLIndex must be an 'FLIndex' object!")

	##par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		counter <- counter + 1
		if(counter == 5 || counter == 9) {
			##windows(width=7, height=6)
			##par(mfrow=c(2,2))
		}

		yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
		arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
		index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
		index <- log(index)

		t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
		t.      <- cbind(t., t.[,1]-t.[,2])

		suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Cohort",ylab="Log-mean-standardised index"))
		for(j in sort(unique(t.[,2]))) {
			lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
			text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
		}
		title(main=cpue[[i]]@name)
	}
	#par(def.par)
}
#################################################################################################

##------------cpue_plot function--------------------------------
single_cpue_plot_year <- function(cpue) {
##	program written by rds, but slightly adapted by jdo
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("FLIndex must be an 'FLIndex' object!")

	##par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		counter <- counter + 1
		if(counter == 5 || counter == 9) {
			##windows(width=7, height=6)
			##par(mfrow=c(2,2))
		}

		yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
		arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
		index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
		index <- log(index)

		t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
		t.      <- cbind(t., t.[,1])

		suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Year",ylab="Log-mean-standardised index"))
		for(j in sort(unique(t.[,2]))) {
			lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
			text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
		}
		title(main=cpue[[i]]@name)
	}
	#par(def.par)
}


#################################################################################################
#Cod34a.tun <- read.FLIndices(paste(path, "EGFS_2006.tun", sep=""))
#cpue_plot(Cod34a.tun)

##------------catchcurve_index function--------------------------
single_catchcurve_index <- function(cpue) {
##	plots catch curves for survey indices
##	rds 27/01/04
##	modified jdo to work with FLCore 16/08/05
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("Input must be an 'FLIndices' object!")

	##par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		dims  <- dimnames(cpue[[i]]@index)
		nages <- length(as.vector(dims$age))
		if(nages > 2) {
			counter <- counter + 1
			if(counter == 5 || counter == 9) {
				##windows(width=7, height=6)
				##par(mfrow=c(2,2))
			}
			yrange <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
			arange <- cpue[[i]]@range["min"] :cpue[[i]]@range["max"]
			catch  <- as.vector(cpue[[i]]@catch.n)
			effort <- as.vector(cpue[[i]]@effort)

			t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch, rep(effort,each=max(arange)-min(arange)+1))
			cohort <- t.[,1]-t.[,2]
			t.     <- cbind(t.,cohort)
			tdf    <- as.data.frame(t.)

			#suppressWarnings(plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="year of survey",ylab="log-abundance index"))
			plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="Year",ylab="log-abundance index")
			uc     <- sort(unique(tdf[,5]))
			for(j in 1:length(uc)) {
				v. <- tdf[tdf[,5]==uc[j],]
				if(length(v.)>1) {
					lines(v.[,1],log(v.[,3]/v.[,4]),col="red")
					text(v.[,1],log(v.[,3]/v.[,4]),as.character(v.[,2]),col="black",cex=0.6)
				}
			}
			title(main=cpue[[i]]@name)
		}
	}
	#par(def.par)
}

#path <- "D:\\2006\\north sea\\cod\\R\\data\\"
#Cod34a.tun <- read.FLIndices(paste(path, "Cod347_2006R.tun", sep=""))
#catchcurve_index(Cod34a.tun)


##------------catchcurvegrad_index function--------------------------
single_catchcurvegrad_index <- function(cpue,agerange) {
##	returns a vector of catch curve gradients for each cohort
##	rds 26/07/04
##	modified jdo to work with FLCore & FLIndex (instead of FLStock) 16/08/05
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
	if(length(cpue) != length(agerange) || !is.list(agerange))
		stop("Second input must be list of same length as first input")

	#par(mfrow=c(2,2))
	counter <- 0
	reslist <- list(1)
	if(length(cpue)>1) reslist <- c(reslist,c(1:(length(cpue)-1)))
	mina<-agerange

	for(i in 1:length(cpue)) {

		dims  <- dimnames(cpue[[i]]@index)
		nages <- length(as.vector(dims$age))
		mina[[i]] <- agerange[[i]][2]-agerange[[i]][1]+1
		if(nages > 2) {
			counter <- counter + 1
			if(counter == 5 || counter == 9) {
				#windows(width=7, height=6)
				#par(mfrow=c(2,2))
			}
			yrange    <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
			arange    <- cpue[[i]]@range["min"]:cpue[[i]]@range["max"]
			cpu       <- cpue[[i]]@index[as.character(arange),]
			logindrat <- cpu
			logindrat[1:nrow(cpu)-1,1:ncol(cpu)-1] <- -log(cpu[1:nrow(cpu)-1,1:ncol(cpu)-1]/cpu[-1,-1])
			logindrat[nrow(cpu),] <- logindrat[nrow(cpu)-1,]
			logindrat[,ncol(cpu)] <- logindrat[,ncol(cpu)-1]
			cpu       <- as.vector(cpu)
			cpu       <- ifelse(cpu==0,NA,cpu)
			logindrat <- as.vector(logindrat)

			t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu)
			t..    <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, logindrat)
			tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
			tdf.   <- as.data.frame(cbind(t.., t..[,1]-t..[,2]))
			colnames(tdf) <- list("year", "age", "cpue", "cohort")
			colnames(tdf.) <- list("year", "age", "cpue", "cohort")

			tdg    <- tdf[tdf$age>=agerange[[i]][1] & tdf$age<=agerange[[i]][2],]
			tdg.    <- tdf.[tdf.$age>=agerange[[i]][1] & tdf.$age<=agerange[[i]][2],]

			uc     <- sort(unique(tdg[,4]))
			res    <- c(0,0)
			suppressWarnings(plot(tdg.[,4],-tdg.[,3],type="null",xlab="cohort",ylab="negative gradient"))
			for(j in 1:length(uc)) {
				dat  <- tdg[tdg$cohort==uc[j],]
				if(length(dat$year)==mina[[i]]) {
					grad <- lm(log(cpue)~year, data  = dat)
					rss  <- c(uc[j],-as.numeric(grad$coefficients[2]))
					res  <- rbind(res, rss)
				}
			}
			res <- res[-1,]
			colnames(res) <- c("cohort", "negative gradient")
			lines(res)
			points(res)
			title(main=paste(cpue[[i]]@name,"- ages",agerange[[i]][1],"to",agerange[[i]][2],sep=" "))
			reslist[[i]]<-res
		}
	}
	#par(def.par)
	return(reslist)
}

#path <- "D:\\2006\\north sea\\cod\\R\\data\\"
#Cod34a.tun <- read.FLIndices(paste(path, "Cod347_2006R.tun", sep=""))
# 1 c() for esch fleet max and min age range
#catchcurvegrad_index(Cod34a.tun,list(c(2,4),c(2,4),c(2,4)))
######################################################################################################

##------------cpue_plot function--------------------------------
multi_cpue_plot_cohort <- function(cpue) {
##	program written by rds, but slightly adapted by jdo
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("FLIndex must be an 'FLIndex' object!")

	par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		counter <- counter + 1
		if(counter == 5 || counter == 9) {
			windows(width=7, height=6)
			par(mfrow=c(2,2))
		}

		yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
		arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
		index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
		index <- log(index)

		t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
		t.      <- cbind(t., t.[,1]-t.[,2])

		suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Cohort",ylab="Log-mean-standardised index"))
		for(j in sort(unique(t.[,2]))) {
			lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
			text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
		}
		title(main=cpue[[i]]@name)
	}
	par(def.par)
}
#################################################################################################

##------------cpue_plot function--------------------------------
multi_cpue_plot_year <- function(cpue) {
##	program written by rds, but slightly adapted by jdo
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("FLIndex must be an 'FLIndex' object!")

	par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		counter <- counter + 1
		if(counter == 5 || counter == 9) {
			windows(width=7, height=6)
			par(mfrow=c(2,2))
		}

		yrange  <- cpue[[i]]@range[4]:cpue[[i]]@range[5]
		arange  <- cpue[[i]]@range[1] :cpue[[i]]@range[2]
		index   <- as.vector(cpue[[i]]@index/rowMeans(cpue[[i]]@index))
		index <- log(index)

		t.      <- cbind(rep(yrange,each=length(arange)),arange, index)
		t.      <- cbind(t., t.[,1])
    
		suppressWarnings(plot(t.[,4],t.[,3],type="null",xlab="Year",ylab="Log-mean-standardised index"))
		for(j in sort(unique(t.[,2]))) {
			lines(t.[t.[,2]==j,4],t.[t.[,2]==j,3],col="red")
			text (t.[t.[,2]==j,4],t.[t.[,2]==j,3],as.character(t.[t.[,2]==j,2]),col="black", cex=0.6)
		}
		title(main=cpue[[i]]@name)
	}
	par(def.par)
}


#################################################################################################
#Cod34a.tun <- read.FLIndices(paste(path, "EGFS_2006.tun", sep=""))
#cpue_plot(Cod34a.tun)

##------------catchcurve_index function--------------------------
multi_catchcurve_index <- function(cpue) {
##	plots catch curves for survey indices
##	rds 27/01/04
##	modified jdo to work with FLCore 16/08/05
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("Input must be an 'FLIndices' object!")

	par(mfrow=c(2,2))
	counter <- 0

	for(i in 1:length(cpue)) {

		dims  <- dimnames(cpue[[i]]@index)
		nages <- length(as.vector(dims$age))
		if(nages > 2) {
			counter <- counter + 1
			if(counter == 5 || counter == 9) {
				windows(width=7, height=6)
				par(mfrow=c(2,2))
			}
			yrange <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
			arange <- cpue[[i]]@range["min"] :cpue[[i]]@range["max"]
			catch  <- as.vector(cpue[[i]]@catch.n)
			effort <- as.vector(cpue[[i]]@effort)

			t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch, rep(effort,each=max(arange)-min(arange)+1))
			cohort <- t.[,1]-t.[,2]
			t.     <- cbind(t.,cohort)
			tdf    <- as.data.frame(t.)

			#suppressWarnings(plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="year of survey",ylab="log-abundance index"))
			plot(tdf[,1],log(tdf[,3]/tdf[,4]),type="n",xlab="Year",ylab="log-abundance index")
			uc     <- sort(unique(tdf[,5]))
			for(j in 1:length(uc)) {
				v. <- tdf[tdf[,5]==uc[j],]
				if(length(v.)>1) {
					lines(v.[,1],log(v.[,3]/v.[,4]),col="red")
					text(v.[,1],log(v.[,3]/v.[,4]),as.character(v.[,2]),col="black",cex=0.6)
				}
			}
			title(main=cpue[[i]]@name)
		}
	}
	par(def.par)
}

##------------catchcurvegrad_index function--------------------------
multi_catchcurvegrad_index <- function(cpue,agerange) {
##	returns a vector of catch curve gradients for each cohort
##	rds 26/07/04
##	modified jdo to work with FLCore & FLIndex (instead of FLStock) 16/08/05
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
	if(length(cpue) != length(agerange) || !is.list(agerange))
		stop("Second input must be list of same length as first input")

	par(mfrow=c(2,2))
	counter <- 0
	reslist <- list(1)
	if(length(cpue)>1) reslist <- c(reslist,c(1:(length(cpue)-1)))
	mina<-agerange

	for(i in 1:length(cpue)) {

		dims  <- dimnames(cpue[[i]]@index)
		nages <- length(as.vector(dims$age))
		mina[[i]] <- agerange[[i]][2]-agerange[[i]][1]+1
		if(nages > 2) {
			counter <- counter + 1
			if(counter == 5 || counter == 9) {
				windows(width=7, height=6)
				par(mfrow=c(2,2))
			}
			yrange    <- cpue[[i]]@range["minyear"]:cpue[[i]]@range["maxyear"]
			arange    <- cpue[[i]]@range["min"]:cpue[[i]]@range["max"]
			cpu       <- cpue[[i]]@index[as.character(arange),]
			logindrat <- cpu
			logindrat[1:nrow(cpu)-1,1:ncol(cpu)-1] <- -log(cpu[1:nrow(cpu)-1,1:ncol(cpu)-1]/cpu[-1,-1])
			logindrat[nrow(cpu),] <- logindrat[nrow(cpu)-1,]
			logindrat[,ncol(cpu)] <- logindrat[,ncol(cpu)-1]
			cpu       <- as.vector(cpu)
			cpu       <- ifelse(cpu==0,NA,cpu)
			logindrat <- as.vector(logindrat)

			t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu)
			t..    <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, logindrat)
			tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
			tdf.   <- as.data.frame(cbind(t.., t..[,1]-t..[,2]))
			colnames(tdf) <- list("year", "age", "cpue", "cohort")
			colnames(tdf.) <- list("year", "age", "cpue", "cohort")

			tdg    <- tdf[tdf$age>=agerange[[i]][1] & tdf$age<=agerange[[i]][2],]
			tdg.    <- tdf.[tdf.$age>=agerange[[i]][1] & tdf.$age<=agerange[[i]][2],]

			uc     <- sort(unique(tdg[,4]))
			res    <- c(0,0)
			suppressWarnings(plot(tdg.[,4],-tdg.[,3],type="null",xlab="cohort",ylab="negative gradient"))
			for(j in 1:length(uc)) {
				dat  <- tdg[tdg$cohort==uc[j],]
				if(length(dat$year)==mina[[i]]) {
					grad <- lm(log(cpue)~year, data  = dat)
					rss  <- c(uc[j],-as.numeric(grad$coefficients[2]))
					res  <- rbind(res, rss)
				}
			}
			res <- res[-1,]
			colnames(res) <- c("cohort", "negative gradient")
			lines(res)
			points(res)
			title(main=paste(cpue[[i]]@name,"- ages",agerange[[i]][1],"to",agerange[[i]][2],sep=" "))
			reslist[[i]]<-res
		}
	}
	par(def.par)
	return(reslist)
}
##############################################################################################
##------------corplotswithin_index function--------------------------
corplotswithin_index <- function(cpue,fl) {
#Written by jdo [22/08/05], based loosely on functions by rds
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
  win.graph() # windows(width=7, height=6)
	par(mfrow=c(2,2))
	counter <- 0
	dims  <- dimnames(cpue[[fl]]@index)
	nages <- length(as.vector(dims$age))
	yrange    <- cpue[[fl]]@range["minyear"]:cpue[[fl]]@range["maxyear"]
	arange    <- cpue[[fl]]@range["min"]:cpue[[fl]]@range["max"]
	cpu       <- cpue[[fl]]@index[as.character(arange),]
	print(cpu)
  cpu <- log(cpu)
	print(cpu)
	t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu)
	tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
	colnames(tdf) <- list("year", "age", "cpue", "cohort")
	tdf$cohort<-substr(as.character(tdf$cohort),3,4)
	#print(tdf)
	if(nages > 1) {
		for(i in arange[1:(length(arange)-1)]) {
			counter <- counter + 1
			if(counter == 5) {
        win.graph() # windows(width=7, height=6)
				par(mfrow=c(2,2))
				counter <- 0
			}
			d1<-tdf[tdf$age==i,]$cpue[1:(length(tdf[tdf$age==i,]$cpue)-1)]
      d1t<-tdf[tdf$age==i,]$cpue[1:(length(tdf[tdf$age==i,]$cpue)-1)]
			#print("d1")
      #print(d1)
      d2<-tdf[tdf$age==i+1,]$cpue[-1]
      d2t<-tdf[tdf$age==i+1,]$cpue[-1]
			#print("d2")
      #print(d2)
			d3<-tdf[tdf$age==i,]$cohort[1:(length(tdf[tdf$age==i,]$cohort)-1)]
			#print("d3")
      #print(d3)
      d1.<-d1[is.finite(d1)&is.finite(d2)]
			d2.<-d2[is.finite(d1)&is.finite(d2)]
			d3.<-d3[is.finite(d1)&is.finite(d2)]
			d4.<-min(d1.)+0.2*(max(d1.)-min(d1.))
			d<-matrix(c(d1.,d2.),length(d1.),2)
			dtt<-matrix(c(d1t,d2t),length(d1t),2)
			print(d)
			print(dtt)
			newf <- data.frame(x = seq(min(d2.), max(d2.), (max(d2.) - min(d2.))/(length(d2.)-1)))
			predict(lm(d2.~d1.), newf, se.fit = TRUE)
      pred.w.plim <- predict(rlm(d2.~d1.), newf, interval="prediction")
      pred.w.clim <- predict(rlm(d2.~d1.), newf, interval="confidence")
      #matplot(new$x,cbind(pred.w.clim, pred.w.plim[,-1]),
      #       lty=c(1,2,2,3,3), type="l", ,xlab=paste("Log-numbers at age",i),
      #       ylab=paste("Log-numbers at age",i+1))			
			xl <- c(min(d[,1])-1,max(d[,1])+1)
			yl <- c(min(d[,2])-1,max(d[,2])+1) #c(min(pred.w.clim[,2],d2.)-1,max(pred.w.clim[,3],d2.)+1)
			plot(d,type="n", xlim=xl, ylim=yl, 
           xlab=paste("Log-numbers at age",i),
           ylab=paste("Log-numbers at age",i+1))
  		title(main=cpue[[fl]]@name)
			text(d1.,d2.,d3.,col="black",cex=0.6)
			abline(lm(d2.~d1.), xlim=xl, ylim=yl)
			abline(rlm(d2.~d1.,maxit=50),lty=3, xlim=xl, ylim=yl)
 			#print(length(d[,1]))
 			#print()
 			if(is.finite(dtt[length(dtt[,1]),1])& is.finite(dtt[length(dtt[,2]),2]))
      {
      text(d[length(d[,1]),1],d[length(d[,2]),2],label="[  ]",font=2,col="black",cex=1.0)
			}
#			print(d)
#			print(pred.w.clim[order(pred.w.clim$fit,pred.w.clim$lwr,pred.w.clim$upr),])
			#lines(sort(d[,1]),sort(pred.w.clim[,3]), xlim=xl, ylim=yl,lty=2)
			#lines(sort(d[,1]),sort(pred.w.clim[,2]), xlim=xl, ylim=yl,lty=2)

			lines(sort(d[,1]),sort(pred.w.plim[,3]), xlim=xl, ylim=yl,lty=3)
			lines(sort(d[,1]),sort(pred.w.plim[,2]), xlim=xl, ylim=yl,lty=3)
			text(d4.,y=max(yl)-0.5,label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
		}
	}
	par(def.par)
}

#path <- "D:\\2006\\north sea\\cod\\R\\data\\"
#Cod34a.tun <- read.FLIndices(paste(path, "Cod347_2006R.tun", sep=""))
#corplotswithin_index(Cod34a.tun,3)

# par(mfrow=c(2, 2)) #, mai=c(0.5, 0.5, 0.5, 0.5))
# corplotswithin_index(Cod34a.tun,1)
# corplotswithin_index(Cod34a.tun,2)
# par(mfrow=c(3, 2)) #, mai=c(0.5, 0.5, 0.5, 0.5))
# corplotswithin_index(Cod34a.tun,2)
# corplotswithin_index(Cod34a.tun,3)
# par(mfrow=c(2, 2)) #, mai=c(0.5, 0.5, 0.5, 0.5))
# corplotswithin_index(Cod34a.tun,3)

##------------corplotsbetween_index function--------------------------
corplotsbetween_index <- function(cpue,fl1,fl2) {
#Written by jdo [22/08/05], based loosely on functions by rds
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
  win.graph() # windows(width=7, height=6)
	par(mfrow=c(2,2))
	counter <- 0

	minage<-max(cpue[[fl1]]@range["min"],cpue[[fl2]]@range["min"])
	maxage<-min(cpue[[fl1]]@range["max"],cpue[[fl2]]@range["max"])
	minyr<-max(cpue[[fl1]]@range["minyear"],cpue[[fl2]]@range["minyear"])
	maxyr<-min(cpue[[fl1]]@range["maxyear"],cpue[[fl2]]@range["maxyear"])

	arange <- minage:maxage
	yrange <- minyr:maxyr
	cpu1 <- cpue[[fl1]]@index[as.character(arange),as.character(yrange)]
	cpu1 <- log(cpu1)
	cpu2 <- cpue[[fl2]]@index[as.character(arange),as.character(yrange)]
	cpu2 <- log(cpu2)
	t1     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu1)
	t2     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu2)
	tdf1    <- as.data.frame(cbind(t1, t1[,1]-t1[,2]))
	tdf2    <- as.data.frame(cbind(t2, t2[,1]-t2[,2]))
	colnames(tdf1) <- list("year", "age", "cpue", "cohort")
	colnames(tdf2) <- list("year", "age", "cpue", "cohort")
	tdf1$cohort<-substr(as.character(tdf1$cohort),3,4)
	tdf2$cohort<-substr(as.character(tdf2$cohort),3,4)
	print("data")
	
	for(i in arange[1:length(arange)]) {
		counter <- counter + 1
		if(counter == 5) {
		  counter <- 0
      win.graph() # windows(width=7, height=6)
	    par(mfrow=c(2,2))
		}
		d1<-tdf1[tdf1$age==i,]$cpue
		d2<-tdf2[tdf2$age==i,]$cpue
		d3<-tdf1[tdf1$age==i,]$cohort
		d1.<-d1[is.finite(d1)&is.finite(d2)]
		d2.<-d2[is.finite(d1)&is.finite(d2)]
		d3.<-d3[is.finite(d1)&is.finite(d2)]
		d4.<-min(d1.)+0.2*(max(d1.)-min(d1.))
		d<-cbind(d1.,d2.)

		#plot(d,type="n",xlab=paste("Log-numbers:",cpue[[fl1]]@name),ylab=paste("Log-numbers:",cpue[[fl2]]@name))
		#title(main=paste("Age",i))
		#text(d1.,d2.,d3.,col="black",cex=0.6)
		#abline(lm(d2.~d1.))
		#abline(rlm(d2.~d1.,maxit=50),lty=3)
		#text(d4.,y=max(d2.),label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)

		newf <- data.frame(x = seq(min(d2.), max(d2.), (max(d2.) - min(d2.))/(length(d2.)-1)))
		predict(lm(d2.~d1.), newf, se.fit = TRUE)
    pred.w.plim <- predict(lm(d2.~d1.), newf, interval="prediction")
    pred.w.clim <- predict(lm(d2.~d1.), newf, interval="confidence")
    #matplot(new$x,cbind(pred.w.clim, pred.w.plim[,-1]),
    #       lty=c(1,2,2,3,3), type="l", ,xlab=paste("Log-numbers at age",i),
    #       ylab=paste("Log-numbers at age",i+1))			
 	  xl <- c(min(d[,1])-1,max(d[,1])+1)
		yl <- c(min(d[,2])-1,max(d[,2])+1) #c(min(pred.w.clim[,2],d2.)-1,max(pred.w.clim[,3],d2.)+1)
		plot(d,type="n", xlim=xl, ylim=yl, 
           xlab=paste("Log-numbers:",cpue[[fl1]]@name),
           ylab=paste("Log-numbers:",cpue[[fl2]]@name))
    title(main=paste("Age",i))       
		text(d1.,d2.,d3.,col="black",cex=0.6)
		abline(lm(d2.~d1.), xlim=xl, ylim=yl)
		abline(rlm(d2.~d1.,maxit=50),lty=3, xlim=xl, ylim=yl)
    text(d[length(d[,1]),1],d[length(d[,2]),2],label="[  ]",font=2,col="black",cex=1.0)
    
#		print(d)
#		print(pred.w.clim[order(pred.w.clim$fit,pred.w.clim$lwr,pred.w.clim$upr),])
#   lines(sort(d[,1]),sort(pred.w.clim[,3]), xlim=xl, ylim=yl,lty=2)
#   lines(sort(d[,1]),sort(pred.w.clim[,2]), xlim=xl, ylim=yl,lty=2)
		lines(sort(d[,1]),sort(pred.w.plim[,3]), xlim=xl, ylim=yl,lty=3)
		lines(sort(d[,1]),sort(pred.w.plim[,2]), xlim=xl, ylim=yl,lty=3)
		text(d4.,y=max(yl),label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)


	}
	par(def.par)
}

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,1,3)

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,2,3)

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,1,3)
 
#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,1,4)

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,2,3)

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,2,4)

#par(mfrow=c(3,2))
#corplotsbetween_index(Cod34a.tun,3,4)
 


##------------catchcurve function--------------------------
catchcurve <- function(stock, agerange=c(NULL, NULL)) {
#Taken without change from rds
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(stock, "FLStock"))
		stop("First input must be an 'FLStock' object!")
  #windows(width=7, height=6)
	#par(mfrow=c(1,1))
	dims   <- dimnames(stock@catch.n)
	arange <- stock@range["min"] :stock@range["max"]
  nages  <- length(as.vector(dims$age))

	if(nages > 2) {
		yrange <- stock@range["minyear"]:stock@range["maxyear"]
		catch  <- as.vector(stock@catch.n)

		t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch)
		tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))

		suppressWarnings(plot(tdf[,1],log(tdf[,3]),type="null",xlab="year of catch",ylab="log-catch"))
		uc     <- sort(unique(tdf[tdf[,2]==arange,4]))
		for(j in 1:length(uc)) {
      if(!is.null(agerange))
          arange <- agerange[1]:agerange[2]
			v. <- tdf[tdf[,4]==uc[j] ,]
			v. <- v.[is.element(v.[,2], arange),]
			if(length(v.[,1])>1) {

		  lines(v.[,1],log(v.[,3]),col="red")
		  text(v.[,1],log(v.[,3]),as.character(v.[,2]),col="black",cex=0.6)

			}
		}
	}
	#par(def.par)
}

#Cod34a    <- no.discards(read.FLStock(paste(path, "Cod347.idx", sep="")))
#catchcurve(Cod34a,c(0,7))


##------------catchcurvegrad function--------------------------
catchcurvegrad <- function(stock, agerange=c(NULL, NULL)) {
#Adapted from rds by jdo [23/08/05] with slight modification (age range, lm and plotting)
	#def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(stock, "FLStock"))
		stop("First input must be an 'FLStock' object!")

	#par(mfrow=c(1,1))
	dims   <- dimnames(stock@catch.n)
	arange <- stock@range["min"] :stock@range["max"]
	nages  <- length(as.vector(dims$age))

	if(nages > 2) {
		yrange <- stock@range["minyear"]:stock@range["maxyear"]
		catch  <- as.vector(stock@catch.n)

		t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, catch)
		tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
		tdg    <- tdf[tdf$arange>=agerange[1] & tdf$arange<=agerange[2],]
		colnames(tdg) <- list("year", "age", "catch", "cohort")

		res <- c(0,0)
		for(cht in sort(unique(tdg$cohort))) {
			tdc <- tdg[tdg$cohort==cht,]
			colnames(tdc) <- list("year", "age", "catch", "cohort")
			if(length(tdc$year)==agerange[2]-agerange[1]+1) {
				grad <- lm(log(catch)~age, data  = tdc)
				rss  <- c(cht,-as.numeric(grad$coefficients[2]))
				res  <- rbind(res, rss)
			}
		}
		res<-res[-1,]
		colnames(res) <- c("cohort", "negative gradient")
		plot(res)
		lines(res)
		title(paste("Ages",agerange[1],"to",agerange[2]))
	}
	#par(def.par)
	return(res)
}

#par(mfrow=c(2,2))
#catchcurvegrad(Cod34a,c(2,4))
#catchcurvegrad(Cod34a,c(2,5))
#catchcurvegrad(Cod34a,c(2,6))
#catchcurvegrad(Cod34a,c(2,7))

##------------corplotswithin function--------------------------
corplotswithin <- function(stock) {
#Written by jdo [30/08/05], based loosely on functions by rds
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(stock, "FLStock"))
		stop("Input must be an 'FLStock' object!")
  windows()
	par(mfrow=c(3,3))
	counter <- 0
	dims  <- dimnames(stock@catch.n)
	nages <- length(as.vector(dims$age))
	yrange <- stock@range["minyear"]:stock@range["maxyear"]
	arange <- stock@range["min"]:stock@range["max"]
	catch  <- stock@catch.n
	logcatch  <- log(catch)
	t.     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, logcatch)
	tdf    <- as.data.frame(cbind(t., t.[,1]-t.[,2]))
	colnames(tdf) <- list("year", "age", "catch", "cohort")
	tdf$cohort<-substr(as.character(tdf$cohort),3,4)
	if(nages > 1) {
		for(i in arange[1:(length(arange)-1)]) {
			counter <- counter + 1
			if(counter == 10 || counter == 19) {
				windows()
				par(mfrow=c(3,3))
			}
			d1<-tdf[tdf$age==i,]$catch[1:(length(tdf[tdf$age==i,]$catch)-1)]
			d2<-tdf[tdf$age==i+1,]$catch[-1]
			d3<-tdf[tdf$age==i,]$cohort[1:(length(tdf[tdf$age==i,]$cohort)-1)]
			d1.<-d1[is.finite(d1)&is.finite(d2)]
			d2.<-d2[is.finite(d1)&is.finite(d2)]
			d3.<-d3[is.finite(d1)&is.finite(d2)]
			d4.<-min(d1.)+0.2*(max(d1.)-min(d1.))
			d<-matrix(c(d1.,d2.),length(d1.),2)
			plot(d,type="n",xlab=paste("Log-numbers at age",i),ylab=paste("Log-numbers at age",i+1))
			text(d1.,d2.,d3.,col="black",cex=0.6)
			abline(lm(d2.~d1.))
			abline(rlm(d2.~d1.,maxit=50),lty=3)
			text(d4.,y=max(d2.),label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
		}
	}
	par(def.par)
}

#######################################################################
##------------corplotsbetween_index function--------------------------
cor_wts <- function(cpue,fl1,fl2) {
  #print(fl2)
 	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
  windows(width=7, height=6)
	par(mfrow=c(2,2))
	counter <- 0

	minage<-max(cpue[[fl1]]@range["min"],cpue[[fl2]]@range["min"])
	maxage<-min(cpue[[fl1]]@range["max"],cpue[[fl2]]@range["max"])
	minyr<-max(cpue[[fl1]]@range["minyear"],cpue[[fl2]]@range["minyear"])
	maxyr<-min(cpue[[fl1]]@range["maxyear"],cpue[[fl2]]@range["maxyear"])

	arange <- minage:maxage
	print(arange)
	yrange <- minyr:maxyr
	print(yrange)
	
	cpu1 <- cpue[[fl1]]@index[as.character(arange),as.character(yrange)] 
	cpu2 <- cpue[[fl2]]@index[as.character(arange),as.character(yrange)] 
	
	t1     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu1)
	t2     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu2)
	
  colnames(t1) <- list("year", "age", "cpue")
	colnames(t2) <- list("year", "age", "cpue")

  #print(t1)
  #print(t2)
 
 	for(i in arange[1:length(arange)]) {
		counter <- counter + 1
		if(counter == 5) {
		  counter <- 0
			windows(width=7, height=6)
			par(mfrow=c(2,2))
		}
		print
		d1<-t1[t1$age==i,]$cpue
		print(d1)
		d2<-t2[t2$age==i,]$cpue
    print(d2)
 		#d1.<-d1[is.finite(d1)&is.finite(d2)]
		#d2.<-d2[is.finite(d1)&is.finite(d2)]

		d<-cbind(d1,d2)
    #print(d1)
    #print(d2)
    
#		newf <- data.frame(x = seq(min(d2.), max(d2.), (max(d2.) - min(d2.))/(length(d2.)-1)))
#		predict(lm(d2~d1), newf, se.fit = TRUE)
#    pred.w.plim <- predict(lm(d2.~d1.), newf, interval="prediction")
#    pred.w.clim <- predict(lm(d2.~d1.), newf, interval="confidence")
#
# 	  xl <- c(min(d[,1])-1,max(d[,1])+1)
#		yl <- c(min(d[,2])-1,max(d[,2])+1) 
#		
#		plot(d,type="n", xlim=xl, ylim=yl, 
#           xlab=paste("Log-numbers:",cpue[[fl1]]@name),
#           ylab=paste("Log-numbers:",cpue[[fl2]]@name))
#    title(main=paste("Age",i))       
#		text(d1.,d2.,d3.,col="black",cex=0.6)
#		abline(lm(d2.~d1.), xlim=xl, ylim=yl)
#		abline(rlm(d2.~d1.,maxit=50),lty=3, xlim=xl, ylim=yl)
#    text(d[length(d[,1]),1],d[length(d[,2]),2],label="[  ]",font=2,col="black",cex=1.0)
#    
#		lines(sort(d[,1]),sort(pred.w.plim[,3]), xlim=xl, ylim=yl,lty=3)
#		lines(sort(d[,1]),sort(pred.w.plim[,2]), xlim=xl, ylim=yl,lty=3)
#		text(d4.,y=max(yl)+1,label=paste("cor =",substr(as.character(corr(d)),1,5)),col=1)
#
#
	}
	par(def.par)
}

##################################################################################

plotsbetween_index <- function(cpue,fl1,fl2) {
	def.par <- par(no.readonly = TRUE)# save default, for resetting...
	if(!inherits(cpue, "FLIndices"))
		stop("First input must be an 'FLIndices' object!")
  windows(width=7, height=6)
	par(mfrow=c(3,2))
	counter <- 0

	minage<-max(cpue[[fl1]]@range["min"],cpue[[fl2]]@range["min"])
	maxage<-min(cpue[[fl1]]@range["max"],cpue[[fl2]]@range["max"])
	minyr<-max(cpue[[fl1]]@range["minyear"],cpue[[fl2]]@range["minyear"])
	maxyr<-min(cpue[[fl1]]@range["maxyear"],cpue[[fl2]]@range["maxyear"])

	arange <- minage:maxage
	yrange <- minyr:maxyr
	cpu1 <- cpue[[fl1]]@index[as.character(arange),as.character(yrange)] #/rep(cpue[[fl1]]@effort[,as.character(yrange)],each=max(arange)-min(arange)+1)
	cpu1 <- log(cpu1)
	cpu2 <- cpue[[fl2]]@index[as.character(arange),as.character(yrange)] #/rep(cpue[[fl2]]@effort[,as.character(yrange)],each=max(arange)-min(arange)+1)
	cpu2 <- log(cpu2)
	t1     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu1)
	t2     <- cbind(rep(yrange,each=max(arange)-min(arange)+1), arange, cpu2)

  colnames(t1) <- list("year", "age", "cpue")
	colnames(t2) <- list("year", "age", "cpue")
  #print(t1)
  #print(t2)

	t1    <- as.data.frame(t1)
	t2    <- as.data.frame(t2)
    
	#tdf1    <- as.data.frame(cbind(t1, t1[,1]-t1[,2]))
	#tdf2    <- as.data.frame(cbind(t2, t2[,1]-t2[,2]))
	#colnames(tdf1) <- list("year", "age", "cpue", "cohort")
	#colnames(tdf2) <- list("year", "age", "cpue", "cohort")
	#tdf1$cohort<-substr(as.character(tdf1$cohort),3,4)
	#tdf2$cohort<-substr(as.character(tdf2$cohort),3,4)
 
 	for(i in arange[1:length(arange)]) {
		counter <- counter + 1
		if(counter == 7) {
		  counter <- 0
			windows(width=7, height=6)
			par(mfrow=c(3,2))
		}
		d1<-t1[t1$age==i,]$cpue
		print(d1)
		d2<-t2[t2$age==i,]$cpue
		print(d2)
		d3<-t1[t1$age==i,]$year
		print(d3)
		#d1.<-d1[is.finite(d1)&is.finite(d2)]
		#d2.<-d2[is.finite(d1)&is.finite(d2)]
		#d3.<-d3[is.finite(d1)&is.finite(d2)]
		#d<-cbind(d3.,d1.)
    #dd<-cbind(d3.,d2.)
		#plot(d,type="n",xlab="Year",ylab="Log-numbers",ylim=c(min(d1.,d2.),max(d1.,d2.)))

		d1.<-d1[is.finite(d1)]
		d3.<-d3[is.finite(d1)]
    d<-cbind(d3.,d1.)
    print(d)
		d2.<-d2[is.finite(d2)]
		d4.<-d3[is.finite(d2)]
    dd<-cbind(d4.,d2.)
    print(dd)
		plot(d,type="n",xlab="Year",ylab="Log-numbers",ylim=c(min(d1.,d2.),max(d1.,d2.)),xlim=c(min(d3.,d4.),max(d3.,d4.)))
    
    points(d)
		points(dd)

		lines(d)
		lines(dd)
		
	}
	par(def.par)
}
