# plot catch correlations - C.P.Millar, 2007.
## FIGURE 7: plot catch correlations
empty <-
function(x,y,groups,subscripts, panel.number, packet.number) {
  # do nufink
}

diag.panel.cm <-
function(draw,axis.line.col,...)
{
  diag.panel.splom(draw=F, axis.line.col=axis.line.col,...)
}

main.panel.cm <-
function(x,y,groups,subscripts, panel.number, packet.number, .tol=0.05) {
  # do sumfink
  panel.xyplot(x,y,pch=19, col=grey(0.35), cex=0.2)

  # fit a linear model to the data
  lm1 <- lm(y ~ x)
  x1<-0:20/20
  fit <- suppressWarnings(predict.lm(lm1, newdata=data.frame(x=x1), se.fit=T))
  y1 <- fit$fit
  yu <- y1 + 2*fit$se
  yl <- y1 - 2*fit$se

  sig <- identical(anova(lm1)$"Pr(>F)"[1]<.tol,T)

  if (sig) {
    line.f <- list(lwd=3, lty=1, col="#000000")
    line.ci <- list(lwd=2, lty=1, col="red4")
  } else {
    line.f <- list(lwd=1, lty=1, col="#0000FF")
    line.ci <- list(lwd=1, lty=2, col="#0000FF")
  }

  panel.lines(x1,y1, lwd=line.f$lwd, lty=line.f$lty, col=line.f$col)
  panel.lines(x1,yu, lwd=line.ci$lwd, lty=line.ci$lty, col=line.ci$col)
  panel.lines(x1,yl, lwd=line.ci$lwd, lty=line.ci$lty, col=line.ci$col)

  # draw in axes
  grid.draw(linesGrob(x=c(-0.1, 1, 1),
                      y=c( 0  , 0, 1.1),
                      gp = gpar(lwd = 2, col=grey(0.5)),
                      default.units = "npc"))
}

panel.pairs.cm <-
function(z, subscripts, panel.subscripts) {
  panel.pairs(z,
    panel.subscripts = panel.subscripts,
    subscripts       = subscripts,
    lower.panel      = empty,
    diag.panel       = diag.panel.cm,
    upper.panel      = main.panel.cm,
    axis.line.col    = "#FFFFFF",
    )
}

centre.log <-
function(mat) {
  mat[mat<=0] <- NA
  mat <- log(mat)
  apply(mat,2,function(x) (x-min(x, na.rm=T))/diff(range(x,na.rm=T)))
}

plot.index.corr <-
function(object, wndows=T) {
  par(bty="n")
  trellis.par.set(box.rectangle=list(col="white"))
  for (i in seq(length(object))) {
    #select one tuning fleet
    tune.mat <- t(object[[i]]@catch.n@.Data[,,1,1,1])
    # make cohort matrix
	  n <- dim(tune.mat)[2]
	  cohort.mat <- matrix(NA, ncol=n, nrow=dim(tune.mat)[1]+n-1)
	  colnames(cohort.mat) <- colnames(tune.mat)
	  for (j in 1:n) {
		  cohort.mat[,j] <- c(rep(NA,n-j),tune.mat[,j],rep(NA,j-1))
    }
    main <- object[[i]]@name
    if (wndows) windows()
    print(splom(~centre.log(cohort.mat), superpanel=panel.pairs.cm, xlab=main, col="white"))
  }
}





