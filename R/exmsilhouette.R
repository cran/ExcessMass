exmsilhouette <-
function(xdata, M=1, CutOff=1, steps=50, rug=TRUE , Lambda=NULL, col=FALSE ,rdata=FALSE, label=TRUE){
  if(all(is.na(xdata))){stop('Input vector is empty.')}else{
  if(is.null(Lambda) == TRUE){
    mLambda <- searchMaxLambda(xdata,4*CutOff)
    if(label==TRUE)
    {
      plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*mLambda), 
           xlab=deparse(substitute(xdata)), ylab="Lambda", main='Excess Mass Silhouette', cex.lab=1.1, cex.main=1.65)
      if(rug)
      {rug(xdata, lwd=.75, col="dimgrey")}
    }else{
      plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*mLambda), 
           xlab="", ylab="", frame.plot=FALSE)
      if(rug)
      {rug(xdata, lwd=.75, col="dimgrey")}
    }
    r <- array(list(), c(steps,3))
    if(col==FALSE)
    {
      for(i in 1:steps)
      {
        res <- excessm(xdata, (i* (mLambda/steps)), M=M)
        r[i,] <- c(i* (mLambda/steps), res)
        for(j in 1:(dim((res$intervals))[1]))
        {
          lines(x=c((res$intervals[j,3]),(res$intervals[j,4])),
                y=c((i* mLambda/steps), (i* mLambda/steps)), lwd=1.7)
        }
      }
    }else{
      for(i in 1:steps)
      {
        res <- excessm(xdata, (i* (mLambda/steps)), M=M)
        r[i,] <- c(i* (mLambda/steps), res)
        len <- length(res[[2]])
        color <- colorizing(res[[2]][len])                      
        for(j in 1:(dim((res$intervals))[1]))
        {
          lines(x=c((res$intervals[j,3]),(res$intervals[j,4])),
                y=c((i* mLambda/steps), (i* mLambda/steps)), col=color, lwd=2)
        }
      }
    }
  }else{   
    if(label==TRUE)
    {
      plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*max(Lambda)), 
           xlab=deparse(substitute(xdata)), ylab="Lambda", main='Excess Mass Silhouette', cex.lab=1.1, cex.main=1.65)
      if(rug)
      {rug(xdata, lwd=.75, col="dimgrey")}
    }else{
      plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*max(Lambda)), 
           xlab="", ylab="", frame.plot=FALSE)
      if(rug)
      {rug(xdata, lwd=.75, col="dimgrey")}
    }
    r <- array(list(), c(length(Lambda),2))
    if(col==FALSE){
      for(i in 1:length(Lambda))
      {
        res <- excessm(xdata, Lambda[i], M=M)
        r[i,] <- res
        for(j in 1:(dim((res$intervals))[1]))
        {
          lines(x=c((res$intervals[j,3]),(res$intervals[j,4])),
                y=c(Lambda[i], Lambda[i]), lwd=1.7)
        }
      }
    }else{
      for(i in 1:length(Lambda))
      {
        res <- excessm(xdata, Lambda[i], M=M)
        r[i,] <- res
        len <- length(res[[2]])
        color <- colorizing(res[[2]][len])  
        for(j in 1:(dim((res$intervals))[1]))
        {
          lines(x=c((res$intervals[j,3]),(res$intervals[j,4])),
                y=c(Lambda[i], Lambda[i]), col=color, lwd=2)
        }
      }
    }
  }
  if(rdata==TRUE)
  return(r)
  }
}
