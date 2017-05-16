mexmsilhouette <-
function(xdata, M=1:3, CutOff=c(1,2,5), steps=30, Lambda=NULL, col=FALSE, rug=TRUE, rdata=FALSE){
  if(all(is.na(xdata))){stop('Input vector is empty.')}else{
  if(is.null(Lambda))
  {
    if((length(M)==1) & (length(CutOff)==1))
    {
      r <- exmsilhouette(xdata, M=M, CutOff=CutOff, steps=steps, label=TRUE, Lambda=Lambda, col=col, rdata=TRUE)
    }else{
      opar <- par(mfrow=c(length(CutOff),length(M)), mar=c(2.5,2.2,0.5,0.5), oma=c(2,2,3,4.5))
      M <- sort(M)
      CutOff <- sort(CutOff)
      max_mod <- max(M)
      r <- array(list(), c(length(CutOff), length(M),steps,3))
      for(i in 1:(length(CutOff)))
      {
        mLambda <- searchMaxLambda(xdata,4*CutOff[i])
        for(j in 1:steps)
        {
          res <- excessm(xdata, (j* (mLambda/steps)), max_mod, UpToM=TRUE)
          if(length(res$intervals) == max_mod){
            for(k in 1:length(M))
            {
              r[i,k,j,2] <- res$intervals[M[k]]
              r[[i,k,j,1]] <- j* (mLambda/steps)
              r[[i,k,j,3]] <- res$excess_mass
            }
          }else{
            t <- 1
            for(k in 1:max_mod)
            {
              if(k==M[t]){
                if(k <= length(res$intervals)){
                  r[i,t,j,2] <- res$intervals[M[t]]
                  r[[i,t,j,1]] <- j* (mLambda/steps)
                  r[[i,t,j,3]] <- res$excess_mass
                }else{
                  r[i,t,j,2] <- res$intervals[length(res$intervals)]
                  r[[i,t,j,1]] <- j* (mLambda/steps)
                  r[[i,t,j,3]] <- res$excess_mass
                }
                t <- t+1
              }
            }
          }
        }
        for(k in 1:length(M)){
          plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*mLambda), 
               xlab="", ylab="", frame.plot=FALSE)
          if(rug)
          {rug(xdata, lwd=.75, col="dimgrey")}
          if(col==FALSE){
            for(j in 1:steps){
              for(n in (1:dim(r[[i,k,j,2]])[1])){
                lines(x=c((r[[i,k,j,2]])[n,3],(r[[i,k,j,2]])[n,4]), 
                      y=c(j* mLambda/steps,j* mLambda/steps), lwd=1.7)
              }
            }
          }else{
            for(j in 1:steps){
              len <- length(r[[i,k,j,3]])
              color <- colorizing(r[[i,k,j,3]][len])           
              for(n in (1:dim(r[[i,k,j,2]])[1])){
                lines(x=c((r[[i,k,j,2]])[n,3],(r[[i,k,j,2]])[n,4]), 
                      y=c(j* mLambda/steps,j* mLambda/steps), col=color, lwd=2)
              }            
            }
          }
          if(i==1)
          {
            mtext(bquote(M == .(M[k])), 3)
          }
        }
        mtext(bquote(mu == .(CutOff[i])), 4, las=1)
      }        
      title(main= "Excess Mass Silhouettes", outer=TRUE, cex.main=1.65)
      mtext("Lambda", 2, outer=TRUE, cex=1.1)
      mtext(deparse(substitute(x)), 1, outer=TRUE, cex=1.1)
      par(opar)
    }
  }else{
    if((length(M)==1) & (length(CutOff)==1))
    {
      r <- exmsilhouette(xdata, M=M, CutOff=CutOff, steps=steps, label=TRUE, Lambda=Lambda, col=col, rdata=TRUE)
    }else{
      opar <- par(mfrow=c(1,length(M)), mar=c(2.5,2.2,0.5,0.5), oma=c(2,2,3,3.2))
      M <- sort(M)
      max_mod <- max(M)
      r <- array(list(), c(length(M), length(Lambda),2))
      for(j in 1:length(Lambda))
      {
        res <- excessm(xdata, Lambda[j], max_mod, UpToM=TRUE)
        if(length(res$intervals) == max_mod){
          for(k in 1:length(M))
          {
            r[k,j,1] <- res$intervals[M[k]]
            r[[k,j,2]] <- res$excess_mass
          }
        }else{
          t <- 1
          for(k in 1:max_mod)
          {
            if(k==M[t]){
              if(k <= length(res$intervals)){
                r[t,j,1] <- res$intervals[M[t]]
                r[[t,j,2]] <- res$excess_mass
              }else{
                r[t,j,1] <- res$intervals[length(res$intervals)]
                r[[t,j,2]] <- res$excess_mass
              }
              t <- t+1
            }
          }
        }
      }
      for(k in 1:length(M)){
        plot(1, type="n", xlim=range(xdata), ylim=c(0,1.05*max(Lambda)), 
             xlab="", ylab="", frame.plot=FALSE)  
        if(rug)
        {rug(xdata, lwd=.75, col="dimgrey")}
        if(col==FALSE){
          for(j in 1:length(Lambda)){
            for(n in (1:dim(r[[k,j,1]])[1])){
              lines(x=c((r[[k,j,1]])[n,3],(r[[k,j,1]])[n,4]), 
                    y=c(Lambda[j],Lambda[j]), lwd=1.7)
            }
          }
        }else{
          for(j in 1:length(Lambda)){
            len <- length(r[[k,j,2]])
            color <- colorizing(r[[k,j,2]][len]) 
            for(n in (1:dim(r[[k,j,1]])[1])){
                lines(x=c((r[[k,j,1]])[n,3],(r[[k,j,1]])[n,4]), 
                      y=c(Lambda[j],Lambda[j]), col=color, lwd=2)
            }
          }
        }
        mtext(bquote(M == .(M[k])), 3)
      }        
      title(main= "Excess Mass Silhouettes", outer=TRUE, cex.main=1.65)
      mtext("Lambda", 2, outer=TRUE, cex=1.1)
      mtext(deparse(substitute(xdata)), 1, outer=TRUE, cex=1.1)
      par(opar)
    }
  }
  if(rdata==TRUE)
  return(r)
  }
}
