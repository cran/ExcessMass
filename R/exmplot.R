exmplot <- function(xdata, M=1, CutOff=1, steps=50, Lambda=NULL){
if(all(is.na(xdata))){stop('Input vector is empty.')}else{
  options(warn=-1)
  if(is.null(Lambda) == TRUE){
    mLambda <- searchMaxLambda(xdata,4*CutOff)
    r <- array(list(), c(steps,2))
    for(i in 1:steps)
    {
      r[i,] <- excessm(xdata, (i* (mLambda/steps)), M=M)
    }
    #because of monotonicity smallest excess mass is achied under M=1, largest Lambda
    #because of monotonicity largest excess mass is achied under M=M, smallest Lambda
    plot(1, type="n", xlim=c(r[[steps,2]][[1]],1.025*r[[1,2]][[length(r[1,2])]]), ylim=c(mLambda/steps,mLambda), 
         ylab="Lambda", xlab="Excess Mass", main='Excess Mass Lambda Plot', cex.lab=1.1, cex.main=1.65)
    for(i in 1:M){
      for(j in 1:(steps-1)){
        l1 = length(r[[j,2]])
        l2 = length(r[[j+1,2]])
        if(l1 >= i ){
          if(l2 >= i){
            lines(x=c(r[[j,2]][[i]], r[[j+1,2]][[i]]),y=c(j* mLambda/steps, (j+1)* mLambda/steps), lwd=1.7)
          }else{
            lines(x=c(r[[j,2]][[i]], r[[j+1,2]][[l2]]),y=c(j* mLambda/steps, (j+1)* mLambda/steps), lwd=1.7)
          }          
        }else{
          if(l2 >= i){
            lines(x=c(r[[j,2]][[l1]], r[[j+1,2]][[i]]),y=c(j* mLambda/steps, (j+1)* mLambda/steps), lwd=1.7)
          }else{
            lines(x=c(r[[j,2]][[l1]], r[[j+1,2]][[l2]]),y=c(j* mLambda/steps, (j+1)* mLambda/steps), lwd=1.7)
          }
        }
      }
    }
    max_dist <- vector(length=M-1)
    max_dist_Lambda <- vector(length=M-1)
    for(i in 2:M){
      max_val <- 0
      max_ind <- 0
      for(j in 1:steps){
        if(length(r[[j,2]]) >= i){ #otherwise the excess mass for i and i+1 are equal
          curval <-  r[[j,2]][[i]] - r[[j,2]][[i-1]]
          if(curval > max_val){
            max_val <- curval
            max_ind <- j
          }
        }
      }
      max_dist[i-1] <- max_val
      max_dist_Lambda[i-1] <- max_ind* mLambda/steps #Lambda for which the maximal distance is achieved 
    }
       
  }else{
    Lambda <- sort(Lambda)
    r <- array(list(), c(length(Lambda),2))
    for(i in 1:length(Lambda))
    {
      r[i,] <- excessm(xdata, Lambda[i], M=M)
    }
    plot(1, type="n", xlim=c(r[[length(Lambda),2]][[1]],1.025*r[[1,2]][[length(r[1,2])]]), ylim=c(Lambda[1],Lambda[length(Lambda)]), 
         ylab="Lambda", xlab="Excess Mass", main='Excess Mass Lambda Plot', cex.lab=1.1, cex.main=1.65)
    for(i in 1:M){
      for(j in 1:(length(Lambda)-1)){
        l1 = length(r[[j,2]])
        l2 = length(r[[j+1,2]])
        if(l1 >= i ){
          if(l2 >=i){
            lines(x=c(r[[j,2]][[i]], r[[j+1,2]][[i]]),y=c(Lambda[j], Lambda[j+1]), lwd=1.7)
          }else{
            lines(x=c(r[[j,2]][[i]], r[[j+1,2]][[l2]]),y=c(Lambda[j], Lambda[j+1]), lwd=1.7)
          }          
        }else{
          if(l2 >= i){
            lines(x=c(r[[j,2]][[l1]], r[[j+1,2]][[i]]),y=c(Lambda[j], Lambda[j+1]), lwd=1.7)
          }else{
            lines(x=c(r[[j,2]][[l1]], r[[j+1,2]][[l2]]),y=c(Lambda[j], Lambda[j+1]), lwd=1.7)
          }
        }
      }
    }
    max_dist <- vector(length=M-1)
    max_dist_Lambda <- vector(length=M-1)
    for(i in 2:M){
      max_val <- 0
      max_ind <- 0
      for(j in 1:length(Lambda)){
        if(length(r[[j,2]]) >= i){ 
          curval <- r[[j,2]][[i]] - r[[j,2]][[i-1]]
          if(curval > max_val){
            max_val <- curval
            max_ind <- j
          }
        }
      }
      max_dist[i-1] <- max_val
      max_dist_Lambda[i-1] <- Lambda[max_ind] #Lambda for which the maximal distance is achieved 
    }
    
  }  

  u <- list("max_dist"=max_dist, "max_dist_Lambda"=max_dist_Lambda)
  return(u)
}
}
