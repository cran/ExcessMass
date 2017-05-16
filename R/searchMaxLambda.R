searchMaxLambda <-
function(x, limcount = 4, step = 1.05, trylambda = 0.01){  
  if(all(is.na(x))){stop('Input vector is empty.')}else{
  count=length(x)                                    
  searchlim <- limcount / sqrt(count)
  x <- sort(x)
  L <- lambdaweight(x, trylambda)
  localmaxres <- localmax(L)
  if (localmaxres[[3]] > searchlim ){
    repeat {
      trylambda <- trylambda * step
      L <- lambdaweight(x,trylambda)   
      localmaxres <- localmax(L)                       
      if ((localmaxres[[3]] <= searchlim) | (localmaxres[[1]]==localmaxres[[2]])){break}    
    }
  } else {
    repeat {
      trylambda <- trylambda / step   
      L <- lambdaweight(x, trylambda)   
      localmaxres <- localmax(L)  
      if((localmaxres[[3]] > searchlim)|(localmaxres[[1]]==1)&(localmaxres[[2]]==length(x))){break}
    }
  }
  return(trylambda)
  }
}
