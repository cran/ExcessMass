localmax <-
function(x,parFrom=1,parTo=length(x),formax=TRUE){
  n <- length(x)                             
  curTo <- parFrom
  curFrom <- parFrom 
  minindex <- parFrom
  if(parFrom==parTo)
  {warning('Start point equals end point')
   localmaxres <- data.frame(parFrom, parFrom, 0)
  }else{if(formax) {                         
    maxVariation <- 0   
    minval <- x[parFrom]   
    maxkrit <- minval                   	   
    for (j in (parFrom + 1) : parTo){
      curval <- x[j]                         
      if (curval < minval){		               
        minval <- curval  		               
        minindex <- j 		                   
        maxkrit <- maxVariation + minval     
      }else{
        if (curval > maxkrit){	             
          maxVariation <- curval - minval    
          curTo <- j  		                   
          curFrom <- minindex   
          maxkrit <- curval	                 
        }
      }       
    }                                        
    if(curFrom==curTo) 
    {
      variation <- 0
    }else{
      variation <- (maxVariation + 1) / n      
    }
    
  } else {                                
    minVariation <- 0
    maxval <- x[parFrom]      
    minkrit <- minVariation + maxval  	     
    for (j in (parFrom + 1) : parTo) {
      curval <- x[j]   
      if (curval < minkrit){		             
        minkrit <- curval
        curTo <- j
        curFrom <- minindex
        minVariation <- curval - maxval
      }else{
        if (curval >= maxval) {		             
          minkrit <- curval + minVariation 
          maxval <- curval
          minindex <- j
        }
      }
    }                                       
    if(curFrom==curTo)                     
    {
      variation <- 0
    }else{
      variation <- (minVariation -1)/n
    }
  }
        IFrom=curFrom
        ITo=curTo      
        localmaxres<-data.frame(IFrom, ITo, variation)
  }
  return(localmaxres)
}
