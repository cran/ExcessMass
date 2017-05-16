excessm <-function(x, lambda, M=1, UpToM=FALSE){
if(all(is.na(x))){stop('Input vector is empty.')}else{
  x <- sort(x)
  L <- lambdaweight(x,lambda)
  localmaxres <- localmax(L)
  if(localmaxres[[1]]!=1 & localmaxres[[2]]!=length(x))
  {
    m=matrix(c(1,localmaxres[[1]],1,localmaxres[[1]],localmaxres[[2]],0,localmaxres[[2]],
               length(x),1),3,3, byrow = TRUE)
  }else{if((localmaxres[[1]]==1) & (localmaxres[[2]]==length(x)))
  {
    m=matrix(c(1,length(x),0),1,3, byrow = TRUE)
  }else{
    if(localmaxres[[1]]==localmaxres[[2]])
    {stop('Interval which maximizes excess mass is empty set')
    }else{if(localmaxres[[1]]==1){
      m=matrix(c(1,localmaxres[[2]],0,localmaxres[[2]],length(x),1),2,3, byrow = TRUE)
    }else{ 
      m=matrix(c(1,localmaxres[[1]],1,localmaxres[[1]],length(x),0),2,3, byrow = TRUE)
    }     
    }
  }         
  }
  if(UpToM==TRUE)
  {
    u <- matrix(nrow=(dim(m)[1]-sum(m[,3])), ncol=4,dimnames=list(1:(dim(m)[1]-sum(m[,3])),
                                                                  c("Start Index","End Index","Start Point","End Point")))
    for(i in 1:(dim(m)[1]))
    {
      if(m[i,3]==0)
      {
        u[1,1] <- m[i,1]
        u[1,2] <- m[i,2]
        u[1,3] <- x[m[i,1]]
        u[1,4] <- x[m[i,2]]
      }
    }
    r <- list(u)
  }
  variation <- localmaxres[[3]]
  if(M>1)
  {
    for(i in 2:M)
    {
      n=matrix(nrow=(dim(m)[1]),ncol=3)
      for(j in (1:(dim(m)[1]))) 
      {
        localmaxres <- localmax(L,m[[j,1]],m[[j,2]],m[[j,3]])
        n[j,1] <- localmaxres[[1]]
        n[j,2] <- localmaxres[[2]]
        n[j,3] <- localmaxres[[3]]
      }
      maxemi <- which.max(abs(n[,3])) 
      
      if(n[maxemi,3] == 0){warning('Number of intervals is smaller than M')
                           break}
      
      if((n[maxemi,1] != m[maxemi,1]) & (n[maxemi,2] != m[maxemi,2])) 
      {
        helpm <- matrix(nrow = ((dim(m)[1])+2), ncol=3)
        helpm[maxemi,1] <- m[[maxemi,1]]
        helpm[maxemi,2] <- n[maxemi,1]
        helpm[maxemi+1,1] <- n[maxemi,1]
        helpm[maxemi+1,2] <- n[maxemi,2]
        helpm[maxemi+2,1] <- n[maxemi,2]
        helpm[maxemi+2,2] <- m[[maxemi,2]]
        if((n[maxemi,3])>0)
        {
          helpm[maxemi,3] <- TRUE
          helpm[maxemi+1,3] <- FALSE
          helpm[maxemi+2,3] <- TRUE
        }else{  
          helpm[maxemi,3] <- FALSE
          helpm[maxemi+1,3] <- TRUE
          helpm[maxemi+2,3] <- FALSE
        }
        if(maxemi < (dim(m)[1])){
            helpm[(maxemi+1):(dim(m)[1])+2,]=m[(maxemi+1):(dim(m)[1]),]
        }
      }else{      
        helpm <- matrix(nrow = ((dim(m)[1])+1), ncol=3)
        if((n[maxemi,1]) == 1) 
        {
          helpm[maxemi,1] <- 1 
          helpm[maxemi,2] <- n[maxemi,2]
          helpm[maxemi+1,1] <- n[maxemi,2]
          helpm[maxemi+1,2] <- m[[maxemi,2]]        
          helpm[maxemi,3] <- FALSE
          helpm[maxemi+1,3] <- TRUE
        }else{                        
          helpm[maxemi,1] <- m[[maxemi,1]]
          helpm[maxemi,2] <- n[maxemi,1]
          helpm[maxemi+1,1] <- n[maxemi,1]
          helpm[maxemi+1,2] <- length(x)
          helpm[maxemi,3] <- TRUE
          helpm[maxemi+1,3] <- FALSE
        }
        if(maxemi < (dim(m)[1])){    
            helpm[((maxemi+1):(dim(m)[1])+1),]=m[(maxemi+1):(dim(m)[1]),]
        }
      }
      if(maxemi > 1){ 
        helpm[1:(maxemi-1),]=m[1:(maxemi-1),]
      }  
      m <- helpm
      variation[i] <- variation[i-1] + abs(n[maxemi,3])
      if(UpToM==TRUE)
      {
        j <- 1
        u <- matrix(nrow=(dim(m)[1]-sum(m[,3])), ncol=4,dimnames=list(1:(dim(m)[1]-sum(m[,3])),
                                                                      c("Start Index","End Index","Start Point","End Point")))
        for(i in 1:(dim(m)[1]))
        {
          if(m[i,3]==0)
          {
            u[j,1] <- m[i,1]
            u[j,2] <- m[i,2]
            u[j,3] <- x[m[i,1]]
            u[j,4] <- x[m[i,2]]
            j <- j+1
          }
        }
        r[[length(r)+1]] <- u
      }
    }
  }
  if(UpToM==TRUE)
  {
    r <- list("intervals"=r, "excess_mass"=variation)
  }else{
    j <- 1
    u <- matrix(nrow=(dim(m)[1]-sum(m[,3])), ncol=4,dimnames=list(1:(dim(m)[1]-sum(m[,3])),
                                                                  c("Start Index","End Index","Start Point","End Point")))
    for(i in 1:(dim(m)[1]))
    {
      if(m[i,3]==0)
      {
        u[j,1] <- m[i,1]
        u[j,2] <- m[i,2]
        u[j,3] <- x[m[i,1]]
        u[j,4] <- x[m[i,2]]
        j <- j+1
      }
    }  
    r <- list("intervals"=u,"excess_mass"=variation)
  }
  return(r)
  }
}