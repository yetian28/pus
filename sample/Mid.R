#helper function of EPUS
source("Fat.R")
source("PUS.R")
Mid <- function(x,n)
{
  #x is the original dataset 
  #n is the size of the subdata
  #output is the subdata selected 
  N <- nrow(x)
  p <- ncol(x)
  A <- combn(p,2)
  m <- ncol(A)
  k <- n/m
  if(k<=4)
  {
    x.sub <- PUS(x,n)
  }
  else 
  {
    v <- n%/%(8*m)
    x.sub <- matrix(data = 0,nrow = 0,ncol = p) #创造一个矩阵
    if(v>0)
    {
      
        for(j in 1:m)
        {
          a <- Fat(x[,A[,j]],v)
          x.sub <- rbind(x.sub,x[a,])
          x <- x[-a,]
        }
      
      N1 <- n-nrow(x.sub)
      x.sub1 <- matrix(data = 0,nrow = 0,ncol = p) #创造一个矩阵
      while(N1>0)
      {
        for(t in 1:m)
        {
          b <- Fat(x[,A[,t]],1)
          x.sub1 <- rbind(x.sub1,x[b,])
          x <- x[-b,]
          
        }
        N1 <- N1-nrow(x.sub1)
      }
      c <- sample(nrow(x.sub1),(n-nrow(x.sub)))
      x.sub <- rbind(x.sub,x.sub1[c,])
    }
    else
    {
      N1 <- n
      while(N1>0)
      {
        for(i in 1:m)
        {
          a <- Fat(x[,A[,i]],1)
          x.sub <- rbind(x.sub,x[a,])
          x <- x[-a,]
        }
        N1 <- N1-nrow(x.sub)
      }
      c <- sample(nrow(x.sub),n)
      x.sub <- x.sub[c,]
    }
  }
  
  return(x.sub)
}
