#PUS Algorithm
#helper function of EPUS
source("TaT.R")
PUS <- function(x,n)
{ #x is the original dataset 
  #n is the size of the subdata
  #output is the subdata selected 
  N <- nrow(x)
  p <- ncol(x)
  #将矩阵按列对分类排序
  A <- combn(p,2)
  m <- ncol(A)
  #计算可以重复取的整数倍
  k <- n%/%(4*m)
  x.sub <- matrix(data = 0,nrow = 0,ncol = p) #创造一个矩阵
  if(k>0)
  {
      for(j in 1:m)
      {
        a <- TaT(x[,A[,j]],k)
        x.sub <- rbind(x.sub,x[a,])
        x <- x[-a,]
      }    #获取第一轮数据
    N1 <- n-nrow(x.sub) #未获取数据
    x.sub1 <- matrix(data = 0,nrow = 0,ncol = p) #创造一个矩阵
    while(N1>0)
    {
      for(t in 1:m)
      {
        b <- TaT(x[,A[,t]],1)
        x.sub1 <- rbind(x.sub1,x[b,])
        x <- x[-b,]
      }
      N1 <- N1-nrow(x.sub1)
    }
    c <- sample(nrow(x.sub1),(n-nrow(x.sub))) #随机抽取未获得数据
    
    x.sub <- rbind(x.sub,x.sub1[c,])
  }
  else
  {
    N1 <- n
    while(N1>0)
    {
      for(i in 1:m)
      {
        a <- TaT(x[,A[,i]],1)
        x.sub <- rbind(x.sub,x[a,])
        x <- x[-a,]
      }
      N1 <- N1-nrow(x.sub)
    }
    c <- sample(nrow(x.sub),n)
    x.sub <- x.sub[c,]
  }
  
  return(x.sub)
}
