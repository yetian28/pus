#helper function of PUS
TaT <- function(x,k)
{
  #x is N x 2 matrix
  #k is the floor of M in PUS algorithm
  #the output is the sample index selected from the original dataset in this projection
  N <- nrow(x)
  p <- ncol(x)                           #即p==2成立
  h <- matrix(data = 0,nrow = N,ncol = p)
  for( i in 1:p)
  {
    h[,i] <- (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))#将协变量矩阵按列放缩到[0,1]之间
  }
  #分成四个格子
  group1 <- which(h[,1] <= 0.5 & h[,2] <= 0.5,  arr.ind = TRUE)
  group2 <- which(h[,1] <= 0.5 & h[,2] > 0.5, arr.ind = TRUE)
  group3 <- which(h[,1] > 0.5 & h[,2] > 0.5, arr.ind = TRUE)
  group4 <- which(h[,1] > 0.5 & h[,2] <= 0.5, arr.ind = TRUE)
  
  if(length(group1)>k)
  {
    #在group1中取离（0.25，0.25）点最近的点
    n1 <- length(group1)
    v1 <- vector(length=n1)
    a1 <- c(0.25,0.25)
    for(j1 in 1:n1)
    {
      v1[j1] <- sum(abs(h[group1[j1],]-a1))
    }
    z1 <- v1[1:k]
    for(i1 in (k+1):n1)
    {
      if(v1[i1]< max(z1))
          {
            z1[which.max(z1)] <- v1[i1] 
      }
    }
    index1 <- group1[which(v1 %in% z1)] #定位位置
  }
  else if(length(group1)>0)
  {
    index1 <- group1
  }
  else
  {
    index1 <- 0
  }
  if(length(group2)>k)
  {
    #在group2中取离（0.25，0.75）点最近的点
    n2 <- length(group2)
    v2 <- vector(length=n2)
    a2 <- c(0.25,0.75)
    for(j2 in 1:n2)
    {
      v2[j2] <- sum(abs(h[group2[j2],]-a2))
    }
    z2 <- v2[1:k]
    for(i2 in (k+1):n2)
    {
      if(v2[i2]< max(z2))
      {
        z2[which.max(z2)] <- v2[i2]
      }
    }
    index2 <- group2[which(v2 %in% z2)]
  }
  else if(length(group2)>0)
  {
    index2 <- group2
  }
  else
  {
    index2 <- 0
  }
  if(length(group3)>k)
  {
    #在group3中取离（0.75，0.75）点最近的点
    n3 <- length(group3)
    v3 <- vector(length=n3)
    a3 <- c(0.75,0.75)
    for(j3 in 1:n3)
    {
      v3[j3] <- sum(abs(h[group3[j3],]-a3))
    }
    z3 <- v3[1:k]
    for(i3 in (k+1):n3)
    {
      if(v3[i3]< max(z3))
      {
        z3[which.max(z3)] <- v3[i3]
      }
    }
    index3 <- group3[which(v3 %in% z3)]
  }
  else if(length(group3)>0)
  {
    index3 <- group3
  }
  else
  {
    index3 <- 0
  }
  if(length(group4)>k)
  {
    #在group4中取离（0.75，0.25）点最近的点
    n4 <- length(group4)
    v4 <- vector(length=n4)
    a4 <- c(0.75,0.25)
    for(j4 in 1:n4)
    {
      v4[j4] <- sum(abs(h[group4[j4],]-a4))
    }
    z4 <- v4[1:k]
    for(i4 in (k+1):n4)
    {
      if(v4[i4]< max(z4))
      {
        z4[which.max(z4)] <- v4[i4]
      }
    }
    index4 <- group4[which(v4 %in% z4)]
  }
  else if(length(group4)>0)
  {
    index4 <- group4
  }
  else
  {
    index4 <- 0
  }
  index <- c(index1,index2,index3,index4)  
  index <- index[index != 0] 
  return(index)
}