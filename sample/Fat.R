#helper function of EPUS
Fat <- function(x,k)
{ #x is N x 2 matrix
  #k is the floor of M/2 in EPUS algorithm
  #the output is the sample index selected from the original dataset in this projection
  N <- nrow(x)
  p <- ncol(x)                           #即p==2成立
  h <- matrix(data = 0,nrow = N,ncol = p)
  for( i in 1:p)
  {
    h[,i] <- (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))#将协变量矩阵按列放缩到[0,1]之间
  }
  t <- sample(c(0,1),1) #格子随机划分
  if(t==0)
  {
    a <- c(1,2)
  }
  else
  {
    a <- c(2,1)
  }
  #分成八(4×2)个格子
  group1 <- which(h[,a[2]] <= 0.5 & h[,a[1]] <= 0.5,  arr.ind = TRUE)
  group2 <- which(h[,a[2]] <= 0.5 & h[,a[1]] > 1/4 & h[,a[1]] <= 1/2, arr.ind = TRUE)
  group3 <- which(h[,a[2]] <= 0.5 & h[,a[1]] > 1/2 & h[,a[1]] <= 3/4, arr.ind = TRUE)
  group4 <- which(h[,a[2]] <= 0.5 & h[,a[1]] > 3/4 & h[,a[1]] <= 1,  arr.ind = TRUE)
  group5 <- which(h[,a[2]] > 1/2 & h[,a[2]] <= 1 & h[,a[1]] <= 0.5, arr.ind = TRUE)
  group6 <- which(h[,a[2]] > 1/2 & h[,a[2]] <= 1 & h[,a[1]] > 1/4 & h[,a[1]] <= 1/2, arr.ind = TRUE)
  group7 <- which(h[,a[2]] > 1/2 & h[,a[2]] <= 1 & h[,a[1]] > 1/2 & h[,a[1]] <= 3/4, arr.ind = TRUE)
  group8 <- which(h[,a[2]] > 1/2 & h[,a[2]] <= 1 & h[,a[1]] > 3/4 & h[,a[1]] <= 1,  arr.ind = TRUE)
  
  if(length(group1)>k)
  {
    #在group1中取离（1/8,1/4）点最近的点
    n1 <- length(group1)
    v1 <- vector(length=n1)
    a1 <- c(1/8,1/4)
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
    index1 <- group1[which(v1 %in% z1)]
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
    #在group2中取离（3/8,1/4）点最近的点
    n2 <- length(group2)
    v2 <- vector(length=n2)
    a2 <- c(3/8,1/4)
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
    #在group3中取离（5/8,1/4）点最近的点
    n3 <- length(group3)
    v3 <- vector(length=n3)
    a3 <- c(5/8,1/4)
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
    #在group4中取离（7/8,1/4）点最近的点
    n4 <- length(group4)
    v4 <- vector(length=n4)
    a4 <- c(7/8,1/4)
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
  if(length(group5)>k)
  {
    #在group5中取离（1/8,3/4）点最近的点
    n5 <- length(group5)
    v5 <- vector(length=n5)
    a5 <- c(1/8,3/4)
    for(j5 in 1:n5)
    {
      v5[j5] <- sum(abs(h[group5[j5],]-a5))
    }
    z5 <- v5[1:k]
    for(i5 in (k+1):n5)
    {
      if(v5[i5]< max(z5))
      {
        z5[which.max(z5)] <- v5[i5]
      }
    }
    index5 <- group5[which(v5 %in% z5)]
  }
  
  else if(length(group5)>0)
  {
    index5 <- group5
  }
  else
  {
    index5 <- 0
  }
  if(length(group6)>k)
  {
    #在group6中取离（3/8,3/4）点最近的点
    n6 <- length(group6)
    v6 <- vector(length=n6)
    a6 <- c(3/8,3/4)
    for(j6 in 1:n6)
    {
      v6[j6] <- sum(abs(h[group6[j6],]-a6))
    }
    z6 <- v6[1:k]
    for(i6 in (k+1):n6)
    {
      if(v6[i6]< max(z6))
      {
        z6[which.max(z6)] <- v6[i6]
      }
    }
    index6 <- group6[which(v6 %in% z6)]
  }
  else if(length(group6)>0)
  {
    index6 <- group6
  }
  else
  {
    index6 <- 0
  }
  if(length(group7)>k)
  {
    #在group1中取离（5/8,3/4）点最近的点
    n7 <- length(group7)
    v7 <- vector(length=n7)
    a7 <- c(5/8,3/4)
    for(j7 in 1:n7)
    {
      v7[j7] <- sum(abs(h[group7[j7],]-a7))
    }
    z7 <- v7[1:k]
    for(i7 in (k+1):n7)
    {
      if(v7[i7]< max(z7))
      {
        z7[which.max(z7)] <- v7[i7]
      }
    }
    index7 <- group7[which(v7 %in% z7)]
  }
  else if(length(group7)>0)
  {
    index7 <- group7
  }
  else
  {
    index7 <- 0
  }
  if(length(group8)>k)
  {
    #在group8中取离（7/8,3/4）点最近的点
    n8 <- length(group8)
    v8 <- vector(length=n8)
    a8 <- c(7/8,3/4)
    for(j8 in 1:n8)
    {
      v8[j8] <- sum(abs(h[group8[j8],]-a8))
    }
    z8 <- v8[1:k]
    for(i8 in (k+1):n8)
    {
      if(v8[i8]< max(z8))
      {
        z8[which.max(z8)] <- v8[i8]
      }
    }
    index8 <- group8[which(v8 %in% z8)]
  }
  else if(length(group8)>0)
  {
    index8 <- group8
  }
  else
  {
    index8 <- 0
  }
  index <- c(index1,index2,index3,index4,index5,index6,index7,index8)
  index <- index[index != 0]
  return(index)
}