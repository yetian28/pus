#helper function of EPUS
FaF <- function(x,k)
{ #x is N x 2 matrix
  #k is the floor of M/4 in EPUS algorithm
  #the output is the sample index selected from the original dataset in this projection
  N <- nrow(x)
  p <- ncol(x)                           #即p==2成立
  h <- matrix(data = 0,nrow = N,ncol = p)
  for( i in 1:p)
  {
    h[,i] <- (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))#将协变量矩阵按列放缩到[0,1]之间
  }
  #分成十六个格子
  group1 <- which(h[,1] <= 1/4 & h[,2] <= 1/4,  arr.ind = TRUE)
  group2 <- which(h[,1] <= 1/4 & h[,2] <= 1/2 & h[,2] > 1/4, arr.ind = TRUE)
  group3 <- which(h[,1] <= 1/4 & h[,2] <= 3/4 & h[,2] > 1/2, arr.ind = TRUE)
  group4 <- which(h[,1] <= 1/4 & h[,2] <= 1 & h[,2] > 3/4, arr.ind = TRUE)
  
  group5 <- which(h[,1] <= 1/2 & h[,1] > 1/4 & h[,2] <= 1/4,  arr.ind = TRUE)
  group6 <- which(h[,1] <= 1/2 & h[,1] > 1/4 & h[,2] <= 1/2 & h[,2] > 1/4, arr.ind = TRUE)
  group7 <- which(h[,1] <= 1/2 & h[,1] > 1/4 & h[,2] <= 3/4 & h[,2] > 1/2, arr.ind = TRUE)
  group8 <- which(h[,1] <= 1/2 & h[,1] > 1/4 & h[,2] <= 1 & h[,2] > 3/4, arr.ind = TRUE)
  
  group9 <- which(h[,1] <= 3/4 & h[,1] > 1/2 & h[,2] <= 1/4,  arr.ind = TRUE)
  group10 <- which(h[,1] <= 3/4 & h[,1] > 1/2 & h[,2] <= 1/2 & h[,2] > 1/4, arr.ind = TRUE)
  group11 <- which(h[,1] <= 3/4 & h[,1] > 1/2 & h[,2] <= 3/4 & h[,2] > 1/2, arr.ind = TRUE)
  group12 <- which(h[,1] <= 3/4 & h[,1] > 1/2 & h[,2] <= 1 & h[,2] > 3/4, arr.ind = TRUE)
  
  group13 <- which(h[,1] <= 1 & h[,1] > 3/4 & h[,2] <= 1/4,  arr.ind = TRUE)
  group14 <- which(h[,1] <= 1 & h[,1] > 3/4 & h[,2] <= 1/2 & h[,2] > 1/4, arr.ind = TRUE)
  group15 <- which(h[,1] <= 1 & h[,1] > 3/4 & h[,2] <= 3/4 & h[,2] > 1/2, arr.ind = TRUE)
  group16 <- which(h[,1] <= 1 & h[,1] > 3/4 & h[,2] <= 1 & h[,2] > 3/4, arr.ind = TRUE)
  
  if(length(group1)>k)
  {
    #在group1中取离（1/8，1/8）点最近的点
    n1 <- length(group1)
    v1 <- vector(length=n1)
    a1 <- c(1/8,1/8)
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
    #在group2中取离（1/8，3/8）点最近的点
    n2 <- length(group2)
    v2 <- vector(length=n2)
    a2 <- c(1/8,3/8)
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
    #在group3中取离（1/8，5/8）点最近的点
    n3 <- length(group3)
    v3 <- vector(length=n3)
    a3 <- c(1/8,5/8)
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
    #在group4中取离（1/8，7/8）点最近的点
    n4 <- length(group4)
    v4 <- vector(length=n4)
    a4 <- c(1/8,7/8)
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
    #在group5中取离（3/8，1/8）点最近的点
    n5 <- length(group5)
    v5 <- vector(length=n5)
    a5 <- c(3/8,1/8)
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
    #在group6中取离（3/8，3/8）点最近的点
    n6 <- length(group6)
    v6 <- vector(length=n6)
    a6 <- c(3/8,3/8)
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
    #在group7中取离（3/8，5/8）点最近的点
    n7 <- length(group7)
    v7 <- vector(length=n7)
    a7 <- c(3/8,5/8)
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
    #在group8中取离（3/8，7/8）点最近的点
    n8 <- length(group8)
    v8 <- vector(length=n8)
    a8 <- c(3/8,7/8)
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
  if(length(group9)>k)
  {
    #在group9中取离（5/8，1/8）点最近的点
    n9 <- length(group9)
    v9 <- vector(length=n9)
    a9 <- c(5/8,1/8)
    for(j9 in 1:n9)
    {
      v9[j9] <- sum(abs(h[group9[j9],]-a9))
    }
    z9 <- v9[1:k]
    for(i9 in (k+1):n9)
    {
      if(v9[i9]< max(z9))
      {
        z9[which.max(z9)] <- v9[i9]
      }
    }
    index9 <- group9[which(v9 %in% z9)]
  }
  else if(length(group9)>0)
  {
    index9 <- group9
  }
  else
  {
    index9 <- 0
  }
  if(length(group10)>k)
  {
    #在group10中取离（5/8，3/8）点最近的点
    n10 <- length(group10)
    v10 <- vector(length=n10)
    a10 <- c(5/8,3/8)
    for(j10 in 1:n10)
    {
      v10[j10] <- sum(abs(h[group10[j10],]-a10))
    }
    z10 <- v10[1:k]
    for(i10 in (k+1):n10)
    {
      if(v10[i10]< max(z10))
      {
        z10[which.max(z10)] <- v10[i10]
      }
    }
    index10 <- group10[which(v10 %in% z10)]
  }
  else if(length(group10)>0)
  {
    index10 <- group10
  }
  else
  {
    index10 <- 0
  }
  if(length(group11)>k)
  {
    #在group11中取离（5/8，5/8）点最近的点
    n11 <- length(group11)
    v11 <- vector(length=n11)
    a11 <- c(5/8,5/8)
    for(j11 in 1:n11)
    {
      v11[j11] <- sum(abs(h[group11[j11],]-a11))
    }
    z11 <- v11[1:k]
    for(i11 in (k+1):n11)
    {
      if(v11[i11]< max(z11))
      {
        z11[which.max(z11)] <- v11[i11]
      }
    }
    index11 <- group11[which(v11 %in% z11)]
  }
  else if(length(group11)>0)
  {
    index11 <- group11
  }
  else
  {
    index11 <- 0
  }
  if(length(group12)>k)
  {
    #在group12中取离（5/8，7/8）点最近的点
    n12 <- length(group12)
    v12 <- vector(length=n12)
    a12 <- c(5/8,7/8)
    for(j12 in 1:n12)
    {
      v12[j12] <- sum(abs(h[group12[j12],]-a12))
    }
    z12 <- v12[1:k]
    for(i12 in (k+1):n12)
    {
      if(v12[i12]< max(z12))
      {
        z12[which.max(z12)] <- v12[i12]
      }
    }
    index12 <- group12[which(v12 %in% z12)]
  }
  else if(length(group12)>0)
  {
    index12 <- group12
  }
  else
  {
    index12 <- 0
  }
  if(length(group13)>k)
  {
    #在group13中取离（7/8，1/8）点最近的点
    n13 <- length(group13)
    v13 <- vector(length=n13)
    a13 <- c(7/8,1/8)
    for(j13 in 1:n13)
    {
      v13[j13] <- sum(abs(h[group13[j13],]-a13))
    }
    z13 <- v13[1:k]
    for(i13 in (k+1):n13)
    {
      if(v13[i13]< max(z13))
      {
        z13[which.max(z13)] <- v13[i13]
      }
    }
    index13 <- group13[which(v13 %in% z13)]
  }
  else if(length(group13)>0)
  {
    index13 <- group13
  }
  else
  {
    index13 <- 0
  }
  if(length(group14)>k)
  {
    #在group14中取离（7/8，3/8）点最近的点
    n14 <- length(group14)
    v14 <- vector(length=n14)
    a14 <- c(7/8,3/8)
    for(j14 in 1:n14)
    {
      v14[j14] <- sum(abs(h[group14[j14],]-a14))
    }
    z14 <- v14[1:k]
    for(i14 in (k+1):n14)
    {
      if(v14[i14]< max(z14))
      {
        z14[which.max(z14)] <- v14[i14]
      }
    }
    index14 <- group14[which(v14 %in% z14)]
  }
  else if(length(group14)>0)
  {
    index14 <- group14
  }
  else
  {
    index14 <- 0
  }
  if(length(group15)>k)
  {
    #在group15中取离（7/8，5/8）点最近的点
    n15 <- length(group15)
    v15 <- vector(length=n15)
    a15 <- c(7/8,5/8)
    for(j15 in 1:n15)
    {
      v15[j15] <- sum(abs(h[group15[j15],]-a15))
    }
    z15 <- v15[1:k]
    for(i15 in (k+1):n15)
    {
      if(v15[i15]< max(z15))
      {
        z15[which.max(z15)] <- v15[i15]
      }
    }
    index15 <- group15[which(v15 %in% z15)]
  }
  else if(length(group15)>0)
  {
    index15 <- group15
  }
  else
  {
    index15 <- 0
  }
  if(length(group16)>k)
  {
    #在group16中取离（7/8，7/8）点最近的点
    n16 <- length(group16)
    v16 <- vector(length=n16)
    a16 <- c(7/8,7/8)
    for(j16 in 1:n16)
    {
      v16[j16] <- sum(abs(h[group16[j16],]-a16))
    }
    z16 <- v16[1:k]
    for(i16 in (k+1):n16)
    {
      if(v16[i16]< max(z16))
      {
        z16[which.max(z16)] <- v16[i16]
      }
    }
    index16 <- group16[which(v16 %in% z16)]
  }
  else if(length(group16)>0)
  {
    index16 <- group16
  }
  else
  {
    index16 <- 0
  }
  
  index <- c(index1,index2,index3,index4,index5,index6,index7,index8,index9,index10,index11,index12,index13,index14,index15,index16)
  index <- index[index != 0]
  return(index)
  
}