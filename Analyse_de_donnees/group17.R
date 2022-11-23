# 22210663
# 22210574
# 22010058


#Imports 
library(matlib)

#Functions (✿◡‿◡)

#1b
disj_tab <- function(X){
  #doesn't handle NA values, i'm not sure how afcp handles them
  U <- data.frame()
  for (variable in names(X)){
    for (i in 1:nrow(X)) {
      value <- X[i,variable]
      name <- paste(as.character(variable),as.character(value), sep = "_")
      U[i,name] = 1
    }
  }
  #substitues NA values with 0
  U[is.na(U)] = 0
  return(U)
}


Qprod <- function(u, v,Q){ #this is a scalar product
  return(t(u)%*%Q%*%v)
}

Qnorm <- function(u,Q){
  return(sqrt(Qprod(u,u,Q)))
}

#Exercise 1 (Putting all together)
AFCM <- function(X) {
  n0 <- nrow(X)
  p <- ncol(X)
  
  #1.b
  U <- data.matrix(disj_tab(X))
  n <- nrow(U)
  K <- ncol(U)
  D <- diag(n) 
  D <- (1/(n)) * D
  Dsum <- diag(K)
  for(i in 1:K){
    Dsum[i,i]<-sum(U[,i])
  }
  Q <-(1/n*K)*Dsum 
    
  #1.c
  X1 <- n*U %*% inv(Dsum) -1
  
  #1.d
  VQ <- t(X1) %*% D %*% X1%*% Q
  ev <- eigen(VQ) #eigenvalues accessed by $values vectors by $vectors
  A <- ev$vectors
  
  #Small eigen values (which should be 0 but approximation errors end up imaginary and then stuff gets messed up)
  for (i in 1:ncol(A)){
    #(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
    
    A[,i] <- A[,i]/Qnorm(A[,i], Q)
    
    #print(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
  }
    print(A)
    
  #1.e
  C = diag(n)
  
  for (i in 1:ncol(A)){
    C[,i] <- 1/sqrt(ev$values[i]) * X1 %*% A[,i]
    

  }
    
  print(C)
  C1 <- X1 %*% D %*% t(X1) %*% C
  print(C1)
  print("divisione componente per componente")
  print(C1/C)

  #Test: the are eigenvectors of XDX'
  
}

#Test

donnees <- data.frame(question1 = c(1,2,3,2,1,3,1,2,1,1), question2 = c(1,1,1,2,3,4,1,2,3,1))


#test of AFCM function:

X <- donnees

#1.b
U <- data.matrix(disj_tab(X))
n <- nrow(U)
K <- ncol(U)
D <- diag(n) 
D <- (1/(n)) * D
Dsum <- diag(K)
for(i in 1:K){
  Dsum[i,i]<-sum(U[,i])
}
Q <-(1/n*K)*Dsum 

#1.c
X1 <- n*U %*% inv(Dsum) -1

#1.d
VQ <- t(X1) %*% D %*% X1%*% Q
ev <- eigen(VQ) #eigenvalues accessed by $values vectors by $vectors
A <- ev$vectors

#Imaginary eigenvalues popping up because of approximation errors
#Also by construction eigen values are alway real, should convert them
#Also consier them null under a certain value
#We should ask prof hot to handle this cases i think
for (i in 1:ncol(A)){
  #(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
  
  A[,i] <- A[,i]/Qnorm(A[,i], Q)
  
  #print(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
}
print(A)

#1.e
C = diag(n)

for (i in 1:ncol(A)){
  print(1/sqrt(ev$values[i]) * X1 %*% A[,i])
  C[,i] <- 1/sqrt(ev$values[i]) * X1 %*% A[,i]
  
  
}

print(C)
C1 <- X1 %*% t(X1) %*% D %*% C
print(C1)
print("divisione componente per componente")
print(C1/C)

#Test: the are eigenvectors of XDX'
