# 22210663
# 22210574
# 22010058
#NOTES, maybe to ask prof
#1) VQ is not necesserly symmetric -> it may not be diagonizable (or have complex eigenvectors)... is that right?

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
  return(t(u)%*% Q %*% v)
}

Qnorm <- function(u,Q){
  return(sqrt(Qprod(u,u,Q)))
}

#calculates the coordinates of the v projection on the span(cols of A)
#in respect to the columns of A
A_coords <- function(A,v) {
  
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
  X <- n*U %*% inv(Dsum) -1
  
  #1.d
  VQ <- t(X) %*% D %*% X%*% Q
  ev <- eigen(VQ, symmetric = TRUE) #eigenvalues accessed by $values vectors by $vectors
  A <- c() 
  
  #Small eigen values (which should be 0 but approximation errors end up imaginary and then stuff gets messed up)
  for (i in 1:length(ev$values)){
    #(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
    if (Re(ev$values[i]) > 0){ #We are only interested in the axis with positive eigenvalues (inertia)
      A <- cbind(A,ev$vectors[,i]/Qnorm(ev$vectors[,i], Q))
    }
    
    
    #print(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
  }
    print(A)
    
  #1.e
  C = c()
    
  for (i in 1:ncol(A)){  
      C <- cbind(C,1/sqrt(Re(ev$values[i])) * X %*% Q %*% A[,i])
  }
  
  #1.f 
  percents_intertie = ev$values[1:ncol(C)]/sum(ev$values[1:ncol(C)])
  
  #1.g For every lin find the coordinates in respect to the columns of A
  
  
}

#Test

donnees <- data.frame(question1 = c(1,2,3,2,1,3,1,2,1,1), question2 = c(1,1,1,2,3,4,1,2,3,1))

#AFCM test

X <- disj_tab(donnees)
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

X <- n * U %*% inv(Dsum) -1

#1.d
V <- t(X) %*% D %*% X
Q_sqrt <- sqrt(Q)
ev <- eigen(Q_sqrt %*% V %*% Q_sqrt) #eigenvalues accessed by $values vectors by $vectors
ord <- order(Re(ev$values), decreasing = TRUE)
ev$vectors <- inv(Q) %*% ev$vectors[,ord] #obtaines eigenvectors for VQ (look fischier pdf "answer" added)
ev$values <- ev$values[ord]

#Small eigen values (which should be 0 but approximation errors end up imaginary and then stuff gets messed up)
A <- c() 
for (i in 1:K){

  if (Re(ev$values[i]) > 0){ #We are only interested in the axis with positive eigenvalues (inertia)
    A <- cbind(A, ev$vectors[,i]/Qnorm(ev$vectors[,i], Q))
    print(sprintf("The ith vector is  with eigval%s", as.character(ev$values[i])))
  }
   #print(sprintf("The ith vector is %s with Qnorm %f",paste(A[,i], collapse = " "),Qnorm(A[,i],Q)))
}
print(A)

#1.e
C = c()

for (i in 1:ncol(A)){  
  C <- cbind(C,1/sqrt(Re(ev$values[i])) * X %*% Q %*% A[,i])
}


W <- X %*% Q %*% t(X) %*% D
C1 <- W %*% C
test_eigenvectoriality <- C1/C #the coulums are constante ---> C is made of the eigenvectors of W
#1.f 
percents_intertie = ev$values[1:ncol(C)]/sum(ev$values[1:ncol(C)])


