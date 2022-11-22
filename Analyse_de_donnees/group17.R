# 22210663
# 22210574
# 22010058


#Imports 

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


Qnorm <- function(u, v,Q){
  return(t(u)%*%Q%*%v)
}

DandQ <- function(X){
  n <- nrow(X)
  K <- ncol(X)
  D <-diag(n) 
  D <- (1/(n)) * D
  Dsum <- diag(K)
  for(i in 1:K){
    Dsum[i,i]<-sum(X[,i])
  }
  Q <-(1/n*K)*Dsum 
  return(c(Q,D,Dsum))
}

DandQ(matrix(nrow = 4))





    


#Exercise 1 (Putting all together)
AFCM <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  #1.b
  U <- data.matrix(disj_tab(X))
  DE <- diag(p) #substitute diag(p) with DE
  D <- diag(n)/n
  Q <- diag(p) #sobsitute diag(p) with Q 
    
  #1.c
  X1 <- n*U %*% inv(DE) -1
  
  #1.d
  VQ <- t(X1) %*% D %*% X1%*%
  ev <- eigen(VQ) #eigenvalues accessed by $values vectors by $vectors
  
    
}

#Test

donnees <- data.frame(question1 = c(1,2,3,2,1,3), question2 = c(1,1,1,2,3,4))
U <- disj_tab(donnees)


#Ex 1

#1a)

D <- diag()

#Ex 2
#Ex 3