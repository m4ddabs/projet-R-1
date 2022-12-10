# 22210663
# 22210574
# 22010058
#NOTES
#!9 fix commentaire

#Imports 
library(matlib)
library(ggplot2)

#Functions (✿◡‿◡)

#1b
disj_tab <- function(X){
  U <- data.frame()
  for (variable in names(X)){
    for (i in 1:nrow(X)) {
      value <- X[i,variable]
      name <- as.character(value) 
      U[i,name] = 1
    }
  }
  U[is.na(U)] = 0
  return(U)
}


Qprod <- function(u, v,Q){
  return(t(u)%*% Q %*% v)
}

Qnorm <- function(u,Q){
  return(sqrt(Qprod(u,u,Q)))
}


#Exercise 1 
AFCM <- function(X) {
  n0 <- nrow(X)
  p <- ncol(X)
  ind_names <- row.names(X)
  
  #1.b
  disjoint_table <- disj_tab(X)
  U <- data.matrix(disjoint_table)
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
  V <- t(X) %*% D %*% X
  Q_sqrt <- sqrt(Q)
  VQ <- V %*% Q
  ev <- eigen(Q_sqrt %*% V %*% Q_sqrt) 
  ord <- order(Re(ev$values), decreasing = TRUE)
  ev$vectors <- inv(Q_sqrt) %*% ev$vectors[,ord] #obtaines eigenvectors for VQ 
  ev$values <- ev$values[ord]
  
  A <- c()
   for (i in 1:length(ev$values)){
    if (Re(ev$values[i]) > 0){ 
      A <- cbind(A,ev$vectors[,i]/Qnorm(ev$vectors[,i], Q))
    }
  }

  #1.e
  C = c()
  for (i in 1:ncol(A)){  
      C <- cbind(C,1/sqrt(Re(ev$values[i])) * X %*% Q %*% A[,i])
  }
  
  #1.f 
  percentages_intertie = ev$values[1:ncol(C)]/sum(ev$values[1:ncol(C)])
  
  #1.g For every lin find the coordinates in respect to the columns of A
  A_tilde <- A #ith line = coordinates of ith-variables on c(k) basis
  C_tilde <- C #ith line =coordinates of ith-individue on a(k) basis
  for (i in 1:ncol(A)) {
    A_tilde[,i] <- A_tilde[,i] * sqrt(ev$values[i])
    C_tilde[,i] <- C_tilde[,i] * sqrt(ev$values[i])
  }
  
  #1.h Renvoyez la liste résultat comportant les pourcentages d’inertie Λ et les matrices A, C, A˜ et C˜.
  A_tilde <- as.data.frame(A_tilde, row.names = names(disjoint_table))
  C_tilde <- as.data.frame(C_tilde, row.names = ind_names)
  return(list(percentages_intertie, A, C, A_tilde, C_tilde))
}

#Exercise 2
plot_individues <- function(afcm) {
  C_tilde <- afcm[[5]]
  ggplot(C_tilde, aes(x = C_tilde[,1], y = C_tilde[,2], label = row.names(C_tilde))) +
    geom_point() +
    geom_text(position=position_jitter(width=1/2,height=1/2) )
}

#Exercise 3
plot_variables <-  function(afcm) {
  A_tilde <- afcm[[4]]
  ggplot(A_tilde, aes(x = A_tilde[,1], y = A_tilde[,2], label = row.names(A_tilde))) +
    geom_point() +
    geom_text(position=position_jitter(width=1/3,height=1/3) )
}

#Exercise 4
exacm <- read.csv("exacm.csv", header = TRUE, row.names = 1, sep = ";")

AFCM_exacm <- AFCM(exacm)

#selection axis
barplot(unlist(AFCM_exacm[1]))

par(mfrow = c(1,3))
plot_individues(AFCM_exacm)
plot_variables(AFCM_exacm)
     






