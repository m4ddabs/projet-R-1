# 22210663
# 22210574
# 22010058
#NOTES
#!9 fix commentaire

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
  ev <- eigen(Q_sqrt %*% V %*% Q_sqrt) #eigenvalues accessed by $values vectors by $vectors
  ord <- order(Re(ev$values), decreasing = TRUE)
  ev$vectors <- inv(Q_sqrt) %*% ev$vectors[,ord] #obtaines eigenvectors for VQ (look fischier pdf "answer" added)
  ev$values <- ev$values[ord]
  
  
  #Small eigen values (which should be 0 but approximation errors end up imaginary and then stuff gets messed up)
  A <- c()
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
  percentages_intertie = ev$values[1:ncol(C)]/sum(ev$values[1:ncol(C)])
  
  #1.g For every lin find the coordinates in respect to the columns of A
  A_tilde <- A #ith line = coordinates of ith-variables on c(k) basis
  C_tilde <- C #ith line =coordinates of ith-individue on a(k) basis
  for (i in 1:ncol(A)) {
    A_tilde[,i] <- A_tilde[,i] * ev$values[i]
    C_tilde[,i] <- C_tilde[,i] * ev$values[i]
  }
  
  #1.h Renvoyez la liste résultat comportant les pourcentages d’inertie Λ et les matrices A, C, A˜ et C˜.
  return(list(percentages_intertie, A, C, A_tilde, C_tilde, names(disjoint_table)))
}

plot_individues <- function(afcm, n_axes) {
  
}

#Test
exacm <- read.csv("exacm.csv", header = TRUE, row.names = 1, sep = ";")

AFCM_exacm <- AFCM(exacm)

#Selection des axes
barplot(unlist(AFCM_donnees[1]))
C_tilde <- as.data.frame(AFCM_exacm[5])
A_tilde <- as.data.frame(AFCM_exacm[4])

variable_names <- AFCM_exacm[6]
plot(C_tilde[,1], C_tilde[,2])
text(C_tilde[,1], C_tilde[,2],  labels = variable_names)
     






