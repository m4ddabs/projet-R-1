# 22210663
# 22210574
# 22010058


#Imports 

#Functions (✿◡‿◡)


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
    
#Test

donnees <- data.frame(question1 = c(1,2,3,2,1,3), question2 = c(1,1,1,2,3,4))
disj_tab(donnees)


#Ex 1
#Ex 2
#Ex 3