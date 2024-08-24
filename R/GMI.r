#'Internal function for rdacca.hp()
GMI<-function(X,tol=0.000001){ # GMI stands for Generalized.Matrix.Inversion
  # X is a square matrix
  # this matrix inversion procedure invert matrices containing highly colinear predictors
  # that cannot be inverted by the standard R inversion function "solve"
  # if solve can't invert the matrix, predictors are replaced by their principal components
  # associated to eigenvalues greater than 0; otherwise the matrix is inverted by solve
  SVD<-svd(X)
  d<-SVD$d
  if (length(which(SVD$d<tol))>0){
    d<-SVD$d[-which(SVD$d<tol)]
  }
  n_effectiveColumns<-length(d)
  n.columns <- dim(X)[2]
  if (n_effectiveColumns<n.columns){
    if (n_effectiveColumns > 1){
      NewInv<-(SVD$u[,1:n_effectiveColumns])%*%diag(1/d)%*%t(SVD$v[,1:n_effectiveColumns])
    } else {
      NewInv<-1/d # this is in the extreme case in which all variables are perfectly correlated - not likely but allows the function to run
    }
  } else {
    NewInv<-solve(X)
  }
  return(NewInv)
}
