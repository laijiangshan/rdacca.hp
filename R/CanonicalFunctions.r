#'Internal function for rdacca.hp() to calculate the Rsquared, adjusted Rsquared, Pseudo-F value, AIC and AICc for canonical analysis (RDA, db-RDA and CCA)
#' @param  dv Response variables. if method="dbRDA", dv is the "dist" matrix.
#' @param  iv Explanatory variables, typically of environmental variables.
#' @param  method The type of canonical analysis: RDA, dbRDA or CCA, the default is "RDA".
#' @param  n.perm Number of permutations to use when computing the adjusted R-squared for a CCA.
#' @keywords internal
Canonical.Rsq <- function(dv,iv,method="RDA",n.perm=1000){ 
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
  
standardize_w <- function(X,w){
    ones <- rep(1,length(w))
    Xc <- X - ones %*% t(w)%*% X
    Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w)) 
  } 
  n <- nrow(iv)
  n.iv <- qr(iv)$rank
  if (method=="RDA"){
    SSY <- sum(diag(t(dv)%*%dv))
    Yhat<-iv%*%GMI(t(iv)%*%iv)%*%t(iv)%*%dv
    SSYhat <- sum(diag(t(Yhat)%*%Yhat))
    e <- dv-Yhat
  }
  if (method=="dbRDA"){
    D <- as.matrix(dv)
    G <- -0.5 * (diag(n) - matrix(1,n,n)/n) %*% (D^2) %*% (diag(n) - matrix(1,n,n)/n)
    SSY <- sum(diag(G))
    H <- iv %*% GMI(t(iv) %*% iv) %*% t(iv)
    SSYhat <- sum(diag(H %*% G %*% H))
    Yhat <- (H %*% G %*% H)
    e <- G-Yhat
  }
  if (method=="CCA"){
    RSq.CCA <- function(dv,iv){
      n.dv <- ncol(dv)
      TotalSum <- sum(dv)
      Pij <- as.matrix(dv/TotalSum)
      SumCols <- apply(dv,2,sum)/TotalSum
      SumRows <- apply(dv,1,sum)/TotalSum
      Q <- diag(1/sqrt(SumRows))%*%Pij%*%diag(1./sqrt(SumCols))-diag(sqrt(SumRows))%*%matrix(1,n,n.dv)%*%diag(sqrt(SumCols))
      SSY <- sum(diag(t(Q) %*% Q))
      iv.std.w <- standardize_w(iv,SumRows)
      B <- GMI(t(iv.std.w) %*% diag(SumRows) %*% iv.std.w) %*% t(iv.std.w) %*% diag(sqrt(SumRows)) %*% Q
      Yhat <- diag(sqrt(SumRows))%*%iv.std.w%*%B
      SSYhat <- sum(diag(t(Yhat) %*% Yhat))
      RSq <- SSYhat/SSY
      e <- Q-Yhat
      F <- (SSYhat/(SSY-SSYhat))*(n-n.iv-1)/n.iv
      return(list(RSq=RSq,F=F,e=e,Yhat=Yhat))
    }
    res.CCA <- RSq.CCA(dv,iv)
    RSq <- res.CCA$RSq
    F <- res.CCA$F
    RSq.perm <- matrix(0,n.perm,1)
    for (i in 1:n.perm){
      RSq.perm[i] <- RSq.CCA(dv,iv[sample(n),])$RSq
    }
    RSq.perm <- mean(RSq.perm)
    RSqAdj <- 1-1/(1-RSq.perm)*(1-RSq)
    e <- res.CCA$e
	Yhat<- res.CCA$Yhat
  } #CCA

  if ((method=="RDA") | (method=="dbRDA")){ 
    RSq <- SSYhat/SSY 
    RSqAdj <- 1-((1-RSq)*(n-1)/(n-1-n.iv))
    F <- (SSYhat/(SSY-SSYhat))*(n-n.iv-1)/n.iv}
	
	AIC<-2*n.iv+n*log((1-RSq)/n)#add AIC and AICc by Jiangshan Lai
	AICc<-(AIC+2*n.iv*(n.iv+1))/(n-n.iv-1)
  Stats <- list(unadj=RSq,adj=RSqAdj,F=F,AIC=AIC,AICc=AICc,e=e,Yhat=Yhat)
  return(Stats)
}
