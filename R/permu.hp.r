#'Permutation Test for Hierarchical Partitioning for Redundancy Analysis and Canonical Correspondence Analysis 
#'
#' This function performs permutation test for the individual contribution of each environmental variable for Redundancy Analysis and Canonical Correspondence Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#' @param  Y Community data matrix for RDA and CCA,or a vector for general linear regression model.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @param  permutations the number of permutations required

#' @details This function perform significance test for individual R2 of each independent variable via random permutation of the
#' rows of the Y data matrix. the test statistic is defined as individual R2. 
#' At this stage, permu.hp  the defined permutation scheme


#' @return a list containing
#' @return \item{R2}{unadjusted R-squared for RDA or CCA  for overall model.}
#' @return \item{hp.R2}{a dataframe for unadjusted R-squared for individual environmental variables and p-value.}
#' @return \item{adj.R2}{adjusted R-squared for RDA or CCA for overall model.}
#' @return \item{hp.adjR2}{a dataframe for adjusted R-squared for individual environmental variables and p-value.}

#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96
#' Chris Walsh and Ralph Mac Nally 2013. hier.part: Hierarchical Partitioning. R package version 1.0-4.https://CRAN.R-project.org/package=hier.part

#' @examples
#'require(vegan)
#'data(varespec)
#'data(varechem)
#'rdacca.hp(varespec,varechem[,c("Al","K","N","Baresoil")],pieplot = "tv",type="RDA")
#'permu.hp(varespec,varechem[,c("Al","K","N","Baresoil")],type="RDA",permutations=999)
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="CCA")
#'permu.hp(varespec,varechem[,c("Al","P","K")],type="CCA",permutations=999)

permu.hp=function(Y,X,type="RDA",permutations=999)
{require(permute)
obs=rdacca.hp(Y,X,type=type,pieplot ="")
Y=data.frame(Y)
n=dim(Y)[1]
r2q=obs$hp.R2
ar2q=obs$hp.adjR2
for(i in 1:permutations)
{newy=Y[shuffle(n),]
 newx=X[shuffle(n),]
 simu=rdacca.hp(newy,X,type=type,pieplot ="")
 r2q=cbind(r2q,simu$hp.R2)
 ar2q=cbind(ar2q,simu$hp.adjR2)
}

Signi=function(x)
{1-ecdf(x)(x[1])+1/(permutations+1)}

p.R2=apply(r2q,1,Signi)
p.adjR2=apply(ar2q,1,Signi)
return(list(R2=obs$R2,hp.R2=data.frame(obs$hp.R2,Pr=p.R2),adj.R2=obs$adj.R2,hp.adjR2=data.frame(adjR2=obs$hp.adjR2,Pr=p.adjR2)))
}

