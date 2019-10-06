#'Permutation Test for Hierarchical Partitioning for Redundancy Analysis and Canonical Correspondence Analysis 
#'
#' This function performs permutation test for the individual explain percentage of each environmental variable for Redundancy Analysis and Canonical Correspondence Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#' @param  Y Community data matrix for RDA and CCA,or a vector for general linear regression model.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @param  permutations the number of permutations required
#' @return a list containing
#' @return \item{R2}{a dataframe for unadjusted R-squared for individual environmental variables and p-value.}
#' @return \item{adjR2}{a dataframe for adjusted R-squared for individual environmental variables and p-value.}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96
#' @examples
#'require(vegan)
#'data(varespec)
#'data(varechem)
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="RDA")
#'permu.hp(varespec,varechem[,c("Al","P","K")],type="RDA",permutations=999)
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
 simu=rdacca.hp(newy,newx,type=type,pieplot ="")
 r2q=cbind(r2q,simu$hp.R2)
 ar2q=cbind(ar2q,simu$hp.adjR2)
}

Signi=function(x)
{1-ecdf(x)(x[1])+1/(permutations+1)}

p.R2=apply(r2q,1,Signi)
p.adjR2=apply(ar2q,1,Signi)
return(list(R2=data.frame(R2=obs$hp.R2,Pr=p.R2),adjR2=data.frame(adjR2=obs$hp.adjR2,Pr=p.adjR2)))
}

