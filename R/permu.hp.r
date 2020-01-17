#'Permutation Test for Hierarchical Partitioning for Redundancy Analysis and Canonical Correspondence Analysis 
#'
#'Randomizes row in Y and X and recalculates rdacca.hp permuation times

#' @param  Y Response variables,typically of community data matrix.
#' @param  X Explanatory variables, typically of environmental variables.
#' @param  type The type of constrained ordination: RDA or CCA, the default is "RDA".
#' @param  permutations The number of permutations required.

#' @details 
#' This function is a permutation routine for the rdacca.hp function which returns a matrix of I values (the independent contribution 
#' towards explained variation in RDA or CCA) for each explanatory variables. For each permutation, the values in each row (both Y and X) 
#' are randomized independently, and rdacca.hp is run on the randomized Y and X.  
#' The function returns a summary table listing the observed I values, the pvalue for each I value.

#' @return a list containing
#' @return \item{R2}{Unadjusted R-squared of RDA or CCA for overall model.}
#' @return \item{hp.R2}{The independent contribution for each explanatory variable (based on unadjusted  R-squared) and p-value.}
#' @return \item{adj.R2}{Adjusted  R-squared of RDA or CCA for overall model.}
#' @return \item{hp.adjR2}{The independent contribution for each explanatory variable (based on adjusted  R-squared) and p-value.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan Albert and Sutherland Michael. 1991. Hierarchical Partitioning. The American Statistician.45:90-96
#' @references
#' Chris Walsh and Ralph Mac Nally 2013. hier.part: Hierarchical Partitioning. R package version 1.0-4.https://CRAN.R-project.org/package=hier.part

#' @examples
#'require(ade4)
#'data(doubs)
#'spe<-doubs$fish
#'env<-doubs$env
#'#Remove empty site 8
#'spe<-spe[-8,]
#'env<-env[-8,]
#'#Hellinger-transform the species dataset for RDA
#'spe.hel <- decostand(spe, "hellinger")
#'#select three variables: alt,oxy,and bdo as main explanatory set, via forward selection.
#'permu.hp(spe.hel,env[,c("alt","oxy","bdo")],type="RDA",permutations=99)
#'permu.hp(spe,env[,c("alt","oxy","bdo")],type="CCA",permutations=99)

permu.hp=function(Y,X,type="RDA",permutations=999)
{
cat("\nPlease wait: running", permutations, "permutations \n")
#require(permute)
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
{pval=1-ecdf(x)(x[1])+1/(permutations+1)
if (pval <= 0.001) {
		return(noquote(paste(pval,"***",sep=" ")))
		}
	if (pval <= 0.01) {
		return(noquote(paste(pval," **",sep=" ")))
		}		
	if (pval <= 0.05) {
		return(noquote(paste(pval,"  *",sep=" ")))
		}
	else {return(noquote(paste(pval,"   ",sep=" ")))
		}
		}
		

p.R2=apply(r2q,1,Signi)
p.adjR2=apply(ar2q,1,Signi)
hp.R2=data.frame(obs$hp.R2,Pr=p.R2)
colnames(hp.R2)[2]="Pr(>I)"
hp.adjR2=data.frame(adjR2=obs$hp.adjR2,Pr=p.adjR2)
colnames(hp.adjR2)[2]="Pr(>I)"
return(list(R2=obs$R2,hp.R2=hp.R2,adj.R2=obs$adj.R2,hp.adjR2=hp.adjR2))
}

