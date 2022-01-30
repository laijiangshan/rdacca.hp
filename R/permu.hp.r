#' Permutation Test of Hierarchical Partitioning for Canonical Analysis

#' @param  dv  Response variable, either a numeric vector or a matrix. If method="dbRDA", dv should be a "dist" matrix.
#' @param  iv Predictors (explanatory variable), either a data frame or a list of data frames. If it is a data frame, the relative importance of each column of the data frame will be evaluated; if it is a list, the relative importance of each element (matrix) will be evaluated.
#' @param  method Type of canonical analysis used for variation partitioning, should be a character string, either "RDA", "dbRDA" or "CCA", the default is "RDA". If the response variable (dv) is a numerical vector and method="RDA", the hierarchical and variation partitioning for the classical multiple regression is implemented.
#' @param  type The type of total explained variation, either "R2" or "adjR2", in which "R2" is unadjusted R-square and "adjR2" is adjusted R-square, the default is "adjR2". The adjusted R-square is calculated using Ezekiel's formula (Ezekiel 1930) for RDA and dbRDA, while permutation procedure is used for CCA (Peres-Neto et al. 2006). 
#' @param  scale Logical; If the columns of dv should be standardized to unit variance when method="RDA" is applied.
#' @param  add Logical; If a constant should be added to the non-diagonal values to euclidify dissimilarities (see dbrda function in vegan for details). Choice "lingoes" (or TRUE) uses the recommended method of Legendre & Anderson (1999: "method 1") and "cailliez" uses their "method 2". The argument has an effect only when method="dbRDA".
#' @param  sqrt.dist Logical, If the square root of dissimilarities should be taken. This often euclidifies dissimilarities. The argument has an effect only when method="dbRDA"(see dbrda function in vegan for details).
#' @param  n.perm An integer; Number of permutations for computing the adjusted R-square for CCA. The argument has an effect only when method="CCA".
#' @param  permutations An integer; Number of permutations for computing p value of individual contribution for the randomized dataset.

#' @details This function is a permutation test of hierarchical partitioning for canonical analysis. It returns a matrix of I values (the individual contribution towards total explained variation) for all values from permutations randomizations. For each permutation, the values in each variable (i.e each column of iv) are randomized independently, and rdacca.hp is run on the randomized iv. As well as the randomized I matrix, the function returns a summary table listing the observed I values, the p value of I for the randomized dataset.

#' @return a data.frame containing a summary table listing the observed individual contribution, the p value of individual contribution for the randomized dataset

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}


#' @export
#'@examples
#'library(vegan)
#'data(mite)
#'data(mite.env)
#'#Hellinger-transform the species dataset for RDA
#'mite.hel <- decostand(mite, "hellinger")
#'permu.hp(mite.hel,mite.env,method="RDA",type="adjR2",permutations=10)


permu.hp=function(dv, iv, method=c("RDA","dbRDA","CCA") ,type=c("adjR2","R2"),scale=FALSE,add = FALSE, sqrt.dist = FALSE,n.perm=1000,permutations=1000)
{
cat("\nPlease wait: running", permutations-1, "permutations \n") 
obs <- rdacca.hp(dv,iv,method=method,type=type,scale=scale,add = add, sqrt.dist = sqrt.dist,n.perm=n.perm)
r2q <- obs$Hier.part[,3]

if(is.data.frame(iv))
{
n <- dim(iv)[1]
nvar <- dim(iv)[2]
for(i in 1:permutations-1)
{
newiv<-iv
for(j in 1:nvar)
{perms <- sample(1:n,n)
 newiv[,j] <-iv[,j][perms]
 }

 row.names(newiv)<-1:n 
 simu=rdacca.hp(dv,newiv,method=method,type=type,scale=scale,add = add, sqrt.dist = sqrt.dist,n.perm=n.perm)
 r2q=cbind(r2q,simu$Hier.part[,3])
}
}

else
{
n <- dim(iv[[1]])[1]
nvar  <-  length(iv)
for(i in 1:permutations-1)
{perms <- sample(1:n,n)
 newiv <- list()
 for(j in 1:nvar)
 {
 newiv[[j]] <- data.frame(iv[[j]][perms, ])
 row.names(newiv[[j]]) <- 1:n 
}
 names(newiv) <- names(iv)
 simu=rdacca.hp(dv,newiv,method=method,type=type,scale=scale,add = add, sqrt.dist = sqrt.dist,n.perm=n.perm)
 r2q=cbind(r2q,simu$Hier.part[,3])
}
}

Signi=function(x)
{pval=round(1-ecdf(x)(x[1])+1/(permutations+1),nchar(permutations))
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
result=data.frame(obs$Hier.part[,3],Pr=p.R2)
colnames(result)<-c("Individual","Pr(>I)")
return(result)
}

