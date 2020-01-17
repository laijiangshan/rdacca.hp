#' Internal functions for rdacca.hp
#'
#' Internal functions for rdacca.hp
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}

#' @param  Y Response variables, typically of community data matrix.
#' @param  current.comb The orders of focus explanatory variables for r-squared
#' @param  X Explanatory variables, typically of environmental variables.
#' @param  type The type of constrained ordination: RDA or CCA, the default is "RDA".

R2current=function(Y, current.comb, X,type)
{#require(vegan)
    data <- data.frame(X[, current.comb])
	if(type=="RDA")
    gf <- vegan::RsquareAdj(rda(Y~., data = data))
	if(type=="CCA")
	gf <- vegan::RsquareAdj(cca(Y~., data = data))
    gf
}
