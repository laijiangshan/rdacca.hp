#' Internal functions for rdacca.hp
#'
#' Internal functions for rdacca.hp
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}

R2current=function(Y, current.comb, X,type)
{require(vegan)
    data <- data.frame(X[, current.comb])
	if(type=="RDA")
    gf <- RsquareAdj(rda(Y~., data = data))
	if(type=="CCA")
	gf <- RsquareAdj(cca(Y~., data = data))
    gf
}
