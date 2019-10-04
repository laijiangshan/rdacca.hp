R2current=function(Y, current.comb, X,type)
{require(vegan)
    data <- data.frame(X[, current.comb])
	if(type=="RDA")
    gf <- RsquareAdj(rda(Y~., data = data))$r.squared
	if(type=="CCA")
	gf <- RsquareAdj(cca(Y~., data = data))$r.squared
    gf
}
