#' Hierarchical Partitioning for Canonical Correspondence Analysis and Redundancy Analysis
#'
#' This function calculates the individual explain percentage of each environmental variable for Canonical Correspondence Analysis and Redundancy Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#'
#' @param  Y Community data matrix.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @param  pieplot a pieplot each variable is plotted expressed as percentage of total variation (pieplot="tv") or total explained variation (pieplot="tev").
#' @return a list containing
#' @return \item{R2}{unadjusted R-squared for RDA or CCA with all environmental variables.}
#' @return \item{hp.R2}{the individual contribution for each variable (based on unadjusted R-squared).}
#' @return \item{adj.R2}{adjusted R-squared for RDA or CCA with all environmental variables.}
#' @return \item{hp.adjR2}{the individual contribution for each variable (based on adjusted R-squared).}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96
#' @examples
#'require(vegan)
#'data(varespec)
#'data(varechem)
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="RDA")
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tev",type="RDA")
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tv",type="CCA")
#'rdacca.hp(varespec,varechem[,c("Al","P","K")],pieplot = "tev",type="CCA")


rdacca.hp=function (Y, X,type="RDA", pieplot = "tv")
{
    Env.num <- dim(X)[2]
    if (Env.num > 13)
        stop("Number of variables must be < 13 for current implementation",call. = FALSE)
    else {
	    hpmodel<- allR2(Y,X,type)
        gfs <- hpmodel$gfs
        HP <- partition.rda(gfs, Env.num, var.names = names(data.frame(X)))
        gfsa <- hpmodel$gfsa
        HPa <- partition.rda(gfsa, Env.num, var.names = names(data.frame(X)))
       par(mfrow=c(1,2))
	   par(mar=c(0.2,0.2,3,0.8))
	   if(type=="RDA")
            rsq <- RsquareAdj(rda(Y~., data = X))
       if(type=="CCA")
             rsq <- RsquareAdj(cca(Y~., data = X))
	   if (pieplot=="tv") {
           lbls<- c("unexplained",row.names(HP$I.perc)) 
		   pct <- round(c(1-sum(HP$IJ$I),HP$IJ$I)*100,1)
           lbls <- paste(lbls, pct) # add percents to labels
           lbls <- paste(lbls,"%",sep="") # ad % to labels
		   pie(pct,labels = lbls, main="% on Total variation\n (R2)")
           lblsa<- c("unexplained",row.names(HPa$I.perc)) 
		   pcta <- round(c(1-sum(HPa$IJ$I),HPa$IJ$I)*100,1)
           lblsa <- paste(lblsa, pcta) # add percents to labels
           lblsa <- paste(lblsa,"%",sep="") # ad % to labels
		   pie(pcta,labels = lblsa, main="% on Total variation\n (adj.R2)")
          
		}
		if (pieplot=="tev") {
           lbls<- row.names(HP$I.perc) 
		   pct <- round(HP$I.perc$I,1)
           lbls <- paste(lbls, pct) # add percents to labels
           lbls <- paste(lbls,"%",sep="") # ad % to labels
           pie(pct,labels = lbls,col=rainbow(length(lbls)), main="Total explained variation\n (R2)")
           lblsa<- row.names(HPa$I.perc) 
		   pcta <- round(HPa$I.perc$I,1)
           lblsa <- paste(lblsa, pcta) # add percents to labels
           lblsa <- paste(lblsa,"%",sep="") # ad % to labels
		 pie(pcta,labels = lblsa,col=rainbow(length(lblsa)), main="Total explained variation\n (adj.R2)")
		}
		
         list(R2=rsq$r.squared, hp.R2 =HP$IJ["I"],adj.R2=rsq$adj.r.squared,hp.adjR2=HPa$IJ["I"])
    }
}

