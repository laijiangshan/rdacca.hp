#' Hierarchical Partitioning for Redundancy Analysis and Canonical Correspondence Analysis
#'
#' Partitions variation for each explanatory variable in Redundancy Analysis and Canonical Correspondence Analysis
 
#' @param  Y Response variables, typically of community data matrix.
#' @param  X Explanatory variables, typically of environmental variables.
#' @param  type The type of constrained ordination: RDA or CCA, the default is "RDA".
#' @param  pieplot A pieplot where each variable is plotted expressed as the percentage of its independent R-squared in the total variation (pieplot="tv") or the total explained variation (pieplot="tev").
#' @param  trace If TRUE,R-squared of all combinations (from all.R2 function),joint contribution of each explanatory variable are printed. The default is FALSE.
#' @details This function calculates the independent contribution of each explanatory variable to explained variation (R-squared) on RDA and CCA,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991). The algorithm is that all joint R-squared will be decomposed into equal fractions by number 
#' of involved explanatory variables and average assigned to these variables. Independent R-squared of each variable will be the sum of assigned R-squared from joint R-squared and unique R-squared. 
#' All combinations of N explanatory variable use the function allr2().
#' It takes the list of R-squared and, using the partition function,to return a simple table listing each variable, its independent contribution (I). 
#' At this stage, the partition routine will not run for more than 9 explanatory variables, due to limitation of computation power and a rounding error for analyses. This code of function is dependent on heir.part packages (Chris and Ralph 2013).

#' @return a list containing
#' @return \item{R2}{Unadjusted R-squared of RDA or CCA for overall model.}
#' @return \item{all.R2}{If trace=TRUE,a vector listing the corresponding R-squared for the model using all combinations of each explanatory variables in ascending order.}
#' @return \item{hp.R2}{The independent contribution for each explanatory variable (based on unadjusted R-squared).If trace=TRUE, the joint and total R-squared of each explanatory variables is listed.}
#' @return \item{adj.R2}{Adjusted R-squared of RDA or CCA for overall model.}
#' @return \item{all.adjR2}{If trace=TRUE, a vector listing the corresponding adjusted R-squared for the model using all combinations of each explanatory variables in ascending order.}
#' @return \item{hp.adjR2}{The independent contribution for each explanatory variable (based on adjusted R-squared).If trace=TRUE, the joint and total adjusted R-squared of each explanatory variables is listed.}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan Albert and Sutherland Michael. 1991. Hierarchical Partitioning. The American Statistician.45:90-96
#' @references
#' Chris Walsh and Ralph Mac Nally 2013. hier.part: Hierarchical Partitioning. R package version 1.0-4.https://CRAN.R-project.org/package=hier.part

#'@examples
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
#'rdacca.hp(spe.hel,env[,c("alt","oxy","bdo")],pieplot = "tv",type="RDA",trace=TRUE)
#'rdacca.hp(spe,env[,c("alt","oxy","bdo")],pieplot = "tev",type="CCA")

#'data(mite)
#'data(mite.env)
#'#Hellinger-transform the species dataset for RDA
#'mite.hel <- decostand(mite, "hellinger")
#'rdacca.hp(mite.hel,mite.env,pieplot = "tv",type= "RDA",trace=TRUE)
#'rdacca.hp(mite,mite.env,pieplot = "tev",type= "CCA")

rdacca.hp=function (Y, X, pieplot = "tv", type = "RDA", trace = FALSE) 
{
  Env.num <- dim(X)[2]
  if (Env.num > 10)
    stop("Number of explanatory variables must be < 10 for current implementation",call. = FALSE)
  else {
    hpmodel<- allR2(Y,X,type)
    gfs <- hpmodel$gfs
    HP <- partition.rda(gfs, Env.num, var.names = names(data.frame(X)))
    gfsa <- hpmodel$gfsa
    gfsa[gfsa<0]=0
    HPa <- partition.rda(gfsa, Env.num, var.names = names(data.frame(X)))
    HPa$IJ$I[HPa$IJ$I<0]=0
    HPa$I.perc[HPa$I.perc<0]=0
    
    op <- par(no.readonly = TRUE) 
    
    if (pieplot=="tv") {
      par(mfrow=c(2,1))
      par(mar=c(0.5,1.5,2.5,1.5)) 
      lbls<- c("unexplained",row.names(HP$I.perc)) 
      pct <- round(c(1-sum(HP$IJ$I),HP$IJ$I)*100,1)
      lbls <- paste(lbls, pct) # add percents to labels
      lbls <- paste(lbls,"%",sep="") # ad % to labels
      pie(pct,labels = lbls, main="% on total variation (R2)",col=hue_pal()(length(lbls)),border = "white")
      lblsa<- c("unexplained",row.names(HPa$I.perc)) 
      pcta <- round(c(1-sum(HPa$IJ$I),HPa$IJ$I)*100,1)
      lblsa <- paste(lblsa, pcta) # add percents to labels
      lblsa <- paste(lblsa,"%",sep="") # ad % to labels
      pie(pcta,labels = lblsa, main="% on total variation (adj.R2)",col=hue_pal()(length(lblsa)),border = "white")
      
    }
    if (pieplot=="tev") {
      par(mfrow=c(2,1))
      par(mar=c(0.5,1.5,2.5,1.5)) 
      lbls<- row.names(HP$I.perc) 
      pct <- round(HP$I.perc$I,1)
      lbls <- paste(lbls, pct) # add percents to labels
      lbls <- paste(lbls,"%",sep="") # ad % to labels
      pie(pct,labels = lbls,col=hue_pal()(length(lbls)), main="% on total explained variation (R2)",border = "white")
      lblsa<- row.names(HPa$I.perc) 
      pcta <- round(HPa$I.perc$I,1)
      lblsa <- paste(lblsa, pcta) # add percents to labels
      lblsa <- paste(lblsa,"%",sep="") # ad % to labels
      pie(pcta,labels = lblsa,col=hue_pal()(length(lblsa)), main="% on total explained variation (adj.R2)",border = "white")
    }
    
    par(op)
    if (trace == FALSE) 
    return(list(R2 = sum(HP$IJ[, "I"]), hp.R2 = HP$IJ["I"], 
                adj.R2 = sum(HPa$IJ[, "I"]), hp.adjR2 = HPa$IJ["I"]))
    if (trace == TRUE) 
            return(list(R2 = sum(HP$IJ[, "I"]), all.R2 = gfs, hp.R2 = HP$IJ, 
                adj.R2 = sum(HPa$IJ[, "I"]), all.adjR2 = gfsa, hp.adjR2 = HPa$IJ))
    }
}
