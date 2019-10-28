#' Hierarchical Partitioning for Redundancy Analysis and Canonical Correspondence Analysis
#'
#' This function calculates the individual contribution of each environmental variable for Redundancy Analysis and Canonical Correspondence Analysis,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' 
#'
#' @param  Y Community data matrix.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @param  pieplot a pieplot each variable is plotted expressed as percentage of total variation (pieplot="tv") or total explained variation (pieplot="tev").
#' @details This function calculates the individual contribution of each independent variable to goodness of fit measures on RDA and CCA,
#' applying the hierarchy algorithm of Chevan and Sutherland (1991) .
#' all combinations of N independent variable use the function allr2.
#' It takes the list of goodness of fit measures and, using the partition function,to return a simple table listing each variable, its individuals contribution (I). 
#' At this stage, the partition routine will not run for more than 12 independent variables. This function requires the vegan,hier.part package.

#' @return a list containing
#' @return \item{R2}{unadjusted R-squared for RDA or CCA  for overall model.}
#' @return \item{hp.R2}{the individual contribution for each variable (based on unadjusted R-squared).}
#' @return \item{adj.R2}{adjusted R-squared for RDA or CCA for overall model.}
#' @return \item{hp.adjR2}{the individual contribution for each variable (based on adjusted R-squared).}
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician 45:90~96

#' Chris Walsh and Ralph Mac Nally 2013. hier.part: Hierarchical Partitioning. R package version 1.0-4.https://CRAN.R-project.org/package=hier.part

#' @examples
#'require(ade4)
#'data(doubs)
#'spe<-doubs$fish
#'env<-doubs$env
#'#Remove empty site 8
#'spe<-spe[-8,]
#'env<-env[-8,]
#'#Remove the 'dfs' variable from the env data frame
#'env <- env[, -1]
#'# Hellinger-transform the species dataset
#'spe.hel <- decostand(spe, "hellinger")
#'#Forward selection for RDA using vegan's ordiR2step()
#'mod0 <- rda(spe.hel~ 1, env)  # Model with intercept only
#'mod1 <- rda(spe.hel~ ., env)
#'ordiR2step(mod0, mod1,direction = "forward")
#'#the remaining variables: alt,oxy,bdo
#'rdacca.hp(spe.hel,env[,c("alt","oxy","bdo")],pieplot = "tv",type="RDA")
#'rdacca.hp(spe.hel,env[,c("alt","oxy","bdo")],pieplot = "tev",type="RDA")
#'#Forward selection for CCA using vegan's ordistep()
#'mod0 <- cca(spe ~ 1, env)  # Model with intercept only
#'mod1 <- cca(spe ~., env)
#'ordistep(mod0, mod1,direction = "forward")
#'#the remaining variables: alt, oxy,bdo
#'rdacca.hp(spe,env[,c("alt","oxy","bdo")],pieplot = "tv",type="CCA")
#'rdacca.hp(spe,env[,c("alt","oxy","bdo")],pieplot = "tev",type="CCA")

#'data(mite)
#'data(mite.env)
#'mite.hel <- decostand(mite, "hellinger")
#Forward selection using vegan's ordiR2step()
#'mod0 <- rda(mite.hel ~ 1, mite.env)  # Model with intercept only
#'mod1 <- rda(mite.hel~ ., mite.env)
#'ordiR2step(mod0, mod1,direction = "forward")
#'#the remaining variables: WatrCont,Shrub, Substrate,Topo
#'rdacca.hp(mite.hel,mite.env[,-1],pieplot = "tv",type="RDA")
#'rdacca.hp(mite.hel,mite.env[,-1],pieplot = "tev",type="RDA")

#'#Forward selection using vegan's ordistep()
#'mod0 <- cca(mite ~ 1, mite.env)  # Model with intercept only
#'mod1 <- cca(mite ~., mite.env)
#'ordistep(mod0, mod1,direction = "forward")
#'#all variables is remaining
#'rdacca.hp(mite,mite.env,pieplot = "tv",type="CCA")
#'rdacca.hp(mite,mite.env,pieplot = "tev",type="CCA")



rdacca.hp=function (Y, X, pieplot="tv", type="RDA")
{
  require(scales)
  Env.num <- dim(X)[2]
  if (Env.num > 13)
    stop("Number of variables must be < 13 for current implementation",call. = FALSE)
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
      par(mfrow=c(1,2))
      par(mar=c(0.5,2.5,2.5,2.5)) 
      lbls<- c("unexplained",row.names(HP$I.perc)) 
      pct <- round(c(1-sum(HP$IJ$I),HP$IJ$I)*100,1)
      lbls <- paste(lbls, pct) # add percents to labels
      lbls <- paste(lbls,"%",sep="") # ad % to labels
      pie(pct,labels = lbls, main="% on Total variation\n (R2)",col=hue_pal()(length(lbls)),border = "white")
      lblsa<- c("unexplained",row.names(HPa$I.perc)) 
      pcta <- round(c(1-sum(HPa$IJ$I),HPa$IJ$I)*100,1)
      lblsa <- paste(lblsa, pcta) # add percents to labels
      lblsa <- paste(lblsa,"%",sep="") # ad % to labels
      pie(pcta,labels = lblsa, main="% on Total variation\n (adj.R2)",col=hue_pal()(length(lblsa)),border = "white")
      
    }
    if (pieplot=="tev") {
      par(mfrow=c(1,2))
      par(mar=c(0.5,2.5,2.5,2.5)) 
      lbls<- row.names(HP$I.perc) 
      pct <- round(HP$I.perc$I,1)
      lbls <- paste(lbls, pct) # add percents to labels
      lbls <- paste(lbls,"%",sep="") # ad % to labels
      pie(pct,labels = lbls,col=hue_pal()(length(lbls)), main="% on Total explained variation\n (R2)",border = "white")
      lblsa<- row.names(HPa$I.perc) 
      pcta <- round(HPa$I.perc$I,1)
      lblsa <- paste(lblsa, pcta) # add percents to labels
      lblsa <- paste(lblsa,"%",sep="") # ad % to labels
      pie(pcta,labels = lblsa,col=hue_pal()(length(lblsa)), main="% on Total explained variation\n (adj.R2)",border = "white")
    }
    
    par(op)
    
    list(R2=sum(HP$IJ[,"I"]), hp.R2 =HP$IJ["I"],adj.R2=sum(HPa$IJ[,"I"]),hp.adjR2=HPa$IJ["I"])
  }
}
