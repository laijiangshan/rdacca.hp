#' R2 for RDA or CCA Hierarchy
#'
#' Internal functions for rdacca.hp'. Calculates R2 for RDA or CCA of a single explanatory variable to all combinations of N explanatory variables.
#' @param  Y Response variables, typically of community data matrix.
#' @param  X Explanatory variables, typically of environmental variables.
#' @param  type The type of constrained ordination: RDA or CCA, the default is "RDA".

#' @details 
#' This function calculates R2 and adjusted R2 for the entire hierarchy of RDA or CCA using all combinations of N 
#' explanatory variables, and returns them as an ordered list ready for input into the function partition.rda. 
#' This function requires the gtools package in the gregmisc bundle.

#' @return a list containing
#' @return \code{gfs} A data frame listing all combinations of explanatory variables in the first column in ascending order, and the corresponding R2 for the model using those variables.
#' @return \code{gfsa} A data frame listing all combinations of explanatory variables in the first column in ascending order, and the corresponding adjusted R2 for the model using those variables.


#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan Albert and Sutherland Michael. 1991. Hierarchical Partitioning. The American Statistician.45:90~96
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
#'allR2(spe.hel,env[,c("alt","oxy","bdo")],type="RDA")
#'allR2(spe,env[,c("alt","oxy","bdo")],type="CCA")



allR2=function (Y,X,type)
{
    Env.num <- dim(X)[2]
    n <- (2^Env.num) - 1
    combs <- hpmatrix(Env.num)
    gfs <- 0
    gfsa<- 0
    for (i in 1:n) {
        current.comb <- as.vector(combs[i, ][combs[i, ] > 0])
        combn <- paste(names(data.frame(X)[current.comb]),
            "", collapse = "")
        new.line <- R2current(Y,current.comb,X,type)
        gfs <- c(gfs, new.line$r.squared)
		gfsa <- c(gfsa, new.line$adj.r.squared)
    }
    list(gfs=gfs,gfsa=gfsa)
}

