#' All Combinations of a Hierarchy of Models of n Variables
#'
#' Internal functions for rdacca.hp. Lists a matrix of combinations of 1 to n variables in ascending order
#' @param  n An integer: the number of variable.

#' @details 
#' Lists hierarchy of all possible combinations of n variables in ascending order, 
#' starting with 1 variable, then all combinations of 2 variables, and so on until 
#' the one combination with all n variables. This function is used by allR2 to 
#' structure the models required for hierarchical partioning. This function requires 
#' the gtools package in the gregmisc bundle.

#' @return  A matrix with all combinations of a hierarchy of models of n variables
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @references
#' Chevan, A. and Sutherland, M. 1991. Hierarchical Partitioning. The American Statistician.45:90-96
#' @references
#' Chris Walsh and Ralph Mac Nally 2013. hier.part: Hierarchical Partitioning. R package version 1.0-4.https://CRAN.R-project.org/package=hier.part


hpmatrix=function(n) 
{   #require(gtools)
    x <- cbind(1:n, matrix(0, n, n - 1))
    for (i in 2:n) {
        nc <- dim(combinations(n, i, 1:n))[1]
        x <- rbind(x, cbind(combinations(n, i, 1:n), matrix(0,nc, n - i)))
    }
x
}

