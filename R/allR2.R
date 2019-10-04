##
#' Goodness of Fit Measures for RDA or CCA Hierarchy
#' Calculates goodness of fit measures for RDA or CCA of a single dependent variable to all combinations of N independent variable
#' @param  Y Community data matrix.
#' @param  X Constraining matrix less than 12 columns, typically of environmental variables.
#' @param  type the Constrained ordination: RDA or CCA, default "RDA"
#' @return \code{gfs} a data frame listing all combinations of independent variables in the first column in ascending order, and the corresponding goodness of fit measure for the model using those variables
#' @export

allR2=function (Y,X,type)
{
    Env.num <- dim(X)[2]
    n <- (2^Env.num) - 1
    combs <- hpmatrix(Env.num)
    gfs <- 0

    for (i in 1:n) {
        current.comb <- as.vector(combs[i, ][combs[i, ] > 0])
        combn <- paste(names(data.frame(X)[current.comb]),
            "", collapse = "")
        new.line <- R2current(Y, current.comb, X,type)
        gfs <- c(gfs, new.line)
    }
    gfs
}
