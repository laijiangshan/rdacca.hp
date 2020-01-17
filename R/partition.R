#' Internal functions for rdacca.hp
#'
#' Internal functions for rdacca.hp
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}

#' @param  gfs a vector of r-squared from a hierarchy of RDA or CCA based on pcan variables in ascending order (as produced by function hpmatrix, but also including the null model as the first element)
#' @param  pcan the number of variables from which the hierarchy was constructed (maximum =9)
#' @param  var.names an array of pcan variable names, if required



partition.rda=function (gfs, pcan, var.names = NULL)
{#require(hier.part)
    if (pcan > 9)
        stop("Number of explanatory variables must be < 10 for current implementation",call. = FALSE)
    else if (pcan > 9)
        warning("rdacca.hp produces a rounding error if number of explanatory variables >9",
            call. = FALSE)
    {
        n <- 2^pcan
        theta <- gfs[1]
        fin <- gfs[2:n]
      
        len <- length(fin)
        IJ <- vector("numeric", pcan * 2)
        storage.mode(pcan) <- "integer"
        storage.mode(len) <- "integer"
        storage.mode(theta) <- "double"
        storage.mode(fin) <- "double"
        storage.mode(IJ) <- "double"
        IJ <- .C("hierpart", pcan, len, theta, fin, IJ = IJ,
            PACKAGE = "hier.part")$IJ
        IJ <- array(IJ, dim = c(pcan, 2))
        IJ <- data.frame(t(data.frame(t(IJ), row.names = c("I",
            "J"))), row.names = var.names)
        IJ.perc <- IJ * 100/sum(IJ)
        I <- data.frame(I = IJ[, 1], row.names = var.names)
        I.perc <- I * 100/sum(I)
        IJ <- data.frame(Total = IJ$I + IJ$J,J=IJ$J,I=IJ$I,row.names = var.names)
        list(gfs = gfs, IJ = IJ, I.perc = I.perc)
    }
}


