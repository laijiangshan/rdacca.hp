partition.rda=function (gfs, pcan, var.names = NULL)
{require(hier.part)
    if (pcan > 12)
        stop("Number of variables must be < 13 for current implementation",call. = FALSE)
    else if (pcan > 9)
        warning("rdaenvpart produces a rounding error if number of variables >9\nSee documentation.",
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
        IJ <- cbind(IJ, Total = IJ$I + IJ$J)
        list(gfs = gfs, IJ = IJ, I.perc = I.perc)
    }
}


