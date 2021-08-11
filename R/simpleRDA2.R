#' An internal function used in rdacca.hp(): Returns only the raw Rsquare and the rank of constraints in RDA, modify from vegan package
#' @param  Y Response matrix.
#' @param  X Explanatory matrix
#' @param  SS.Y Sum of square of response matrix (Y).
#' @keywords internal

`simpleRDA2` <-
    function (Y, X, SS.Y, ...)
{
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    if (missing(SS.Y)) SS.Y <- sum(Y^2)
    Rsquare <- SS/SS.Y
    R2adj <- RsquareAdj(Rsquare, nrow(Y), Q$rank)
    list(Rsquare = Rsquare, RsquareAdj = R2adj, m = Q$rank)
}

