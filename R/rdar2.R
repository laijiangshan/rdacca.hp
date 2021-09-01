#' An internal function used in rdacca.hp(): Returns only the raw Rsquare and the rank of constraints in RDA, modify "simpleRDA2" from vegan package
#' @param  dv Response matrix.
#' @param  iv Explanatory matrix
#' @keywords internal

rdar2 <-
    function (dv, iv)
{
    Q <- qr(iv, tol=1e-6)
    Yfit.X <- qr.fitted(Q, dv)
    SS <- sum(Yfit.X^2)
    SS.Y <- sum(dv^2)
    R2 <- SS/SS.Y
    adjR2 <- RsquareAdj(R2, nrow(dv), Q$rank)
    list(Rsquare = R2, RsquareAdj = adjR2)
}

