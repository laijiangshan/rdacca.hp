% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CanonicalFunctions.r
\name{Canonical.Rsq}
\alias{Canonical.Rsq}
\title{Internal function for rdacca.hp() to calculate the Rsquared, adjusted Rsquared, Pseudo-F value, AIC and AICc for canonical analysis (RDA, db-RDA and CCA)}
\usage{
Canonical.Rsq(dv, iv, method = "RDA", n.perm = 1000)
}
\arguments{
\item{dv}{Response variables. if method="dbRDA", dv is the "dist" matrix.}

\item{iv}{Explanatory variables, typically of environmental variables.}

\item{method}{The type of canonical analysis: RDA, dbRDA or CCA, the default is "RDA".}

\item{n.perm}{Number of permutations to use when computing the adjusted R-squared for a CCA.}
}
\description{
Internal function for rdacca.hp() to calculate the Rsquared, adjusted Rsquared, Pseudo-F value, AIC and AICc for canonical analysis (RDA, db-RDA and CCA)
}
\keyword{internal}
