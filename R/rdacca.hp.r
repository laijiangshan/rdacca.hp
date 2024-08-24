#' Hierarchical Partitioning for Canonical Analysis

#' @param  dv Response variables. if method="dbRDA", dv is the "dist" matrix.
#' @param  iv Explanatory variables, typically of environmental variables.
#' @param  method The type of canonical analysis: RDA, dbRDA or CCA, the default is "RDA".
#' @param  type The type of total explained variation: "adjR2" is adjusted R-squared and "R2" for unadjusted R-squared, the default is "adjR2".
#' @param  trace logical value, if TRUE, the vaules of commonality (2^N-1for N explanatory variables) are outputed,the default is FALSE.

#' @details This function calculates the independent contribution of each explanatory variable to explained variation (R-squared) on canonical analysis (RDA,CCA and dbRDA),
#' applying the hierarchy algorithm of Chevan, A. and Sutherland, M. 1991 Hierarchical Partitioning.The American Statistician, 90--96. DOI: 10.1080/00031305.1991.10475776. 
#' Under the idea of hierarchy algorithm, the shared R2 can be divided into equal components by number of involved variables, and then allocated equally to these variables as joint effects.
#' The independent contribution of each explanatory variable is the sum of all its allocated common R2 and its unique R2. The order of importance of explanatory variables are determined by their independent contributions. 

#' @return a list containing
#' @return \item{Method_Type}{The type of canonical analysis and the type of total explained variation.}
#' @return \item{R.squared}{Total explained variation.}
#' @return \item{Commonality}{If trace=TRUE,a mtrix listing tha value and percentage of all commonality (2^N-1 for N explanatory variables).}
#' @return \item{Var.part}{A matrix listing independent effect and its percentage to total explained variation for each explanatory variable.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Pedro Peres-Neto} \email{pedro.peres-neto@concordia.ca}

#' @export
#'@examples
#'library(vegan)
#'data(mite)
#'data(mite.env)
#'#Hellinger-transform the species dataset for RDA to deal with the "double zero" problem
#'mite.hel <- decostand(mite, "hellinger")
#'rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#'rdacca.hp(vegdist(mite),mite.env,method="dbRDA",type="adjR2")
#'rdacca.hp(mite,mite.env,method="CCA",type="adjR2")


rdacca.hp <- function (dv,iv,method=c("RDA","dbRDA","CCA"),type=c("adjR2","R2"),trace = FALSE)
{
  if(sum(is.na(dv))>=1|sum(is.na(iv))>=1)
  {cat("Error: NA is not allowed in this analysis")}
  else
  {method <- method[1]
  type <- type[1]
  if(inherits(dv, "dist"))
  {method <- "dbRDA"}
  if(method=="dbRDA"){
    if(!inherits(dv, "dist"))
      return("dv should be a 'dist' matrix for dbRDA method")
  }

  if(method=="RDA")
  {dv<-scale(dv,scale=F)}

  iv <- data.frame(iv)
  ivname <- colnames(iv)
  iv.name <- ivname
  nvar <- dim(iv)[2]
  if (nvar < 2)
    stop("Analysis not conducted. Insufficient number of explanatory variables.")

  totalN <- 2^nvar - 1
  binarymx <- matrix(0, nvar, totalN)
  for (i in 1:totalN) {
    binarymx <- creatbin(i, binarymx)
  }

  commonM <- matrix(nrow = totalN, ncol = 3)
  for (i in 1:totalN) {
    tmp <- iv[as.logical(binarymx[, i])]
    tmp.design <- as.matrix(stats::model.matrix(~ ., tmp))
    tmp.design.ct <- as.matrix(data.frame(scale(tmp.design,scale=FALSE))[-1])
    if(method=="CCA")
    {gfa=vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)))
    if(type=="R2")commonM[i, 2]=gfa$r.squared
    if(type=="adjR2")commonM[i, 2]=gfa$adj.r.squared
    }
    if(method=="RDA"||method=="dbRDA")
    {gfa=Canonical.Rsq(dv,tmp.design.ct,method=method)
    if(type=="R2")commonM[i, 2]=gfa$unadj
    if(type=="adjR2")commonM[i, 2]=gfa$adj
    }
  }

  commonlist <- vector("list", totalN)

  seqID <- vector()
  for (i in 1:nvar) {
    seqID[i] = 2^(i-1)
  }


  for (i in 1:totalN) {
    bit <- binarymx[1, i]
    if (bit == 1)
      ivname <- c(0, -seqID[1])
    else ivname <- seqID[1]
    for (j in 2:nvar) {
      bit <- binarymx[j, i]
      if (bit == 1) {
        alist <- ivname
        blist <- genList(ivname, -seqID[j])
        ivname <- c(alist, blist)
      }
      else ivname <- genList(ivname, seqID[j])
    }
    ivname <- ivname * -1
    commonlist[[i]] <- ivname
  }

  for (i in 1:totalN) {
    r2list <- unlist(commonlist[i])
    numlist  <-  length(r2list)
    ccsum  <-  0
    for (j in 1:numlist) {
      indexs  <-  r2list[[j]]
      indexu  <-  abs(indexs)
      if (indexu != 0) {
        ccvalue  <-  commonM[indexu, 2]
        if (indexs < 0)
          ccvalue  <-  ccvalue * -1
        ccsum  <-  ccsum + ccvalue
      }
    }
    commonM[i, 3]  <-  ccsum
  }

  orderList <- vector("list", totalN)
  index  <-  0
  for (i in 1:nvar) {
    for (j in 1:totalN) {
      nbits  <-  sum(binarymx[, j])
      if (nbits == i) {
        index  <-  index + 1
        commonM[index, 1] <- j
      }
    }
  }

  outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
  totalRSquare <- sum(commonM[, 3])
  for (i in 1:totalN) {
    outputcommonM[i, 1] <- round(commonM[commonM[i,
                                                 1], 3], digits = 4)
    outputcommonM[i, 2] <- round((commonM[commonM[i,
                                                  1], 3]/totalRSquare) * 100, digits = 2)
  }
  outputcommonM[totalN + 1, 1] <- round(totalRSquare,
                                        digits = 4)
  outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
  rowNames  <-  NULL
  for (i in 1:totalN) {
    ii  <-  commonM[i, 1]
    nbits  <-  sum(binarymx[, ii])
    cbits  <-  0
    if (nbits == 1)
      rowName  <-  "Unique to "
    else rowName  <-  "Common to "
    for (j in 1:nvar) {
      if (binarymx[j, ii] == 1) {
        if (nbits == 1)
          rowName  <-  paste(rowName, iv.name[j], sep = "")
        else {
          cbits  <-  cbits + 1
          if (cbits == nbits) {
            rowName  <-  paste(rowName, "and ", sep = "")
            rowName  <-  paste(rowName, iv.name[j], sep = "")
          }
          else {
            rowName  <-  paste(rowName, iv.name[j], sep = "")
            rowName  <-  paste(rowName, ", ", sep = "")
          }
        }
      }
    }
    rowNames  <-  c(rowNames, rowName)
  }
  rowNames  <-  c(rowNames, "Total")
  rowNames <- format.default(rowNames, justify = "left")
  colNames <- format.default(c("Fractions", " % Total"),
                             justify = "right")
  dimnames(outputcommonM) <- list(rowNames, colNames)

  VariableImportance <- matrix(nrow = nvar, ncol = 2)
  for (i in 1:nvar) {
    VariableImportance[i, 1] <-  round(sum(binarymx[i, ] * (commonM[,
                                                                    3]/apply(binarymx,2,sum))), digits = 4)
  }
  total=round(sum(VariableImportance[,1]),4)
  VariableImportance[, 2] <- round(100*VariableImportance[, 1]/total,2)

  dimnames(VariableImportance) <- list(iv.name, c("Independent","I.perc(%)"))




if(trace)
{outputList <- list(Method_Type=c(method,type),R.squared=total,Commonality = outputcommonM, Var.part = VariableImportance)}
else
{outputList<-list(Method_Type=c(method,type),R.squared=total,Var.part= VariableImportance)}
class(outputList) <- "rdaccahp" # Class definition

return(outputList)
}
}
