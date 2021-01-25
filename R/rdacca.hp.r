#' Hierarchical Partitioning for Canonical Analysis via Commonality Analysis

#' @param  dv Response variables. if method="dbRDA", dv is the "dist" matrix.
#' @param  iv Explanatory variables, both data frame and list are allowed, data frame is for accessing each variable and list is for accessing each explanatory matrix.
#' @param  method The type of canonical analysis: RDA, dbRDA or CCA, the default is "RDA".
#' @param  type The type of total explained variation: "adjR2" is adjusted R-squared and "R2" for unadjusted R-squared, the default is "adjR2".
#' @param  trace logical value, if TRUE, the vaules of commonality (2^N-1 fractions for N explanatory variables or groups) are outputed,the default is FALSE.
#' @param  plot.perc logical value, if TRUE, the bar plot (based on ggplot2) of the percentage to independent effects of variables to total Rsquared, the default is FALSE to show plot with original independent effects.

#' @details This function calculates the independent contribution of each explanatory variable or group to explained variation (R-squared) on canonical analysis (RDA,CCA and dbRDA),
#' applying the hierarchy algorithm of Chevan and Sutherland (1991). The algorithm is that all joint R-squared will be decomposed into equal fractions by number
#' of involved explanatory variables and average assigned to these variables. Independent R-squared of each variable or group will be the sum of assigned R-squared from joint R-squared and unique R-squared.

#' @return a list containing
#' @return \item{Method_Type}{The type of canonical analysis and the type of total explained variation.}
#' @return \item{R.squared}{The explained variation for global model.}
#' @return \item{Commonality}{If trace=TRUE,a mtrix listing tha value and percentage of all commonality (2^N-1 for N explanatory variables or groups).}
#' @return \item{Var.part}{A matrix listing independent effect and its percentage to total explained variation for each explanatory variable or groups.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Pedro Peres-Neto} \email{pedro.peres-neto@concordia.ca}

#' @export
#'@examples
#'library(vegan)
#'data(mite)
#'data(mite.env)
#'data(mite.xy)
#'data(mite.pcnm)
#'#Hellinger-transform the species dataset for RDA to deal with the "double zero" problem
#'mite.hel <- decostand(mite, "hellinger")
#'rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#'rdacca.hp(vegdist(mite),mite.env,method="dbRDA",type="adjR2")
#'rdacca.hp(mite,mite.env,method="CCA",type="adjR2")
#'iv <- list(env=mite.env,xy=mite.xy,pcnm=mite.pcnm)
#'rdacca.hp(mite.hel,iv,method="RDA",trace = TRUE,plot.perc = FALSE)
#'rdacca.hp(vegdist(mite),iv,method="dbRDA",trace = TRUE,plot.perc = FALSE)
#'rdacca.hp(mite,iv,method="CCA",trace = TRUE,plot.perc = FALSE)


rdacca.hp <- function (dv,iv,method=c("RDA","dbRDA","CCA"),type=c("adjR2","R2"),trace = FALSE,plot.perc = FALSE) 
{
if(is.data.frame(iv))
#{rdacca.hp(dv, iv,method,type,trace,plot.perc)} 
#if(is.list(iv))
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
#class(outputList) <- "rdacca.hp" # Class definition

if(plot.perc)
{tips3=data.frame(variable=rownames(outputList$Var.part), value=as.numeric(outputList$Var.part[,"I.perc(%)"]))
gg=ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="Independent effect to Rsquared (%)")}
else
{tips2=data.frame(variable=rownames(outputList$Var.part), value=as.numeric(outputList$Var.part[,"Independent"]))
gg=ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="Independent effect")}

print(gg)
#class(outputList) <- "rdaccahp" # Class definition

return(outputList)
}
}
#
else
{rdacca.mhp(dv, iv,method,type,trace,plot.perc)}
}

