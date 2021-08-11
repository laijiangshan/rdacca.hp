#' Hierarchical and Variation Partitioning for Canonical Analysis Without Limiting the Number of Predictors (Matrices)

#' @param  dv  Response variable, either a numeric vector or a matrix. If method="dbRDA", dv should be a "dist" matrix.
#' @param  iv Predictors (explanatory variable), either a data frame or a list of data frames. If it is a data frame, the relative importance of each column of the data frame will be evaluated; if it is a list, the relative importance of each element (matrix) will be evaluated.
#' @param  method Type of canonical analysis used for variation partitioning, should be a character string, either "RDA", "dbRDA" or "CCA", the default is "RDA". If the response variable (dv) is a numerical vector and method="RDA", the hierarchical and variation partitioning for the classical multiple regression is implemented.
#' @param  type The type of total explained variation, either "R2" or "adjR2", in which "R2" is unadjusted R-square and "adjR2" is adjusted R-square, the default is "adjR2". The adjusted R-square is calculated using Ezekiel's formula (Ezekiel 1930) for RDA and dbRDA, while permutation procedure is used for CCA (Peres-Neto et al. 2006). 
#' @param  scale Logical; If the columns of dv should be standardized to unit variance when method="RDA" is applied.
#' @param  add Logical; If a constant should be added to the non-diagonal values to euclidify dissimilarities (see dbrda function in vegan for details). Choice "lingoes" (or TRUE) uses the recommended method of Legendre & Anderson (1999: "method 1") and "cailliez" uses their "method 2". The argument has an effect only when method="dbRDA".
#' @param  sqrt.dist Logical, If the square root of dissimilarities should be taken. This often euclidifies dissimilarities. The argument has an effect only when method="dbRDA"(see dbrda function in vegan for details).
#' @param  n.perm An integer; Number of permutations for computing the adjusted R-square for CCA.
#' @param  var.part Logical; If TRUE, the result of variation partitioning (2^N-1 fractions for N predictors or matrices) are shown, the default is FALSE.

#' @details This function conducts variation partitioning and hierarchical partitioning to calculate the unique, average shared (referred as to "common") and individual contributions of each predictor (or matrix) towards explained variation (R-square) on canonical analysis (RDA,CCA and dbRDA).
#' Variation partitioning should be conducted before hierarchical partitioning. The former emphasizes unique and common variation among predictors, the latter emphasizes the overall importance of each predictor (or group of predictors). This function simultaneously implements variation and hierarchical partitioning for single- and multiple-response models without limiting in the number of predictors / matrices of predictors. 
#' This package was inspired by Chevan & Sutherland (1991)’s paper, the "yhat" (Nimon, Oswald & Roberts 2013) and "hier.part"(Walsh & Mac Nally 2013) R packages. 

#' @return a list containing
#' @return \item{Method_Type}{Type of canonical analysis and the type of total explained variation.}
#' @return \item{Total_explained_variation}{The explained variation for the global model.}
#' @return \item{Var.part}{If var.part=TRUE, a matrix containing the value and percentage of all commonality (2^N-1 for N predictors or matrices).}
#' @return \item{Hier.part}{A matrix containing unique, average shared, individual effects and percentage of individual effect towards total explained variation for each predictor or matrix.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Pedro Peres-Neto} \email{pedro.peres-neto@concordia.ca}

#' @references
#' \itemize{
#' \item Lai J.,Zou Y., Zhang J.,Peres-Neto P.(2021) rdacca.hp: an R package for generalizing hierarchical and variation partitioning in multiple regression and canonical analysis.bioRxiv 2021.03.09.434308,<DOI:10.1101/2021.03.09.434308>
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning. American Statistician, 45, 90-96. doi:10.1080/00031305.1991.10475776
#' \item Nimon, K., Oswald, F.L. & Roberts, J.K. (2013). Yhat: Interpreting regression effects. R package version 2.0.0.
#' \item Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. Ecological Monographs, 69, 1–24.<DOI:10.1890/0012-9615(1999)069[0001:dbratm]2.0.co;2>
#' \item Walsh, C.J. & Mac Nally, R. (2013) hier.part: Hierarchical Partitioning. R package version 1.0-4.
#' \item Peres-Neto, P.R., Legendre, P., Dray, S. & Borcard, D. (2006) Variation partitioning of species data matrices: Estimation and comparison of fractions. Ecology, 87, 2614-2625.<DOI: doi.org/10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2>
#' \item Ezekiel, M. (1930) Methods of Correlational Analysis. Wiley, New York
#' }


#'@export
#'@examples
#'library(vegan)
#'data(mite)
#'data(mite.env)
#'data(mite.xy)
#'data(mite.pcnm)
#'#Hellinger-transform the species dataset for RDA
#'mite.hel <- decostand(mite, "hellinger")
#'rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#'rdacca.hp(vegdist(mite),mite.env,method="dbRDA",type="adjR2")
#'rdacca.hp(mite,mite.env,method="CCA",type="adjR2")
#'iv <- list(env=mite.env,xy=mite.xy,pcnm=mite.pcnm)
#'rdacca.hp(mite.hel,iv,method="RDA",var.part = TRUE)
#'rdacca.hp(vegdist(mite),iv,method="dbRDA",var.part = TRUE)
#'rdacca.hp(mite,iv,method="CCA",var.part = TRUE)


rdacca.hp <- function (dv,iv,method=c("RDA","dbRDA","CCA"),type=c("adjR2","R2"),scale=FALSE,add = FALSE, sqrt.dist = FALSE,n.perm=1000,var.part = FALSE) 
{

if(is.data.frame(iv))
{
  if(sum(is.na(dv))>=1|sum(is.na(iv))>=1)
  {stop("NA/NaN/Inf is not allowed in this analysis")}
  if(nrow(iv)<=ncol(iv))
  {stop("sample size (row) is less than the number of predictors")}

  else
  {method <- method[1]
  type <- type[1]
  if(inherits(dv, "dist"))
  {method <- "dbRDA"}
  if(method=="dbRDA"||method=="dbrda"||method=="DBRDA")
  {
    if(!inherits(dv, "dist"))
      return("response variables should be a 'dist' matrix for dbRDA")
  }



  if(method=="RDA"||method=="rda")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  {dv<-scale(dv,scale=scale)}

  iv <- data.frame(iv)
  ivname <- colnames(iv)
  iv.name <- ivname
  nvar <- dim(iv)[2]
  if (nvar < 2)
    stop("Analysis not conducted. Insufficient number of predictors.")

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
    
    if(method=="RDA"||method=="rda")
    {gfa <- simpleRDA2(dv,tmp.design.ct)
    if(type=="R2")commonM[i, 2] <- gfa$Rsquare
    if(type=="adjR2")commonM[i, 2] <- gfa$RsquareAdj
    }
	if(method=="CCA"||method=="cca")
    {gfa <- vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)),permutations = n.perm)
    if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
    }
	if(method=="dbRDA"||method=="dbrda"||method=="DBRDA")
    {
	gfa <- vegan::RsquareAdj(vegan::dbrda(dv~.,data=data.frame(tmp.design.ct),add=add,sqrt.dist = sqrt.dist))
	if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
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

  VariableImportance <- matrix(nrow = nvar, ncol = 4)
  for (i in 1:nvar) {
	VariableImportance[i, 3] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
  }
  
  VariableImportance[,1] <- outputcommonM[1:nvar,1]
  VariableImportance[,2] <- VariableImportance[,3]-VariableImportance[,1]
  
  total=round(sum(VariableImportance[,3]),digits = 4)
  VariableImportance[, 4] <- round(100*VariableImportance[, 3]/total,2)

  dimnames(VariableImportance) <- list(iv.name, c("Unique","Average.share","Individual","I.perc(%)"))


if(var.part)
{outputList <- list(Method_Type=c(method,type),Total_explained_variation=total,Var.part = outputcommonM, Hier.part = VariableImportance)}
else
{outputList<-list(Method_Type=c(method,type),Total_explained_variation=total,Hier.part= VariableImportance)}

class(outputList) <- "rdaccahp" # Class definition
outputList

}
}
#
else
{nvar  <-  length(iv)
  if(sum(unlist(lapply(iv,is.data.frame)))<nvar)
  stop("data.frame is required for each group explanatory table")

  if(sum(is.na(dv))>=1|sum(is.na(unlist(iv)))>=1)
  {stop("NA/NaN/Inf is not allowed in this analysis")}
  
  else
  {method <- method[1]
  type <- type[1]
  if(inherits(dv, "dist"))
  {method <- "dbRDA"}
  if(method=="dbRDA"||method=="dbrda"||method=="DBRDA"){
    if(!inherits(dv, "dist"))
      return("dv should be a 'dist' matrix for dbRDA")
  }

  if(method=="RDA"||method=="rda")
  {dv<-scale(dv,scale=scale)}
	
	
	#method <- "RDA"
	#if(CCA){method="CCA"}
	#if(inherits(dv, "dist")){method="dbRDA"}
		
    ilist <- names(iv)
	if(is.null(ilist))
	{names(iv) <- paste("X",1:nvar,sep="")}
	else
	{whichnoname <- which(ilist=="")
    names(iv)[whichnoname] <- paste("X",whichnoname,sep="")}
	
	ilist <- names(iv)
	
	ivlist <- ilist
	iv.name <- ilist
	
    if (nvar < 2) 
    stop("Analysis not conducted. Insufficient number of predictor groups.")
  
   ivID <- matrix(nrow = nvar, ncol = 1)
    for (i in 0:nvar - 1) {
        ivID[i + 1]  <-  2^i
    }

   totalN  <-  2^nvar - 1
    #binarymx <- matrix(0, nvar, totalN)
    #for (i in 1:totalN) {
    #    binarymx <- creatbin(i, binarymx)
    # }

  commonM <- matrix(nrow = totalN, ncol = 3)
  binarymx <- matrix(0, nvar, totalN)
  for (i in 1:totalN) {
    binarymx <- creatbin(i, binarymx)
  }

  commonM <- matrix(nrow = totalN, ncol = 3)
	
    for (i in 1:totalN) {
        ivls <- iv[as.logical(binarymx[, i])]
		N <- length(ivls)
		if(N==1)
		{
		tmp <- ivls[[1]]
		tmp.design <- as.matrix(model.matrix(~ ., tmp))
        tmp.design.ct <- as.matrix(data.frame(scale(tmp.design,scale=FALSE))[-1])
	if(method=="RDA"||method=="rda")
    {gfa <- simpleRDA2(dv,tmp.design.ct)
    if(type=="R2")commonM[i, 2] <- gfa$Rsquare
    if(type=="adjR2")commonM[i, 2] <- gfa$RsquareAdj
    }
	if(method=="CCA"||method=="cca")
    {gfa <- vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)),permutations = n.perm)
    if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
    }

	if(method=="dbRDA"||method=="dbrda"||method=="DBRDA")
    {
	gfa <- vegan::RsquareAdj(vegan::dbrda(dv~.,data=data.frame(tmp.design.ct),add=add,sqrt.dist = sqrt.dist))
	if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
    }	
  }
  
		#commonM[i, 2] <-Canonical.Rsq(as.matrix(dv),tmp.design.ct,method=method)$adj}
		if(N>1)
		{inv <- ivls[[1]]
   for(k in 2:N)
  {inv <- cbind(inv,ivls[[k]])}
		tmp.design <- as.matrix(model.matrix(~ ., inv))
        tmp.design.ct <- as.matrix(data.frame(scale(tmp.design,scale=FALSE))[-1])
	if(method=="RDA"||method=="rda")
    {gfa <- simpleRDA2(dv,tmp.design.ct)
    if(type=="R2")commonM[i, 2] <- gfa$Rsquare
    if(type=="adjR2")commonM[i, 2] <- gfa$RsquareAdj
    }	

	if(method=="CCA"||method=="cca")
    {gfa <- vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)),permutations = n.perm)
    if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
    }
   if(method=="dbRDA"||method=="dbrda"||method=="DBRDA")
    {
	gfa <- vegan::RsquareAdj(vegan::dbrda(dv~.,data=data.frame(tmp.design.ct),add=add,sqrt.dist = sqrt.dist))
	if(type=="R2")commonM[i, 2] <- gfa$r.squared
    if(type=="adjR2")commonM[i, 2] <- gfa$adj.r.squared
    }		
	#commonM[i, 2] <-Canonical.Rsq(as.matrix(dv),tmp.design.ct,method=method)$adj}		
	}
    }
    commonalityList <- vector("list", totalN)
	
    for (i in 1:totalN) {
        bit  <-  binarymx[1, i]
        if (bit == 1) 
            ilist <- c(0, -ivID[1])
        else ilist <- ivID[1]
        for (j in 2:nvar) {
            bit  <-  binarymx[j, i]
            if (bit == 1) {
                alist <- ilist
                blist <- genList(ilist, -ivID[j])
                ilist <- c(alist, blist)
            }
            else ilist <- genList(ilist, ivID[j])
        }
        ilist <- ilist * -1
        commonalityList[[i]] <- ilist
    }

    for (i in 1:totalN) {
        r2list <- unlist(commonalityList[i])
        numlist  <- length(r2list)
        ccsum = 0
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
    rowNames = NULL
    for (i in 1:totalN) {
        ii  <-  commonM[i, 1]
        nbits  <-  sum(binarymx[, ii])
        cbits  <-  0
        if (nbits == 1) 
            rowName  <-  "Unique to "
        else rowName = "Common to "
        for (j in 1:nvar) {
            if (binarymx[j, ii] == 1) {
                if (nbits == 1) 
                  rowName  <-  paste(rowName, ivlist[j], sep = "")
                else {
                  cbits = cbits + 1
                  if (cbits == nbits) {
                    rowName  <-  paste(rowName, "and ", sep = "")
                    rowName  <-  paste(rowName, ivlist[j], sep = "")
                  }
                  else {
                    rowName  <-  paste(rowName, ivlist[j], sep = "")
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
  
   # outputCCbyVar <- matrix(nrow = nvar, ncol = 3)
    #for (i in 1:nvar) {
   #     outputCCbyVar[i, 1]  <-  outputcommonM[i, 1]
    #    outputCCbyVar[i, 3]  <-  round(sum(binarymx[i, ] * (commonM[, 
    #        3])), digits = 4)
    #    outputCCbyVar[i, 2]  <- outputCCbyVar[i, 3] - outputCCbyVar[i, 1]
   # }
    #dimnames(outputCCbyVar) <- list(ivlist, c("Unique", "Common", "Total"))
	
   # outputList <- list(Method=paste("Partition of variance in ",method,sep=""),Partition = outputcommonM, CCTotalbyVar = outputCCbyVar)

  VariableImportance <- matrix(nrow = nvar, ncol = 4)
  for (i in 1:nvar) {
	VariableImportance[i, 3] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
  }
  
  VariableImportance[,1] <- outputcommonM[1:nvar,1]
  VariableImportance[,2] <- VariableImportance[,3]-VariableImportance[,1]
  
  total=round(sum(VariableImportance[,3]),digits = 4)
  VariableImportance[, 4] <- round(100*VariableImportance[, 3]/total,2)

  dimnames(VariableImportance) <- list(iv.name, c("Unique","Average.share","Individual","I.perc(%)"))


if(var.part)
{outputList <- list(Method_Type=c(method,type),Total_explained_variation=total,Var.part = outputcommonM, Hier.part = VariableImportance)}
else
{outputList<-list(Method_Type=c(method,type),Total_explained_variation=total,Hier.part= VariableImportance)}

class(outputList) <- "rdaccahp" # Class definition
outputList
}
}		
}

