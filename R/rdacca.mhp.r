#' Hierarchical Partitioning for Canonical Analysis by Explanatory Matrices


#' @param  dv Data frame or matrix containing the response data table or dissimilarity structure inheriting from dist. In community ecology, that table is often a site-by-species table or a dissimilarity object.
#' @param  iv A list containing N groups explanatory data frame or matrix.
#' @param  method The type of canonical analysis: RDA, dbRDA or CCA, the default is "RDA".
#' @param  type The type of total explained variation: "adjR2" is adjusted R-squared and "R2" for unadjusted R-squared, the default is "adjR2".
#' @param  trace logical value, if TRUE, the vaules of commonality (2^N-1for N group of explanatory variables) are outputed,the default is FALSE.
#' @param  plot.perc logical value, if TRUE, the bar plot (based on ggplot2) of the percentage to independent effects of variables group to total Rsquared, the default is FALSE to show plot with original independent effects.

#' @details This function calculates the independent contribution of each explanatory group to explained variation (R-squared) on canonical analysis (RDA,CCA and dbRDA),
#' applying the hierarchy algorithm of Chevan and Sutherland (1991). The algorithm is that all joint R-squared will be decomposed into equal fractions by number
#' of involved explanatory groups and average assigned to these group. Independent R-squared of each variable will be the sum of assigned R-squared from joint R-squared and unique R-squared.


#' @return a list containing
#' @return \item{Method_Type}{The type of canonical analysis and the type of total explained variation.}
#' @return \item{R.squared}{The explained variation for global model.}
#' @return \item{Commonality}{If trace=TRUE,a mtrix listing tha value and percentage of all commonality (2^N-1 for N explanatory variables groups).}
#' @return \item{Var.part}{A matrix listing independent effect and its percentage to total explained variation for each explanatory variable group.}

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Pedro Peres-Neto} \email{pedro.peres-neto@concordia.ca}

#'@examples
#'require(vegan)
#'data(mite)
#'data(mite.env)
#'data(mite.xy)
#'data(mite.pcnm)
#'mite.hel <- decostand(mite, "hellinger")
#'iv <- list(env=mite.env,xy=mite.xy,pcnm=mite.pcnm)
#'rdacca.mhp(mite.hel,iv,method="RDA",trace = TRUE,plot.perc = FALSE)
#'rdacca.mhp(vegdist(mite),iv,method="dbRDA",trace = TRUE,plot.perc = FALSE)
#'rdacca.mhp(mite,iv,method="CCA",trace = TRUE,plot.perc = FALSE)

rdacca.mhp <- function (dv, iv,method=c("RDA","dbRDA","CCA"),type=c("adjR2","R2"),trace = FALSE,plot.perc = FALSE) 
  { 
  nvar  <-  length(iv)
  if(sum(unlist(lapply(iv,is.data.frame)))<nvar)
  stop("Error: data.frame is required for each group explanatory table")

  if(sum(is.na(dv))>=1|sum(is.na(unlist(iv)))>=1)
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
    stop("Analysis not conducted. Insufficient number of explanatory variables.")
  
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
	if(method=="CCA")
    {gfa=vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)))
    if(type=="R2")commonM[i, 2]=gfa$r.squared
    if(type=="adjR2")commonM[i, 2]=gfa$adj.r.squared
    }
    if(method=="RDA"||method=="dbRDA")
    {gfa=Canonical.Rsq(dv,tmp.design.ct,method=method)
    if(type=="R2")commonM[i, 2]=gfa$unadj
    if(type=="adjR2")commonM[i, 2]=gfa$adj}
	}		
		#commonM[i, 2] <-Canonical.Rsq(as.matrix(dv),tmp.design.ct,method=method)$adj}
		if(N>1)
		{inv <- ivls[[1]]
   for(k in 2:N)
  {inv <- cbind(inv,ivls[[k]])}
		tmp.design <- as.matrix(model.matrix(~ ., inv))
        tmp.design.ct <- as.matrix(data.frame(scale(tmp.design,scale=FALSE))[-1])
		if(method=="CCA")
    {gfa=vegan::RsquareAdj(vegan::cca(dv~.,data=data.frame(tmp.design.ct)))
    if(type=="R2")commonM[i, 2]=gfa$r.squared
    if(type=="adjR2")commonM[i, 2]=gfa$adj.r.squared
    }
    if(method=="RDA"||method=="dbRDA")
    {gfa=Canonical.Rsq(dv,tmp.design.ct,method=method)
    if(type=="R2")commonM[i, 2]=gfa$unadj
    if(type=="adjR2")commonM[i, 2]=gfa$adj}	
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
