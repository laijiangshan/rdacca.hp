#' plot output of rdacca.hp() based on ggplot2
#' @param  outputList  the output of rdacca.hp().
#' @param  plot.perc logical value, if TRUE, the bar plot (based on ggplot2) of the percentage ot independent effects of variables to total Rsquared, the default is FALSE to show plot with original independent effects.

#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @author {Pedro Peres-Neto} \email{pedro.peres-neto@concordia.ca}

#' @export
#'@examples
#'library(vegan)
#'data(mite)
#'data(mite.env)
#'#Hellinger-transform the species dataset for RDA to deal with the "double zero" problem
#'mite.hel <- decostand(mite, "hellinger")
#'rda.hp <- rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#'plot(rda.hp)
#'plot(rda.hp,plot.perc=TRUE)

plot.rdaccahp<-function(outputList,plot.perc=FALSE)
{
  if(class(x)!="rdaccahp")
   {
    stop("x should be the output of rdacca.hp()")
  }
  
if(plot.perc)
{tips3=data.frame(variable=rownames(outputList$Var.part), value=as.numeric(outputList$Var.part[,"I.perc(%)"]))
gg=ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="Independent effect to Rsquared (%)")}
else
{tips2=data.frame(variable=rownames(outputList$Var.part), value=as.numeric(outputList$Var.part[,"Independent"]))
gg=ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="Independent effect")}
print(gg)
}
