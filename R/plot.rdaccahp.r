#' Plot for a \code{\link{rdacca.hp}} object
#'
#' @param object A \code{\link{rdacca.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables or groups towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@ibcas.ac.cn}
#' @export
#' @examples
#' library(vegan)
#' data(mite)
#' data(mite.env)
#' mite.hel <- decostand(mite, "hellinger")
#' avc<-rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#' plot(avc)
#' plot(avc,plot.perc=TRUE)


plot.rdaccahp <- function(object, plot.perc=FALSE,...) 
{
  if(class(object)!="rdaccahp")
   {
    stop("object should be the output of rdacca.hp()")
  }
  
if(plot.perc)
{tips3=data.frame(variable=rownames(object$Hier.part), value=as.numeric(object$Hier.part[,"I.perc(%)"]))
gg=ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="% Individual effect to Rsquare (%I)")+ ggplot2::theme(axis.text = element_text(size = 10))+ ggplot2::theme(axis.title = element_text(size = 13))}
else
{tips2=data.frame(variable=rownames(object$Hier.part), value=as.numeric(object$Hier.part[,"Individual"]))
gg=ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable, -value), y = value)) + ggplot2::geom_bar(stat = "identity")+
  ggplot2::theme_minimal()+ggplot2::labs(x="Variables",y="Individual effect")+ ggplot2::theme(axis.text = element_text(size = 10))+ ggplot2::theme(axis.title = element_text(size = 13))
}

print(gg)

}


