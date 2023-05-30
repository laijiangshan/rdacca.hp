#' Plot for a \code{\link{rdacca.hp}} object
#'
#' @param x A \code{\link{rdacca.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables or groups towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param color Color of variables.
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}
#' @author {Yao Liu} \email{liuyao@stu.xmu.edu.cn}

#' @export
#' @examples
#' library(vegan)
#' data(mite)
#' data(mite.env)
#' mite.hel <- decostand(mite, "hellinger")
#' avc<-rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
#' plot(avc)
#' plot(avc, plot.perc=TRUE)
#' avc<-rdacca.hp(mite.hel,mite.env[,1:3],method="RDA",type="adjR2",var.part=TRUE)
#' plot(avc,color = c("#8DD3C7", "#FFFFB3", "#BEBADA"))


plot.rdaccahp <- function(x, plot.perc = FALSE, color = NULL, ...){
  if (!inherits(x,"rdaccahp")){
    stop("x should be the output of rdacca.hp()")
  }
if(names(x)[3]=="Hier.part")
 { if (plot.perc){
    tips3 = data.frame(variable = rownames(x$Hier.part),
                       value = as.numeric(x$Hier.part[, "I.perc(%)"]))
    gg = ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "% Individual effect to Rsquare (%I)") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))
  } else {
    tips2 = data.frame(variable = rownames(x$Hier.part),
                       value = as.numeric(x$Hier.part[, "Individual"]))
    gg = ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "Individual effect") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))
  }

  return(gg)
}
if(names(x)[3]=="Var.part")
{ 
  Var.part <- as.data.frame(x[[3]])
  #Var.part <- Var.part[which(Var.part$Fractions >=cutoff), ]
  Var.part$Fractions <- round(Var.part$Fractions,3)
  nvar <- nrow(x$Hier.part)
  Constrained <- Var.part$Fractions[2^nvar]
  if (!nvar%in% 2:4)
    stop("Venn diagram supports only 2-4 variables")
  else if (nvar == 2)
	Var <- Var.part$Fractions[c(1, 3, 2)]
  else if (nvar == 3)
    Var <- Var.part$Fractions[c(1:4, 6, 5, 7)]
  else if (nvar == 4)
    Var <- Var.part$Fractions[c(1:5, 7, 6, 8:10, 12, 11, 14, 13, 15)]
    vegan::showvarparts(part = nvar, bg = color, Xnames = rownames(x$Hier.part), labels = as.character(c(Var, 1-Constrained)))
}
}
