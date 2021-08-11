#'Internal function for rdacca.hp() to create diagonal matrix
#' @param  col Imput number.
#' @param  binmatrix Imput empty matrix.
#' @keywords internal
creatbin <-
function(col, binmatrix) {
row<-1
val<-col
	while (val!=0){              
	if (odd(val)) {
		binmatrix[row,col]=1 
	}
	val<-as.integer(val/2)
	row<-row+1
}
##Return matrix
return(binmatrix)
}


#'Internal function for rdacca.hp() to determine whether the odd number
#' @param  val Imput number.
#' @keywords internal
odd <- function (val) 
{
    if (val%%2 == 0) 
     return(FALSE)
	 else
    return(TRUE)
}