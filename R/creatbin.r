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


odd <- function (val) 
{
    if (val%%2 == 0) 
     return(FALSE)
	 else
    return(TRUE)
}