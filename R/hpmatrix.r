hpmatrix=function (n) 
{   require(gtools)
    x <- cbind(1:n, matrix(0, n, n - 1))
    for (i in 2:n) {
        nc <- dim(combinations(n, i, 1:n))[1]
        x <- rbind(x, cbind(combinations(n, i, 1:n), matrix(0,nc, n - i)))
    }
x
}

