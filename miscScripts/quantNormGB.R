######################
#
# 140530 
# The quantNorm function in the Bolstad 'preprocessCore' R package does not appear to retaing the order of 'probes' of equal value.
# So wrote a custom version that does.
#
######################

quantNormGB = function(x) {
	rows <- dim(x)[1]
    cols <- dim(x)[2]
    if (!is.matrix(x)) {
        stop("Matrix expected in normalize.quantiles")
    }
    if (is.integer(x)) {
        x <- matrix(as.double(x), rows, cols)
        copy <- FALSE
    }
	
	y = sortIdxMatrix = c()
	for (i in 1:cols) {
		colx = x[,i]
		zeroIdx = colx==0
		colx[zeroIdx] = -1 * (rows - which(zeroIdx))
		xord = order(colx)
		colx[zeroIdx] = 0
		coly = colx[xord]
		y = cbind(y, coly)
		sortIdxMatrix = cbind(sortIdxMatrix, xord)
	}
	yQuantMeans = rowMeans(y)
	# print(head(yQuantMeans))
	
	xNorm = c()
	for (i in 1:cols) {
		colxNorm = rep(0, rows)
		colxNorm[sortIdxMatrix[,i]] = yQuantMeans
		xNorm = cbind(xNorm, colxNorm)
	}
	
	colnames(xNorm) = colnames(x)
	rownames(xNorm) = rownames(x)
	return(as.matrix(xNorm))
}	
