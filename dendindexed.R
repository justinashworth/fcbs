# borrowed from R source
.memberDend <- function(x) {
	r <- attr(x,"x.member")
	if(is.null(r)) {
		r <- attr(x,"members")
		if(is.null(r)) r <- 1L
	}
	r
}

# node-indexed dendrogram from hclust!
# essential for labeling/picking subdendrograms based on agglomerative order (e.g. for pvclust)
# this should be made the default in R
as.dendrogram.hclust.indexed <- function (object, hang = -1, ...)
{
	stopifnot(length(object$order) > 0L)
	if (is.null(object$labels)) object$labels <- seq_along(object$order)
	z <- list()
	nMerge <- length(oHgt <- object$height)
	if (nMerge != nrow(object$merge)) stop("'merge' and 'height' do not fit!")
	hMax <- oHgt[nMerge]

	for (k in 1L:nMerge) {
		x <- object$merge[k, ]
		if (any(neg <- x < 0)) h0 <- if (hang < 0) 0 else max(0, oHgt[k] - hang * hMax)

		# two leaves
		if (all(neg)) {
			zk <- as.list(-x)
			attr(zk, "members") <- 2L
			attr(zk, "midpoint") <- 0.5
			objlabels <- object$labels[-x]
			attr(zk[[1L]], "label") <- objlabels[1L]
			attr(zk[[2L]], "label") <- objlabels[2L]
			attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- 1L
			attr(zk[[1L]], "height")  <- attr(zk[[2L]], "height") <- h0
			attr(zk[[1L]], "leaf")	<- attr(zk[[2L]], "leaf") <- TRUE
			attr(zk, "index") = k

		# one leaf, one node
		} else if (any(neg)) {
			X <- as.character(x)
			isL <- x[1L] < 0 # is leaf left?
			zk <- if(isL) list(-x[1L], z[[X[2L]]]) else	list(z[[X[1L]]], -x[2L])
			attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + 1L
			attr(zk, "midpoint") <- (.memberDend(zk[[1L]]) + attr(z[[X[1 + isL]]], "midpoint"))/2
			attr(zk, "index") = k
			attr(zk[[2 - isL]], "members") <- 1L
			attr(zk[[2 - isL]], "height") <- h0
			attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - isL]]
			attr(zk[[2 - isL]], "leaf") <- TRUE
			z[[X[1 + isL]]] <- NULL

		# two non-leaf nodes
		} else {
			x <- as.character(x)
			# "merge" the two ('earlier') branches:
			zk <- list(z[[x[1L]]], z[[x[2L]]])
			attr(zk, "members") <- attr(z[[x[1L]]], "members") + attr(z[[x[2L]]], "members")
			attr(zk, "midpoint") <- (attr(z[[x[1L]]], "members") + attr(z[[x[1L]]], "midpoint") + attr(z[[x[2L]]], "midpoint"))/2
			attr(zk, "index") = k
			z[[x[1L]]] <- z[[x[2L]]] <- NULL
		}
		attr(zk, "height") <- oHgt[k]
		z[[as.character(k)]] <- zk
	}

	z <- z[[as.character(k)]]
	class(z) <- "dendrogram"
	z
}

