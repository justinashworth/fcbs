library(pvclust)
source('~/code/dendindexed.R')

# plot a subdendrogram (instead of the whole tree) with pvclust bootstrap p-values at the tree nodes
plot_pvdend = function(dend,hclabelorder,pvclustaxes,desc=NULL,...){
	require(pvclust)
	ids = labels(dend)
	if(!is.null(desc)){
		dend = dendrapply(dend, function(x){
			if(is.leaf(x)) attr(x,'label') = paste(attr(x,'label'),desc[attr(x,'label')])
			x
		})
	}
	inds = which(hclabelorder %in% ids)
	par(mar=c(5,1,3,20))
	plot(dend,horiz=T,edge.root=T,...)

	relev = pvclustaxes$x.axis >= min(inds) & pvclustaxes$x.axis <= max(inds)
	claxes = pvclustaxes[relev,]
	indcoord = pvclustaxes$x.axis[relev]
	heightcoord = pvclustaxes$y.axis[relev]
	indadj = 1
#	indadj = 1 + 0.01*length(inds)
	hrange = range(heightcoord)
	heightadj = 0
#	heightadj = 0 + 0.01*(hrange[2]-hrange[1])
	adj = c(1.3,-0.5)
#	adj = NULL
	pvs = pvc$edges$au[relev]
	pvdisp = sapply(pvs,function(x){if(x==0) 0 else round(100*x)})
	font = sapply(pvdisp,function(x){if(x<pvcut*100) 1 else 2})
	cols = rep('black',length(pvs))
#	cols[pvs<pvcut*100] = 'gray'
	text(heightcoord+heightadj, indcoord+indadj-min(inds), pvdisp, font=font, srt=0, cex=0.7, col=cols, adj=adj)
}
