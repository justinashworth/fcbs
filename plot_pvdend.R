library(pvclust)
#source('dendindexed.R')

# plot a subdendrogram (instead of the whole tree) with pvclust bootstrap p-values at the tree nodes
plot_pvdend = function(dend,labord,pvaxes,desc=NULL,...){
	require(pvclust)
	ids = labels(dend)
	if(!is.null(desc)){
		dend = dendrapply(dend, function(x){
#			if(is.leaf(x)) attr(x,'label') = paste(attr(x,'label'),desc[attr(x,'label')])
			if(is.leaf(x)) attr(x,'label') = desc[attr(x,'label')]
			x
		})
	}
	inds = which(labord %in% ids)
	par(mar=c(5,1,3,24))
	plot(dend,horiz=T,edge.root=T,...)

	relev = pvaxes$x.axis >= min(inds) & pvaxes$x.axis <= max(inds)
#	claxes = pvaxes[relev,]
	xcrd = pvaxes$x.axis[relev]
	hcrd = pvaxes$y.axis[relev]
	xadj = 1
#	xadj = 1 + 0.01*length(inds)
	hrange = range(hcrd)
	hadj = 0
#	hadj = 0 + 0.01*(hrange[2]-hrange[1])
	adj = c(1.3,-0.5)
#	adj = NULL
	pvs = ""
	if(exists('pvc')) {
		pvs = pvc$edges$au[relev]
		pvs = sapply(pvs,function(x){if(x==0) 0 else round(100*x)})
	}
	else cat('no pvc, skipping pvs\n')
	font = 2
#	font = sapply(pvdisp,function(x){if(x<pvcut*100) 1 else 2})
	cols = rep('black',length(pvs))
#	cols[pvs<pvcut*100] = 'gray'
	text(hcrd+hadj, xcrd+xadj-min(inds), pvs, font=font, srt=0, cex=0.7, col=cols, adj=adj)
}
