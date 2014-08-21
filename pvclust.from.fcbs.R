source('~/xcode/fcbs/fcbs/hclust.from.fcbs.R')

# example usage
#system('ls */nodecounts > countfiles')
#pvc = pvclust_from_cpp('1/hc','countfiles',1000)
#pdf('pvc.pdf')
#plot(pvc)
#pvrect(pvc)
#pvrect(pvc,lty=2,alpha=0.9)
#dev.off()

pvclust_from_cpp = function(hcfile,countsfiles,nboot){
	require(pvclust)
	hc = hc_from_fcbs(hcfile)
	countsfiles = readLines(countsfiles)
	rs = as.numeric( gsub('/.*','',countsfiles) )
	counts = do.call(cbind, lapply(countsfiles,function(x){as.integer(readLines(x))}))

	mboot = list()
	for(i in 1:ncol(counts)){
		mboot[[i]] = list(edges.cnt=counts[,i], nboot=nboot, r=rs[i])
	}

	# requires modified version of pvclust to expose merge function (and in the case of big trees to avoid some unnecessarily slow bottlenecks in the code)
	pvclust.merge(data=NULL, object.hclust=hc, mboot=mboot)
}

save(pvc,file='pvc.RData')

# plot a subdendrogram (instead of the whole tree) with pvclust bootstrap p-values at the tree nodes
# work in progress
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
