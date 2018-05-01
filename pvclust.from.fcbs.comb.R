source('hclust.from.fcbs.R')

# example usage
#system('ls */nodecounts > countfiles')
#pvc = pvclust_from_cpp('1/hc',readLines('countfiles'))
#pdf('pvc.pdf')
#plot(pvc)
#pvrect(pvc)
#pvrect(pvc,lty=2,alpha=0.9)
#dev.off()
#save(pvc,file='pvc.RData')

pvclust_from_cpp = function(hcfile,countfiles=NULL){
	require(pvclust)
	hc = hc_from_fcbs(hcfile)
#	if(is.null(countfiles)) countfiles = dir('.','curr',recursive=T)
	if(is.null(countfiles)) countfiles = dir('.','counts.combined',recursive=T)
	cat(countfiles, '\n')
#	countfiles = readLines(countfiles)
	counts = do.call(cbind, lapply(countfiles,function(x){as.integer(readLines(x))}))

	# fcbs writes of the number of interations it actually finished (reading this is especially useful for incomplete numbers of iterations)
	rs = gsub('/.*','',countfiles)
	cat(rs, '\n')
#	iters = sapply(rs, function(x){as.integer(readLines(sprintf('%s/lastiter',x))[1])})
	iters = sapply(rs, function(x){as.integer(readLines(sprintf('%s/iters.combined',x))[1])})
	cat(iters, '\n')

	rs = as.numeric(rs)
	mboot = list()
	for(i in 1:ncol(counts)){
		mboot[[i]] = list(edges.cnt=counts[,i], nboot=iters[i], r=rs[i])
	}

	for(mb in mboot){
		cat(mb$r,mb$iters,'\n')
	}

	# requires modified version of pvclust to expose merge function (and in the case of big trees to avoid some unnecessarily slow bottlenecks in the code)
	pvclust.merge(data=NULL, object.hclust=hc, mboot=mboot)
}
