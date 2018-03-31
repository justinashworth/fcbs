# simple function to load up an hclust object that was written to file by fcbs (not bothering with a real R package interface for now)
hc_from_fcbs = function(hcfile='hc',dist.method='',hclust.method='fastcluster'){
	# convert output files from c++ fastcluster mod into R hclust object
	hc = list()
	class(hc) = 'hclust'
	hc$merge =  as.integer( unlist( strsplit(readLines( sprintf('%s.merge',hcfile)),' ')))
	hc$height = as.numeric( unlist( strsplit(readLines( sprintf('%s.height',hcfile)),' ')))
	hc$order =  as.integer( unlist( strsplit(readLines( sprintf('%s.order',hcfile)),' ')))
	hc$labels =             unlist( strsplit(readLines( sprintf('%s.labels',hcfile)),' '))

	dim(hc$merge) = c(length(hc$merge)/2, 2)
	if(nrow(hc$merge) != length(hc$height)){
		cat('ERROR: merges/heights length mismatch!!\n')
		return(NULL)
	}
	hc$labels = gsub('\"','',hc$labels)
	hc$method = hclust.method
	hc$call = match.call()
	hc$dist.method = dist.method
	save(hc,file='hc.RData')
	hc
}
