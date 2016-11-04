get_hpk = function(hc){
	if(!exists('hpk')){
		library(dendextendRcpp,quietly=T)
		# (to do: can probably re-implement the height-based tree cutting in c++ for bootstrapping/empirical distributions,
		# and/or adapt code from source of dendextendRcpp's Rcpp_cut_lower function)
		cat('Heights per k...\n')
		hpk = dendextendRcpp_heights_per_k.dendrogram(as.dendrogram(hc))
		save(hpk,file='hpk.RData')
	}
	hpk
}

cut_tree = function(hc,hs){
	library(dendextendRcpp,quietly=T)
	clusters = list()
	for(h in hs){
		if(is.na(h))next
		cat('Cutting at height',h,':')
		cls = Rcpp_cut_lower(as.dendrogram(hc), h=h)
		l = length(cls)
		cat(l,'clusters\n')
		clusters[[as.character(l)]] = cls
	}
#	save(clusters,file='clusters.RData')
	clusters
}

plot_exp_and_dendro = function(mat,dendro,expmain='Expression cluster',dendmain='Dendrogram'){
	ng = 1
	if(is.matrix(mat)) ng = nrow(mat)
	layout(matrix( c(1,2), nrow=2, ncol=1 ), heights=c(6,max(3,ng/3)), T)
	matplot(t(mat), type='l', lwd=2, lty=1, col=rgb(0,0,0,0.5), main=expmain, ylab='log2 expression changes')
	mtext(sprintf('%i genes',ng),side=3,line=-2)
	par(mar=c(4,2,2,10))
	plot(dendro, horiz=T,main=dendmain,edge.root=T)
}

plot_clusters_for_ids = function(ratios,clusters_,qids,maxplot=200){
	for(qid in qids){
		cat('plotting clusters for',qid,'\n')
		hasit = sapply(clusters_, function(cls){ which( sapply(cls, function(x){qid %in% labels(x)})) })
		names(hasit) = names(clusters_)
		if(length(hasit)==0) next

		# plot expression and subdendrogram
		for(i in 1:length(hasit)){
			fname = sprintf('cluster.%05i.%s.pdf',as.integer(names(hasit)[i]),qid)
			cat(fname,'\n')
#			if(file.exists(fname)){ cat('plot',fname,'exists, skipping\n'); next }
			if(length(hasit[[i]])==0) next
			if(is.na(hasit[[i]] | is.null(hasit[[i]]))) next
			cl = clusters_[[i]][[hasit[[i]]]]
			ids = as.character(labels(cl))
			ng = length(ids)
			if(ng>maxplot) next
#			labels(cl) = paste( labels(cl), desc[ labels(cl) ] )

			baseh = 6
			pdf(fname,useDingbats=F,height=baseh+max(baseh/2,ng/3),width=baseh*1.5)
			main=sprintf('cluster %i; k=%s',hasit[[i]],names(clusters_)[i])
			plot_exp_and_dendro(rr[ids,], cl, expmain=fname, dendmain=main)
			dev.off()

		}
	}
}

plot_sizedists = function(clusterids){
	for(k in names(clusterids)){
		clids = clusterids[[k]]
		if(length(clids)<=1)next
		fname=sprintf('sizedist.%s.pdf',k)
#		if(file.exists(fname)){ cat('plot',fname,'exists, skipping\n'); next }
		pdf(fname,useDingbats=F)
		plot(density(sapply(clids,length)),main=sprintf('cluster size distribution: k=%s',k))
		dev.off()
	}
}

plot_dendrograms = function(clusters_,dirn='dendro',maxdend=400,descs=NULL){
	cat('dendrograms\n')
	dir.create(dirn,recursive=T)
	for(i in 1:length(clusters_)){
		cat('.')
		fname = sprintf('%s/dendro.%i.png',dirn,i)
		cluster = clusters_[[i]]
		ids = labels(cluster)
		sz = length(ids)
		if(sz > maxdend) next
		if(!is.null(descs)){
			cluster = dendrapply(cluster, function(x){
				if(is.leaf(x)){
					id = attr(x,'label')
					if(id %in% names(descs)) attr(x,'label') = paste(id, descs[id] )
				}
				x
			})
		}

		baseh = 480
		png(fname,height=max(baseh*0.5,sz*15),width=baseh)
		par(mar=c(2,1,2,20))
		plot(cluster, horiz=T, main=sprintf('cluster %i',i))
		dev.off()
	}
	cat('\n')
}
