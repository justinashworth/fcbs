source('~/code/hclust.lib.R')
source('~/code/fcbs/pvclust.from.fcbs.R')
source('~/code/fcbs/plot_pvdend.R')
require(pvclust)

if(!exists('dd')){
	rfile = 'data'
	cat('reading data',rfile,'\n')
	dd = read.table(rfile)
}

cache = c('pvc','hpk')
cache = c(cache, 'uncontained.clusters')

for(f in cache){
	fd = sprintf('%s.RData',f)
	if(file.exists(fd)){
		cat('reading cache file',fd,'\n')
		load(fd)
	}
}

if(!exists('pvc')){
	pvc = pvclust_from_cpp('1.0/hc')
	save(pvc,file='pvc.RData')
}

if(!exists('clusters')){
	pt='au'
	pvcut=0.90
	cat(sprintf('pvpick with P[%s]>=%f\n',pt,pvcut))
	clpick = pvpick(pvc,pv=pt,alpha=pvcut,max.only=T)
	save(clpick,file='pvpick.RData')
	clusters = clpick$clusters
}

plots = T
if(plots){

if(!exists('hpk')){
	hpk=get_hpk(pvc$hclust)
	save(hpk,file='hpk.RData')
}

if(!exists('cls')){
library(dendextendRcpp)
ncl=500
if(!ncl %in% names(hpk)){
	ncls = as.numeric(names(hpk))
	ncl = ncls[ ncls>ncl ][1]
}
cat('Rcpp_cut_lower for',ncl,'clusters\n')
cls = dendextendRcpp::Rcpp_cut_lower(as.dendrogram(pvc$hclust),height=hpk[[as.character(ncl)]])
save(cls,file='cls.RData')
}

plot_dends=T
if(plot_dends){
	library(fastcluster)
	baseh = 480
	maxdend = 200
	denddir = 'dendro'
	dir.create(denddir)
	# requires modified version of pvclust to expose functions and work properly with externally produced bootstraps
	axes = hc2axes(pvc$hclust)
	ordlabs = pvc$hclust$labels[ pvc$hclust$order ]

	qsi = sapply(qs, function(x){ which(sapply(clusters, function(y){x %in% y})) } )
	for(qq in qs){
		for(i in qsi[[qq]]){
			cat('plots for',qq,'\n')
			cl = which(sapply(cls,function(x){qq %in% labels(x)}))
			ids = labels(cls[[cl]])
			sz = length(ids)
			if(sz>maxdend){
				cat('cluster',i,'too big to plot dendrogram (',sz,'genes)\n')
				next
			}
#			png(paste(denddir,'/',qq,'.',i,'.dend.png',sep=''),height=max(baseh/2,length(ids)*15),width=baseh)
			pdf(paste(denddir,'/',qq,'.',i,'.dend.pdf',sep=''),height=max(4,length(ids)/4))
			plot_pvdend(cls[[cl]], ordlabs, axes, desc)
			dev.off()

			pdf(paste(qq,'.',i,'.pdf',sep=''))
			matplot( t(dd[ids,]), main=sprintf('cluster %s: %i genes',i,length(ids)), type='l', col=rgb(0,0,0,0.6), lty=1, lwd=2)
			mtext(side=3, line=-1*(1:length(ids)),text=ids[ order(as.numeric(ids)) ], adj=0.9)
			dev.off()
		}
	}
#} else {
	for(i in 1:length(cls)){
		ids = labels(cls[[i]])
		sz = length(ids)
		cat('dendro',i,sz,'transcripts\n')
		if(sz>maxdend){
			cat('cluster',i,'too big to plot dendrogram (',sz,'genes)\n')
			next
		}
		if(sz<2){
			cat('cluster',i,'too small to plot dendrogram (',sz,'genes)\n')
			next
		}
		fname=paste(denddir,'/pvdend.',i,'.png',sep='')
#		if(file.exists(fname)){
#			cat('skipping',fname,'\n')
#			next
#		}
		png(fname,height=max(baseh*0.8,sz*15),width=baseh)
#		pdf(paste(denddir,'/pvdend.',i,'.pdf',sep=''),height=max(4,length(ids)/4))
		main = sprintf('Cluster %i: dendrogram with pvclust AU p-values',i)
		plot_pvdend(cls[[i]], ordlabs, axes, desc)
		dev.off()
	}
}

if(!exists('uncontained')){
# figure out which pvpicked clusters are not contained within height-based clusters, and create output for them
cat('adding significant pvclust clusters not already contained within height-based clusters\n')
uncontained = which(sapply(clusters, function(cluster){
	cat('.')
	contained = F
	for(cl in cls){
		if(all(cluster %in% labels(cl))){
			contained = T
			break
		}
	}
	!contained
}))
uncontained = clusters[uncontained]
save(uncontained,file='uncontained.clusters.RData')
cat('\n')
}

allcls = lapply(cls, function(x){labels(x)})
# add any uncontained pvclust clusters to the list of clusters
#allcls = c(allcls, uncontained)

#cluster_stats(dd,allcls)

}
