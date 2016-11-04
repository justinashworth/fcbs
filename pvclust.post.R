source('hclust.lib.R')
source('pvclust.from.fcbs.R')
source('plot_pvdend.R')
require(pvclust)

if(!exists('dd')){
	rfile = 'data'
	cat('reading data',rfile,'\n')
	dd = read.table(rfile)
}

cache = c('pvc','hpk','cls')
cache = c(cache, 'uncontained.clusters')

for(f in cache){
	fd = sprintf('%s.RData',f)
	if(file.exists(fd)){
		cat('reading cache file',fd,'\n')
		load(fd)
	}
}

do_pvclust = F

if(do_pvclust){
if(!exists('pvc')){
	pvc = pvclust_from_cpp('1.0/hc')
	save(pvc,file='pvc.RData')
}

hc = pvc$hclust

if(!exists('clusters')){
	pt='au'
	pvcut=0.90
	cat(sprintf('pvpick with P[%s]>=%f\n',pt,pvcut))
	clpick = pvpick(pvc,pv=pt,alpha=pvcut,max.only=T)
	save(clpick,file='pvpick.RData')
	clusters = clpick$clusters
}

} else {
	# no pvclust
	hc = hc_from_fcbs()
}

plots = T
if(plots){

if(!exists('hpk')){
	hpk=get_hpk(hc)
	save(hpk,file='hpk.RData')
}

if(!exists('cls')){
library(dendextendRcpp)
ncl=200
if(!ncl %in% names(hpk)){
	ncls = as.numeric(names(hpk))
	ncl = ncls[ ncls>ncl ][1]
}
cat('Rcpp_cut_lower for',ncl,'clusters\n')
cls = dendextendRcpp::Rcpp_cut_lower(as.dendrogram(hc),height=hpk[[as.character(ncl)]])
save(cls,file='cls.RData')
}

plot_dends=T
if(plot_dends){
	#qs = readLines('qs')
	qs =c()
	desctab = read.delim('desc')
	desc = as.character(desctab$desc)
	names(desc) = as.character(desctab$id)

	library(fastcluster)
	baseh = 480
	maxdend = 2000
	denddir = 'dendro'
	dir.create(denddir)
	expdir = 'exp'
	dir.create(expdir)
	# requires modified version of pvclust to expose functions and work properly with externally produced bootstraps
	axes = hc2axes(hc)
	ordlabs = hc$labels[ hc$order ]

	qsi = sapply(qs, function(x){ which(sapply(cls, function(y){x %in% labels(y)})) } )
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

	for(i in 1:length(cls)){
		ids = labels(cls[[i]])
		sz = length(ids)
		cat('cluster',i,sz,'transcripts\n')
		if(sz>maxdend){
			cat('cluster',i,'too big to plot (',sz,'genes)\n')
			next
		}
		if(sz<2){
			cat('cluster',i,'too small to plot (',sz,'genes)\n')
			next
		}
#		fname=paste(denddir,'/pvdend.',i,'.png',sep='')
#		if(file.exists(fname)){
#			cat('skipping',fname,'\n')
#			next
#		}
#		png(fname,height=max(baseh*0.8,sz*15),width=baseh)
		pdf(sprintf('%s/pvdend.%04d.pdf',denddir,i),height=max(4,length(ids)/4),width=10)
		main = sprintf('Cluster %i: dendrogram with pvclust AU p-values',i)
		plot_pvdend(cls[[i]], ordlabs, axes, desc)
		dev.off()

		# plot expression lineplot
		pdf(sprintf('%s/exp.%04d.pdf',expdir,i),height=8,width=10)
		main = sprintf('Cluster %i: expression values',i)
		matplot( t(dd[ids,]), type='l', lty=1, lwd=2)
		dev.off()

	}
}

do_uncontained = F
if(do_uncontained & !exists('uncontained')){
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
