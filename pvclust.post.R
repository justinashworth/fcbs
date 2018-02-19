source('hclust.lib.R')
source('pvclust.from.fcbs.R')
source('plot_pvdend.R')
source('dendindexed.R')
require(pvclust)
require(fastcluster)

mindend = 2
maxdend = 2000

expylab = "Transcript level (CPM)"

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

do_pvclust = T

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

# for plotting dendrograms with pvclust information
if(exists('clusters')){
	# requires modified version of pvclust to expose functions and work properly with externally produced bootstraps
	axes = hc2axes(hc)
	ordlabs = hc$labels[ hc$order ]
}

plots = T
if(plots){

if(!exists('hpk')){
	hpk=get_hpk(hc,0.8)
	save(hpk,file='hpk.RData')
}

if(!exists('cls')){
require(dendextendRcpp)
ncl=1000
if(!ncl %in% names(hpk)){
	ncls = as.numeric(names(hpk))
	ncl = ncls[ ncls>ncl ][1]
}
cat('Rcpp_cut_lower for',ncl,'clusters\n')
cls = dendextendRcpp::Rcpp_cut_lower(as.dendrogram(hc),height=hpk[[as.character(ncl)]])
save(cls,file='cls.RData')
}

clids = sapply(cls, function(x){labels(x)})

exp_plot = function(expvals, desc=NULL, ...){
	log='y'
	par(mar=c(10,4,6,3))
	ids = rownames(expvals)
	cols = rainbow(length(ids))
	matplot( t(expvals), type='n', xaxt='n', xaxs='i', log=log, ...)
	matlines( t(expvals), lty=1, lwd=2, col=cols)
	axis(1, las=2, at=1:ncol(expvals), labels=colnames(expvals), cex.axis=0.5)
	if(!is.null(desc)){
		for(i in 1:length(ids)){
			mtext(ids[i], line=-1*i, col=cols[i])
		}
	}
}

plot_dends=T
if(plot_dends){
	desctab = read.delim('desc')
	desc = paste(desctab$id,desctab$desc)
	names(desc) = as.character(desctab$id)

	baseh = 480
	denddir = 'dendro'
	dir.create(denddir)
	expdir = 'exp'
	dir.create(expdir)

	# looping over non-pvc height-based hc clusters here:
	for(i in 1:length(cls)){
		ids = labels(cls[[i]])
		sz = length(ids)
		cat('cluster',i,sz,'transcripts\n')
		if(sz>maxdend){
			cat('cluster',i,'too big to plot (',sz,'genes)\n')
			next
		}
		if(sz<mindend){
			cat('cluster',i,'too small to plot (',sz,'genes)\n')
			next
		}

		# dendrogram for height-based hc clusters, including pvclust information at branchpoints
		pdf(sprintf('%s/pvdend.%04d.pdf',denddir,i),height=max(4,length(ids)/4),width=10,colormodel='cmyk')
		main = sprintf('Cluster %i: dendrogram with pvclust AU p-values',i)
		plot_pvdend(cls[[i]], ordlabs, axes, desc)
		dev.off()

		# plot expression lineplot
		pdf(sprintf('%s/exp.%04d.pdf',expdir,i),height=8,width=20,colormodel='cmyk')
		main = sprintf('Cluster %i: expression values',i)
		exp_plot(dd[ids,],desc,main=main,ylab=expylab)
		dev.off()
	}
}

# plot out all of the actual pvclust clusters themselves, agnostic of the hc height-based clustering
plot_all_pvclusters=T
if(plot_all_pvclusters){
cat('Plotting all PVClust clusters\n')

pvcdenddir = 'pvc_dendro'
pvcexpdir = 'pvc_exp'

dir.create(pvcdenddir)
dir.create(pvcexpdir)

dnd = as.dendrogram.hclust.indexed(pvc$hclust)
invisible(dendrapply(dnd, function(x){
	index=attr(x,'index')
	if(!is.leaf(x)){
		if(index %in% clpick$edges){

			ids = labels(x)
			sz = length(ids)
			pickindex = which(clpick$edges==index)
			cat('dendro',pickindex,sz,'transcripts\n')
			if(sz<mindend){cat('cluster',pickindex,'too small (',sz,'genes)\n'); return()}
			if(sz>maxdend){cat('cluster',pickindex,'too big (',sz,'genes)\n'); return()}
			main = sprintf('Cluster %i: dendrogram with pvclust AU p-values',pickindex)

			# dendrogram
			pdf(sprintf('%s/pvcdend.%04d.pdf',pvcdenddir,pickindex),width=10,height=max(4,sz/4),colormodel='cmyk')
			plot_pvdend(x, ordlabs, axes, desc, main=main)
			dev.off()

			# exp lineplot
			pdf(sprintf('%s/pvcexp.%04d.pdf',pvcexpdir,pickindex),height=8,width=20,colormodel='cmyk')
			main = sprintf('PVClust %i: expression values',pickindex)
			exp_plot(dd[ids,],desc,main=main,ylab=expylab)
			dev.off()
		}
	}
}))
}

do_uncontained = T
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

} # plots
