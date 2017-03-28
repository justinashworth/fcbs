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

do_pvclust = T

if(do_pvclust){
if(!exists('pvc')){
	pvc = pvclust_from_cpp('1.0/hc')
	save(pvc,file='pvc.RData')
}

hc = pvc$hclust

if(!exists('clusters')){
	pt='au'
	pvcut=0.80
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
	hpk=get_hpk(hc,0.8)
	save(hpk,file='hpk.RData')
}

if(!exists('cls')){
library(dendextendRcpp)
ncl=750
if(!ncl %in% names(hpk)){
	ncls = as.numeric(names(hpk))
	ncl = ncls[ ncls>ncl ][1]
}
cat('Rcpp_cut_lower for',ncl,'clusters\n')
cls = dendextendRcpp::Rcpp_cut_lower(as.dendrogram(hc),height=hpk[[as.character(ncl)]])
save(cls,file='cls.RData')
}

desctab = read.delim('desc')
clstab = as.data.frame( do.call(rbind, lapply(1:length(cls), function(i){ cbind(i,labels(cls[[i]]))})))
names(clstab) = c('cl','id')
pvctab = as.data.frame( do.call(rbind, lapply(1:length(clpick$clusters), function(i){ cbind(i,clpick$clusters[[i]])})))
names(pvctab) = c('pvc','id')
m = merge(clstab, pvctab, by='id', all.x=T, sort=F)
m = merge(m, desctab, by='id', all.x=T, sort=F)
hc_for_id = as.numeric(as.character(m$cl))
names(hc_for_id) = as.character(m$id)
m = m[order(m$cl,m$pvc,m$id),]
tsv(m,'cls.tsv')

#warning: tapply sorts by as.character(index) here:
clids = tapply(m$id,m$cl,list)
# sort numerically so list index == cluster index
clids = clids[ order(as.numeric(names(clids))) ]

# convenient lists of image files corresponding to larger bootstrap clusters
pvclens = sapply(clpick$clusters, length)
minlen = 9
cat( sprintf('pvc/exp.%04d.pdf', which( pvclens > minlen )), sep='\n', file='pvc_big_exp')
hc_for_pvc = m$cl
names(hc_for_pvc) = m$pvc
hc_big_pvcs = unique( hc_for_pvc[ pvclens > minlen ] )
cat( sprintf('dendro/dendro.%04d.pdf', hc_big_pvcs), sep='\n', file='hc_big_pvcs_dendro')
cat( sprintf('exp/exp.%04d.pdf', hc_big_pvcs), sep='\n', file='hc_big_pvcs_exp')

plot_exp = function(values,fname,log=T,...){
	values = values[ order(rownames(values)), ]
	pdf(fname,height=8,width=16,colormodel='cmyk')
	if(log) values = log(values,10)
	par(mar=c(8,4,4,5))
	cols = rainbow(nrow(values))
	yrange=range(values,na.rm=T)
	if(all(yrange==0)){
		cat('all-zero cluster...\n')
	} else {
		matplot( t(values), type='n', xaxt='n',, xaxs='i', yaxt='n',...)
		abline(v=1:ncol(values), lty=1, lwd=0.5, col=rgb(0.9,0.9,0.9))
		matlines( t(values), lty=1, lwd=2, col=cols)
		axis(1,at=1:ncol(values),labels=NA)
		text(1:ncol(values), par("usr")[3]-0.15, srt=60, adj=1, xpd=TRUE, labels=colnames(values), cex=0.4)
		ylv = c(0.01,0.1,1,10,100,1000,10000)
		ylab=as.character(ylv)
		if(log) ylv=log(ylv,10)
		axis(2,at=ylv,labels=ylab,las=2)
		for(i in 1:nrow(values)){
			mtext(rownames(values)[i], line=(-0.7*i)-2, col=cols[i],adj=0.97,outer=T,cex=0.7)
		}
	}
	dev.off()
}

plot_dends=T
if(plot_dends){
	qs =c()
	if(file.exists('qs')) qs = readLines('qs')
	desc = paste(desctab$id,desctab$desc)
	names(desc) = as.character(desctab$id)

	library(fastcluster)
	baseh = 480
	maxdend = 2000
	denddir = 'dendro'
	dir.create(denddir)
	expylab = "transcript expression"
	expdir = 'exp'
	dir.create(expdir)
	dir.create('pvc')
	# requires modified version of pvclust to expose functions and work properly with externally produced bootstraps
	axes = hc2axes(hc)
	ordlabs = hc$labels[ hc$order ]

	# make extra plots for specific queries
	qsi = sapply(qs, function(x){ which(sapply(cls, function(y){x %in% labels(y)})) } )
	pvi = sapply(qs, function(x){ which(sapply(clpick$clusters, function(y){x %in% y})) } )
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

		pdf(sprintf('%s/pvdend.%04d.pdf',denddir,i),height=max(4,length(ids)/4),width=12)
		main = sprintf('Cluster %i: dendrogram with pvclust AU p-values',i)
		plot_pvdend(cls[[i]], ordlabs, axes, desc)
		dev.off()

		# plot expression lineplot
		main = sprintf('Cluster %i: expression values',i)
		values = as.matrix(dd[ids,])
		rownames(values)=ids
		fname = sprintf('%s/exp.%04d.pdf',expdir,i)

		plot_exp(values,fname,main=main,ylab=expylab,xlab='')

	}

	# now pvclust pvpicked cluster expression plots
	for(i in 1:length(clpick$clusters)){
		cat( sprintf('pvc cluster %i\n',i) )

		fname = sprintf('pvc/exp.%04d.pdf',i)
		ids = clpick$clusters[[i]]
		values = as.matrix(dd[ids,])
		rownames(values)=ids
		main = sprintf('Bootstrap Cluster %i',i)
		plot_exp(values,fname,main=main,ylab=expylab,xlab='')
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
