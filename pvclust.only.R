load('pvc.RData')
load('pvpick.RData')

baseh = 480
mindend = 15
maxdend = 200
denddir = 'dendro'
desc=NULL
dir.create(denddir)
# requires modified version of pvclust to expose functions and work properly with externally produced bootstraps
axes = hc2axes(pvc$hclust)
ordlabs = pvc$hclust$labels[ pvc$hclust$order ]

require('dendindexed.R')
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

			#png(paste(denddir,'/pvdend.',pickindex,'.png',sep=''),width=baseh,height=max(baseh*0.8,sz*15))
			#plot_pvdend(x, ordlabs, axes, desc, main=main)
			#dev.off()

			pdf(paste(denddir,'/pvdend.',pickindex,'.pdf',sep=''),width=10,height=max(4,sz/4))
			plot_pvdend(x, ordlabs, axes, desc, main=main)
			dev.off()

		}
	}
}))
