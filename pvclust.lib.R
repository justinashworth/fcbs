source('~/code/fcbs/pvclust.from.fcbs.R')
pvc = pvclust_from_cpp('1.0/hc')
save(pvc,file='pvc.RData')
