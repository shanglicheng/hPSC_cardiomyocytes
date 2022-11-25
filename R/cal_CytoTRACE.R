teches = c('cmm', 'cas', 'cts', 'coc')
list.counts = list()
list.inf = list()
list.batch = list()
for(i in c(1:4)) {
  list.counts[[i]] = as.matrix(cardiomyocyte.cell.counts)[, cell.cluster$tech == teches[i]]
  list.inf[[i]] = as.matrix(cell.cluster)[cell.cluster$tech == teches[i],]
  list.batch[[i]] = as.matrix(cell.cluster)[cell.cluster$tech == teches[i], 'batch']
}

library(Linnorm)
list.counts.norm = lapply(list.counts, function(x) Linnorm.Norm(datamatrix = x))

library(plyr)
counts.norm = do.call(cbind, list.counts.norm)
inf.norm = cell.cluster[colnames(counts.norm),]
plot(c(1:nrow(inf.norm)), match(rownames(inf.norm), colnames(counts.norm)))

library(limma)
list.count.edata = lapply(list.counts.norm, function(x) log2(x + 1))

list.count.edata.rbe = mapply(FUN = function(x,y) {removeBatchEffect(x = x, batch = y)}, x = list.count.edata, y = list.batch)

library(plyr)
edata_rbe = do.call(cbind, list.count.edata.rbe)
edata_inf = do.call(rbind, list.inf)

edata_rbe[edata_rbe < 0] = 0
inf.norm = cell.cluster[colnames(edata_rbe),]

#####
library(CytoTRACE)
results2 <- CytoTRACE(edata_rbe, batch = cell.cluster$batch)
