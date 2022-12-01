mixed.cell.counts.cmm.01 = readRDS("~/cmm01_allcounts.Rds")#b1
mixed.cell.counts.cmm.02 = readRDS("~/cmm02_allcounts.Rds")#b2

mixed.cell.counts.cas.01 = readRDS("~/cas01_allcounts.Rds")#b1
mixed.cell.counts.cas.02 = readRDS("~/cas02_allcounts.Rds")#b3

mixed.cell.counts.cts.01 = readRDS("~/cts01_allcounts.Rds")#b1
mixed.cell.counts.cts.02 = readRDS("~/cts02_allcounts.Rds")#b2
mixed.cell.counts.cts.03 = readRDS("~/OneDrive - KI.SE/R1/0900_allcounts/cts03_allcounts.Rds")#b2

mixed.cell.counts.coc.01 = readRDS("~/coc01_allcounts.Rds")#b1
mixed.cell.counts.coc.02 = readRDS("~/coc02_allcounts.Rds")#b1
mixed.cell.counts.coc.03 = readRDS("~/coc03_allcounts.Rds")#b3

list.batch.all.umi.obj.matrix = list(mixed.cell.counts.cmm.01,
                                     mixed.cell.counts.cmm.02,
                                     
                                     mixed.cell.counts.cas.01,
                                     mixed.cell.counts.cas.02,
                                     
                                     mixed.cell.counts.cts.01,
                                     mixed.cell.counts.cts.02,
                                     mixed.cell.counts.cts.03,
                                     
                                     mixed.cell.counts.coc.01,
                                     mixed.cell.counts.coc.02,
                                     mixed.cell.counts.coc.03)

count.matrix <- as.matrix(list.batch.all.umi.obj.matrix[[1]])
for (i in 2:length(list.batch.all.umi.obj.matrix)) {
  count.matrix <- merge(count.matrix, as.matrix(list.batch.all.umi.obj.matrix[[i]]), by = "row.names", all=T)
  rownames(count.matrix) <- count.matrix$Row.names
  count.matrix <- count.matrix[ , !(names(count.matrix) %in% "Row.names")]  
}

dim(count.matrix)
count.matrix[is.na(count.matrix)] <- 0

mixed.cell.counts = count.matrix

subject = colnames(mixed.cell.counts)
subject = gsub(pattern = '.*_', replacement = '', x = subject)
subject = gsub(pattern = 'data.', replacement = '', x = subject)
subject = gsub(pattern = '.f.hes.', replacement = '', x = subject)

tech = gsub(pattern = '\\..*', replacement = '', x = subject)
cell.cluster = data.frame(subject,
                          tech,
                          row.names = colnames(mixed.cell.counts),
                          stringsAsFactors = F)

cell.cluster$batch = cell.cluster$subject
cell.cluster$batch = gsub(pattern = 'cmm.01', replacement = 'b1', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'cmm.02', replacement = 'b2', x = cell.cluster$batch)

cell.cluster$batch = gsub(pattern = 'cas.01', replacement = 'b1', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'cas.02', replacement = 'b3', x = cell.cluster$batch)

cell.cluster$batch = gsub(pattern = 'cts.01', replacement = 'b1', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'cts.02', replacement = 'b2', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'cts.03', replacement = 'b2', x = cell.cluster$batch)

cell.cluster$batch = gsub(pattern = 'coc.01', replacement = 'b1', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'coc.02', replacement = 'b1', x = cell.cluster$batch)
cell.cluster$batch = gsub(pattern = 'coc.03', replacement = 'b3', x = cell.cluster$batch)

###
geneY = read.table(file = '~/ygene.txt', header = T, row.names = NULL, sep = "\t", as.is = T)
geneY = geneY[!duplicated(geneY$Symbol),]

geneY = geneY[geneY$Symbol %in% rownames(mixed.cell.counts),]

df.y = data.frame(expY = apply(mixed.cell.counts[geneY$Symbol,], 2, sum), 
                  cell.cluster)
df.y$female = rep('F', nrow(df.y))
df.y$female[df.y$expY > 1] = 'M'

table(df.y$female, df.y$tech)


df.y = df.y[order(df.y$female, decreasing = T),]
df.y$tech = factor(x = df.y$tech,
                   levels = c('cmm', 'cas', 'cts', 'coc'),
                   labels = c('CS.Mixed', 'hvCAS.Mixed', 'hvCTS.Mixed', 'hvCOC.Mixed'))

library(ggplot2)
library(RColorBrewer)
library(scales)
g = ggplot(df.y, mapping = aes(x = c(1:nrow(df.y)), y = expY, color = tech, fill = tech))
g = g + geom_bar(stat = 'identity') +
  scale_color_manual(values = brewer.pal(n = 4, name = 'Dark2')) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'Dark2')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) + 
  theme(legend.position="bottom",
        legend.title = element_blank()) +
  xlab('Mixed cells') + 
  ylab('Sum counts of \n Y-chromosome genes') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
g
ggsave(filename = paste0(outIndex, "ExpOfYchromosomeGenesInMixedCells.pdf"), plot = g, width = 6, height = 3, useDingbats=FALSE)


df.y.table = as.data.frame(table(df.y$female, df.y$tech))
colnames(df.y.table) = c('sex', 'tech', 'freq')

g = ggplot(df.y.table, mapping = aes(x = tech, y = sex, fill = log10(freq+1)))
g = g + geom_tile(color = 'white') +
  scale_fill_gradient2(low = "blue", mid = "white" ,high = "red", 
                       midpoint = log10(max(df.y.table$freq)+1)/2, limit = c(0,log10(max(df.y.table$freq)+1)), space = "Lab", 
                       name="Number of\nMixed cells")
g = g + geom_text(aes(tech, sex, label = freq), color = "black", size = 6)
g = g +  theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               panel.background = element_blank()) + 
  theme(legend.position="bottom",
        legend.text = element_blank()) +
  xlab('Mixed cells from \n different technologies') + 
  ylab('Cell sex')
g
ggsave(filename = paste0(outIndex, "NumberOfFemaleMixedCells.pdf"), plot = g, width = 5, height = 3.5, useDingbats=FALSE)



df.y.table = as.data.frame(table(df.y$female, df.y$subject))
colnames(df.y.table) = c('sex', 'subject', 'freq')

g = ggplot(df.y.table, mapping = aes(x = subject, y = sex, fill = log10(freq+1)))
g = g + geom_tile(color = 'white') +
  scale_fill_gradient2(low = "blue", mid = "white" ,high = "red", 
                       midpoint = log10(max(df.y.table$freq)+1)/2, limit = c(0,log10(max(df.y.table$freq)+1)), space = "Lab", 
                       name="Number of\nMixed cells")
g = g + geom_text(aes(subject, sex, label = freq), color = "black", size = 6)
g = g +  theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               panel.background = element_blank()) + 
  theme(legend.position="bottom",
        legend.text = element_blank()) +
  xlab('Mixed cells from \n different samples') + 
  ylab('Cell sex')
g
ggsave(filename = paste0(outIndex, "NumberOfFemaleMixedCells_subject.pdf"), plot = g, width = 10, height = 3.5, useDingbats=FALSE)


female.cells.fib = read.table(file = '~/df.y.txt', header = T, row.names = 1, sep = '\t')
female.cells.car = read.table(file = '~/df.plot.txt', header = T, row.names = 1, sep = '\t')

table(rownames(female.cells.fib) %in% rownames(female.cells.car))

table(female.cells.fib$female)

df.y$type = rep('Other', nrow(df.y))
df.y$type[rownames(df.y) %in% rownames(female.cells.fib[female.cells.fib$female == 'M',])] = 'Fibroblasts.Male'
df.y$type[rownames(df.y) %in% rownames(female.cells.fib[female.cells.fib$female == 'F',])] = 'Fibroblasts.Female'

df.y$type[rownames(df.y) %in% rownames(female.cells.car)] = 'Cardiomyocytes'

table(df.y$type)
df.y$type = factor(x = df.y$type,
                   levels = c('Cardiomyocytes',    'Fibroblasts.Female', 'Fibroblasts.Male',          'Other'),
                   labels = c('Cardiomyocytes',    'Fibroblasts.Female', 'Fibroblasts.Male',          'Other'))

write.table(x = df.y, file = '~/df.y.txt', sep = '\t')

library(ggplot2)
library(RColorBrewer)
library(scales)
g = ggplot(df.y, mapping = aes(x = c(1:nrow(df.y)), y = expY, color = tech, fill = tech))
g = g + geom_bar(stat = 'identity') +
  scale_color_manual(values = brewer.pal(n = 4, name = 'Dark2')) +
  scale_fill_manual(values = brewer.pal(n = 4, name = 'Dark2')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) + 
  theme(legend.position="bottom",
        legend.title = element_blank()) +
  xlab('Mixed cells') + 
  ylab('Sum counts of \n Y-chromosome genes') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
g = g + facet_wrap(~ type, nrow = 4)
g
ggsave(filename = paste0(outIndex, "ExpOfYchromosomeGenesInMixedCells_Type.pdf"), plot = g, width = 8, height = 10, useDingbats=FALSE)



###Female cells

plot(c(1:nrow(cell.cluster)), match(rownames(df.y), colnames(mixed.cell.counts)))
df.y = df.y[colnames(mixed.cell.counts),]

mixed.cell.counts = mixed.cell.counts[,rownames(df.y[df.y$female == 'F',])]
cell.cluster = cell.cluster[colnames(mixed.cell.counts),]

###
df.y.female = df.y[rownames(cell.cluster),]
table(df.y.female$female)

cell.cluster$sex = df.y.female$female
cell.cluster$celltype = df.y.female$type

table(cell.cluster$tech, cell.cluster$subject)

data.cmm.f.hes2.01 = mixed.cell.counts[,cell.cluster$subject == 'cmm.01']
data.cmm.f.hes2.02 = mixed.cell.counts[,cell.cluster$subject == 'cmm.02']

data.cas.f.hes2.01 = mixed.cell.counts[,cell.cluster$subject == 'cas.01']
data.cas.f.hes2.02 = mixed.cell.counts[,cell.cluster$subject == 'cas.02']

data.cts.f.hes2.01 = mixed.cell.counts[,cell.cluster$subject == 'cts.01']
data.cts.f.hes2.02 = mixed.cell.counts[,cell.cluster$subject == 'cts.02']
data.cts.f.hes2.03 = mixed.cell.counts[,cell.cluster$subject == 'cts.03']

data.coc.f.hes2.01 = mixed.cell.counts[,cell.cluster$subject == 'coc.01']
data.coc.f.hes2.02 = mixed.cell.counts[,cell.cluster$subject == 'coc.02']
data.coc.f.hes2.03 = mixed.cell.counts[,cell.cluster$subject == 'coc.03']


list.data = list(data.cmm.f.hes2.01 = data.cmm.f.hes2.01,
                 data.cmm.f.hes2.02 = data.cmm.f.hes2.02,
                 
                 data.cas.f.hes2.01 = data.cas.f.hes2.01,
                 data.cas.f.hes2.02 = data.cas.f.hes2.02,
                 
                 data.cts.f.hes2.01 = data.cts.f.hes2.01,
                 data.cts.f.hes2.02 = data.cts.f.hes2.02,
                 data.cts.f.hes2.03 = data.cts.f.hes2.03,
                 
                 data.coc.f.hes2.01 = data.coc.f.hes2.01,
                 data.coc.f.hes2.02 = data.coc.f.hes2.02,
                 data.coc.f.hes2.03 = data.coc.f.hes2.03)


###
myCreateSeuratObjectFromReadMatrix = function(ReadMatrix, projectName) {
  ctrl.data = ReadMatrix
  colnames(ctrl.data) = paste0(colnames(ctrl.data))
  ctrl <- CreateSeuratObject(counts = ctrl.data, project = projectName, min.cells = 0, min.features = 0)
  
  ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
  
  ctrl <- NormalizeData(ctrl)
  ctrl <- ScaleData(ctrl, display.progress = F)
  return(ctrl)
}
list.data.obj = mapply(function(x, y) myCreateSeuratObjectFromReadMatrix(x, y), list.data, names(list.data))

immune.anchors <- FindIntegrationAnchors(object.list = list.data.obj, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

immune.combined.tmp = immune.combined

##
immune.combined = immune.combined.tmp
DefaultAssay(immune.combined) = "integrated"

# Run the standard workflow for visualization and clustering
immune.combined = ScaleData(immune.combined, verbose = FALSE)
immune.combined = RunPCA(immune.combined, npcs = 50, verbose = FALSE)

immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20, seed.use = 0)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.1)

tech = as.vector(immune.combined$orig.ident)
tech = gsub(pattern = 'data\\.', replacement = '', x = tech)
tech = gsub(pattern = 'f\\.', replacement = '', x = tech)
tech = gsub(pattern = '\\..*', replacement = '', x = tech)
table(tech)

subject = as.vector(immune.combined$orig.ident)
subject = gsub(pattern = '.*_', replacement = '', x = subject)
subject = gsub(pattern = 'data.', replacement = '', x = subject)
subject = gsub(pattern = '.f.hes.', replacement = '', x = subject)
table(subject)

emb.pca = immune.combined@reductions$pca@cell.embeddings
emb.tsne = immune.combined@reductions$tsne@cell.embeddings
emb.umap = immune.combined@reductions$umap@cell.embeddings
