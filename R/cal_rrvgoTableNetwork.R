
load('tmp_rrvgo.Rd')

GOTermVector = unique(go.term.input$term_id[go.term.input$source == 'GO:BP'])

library(rrvgo)
library(org.Hs.eg.db)
simMatrix <- calculateSimMatrix(GOTermVector,
                                orgdb = "org.Hs.eg.db",
                                ont = 'BP',
                                method="Rel")
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold = 0.7,
                                orgdb = "org.Hs.eg.db")
###
df.rrgo = cbind(go.term.input,
                reducedTerms[match(go.term.input$term_id, reducedTerms$go),])
as.character(as.data.frame(table(df.rrgo$parentTerm))$Var1)[as.data.frame(table(df.rrgo$parentTerm))$Freq ==1]

###
reducedTerms.used = reducedTerms[!(reducedTerms$term %in% as.character(
  as.data.frame(table(df.rrgo$parentTerm))$Var1)[as.data.frame(table(df.rrgo$parentTerm))$Freq ==1]
),]
df.rrgo = cbind(go.term.input,
                reducedTerms.used[match(go.term.input$term_id, reducedTerms.used$go),])
df.rrgo = df.rrgo[!is.na(df.rrgo$go),]

###
cmm.cas.rrgo = df.rrgo[grep("^CS vs hvCAS$", df.rrgo$Comparisons), 10:17]
cas.cmm.rrgo = df.rrgo[grep("^hvCAS vs CS$", df.rrgo$Comparisons), 10:17]
cmm.cts.rrgo = df.rrgo[grep("^CS vs hvCTS$", df.rrgo$Comparisons), 10:17]
cts.cmm.rrgo = df.rrgo[grep("^hvCTS vs CS$", df.rrgo$Comparisons), 10:17]
cmm.coc.rrgo = df.rrgo[grep("^CS vs hvCOC$", df.rrgo$Comparisons), 10:17]
coc.cmm.rrgo = df.rrgo[grep("^hvCOC vs CS$", df.rrgo$Comparisons), 10:17]
cas.cts.rrgo = df.rrgo[grep("^hvCAS vs hvCTS$", df.rrgo$Comparisons), 10:17]
cts.cas.rrgo = df.rrgo[grep("^hvCTS vs hvCAS$", df.rrgo$Comparisons), 10:17]
cas.coc.rrgo = df.rrgo[grep("^hvCAS vs hvCOC$", df.rrgo$Comparisons), 10:17]
coc.cas.rrgo = df.rrgo[grep("^hvCOC vs hvCAS$", df.rrgo$Comparisons), 10:17]
cts.coc.rrgo = df.rrgo[grep("^hvCTS vs hvCOC$", df.rrgo$Comparisons), 10:17]
coc.cts.rrgo = df.rrgo[grep("^hvCOC vs hvCTS$", df.rrgo$Comparisons), 10:17]

network.input = data.frame(term = c(cmm.cas.rrgo$parentTerm,
                                    cas.cmm.rrgo$parentTerm,
                                    cmm.cts.rrgo$parentTerm,
                                    cts.cmm.rrgo$parentTerm,
                                    cmm.coc.rrgo$parentTerm,
                                    coc.cmm.rrgo$parentTerm,
                                    cas.cts.rrgo$parentTerm,
                                    cts.cas.rrgo$parentTerm,
                                    cas.coc.rrgo$parentTerm,
                                    coc.cas.rrgo$parentTerm,
                                    cts.coc.rrgo$parentTerm,
                                    coc.cts.rrgo$parentTerm),
                           box = c(rep('CS',    length(cmm.cas.rrgo$parentTerm)),
                                   rep('hvCAS', length(cas.cmm.rrgo$parentTerm)),
                                   rep('CS',    length(cmm.cts.rrgo$parentTerm)),
                                   rep('hvCTS', length(cts.cmm.rrgo$parentTerm)),
                                   rep('CS',    length(cmm.coc.rrgo$parentTerm)),
                                   rep('hvCOC', length(coc.cmm.rrgo$parentTerm)),
                                   rep('hvCAS', length(cas.cts.rrgo$parentTerm)),
                                   rep('hvCTS', length(cts.cas.rrgo$parentTerm)),
                                   rep('hvCAS', length(cas.coc.rrgo$parentTerm)),
                                   rep('hvCOC', length(coc.cas.rrgo$parentTerm)),
                                   rep('hvCTS', length(cts.coc.rrgo$parentTerm)),
                                   rep('hvCOC', length(coc.cts.rrgo$parentTerm))),
                           tech = c(rep('CS',    length(cmm.cas.rrgo$parentTerm)),
                                    rep('hvCAS', length(cas.cmm.rrgo$parentTerm)),
                                    rep('CS',    length(cmm.cts.rrgo$parentTerm)),
                                    rep('hvCTS', length(cts.cmm.rrgo$parentTerm)),
                                    rep('CS',    length(cmm.coc.rrgo$parentTerm)),
                                    rep('hvCOC', length(coc.cmm.rrgo$parentTerm)),
                                    rep('hvCAS', length(cas.cts.rrgo$parentTerm)),
                                    rep('hvCTS', length(cts.cas.rrgo$parentTerm)),
                                    rep('hvCAS', length(cas.coc.rrgo$parentTerm)),
                                    rep('hvCOC', length(coc.cas.rrgo$parentTerm)),
                                    rep('hvCTS', length(cts.coc.rrgo$parentTerm)),
                                    rep('hvCOC', length(coc.cts.rrgo$parentTerm))))

dup.column = paste0(network.input$term, "_", network.input$box)
dup.column.se = duplicated(dup.column)
network.input = network.input[!dup.column.se,]

###
network.input.term = data.frame(pterm = c(cmm.cas.rrgo$parentTerm,
                                          cas.cmm.rrgo$parentTerm,
                                          cmm.cts.rrgo$parentTerm,
                                          cts.cmm.rrgo$parentTerm,
                                          cmm.coc.rrgo$parentTerm,
                                          coc.cmm.rrgo$parentTerm,
                                          cas.cts.rrgo$parentTerm,
                                          cts.cas.rrgo$parentTerm,
                                          cas.coc.rrgo$parentTerm,
                                          coc.cas.rrgo$parentTerm,
                                          cts.coc.rrgo$parentTerm,
                                          coc.cts.rrgo$parentTerm),
                                term = c(cmm.cas.rrgo$term,
                                         cas.cmm.rrgo$term,
                                         cmm.cts.rrgo$term,
                                         cts.cmm.rrgo$term,
                                         cmm.coc.rrgo$term,
                                         coc.cmm.rrgo$term,
                                         cas.cts.rrgo$term,
                                         cts.cas.rrgo$term,
                                         cas.coc.rrgo$term,
                                         coc.cas.rrgo$term,
                                         cts.coc.rrgo$term,
                                         coc.cts.rrgo$term),
                                box = c(rep('CS',    length(cmm.cas.rrgo$parentTerm)),
                                        rep('hvCAS', length(cas.cmm.rrgo$parentTerm)),
                                        rep('CS',    length(cmm.cts.rrgo$parentTerm)),
                                        rep('hvCTS', length(cts.cmm.rrgo$parentTerm)),
                                        rep('CS',    length(cmm.coc.rrgo$parentTerm)),
                                        rep('hvCOC', length(coc.cmm.rrgo$parentTerm)),
                                        rep('hvCAS', length(cas.cts.rrgo$parentTerm)),
                                        rep('hvCTS', length(cts.cas.rrgo$parentTerm)),
                                        rep('hvCAS', length(cas.coc.rrgo$parentTerm)),
                                        rep('hvCOC', length(coc.cas.rrgo$parentTerm)),
                                        rep('hvCTS', length(cts.coc.rrgo$parentTerm)),
                                        rep('hvCOC', length(coc.cts.rrgo$parentTerm))),
                                tech = c(rep('CS',    length(cmm.cas.rrgo$parentTerm)),
                                         rep('hvCAS', length(cas.cmm.rrgo$parentTerm)),
                                         rep('CS',    length(cmm.cts.rrgo$parentTerm)),
                                         rep('hvCTS', length(cts.cmm.rrgo$parentTerm)),
                                         rep('CS',    length(cmm.coc.rrgo$parentTerm)),
                                         rep('hvCOC', length(coc.cmm.rrgo$parentTerm)),
                                         rep('hvCAS', length(cas.cts.rrgo$parentTerm)),
                                         rep('hvCTS', length(cts.cas.rrgo$parentTerm)),
                                         rep('hvCAS', length(cas.coc.rrgo$parentTerm)),
                                         rep('hvCOC', length(coc.cas.rrgo$parentTerm)),
                                         rep('hvCTS', length(cts.coc.rrgo$parentTerm)),
                                         rep('hvCOC', length(coc.cts.rrgo$parentTerm))))
library(igraph)
g = graph.data.frame(d = network.input, directed = F)

###pie
pie.list = list()
for(i in c(1:length(names(V(g))))) {
  term.tmp = names(V(g))[i]
  pie.tmp = c(0,0,0,0)
  for(j in c(1:nrow(network.input.term))) {
    if(network.input.term$pterm[j] == term.tmp & network.input.term$tech[j] == 'CS')    { pie.tmp[1] = pie.tmp[1] + 1 }
    if(network.input.term$pterm[j] == term.tmp & network.input.term$tech[j] == 'hvCAS') { pie.tmp[2] = pie.tmp[2] + 1 }
    if(network.input.term$pterm[j] == term.tmp & network.input.term$tech[j] == 'hvCTS') { pie.tmp[3] = pie.tmp[3] + 1 }
    if(network.input.term$pterm[j] == term.tmp & network.input.term$tech[j] == 'hvCOC') { pie.tmp[4] = pie.tmp[4] + 1 }
    if(network.input.term$box[j] == term.tmp & network.input.term$tech[j] == 'CS')    { pie.tmp[1] = pie.tmp[1] + 1 }
    if(network.input.term$box[j] == term.tmp & network.input.term$tech[j] == 'hvCAS') { pie.tmp[2] = pie.tmp[2] + 1 }
    if(network.input.term$box[j] == term.tmp & network.input.term$tech[j] == 'hvCTS') { pie.tmp[3] = pie.tmp[3] + 1 }
    if(network.input.term$box[j] == term.tmp & network.input.term$tech[j] == 'hvCOC') { pie.tmp[4] = pie.tmp[4] + 1 }
  }
  pie.list[[i]] = pie.tmp
}

library(RColorBrewer)
V(g)$pie.color=list(c(brewer.pal(n = 4, name = 'Dark2')))
V(g)$size = rep(8, length(V(g)))
V(g)$size[(length(V(g))-3):length(V(g))] = 10
V(g)$shape = rep('pie', length(V(g)))
V(g)$shape[(length(V(g))-3):length(V(g))] = 'square'

m = "layout_nicely"
set.seed(0)
l <- do.call(m, list(g)) 
plot(g,
     vertex.pie = pie.list,
     vertex.label.dist = 0.8,
     vertex.label.cex = 0.6, 
     vertex.label.degree = pi,
     layout = l)
