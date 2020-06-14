#Step 1: Pre-processing the scRNA-seq data
data_magic <- magic(GSE114397, genes="all_genes")
data_smooth <- as.matrix(data_magic)
GSE114397.mRNA = data_smooth
data_filtered <- data_smooth[rowMeans(data_smooth) > 0.3 & rowMeans(data_smooth > 0) > 0.2,]

#Step 2: Estimating EMT Pseudotime for individual cells
#
vimTime = getPseudotime(GSE114397.mRNA["VIM",])
# 
# EMTtime is calculated by Wanderlust 
#
#Or load pseudotime
load("./data/pseudoTime.rda")
# 
# Check the relaionship of EMT markers and EMTtime
VIM = GSE114397.mRNA["VIM",]
ZEB1 = GSE114397.mRNA["ZEB1",]
KRT19 = GSE114397.mRNA["KRT19",]
KRT7 = GSE114397.mRNA["KRT7",]
cor(emtTime, VIM) #0.4572169
cor(emtTime, ZEB1) #0.776441
cor(emtTime, KRT19) #-0.5158946
cor(emtTime, KRT7) #-0.3828002


#Step 3: Ranking genes by an integrative model
#
###3.1 MAD
mads = apply(data_filtered,1,mad)
index = order(mads, decreasing = T)
mad.ranking= mads[index]
mad.ranking=abs(mad.ranking)/sum(abs(mad.ranking))

###3.2 SDE
sdes <- switchde(data_filtered, as.numeric(vimTime))
sde.filtered = filter(sdes, qval < 0.05)
index = order(abs(sde.filtered$k), decreasing = T)
vim.sdes.rank = sde.filtered[index,]
vim.sdes.ranking = vim.sdes.rank$k
names(vim.sdes.ranking) = vim.sdes.rank$gene
vim.sdes.ranking=abs(vim.sdes.ranking)/sum(abs(vim.sdes.ranking))

sdes <- switchde(data_filtered, as.numeric(emtTime))
sde.filtered = filter(sdes, qval < 0.05)
index = order(abs(sde.filtered$k), decreasing = T)
emt.sdes.rank = sde.filtered[index,]
emt.sdes.ranking = emt.sdes.rank$k
names(emt.sdes.ranking) = emt.sdes.rank$gene
emt.sdes.ranking=abs(emt.sdes.ranking)/sum(abs(emt.sdes.ranking))

###3.3 co-expression network
X_order = data_filtered[,order(vimTime)] 
vim_net = netLeap(data = X_order, cutoff = 0.9) 
g <- graph.data.frame(vim_net, directed=TRUE)
vim.net.ranking = page_rank(g)$vector
index = order(vim.net.ranking, decreasing = T)
vim.net.ranking= vim.net.ranking[index]
vim.net.ranking=vim.net.ranking/sum(vim.net.ranking)

X_order = data_filtered[,order(emtTime)] 
emt_net = netLeap(data = X_order, cutoff = 0.9) 
g <- graph.data.frame(emt_net, directed=TRUE)
emt.net.ranking = page_rank(g)$vector
index = order(emt.net.ranking, decreasing = T)
emt.net.ranking= emt.net.ranking[index]
emt.net.ranking=emt.net.ranking/sum(emt.net.ranking) 

# Step 4. integrative method
vim.ranking = geneRank(mad.ranking,vim.sdes.ranking,vim.net.ranking,a1 = 0.2,a2 = 0.5,a3 = 0.3)
emt.ranking = geneRank(mad.ranking,emt.sdes.ranking,emt.net.ranking,a1 = 0.4,a2 = 0.2,a3 = 0.2)

# Step 5. Cancer prognosis
#
####5.1 test the top 50 genes in single method
res = benchdb(mad.ranking)
res = cbind(res, benchdb(vim.sdes.ranking))
res = cbind(res, benchdb(vim.net.ranking))
res = cbind(res, benchdb(emt.sdes.ranking))
res = cbind(res, benchdb(emt.net.ranking))
write.csv(res, file = "res.csv")

###5.2 Optimized parameters by grid search

vim.grid.search = vector("list", length = 60)
for (gsi in 1:60){
  vim.grid.search[[gsi]] = grid.search(mads,vim.sdes.ranking,vim.net.ranking,gsi)
}
save(vim.grid.search, file = "vim.grid.search.rda")

emt.grid.search = vector("list", length = 60)
for (gsi in 1:60){
  emt.grid.search[[gsi]] = grid.search(mads,emt.sdes.ranking,emt.net.ranking,gsi)
}
save(emt.grid.search, file = "emt.grid.search.rda")

###5.3 Independent test based on the breast cancer signatures
##The identified signatures based on METABRIC dataset
indepTest1(names(emt.ranking)[1:10])

##The benchmark methods
indepTest1(PAM50)
indepTest1(Mamma)
indepTest1(RS)
indepTest1(GGI97)
indepTest1(Endo)
indepTest1(LM)


###Step 5.4 Downsteam analysis
source("./R/DownstreamAnalysis.R")
