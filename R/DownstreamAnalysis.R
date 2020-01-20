library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

test1 = bitr(vim.cor.os.scde.rank.2[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                    universe = names(vim.cor.os.scde.rank.2), 
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05,
                    readable = TRUE) 
ego_ALL1 <- setReadable(ego_ALL, OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'hsa', 
                 pvalueCutoff = 1)
head(kk,2)


test2 = bitr(vim.cor.rf.scde.rank.2[1:50], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
marker_go <- enrichGO(gene = test1$ENTREZID, 
                    universe = names(rownames(GSE114397.mRNA)), 
                    OrgDb = org.Hs.eg.db, 
                    #keytype = 'ENSEMBL',
                    ont = "ALL", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05,
                    readable = TRUE) 
marker_go1 <- setReadable(marker_go, OrgDb = org.Hs.eg.db)

