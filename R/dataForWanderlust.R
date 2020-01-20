library(CancerSubtypes)
library(dplyr)

Mesg = filter(EMTgenelist, EMTtype == "Mes")
index = which(Mesg$gene %in% rownames(GSE114397.mRNA))
temp = GSE114397.mRNA[index,]
dmes = CancerSubtypes::FSbyMAD(temp,cut.type = "topk", 10)
smes = colSums(dmes)
osmesi = smes[order(smes)]

Epi = filter(EMTgenelist, EMTtype == "Epi")
index = which(Epi$gene %in% rownames(GSE114397.mRNA))
temp = GSE114397.mRNA[index,]
depi = CancerSubtypes::FSbyMAD(temp,cut.type = "topk", 10)
sepi = colSums(depi)
osepid = sepi[order(sepi,decreasing = T)]

startPiont = intersect(names(osmesi)[1:3000],names(osepid)[1:3000])
endPiont = intersect(names(osmesi)[4200:7523],names(osepid)[4200:7523])
print(startPiont)

data = rbind(dmes,depi)
write.table(t(data), file = "EMT20.csv",sep = ",",row.names = FALSE,col.names = FALSE)

##run wanderlust in matlab
