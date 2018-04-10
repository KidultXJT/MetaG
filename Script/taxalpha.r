Args <- commandArgs()

# Necessary
intable      <- Args[6]
outDir       <- Args[7]
samples      <- Args[8]
groups       <- Args[9]
prefix       <- Args[10] # Taxonomy

library(Kitult)
library(lance)
require(ape)
require(vegan)
require(sprex)
require(pheatmap)
require(RColorBrewer)
require(corrgram)
require(ggplot2)
require(methods)
require(WriteXLS)

# samples = "hsm2-1,ksm2-1,msm2-1,ssm2-1,hsm2-2,ksm2-2,msm2-2,ssm2-2,hsm2-3,ksm2-3,msm2-3,ssm2-3,hsm3-1,ksm3-1,msm3-1,ssm3-1,hsm3-2,ksm3-2,msm3-2,ssm3-2,hsm3-3,ksm3-3,msm3-3,ssm3-3"
# groups  = "hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm"
samples = gsub("-", ".", samples)
groups  = gsub("-", ".", groups)
smplst = strsplit(samples,split = ",")[[1]]
smplst.ord = smplst[order(smplst)]
grplst = strsplit(groups,split = ",")[[1]]
grplst.ord = grplst[order(smplst)]

annotation_row = data.frame(Group = grplst.ord)
rownames(annotation_row) = smplst.ord
map = data.frame(Sample = smplst.ord,
                 Group = grplst.ord)

#intable = "/Bio2/Project/SGM/SGM1052-Pig-Gut-24-metaGenome-10G.5GAssembly/BinGenome_Result3/TaxAnalysis/Reads/Family.xls"
table = data.frame(t(read.table(intable,sep = "\t",header = T,row.names = 1,quote = '')))
table = table[order(rownames(table)),]

## Alpha Diversity
(obs_sp     <- specnumber(table)) # Species ~= Obs
(evenness   <- diversity(x = table)/(specnumber(table)))
#(richness   <- specaccum(table))
#(fishalpha  <- fisher.alpha(table))
(shannon    <- diversity(x = table,index = "shannon",MARGIN = 1))
(simpson    <- diversity(x = table,index = "simpson",MARGIN = 1))
(invsimpson <- diversity(x = table,index = "invsimpson",MARGIN = 1))
# Alpha INDEX
(INDEX = data.frame(Obs_sp = obs_sp,
                    Shannon = shannon,
                    Simpson = simpson))
(richness.s <- with(map, specpool(table, Sample, smallsample = F)))
aceTable <- estimateR(table)
(AlphaINDEX = data.frame(Species=INDEX[,1],Chao1=richness.s$chao,INDEX[,2:3]))

## Rarecurve
c = c("#001F3F","#FF5722","#51710A","#4E1184","#F87D09","#C5D200",
      "#932B77","#0D63A5","#F9C535","#9BDF46","#FD367E","#228896",
      "#F8DA5B","#83CC61","#E84A5F","#3D84A8","#F8C957","#B7E576",
      "#FF847C","#46CDCF","#99CDA9","#FECEA8","#ABEDD8","#D1E9D2",
      "#FFD6A4","#A7CDCC","#E5F4E7","#FDE9DF","#D6E6F2","#F1FDF3",
      "#43496E","#544D7E","#65589C","#3D3551","#884EA2","#3EC280",
      "#D15400","#EDEDEE","#E16A6B","#D9DEE1","#EA9432","#36D6B5",
      "#68C2A1","#AAB6B6","#2474A8","#3597DA","#E57E22","#81CEE0",
      "#D54542","#D24D58","#51B2D7","#F5AA36","#AEA7D1","#00B069",
      "#EE4836","#E77E04","#F6C917","#1F3A92","#2ABA99","#169F85",
      "#F5D66D","#BDC2C6","#F17835","#1E8AC2","#F44746","#A0DDCE",
      "#66CB99","#59AAE2","#903D87","#D81E17","#8C44AC","#4ECCC3",
      "#F22613","#DA0A5B","#9959B4","#95281B","#03C8A7","#F1F0EE",
      "#BEBEBE","#F62459","#94A4A4","#C7F5C4","#D1D6D3","#F89306",
      "#EA964E","#65C5BA","#5C96BE","#23A6F0")
C = 1

for (i in 1:length((grplst.ord))){
  r = rarecurve(table[i,], start = 0 ,step = 1000000, sample = raremax, col = c[C])
  par(new=T)
}

# RankAbundance
plot(radfit(table[3,]))
