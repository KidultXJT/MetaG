Args <- commandArgs()

# Necessary
intable      <- Args[6]
outDir       <- Args[7]
prefix       <- Args[8]

library(Kitult)
require(pheatmap)
require(ggplot2)
require(methods)

#intable = "/Bio2/Project/SGM/SGM1059-Environment-Soil_Feces-12-metaGenome-10G/Taxonomy_Result/TaxAnalysis/AvgDepth/RelativeAbundance/Kingdom.RelativeAbundance.xls"
table = read.table(intable,sep = "\t",header = T,row.names = 1,quote = '')

# Heatmap
hmtable = table
#ind <- apply(hmtable, 1, sd) == 0
#hmtable <- hmtable[!ind,]
#ind <- apply(hmtable, 2, sd) == 0
#hmtable <- hmtable[,!ind]
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "none",
            filename = paste(outDir,paste(prefix,"_Heatmap.png",sep=""),sep="/"))
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "column",
            filename = paste(outDir,paste(prefix,"_cHeatmap.png",sep=""),sep="/"))
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "row",
            filename = paste(outDir,paste(prefix,"_rHeatmap.png",sep=""),sep="/"))
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "none",
            filename = paste(outDir,paste(prefix,"_Heatmap.pdf",sep=""),sep="/"))
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "column",
            filename = paste(outDir,paste(prefix,"_cHeatmap.pdf",sep=""),sep="/"))
Kitult::phm(hmtable,brwCol = "RdYlGn",scale = "row",
            filename = paste(outDir,paste(prefix,"_rHeatmap.pdf",sep=""),sep="/"))

# Dominant
DominantTable = table[apply(table,1,max) > 0.01,]
outTable = cbind(SampleID = rownames(DominantTable),data.frame(DominantTable))
write.table(outTable,paste(outDir,paste(prefix,"_Dominant.xls",sep=""),sep="/"),sep = "\t",quote = FALSE,row.names = F)

hmDominantTable = DominantTable[rownames(DominantTable) !="",]
if (nrow(hmDominantTable) >= 2) {
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none",
              filename = paste(outDir,paste(prefix,"_Dominant_Heatmap.png",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "column",
              filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.png",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "row",
              filename = paste(outDir,paste(prefix,"_Dominant_rHeatmap.png",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none",
              filename = paste(outDir,paste(prefix,"_Dominant_Heatmap.pdf",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "column",
              filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.pdf",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "row",
              filename = paste(outDir,paste(prefix,"_Dominant_rHeatmap.pdf",sep=""),sep="/"))
}else{
  if (nrow(hmDominantTable) < 2) {
    Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none",cluster_rows = F,main="Scale::None",
                filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.png",sep=""),sep="/"))
    Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none",cluster_rows = F,main="Scale::None",
                filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.pdf",sep=""),sep="/"))
  }else{
    Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "column",cluster_rows = F,
                filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.png",sep=""),sep="/"))
    Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "column",cluster_rows = F,
                filename = paste(outDir,paste(prefix,"_Dominant_cHeatmap.pdf",sep=""),sep="/"))
  }
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none", cluster_rows = F,
              filename = paste(outDir,paste(prefix,"_Dominant_Heatmap.png",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "row",cluster_rows = F,
              filename = paste(outDir,paste(prefix,"_Dominant_rHeatmap.png",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "none", cluster_rows = F,
              filename = paste(outDir,paste(prefix,"_Dominant_Heatmap.pdf",sep=""),sep="/"))
  Kitult::phm(hmDominantTable,brwCol = "RdYlGn",scale = "row",cluster_rows = F,
              filename = paste(outDir,paste(prefix,"_Dominant_rHeatmap.pdf",sep=""),sep="/"))
}

# 
Others = t(data.frame(1 - colSums(hmDominantTable)))
rownames(Others) = "Others"
stackDominantTable = rbind(Others,hmDominantTable)
DF = data.frame(Name=rownames(stackDominantTable),stack(stackDominantTable))
DF$Name = factor(DF$Name, levels=unique(rownames(stackDominantTable)))
stackplot = ggplot(DF,aes(x=ind,y=values,color=Name,fill=Name)) + 
  Kitult::ggBasicStack() +
  Kitult::fillScale(nrow(stackDominantTable)) +
  Kitult::colScale(nrow(stackDominantTable)) + 
  Kitult::reportStyle() + 
  theme(axis.text.x = element_text(angle=50,hjust=1)) + 
  labs(x="Samples",y="Relative Abundance") +
  guides(fill = guide_legend(ncol = 1))
ggsave(stackplot,
       filename=paste(outDir,paste(prefix,"_Stack.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=15, units = "cm")
ggsave(stackplot,
       filename=paste(outDir,paste(prefix,"_Stack.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=15, units = "cm")

DF = DF[DF$Name != "Others",]
boxplot = ggplot(DF,aes(Name,values)) + 
  Kitult::ggBasicBox() +
  Kitult::fillScale(nrow(stackDominantTable)) +
  Kitult::colScale(nrow(stackDominantTable)) + 
  Kitult::reportStyle() + 
  theme(axis.text.x = element_text(angle=50,hjust=1)) + 
  labs(x="",y="")
ggsave(boxplot,
       filename=paste(outDir,paste(prefix,"_Box.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=15, units = "cm")
ggsave(boxplot,
       filename=paste(outDir,paste(prefix,"_Box.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=15, units = "cm")

