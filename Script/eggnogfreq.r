Args <- commandArgs()

# Necessary
intable      <- Args[6]
outDir       <- Args[7]
prefix       <- Args[8]

library(Kitult)
require(pheatmap)
require(ggplot2)
require(methods)
require(dplyr)

#intable = "/Bio2/Project/SGM/SGM1052-Pig-Gut-24-metaGenome-10G.5GAssembly/FunAnno_Result3/06.BasicFun/COG/SG.COG.Anno.txt"
eggnog = read.table(intable,
                    sep = "\t",header = T,
                    quote = "", 
                    #row.names = NULL, 
                    stringsAsFactors = FALSE)

##########################################################################################################
df = data.frame(eggnog)

all <- df %>% count(COG.cat,eggNOG.annot,OGs) %>% mutate(prop = prop.table(n))
eggNOG <- df %>% count(COG.cat,eggNOG.annot) %>% mutate(prop = prop.table(n))
COG <- df %>% count(COG.cat) %>% mutate(prop = prop.table(n))

all = all[order(all$n,decreasing = T),1:4]
write.table(all,file = paste(outDir,paste(prefix,".EggNOG.xls",sep=""),sep="/"),
            row.names=F,sep = "\t",quote=F)


# Sort by 
COG <- COG[order(COG$COG.cat,decreasing = T),]
COG$COG.cat <- factor(COG$COG.cat, levels=unique(COG$COG.cat))

(plot = ggplot(COG,aes(x=COG.cat,y=n,label=n)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    labs(y = "Frequency",
         x = "Ortholog Group"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_COG_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=35, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_COG_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=35, units = "cm")

####################################################################################
# Sort by 
eggNOG <- eggNOG[order(eggNOG$n,decreasing = T),]
eggNOG <- eggNOG[1:50,]
eggNOG <- eggNOG[!is.na(eggNOG$eggNOG.annot),]

eggNOG <- eggNOG[order(eggNOG$COG.cat,decreasing = T),]
eggNOG$eggNOG.annot <- substr(eggNOG$eggNOG.annot, start = 0,stop = 50)
eggNOG$eggNOG.annot <- factor(eggNOG$eggNOG.annot, levels=unique(eggNOG$eggNOG.annot))

(plot = ggplot(eggNOG,aes(x=eggNOG.annot,y=n,fill=COG.cat,label=n)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    labs(y = "Frequency",
         x = "COG"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_EggNOG_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=25, height=35, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_EggNOG_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=25, height=50, units = "cm")

