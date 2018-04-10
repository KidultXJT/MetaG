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

#intable = "/Bio2/Project/SGM/SGM1052-Pig-Gut-24-metaGenome-10G.5GAssembly/FunAnno_Result3/07.EnhanceFun/CAZy/SG.CAZy.xls"
cazy = read.table(intable,
                  sep = "\t",
                  header = F,
                  #row.names = NULL, 
                  stringsAsFactors = FALSE)

##########################################################################################################
df = data.frame(cazy)

all <- df %>% count(V2,V5,V6,V7) %>% mutate(prop = prop.table(n))
family <- df %>% count(V2,V5) %>% mutate(prop = prop.table(n))
class <- df %>% count(V5) %>% mutate(prop = prop.table(n))

all = all[order(all$n,decreasing = T),1:4]
colnames(all) = c("CAZy Class","CAZy Family","ECs","Freq")
write.table(all,file = paste(outDir,paste(prefix,".Freq.xls",sep=""),sep="/"),
            row.names=F,sep = "\t",quote=F)


# Sort by 
class <- class[order(class$V5,decreasing = T),]
class$V5 <- factor(class$V5, levels=unique(class$V5))

(plot = ggplot(class,aes(x=V5,y=n,label=n,fill=V5)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    labs(y = "Frequency",
         x = "CAZy Class(Level1)",
         fill = "Level1")+
    theme(legend.position = "none"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Level1_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=14, height=16, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Level1_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=35, units = "cm")

class <- class[order(class$V5,decreasing = T),1:2]
colnames(class) = c("CAZy Class","Freq")
write.table(class,file = paste(outDir,paste(prefix,"_Level1_Freq.xls",sep=""),sep="/"),
            row.names=F,sep = "\t",quote=F)

####################################################################################
# Sort by 
family <- family[order(family$n,decreasing = T),]
family <- family[1:50,]
family <- family[!is.na(family$V2),]
#family$V6 <- substr(family$V2,start = 0,stop = 50)

family <- family[order(family$V5,decreasing = T),]
family$V2 <- factor(family$V2, levels=unique(family$V2))

family = family[family$V5 != "-",]
family = family[family$V2 != "-",]

(plot = ggplot(family,aes(x=V2,y=n,fill=V5,label=n)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    labs(title = "Top 50",
         fill = "Level2",
         y = "Frequency",
         x = "CAZy Family(Level2)"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Level2_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=16, height=28, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Level2_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=35, height=55, units = "cm")

family <- family[order(family$V2,decreasing = T),1:3]
colnames(family) = c("CAZy Class","CAZy Family","Freq")
write.table(family,file = paste(outDir,paste(prefix,"_Level2_Freq.xls",sep=""),sep="/"),
            row.names=F,sep = "\t",quote=F)