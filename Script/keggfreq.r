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

theme = theme(strip.background     = element_blank(),
              strip.text           = element_blank())

kegg_path = read.table(intable,
                       sep = "\t",header = T)
##########################################################################################################
df = data.frame(kegg_path)
by_A <- df[,c(1,5)] %>% group_by(KEGG_A_class)
dfA <- by_A %>% summarise(sum(Count))

dfA <- dfA[order(dfA$`sum(Count)`),]
dfA$KEGG_A_class <- factor(dfA$KEGG_A_class, levels=unique(dfA$KEGG_A_class))

(plot = ggplot(dfA,aes(x=KEGG_A_class,y=`sum(Count)`,fill=KEGG_A_class,label=`sum(Count)`)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    labs(y = "Frequency",
         x = "Level A",
         fill = "Level A"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_LevelA_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=35, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_LevelA_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=35, units = "cm")

##########################################################################################################
df = data.frame(kegg_path)
by_AB <- df[,c(1,2,5)] %>% group_by(KEGG_A_class,KEGG_B_class)
dfB <- by_AB %>% summarise(sum(Count))

dfB <- dfB[order(dfB$KEGG_A_class,dfB$`sum(Count)`),]
dfB$KEGG_B_class <- factor(dfB$KEGG_B_class, levels=unique(dfB$KEGG_B_class))

(plot = ggplot(dfB,aes(x=KEGG_B_class,y=`sum(Count)`,fill=KEGG_A_class,label=`sum(Count)`)) + 
  Kitult::ggBasicBar() +
  Kitult::ggBasicText(size = 2) + 
  Kitult::reportStyle() +
  coord_flip() + 
  #facet_grid(KEGG_A_class ~ .,scales="free_y", space="free_y",as.table=T) + 
  #theme +
  labs(y = "Frequency",
       x = "Level B",
       fill = "Level A"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_LevelB_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=35, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_LevelB_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=35, units = "cm")
#######################################################################################################
df = data.frame(kegg_path)
by_AC <- df[,c(1,3,5)] %>% group_by(KEGG_A_class,Pathway)
dfC <- by_AC %>% summarise(sum(Count))

dfC <- dfC[order(dfC$KEGG_A_class,dfC$`sum(Count)`),]
dfC$Pathway <- factor(dfC$Pathway, levels=unique(dfC$Pathway))
dfC <- dfC[order(dfC$`sum(Count)`,decreasing = T),]
dfC <- dfC[1:50,]

(plot = ggplot(dfC,aes(x=Pathway,y=`sum(Count)`,fill=KEGG_A_class,label=`sum(Count)`)) + 
    Kitult::ggBasicBar() +
    Kitult::ggBasicText(size = 2) + 
    Kitult::reportStyle() +
    coord_flip() +
    #theme + 
    #facet_grid(KEGG_A_class ~ .,scales="free_y", space="free_y",as.table=T) + 
    labs(y = "Frequency",
         x = "Pathway",
         title = "Top 50 Pathway",
         fill = "Level A"))

ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Pathway_FreqBar.png",sep=""),sep="/"),
       device = "png", dpi = 100, width=28, height=35, units = "cm")
ggsave(plot,
       filename=paste(outDir,paste(prefix,"_Pathway_FreqBar.tiff",sep=""),sep="/"),
       device = "tiff", width=28, height=35, units = "cm")