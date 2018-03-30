# function
metapca <- function(
  m,
  main = NULL,
  group = NULL,
  group_colour = NULL,
  select = NULL,
  point_size  = 8,
  point_alpha = 1,
  point_text_size = 4
){
  pca <- vegan::rda(m)
  sum = summary(pca)
  sites_df <-  data.frame(sum$sites)
  sum_df <- data.frame(sum$cont$importance)
  # Eigenvalue
  (eigv_pca_1 <- as.numeric(sprintf("%.3f",sum_df[2,1]))*100)
  (eigv_pca_2 <- as.numeric(sprintf("%.3f",sum_df[2,2]))*100)
  name <- row.names(sites_df)
  sites_df$name = name
  if(!is.null(group)) {
    groups <- unlist(strsplit(group,","))
    sites_df$Groups = groups
  } else {
    groups = group
    sites_df$Groups = groups
  }
  xlab=paste("PC1(",eigv_pca_1,"%)",sep="")
  ylab=paste("PC2(",eigv_pca_2,"%)",sep="")
  
  ggpca = ggplot(sites_df,aes(PC1,PC2))
  if(!is.null(main)) {
    ggpca <- ggpca + labs(title = main)
    ggpca = ggpca + xlab(xlab)
    ggpca = ggpca + ylab(ylab)
    ggpca = ggpca +
      geom_hline(yintercept=0,linetype=2,color="grey",size=.5) + 
      geom_vline(xintercept=0,linetype=2,color="grey",size=.5)
  }
  if(!is.null(group)) {
    ggpca <- ggpca +
      geom_point(data  = sites_df,
                 size  =  point_size,
                 alpha =  point_alpha, 
                 aes(color=Groups)) + 
      geom_text(data  = sites_df,
                size  = point_text_size,
                aes(label = name),
                color = "grey56")
  }else{
    ggpca <- ggpca +
      geom_point(data  = sites_df,
                 size  = point_size,
                 alpha = point_alpha) + 
      geom_text(data  = sites_df,
                size  = point_text_size,
                aes(label = name),
                color = "grey56")
  }
  colours <- c("#F4BD6D","#87CCC5","olivedrab3","brown2","goldenrod1","deepskyblue3","maroon2","mediumpurple4",
               "cadetblue4","darkblue","hotpink3","indianred1","gray59",
               "darkorange1","slategray3","grey95","mediumorchid3","paleturquoise1",
               "orangered4","grey64","gray62","royalblue4","antiquewhite3")
  
  cols <- colours[1:dim(as.data.frame(table(sites_df$Groups)))[1]]
  ggpca <- ggpca + 
    scale_colour_manual(values = cols)
  if(!is.null(group_colour)) {
    col <- unlist(strsplit(group_colour,","))
    ggpca = ggpca + scale_colour_manual(values = col)
  }
  ggpca <- ggpca + theme_bw()
  ggpca <- ggpca + theme(
    legend.title = element_text(angle = 0,  
                                hjust =.95, 
                                vjust =.95,
                                size  = 12),
    legend.position   = "right",
    legend.background = element_rect(
      colour = "white", 
      fill = "white", 
      size = 1, 
      linetype = 1),#'dashed'
    legend.key.size   = unit(.5, "cm"),
    legend.key.height = unit(.5, "cm"),
    legend.key.width  = unit(.5, "cm"),
    legend.text       = element_text(angle = 0,  
                                     hjust =.95, 
                                     vjust =.95,
                                     size  = 10),
    legend.justification = "top",
    panel.background = element_rect(fill = "white",colour = "white",size = 1.5))
  
  return(list(plot = ggpca,
              sites= sites_df,
              sum  = sum_df))
}

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
require(pheatmap)
require(RColorBrewer)
require(ggplot2)
require(methods)
require(WriteXLS)

samples = gsub("-", ".", samples)
groups  = gsub("-", ".", groups)
smplst = strsplit(samples,split = ",")[[1]]
smplst.ord = smplst[order(smplst)]
grplst = strsplit(groups,split = ",")[[1]]
grplst.ord = grplst[order(smplst)]

annotation_row = data.frame(Group = grplst.ord)
rownames(annotation_row) = smplst.ord

# intable = "/Bio2/Project/SGM/SGM1052-Pig-Gut-24-metaGenome-10G.5GAssembly/BinGenome_Result/05.Taxonomy/RPKM/Order.xls"
table = t(read.table(intable,sep = "\t",header = T,row.names = 1,quote = ''))
table = table[order(rownames(table)),]

if (is.null(dim(table))){
  print("Do Nothing !! Because We Can Calculate the Div ~ ")
}else{
  eupc = metapca(table,main = "Euclidean",group = paste(grplst[order(smplst)],collapse = ","))
  write.table(eupc$sites,paste(outDir,paste(prefix,"_Euclidean_pca.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  ggsave(eupc$plot,
         filename = paste(outDir,paste(prefix,"_Euclidean_pca.png",sep = ""),sep="/"),
         device = "png",
         units = "cm",
         width = 20,
         height = 17)
  ggsave(eupc$plot,
         filename = paste(outDir,paste(prefix,"_Euclidean_pca.png",sep = ""),sep="/"),
         device = "png",
         width = 10,
         height = 8.5)
}