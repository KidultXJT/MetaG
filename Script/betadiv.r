# function
metapcoa <- function(
  div,
  main = NULL,
  group = NULL,
  group_colour = NULL,
  select = NULL,
  point_size  = 8,
  point_alpha = 1,
  point_text_size = 4
  ){
  pcoa <- ape::pcoa(div)
  sites_df <-  data.frame(pcoa$vectors)
  sum_df <- data.frame(pcoa$values)
  # Eigenvalue
  (eigv_pcoa_1 <- as.numeric(sprintf("%.3f",pcoa$values$Rel_corr_eig[1]))*100)
  (eigv_pcoa_2 <- as.numeric(sprintf("%.3f",pcoa$values$Rel_corr_eig[2]))*100)
  name <- row.names(sites_df)
  sites_df$name = name
  if(!is.null(group)) {
    groups <- unlist(strsplit(group,","))
    sites_df$Groups = groups
  } else {
    groups = group
    sites_df$Groups = groups
  }
  xlab=paste("PC1(",eigv_pcoa_1,"%)",sep="")
  ylab=paste("PC2(",eigv_pcoa_2,"%)",sep="")
  
  ggpcoa = ggplot(sites_df,aes(Axis.1,Axis.2))
  if(!is.null(main)) {
    ggpcoa <- ggpcoa + labs(title = main)
    ggpcoa = ggpcoa +
      xlab(xlab)
    ggpcoa = ggpcoa +
      ylab(ylab)
    ggpcoa = ggpcoa +
      geom_hline(yintercept=0,linetype=2,color="grey",size=.5) + 
      geom_vline(xintercept=0,linetype=2,color="grey",size=.5)
  }
  if(!is.null(group)) {
    ggpcoa <- ggpcoa +
      geom_point(data  = sites_df,
                 size  =  point_size,
                 alpha =  point_alpha, 
                 aes(color=Groups)) + 
      geom_text(data  = sites_df,
                size  = point_text_size,
                aes(label = name),
                color = "grey56")
  }else{
    ggpcoa <- ggpcoa +
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
  ggpcoa <- ggpcoa + 
    scale_colour_manual(values = cols)
  if(!is.null(group_colour)) {
    col <- unlist(strsplit(group_colour,","))
    ggpcoa = ggpcoa + scale_colour_manual(values = col)
  }
  ggpcoa <- ggpcoa + theme_bw()
  ggpcoa <- ggpcoa + theme(
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
  
  return(list(plot = ggpcoa,
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

#intable = "/Bio/User/xiejunting/EzMetaG/EzMetaAs/BinGenome/Example/04.Taxonomy/AvgDepth/RelativeAbundance/Phylum_Dominant.xls"
table = t(read.table(intable,sep = "\t",header = T,row.names = 1,quote = ''))
table = table[order(rownames(table)),]

if (is.null(dim(table))){
  print("Do Nothing !! Because We Can Calculate the Div ~ ")
}else{
  
  bcdiv = as.matrix(vegan::vegdist(table,method = "bray",diag = T,upper = T))
  eudiv = as.matrix(vegan::vegdist(table,method = "euclidean",diag = T,upper = T))
  mhdiv = as.matrix(vegan::vegdist(table,method = "mahalanobis",diag = T,upper = T))
  jcdiv = as.matrix(vegan::vegdist(table,method = "jaccard",diag = T,upper = T))
  
  outTable = cbind(SampleID = rownames(bcdiv),data.frame(bcdiv))
  write.table(outTable,paste(outDir,paste(prefix,"_BrayCurtis_div.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = F)
  outTable = cbind(SampleID = rownames(eudiv),data.frame(eudiv))
  write.table(outTable,paste(outDir,paste(prefix,"_Euclidean_div.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = F)
  outTable = cbind(SampleID = rownames(mhdiv),data.frame(mhdiv))
  write.table(outTable,paste(outDir,paste(prefix,"_Mahalanobis_div.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = F)
  outTable = cbind(SampleID = rownames(jcdiv),data.frame(jcdiv))
  write.table(outTable,paste(outDir,paste(prefix,"_Jaccard_div.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = F)
  
  # hclust
  # tiff
  png(filename = paste(outDir,paste(prefix,"_BrayCurtis","_ClusterDendrogram.png",sep=""),sep="/"),
      width  = 850, 
      height = 400, 
      units = "px")
  plot(hclust(vegan::vegdist(table,method = "bray"),method = "average"),main = "Bray Curtis",sub = "",ylab = "", xlab = "")
  dev.off()
  png(filename = paste(outDir,paste(prefix,"_Euclidean","_ClusterDendrogram.png",sep=""),sep="/"),
      width  = 850, 
      height = 400, 
      units = "px")
  plot(hclust(vegan::vegdist(table,method = "euclidean"),method = "average"),main = "Euclidean",sub = "",ylab = "", xlab = "")
  dev.off()
  png(filename = paste(outDir,paste(prefix,"_Mahalanobis","_ClusterDendrogram.png",sep=""),sep="/"),
      width  = 850, 
      height = 400, 
      units = "px")
  plot(hclust(vegan::vegdist(table,method = "mahalanobis"),method = "average"),main = "Mahalanobis",sub = "",ylab = "", xlab = "")
  dev.off()
  png(filename = paste(outDir,paste(prefix,"_Jaccard","_ClusterDendrogram.png",sep=""),sep="/"),
      width  = 850, 
      height = 400, 
      units = "px")
  plot(hclust(vegan::vegdist(table,method = "jaccard"),method = "average"),main = "Jaccard",sub = "",ylab = "", xlab = "")
  dev.off()
  
  # tiff
  tiff(filename = paste(outDir,paste(prefix,"_BrayCurtis","_ClusterDendrogram.tiff",sep=""),sep="/"),
       width  = 850, 
       height = 400, 
       units = "px")
  plot(hclust(vegan::vegdist(table,method = "bray"),method = "average"),main = "Bray Curtis",sub = "",ylab = "", xlab = "")
  dev.off()
  tiff(filename = paste(outDir,paste(prefix,"_Euclidean","_ClusterDendrogram.tiff",sep=""),sep="/"),
       width  = 850, 
       height = 400, 
       units = "px")
  plot(hclust(vegan::vegdist(table,method = "euclidean"),method = "average"),main = "Euclidean",sub = "",ylab = "", xlab = "")
  dev.off()
  tiff(filename = paste(outDir,paste(prefix,"_Mahalanobis","_ClusterDendrogram.tiff",sep=""),sep="/"),
       width  = 850, 
       height = 400, 
       units = "px")
  plot(hclust(vegan::vegdist(table,method = "mahalanobis"),method = "average"),main = "Mahalanobis",sub = "",ylab = "", xlab = "")
  dev.off()
  tiff(filename = paste(outDir,paste(prefix,"_Jaccard","_ClusterDendrogram.tiff",sep=""),sep="/"),
       width  = 850, 
       height = 400, 
       units = "px")
  plot(hclust(vegan::vegdist(table,method = "jaccard"),method = "average"),main = "Jaccard",sub = "",ylab = "", xlab = "")
  dev.off()
  
  # PCoA
  bcpc = metapcoa(bcdiv,main = "Bray-Curtis",group = paste(grplst[order(smplst)],collapse = ","))
  #eupc = metapcoa(eudiv,main = "Bray-Curtis",group = paste(grplst[order(smplst)],collapse = ","))
  mhpc = metapcoa(mhdiv,main = "Mahalanobis",group = paste(grplst[order(smplst)],collapse = ","))
  jcpc = metapcoa(jcdiv,main = "Jaccard",group = paste(grplst[order(smplst)],collapse = ","))
  
  write.table(bcpc$sites,paste(outDir,paste(prefix,"_BrayCurtis_pcoa.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  #write.table(eupc$sites,paste(outDir,paste(prefix,"_Euclidean_pcoa.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  write.table(mhpc$sites,paste(outDir,paste(prefix,"_Mahalanobis_pcoa.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  write.table(jcpc$sites,paste(outDir,paste(prefix,"_Jaccard_pcoa.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  
  write.table(bcpc$sum,paste(outDir,paste(prefix,"_BrayCurtis_pcoaSUM.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  #write.table(eupc$sum,paste(outDir,paste(prefix,"_Euclidean_pcoaSUM.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  write.table(mhpc$sum,paste(outDir,paste(prefix,"_Mahalanobis_pcoaSUM.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  write.table(jcpc$sum,paste(outDir,paste(prefix,"_Jaccard_pcoaSUM.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE)
  
  ggsave(bcpc$plot,
         filename = paste(outDir,paste(prefix,"_BrayCurtis_pcoa.png",sep = ""),sep="/"),
         device = "png",
         units = "cm",
         width = 20,
         height = 17)
  #ggsave(eupc$plot,
  #       filename = paste(outDir,paste(prefix,"_Euclidean_pcoa.png",sep = ""),sep="/"),
  #       device = "png",
  #       units = "cm",
  #       width = 20,
  #       height = 17)
  ggsave(mhpc$plot,
         filename = paste(outDir,paste(prefix,"_Mahalanobis_pcoa.png",sep = ""),sep="/"),
         device = "png",
         units = "cm",
         width = 20,
         height = 17)
  ggsave(jcpc$plot,
         filename = paste(outDir,paste(prefix,"_Jaccard_pcoa.png",sep = ""),sep="/"),
         device = "png",
         units = "cm",
         width = 20,
         height = 17)
  ggsave(bcpc$plot,
         filename = paste(outDir,paste(prefix,"_BrayCurtis_pcoa.pdf",sep = ""),sep="/"),
         device = "pdf",
         width = 10,
         height = 8.5)
  #ggsave(eupc$plot,
  #       filename = paste(outDir,paste(prefix,"_Euclidean_pcoa.pdf",sep = ""),sep="/"),
  #       device = "pdf",
  #       width = 10,
  #       height = 8.5)
  ggsave(mhpc$plot,
         filename = paste(outDir,paste(prefix,"_Mahalanobis_pcoa.pdf",sep = ""),sep="/"),
         device = "pdf",
         width = 10,
         height = 8.5)
  ggsave(jcpc$plot,
         filename = paste(outDir,paste(prefix,"_Jaccard_pcoa.pdf",sep = ""),sep="/"),
         device = "pdf",
         width = 10,
         height = 8.5)
  
  
  
  # heatmap
  #png
  pheatmap(bcdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_BrayCurtis_div.png",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(eudiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           #display_numbers = T,
           #number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Euclidean_div.png",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(mhdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Mahalanobis_div.png",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(jcdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Jaccard_div.png",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  
  # PDF
  pheatmap(bcdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_BrayCurtis_div.pdf",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(eudiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           #display_numbers = T,
           #number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Euclidean_div.pdf",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(mhdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Mahalanobis_div.pdf",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(jcdiv,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(100),
           border_color = "white",
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_Jaccard_div.pdf",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  
  # Samples
  spcor = cor(t(table),method = "spearman")
  outTable = cbind(SampleID = rownames(spcor),data.frame(spcor))
  write.table(outTable,paste(outDir,paste(prefix,"_CorSpearman_div.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = FALSE)
  if (nrow(data.frame(table(spcor)))==1){
    pheatmap(spcor,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(2),
             border_color = "white",
             annotation_row = annotation_row,
             annotation_col = annotation_row,
             annotation_legend = F,
             display_numbers = T,
             number_color = "white",
             cluster_rows = F,
             cluster_cols = F,
             legend = F,
             filename = paste(outDir,paste(prefix,"_CorSpearman_div.png",sep = ""),sep="/"),
             width = 12,
             height = 12
    )
    pheatmap(spcor,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(2),
             border_color = "white",
             annotation_row = annotation_row,
             annotation_col = annotation_row,
             annotation_legend = F,
             display_numbers = T,
             number_color = "white",
             cluster_rows = F,
             cluster_cols = F,
             legend = F,
             filename = paste(outDir,paste(prefix,"_CorSpearman_div.pdf",sep = ""),sep="/"),
             width = 12,
             height = 12
    )
  }else {
    pheatmap(spcor,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100),
             border_color = "white",
             annotation_row = annotation_row,
             annotation_col = annotation_row,
             annotation_legend = F,
             display_numbers = T,
             number_color = "white",
             cluster_rows = F,
             cluster_cols = F,
             legend = F,
             filename = paste(outDir,paste(prefix,"_CorSpearman_div.png",sep = ""),sep="/"),
             width = 12,
             height = 12
    )
    pheatmap(spcor,
             color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100),
             border_color = "white",
             annotation_row = annotation_row,
             annotation_col = annotation_row,
             annotation_legend = F,
             display_numbers = T,
             number_color = "white",
             cluster_rows = F,
             cluster_cols = F,
             legend = F,
             filename = paste(outDir,paste(prefix,"_CorSpearman_div.pdf",sep = ""),sep="/"),
             width = 12,
             height = 12
    )
  }
  
  
  
  
  
  # vars
  spcor = cor(table,method = "spearman")
  outTable = cbind(SampleID = rownames(spcor),data.frame(spcor))
  write.table(outTable,paste(outDir,paste(prefix,"_CorSpearman.xls",sep = ""),sep="/"),sep = "\t",quote = FALSE,row.names = FALSE)
  pheatmap(spcor,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100),
           border_color = "white",
           #annotation_row = annotation_row,
           #annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_CorSpearman.png",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
  pheatmap(spcor,
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="Reds")))(100),
           border_color = "white",
           #annotation_row = annotation_row,
           #annotation_col = annotation_row,
           annotation_legend = F,
           display_numbers = T,
           number_color = "white",
           cluster_rows = F,
           cluster_cols = F,
           legend = F,
           filename = paste(outDir,paste(prefix,"_CorSpearman.pdf",sep = ""),sep="/"),
           width = 12,
           height = 12
  )
}


