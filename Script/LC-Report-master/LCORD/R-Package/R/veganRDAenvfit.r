## --------------------------- ##
##   Title: veganRDA::Points   ##
##  Author: Kidult             ##
##    Date: 2018/1/16          ##
## --------------------------- ##
veganRDAenvfit <- function(
  X, # Table X::sp, with more var
  Y=NULL, # Table Y::env
  outDir = getwd(),
  Group  = NULL # Like A,A,A,B,B,B
){
  # Description:
  #
  # Function cca performs correspondence analysis, or optionally 
  # constrained correspondence analysis (a.k.a. canonical 
  # correspondence analysis), or optionally partial constrained 
  # correspondence analysis. Function rda performs redundancy 
  # analysis, or optionally principal components analysis. 
  # These are all very popular ordination techniques in community ecology.
  # 
  # Package:
  require(vegan)
  require(ggplot2)
  
  if(!is.null(Y)){prefix = "RDA"}else{prefix = "PCA"}
  
  rda = rda(X,Y)
  
  sum = sumRDPenvfit(X,Y,outDir = outDir)
  sumTable  = sum$sumTable
  siteTable = sum$siteTable
  spTable   = sum$spTable
  
  if(!is.null(Y)){
    envTable = sum$envTable
    waTable  = sum$WeightAvg
  }
  
  
  # Sites DataFrame
  if(is.null(Group)){
    
    Col = "#009944"
    
    sidf = data.frame(X = siteTable[,1],
                      Y = siteTable[,2],
                      Name = rownames(siteTable),
                      Col = rep(Col,nrow(siteTable)))
    
  }else{
    
    Groups = c(strsplit(Group,split = ",")[[1]])
    
    
    
    if(length(Groups) == length(row.names(siteTable))){
      
      sidf = data.frame(X = siteTable[,1],
                        Y = siteTable[,2],
                        Name = rownames(siteTable),
                        Col = Groups)
      
    }else{
      
      Col = "#009944"
      
      sidf = data.frame(X = siteTable[,1],
                        Y = siteTable[,2],
                        Name = rownames(siteTable),
                        Col = rep(Col,nrow(siteTable)))
      
    }
  }
  
  # Species DataFrame
  spdf = data.frame(X = spTable[,1],
                    Y = spTable[,2],
                    Name = rownames(spTable),
                    Col  = rep("Species",nrow(spTable)))
  
  if (!is.null(Y)){
    axisName1 = "RDA1("
    axisName2 = "RDA2("
  }else{
    axisName1 = "PCA1("
    axisName2 = "PCA2("
  }
  
  XTitle = paste(axisName1,as.character(round(sumTable[2,1]*100,2)),"%"," ",as.character(round(sumTable[1,1],2)),")",sep = "")
  YTitle = paste(axisName2,as.character(round(sumTable[2,2]*100,2)),"%"," ",as.character(round(sumTable[1,2],2)),")",sep = "")
  
  
  # Basic with Sites
  p1 = ggplot(data = sidf, 
              aes(x = X,
                  y = Y,
                  color = Col,
                  label = Name)) +
    # Sites
    ggBasicPoint()+
    ggBasicText(Col = "grey50")+
    # Theme::reportStyle
    reportStyle() +
    geom_vline(xintercept=0,linetype=2,color="grey",size=.5)+
    geom_hline(yintercept=0,linetype=2,color="grey",size=.5)+
    labs(x=XTitle,y=YTitle)
  
  ggsave(filename = paste(path.expand(outDir),"/",prefix,".","site.pdf",sep = ""),
         p1,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/",prefix,".","site.png",sep = ""),
         p1,
         width  = 9,
         height = 8)
  
  # With Sites + Species
  p2 = p1 +
    # Sp
    geom_point(data = spdf,
               aes(x = spdf[,1],
                   y = spdf[,2]),
               size = 2,
               shape = 2,
               color = "#009944")
  
  
  p3 = p1 +
    # Sp
    geom_text(data = spdf,
              aes(x = spdf[,1],
                  y = spdf[,2],
                  color=spdf[,4]),
              size = 3,
              color = "#009944")
  ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","sp1.pdf",sep = ""),
         p2,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","sp1.png",sep = ""),
         p2,
         width  = 9,
         height = 8)
  ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","sp2.pdf",sep = ""),
         p3,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","sp2.png",sep = ""),
         p3,
         width  = 9,
         height = 8)
  
  
  if(!is.null(Y)){
    # Environment DataFrame
    envdf = data.frame(X = envTable[,1],
                       Y = envTable[,2],
                       Name = rownames(envTable),
                       Col  = rep("Env",nrow(envTable)))
    
    # With Sites + Species + Env
    p4 = p2 + 
      geom_segment(
        data  = envdf,
        aes(x = 0, y = 0, xend = X*3.5, yend =Y*3.5),
        arrow = arrow(length = unit(0.01, "npc")),
        size  = .2,
        color = "grey10") + 
      geom_text(
        data  = envdf,
        aes(x = X*3.8, y = Y*3.8, label = Name),
        color = "grey10"
      )
    
    #ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","env.pdf",sep = ""),
    #       p4,
    #       width  = 18,
    #       height = 16)
    #ggsave(filename = paste(path.expand(outDir),"/",prefix,"envfit.","env.png",sep = ""),
    #       p4,
    #       width  = 9,
    #       height = 8)
    png(filename = paste(path.expand(outDir),"/",prefix,".","FLTenvfitLC.png",sep = ""), # filtered by Maximum estimated P value
        width  = 480, 
        height = 480, 
        units = "px", 
        pointsize = 12,
        bg = "white")
    plot(sum$rda,display = "lc")
    plot(sum$LC,col = "grey")
    plot(sum$LC,col = "red",p.max=pthreshold)
    dev.off()
    pdf(file = paste(path.expand(outDir),"/",prefix,".","FLTenvfitLC.pdf",sep = ""),
        width  = 6.4,
        height = 6.4)
    plot(sum$rda,display = "lc")
    plot(sum$LC,col = "grey")
    plot(sum$LC,col = "red",p.max=pthreshold)
    dev.off()
    png(filename = paste(path.expand(outDir),"/",prefix,".","FLTenvfitWA.png",sep = ""), # filtered by Maximum estimated P value
        width  = 480, 
        height = 480, 
        units = "px", 
        pointsize = 12,
        bg = "white")
    plot(sum$rda,display = "wa")
    plot(sum$LC,col = "grey")
    plot(sum$LC,col = "red",p.max=pthreshold)
    dev.off()
    pdf(file = paste(path.expand(outDir),"/",prefix,".","FLTenvfitWA.pdf",sep = ""),
        width  = 6.4,
        height = 6.4)
    plot(sum$rda,display = "wa")
    plot(sum$LC,col = "grey")
    plot(sum$LC,col = "red",p.max=pthreshold)
    dev.off()
  }
  return(sum)
}
