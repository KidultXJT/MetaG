## --------------------------- ##
##   Title: veganCCA::Points   ##
##  Author: Kidult             ##
##    Date: 2017/9/27          ##
## --------------------------- ##
veganCCAenvfit <- function(
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
  
  cca = cca(X,Y)
  
  sum = sumCCAenvfit(X,Y,outDir = outDir)
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
  
  
  XTitle = paste("CCA1(",as.character(round(sumTable[2,1]*100,2)),"%"," ",as.character(round(sumTable[1,1],2)),")",sep = "")
  YTitle = paste("CCA2(",as.character(round(sumTable[2,2]*100,2)),"%"," ",as.character(round(sumTable[1,2],2)),")",sep = "")
  
  
  # Basic with Sites
  p1 = ggplot(data = sidf, 
              aes(x = X,
                  y = Y,
                  color = Col,
                  label = Name)) +
    # Sites
    BasicPoint()+
    BasicText(Col = "grey50")+
    # Theme::reportStyle
    reportStyle() +
    geom_vline(xintercept=0,linetype=2,color="grey",size=.5)+
    geom_hline(yintercept=0,linetype=2,color="grey",size=.5)+
    labs(x=XTitle,y=YTitle)
  
  ggsave(filename = paste(path.expand(outDir),"/CCA.","site.pdf",sep = ""),
         p1,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/CCA.","site.png",sep = ""),
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
  ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","sp1.pdf",sep = ""),
         p2,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","sp1.png",sep = ""),
         p2,
         width  = 9,
         height = 8)
  ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","sp2.pdf",sep = ""),
         p3,
         width  = 18,
         height = 16)
  ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","sp2.png",sep = ""),
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
        aes(x = 0, y = 0, xend = X*2, yend =Y*2),
        arrow = arrow(length = unit(0.01, "npc")),
        size  = .2,
        color = "grey10") + 
      geom_text(
        data  = envdf,
        aes(x = X*2.2, y = Y*2.2, label = Name),
        color = "grey10"
      )
    
    ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","env.pdf",sep = ""),
           p4,
           width  = 18,
           height = 16)
    ggsave(filename = paste(path.expand(outDir),"/CCAenvfit.","env.png",sep = ""),
           p4,
           width  = 9,
           height = 8)
  }
  
  #df = rbind(sidf,spdf,envdf)
  
  
  return(sum)
  
}
