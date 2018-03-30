## --------------------------- ##
##   Title: cor::Heatmap       ##
##  Author: Kidult             ##
##    Date: 2017/9/27          ##
## --------------------------- ##
corHeatmap <- function(
  X, # Table X::sp, with more var
  Y, # Table Y::env
  outDir = getwd(),   # outDir
  method = "pearson", # Correlation Method::pearson or kandall or spearman. Can Be a list::"pearson,kendall,spearman"
  brwCol = "RdYlGn",  # select From RColorBrewer::Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd 
                                                # Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr
                                                # RdBu RdGy RdYlBu RdYlGn Spectral
  width  = 30,
  height = 30,
  free   = FALSE
  ){
  # Description:
  #
  # cor compute the variance of x and the covariance or correlation of x 
  # and y if these are vectors. If x and y are matrices then the 
  # correlations between the columns of x and the columns of y are computed.
  # 
  # Package:
  require(pheatmap)
  require(RColorBrewer)

  methods = c(strsplit(method,split = ",")[[1]])
  for(i in methods){
    corTable = cor(X,Y,method = paste(i))
    write.table(corTable,paste(path.expand(outDir),"/cor.",i,".xls",sep = ""),sep = "\t",quote = F)
    pheatmap(corTable,
             filename=paste(path.expand(outDir),"/cor.",i,".heatmap.png",sep = ""),
             width=8,
             height=8,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
             )
    pheatmap(corTable,
             filename=paste(path.expand(outDir),"/cor.",i,".heatmap.pdf",sep = ""),
             width=width,
             height=height,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
             )
    
    if(free){
      print("We Don't Like the Free Project !!!")
    }else{
      require(corrgram)
      png(filename = paste(path.expand(outDir),"/cor.",i,".halfheatmap.png",sep = ""),
          width = width*20,
          height = height*20,
          units = "px")
      corrgram(corTable,
               main="Correlation Matrix",
               order=NULL, 
               lower.panel=panel.shade,
               upper.panel=NULL, 
               text.panel=panel.txt,
               col.regions = colorRampPalette(c(rev(brewer.pal(n = 7, name = brwCol)))),
               cor.method = paste(i)
               )
      dev.off()
    }
  }
  return(corTable)
}

corTESTHeatmap <- function(
  X, # Table X::sp, with more var
  Y, # Table Y::env
  outDir = getwd(),   # outDir
  method = "pearson", # Correlation Method::pearson or kandall or spearman. Can Be a list::"pearson,kendall,spearman"
  brwCol = "RdYlGn",  # select From RColorBrewer::Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd 
  # Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr
  # RdBu RdGy RdYlBu RdYlGn Spectral
  width  = 30,
  height = 30,
  free   = TRUE
){
  # Description:
  #
  # cor compute the variance of x and the covariance or correlation of x 
  # and y if these are vectors. If x and y are matrices then the 
  # correlations between the columns of x and the columns of y are computed.
  # 
  # Package:
  require(pheatmap)
  require(RColorBrewer)
  
  methods = c(strsplit(method,split = ",")[[1]])
  for(i in methods){
    rlst = c()
    plst = c()
    for(x in colnames(X)){
      assign(paste(x), X[,colnames(X) == x])
      for(y in colnames(Y)){
        assign(paste(y), Y[,colnames(Y) == y])
        ct <- cor.test(get(paste(x)),get(paste(y)),method = i)
        r          <- ct$estimate
        rlst       <- c(rlst,r)
        p          <- ct$p.value
        plst       <- c(plst,p)
      }
    }
    r.m <- matrix(rlst,ncol(Y)) 
    rownames(r.m) = colnames(Y)
    colnames(r.m) = colnames(X)
    
    p.m <- matrix(plst,ncol(Y))
    rownames(p.m) = colnames(Y)
    colnames(p.m) = colnames(X)
    
    pheatmap(r.m,
             filename=paste(path.expand(outDir),"/cor.",i,".heatmap.png",sep = ""),
             width=8,
             height=8,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
    )
    pheatmap(r.m,
             filename=paste(path.expand(outDir),"/cor.",i,".heatmap.pdf",sep = ""),
             width=width,
             height=height,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
    )
    pheatmap(p.m,
             filename=paste(path.expand(outDir),"/sig.",i,".heatmap.png",sep = ""),
             width=8,
             height=8,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
    )
    pheatmap(p.m,
             filename=paste(path.expand(outDir),"/sig.",i,".heatmap.pdf",sep = ""),
             width=width,
             height=height,
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white",
             cluster_rows = F,
             cluster_cols = F
    )
    
    write.table(r.m,paste(path.expand(outDir),"/cor.",i,".xls",sep = ""),sep = "\t",quote = F)
    write.table(p.m,paste(path.expand(outDir),"/sig.",i,".xls",sep = ""),sep = "\t",quote = F)
  }
  
  return(list(statistic=r.m,pvalue=p.m))
}
