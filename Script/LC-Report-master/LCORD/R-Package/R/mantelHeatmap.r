## ------------------------- ##
##   Title: mantel::Heatmap  ##
##  Author: Kidult           ##
##    Date: 2017/9/28        ##
## ------------------------- ##
mantelHeatmap <- function(
  X, # Table X::sp, with more var
  Y, # Table Y::env
  Z=NULL, # Table Z
  outDir = getwd(),    # outDir
  dist   = "euclidean",# distance Method ::  "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", 
                       # "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", 
                       # "mahalanobis"
                       # SP  :: bray
                       # Env :: euclidean
  method = "pearson",  # Correlation Method::pearson or kandall or spearman. Can Be a list::"pearson,kendall,spearman"
  brwCol = "RdYlGn",   # select From RColorBrewer::Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd 
                                                 # Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr
                                                 # RdBu RdGy RdYlBu RdYlGn Spectral
  width  = 30,
  height = 30
  ){
  # Description:
  #
  # Before You Calculate the Mantel Test Correlation stat :: r, You should Distance Your Table X and Table Y First.
  # 
  # Distance(Dissimilarity Indices)::vegdist
  # The function computes dissimilarity indices that are useful for or popular with COMMUNITY ECOLOGISTS. All indices
  # use quantitative data, although they would be named by the corresponding binary index, but you can calculate the 
  # binary index using an appropriate argument.
  # 
  # Mantel Test
  # 
  # Package:
  require(vegan)
  require(pheatmap)
  require(RColorBrewer)
  
  # Distance
  # mantel
  # mantel partial
  # function mantel.partial finds the partial Mantel statistic as the partial matrix correlation between THREE(X,Y,Z) 
  # dissimilarity matrices. The significance of the statistic is evaluated by permuting rows and columns of the first
  # dissimilarity matrix.
  X.dist = vegdist(X,method=dist)
  Y.dist = vegdist(Y,method=dist)
  if(is.null(Z)){
    Z = NULL
    methods = c(strsplit(method,split = ",")[[1]])
    for(i in methods){
      mantel.XY = mantel(X.dist,Y.dist, method = i)
      mantel.XY.df = data.frame(dist.method = dist,
                                cor.method  = i,
                                stat.r = mantel.XY$statistic,
                                sig.p  = mantel.XY$signif)
      write.table(mantel.XY.df,paste(path.expand(outDir),"/mantel.",i,".",dist,".xls",sep = ""),sep = "\t",quote = F)
    }
  }else{
    Z.dist = vegdist(Z,method=dist)
    methods = c(strsplit(method,split = ",")[[1]])
    for(i in methods){
      mantel.XYZ = mantel.partial(X.dist,Y.dist,Z.dist,method = i)
      mantel.XYZ.df = data.frame(dist.method = dist,
                                 cor.method  = i,
                                 stat.r = mantel.XY$statistic,
                                 sig.p  = mantel.XY$signif)
      write.table(mantel.XYZ.df,paste(path.expand(outDir),"/mantelpartial.",i,".xls",sep = ""),sep = "\t",quote = F)
    }
  }
  
  #methods = c(strsplit(method,split = ",")[[1]])
  #for(i in methods){
  #mantel.XY = mantel(X.dist,Y.dist, method = i)
  #mantel.XY.df = data.frame(dist.method = dist,
  #                          cor.method  = i,
  #                          stat.r = mantel.XY$statistic,
  #                          sig.p  = mantel.XY$signif)
  #write.table(mantel.XY.df,paste(path.expand(outDir),"/mantel.",i,".xls",sep = ""),sep = "\t")
  #}
  
  methods = c(strsplit(method,split = ",")[[1]])
  for(i in methods){
    rlst = c()
    plst = c()
    for(x in colnames(X)){
      assign(paste("dist",x,sep = "."), dist(X[,colnames(X) == x],
                                             method = dist))
      for(y in colnames(Y)){
        assign(paste("dist",y,sep = "."), dist(Y[,colnames(Y) == y],
                                               method = dist))
        m <- mantel(get(paste("dist",x,sep = ".")),
                    get(paste("dist",y,sep = ".")))
        r          <- m$statistic
        rlst       <- c(rlst,r)
        p          <- m$signif
        plst       <- c(plst,p)
      }
    }
    r.m <- matrix(rlst,ncol(Y)) 
    rownames(r.m) = colnames(Y)
    colnames(r.m) = colnames(X)
    
    p.m <- matrix(plst,ncol(Y))
    rownames(p.m) = colnames(Y)
    colnames(p.m) = colnames(X)
    
    
    write.table(r.m,paste(path.expand(outDir),"/mantel.cor.",i,".",dist,".xls",sep = ""),sep = "\t",quote = F)
    write.table(p.m,paste(path.expand(outDir),"/mantel.sig.",i,".",dist,".xls",sep = ""),sep = "\t",quote = F)
    
    pheatmap(r.m,
             filename=paste(path.expand(outDir),"/mentel.cor.",i,".",dist,".heatmap.png",sep = ""),
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white", # Kidult Style
             width=width, 
             height=height
    )
    pheatmap(p.m,
             filename=paste(path.expand(outDir),"/mantel.sig.",i,".",dist,".heatmap.png",sep = ""),
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white", # Kidult Style
             width=width, 
             height=height
    )
    pheatmap(r.m,
             filename=paste(path.expand(outDir),"/mantel.cor.",i,".",dist,".heatmap.pdf",sep = ""),
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white", # Kidult Style
             width=width, 
             height=height
    )
    pheatmap(p.m,
             filename=paste(path.expand(outDir),"/mantel.sig.",i,".",dist,".heatmap.pdf",sep = ""),
             silent=FALSE,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
             border_color = "white", # Kidult Style
             width=width, 
             height=height
    )
  }
  return(list(Sig  = p.m,
              Stat = r.m,
              mantel = mantel.XY.df))
}
