## --------------------------- ##
##   Title: veganCCA::Points   ##
##  Author: Kidult             ##
##    Date: 2017/9/27          ##
## --------------------------- ##
sumCCA <- function(
  X,      # Table X::sp, with more var
  Y=NULL, # Table Y::env
  outDir=NULL
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
  
  cca = cca(X,Y)
  
  cca.sum = summary(cca)
  constraintTable = cca.sum$constraints
  spTable         = cca.sum$species
  siteTable       = cca.sum$sites
  
  if(is.null(Y)){
    sumTable = cca.sum$cont$importance
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable)
  }else{
    sumTable = cca.sum$concont$importance
    envTable = cca.sum$biplot
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/ccaY",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envTable,
                   WeightAvg = cca$CCA$wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/ccaSum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/ccaX",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/ccaSite",".xls",sep=""),sep="\t")
  }
  return(sumlist)
}

sumCCAenvfit <- function(
  X,     # Table X::sp, with more var
  Y,     # Table Y::env
  outDir=NULL
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
  
  cca = cca(X,Y)
  
  cca.sum = summary(cca)
  constraintTable = cca.sum$constraints
  spTable         = cca.sum$species
  siteTable       = cca.sum$sites
  
  if(is.null(Y)){
    sumTable = cca.sum$cont$importance
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable)
  }else{
    sumTable = cca.sum$concont$importance
    envTable = cca.sum$biplot
    envAll.lc <- envfit(cca, Y, display = "lc")
    envAll.wa <- envfit(cca, Y, display = "wa")
    
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/ccaY",".xls",sep=""),sep="\t")
      write.table(envAll.lc$vectors$arrows,paste(path.expand(outDir),"/ccaY.envfit.lc",".xls",sep=""),sep="\t")
      write.table(envAll.wa$vectors$arrows,paste(path.expand(outDir),"/ccaY.envfit.wa",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envAll.lc$vectors$arrows,
                   WeightAvg = cca$CCA$wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/ccaSum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/ccaX",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/ccaSite",".xls",sep=""),sep="\t")
  }
  return(sumlist)
}

sumRDA <- function(
  X,     # Table X::sp, with more var
  Y=NULL # Table Y::env
){
  # Description:
  #
  # Function cca performs correspondence analysis, or optionally constrained 
  # correspondence analysis (a.k.a. canonical correspondence analysis), or optionally
  # partial constrained correspondence analysis. Function rda performs redundancy
  # analysis, or optionally principal components analysis. These are all very popular
  # ordination techniques in community ecology.
  # 
  # Package:
  require(vegan)
  
  rda = rda(X,Y)
  return(summary(rda))
}

