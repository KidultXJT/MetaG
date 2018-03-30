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
    prefix = "ca"
  }else{
    sumTable = cca.sum$concont$importance
    envTable = cca.sum$biplot
    prefix = "cca"
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/",prefix,"Y",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envTable,
                   cca       = cca,
                   WeightAvg = cca$CCA$wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/",prefix,"Sum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/",prefix,"X",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/",prefix,"Site",".xls",sep=""),sep="\t")
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
    prefix = "ca"
  }else{
    sumTable = cca.sum$concont$importance
    envTable = cca.sum$biplot
    prefix = "cca"
    envAll.lc <- envfit(cca, Y, display = "lc")
    envAll.wa <- envfit(cca, Y, display = "wa")
    
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/",prefix,"Y",".xls",sep=""),sep="\t")
      write.table(envAll.lc$vectors$arrows,paste(path.expand(outDir),"/",prefix,"Y.envfit.lc",".xls",sep=""),sep="\t")
      write.table(envAll.wa$vectors$arrows,paste(path.expand(outDir),"/",prefix,"Y.envfit.wa",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envAll.lc$vectors$arrows,
                   WeightAvg = cca$CCA$wa,
                   cca       = cca,
                   LC        = envAll.lc,
                   WA        = envAll.wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/",prefix,"Sum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/",prefix,"X",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/",prefix,"Site",".xls",sep=""),sep="\t")
  }
  return(sumlist)
}

sumRDA <- function(
  X,     # Table X::sp, with more var
  Y=NULL,# Table Y::env
  outDir=NULL
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
  
  rda.sum = summary(rda)
  constraintTable = rda.sum$constraints
  spTable         = rda.sum$species
  siteTable       = rda.sum$sites
  
  if(is.null(Y)){
    sumTable = rda.sum$cont$importance
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable)
    prefix = "pca"
  }else{
    sumTable = rda.sum$concont$importance
    envTable = rda.sum$biplot
    prefix = "rda"
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/",prefix,"Y",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envTable,
                   rda       = rda,
                   WeightAvg = rda$CCA$wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/",prefix,"Sum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/",prefix,"X",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/",prefix,"Site",".xls",sep=""),sep="\t")
  }
  return(sumlist)
}

sumRDPenvfit <- function(
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
  
  rda = rda(X,Y)
  
  rda.sum = summary(rda)
  constraintTable = rda.sum$constraints
  spTable         = rda.sum$species
  siteTable       = rda.sum$sites
  
  if(is.null(Y)){
    sumTable = rda.sum$cont$importance
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable)
    prefix = "pca"
  }else{
    sumTable = rda.sum$concont$importance
    envTable = rda.sum$biplot
    prefix = "rda"
    envAll.lc <- envfit(rda, Y, display = "lc")
    envAll.wa <- envfit(rda, Y, display = "wa")
    
    if(!is.null(outDir)){
      write.table(envTable,paste(path.expand(outDir),"/",prefix,"Y",".xls",sep=""),sep="\t")
      write.table(envAll.lc$vectors$arrows,paste(path.expand(outDir),"/",prefix,"Y.envfit.lc",".xls",sep=""),sep="\t")
      write.table(envAll.wa$vectors$arrows,paste(path.expand(outDir),"/",prefix,"Y.envfit.wa",".xls",sep=""),sep="\t")
    }
    sumlist = list(sumTable  = sumTable,
                   siteTable = siteTable,
                   spTable   = spTable,
                   envTable  = envAll.lc$vectors$arrows,
                   WeightAvg = rda$CCA$wa,
                   rda       = rda,
                   LC        = envAll.lc,
                   WA        = envAll.wa)
  }
  if(!is.null(outDir)){
    write.table(sumTable,paste(path.expand(outDir),"/",prefix,"Sum",".xls",sep=""),sep="\t")
    write.table(spTable,paste(path.expand(outDir),"/",prefix,"X",".xls",sep=""),sep="\t")
    write.table(siteTable,paste(path.expand(outDir),"/",prefix,"Site",".xls",sep=""),sep="\t")
  }
  return(sumlist)
}