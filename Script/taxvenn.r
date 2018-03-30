Args <- commandArgs()

# Necessary
intable      <- Args[6]
samples      <- Args[7] # samples = "hsm2-1,ksm2-1,msm2-1,ssm2-1,hsm2-2,ksm2-2,msm2-2,ssm2-2,hsm2-3,ksm2-3,msm2-3,ssm2-3,hsm3-1,ksm3-1,msm3-1,ssm3-1,hsm3-2,ksm3-2,msm3-2,ssm3-2,hsm3-3,ksm3-3,msm3-3,ssm3-3"
groups       <- Args[8] # groups  = "hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm,hsm,ksm,msm,ssm"
compare      <- Args[9] # compare = "ksm:hsm:msm:ssm"
outDir       <- Args[10]
prefix        <- Args[11] # Taxonomy

library(Kitult)
library(lance)
require(pheatmap)
require(ggplot2)
require(methods)
require(WriteXLS)

#intable = "/Bio2/Project/SGM/SGM1059-Environment-Soil_Feces-12-metaGenome-10G/Taxonomy_Result2/TaxAnalysis/Reads/Phylum.xls"
table = read.table(intable,sep = "\t",header = T,row.names = 1,quote = '')
samples = gsub("-", ".", samples)
groups  = gsub("-", ".", groups)
compare = gsub("-", ".", compare)
smplst = strsplit(samples,split = ",")[[1]]
grplst = strsplit(groups,split = ",")[[1]]
cmplst = strsplit(compare,split = ",")[[1]]

# MultiVenn(<=5)
for(c in cmplst){
  Clst = strsplit(c,split = ":")[[1]]
  gTABLE = table[,order(smplst)]
  gTABLE = gTABLE[,smplst[order(smplst)][grplst[order(smplst)] %in% Clst]]
  grplst.ord = paste(grplst[order(smplst)][grplst[order(smplst)] %in% Clst],collapse = ",")
  gv = lance::lc.vennFromAnnot(df=gTABLE,groups=grplst.ord)
  vstring =  gsub(":", "-V.S-", cmplst)
  png(filename = paste(outDir,paste(prefix,"_",vstring,"_Venn.png",sep=""),sep="/"),
      width = 450, height = 450, units = "px")
  grid.draw(gv$plot)
  dev.off()
  tiff(filename = paste(outDir,paste(prefix,"_",vstring,"_Venn.tiff",sep=""),sep="/"),
       width = 450, height = 450, units = "px")
  grid.draw(gv$plot)
  dev.off()
  # save Table
  sets = gv$set
  names = as.character()
  for (i in seq(length(sets))){
    name = names(sets)[[i]]
    data = sets[[i]]$data
    data = cbind(rownames(data), data)
    colnames(data)[1] = prefix
    cmd = paste0("`",name,"`","=data")
    eval(parse(text = cmd))
    names = c(names, name)
  }
  WriteXLS(names, ExcelFileName = paste(outDir,paste(prefix,"_",vstring,"_Venn.xls",sep=""),sep="/"), SheetNames = names)
}


## Two Groups
#for(i in unique(grplst[order(smplst)])){
#  print(i)
#  for(j in unique(grplst[order(smplst)])){
#    if(j != i){
#      print(j)
#      print(paste(i,j,sep = " - "))
#      g2TABLE = gTABLE[,grplst[order(smplst)] %in% c(i,j)]
#      grp2 = grplst[order(smplst)][grplst[order(smplst)] %in% c(i,j)]
#      gv = lance::lc.vennFromAnnot(df=g2TABLE,groups=grp2,label.dist = -15)
#      vstring = paste(i,j,sep = "-V.S-")
#      print(vstring)
#      png(filename = paste(outDir,paste(vstring,"_2Groups_Venn.png",sep=""),sep="/"),
#          width = 450, height = 450, units = "px")
#      grid.draw(gv$plot)
#      dev.off()
#      tiff(filename = paste(outDir,paste(vstring,"_2Groups_Venn.tiff",sep=""),sep="/"),
#           width = 450, height = 450, units = "px")
#      grid.draw(gv$plot)
#      dev.off()
#      # save Table
#      sets = gv$set
#      names = as.character()
#      for (i in seq(length(sets))){
#        name = names(sets)[[i]]
#        data = sets[[i]]$data
#        data = cbind(rownames(data), data)
#        colnames(data)[1] = prefix
#        cmd = paste0("`",name,"`","=data")
#        eval(parse(text = cmd))
#        names = c(names, name)
#      }
#      WriteXLS(names, ExcelFileName = paste(outDir,paste(vstring,"_2Groups_Venn.xls",sep=""),sep="/"), SheetNames = names)
#    }
#  }
#}
#
### Two Samples
#for(i in smplst[order(smplst)]){
#  for(j in smplst[order(smplst)]){
#    if(j != i){
#      print(paste(i,j,sep = " - "))
#      s2TABLE = gTABLE[,c(i,j)]
#      sv = lance::lc.vennFromAnnot(df=s2TABLE,groups=paste(i,j,sep = ","),label.dist = -15)
#      vstring = paste(i,j,sep = "-V.S-")
#      png(filename = paste(outDir,paste(vstring,"_2Samples_Venn.png",sep=""),sep="/"),
#          width = 450, height = 450, units = "px")
#      grid.draw(sv$plot)
#      dev.off()
#      tiff(filename = paste(outDir,paste(vstring,"_2SampleGroups_Venn.tiff",sep=""),sep="/"),
#           width = 450, height = 450, units = "px")
#      grid.draw(sv$plot)
#      dev.off()
#      # save Table
#      sets = sv$set
#      names = as.character()
#      for (i in seq(length(sets))){
#        name = names(sets)[[i]]
#        data = sets[[i]]$data
#        data = cbind(rownames(data), data)
#        colnames(data)[1] = prefix
#        cmd = paste0("`",name,"`","=data")
#        eval(parse(text = cmd))
#        names = c(names, name)
#      }
#      WriteXLS(names, ExcelFileName = paste(outDir,paste(vstring,"_2Groups_Venn.xls",sep=""),sep="/"), SheetNames = names)
#    }
#  }
#}
#
## All Samples Core-Pan Graph

