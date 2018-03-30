import pandas as pd
import sys
import os
# personal library
import forSum
import TableDealer

tax     = sys.argv[1] # Taxonomy From Kraken
covs    = sys.argv[2] # A1.cov.txt,A2.cov.txt,B.cov.txt,C.cov.txt...
delimit = sys.argv[3] # when k__|p__|c__|o__|f__|g__|s__ , then delimit == "|" 
smps    = sys.argv[4] # A1,A2,B,C
grps    = sys.argv[5] # A,A,B,C
values  = sys.argv[6] # Values Column
key     = sys.argv[7] # Key Column # NOT YET !!
header  = sys.argv[8] # No Header(0) | HAVE Header(1) # NOT YET !!
outDir  = sys.argv[9]
prefix  = sys.argv[10]

All,Set,lst = TableDealer.TableCombine(Anno=tax,vfilelist=covs,ColumNames=smps,Values = int(values))
Set.index = Set["Description"]
TABLE = Set.iloc[:,:-1]
G = pd.DataFrame({"Groups":grps.split(","),"ID":smps.split(",")}).T
G.columns = smps.split(",")
OUT = pd.concat([G,TABLE])
OUT.to_csv(outDir+"/"+prefix+".LEfSe.txt",sep='\t', encoding='utf-8',header=None,index=True)


