import pandas as pd
import sys
import re
import BASICGraph
import forSum

blastout = sys.argv[1]
smp = blastout.split('.')[0]

#### blast results

df = pd.read_table(blastout,header=None,dtype=str)
df[13] = df[12].apply(
    lambda x:re.findall('(?<= )(.*?)(?=OS)',x)[0] \
if re.findall('(?<= )(.*?)(?=OS)',x) else '-')
df[14] = df[12].apply(
    lambda x:re.findall('(?<=GN=)(.*?)(?= )',x)[0] \
if re.findall('(?<=GN=)(.*?)(?= )',x) else '-')
df[[0,1,2,3,10,11,14,13]].to_csv('%s.SwissProt.Anno.xls'%smp,index=None,sep='\t',quoting=3,
                              header=['query','subject','identity(%)','align_len',
                                      'Evalue','score','Gene_name','Protein_name'])
out = pd.read_table("%s.SwissProt.Anno.xls"%smp,sep="\t")
df = pd.DataFrame(out["Protein_name"].value_counts()[:20])
DF = pd.DataFrame({"Protein":df.index,"values":df.iloc[:,0]})
ax,fig,c = BASICGraph.BasicHBar(df=DF,colorsets="pink",Title="SwissProt Protein",xTitle="Freq",yTitle="",width=2.4)
fig.savefig(smp+".SwissProt_Top10.PNG",dpi=100)
fig.savefig(smp+".SwissProt_Top10.PDF")
