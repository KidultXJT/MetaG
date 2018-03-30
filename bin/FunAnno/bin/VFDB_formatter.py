import pandas as pd
import sys
import math
import os
import EzDraw
import EzStat
import forSum

blast6out = sys.argv[1]
mapf = sys.argv[2]
evalue = sys.argv[3]

#try:
#    if 'e' in evalue:
#        base,index = map(int,evalue.split('e'))
#        threshold = base*(10**index)
#    else:
#        threshold = evalue
#except:
#    raise Exception('Wrong evalue format! use 0.05 or 1e-5')

df = pd.read_table(blast6out,
              names=['#query','target','identity','aln_len','mismatches','gaps','query_start','query_end','hit_start','hit_end','Evalue','bitscore'])
df = df[df['Evalue'] < float(evalue)]
df['VFGID']=df['target'].apply(lambda x:x.split('(')[0])


ref = pd.read_table(mapf)
out = pd.merge(df,ref,on='VFGID',how='left')
out = out[out['VFID'].notnull()]

smp = blast6out.split('.')[0]
outname = '%s.VFDB.xls'%smp
del out['mismatches'],out['gaps'],out['query_start'],out['query_end'],out['hit_start'],out['hit_end']
out.to_csv(outname,sep='\t',index=None,na_rep='-')

pp = EzStat.cmplxCol2line(out[out['Keyword'].notnull()],'Keyword',sep=';')
tbl = pp['Keyword'].value_counts().drop("-")[:10]
fig = EzDraw.colCatBar(tbl,"VFDB classification")
fig.savefig('%s.VFDB_Top10.pdf'%smp)
fig.savefig('%s.VFDB_Top10.png'%smp)
data = pp["Keyword"][pp['Keyword'] != "-"]
df,p,fig = forSum.freqHBar(data=data,colorsets="cornflowerblue",Title="VFDB(Set A) Annotation",xTitle="Frequency",yTitle="")
fig.savefig(smp+".VFDB.PNG",dpi=100)
fig.savefig(smp+".VFDB.PDF")

