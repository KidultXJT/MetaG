import pandas as pd
import sys
import math
import os
import EzDraw
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

ref = pd.read_table(mapf)

out = pd.merge(df,ref,on='target',how='left')
out = out[out['resistance_type'].notnull()]

smp = blast6out.split('.')[0]
del out['mismatches'],out['gaps'],out['query_start'],out['query_end'],out['hit_start'],out['hit_end']
outname = '%s.ARDB.xls'%smp
out.to_csv(outname,sep='\t',index=None,na_rep='-')


tbl = out['resistance_class'].value_counts()[:10]#.drop("-")[:10]
fig = EzDraw.colCatBar(tbl,"ARDB resistance classification")
fig.savefig('%s.ARDB_Top10.pdf'%smp)
fig.savefig('%s.ARDB_Top10.png'%smp)
data = out["resistance_class"][out['resistance_class'] != "-"]
df,p,fig = forSum.freqHBar(data=data,colorsets="violet",Title="ARDB Annotation",xTitle="Frequency",yTitle="")
fig.savefig(smp+".ARDB.PNG",dpi=100)
fig.savefig(smp+".ARDB.PDF")
