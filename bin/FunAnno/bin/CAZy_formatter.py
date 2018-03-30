from Bio import SearchIO
import pandas as pd
import sys
import EzDraw
import forSum

hmmtsv = sys.argv[1]
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

tsv = SearchIO.parse(open(hmmtsv),'hmmer3-tab')
lst = []
for domain in tsv:
    for line in domain:
        lst.append([line.id,domain.id,line[0].evalue,line[0].bitscore])

extr=pd.DataFrame(lst,columns=['#seq_name','domain','Evalue','bitscore'])
extr['domain'] = extr['domain'].apply(lambda x:x.replace('.hmm',''))
extr = extr[extr['Evalue']<float(evalue)]


ref = pd.read_table(mapf)
merg = pd.merge(extr,ref,left_on='domain',right_on='#Family',how='left')
out = merg[['#seq_name','domain','Evalue','bitscore','cazy-class','cazy-note','cazy-activities']].drop_duplicates()
out = out.sort_values("#seq_name")
dic = {"GH":"Glycoside Hydrolases (GHs)",
       "GT":"GlycosylTransferases (GTs)",
       "PL":"Polysaccharide Lyases (PLs)",
       "CE":"Carbohydrate Esterases (CEs)",
       "AA":"Auxiliary Activities (AAs)",
       "CBM":"Carbohydrate-Binding Modules (CBMs)"}
out['cazy-class'] = out['cazy-class'].apply(
    lambda x:dic[x] if x in dic.keys() else x)

smp = hmmtsv.split('.')[0]
outname = '%s.CAZy.xls'%smp
out.to_csv(outname,sep='\t',index=None,na_rep='-')

tbl = out['cazy-class'].value_counts()[:10]
fig = EzDraw.colCatBar(tbl,"CAZy classification")
fig.savefig('%s.CAZy_Top10.pdf'%smp)
fig.savefig('%s.CAZy_Top10.png'%smp)

out = pd.read_table('%s.CAZy.xls'%smp,sep="\t")
ser= out["cazy-class"]
ser.index = out["cazy-class"]
ser = ser.drop('-')  
df,p,fig = forSum.freqHBar(data=ser,colorsets="deepskyblue",Title="CAZy Annotation",xTitle="Frequency",yTitle="",width=1.8)
fig.savefig(smp+".CAZy.PNG",dpi=100)
fig.savefig(smp+".CAZy.PDF")
