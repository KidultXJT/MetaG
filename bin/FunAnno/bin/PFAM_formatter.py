from Bio import SearchIO
import pandas as pd
import sys

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
        lst.append([line.id,domain.id,domain.accession,line[0].evalue,line[0].bitscore])

extr=pd.DataFrame(lst,columns=['#seq_name','domain','domain_acc','Evalue','bitscore'])
extr = extr[extr['Evalue']<float(evalue)]

ref = pd.read_table(mapf)
merg = pd.merge(extr,ref,left_on='domain',right_on='PFAM_Name',how='left')
out = merg[['#seq_name','domain','domain_acc','Evalue','bitscore','PFAM_desc']].drop_duplicates()


outname = '%s.PFAM.xls'%(hmmtsv.split('.')[0])
out.to_csv(outname,sep='\t',index=None,na_rep='-')

