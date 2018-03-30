import pandas as pd
import sys
import EzDraw

blast6out = sys.argv[1]
evalue = float(sys.argv[2])
identity = float(sys.argv[3])
smp = blast6out.split('.')[0]

anno=pd.read_table("/Bio/User/yinsijun/database/FuncAnno/PHI/current/for_pipe/phi_anno.xls")

df = pd.read_table(blast6out,
              names=['#query','target','identity','aln_len','mismatches','gaps','query_start','query_end','hit_start','hit_end','Evalue','bitscore'])
df = df.sort_values('#query')
df = df[(df['Evalue']<evalue)&(df['identity']>identity)]
df['Phibase_acc']=df['target'].apply(lambda x:x.split('#')[1])

ndf = pd.merge(df,anno,left_on="Phibase_acc",right_on="PHI_MolConn_ID",how="left")

outname = '%s.PHI.xls'%smp
del ndf['target']
del ndf['aln_len']
del ndf['mismatches']
del ndf['gaps']
del ndf['query_start']
del ndf['query_end']
del ndf['hit_start']
del ndf['hit_end']
del ndf['PHI_MolConn_ID']
#df[['#query','target','identity','aln_len','Evalue','bitscore','protein_id','Phibase_acc','gene','pathogen','pathogen_species','mutant_phenotype']].to_csv(outname,sep='\t',index=None,na_rep='-')
ndf.to_csv(outname,sep='\t',index=None,na_rep='-')
tbl = ndf['Phenotype_of_mutant'].value_counts()[:10]
fig = EzDraw.colCatBar(tbl,"PHI phenotype classification")
fig.savefig('%s.PHI.pdf'%smp)
fig.savefig('%s.PHI.png'%smp)
