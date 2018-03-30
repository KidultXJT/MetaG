from Bio.Blast import NCBIXML
import sys
import pandas as pd

blastxml = sys.argv[1]
evalue = sys.argv[2]
smp = blastxml.split('.')[0]

lst = []
for seq in NCBIXML.parse(open(blastxml)):
    if seq.alignments:
        slst = []
        slst.append(seq.query.split(' ')[0])
        slst.append(seq.descriptions[0].title.split(' ')[0])
        slst.append(seq.alignments[0].hsps[0].identities/float(seq.alignments[0].length))
        slst.append(seq.alignments[0].length)
        slst.append(seq.descriptions[0].e)
        slst.append(seq.descriptions[0].bits)
        slst += seq.descriptions[0].title.split(' ',2)[1:]
        lst.append(slst)

tbl = pd.DataFrame(lst,columns=['#query','target','identity','aln_len','Evalue','bitscore','protein_id','descriptions'])
fltr = tbl[tbl['Evalue']<float(evalue)]
outname = '%s.tbl'%smp
fltr.to_csv(outname,sep='\t',index=None,na_rep='-')

