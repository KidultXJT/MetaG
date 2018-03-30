import sys
from Bio import SeqIO

fasta = sys.argv[1]
smp = fasta.split('.')[0]
handle = SeqIO.parse(open(fasta),'fasta')
name = [i.name for i in handle]
out = open('%s.glst'%smp,'w')
out.write('\n'.join(name))
out.close()
