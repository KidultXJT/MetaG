import pandas as pd
import sys

sys.path.append('/Bio/User/yinsijun/NeoFlow/pylib')
import EzDraw
import EzStat

infile = sys.argv[1]
smp = infile.split('.')[0]

cogdict = {'A':'RNA processing and modification',
       'B':'Chromatin Structure and dynamics',
       'C':'Energy production and conversion',
       'D':'Cell cycle control and mitosis',
       'E':'Amino Acid metabolis and transport',
       'F':'Nucleotide metabolism and transport',
       'G':'Carbohydrate metabolism and transport',
       'H':'Coenzyme metabolis',
       'I':'Lipid metabolism',
       'J':'Tranlsation',
       'K':'Transcription',
       'L':'Replication and repair',
       'M':'Cell wall/membrane/envelop biogenesis',
       'N':'Cell motility',
       'O':'Post-translational modification, protein turnover, chaperone functions',
       'P':'Inorganic ion transport and metabolism',
       'Q':'Secondary Structure',
       'T':'Signal Transduction',
       'U':'Intracellular trafficing and secretion',
       'V':'Defense mechanisms',
       'W':'Extracellular structures',
       'Y':'Nuclear structure',
       'Z':'Cytoskeleton',
       'R':'General Functional Prediction only',
       'S':'Function Unknown'}

f= pd.read_table(infile,comment='#',header=None,dtype=str)

f[10] = f[10].apply(lambda x:x if pd.isnull(x) else
                    " || ".join(map(lambda y:'['+y+']'+cogdict[y],
                                    filter(lambda z:z in cogdict.keys()
                                           ,x.split(',')))))

outname = '%s.COG.Anno.xls'%smp
f[[0,1,2,3,4,7,8,9,10,11]].to_csv(outname,sep='\t',index=None,
         header=['#queryname','seed_eggNOG_ortholog','seed_ortholog_evalue',
                  'seed_ortholog_score','predicted_gene_name','Annotation_tax_scope',
                 'OGs','bestOG|evalue|score','COG cat','eggNOG annot'],na_rep='-')

f[10] = f[10].apply(lambda x:x[:50]+'...' if len(x)>50 else x)
figcog = EzDraw.colCatBar(EzStat.cmplxColPiv(f[10],10),'COG categroy classification')
cogpdf = '%s.COG.pdf'%smp
figcog.savefig(cogpdf)
cogpng = '%s.COG.png'%smp
figcog.savefig(cogpng)



