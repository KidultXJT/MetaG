import sqlite3
import pandas as pd
import sys
import re

sys.path.append('/Bio/User/yinsijun/NeoFlow/pylib')
import EzDraw
import EzStat

handle =sys.argv[1]
db = sys.argv[2]
smp = handle.split('.')[0]

conn = sqlite3.connect(db)
db = conn.cursor()

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



def checkOG(row):
    try:
        tax = re.sub('\[[\d\w]+\]','',row[7])
        OG = re.findall('[^,]+(?=@%s)'%tax,row[8])[0]
        return OG
    except:
        return '-'


inptDF = pd.read_table(handle,comment='#',dtype=str,header=None)
inptDF['usedOG'] = inptDF.apply(checkOG,axis=1)
OGlst = inptDF['usedOG'].tolist()
OGs = ','.join(map(lambda x:'"%s"'%x,OGlst))
cmd='SELECT og.og, description, COG_categories FROM og WHERE og.og IN (%s)'%OGs
ogdic = {}
if db.execute(cmd):
    for og, desc, cat in db.fetchall():
        cat = ','.join(re.findall('(?<=u\')\w(?=\')',cat))
        ogdic[og] = [cat,desc]
inptDF[10] = inptDF['usedOG'].apply(
    lambda x:ogdic[x][0] if x in ogdic.keys() else '-')
inptDF[11] = inptDF['usedOG'].apply(
    lambda x:ogdic[x][1] if x in ogdic.keys() else '-')
inptDF[10] = inptDF[10].apply(
    lambda x:'['+x+']'+cogdict[x] if x in cogdict else '-')

outname = '%s.COG.Anno.xls'%smp
inptDF[[0,1,2,3,4,7,8,10,11]].to_csv(outname,sep='\t',index=None,
         header=['#queryname','seed_eggNOG_ortholog','seed_ortholog_evalue',
                  'seed_ortholog_score','predicted_gene_name','Annotation_tax_scope',
                 'OGs','COG cat','eggNOG annot'],na_rep='-')
inptDF[10] = inptDF[10].apply(lambda x:x[:50]+'...' if len(x)>50 else x)
figcog = EzDraw.colCatBar(EzStat.cmplxColPiv(inptDF[10],10),'COG categroy classification')
cogpdf = '%s.COG.pdf'%smp
figcog.savefig(cogpdf)
cogpng = '%s.COG.png'%smp
figcog.savefig(cogpng)



