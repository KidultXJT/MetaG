import json
import networkx as nx
import sys
import re
import os
from PIL import Image,ImageDraw
import pandas as pd
from jinja2 import Template

blastout = sys.argv[1]
dbpath = sys.argv[2]

smp = blastout.split('.')[0]
PATH = "%s/brite"%dbpath
MAP = "%s/ko_desc.map"%dbpath

##### blastout manipulation

tbl = pd.read_table(blastout,header=None)
tbl[13] = tbl[12].apply(lambda x:x.split(' ',2)[1])
tbl[14] = tbl[12].apply(lambda x:x.split(' ',2)[2])
orglst = set(tbl[1].apply(lambda x:x.split(':')[0]).tolist())
tbl[[0,1,2,3,10,11,13,14]].to_csv('%s.KEGG.Anno.xls'%smp,index=None,sep='\t',
                                  header=['query','subject','identity(%)','align_len',
                                          'Evalue','score','K-id','description'])
gene2K = dict(zip(tbl[0],tbl[13]))
genelst = gene2K.keys()
Klst = gene2K.values()

##### create KEGG hierarchy

def crtNetwork(kegfile):
    G = nx.DiGraph()
    Z = json.load(open(kegfile))
    Zname = 'KEGG PATHWAY'
    G.add_node(Zname)
    for A in Z['children']:
        Aname = A['name']
        for B in A['children']:
            Bname = B['name']
            for C in B['children']:
                Cname = C['name'].split(' ',1)[0]
                if 'children' in C.keys():
                    G.add_edges_from([[Zname,Aname],[Aname,Bname],[Bname,Cname]])
                    Cchildren = '\n'.join(
                        map(lambda x:x['name'],C['children']))
                    Kset = re.findall('(?<=\t)[K\d]+(?= )',Cchildren)
                    G.add_edges_from([[Cname,K] for K in Kset])
    return G

print "generating Chimaera hierarchy network, org number is %i"%len(orglst)
allNet = nx.DiGraph()
for i in orglst:
        allNet=nx.compose(allNet,crtNetwork('%s/%s00001.json'%(PATH,i.replace('\n',''))))
print "network generation done..."

##### prune network

allNet.add_edges_from([[gene2K[g],g] for g in genelst])
print "begin prunning, %i edges exists before processing."%len(allNet.edges)
kolst = []
for K in gene2K.values():
    kolst+=allNet.predecessors(K)
kolst = list(set(kolst))
Kko = allNet.subgraph(kolst+Klst+genelst)
print "prunning done, %i edges left..."%len(Kko.edges)

##### enrichment
print "begin enrichment ..."
ko2desc = pd.read_table(MAP,header=None,dtype=str)
K2D = dict(zip(ko2desc[0],ko2desc[1]))
kolst = [ko for ko in kolst if ko in K2D.keys()]
outlst = []
for ko in kolst:
    mix = nx.descendants(Kko,ko)
    geneS = [g for g in mix if g in genelst]
    geneAll = ','.join(geneS)
    Kall = ','.join([K for K in mix if K in Klst])
    for path in nx.all_simple_paths(allNet,"KEGG PATHWAY",ko):
        Z,A,B,C = path
        outlst.append([A,B,K2D[ko],ko,len(geneS),geneAll,Kall])
outDF = pd.DataFrame(outlst)
outDF = outDF.sort_values(4,ascending=False)
outDF.to_csv('%s.KEGG.Path.xls'%smp,index=None,sep='\t',
             header=['KEGG_A_class','KEGG_B_class','Pathway',
                     'ID','Count(%i)'%len(genelst),'Genes','KOs'])
print "enrichment done ..."

##### generate HTML

print "generating HTML."
os.system("mkdir -p %s_map"%smp)
for ko in kolst:
    mapid = "map%s"%ko
    conf = "%s/map/%s.conf"%(dbpath,mapid)
    mapimg = "%s/map/%s.png"%(dbpath,mapid)
    source_img = Image.open(mapimg).convert("RGBA")
    draw = ImageDraw.Draw(source_img)
    sKlst = Kko.successors(ko)
    sKlst = [i for i in sKlst]
    sKdic = {}
    for line in open(conf).readlines():
        sublst = [i for i in sKlst if i in line]
        if (len(sublst)>0) & line.startswith('rect'):
            rect,url,desc = line.split('\t')
            Kpos = re.findall('[\d]+',rect)
            rect = ','.join(Kpos)
            Kpos = map(lambda x:int(x),Kpos)
            draw.rectangle(Kpos, outline="red")
            url = url.split('?')[1]
            sKdic[rect] = {}
            sKdic[rect]['label'] = []
            for K in sublst:
                sKdic[rect]['label'].append(K+': '+','.join([i for i in Kko.successors(K)]))
            sKdic[rect]['show']='</li><li>'.join(sKdic[rect]['label'])
            sKdic[rect]['pos']= rect
    template = Template(open('map.template').read())
    out = open('%s_map/%s.html'%(smp,mapid),'w')
    rendr = template.render(mapid=mapid,Klst=sKdic.values())
    out.write(rendr)
    out.close()
    source_img.save('%s_map/%s.png'%(smp,mapid), 'PNG')
print "generation done ..."






