import pandas as pd
import sys
import re
import networkx as nx
import itertools

blastout = sys.argv[1]
obofile = sys.argv[2]
smp = blastout.split('.')[0]

def merge(a,b):
    return a+b

##### create GO network

three = ['GO:0005575','GO:0003674','GO:0008150']

dic = {'id':[],'name':'name','parents':[]}
G = nx.DiGraph()

for line in open(obofile).readlines():
    if '[Term]' in line:
        G.add_nodes_from(dic['id'],name=dic['name'])
        G.add_edges_from([i for i in itertools.product(dic['parents'],dic['id'])])
        dic = {'id':[],'name':'-','parents':[]}
    if '[Typedef]' in line:
        break
    if re.findall('(?<=^id:)(.*?)(?=\n)',line):
        dic['id'].append(re.findall('(?<=^id:)(.*?)(?=\n)',line)[0].strip())
    if re.findall('(?<=^alt_id:)(.*?)(?=\n)',line):
        dic['id'].append(re.findall('(?<=^alt_id:)(.*?)(?=\n)',line)[0].strip())
    if re.findall('(?<=^name:)(.*?)(?=\n)',line):
        dic['name']=re.findall('(?<=^name:)(.*?)(?=\n)',line)[0].strip()
    if re.findall('(?<=^is_a:)(.*?)(?=\n)',line):
        dic['parents'].append(re.findall('GO:[\d]+',line)[0].strip())
    if re.findall('(?<=^relationship: part_of)(.*?)(?=\n)',line):
        dic['parents'].append(re.findall('GO:[\d]+',line)[0].strip())

sole = filter(lambda i:False if G.in_edges(i) else True,G.nodes)
G.add_edges_from([i for i in itertools.product(['GO:0003673'],sole)])
N = nx.shortest_path(G)
#nx.write_gml(G,'test.gml')

#### enrich parser

def lineDeal(line):
    lnsplt = line.split('\t')
    gene = lnsplt[0]
    golst = re.findall('(?<=\[)GO:[\d]+',lnsplt[12])
    expand = filter(lambda x:len(x)<=5,[N['GO:0003673'][i] for i in golst])
    expand.append(golst)
    st = []
    for i in expand:
        st+=i
    expand = list(set(st))
    out = [[i,gene,G.nodes[i]['name']] for i in expand]
    return out

out = pd.DataFrame(reduce(merge,
                [lineDeal(line) for line in open(blastout).readlines()]))

out[out[0].apply(lambda x:nx.has_path(G,'GO:0008150',x))].\
    to_csv('%s.P'%smp,index=None,sep='\t',header=None)
out[out[0].apply(lambda x:nx.has_path(G,'GO:0003674',x))].\
    to_csv('%s.F'%smp,index=None,sep='\t',header=None)
out[out[0].apply(lambda x:nx.has_path(G,'GO:0005575',x))].\
    to_csv('%s.C'%smp,index=None,sep='\t',header=None)

#### blast results

df = pd.read_table(blastout,header=None,dtype=str)
df[12] = df[12].apply(
    lambda x:re.findall('(?<=\[)GO:(.*?)(?=\])',x))
df[12] = df[12].apply(
    lambda x:' || '.join(map(lambda y:'['+y+']',x)))
df[[0,1,2,3,10,11,12]].to_csv('%s.GO.Anno.xls'%smp,index=None,sep='\t',quoting=3,
                              header=['query','subject','identity(%)','align_len',
                                      'Evalue','score','GO annotation'])

