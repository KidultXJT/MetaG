#!/usr/bin/env python
# For MetaSSU NF workflow
# Format the output from LEfSe

import sys
import pandas as pd

help = "[Usage] python " + sys.argv[0] + " <input_file>"

if len(sys.argv[1:]) == 0:
	print help
	sys.exit(0)

if len(sys.argv[1:]) >1:
	print "One file each time please."
	print help
	sys.exit(0)

fIn = sys.argv[1]
df = pd.read_csv(fIn, sep = '\t', header = None)
df.columns = ["Taxonomy","Highest Mean","Groups(with Highest Mean)","LDA SCORE (log 10)","P-Value"]
df = df.sort_values(['Groups(with Highest Mean)','P-Value', 'LDA SCORE (log 10)'],ascending=[True,True, False])

df.ix[:,"ind"] = [not i.endswith('.') for i in df.iloc[:,0].tolist()]
df = df[df["ind"]]
df.dropna(inplace = True)
df.drop("ind", axis = 1, inplace = True)

fOut1 = sys.argv[1] + '.formatted'
fOut2 = sys.argv[1] + '.formatted.forplot'
df.to_csv(fOut1, sep = "\t", na_rep="-", header=True, index=False)
df.to_csv(fOut2, sep = "\t", na_rep="-", header=False, index=False)
