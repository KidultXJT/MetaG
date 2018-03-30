#!usr/bin/python
# -*- coding:utf-8 -*-

###########################################################
## Author: Kidult
## Email : junting.xie@sagene.com.cn
## Date  : 
###########################################################
## Your code go here

from optparse import OptionParser
from pandas import Series, DataFrame
import pandas as pd
import os
import sys
import re

from jinja2 import Environment, FileSystemLoader

SOURCEDIR = "/Bio/User/xiejunting/KitulTHings/LC-Report/"

################################# Arguments and Error #################################
usage = """
copyright: 
	Guangzhou Sagene Biotech Co.,Ltd

contact e-mail:
	junting.xie@sagene.com.cn

description:
	

options:
	-h	--help

	-a	--infile1	(necessary): The Infile1 (str)
	-b	--infile2	(necessary): The Infile2 (str)
	-o	--outDir	(optional):  The Output Direction (str, Default = "./PAnalysis")
	-m	--method	(optional):  Correlation Method(str, can be a list, Default = "pearson")
	-d	--dist		(optional):  Distance Method(str, can be a list, Default = "euclidean")
	-P, --sigthrd	(optional):  Correlation TEST Pvalue Result Table (float, Defalult = 0.05)
	-R, --corthrd	(optional):  Correlation Analysis R statistics Result Tabule (float, Default = 0.5)

	-L,	--language	(optional):  Report Language (str, Default = EN)

	Project Information:
	-c	--client	(optional):  Client Name(str, Default = "Kidult")
	-i	--institute	(optional):  Institute Name(str, Default = "Sagene")
	-s	--samples	(optional):  Samples Name(str, Default = "")
	-t	--pjectypes	(optional):  Types of This Project(str, Default = "SGA PAnalysis")
	-p	--project	(optional):  Project ID(str, Default = "SGA0000")

usage e.g.:
	python EzCor.py -h 
	python EzCor.py -a infile1 -b infile2
	python EzCor.py --infile1 infile1 --infile2 infile2

Use -h/--help option to see the usage.
"""
parser = OptionParser(usage)
parser.add_option('-a', '--infile1', type="string", dest="infile1",   help="Input file1")
parser.add_option('-b', '--infile2', type="string", dest="infile2",   help="Input file2")
parser.add_option('-o', '--outDir',  type="string", dest="outDir",    help="Working Direction", default = "./PAnalysis")
parser.add_option('-m', '--method',  type="string", dest="cormethod", help="Correlation Method",default = "pearson")
parser.add_option('-d', '--dist',    type="string", dest="distmethod",help="distance Method",   default = "euclidean")
parser.add_option('-P', '--sigthrd', type="string", dest="pthreshold",help="p threshold",       default = "0.05")
parser.add_option('-R', '--corthrd', type="string", dest="rthreshold",help="r threshold",       default = "0.5")
parser.add_option('-L', '--language',type="string", dest="language",  help="Report Language",   default = "EN")
parser.add_option('-c', '--client',  type="string", dest="client",    help="Client Name",       default = "Kidult")
parser.add_option('-i', '--institu', type="string", dest="institute", help="Institute Name",    default = "Sagene")
parser.add_option('-s', '--samples', type="string", dest="samples",   help="Samples Name",      default = "")
parser.add_option('-t', '--types',   type="string", dest="types",     help="Project Types",     default = "SGA PAnalysis")
parser.add_option('-p', '--project', type="string", dest="project",   help="Project ID",        default = "SGA0000")
parser.add_option('-C', '--contract',type="string", dest="contract",  help="Contract ID",       default = "SG0000000")


(options,args) = parser.parse_args()

if not options.infile1:
	print "-a/--infile1 is necessary!! please check your parameters."
	print usage
	sys.exit()
else:
	infile1 = os.path.abspath(options.infile1)

if not options.infile2:
	print "-a/--infile2 is necessary!! please check your parameters."
	print usage
	sys.exit()
else:
	infile2 = os.path.abspath(options.infile2)

if "-h" in dir():
	print usage
	sys.exit()

outDir      = options.outDir
cormethod   = options.cormethod
distmethod  = options.distmethod
pthreshold  = options.pthreshold
rthreshold  = options.rthreshold

language    = options.language

client    = options.client
institute = options.institute
samples   = options.samples
types     = options.types
project   = options.project
contract  = options.contract

# jinja2
env = Environment(loader=FileSystemLoader(SOURCEDIR+"/LCCor"))
if language == "EN":
	temp = env.get_template("Temp.rmd")
elif language == "CN":
	temp = env.get_template("Temp.rmd")
# Make R markdown to Analysis.
os.system("mkdir -p {outDir}".format(outDir = outDir))

out = open(os.path.abspath(outDir)+"/PAnalysis.rmd","w")
out.write(
		temp.render(contract=contract,
					project=project,
					client=client,
					institute=institute,
					samples=samples,
					types=types,
					cormethod=cormethod,
					distmethod=distmethod,
					pthreshold=pthreshold,
					rthreshold=rthreshold,
					infile1=infile1,
					infile2=infile2,
					outDir=os.path.abspath(outDir)
					)
		)
out.close()
# make R markdown to make Dynamic Result
method = cormethod.split(",")
for m in method: 
	env = Environment(loader=FileSystemLoader(SOURCEDIR+"/LCCor"))
	temp = env.get_template("DynamicTemp.rmd")
	if os.path.exists(outDir):
		outDir = outDir
	else:
		os.system("mkdir {outDir}".format(outDir = outDir))
	out = open(os.path.abspath(outDir)+"/PAnalysis."+m+".DynamicHM.rmd","w") # Dynamic Result
	out.write(
			temp.render(method=m,
				outDir=os.path.abspath(outDir),
				distmethod=distmethod
				)
			)
	out.close()

