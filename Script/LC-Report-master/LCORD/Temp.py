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

SOURCEDIR = "/Bio/User/xiejunting/KitulTHings/miniReport/"

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
	-m	--method	(optional):  Ordination Method(str, can be a list, Default = "cca")
	-P,	--sigthrd	(optional):  Premu TEST Pvalue Result Table (float, Defalult = 0.05)

	-L,	--language	(optional):  Report Language (str, Default = EN)
	
	Project Information:
	-c	--client	(optional):  Client Name(str, Default = "Kidult")
	-i	--institute	(optional):  Institute Name(str, Default = "Sagene")
	-s	--samples	(optional):  Samples Name(str, Default = "")
	-g	--groups	(optional):  Groups Name(str, Default = "")
	-t	--pjectypes	(optional):  Types of This Project(str, Default = "SGA PAnalysis")
	-p	--project	(optional):  Project ID(str, Default = "SGA0000")

usage e.g.:
	python LCORD.py -h 
	python LCORD.py -a infile1 -b infile2
	python LCORD.py --infile1 infile1 --infile2 infile2

Use -h/--help option to see the usage.
"""
parser = OptionParser(usage)
parser.add_option('-a', '--infile1', type="string", dest="infile1",   help="Input file1")
parser.add_option('-b', '--infile2', type="string", dest="infile2",   help="Input file2")
parser.add_option('-o', '--outDir',  type="string", dest="outDir",    help="Working Direction", default = "./PAnalysis")
parser.add_option('-m', '--method',  type="string", dest="ordmethod", help="Correlation Method",default = "cca")
parser.add_option('-f', '--envfit',  type="string", dest="envfit",    help="Fit or NOT",        default = "nofit")
parser.add_option('-P', '--sigthrd', type="string", dest="pthreshold",help="p threshold",       default = "0.05")
parser.add_option('-c', '--client',  type="string", dest="client",    help="Client Name",       default = "Kidult")
parser.add_option('-i', '--institu', type="string", dest="institute", help="Institute Name",    default = "Sagene")
parser.add_option('-s', '--samples', type="string", dest="samples",   help="Samples Name",      default = "")
parser.add_option('-g', '--groups',  type="string", dest="groups",    help="Groups Name",       default = "")
parser.add_option('-L', '--language',type="string", dest="language",  help="Report Language",   default = "EN")
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
	print "-b/--infile2 is not necessary!! But You Don't Have infile2 You can't do the [CCA,RDA] That input should be 2 VarTables."
	infile2 = os.path.abspath(options.infile1)
else:
	infile2 = os.path.abspath(options.infile2)

if "-h" in dir():
	print usage
	sys.exit()

outDir      = options.outDir
ordmethod   = options.ordmethod
envfit      = options.envfit
pthreshold  = options.pthreshold

samples = options.samples
groups = options.groups 

language    = options.language

client    = options.client
institute = options.institute
types     = options.types
project   = options.project
contract  = options.contract

# jinja2
env = Environment(loader=FileSystemLoader(SOURCEDIR+"/LCORD"))
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
					groups=groups,
					types=types,
					ordmethod=ordmethod,
					envfit=envfit,
					pthreshold=pthreshold,
					infile1=infile1,
					infile2=infile2,
					outDir=os.path.abspath(outDir)
					)
		)
out.close()
