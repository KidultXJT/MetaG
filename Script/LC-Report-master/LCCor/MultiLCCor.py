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
	-m	--method	(optional):  Correlation Method(str, can be a list, Default = "pearson,spearman" Can be a List)
	-d	--dist		(optional):  Distance Method(str, can be a list, Default = "euclidean,bray" Can be a List)
	-P,	--sigthrd	(optional):  Correlation TEST Pvalue Result Table (float, Defalult = 0.05)
	-R,	--corthrd	(optional):  Correlation Analysis R statistics Result Tabule (float, Default = 0.5)
	
	-L,	--language	(optional):  Report Language (str, Default = EN)
	
	Project Information:
	-c	--client	(optional):  Client Name(str, Default = "Kidult")
	-i	--institute	(optional):  Institute Name(str, Default = "Sagene")
	-s	--samples	(optional):  Samples Name(str, Default = "")
	-t	--pjectypes	(optional):  Types of This Project(str, Default = "SGA PAnalysis")
	-p	--project	(optional):  Project ID(str, Default = "SGA0000")

usage e.g.:
	python LCCor.py -h 
	python LCCor.py -a infile1 -b infile2
	python LCCor.py --infile1 infile1 --infile2 infile2

Use -h/--help option to see the usage.
"""
parser = OptionParser(usage)
parser.add_option('-a', '--infile1', type="string", dest="infile1",   help="Input file1")
parser.add_option('-b', '--infile2', type="string", dest="infile2",   help="Input file2")
parser.add_option('-o', '--outDir',  type="string", dest="outDir",    help="Working Direction", default = "./PAnalysis")
parser.add_option('-m', '--method',  type="string", dest="cormethod", help="Correlation Method",default = "spearman,pearson")
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
parser.add_option('-k', '--keepRMD', type="string", dest="keeprmd",   help="KeepRMD or Not",    default = "no")


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
keeprmd   = options.keeprmd

# Make R markdown to Analysis.
if os.path.exists(outDir):
	outDir = outDir
	#os.system("echo 'Pay ATTENTION please !! Your PATH {outDir} is Exist !!'".format(outDir = os.path.abspath(outDir)))
else:
	os.system("mkdir -p {outDir}".format(outDir = os.path.abspath(outDir)))

# Run Temp.py to make R markdown Script ~
for m in cormethod.split(","):
	for d in distmethod.split(","):
		os.system("python {SOURCEDIR}/LCCor/Temp.py \
				-a {infile1} \
				-b {infile2} \
				-m {cormethod} \
				-d {distmethod} \
				-P {pthreshold} \
				-R {rthreshold} \
				-o {outDir} \
				-L {language}".format(
					SOURCEDIR  = SOURCEDIR,
					infile1    = infile1,
					infile2    = infile2,
					cormethod  = cormethod,
					distmethod = distmethod,
					pthreshold = pthreshold,
					rthreshold = rthreshold,
					language   = language,
					outDir     = os.path.abspath(outDir)
					)
				)
		# RUN Analysis Script
		try:
			RMD = os.path.abspath(outDir)+"/PAnalysis.rmd"
			outputdir = os.path.abspath(outDir)
			outfile   = "Correlation.html"
			CMD = "Rscript -e " + '"' + "rmarkdown::render('" + RMD +"', output_dir = '" + outputdir  + "', output_file = '" + outfile + "')" + '"'
			os.system(CMD)
			DynameRMD = os.path.abspath(outDir)+"/PAnalysis."+m+".DynamicHM.rmd"
			outputdir = os.path.abspath(outDir) + "/Correlation/" + m
			outfile   = m + ".html"
			DynameCMD = "Rscript -e " + '"' + "rmarkdown::render('" + DynameRMD +"', output_dir = '" + outputdir  + "', output_file = '" + outfile + "')" + '"'
			os.system(DynameCMD)
		except:
			print "##################################################################################################################"
			print
			'[LANCE] Error in generating the report! in Correlation Method :: {cormethod} AND Distance Method :: {distmethod}'.format(cormethod = m, distmethod = d)
			print "##################################################################################################################"
			sys.exit(-1)


if keeprmd == "yes":
	os.system("echo 'Run Script is in {outDir}'".format(outDir=os.path.abspath(outDir)))
else:
	os.system("rm -rf {outDir}/*.rmd".format(outDir=os.path.abspath(outDir)))
	os.system("zip -r {outDir}.zip {outDir}".format(outDir = os.path.abspath(outDir)))


