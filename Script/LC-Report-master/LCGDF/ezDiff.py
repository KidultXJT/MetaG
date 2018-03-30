#!/usr/bin/env python
# -*- coding: UTF-8 -*-

##############################################
#                ezDiff Main                 #
#            ver 0.1   2018.01.08            #
#                                            #
#       Created by Liyan Yang @ Sagene       #
##############################################

import argparse
import os
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import time

SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]
WORK_DIR = os.getcwd()

################# Argument Parse #################

def argparser():

	description = '''
	*** LCDIFF ***
	'''
	parser = argparse.ArgumentParser(description=description, add_help=True)
	
	help = 'Input data table. Should be a tab-separated table in plain-text format. Each row represent a record (eg. a gene, a species or a OTU) and each column represents a sample. [Required]'
	parser.add_argument('-d', metavar="FILE", required = True, help = help)
	
	help = 'Differentiation test type. 1 - parametric tests; 2 - non-parametic tests; 3 - Fisher-Pitman permutation tests. [Required]'
	parser.add_argument('-t', metavar = "INT", required = True, choices = [1, 2, 3], help = help, type = int)

	help = 'Group labels for each of the sample separated by commas. The order of the labels should correspond to the order of the sample name in the input table. [Required]'
	parser.add_argument('-g', metavar = "STR", required = True, help = help)

	help = 'Output a HTML report along with results table and figures.'
	parser.add_argument('--html', action = "store_true", required = False, help = help)
	
	help = "Keep the intermediate RMD file."
	parser.add_argument('--keep_rmd', action = "store_true", required = False, help = help)

	help = "Prefix of the results files. [Default: TEST]"
	parser.add_argument('-p', metavar = "STR", required = False, default = "TEST", help = help)

	help = "Index of the column for the record name. Index starting from 1. [Default: 1]"
	parser.add_argument('--rn', metavar = "INT", required = False, default = 1, type = int, help = help)

	help = "Index of the columns to be removed when doing tests. Index starting from 1."
	parser.add_argument('--c_skip', metavar = "STR", required = False, default = "0", help = help)

	help = "Key word of the record. eg. Gene, OTU, Transcript, etc. [Default: Record]"
	parser.add_argument('-r', metavar = "STR", required = False, default = "Record", help = help)
	
	help = "Name of the measurments in the data table. eg. Abundance, Frequency, Counts, etc. [Default: Measurement]"
	parser.add_argument('-m', metavar = "STR", required = False, default = "Measurement", help = help)

	help = "Subset of the groups to be tests."
	parser.add_argument('--sel', required = False)

	pAdjMethods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"]

	help = "Method for intra-group p.value adjustment. The methods available are: \
			bonferroni: Bonferroni correction (1961); \
			holm: Holm's method (1979); \
			hochberg: Hochberg's method (1988); \
			hommel: Hommel's method (1988); \
			BH: Benjamini & Hochberg's method (1995); \
			BY: Benjamini & Yekutieli's method (2001); \
			fdr: same as BH; \
			none: no adjustment will be done. [Default: BH]"
	parser.add_argument('--padj1', required = False, choices = pAdjMethods, default = "BH", help = help)

	help = "Method for inter-group p.value adjustment. The methods available are the same as those for --padj2. [Default: bonferroni]"
	parser.add_argument('--padj2', required = False, choices = pAdjMethods, default = "bonferroni", help = help)

	help = "Variation equality of the records between groups are considered as same. =auto will do variation test for each of the records and do the differentiation test according to the results. [Default: auto]"
	parser.add_argument('--ev', metavar = "STR", required = False, default = "auto", choices = ["auto", "T", "F"], help = help)

	help = "The confidence interval for the tests and figures. [Default: 0.95]"
	parser.add_argument('--conf', metavar = "FLOAT", required = False, default = 0.95, type = float, help =help)

	help = "Paired measurements between groups."
	parser.add_argument('--pair', required = False, action = "store_true", help = help)

	help = "Filter threashold for group-wise tests results. [Default: 0.05]"
	parser.add_argument('--gfilter', metavar = "FLOAT", required = False, default = 0.05, type = float, help = help)

	help = "Filter threashold for pair-wise tests results. [Default: 0.05]"
	parser.add_argument('--pfilter', metavar = "FLOAT", required = False, default = 0.05, type = float, help = help)
	
	help = "Output directory for the group-wise tests results. [Default: the working dir]"
	parser.add_argument('--gdir', metavar = "DIR", required = False, default = WORK_DIR, help = help)

	help = "Output directory for the pair-wise tests results. [Default: the working dir]"
	parser.add_argument('--pdir', metavar = "DIR", required = False, default = WORK_DIR, help = help)

	help = "Do pair-wise Pitman-Fisher Permutation Test when test type is 3."
	parser.add_argument('--ptest', action = "store_true", help = help)

	help = "Output the qq plot for each variable in each group."
	parser.add_argument('--qq', action = "store_true", required = False, help = help)

	help = "Output directory for the qq plots"
	parser.add_argument('--qdir', metavar = "DIR", required = False, default = WORK_DIR, help = help)

	return parser.parse_args()

################# Main #################

print '[LANCE] Start working...'
print '[LANCE] Checking the arguments...'

args = argparser()
IN_FILE = os.path.realpath(args.d)
if (not os.path.exists(IN_FILE)):
	print '[ERROR] Input data does not exist!'
	sys.exit(-1)

try:
	d = pd.read_csv(IN_FILE, sep = '\t')
except:
	print '[ERROR] Error in input data!'
	sys.exit(-1)

GROUP = args.g
if (d.shape[1] != len(GROUP.split(',')) + 1):
	print '[ERROR] The group list does not match the column in the input table!'
	sys.exit(-1)

ROW_NAME = args.rn
if ROW_NAME > d.shape[1]:
	print '[ERROR] The column for the record name is out of index!'
	sys.exit(-1)

CSKIP = args.c_skip
for i in CSKIP.split(','):
	if int(i) > d.shape[1]:
		print '[ERROR] The column(s) to be skipped is out of index!'
		sys.exit(-1)

TESTS = ["param", "nparam", "perm"]
TESTS_NAME = ["Parametric tests", "Non-parametric tests", "Permutation tests"]
TEST_TYPE = TESTS[args.t-1]
TEST_TMPL = "ezDiff-" + TEST_TYPE + ".rmd.tmpl"
TEST_FILE = SCRIPT_DIR + "/" + TEST_TMPL
if (not os.path.exists(TEST_FILE)):
	print "[ERROR] Template file for " + TESTS_NAME[int(args.t)-1] + " not found!"
	sys.exit(-1)

SEL = args.sel if args.sel else "NULL"
if args.sel:
	for i in SEL.split(','):
		if i not in GROUP.split(','):
			print "Group(s) in the subset not found!"
			sys.exit(-1)
SEL = '"' + SEL + '"' if args.sel else SEL

EV = '"' + args.ev + '"' if args.ev == "auto" else args.ev

RMD = args.p + ".rmd"
HTML = args.p + ".html"

CONF = args.conf
if (CONF <= 0) | (CONF >= 1):
	print '[ERROR] Invalid confidence interval!'
	sys.exit(-1)

PAIR = 'T' if args.pair else 'F'

GFILTER = args.gfilter
PFILTER = args.pfilter
if (GFILTER <= 0) | (GFILTER >= 1):
	print '[ERROR] Invalid filter threashold for group-wise tests.'
	sys.exit(-1)
if (PFILTER <= 0) | (PFILTER >= 1):
	print '[ERROR] Invalid filter threashold for pair-wise tests.'
	sys.exit(-1)

GDIR = args.gdir
os.makedirs(GDIR) if not os.path.exists(args.gdir) else 0
PDIR = args.pdir
os.makedirs(PDIR) if not os.path.exists(args.pdir) else 0
QDIR = args.qdir
os.makedirs(QDIR) if (args.qq and (not os.path.exists(args.qdir))) else 0

OUTQQ = 'T' if args.qq else 'F'
PTEST = 'T' if args.ptest else 'F'

print '[LANCE] Arguments all right!'
print '[LANCE] Processing the test file...'

env = Environment(loader=FileSystemLoader(SCRIPT_DIR))
temp = env.get_template(TEST_TMPL)
out = open(RMD, "w")
out.write(
	temp.render(
		IN_FILE      = IN_FILE,
		RN           = args.rn,
		CSKIP        = CSKIP,
		GROUP        = GROUP,
		ID           = args.r,
		VAL          = args.m,
		SEL          = SEL,
		PADJ1        = args.padj1,
		PADJ2        = args.padj2,
		EV           = EV,
		CONF         = CONF,
		PAIR         = PAIR,
		GFIL         = GFILTER,
		PFIL         = PFILTER,
		OUTPFX       = args.p,
		GDIR         = GDIR,
		PDIR         = PDIR,
		PAIRTEST     = PTEST,
		DATE         = time.strftime('%Y-%m-%d %H:%M',time.localtime(time.time())),
		OUTQQ        = OUTQQ,
		QQDIR        = QDIR
	)
)
out.close()

print '[LANCE] Doing the test and generating the report...'

CMD = "Rscript -e " + '"' + "rmarkdown::render('" + RMD +"')" + '"' + " >/dev/null 2>&1"

try:
	os.system(CMD)
except:
	print '[LANCE] Error in doing the tests. Please contact the author.'
	sys.exit(-1)

if not args.keep_rmd:
	os.remove(RMD)

if not args.html:
	os.remove(HTML)

print '[LANCE] All Done! Thanks for using LANCE!'
