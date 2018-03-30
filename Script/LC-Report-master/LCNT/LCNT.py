#!/usr/bin/env python
# -*- coding: UTF-8 -*-

##############################################
#         LCNT - Lance Noramlity Test        #
#            ver 1.0   2017.11.16            #
#                                            #
#         Created by Jeffrey @ Sagene        #
#       Email: liyan.yang@sagene.com.cn      #
##############################################


import argparse
import os
from jinja2 import Environment, FileSystemLoader
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import os
import os.path
import pandas as pd

scriptDir = os.path.split(os.path.realpath(__file__))[0]
workDir = os.getcwd()

################# Argument Parse #################

def argparser():

	description = '''
This tool is used to generate the normality test report
for statistical tabels with rows representing different
variables and columns representing different observati-
ons. This tool will create a result folder with tables
and plots in it and a html report.

If you have any questions and problems in using this
tool, please contact:

liyan.yang@sagene.com.cn
'''
	parser = argparse.ArgumentParser(description=description, add_help=True)

	help = 'The input tab-separated table for the normality test. The first row should be the header with each row following represents a variable. The first column should be the name of the features and each column following represents a single observation.'
	parser.add_argument('-f', required = True, help = help)

	help = 'Group label for each of the observation in the script separated by commas. '
	parser.add_argument('-g', required = True, help = help)

	help = 'LCNT RMD template. [default: LCNT.rmd.tmpl in the folder of this script]'
	parser.add_argument('-m', default = scriptDir + '/LCNT.rmd.tmpl', help = help)

	help = 'The prefix of the result tabels and plots. If not designated, the default will be the name of the input table'
	parser.add_argument('-o', help = help)

	help = 'Whether to output the P-P and Q-Q plots in the result folder.'
	parser.add_argument('-p', choices = ['none','qq','pp','both'], default = 'both', help = help)

	help = 'Whether to output the test result by test. Each output table will have columns representing variables and rows representing groups.'
	parser.add_argument('-t', action = 'store_true', help = help)

	help = 'Whether to output the test result by variables. Each output table will have columns representing tests types and rows representing groups.'
	parser.add_argument('-v', action = 'store_true', help = help)

	args = parser.parse_args()
	return args

################# Main #################

print '[LANCE] Start working...'
print '[LANCE] Checking the arguments...'

args = argparser()

IN_FILE = os.path.realpath(args.f)
if (not os.path.exists(IN_FILE)):
	print 'Input data does not exist!'
	sys.exit(-1)
try:
	d = pd.read_csv(IN_FILE, sep = '\t')
except:
	print '[ERROR] Error in input data!'
	sys.exit(-1)
IN_FILE = '"' + IN_FILE + '"'

GROUP_LIST = args.g
if (d.shape[1] != len(GROUP_LIST.split(',')) + 1): 
	print '[ERROR] The group list does not match the column in the input table!'
	sys.exit(-1)
else: GROUP_LIST = '"' + GROUP_LIST + '"'

TMPL = os.path.realpath(args.m)
if (not os.path.exists(TMPL)):
	print '[ERROR] RMD Template does not exist!'
	sys.exit(-1)
tmplDir = os.path.split(TMPL)[0]
tmplFile = os.path.split(TMPL)[1]

PREFIX = args.o if (args.o) else os.path.splitext(os.path.split(os.path.realpath(args.f))[1])[0]
htmlFile = PREFIX + '_LCNT.Report.html'
resDir = workDir + '/' + PREFIX + '_LCNT_Report'
RMD = PREFIX + '_LCNT.Report.rmd'

PREFIX = '"' + PREFIX + '"'
OUT_PLOT = '"' + args.p + '"'
OUT_BY_TEST = 'TRUE' if args.t else 'FALSE'
OUT_BY_VAR = 'TRUE' if args.v else 'FALSE'

print '[LANCE] Arguments all right!'
print '[LANCE] Rendering the RMD File...'

env = Environment(loader=FileSystemLoader(tmplDir))
temp = env.get_template(tmplFile)
out = open(RMD, "w")
out.write(
		temp.render(
		OUT_PLOT     = OUT_PLOT,
		GROUP_LIST   = GROUP_LIST,
		PREFIX       = PREFIX,
		OUT_BY_TEST  = OUT_BY_TEST,
		OUT_BY_VAR   = OUT_BY_VAR,
		IN_FILE      = IN_FILE
		)
)
out.close()

print '[LANCE] RMD File Rendered!'
print '[LANCE] Doing the test and generating the report...'
CMD = "Rscript -e " + '"' + "rmarkdown::render('" + RMD +"')" + '"' + " >/dev/null 2>&1"

try:
	os.system(CMD)
except:
	print '[LANCE] Error in generating the report!'
	sys.exit(-1)

print '[LANCE] Report finished!'
print '[LANCE] Packing the results...'

if os.path.exists(resDir):
	os.system('rm -rf ' + resDir)
os.makedirs(resDir)

if args.p == 'pp' or args.p == 'both':
	if not os.path.exists(workDir + '/PP-Plots'):
		print '[WARNING] PP Plots were manually deleted!'
	else: os.system('mv PP-Plots ' + resDir)

if args.p == 'qq' or args.p == 'both':
	if not os.path.exists(workDir + '/QQ-Plots'):
		print '[WARNING] QQ Plots were manually deleted!'
	else: os.system('mv QQ-Plots ' + resDir)

if args.t:
	if not os.path.exists(workDir + '/BY_TEST'):
		print '[WARNING] Tables by tests were mannually deleted!'
	else: os.system('mv BY_TEST ' + resDir)

if args.v:
	if not os.path.exists(workDir + '/BY_VAR'):
		print '[WARNING] Tables by variables were mannually deleted!'
	else: os.system('mv BY_VAR ' + resDir)

if not os.path.exists(htmlFile):
	print '[WARNING] Report HTML was mannually deleted!'
else: os.system('mv ' + htmlFile + ' ' + resDir)

if not os.path.exists(resDir):
	print '[WARNING] Report Folder was mannually deleted!'
else: os.system('tar -jcf ' + os.path.split(resDir)[1] + '.tbz ' + os.path.split(resDir)[1])

print '[LANCE] All Done! Thanks for using LANCE!'

