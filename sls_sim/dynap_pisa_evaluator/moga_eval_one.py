#!/usr/bin/env python

import re, os, numpy
import errno
import sys
from namelist import namelist2dict

common_params_file = "common.in"

#configuration
common_in = namelist2dict(common_params_file)
periodicity = int(common_in['evaluator'][0]['plot_periodicity'])
doLA = common_in['evaluator'][0]['doLA']
doNDP = common_in['evaluator'][0]['doNDP']
doDA = common_in['evaluator'][0]['doDA']
doRA = common_in['evaluator'][0]['doRA']
doTL = common_in['evaluator'][0]['doTL']
doFP = common_in['evaluator'][0]['doFP']
doADTS = common_in['evaluator'][0]['doADTS']
#end configuration

print("Plot periodicity is set to "+str(periodicity))

exePath = '/nfs/acc/user/ehrlichm/openmp_bmad/mpe_bmad/production/bin/'
scriptPath = '/home/ehrlichm/openmp_bmad/mpe_bmad/sls_sim/dynap_pisa_evaluator/'

def try_makedirs(dirName):
	try:
		os.makedirs(dirName)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise e
		pass

def plot():
	if doDA: 
		os.system("gnuplot "+scriptPath+"plot.da.gp")
	if doRA: 
		os.system("gnuplot "+scriptPath+"plot.raster.gp")
		os.system(scriptPath+"plot.fprint.py "+str(periodicity))
	if doADTS: 
		os.system(scriptPath+"plot.adts.py "+str(periodicity))
	if doTL: 
		os.system("gnuplot "+scriptPath+"plot.ma.gp")
	if doFP: 
		os.system(scriptPath+"plot.pzfp.py "+str(periodicity))
		if doADTS:
			os.system(scriptPath+"pretty.footprints.py "+str(periodicity))

def crunch():
	if doLA: 
		try_makedirs('00da_linear')
		os.system("cd 00da_linear"+";"+exePath+"/dynap_linear ../"+common_params_file)
	if doNDP: 
		try_makedirs('00da_ndp')
		os.system("cd 00da_ndp"+";"+exePath+"/dynap_ndp ../"+common_params_file)
	if doDA: 
		try_makedirs('00da')
		os.system("cd 00da"+";"+exePath+"/dynap_slim ../"+common_params_file)
	if doRA: 
		try_makedirs('00da_raster')
		os.system("cd 00da_raster"+";"+exePath+"/dynap_raster ../"+common_params_file)
	if doADTS: 
		try_makedirs('00adts')
		os.system("cd 00adts"+";"+exePath+"/tracker_adts ../"+common_params_file)
	if doTL: 
		try_makedirs('00touschek')
		os.system("cd 00touschek"+";"+exePath+"/aperture_and_lifetime ../"+common_params_file)
	if doFP: 
		try_makedirs('00footprint')
		os.system("cd 00footprint"+";"+exePath+"/footprint ../"+common_params_file)

#crunch()
plot()
