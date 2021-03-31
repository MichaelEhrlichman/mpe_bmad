#!/usr/bin/env python

import re, os, numpy
import multiprocessing
import errno
import sys
from namelist import namelist2dict
import glob

common_in = namelist2dict("../../../common.in")

nProcs=multiprocessing.cpu_count()

#configuration
periodicity = int(common_in['evaluator'][0]['plot_periodicity'])
makeLats = common_in['evaluator'][0]['makeLats']
makeStock = common_in['evaluator'][0]['makeStock']
plotOnly = common_in['evaluator'][0]['plotOnly']
doLA = common_in['evaluator'][0]['doLA']
doNDP = common_in['evaluator'][0]['doNDP']
doDA = common_in['evaluator'][0]['doDA']
doRA = common_in['evaluator'][0]['doRA']
doTL = common_in['evaluator'][0]['doTL']
doFP = common_in['evaluator'][0]['doFP']
doADTS = common_in['evaluator'][0]['doADTS']
#end configuration

#scriptPath = '/gpfs/home/ehrlichman_m/afs_home/afs_bmad/bmad_dist_2014_0429/sls_sim/dynap_pisa_evaluator/'
scriptPath = '/afs/psi.ch/user/e/ehrlichman_m/afs_bmad/bmad_dist_2014_0429/sls_sim/dynap_pisa_evaluator/'

def try_makedirs(dirName):
	try:
		os.makedirs(dirName)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise e
		pass

def crunch(dirName):
	if doLA: 
		try_makedirs(dirName+'/00da_linear')
		os.system("cd "+dirName+"/00da_linear"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_linear ../../../../../common.in ../phase_correction.lat")
	if doNDP: 
		try_makedirs(dirName+'/00da_ndp')
		os.system("cd "+dirName+"/00da_ndp"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_ndp ../../../../../common.in ../phase_correction.lat")
	if doDA: 
		try_makedirs(dirName+'/00da')
		os.system("cd "+dirName+"/00da"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_slim ../../../../../common.in ../phase_correction.lat")
	if doRA: 
		try_makedirs(dirName+'/00da_raster')
		os.system("cd "+dirName+"/00da_raster"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_raster ../../../../../common.in ../phase_correction.lat")
	if doADTS: 
		try_makedirs(dirName+'/00adts')
		os.system("cd "+dirName+"/00adts"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/tracker_adts ../../../../../common.in ../phase_correction.lat")
	if doTL: 
		try_makedirs(dirName+'/00touschek')
		os.system("cd "+dirName+"/00touschek"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/aperture_and_lifetime ../../../../../common.in ../phase_correction.lat")
	if doFP: 
		try_makedirs(dirName+'/00footprint')
		os.system("cd "+dirName+"/00footprint"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/footprint ../../../../../common.in ../phase_correction.lat")

	os.system("cd "+dirName+";rm -f phase_correction.lat.digested*")

def plot(dirName):
	if doDA: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.da.gp")
	if doRA: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.raster.gp")
	if doADTS: 
		os.system("cd "+dirName+"; "+scriptPath+"plot.adts.py")
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.adts.gp")
	if doTL: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.ma.gp")
	if doFP: 
		os.system("cd "+dirName+"; "+scriptPath+"plot.pzfp.py "+str(periodicity))
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.trace.gp")
		if doADTS:
			os.system("cd "+dirName+"; "+scriptPath+"pretty.footprints.py "+str(periodicity))

dirs = glob.glob("seed_*")

pool = multiprocessing.Pool(processes=nProcs)
if not plotOnly:
	pool.map(crunch,dirs)
pool.map(plot,dirs)
pool.close()
pool.join()




