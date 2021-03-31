#!/usr/bin/env python

import re, os, numpy
from jug import TaskGenerator
import errno
import sys
from namelist import namelist2dict

common_params_file = "../common.in"

common_in = namelist2dict(common_params_file)

elementNames = ()
elementProps = ()
for mag in common_in['nl_moga'][0]['mags_in']:
	elementNames = elementNames+(mag[1],)
	elementProps = elementProps+(mag[2],)

baseLat = '../../'+common_in['general'][0]['lat_file']
pop_file = 'moga_picked.out'
nObjs=3

#configuration
periodicity = int(common_in['evaluator'][0]['plot_periodicity'])
makeLats = common_in['evaluator'][0]['makeLats']
plotOnly = common_in['evaluator'][0]['plotOnly']
doLA = common_in['evaluator'][0]['doLA']
doNDP = common_in['evaluator'][0]['doNDP']
doDA = common_in['evaluator'][0]['doDA']
doRA = common_in['evaluator'][0]['doRA']
doTL = common_in['evaluator'][0]['doTL']
doFP = common_in['evaluator'][0]['doFP']
doADTS = common_in['evaluator'][0]['doADTS']
#end configuration

print "Plot periodicity is set to "+str(periodicity)
nElements=len(elementNames)

#scriptPath = '/gpfs/home/ehrlichman_m/afs_home/afs_bmad/bmad_dist_2014_0429/sls_sim/dynap_pisa_evaluator/'
scriptPath = '/afs/psi.ch/user/e/ehrlichman_m/afs_bmad/bmad_local/mpe_bmad/sls_sim/dynap_pisa_evaluator/'

seedNames = []
dirs = []
seeds = []
file = open(pop_file, 'r')
for line in file:
	fileLine = line.split()
	if len(fileLine) == 0 or fileLine[0][0] == "#": 
		continue

	seedNames.append(fileLine[0])
	dirs.append('seed_'+fileLine[0])

	seedVars=["" for x in range(nElements)]
	for i in numpy.arange(1,nElements+1):
		seedVars[i-1] = fileLine[i]
	seeds.append(seedVars)
file.close()

def try_makedirs(dirName):
	try:
		os.makedirs(dirName)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise e
		pass

def plot(n):
	dirName = dirs[n]
	if doDA: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.da.gp")
	if doRA: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.raster.gp")
		os.system("cd "+dirName+"; "+scriptPath+"plot.fprint.py "+str(periodicity))
	if doADTS: 
		os.system("cd "+dirName+"; "+scriptPath+"plot.adts.py "+str(periodicity))
		#os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.adts.gp")
	if doTL: 
		os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.ma.gp")
	if doFP: 
		os.system("cd "+dirName+"; "+scriptPath+"plot.pzfp.py "+str(periodicity))
		#os.system("cd "+dirName+";gnuplot "+scriptPath+"plot.trace.gp")
		if doADTS:
			os.system("cd "+dirName+"; "+scriptPath+"pretty.footprints.py "+str(periodicity))

@TaskGenerator
def crunch(n):
	dirName = dirs[n]
	print "Crunching "+dirName+" ..."
	if makeLats:
		seedName = seedNames[n]
		print "       ... making seed "+seedName
		seed = seeds[n]
		try_makedirs(dirName)

		latLat = open(dirName+'/'+'lat.lat', 'w')
		latLat.write('call, file = "'+baseLat+'"\n')
		latLat.write('call, file = "'+'moga_'+seedName+'.lat'+'"\n')
		latLat.close()

		seedLat = open(dirName+'/'+'moga_'+seedName+'.lat', 'w')
		for i in numpy.arange(nElements):
			seedLat.write(elementNames[i]+'['+elementProps[i]+'] = '+seed[i]+'\n')
		seedLat.close()
	if not plotOnly:
		if doLA: 
			try_makedirs(dirName+'/00da_linear')
			os.system("cd "+dirName+"/00da_linear"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_linear ../../"+common_params_file+" ../lat.lat")
		if doNDP: 
			try_makedirs(dirName+'/00da_ndp')
			os.system("cd "+dirName+"/00da_ndp"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_ndp ../../"+common_params_file+" ../lat.lat")
		if doDA: 
			try_makedirs(dirName+'/00da')
			os.system("cd "+dirName+"/00da"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_slim ../../"+common_params_file+" ../lat.lat")
		if doRA: 
			try_makedirs(dirName+'/00da_raster')
			os.system("cd "+dirName+"/00da_raster"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_raster ../../"+common_params_file+" ../lat.lat")
		if doADTS: 
			try_makedirs(dirName+'/00adts')
			os.system("cd "+dirName+"/00adts"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/tracker_adts ../../"+common_params_file+" ../lat.lat")
		if doTL: 
			try_makedirs(dirName+'/00touschek')
			os.system("cd "+dirName+"/00touschek"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/aperture_and_lifetime ../../"+common_params_file+" ../lat.lat")
		if doFP: 
			try_makedirs(dirName+'/00footprint')
			os.system("cd "+dirName+"/00footprint"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/footprint ../../"+common_params_file+" ../lat.lat")
		os.system("cd "+dirName+";rm -f lat.lat.digested*")
	plot(n)

for n in range(len(dirs)):
	crunch(n)

