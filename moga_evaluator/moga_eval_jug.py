#!/usr/bin/env python

import re, os, numpy
from jug import TaskGenerator
import errno
import sys
import yaml

cfg_file = "moga_eval.cfg"
common_params_file = 'common.in'

with open(cfg_file, 'r') as cfg_file_stream:
	try:
		cfg = yaml.safe_load(cfg_file_stream)
	except yaml.YAMLError as exc:
		print(exc)

elementNames = ()
elementProps = ()
for magnet,prop in zip(cfg['magnet_names'],cfg['magnet_props']):
	elementNames = elementNames+(magnet,)
	elementProps = elementProps+(prop,)

baseLat = '../../'+cfg['lat_file']
pop_file = 'moga_picked.out'
nObjs=4

#configuration
periodicity = cfg['periodicity']
makeLats = cfg['makeLats']
plotOnly = cfg['plotOnly']
includeBase = cfg['includeBase']
doLA = cfg['doLA']
doDA = cfg['doDA']
doRA = cfg['doRA']
doTL = cfg['doTL']
doFP = cfg['doFP']
doADTS = cfg['doADTS']
#end configuration

print("Plot periodicity is set to "+str(periodicity))
nElements=len(elementNames)

scriptPath = '/global/data/mike/bmad_dist_2019_0226/mpe_bmad/moga_evaluator/'

seedNames = []
dirs = []
seeds = []
with open(pop_file, 'r') as file:
	if includeBase:
		dirs.append('seed_0')
		seedNames.append('0')
		seeds.append('')
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

def try_makedirs(dirName):
	try:
		os.makedirs(dirName)
	except OSError as e:
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
	print("Crunching "+dirName+" ...")
	if makeLats:
		seedName = seedNames[n]
		print("       ... making seed "+seedName)
		seed = seeds[n]
		try_makedirs(dirName)

		latLat = open(dirName+'/'+'lat.lat', 'w')
		latLat.write('call, file = "'+baseLat+'"\n')
		if seedName != '0':
			latLat.write('call, file = "'+'moga_'+seedName+'.lat'+'"\n')
		latLat.close()

		if seedName != '0':
			seedLat = open(dirName+'/'+'moga_'+seedName+'.lat', 'w')
			for i in numpy.arange(nElements):
				seedLat.write(elementNames[i]+'['+elementProps[i]+'] = '+seed[i]+'\n')
			seedLat.close()
	if not plotOnly:
		if doLA: 
			try_makedirs(dirName+'/00da_linear')
			os.system("cd "+dirName+"/00da_linear"+";dynap_linear ../../../"+common_params_file+" ../lat.lat")
#		if doDA: 
#			try_makedirs(dirName+'/00da')
#			os.system("cd "+dirName+"/00da"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/dynap_slim ../../"+common_params_file+" ../lat.lat")
		if doRA: 
			try_makedirs(dirName+'/00da_raster')
			os.system("cd "+dirName+"/00da_raster"+";dynap_raster ../../../"+common_params_file+" ../lat.lat")
#		if doADTS: 
#			try_makedirs(dirName+'/00adts')
#			os.system("cd "+dirName+"/00adts"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/tracker_adts ../../"+common_params_file+" ../lat.lat")
#		if doTL: 
#			try_makedirs(dirName+'/00touschek')
#			os.system("cd "+dirName+"/00touschek"+";/afs/psi.ch/user/e/ehrlichman_m/bbin/aperture_and_lifetime ../../"+common_params_file+" ../lat.lat")
		if doFP: 
			try_makedirs(dirName+'/00footprint')
			os.system("cd "+dirName+"/00footprint"+";footprint ../../../"+common_params_file+" ../lat.lat")
		os.system("cd "+dirName+";rm -f lat.lat.digested*")
	plot(n)

for n in range(len(dirs)):
	crunch(n)

