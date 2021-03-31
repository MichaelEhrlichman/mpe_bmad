#!/usr/bin/env python3

import re, os, numpy
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

#configuration
periodicity = cfg['periodicity']
makeLats = cfg['makeLats']
applyErrors = cfg['applyErrors']
includeBase = cfg['includeBase']
doLA = cfg['doLA']
doDA = cfg['doDA']
doRA = cfg['doRA']
doXPZ = cfg['doXPZ']
doTL = cfg['doTL']
doFP = cfg['doFP']
doADTS = cfg['doADTS']
#end configuration

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

def pound(b):
	if b:
		return ''
	else:
		return '# '

for n in range(len(dirs)):
	dirName = dirs[n]
	print("Staging "+dirName+" ...")
	if makeLats:
		seedName = seedNames[n]
		print("       ... making seed "+seedName)
		seed = seeds[n]
		try_makedirs(dirName)

		latLat = open(dirName+'/'+'lat_pre.lat', 'w')
		latLat.write('call, file = "'+baseLat+'"\n')
		if seedName != '0':
			latLat.write('call, file = "'+'moga_'+seedName+'.lat'+'"\n')
		latLat.write('*[x_limit]=0.006 ! 0.03\n')
		latLat.write('*[y_limit]=0.006 ! 0.03\n')
		latLat.write('D11[x_limit]=0.003\n')
		latLat.write('D11[y_limit]=0.003\n')
		latLat.close()

		latLat = open(dirName+'/'+'lat.lat', 'w')
		latLat.write('call, file = "'+baseLat+'"\n')
		if seedName != '0':
			latLat.write('call, file = "'+'moga_'+seedName+'.lat'+'"\n')
		latLat.write('*[x_limit]=0.006 ! 0.03\n')
		latLat.write('*[y_limit]=0.006 ! 0.03\n')
		latLat.write('D11[x_limit]=0.003\n')
		latLat.write('D11[y_limit]=0.003\n')
		if applyErrors:
			latLat.write('call, file = errors.lat\n')
		latLat.close()

		if seedName != '0':
			seedLat = open(dirName+'/'+'moga_'+seedName+'.lat', 'w')
			for i in numpy.arange(nElements):
				seedLat.write(elementNames[i]+'['+elementProps[i]+'] = '+seed[i]+'\n')
			seedLat.close()
	try_makedirs(dirName+'/00da_linear')
	try_makedirs(dirName+'/00da')
	try_makedirs(dirName+'/00da_raster')
	try_makedirs(dirName+'/00xpz_raster')
	try_makedirs(dirName+'/00adts')
	try_makedirs(dirName+'/00touschek')
	try_makedirs(dirName+'/00footprint')
	jobsh = open(dirName+'/'+'evaluate.sh', 'w')
	if applyErrors:
		jobsh.write("cd "+dirName+";error_tool ../../common.in lat_pre.lat"+"\n")
		jobsh.write("basic_stats lat.lat"+"\n")
		useLattice = '../lat.lat'
	else:
		jobsh.write("cd "+dirName+";basic_stats lat.lat"+"\n")
		useLattice = '../lat.lat'
	jobsh.write("cd 00da_linear"+"\n")
	jobsh.write(pound(doLA)+"cd ../00da_linear"+";dynap_linear ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doDA)+"cd ../00da"+";dynap_slim ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doRA)+"cd ../00da_raster"+";dynap_raster ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doXPZ)+"cd ../00xpz_raster"+";dynap_pz ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doADTS)+"cd ../00adts"+";tracker_adts ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doFP)+"cd ../00footprint"+";footprint ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write(pound(doTL)+"cd ../00touschek"+";aperture_and_lifetime ../../../"+common_params_file+" "+useLattice+"\n")
	jobsh.write("cd .."+";rm -f lat.lat.digested*"+"\n")
	jobsh.close()

