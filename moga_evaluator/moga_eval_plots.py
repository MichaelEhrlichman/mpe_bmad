#!/usr/bin/env python3

import re, os, numpy
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
includeBase = cfg['includeBase']
doLA = cfg['doLA']
doDA = cfg['doDA']
doRA = cfg['doRA']
doXPZ = cfg['doXPZ']
doTL = cfg['doTL']
doFP = cfg['doFP']
doADTS = cfg['doADTS']
#end configuration

print("Plot periodicity is set to "+str(periodicity))
nElements=len(elementNames)

#scriptPath = '/global/data/mike/bmad_dist_2019_0226/mpe_bmad/moga_evaluator/'
#scriptPath = '~/bmad_dist_2019_0521/mpe_bmad/moga_evaluator/'
scriptPath = os.environ['ACC_ROOT_DIR']+'/mpe_bmad/moga_evaluator/'

dirs = []
with open(pop_file, 'r') as file:
	if includeBase:
		dirs.append('seed_0')
	for line in file:
		fileLine = line.split()
		if len(fileLine) == 0 or fileLine[0][0] == "#": 
			continue

		dirs.append('seed_'+fileLine[0])

def pound(b):
	if b:
		return ''
	else:
		return '# '

for n in range(len(dirs)):
	dirName = dirs[n]
	jobsh = open(dirName+'/'+'plot.sh', 'w')
	jobsh.write("cd "+dirName+"\n")
	jobsh.write(pound(doDA)+"gnuplot "+scriptPath+"plot.da.gp"+"\n")
	jobsh.write(pound(doRA)+"gnuplot "+scriptPath+"plot.raster.gp"+"\n")
	jobsh.write(pound(doRA)+"gnuplot "+scriptPath+"plot.fma.gp"+"\n")
	jobsh.write(pound(doRA)+scriptPath+"plot.fprint.py "+str(periodicity)+"\n")
	jobsh.write(pound(doXPZ)+"gnuplot "+scriptPath+"plot.raster_pz.gp"+"\n")
	jobsh.write(pound(doXPZ)+"gnuplot "+scriptPath+"plot.fma_pz.gp"+"\n")
	jobsh.write(pound(doADTS)+scriptPath+"plot.adts.py "+str(periodicity)+"\n")
	#jobsh.write(pound(doADTS)+";gnuplot "+scriptPath+"plot.adts.gp"+"\n")
	jobsh.write(pound(doFP)+scriptPath+"plot.pzfp.py "+str(periodicity)+"\n")
	#jobsh.write(pound(doFP)+"gnuplot "+scriptPath+"plot.trace.gp")
	jobsh.write(pound(doADTS)+pound(doDA)+scriptPath+"pretty.footprints.py "+str(periodicity)+"\n")
	jobsh.write(pound(doTL)+"gnuplot "+scriptPath+"plot.ma.gp"+"\n")
	jobsh.close()

