#!/usr/bin/env python

import sys
import shutil as sh
import subprocess
import random

template = sys.argv[1]
nstart = int(sys.argv[2])
nerrors = int(sys.argv[3])

print("template: {}".format(template))
print("nerrors:  {}".format(nerrors))

for i in range(nstart,nerrors+nstart):
	newdir = "{}_{}".format(template,i)
	sh.copytree(template,newdir)
	subprocess.call(["sed -i '/magnet_error_seed/c\  magnet_error_seed = {}' common.in".format(random.randint(1,99999))],shell=True,cwd=newdir)
	p = subprocess.Popen(["error_tool", "common.in", "lat.lat"], cwd=newdir)
	p.wait()
	p = subprocess.Popen(["beta_beat", "lat.lat", "errors.lat"], cwd=newdir)
	p.wait()
	p = subprocess.Popen(["sbatch", "submit.sl"], cwd=newdir+"/00touschek")
	p.wait()
