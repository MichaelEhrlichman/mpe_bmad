#!/bin/bash
#$ -cwd
#$ -N aperture_and_lifetime
#$ -pe orte 16
#$ -q prime_bd.q
#$ -l s_rt=4:39:00,h_rt=4:59:30
#$ -m abe -M michael.ehrlichman@psi.ch

source /opt/psi/config/profile.bash
export -n -f module
module use unstable
module add intel/16.3
module add gcc/4.9.3
module add openmpi/1.10.2

# For informational purposes
cat $PE_HOSTFILE | awk '{print $1" "$2}' > nodes.list

mpirun /afs/psi.ch/user/e/ehrlichman_m/bbin/aperture_and_lifetime ../../../common.in ../lat.lat


