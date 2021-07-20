cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(sls_simple)

set (SRC_DIRS
)

set (EXE_SPECS
  cmake.ps_tune_perturb
# cmake.make_hybrid_lat
# #cmake.show_field
# cmake.toss
 cmake.toss_many
# cmake.toss_mc
# #cmake.toss_envelope
cmake.tunes
cmake.basic_stats
# cmake.make_2fam
# cmake.crm_frontend
 cmake.analytic_h
 cmake.amp_fact
 cmake.amp_fact2
# cmake.srdt_lsq_soln  #now in bsim
# cmake.ptc_via_bmad
# cmake.smat_sizes
cmake.co
# cmake.phase_space
 cmake.rad_int
 cmake.for_dcs
 cmake.linear_optics
# cmake.phase_beating
 cmake.mode3_optics
# cmake.linear_ma
 cmake.off_optics
 cmake.element_list
 cmake.extended_names
cmake.phys_apertures
# cmake.make_opa_eles
# cmake.grid
# cmake.one_turn_mat
# cmake.mat_a_to_b
# cmake.tilts
#cmake.cor_perturb
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
