set (LIBNAME sls_sim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(sls_sim)

set(INC_DIRS 
)

set (SRC_DIRS
  code
  modules
  namelists
)

set (EXE_SPECS
 #cmake.bend_maker
 #cmake.injection
# cmake.inj_dist_tao
# cmake.tracker_train
# cmake.special_chrom_set
 #cmake.write_quads
# cmake.dynap_slim
# cmake.dynap_linear
 cmake.dynap_raster
 cmake.aperture_and_lifetime
 #cmake.divergence_aperture
 #cmake.scan_dp
 #cmake.patch_maker
# cmake.tracker_adts
# cmake.footprint
# cmake.footprint_simple
#cmake.tracker_footprint
cmake.dynap_pz
cmake.tracker_1
 #fix ele_noise_struct_in cmake.tracker_noise_multi
# cmake.tracker_noise_single
# cmake.portrait
 cmake.tracker_simple
# cmake.tracker_dist
# cmake.orm
# cmake.orm_print
 #cmake.prm
# cmake.ma_lat
# cmake.beta_beat
 # cmake.freq_map
)

#set (FFLAGS "-pg")     # For profiling
#set (LINK_FLAGS "-pg")


include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)

#set(EXECUTABLE_OUTPUT_PATH ../../production/bin)




