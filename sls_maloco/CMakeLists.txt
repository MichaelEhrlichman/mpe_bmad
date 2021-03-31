set (LIBNAME sls_sim)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
  ../include
)

set (SRC_DIRS
  code
  modules
  namelists
)

set (EXE_SPECS
 cmake.write_quads
 cmake.dynap_slim
 cmake.dynap_linear
 cmake.dynap_raster
 cmake.aperture_and_lifetime
 cmake.scan_dp
 cmake.patch_maker
 cmake.tracker_adts
 cmake.footprint
 #needs coarray fortran.  cmake.tracker_footprint
 cmake.dynap_pz
 #needs to be updated for new naff interface. cmake.tracker_1
 cmake.tracker_simple
 cmake.orm
 cmake.prm
 cmake.ma_lat
 cmake.beta_beat
 #out of date. cmake.freq_map
)
set(EXECUTABLE_OUTPUT_PATH ../../production/bin)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
