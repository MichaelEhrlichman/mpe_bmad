cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(lslc_sim)

set (SRC_DIRS
)

set (EXE_SPECS
 cmake.lc_toss
 cmake.lc_unit_test
 cmake.und_toss
 cmake.lc_optics
# cmake.bug
# cmake.bug2
# cmake.bug3
 cmake.rmat
# cmake.checker1
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
