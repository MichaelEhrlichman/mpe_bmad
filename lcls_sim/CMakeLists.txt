cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(lslc_sim)

set (SRC_DIRS
)

set (EXE_SPECS
 cmake.lc_toss
 cmake.lc_optics
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
