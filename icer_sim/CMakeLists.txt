cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(icer_sim)

set (SRC_DIRS
)

set (EXE_SPECS
cmake.icer_optics
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
