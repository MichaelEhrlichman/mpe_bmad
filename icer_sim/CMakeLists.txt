cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(icer_sim)

set (SRC_DIRS
)

set (EXE_SPECS
#cmake.icer_optics
cmake.tracker_icer
cmake.tracker_icer_many
#cmake.rad_int
cmake.toss_alphap
cmake.icer_investigator
cmake.profile_simulator
#cmake.diffusion_simulator
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
