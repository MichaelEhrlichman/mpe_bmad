cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
)

set (SRC_DIRS
)

set (EXE_SPECS
 cmake.tracker_train
)

#set (FFLAGS "-pg")     # For profiling
#set (LINK_FLAGS "-pg")


include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)

#set(EXECUTABLE_OUTPUT_PATH ../../production/bin)




