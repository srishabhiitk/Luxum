#Cmake setup
cmake_minimum_required (VERSION 2.8.3)
project (luxum)



#set (CMAKE_CXX_COMPILER mpic++)
get_filename_component (INSTALL_DIR ${CMAKE_BINARY_DIR} PATH)



set (HEADER_DIR ${INSTALL_DIR}/lib)
include_directories(${HEADER_DIR}/field
		            ${HEADER_DIR}/io
		            ${HEADER_DIR}/parallel)

message(STATUS "${HEADER_DIR}")

#OpenMP
find_package(OpenMP)
message(status ${OpenMP_C_FLAGS})




#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -fopenmp")
set( CMAKE_CXX_FLAGS "-Wall -fopenmp -O3" )


ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(src)
