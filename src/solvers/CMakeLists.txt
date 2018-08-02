
INCLUDE_DIRECTORIES(
	${PARENT_DIR}/src/solvers
)


#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -fopenmp")


set (EXECUTABLE_OUTPUT_PATH ${INSTALL_DIR})

ADD_EXECUTABLE(luxum
	main.cc
	two_dim/maxwell_two_dim.cc
	two_dim/maxwell_two_dim_pml.cc
	three_dim/maxwell_three_dim.cc
	three_dim/maxwell_three_dim_pml.cc
	maxwell.cc
	)

#set_target_properties(luxum PROPERTIES COMPILE_FLAGS "-Wall -g -DDEBUG_RUN -DBZ_DEBUG")
set_target_properties(luxum PROPERTIES COMPILE_FLAGS "-O3 -std=c++11 -Wall -pedantic -fopenmp")
#set( CMAKE_CXX_FLAGS "-Wall -fopenmp -O3" )


TARGET_LINK_LIBRARIES(luxum fileio field parallel reader yaml-cpp hdf5)
