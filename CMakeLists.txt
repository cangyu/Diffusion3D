cmake_minimum_required(VERSION 3.10)

project(Diffusion3D VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/usr/local/include/eigen3)

add_library(SOLVER
	src/io.cc
	src/noc.cc
	src/geom.cc)

add_executable(CAVITY
	src/main.cc
	case/cavity/ic.cc
	case/cavity/bc.cc)

add_executable(PIPE
	src/main.cc
	case/pipe/ic.cc
	case/pipe/bc.cc)

target_link_libraries(CAVITY PUBLIC SOLVER)

target_link_libraries(PIPE PUBLIC SOLVER)

install(TARGETS CAVITY RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

install(TARGETS PIPE RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
