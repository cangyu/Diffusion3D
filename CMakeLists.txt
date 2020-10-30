cmake_minimum_required(VERSION 3.10)

project(Diffusion3D VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3.7 REQUIRED)

add_library(SOLVER STATIC
	src/misc.cc
	src/property.cc
	src/diagnose.cc
	src/io.cc
	src/noc.cc
	src/geom.cc
	src/temporal.cc
	src/spatial.cc
	src/gradient.cc
	src/green-gauss.cc
	src/lsq.cc)

target_link_libraries(SOLVER PUBLIC Eigen3::Eigen)

add_executable(CAVITY
	app/main.cc
	case/cavity/ic.cc
	case/cavity/bc.cc)

target_link_libraries(CAVITY PUBLIC SOLVER)
install(TARGETS CAVITY RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)


add_executable(PIPE
	app/main.cc
	case/pipe/ic.cc
	case/pipe/bc.cc)

target_link_libraries(PIPE PUBLIC SOLVER)
install(TARGETS PIPE RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(GRADIENT-GG1 app/benchmark1.cc)
target_link_libraries(GRADIENT-GG1 PUBLIC SOLVER)
