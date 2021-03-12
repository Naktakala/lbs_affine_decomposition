cmake_minimum_required(VERSION 3.2)
 
set(TARGET affine)
project(${TARGET} C CXX)
set(CHI_TECH_DIR "~/Desktop/ChiTech/chi-tech")
include("${CHI_TECH_DIR}/ChiResources/Macros/Downstream.cmake")
 
file (GLOB_RECURSE SOURCES "*.cc")
file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/LinearBoltzmannSolver/*.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES})
add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})
