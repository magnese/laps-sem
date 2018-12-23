project(laps-sem)
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
aux_source_directory(. SRC_LIST)
set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_CXX_FLAGS "-Wall -pedantic -O3" ${CMAKE_CXX_FLAGS})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_PATH})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
