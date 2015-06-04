cmake_minimum_required (VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project (vegas)

set (vegas_VERSION_MAJOR 0)
set (vegas_VERSION_MINOR 1)

set (CMAKE_C_FLAGS "-std=c99")

set (CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -lprofiler")

include_directories("${PROJECT_BINARY_DIR}")

find_package (GSL REQUIRED)

add_library (mcvegas src/mc.c src/system.c src/util.c)
target_link_libraries (mcvegas m ${GSL_LIBRARIES})

include_directories ("${PROJECT_SOURCE_DIR}/include/")

add_executable (vegas src/main.c)
target_link_libraries (vegas mcvegas ${GSL_LIBRARIES})

find_package (CUDA 6.5 REQUIRED)

cuda_add_library (mcvegas_gpu src/mc.c src/system.c src/util.c)
target_link_libraries (mcvegas_gpu m ${GSL_LIBRARIES})

cuda_build_clean_target()

cuda_add_executable (vegas_gpu src/main.cu)
target_link_libraries (vegas_gpu mcvegas_gpu ${GSL_LIBRARIES} ${CUDA_curand_LIBRARY})

install (TARGETS mcvegas DESTINATION lib)
install (TARGETS mcvegas_gpu DESTINATION lib)
install (FILES src/util.h src/system.h src/mc.h DESTINATION include)
install (TARGETS vegas DESTINATION bin)
install (TARGETS vegas_gpu DESTINATION bin)

