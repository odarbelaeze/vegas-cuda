cmake_minimum_required (VERSION 2.8)

project (vegas)

set (vegas_VERSION_MAJOR 0)
set (vegas_VERSION_MINOR 1)

set (CMAKE_C_FLAGS "-std=c99")

set (CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS} -pg")

include_directories("${PROJECT_BINARY_DIR}")

add_library (mcvegas src/mc.c src/system.c src/util.c)

target_link_libraries (mcvegas m gsl cblas)

include_directories ("${PROJECT_SOURCE_DIR}/include/")

add_executable (vegas src/main.c)
target_link_libraries (vegas mcvegas gsl cblas)

install (TARGETS mcvegas DESTINATION lib)
install (FILES src/util.h src/system.h src/mc.h DESTINATION include)
install (TARGETS vegas DESTINATION bin)

