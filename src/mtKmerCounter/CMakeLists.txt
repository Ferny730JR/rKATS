cmake_minimum_required(VERSION 3.7.1)
project(KATSS_KMERCOUNTER
	VERSION 0.9.0
	DESCRIPTION "K-mer counting function utilities."
	LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED TRUE)

# Set cmake variables 
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/KATSS" CACHE PATH "Installation directory for headers")

set(TINYCTHREAD_DISABLE_TESTS ON CACHE BOOL "Disable TinyCThread unit tests")
set(TINYCTHREAD_INSTALL OFF CACHE BOOL "Install a static library for TinyCThread")

option(VERBOSE "Output error/warning messages to stderr")
option(SKIP_INSTALL_ALL "Don't install any targets")
option(SKIP_INSTALL_LIBRARIES "Don't install shared and static libraries")
option(SKIP_INSTALL_STATIC "Don't install static library")
option(SKIP_INSTALL_SHARED "Don't install shared library")
option(SKIP_INSTALL_HEADER "Don't install header files")

add_compile_definitions(KATSS_VERBOSE=${VERBOSE})

#Include dependencies
set(TINYCTHREAD source/tinycthread)
#add_subdirectory(${TINYCTHREAD})
find_package(ZLIB REQUIRED)

# Set optimization flags for release builds
set(CMAKE_C_FLAGS_RELEASE "-O3 -flto")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -flto")

# Set source files for library
set(KATSS_SOURCE_FILES
	"source/hash_functions.c"
	"source/tables.c"
	"source/seqseq.c"
	"source/counter.c"
	"source/uncounter.c"
	"source/enrichments.c"
	"source/helpers/rnafiles.c"
	"source/helpers/memory_utils.c"
	"source/tinycthread/source/tinycthread.c")

# Create objects for source files
add_library(OBJLIB OBJECT "${KATSS_SOURCE_FILES}")

# Shared libraries need PIC
set_property(TARGET OBJLIB PROPERTY POSITION_INDEPENDENT_CODE 1)

# Compiler optimizations
target_compile_options(OBJLIB PRIVATE -Wall -O3)

# Specify include directories
target_include_directories(OBJLIB PRIVATE
	${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/helpers ${TINYCTHREAD}/source)

# Create shared and static library
add_library(katss_shared SHARED $<TARGET_OBJECTS:OBJLIB>)
add_library(katss_static STATIC $<TARGET_OBJECTS:OBJLIB>)

# Link library with zlib
target_link_libraries(katss_shared ZLIB::ZLIB)
target_link_libraries(katss_static ZLIB::ZLIB)

# Set public header/API for library
set(KATSS_PUBLIC_HEADERS
	"include/counter.h"
	"include/hash_functions.h"
	"include/enrichments.h"
	"include/seqseq.h")

# Set library properties
set_target_properties(katss_shared PROPERTIES
	VERSION ${PROJECT_VERSION}
	OUTPUT_NAME katss)
set_target_properties(katss_static PROPERTIES
	VERSION ${PROJECT_VERSION}
	OUTPUT_NAME katss)

# Install the library 
if(NOT SKIP_INSTALL_LIBRARIES AND NOT SKIP_INSTALL_ALL)
	if(NOT SKIP_INSTALL_STATIC)
		install(TARGETS katss_static LIBRARY DESTINATION ${INSTALL_LIB_DIR})
	endif()
	if(NOT SKIP_INSTALL_SHARED)
		install(TARGETS katss_shared LIBRARY DESTINATION ${INSTALL_LIB_DIR})
	endif()
endif()

if(NOT SKIP_INSTALL_HEADER AND NOT SKIP_INSTALL_ALL)
	install(FILES ${KATSS_PUBLIC_HEADERS} DESTINATION ${INSTALL_INC_DIR})
endif()
