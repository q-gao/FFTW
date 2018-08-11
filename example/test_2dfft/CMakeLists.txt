# Comman line agruments to cmake:
#  * -DCMAKE_BUILD_TYPE=Debug instead of set(CMAKE_BUILD_TYPE Debug)
#  * -D_CRT_SECURE_NO_WARNINGS= : 
#
# See https://stackoverflow.com/questions/34438357/cmake-one-build-directory-for-multiple-projects
cmake_minimum_required(VERSION 3.00)
# -------- C/C++ preprocessor macros -----
#  You can use cmake command line to turn them on, e.g.,
#    -DDEBUG_LEVEL_INFO_LOW=ON 
#
if(MSVC)  # MSVC is a built-in macro. It is ON if the generator is MSVC
	#in Visual Studio to prevent compiling errors related to use
	#unsafe functions like sprintf
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)   

if(DEBUG_LEVEL_INFO_MEDIUM)
	add_definitions(-DDEBUG_LEVEL_INFO_MEDIUM)
endif(DEBUG_LEVEL_INFO_MEDIUM)	

if(DEBUG_LEVEL_INFO_LOW)
	add_definitions(-DDEBUG_LEVEL_INFO_LOW)
endif(DEBUG_LEVEL_INFO_LOW)	

#----------------------------------------------
project(2dfft)

if( NOT MSVC )
	include_directories("D:/code/fftw-3.3.8/api")
	link_directories("D:/code/fftw-3.3.8/build_android")

	add_executable(2dfft
	               simple_example.c
	)
	target_link_libraries(2dfft fftw3)
endif()

add_executable(ourfft 
	ourfft.cpp
)

