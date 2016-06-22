cmake_minimum_required(VERSION 2.8)
project(ExamplesNumericalSolvers)
#
# required libraries: libeigen3, libgtest
#
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# check existence of c++11 compiler
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# Eigen3 Headerfiles
# find_package(Eigen3 REQUIRED)
include_directories(./include/eigen/)

# GTest headerfiles
# find_package(GTest REQUIRED)
# include_directories(SYSTEM ${GTEST_INCLUDE_DIRS})
include_directories(./include/gtest/googletest/include)
include_directories(./include/cppoptlib)


SET( EXAMPLE_FILES  L-BFGS-gp)





set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} ) #/bin/examples
foreach( currentfile ${EXAMPLE_FILES} )
    add_executable( ${currentfile} "${currentfile}.cpp"  )
    target_link_libraries(${currentfile} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
endforeach( currentfile ${EXAMPLE_FILES} )
