#CMake list file for compiling mesytec DAQ data processing software 
#library
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(NE255HW5)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11")
set(HEADERS 
      ${PROJECT_SOURCE_DIR}/hw5.h
    )

set(SOURCES 
      ${PROJECT_SOURCE_DIR}/hw5.cpp
    )
include_directories(${PROJECT_SOURCE_DIR}/include 
                    )
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/linkdef)
set(LINKDEF ${PROJECT_BINARY_DIR}/linkdef/LinkDef.h)
set(LD_HEADER
"#ifdef __CINT__

#pragma link off all globals\;
#pragma link off all classes\;
#pragma link off all functions\;
#pragma link C++ nestedclasses\;\n")
file(WRITE ${LINKDEF} ${LD_HEADER})
foreach(HEADER ${HEADERS})
  file(APPEND ${LINKDEF} "#pragma link C++ defined_in " ${HEADER}\;\n)
endforeach()
file(APPEND ${LINKDEF} "\n#endif\n")


add_library(NE255HW5 SHARED ${SOURCES})

target_link_libraries(NE255HW5)

target_include_directories(NE255HW5 
                           PUBLIC ${PROJECT_SOURCE_DIR}/include
                          )
add_executable(hw5_answer.exe hw5_answer.cpp)
target_link_libraries(hw5_answer.exe NE255HW5)
target_include_directories(hw5_answer.exe PUBLIC ${PROJECT_SOURCE_DIR}/include)