#CMake list file for compiling mesytec DAQ data processing software 
#library
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(NE255HW3)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11")
set(HEADERS 
      ${PROJECT_SOURCE_DIR}/ne255_hw3.h
    )

set(SOURCES 
      ${PROJECT_SOURCE_DIR}/ne255_hw3.cpp
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


add_library(NE255HW3 SHARED ${SOURCES})

target_link_libraries(NE255HW3)

target_include_directories(NE255HW3 
                           PUBLIC ${PROJECT_SOURCE_DIR}/include
                          )
add_executable(hw3_answer.exe hw3_answer.cpp)
target_link_libraries(hw3_answer.exe )
target_include_directories(hw3_answer.exe PUBLIC ${PROJECT_SOURCE_DIR}/include)