#CMake list file for compiling mesytec DAQ data processing software 
#library
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(NE255FINAL)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++11")
set(HEADERS 
      ${PROJECT_SOURCE_DIR}/include/readXSec.h
      ${PROJECT_SOURCE_DIR}/include/randNum.h
      ${PROJECT_SOURCE_DIR}/include/dataTypes.h
    )

set(SOURCES 
      ${PROJECT_SOURCE_DIR}/src/readXSec.cpp
      ${PROJECT_SOURCE_DIR}/src/dataTypes.cpp
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


add_library(NE255FINAL SHARED ${SOURCES})

target_link_libraries(NE255FINAL)

target_include_directories(NE255FINAL
                           PUBLIC ${PROJECT_SOURCE_DIR}/include
                          )
add_executable(eventGenerator.exe eventGenerator.cpp)
target_link_libraries(eventGenerator.exe NE255FINAL)
target_include_directories(eventGenerator.exe PUBLIC ${PROJECT_SOURCE_DIR}/include)