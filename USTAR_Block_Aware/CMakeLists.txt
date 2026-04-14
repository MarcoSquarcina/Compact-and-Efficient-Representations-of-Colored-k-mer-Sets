cmake_minimum_required(VERSION 3.16)
project(USTAR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)

# Usa O3 per ottimizzazione massima
add_compile_options(-O3 -W)

add_executable(ustar
    src/ustar.cpp
    src/commons.cpp
    src/DBG.cpp        
    src/SPSS.cpp       
    src/Sorter.cpp     
    src/Encoder.cpp    
    src/ColorGraph.cpp # Keep commented if truly not used
)