cmake_minimum_required(VERSION 3.9)
project(PQsampling)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(PQsampling
    Fig1/main.cpp
    BD/bruteForce.hpp
    burst.hpp
    checks.hpp
    draw.hpp
    greensFunct.hpp
    init.hpp
    parameters.hpp
    print.hpp
    BD/run_BD.hpp
    shell.hpp
    step.hpp
    tools.hpp toolsAnnih.hpp)


find_package(GSL REQUIRED)
target_link_libraries(PQsampling GSL::gsl GSL::gslcblas)