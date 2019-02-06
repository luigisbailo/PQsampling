cmake_minimum_required(VERSION 3.9)
project(PQsampling)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(PQsampling
#    Fig1/main.cpp
#    Fig2/main.cpp
#    Fig2/main_BM.cpp
#    Fig2/main_GF.cpp
#    Fig4/main_GFproj.cpp
#    main.cpp
#    Annihilation/main.cpp
        Annihilation_proj/main_annih_proj.cpp
    Brute_force/bruteForce.hpp
    burst.hpp
    checks.hpp
    draw.hpp
        greensFunct.hpp
        init.hpp
    parameters.hpp
        print.hpp
        Brute_force/run_BM.hpp
        run_hybGF_P.hpp
    run_hybGF_P_proj.hpp
    run_hybGF_PQ.hpp
    run_hybGF_PQ_proj.hpp
    shell.hpp
    step.hpp
    tools.hpp tools_annih.h)


find_package(GSL REQUIRED)
target_link_libraries(PQsampling GSL::gsl GSL::gslcblas)