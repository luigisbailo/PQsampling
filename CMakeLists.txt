cmake_minimum_required(VERSION 3.9)
project(paper_formalization_first_passage)

set(CMAKE_C_STANDARD 11)

add_executable(main
        main.c
        burst.hpp
        checks.hpp
        draw.hpp
        greensFunct.hpp
        init.hpp
        parameters.hpp
        print.hpp
        shell.hpp
        step.hpp
        tools.hpp
        toolsAnnih.hpp
        BD/bruteForce.hpp
        BD/run_annih_BD.hpp
        BD/run_BD.hpp
        Fig1/fig1.hpp
        Fig2/fig2_diff.hpp
        Fig2/fig2_annih.hpp
        Fig2/run_annih_P.hpp
        Fig2/run_annih_PQ.hpp
        Fig2/run_hybGF_P.hpp
        Fig2/run_annih_PQ.hpp
        Fig4/fig4_annih_proj.hpp
        Fig4/fig4_diff_proj.hpp
        Fig4/run_annih_P_proj.hpp
        Fig4/run_annih_PQ_proj.hpp
        Fig4/run_hybGF_P_proj.hpp
        Fig4/run_annih_PQ_proj.hpp
        Fig4/run_hybGF_PQ_proj.hpp
        )

find_package(GSL REQUIRED)
target_link_libraries(main GSL::gsl GSL::gslcblas)
