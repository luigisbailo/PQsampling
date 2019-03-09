cmake_minimum_required(VERSION 3.9)
project(paper_formalization_first_passage)

set(CMAKE_C_STANDARD 11)

add_executable(main
        main.c
        burst.h
        checks.h
        draw.h
        greensFunct.h
        init.h
        parameters.h
        print.h
        shell.h
        step.h
        tools.h
        toolsAnnih.h
        BD/bruteForce.h
        BD/run_annih_BD.h
        BD/run_BD.h
        Fig1/fig1.h
        Fig2/fig2_diff.h
        Fig2/fig2_annih.h
        Fig2/run_annih_P.h
        Fig2/run_annih_PQ.h
        Fig2/run_hybGF_P.h
        Fig2/run_annih_PQ.h
        Fig4/fig4_annih_proj.h
        Fig4/fig4_diff_proj.h
        Fig4/run_annih_P_proj.h
        Fig4/run_annih_PQ_proj.h
        Fig4/run_hybGF_P_proj.h
        Fig4/run_annih_PQ_proj.h
        Fig4/run_hybGF_PQ_proj.h
        )

find_package(GSL REQUIRED)
target_link_libraries(main GSL::gsl GSL::gslcblas)
