#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>


#include "parameters.hpp"
#include "greensFunct.hpp"
#include "draw.hpp"
#include "tools.hpp"
#include "init.hpp"
#include "step.hpp"
#include "shell.hpp"
#include "print.hpp"
#include "burst.hpp"
#include "checks.hpp"
#include "toolsAnnih.hpp"

#include "BD/bruteForce.hpp"
#include "BD/run_BD.hpp"
#include "BD/run_annih_BD.hpp"
#include "Fig1/fig1.hpp"
#include "Fig2/run_hybGF_P.hpp"
#include "Fig2/run_hybGF_PQ.hpp"
#include "Fig2/run_annih_PQ.hpp"
#include "Fig2/run_annih_P.hpp"
#include "Fig2/run_annih_PQ.hpp"
#include "Fig2/fig2_diff.hpp"
#include "Fig2/fig2_annih.hpp"
#include "Fig4/run_hybGF_P_proj.hpp"
#include "Fig4/run_hybGF_PQ_proj.hpp"
#include "Fig4/run_annih_PQ_proj.hpp"
#include "Fig4/run_annih_P_proj.hpp"
#include "Fig4/run_annih_PQ_proj.hpp"
#include "Fig4/fig4_diff_proj.hpp"
#include "Fig4/fig4_annih_proj.hpp"


int main () {

    std::cout << std::endl<< std::endl << "Fig1 is being produced:" << std::endl<< std::endl;

    fig1();

    std::cout << std::endl<< std::endl << "Fig2-MSD is being produced:" << std::endl<< std::endl;

    fig2_diff();

    std::cout << std::endl<< std::endl << "Fig2-kinetics is being produced:" << std::endl<< std::endl;

    fig2_annih();

    std::cout << std::endl<< std::endl << "Fig4-MSD is being produced:" << std::endl<< std::endl;

    fig4_diff_proj();

    std::cout << std::endl<< std::endl << "Fig4-kinetics is being produced:" << std::endl<< std::endl;

    fig4_annih_proj();


}