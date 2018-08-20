//TO COMPILE: g++ -std=c++11 main.cpp -o main -lgsl -lgslcblas -lm

#define print(x) cout << x << endl;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <sstream>

using namespace std::chrono;

#include "parameters.hpp"
#include "tools.hpp"
#include "greensFunct.hpp"
#include "draw.hpp"
#include "init.hpp"
#include "step.hpp"
#include "shell.hpp"
#include "print.hpp"
#include "burst.hpp"
#include "bruteForce.hpp"
#include "checks.hpp"

int main () {

// cout << drawTimeNewt (1.5541445463895798e-08,0.001,0.38271502428688109) << endl;
    std::cout << drawPosNewt (0.1,0.2,0.01,0.2) << std::endl;

}







