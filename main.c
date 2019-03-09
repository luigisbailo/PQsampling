// author luigisbailo


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>


#include "parameters.h"
#include "greensFunct.h"
#include "draw.h"
#include "tools.h"
#include "init.h"
#include "step.h"
#include "shell.h"
#include "print.h"
#include "burst.h"
#include "checks.h"
#include "toolsAnnih.h"
#include "BD/bruteForce.h"
#include "BD/run_BD.h"
#include "BD/run_annih_BD.h"
#include "Fig1/fig1.h"
#include "Fig2/run_hybGF_P.h"
#include "Fig2/run_hybGF_PQ.h"
#include "Fig2/run_annih_P.h"
#include "Fig2/run_annih_PQ.h"
#include "Fig2/fig2_diff.h"
#include "Fig2/fig2_annih.h"
#include "Fig4/run_hybGF_P_proj.h"
#include "Fig4/run_hybGF_PQ_proj.h"
#include "Fig4/run_annih_PQ_proj.h"
#include "Fig4/run_annih_P_proj.h"
#include "Fig4/fig4_diff_proj.h"
#include "Fig4/fig4_annih_proj.h"


int main () {

    printf("\n\nFig1 is being produced:\n\n");

    fig1();

    printf("\n\nFig2-MSD is being produced:\n\n");

    fig2_diff();

    printf("\n\nFig2-kinetiks is being produced:\n\n");

    fig2_annih();

    printf("\n\nFig4-MSD is being produced:\n\n");

    fig4_diff_proj();

    printf("\n\nFig4-kinetiks is being produced:\n\n");

    fig4_annih_proj();


}