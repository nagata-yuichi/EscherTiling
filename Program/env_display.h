/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#define __ENVIRONMENT__

#ifndef __DISPLAY__
#include "display.h"
#endif

#include "Eigen/Core"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>


using namespace Eigen;
using namespace std;

class TEnvironment {
 public:
  TEnvironment(); 
  ~TEnvironment(); 
  void SetInit();   // Initial setup
  void ReadTile();  // Read tile shapes from the specified file
  void DoIt();      // Main procedure
  void Display_top(); // Display tile shpaes and tiling results
  void SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd ); // Sort data

  
  int fN, fNv;        // n, n_v
  int fN_in, fNN;     // n', n+n'
  VectorXd fW;        // w
  VectorXd fW_d;      // w_c
  VectorXd fWW;       // \tilde{w}
  TDisplay *tDisplay; // // a class of drawing the tile shapes
  char *fResultFileName; // input file name 

  int fNum_top;       // the number of top solutions stored in the EST
  double *fEval_top;  // the evaluating values of the top solutions 
  VectorXd *fU_top;   // the top solutions (tile shape)
  VectorXd *fUU_top;  // the top solutions (mesh representation) 
  int *fIndex_top;    // [i] -> index of the i-th best solution
  int **ffn_top;      // [s] -> fi[] of the tile shape fU_top[s]
  int *fNv_top;       // [s] -> fNv of the tile shape fU_top[s]
  int *fIH_top;       // [s] -> IH type of the tile shape fU_top[s]
  int *fSt1_top;      // [s] -> st1 of the tile shape fU_top[s]

  double fTime;       // execution time 
  MatrixXi fKKi;      // Laplacian matrix of the goal mesh

  int inline ix( int i ){ // positoin of the x-coordinate of the i-th point of the goal polygon
    if ( i >= fN )
      i -= fN;
    if ( i < 0 )
      i += fN;
    return 2*i;
  }
  int inline iy( int i ){ // positoin of the y-coordinate of the i-th point of the goal polygon
    if ( i >= fN )
      i -= fN;
    if ( i < 0 )
      i += fN;
    return 2*i+1;
  }
};

#endif

  
