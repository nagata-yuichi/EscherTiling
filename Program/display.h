/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#define __DISPLAY__

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include  "xw.h"
#include  <X11/Xlib.h>
#include <unistd.h>



using namespace Eigen;
using namespace std;

class TDisplay {
 public:
  TDisplay(); 
  ~TDisplay();
  void SetInit( VectorXd w ); // Set data of the goal polygon 
  void SetInit( int N, int N_in, VectorXd ww, MatrixXi kki ); // Set data of the goal mesh
  void Set_range(); // Get the maximum and minimum values of the coordinates of the goal shape
  void Goal();      // Display the goal polygon
  void Goal_mesh(); // Display the goal mesh
  void Tile( VectorXd& u ); // Display the tile polygon
  void Tile_mesh( VectorXd& uu ); // Display the tile mesh
  void Tiling_IH4( VectorXd& u, int *fn, int nv ); // Display the tiling result 
  void Tiling_IH5( VectorXd& u, int *fn, int nv ); // Display the tiling result 
  void Tiling_IH6( VectorXd& u, int *fn, int nv ); // Display the tiling result
  void Tiling_IH1( VectorXd& u, int *fn, int nv ); // Display the tiling result 
  void All_tiles( VectorXd* u );                   // Display all the top tile shapes 
  
  int fNumDisplay;
  int fOpen_flag_tiling;  
  double fWidth, fHeight;
  int fxs, fxe, fys, fye;
  int fxs_t, fxe_t, fys_t, fye_t;
  int fN;
  int fN_in;
  int fNN;
  VectorXd fW;
  VectorXd fW_in;
  VectorXd fWW;
  MatrixXi fKKi;
  
  int inline ix( int i ){ // positoin of the x-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i;
  }

  int inline iy( int i ){ // positoin of the y-coordinate of the i-th point of the goal polygon
    if( i < 0 )
      i += fN;
    else if ( i >= fN )
      i -= fN;
    return 2*i+1;
  }

  int inline iix( int i ){ // positoin of the x-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN;
    else if ( i >= fNN )
      i -= fNN;
    return 2*i;
  }

  int inline iiy( int i ){ // positoin of the y-coordinate of the i-th point of the goal mesh
    if( i < 0 )
      i += fNN;
    else if ( i >= fNN )
      i -= fNN;
    return 2*i+1;
  }
};

#endif
