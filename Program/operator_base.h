/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_BASE__
#define __OPERATOR_BASE__

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;


class TOperator_base {
 public:
  TOperator_base( int N );
  ~TOperator_base();

  MatrixXd Trans_AtoB( const MatrixXd& A ); // return matrix B from matrix A (Eqs.3-4)
  void Orthogonalization_Bd(); // orthonormalize the column vectors of matrix fBd 
  void Check_Orthogonal();     // check the orthonormality of matrices fBd and fBs 

 // set matrix A (Eqs. 2-3)
  void Set_As_Parallel( int ivs, int ive, MatrixXd& As, int& m ); // parallel J edges
  void Set_As_Deg( int iv, int deg, MatrixXd& As, int& m );       // consecutive J edges
  
  // set matrix fBe (before Eq. 22)
  void Set_Be_Jd_Parallel( int ivs, int ive );  // parallel J edges between the ivs and ive-th tiling vertices
  void Set_Be_S_Parallel( int ivs, int ive );   // parallel S edges ...
  void Set_Be_J( int ivs );                     // sigle J edge starting at the ivs-th tiling vertex
  void Set_Be_S( int ivs);                      // single S edge ...
  void Set_Be_Js_Xsymmetry( int ivs, int ive ); // J edges (symmetric w.r.t the x-axis) ...
  void Set_Be_Js_Ysymmetry( int ivs, int ive ); // J edges (symmetric w.r.t the y-axis) ...
  void Set_Be_Jd_Deg( int iv, int deg );        // consecutive J edges with a apecified angle at the iv-th tiling vertex

  // set matrix fBd (Eq.23)
  void Set_Bd_Jd_Parallel( int ivs, int ive );  
  void Set_Bd_S_Parallel( int ivs, int ive );   
  void Set_Bd_S( int ivs );                     
  void Set_Bd_Js_Xsymmetry( int ivs, int ive );
  void Set_Bd_Js_Ysymmetry( int ivs, int ive );
  void Set_Bd_Jd_Deg( int iv, int deg );

  void Make_Bv_Be_IH4(); // Set fNv, fBe, fMd;  Resize fBd, fBs;
  void Make_B_IH4( int k1, int k2, int k3, int k4 ); // Set fBd, fBs, fMs

  void Make_Bv_Be_IH5();
  void Make_B_IH5( int k1, int k2, int k3 );

  void Make_Bv_Be_IH6();
  void Make_B_IH6( int k1, int k2, int k3 );

  void Make_Bv_Be_IH1();
  void Make_B_IH1( int k1, int k2 );

  void Make_Bv_Be_IH2();
  void Make_B_IH2( int k1, int k2 );

  void Make_Bv_Be_IH3();
  void Make_B_IH3( int k1, int k2 );

  void Make_Bv_Be_IH7( int sign ); // sign -> -1:unclockwise, 1:clockwise
  void Make_B_IH7( int k1, int k2, int sign );

  void Make_Bv_Be_IH21( int sign );
  void Make_B_IH21( int k1, int k2, int sign );

  void Make_Bv_Be_IH28( int sign );
  void Make_B_IH28( int k1, int k2, int sign );

  int fN;        // the number of points of the goal polygon (n)
  int fNv;       // the number of tiling vertices of a template (n_v)
  int fIH;       // isohedral type
  MatrixXd fBd;  // dens part of the matrix B (B''_d in Eq. 23)
  int fMd;       // = fBd.cols() (m_d)
  SparseMatrix<double> fBs; //  sparse part of the matrix B (B_s in Eq. 24)
  int fMs;       // = fBs.cols() (m_s)

  MatrixXd fBv;  // (B_v in Eq. 16)
  int fMv;       // = fBv.cols() 
  MatrixXd fBe;  // The row elements store {d^b_s}^T and {d^f_s}^T (before Eq. 22)
  vector<T> fTripletList;  

  int *fi; // fi[v] (v=0, ..., fNv) index of the point of the tile polygon corresponting to the v-th tiling vertex of a template
  int *fn; // fn[v] (v=0, ..., fNv-1) the numbers of points on the tiling edges of a template (k_1, k_2, ... )

  
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

  int inline sx_b( int iv ){ 
    return 4*iv;
  }
  int inline sy_b( int iv ){
    return 4*iv+1;
  }
  int inline sx_f( int iv ){
    return 4*iv+2;
  }
  int inline sy_f( int iv ){
    return 4*iv+3;
  }

};

#endif


// Note:
// The equation numbers corresponds to those in [1].
// The coordinates of the goal (tile) shape are defined as (x0, y0, x1, y1, ...), unlike in [1]. 
// [1] An Efficient Exhaustive Search Algorithm for the Escherization Problem, Algorithmica

