/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_E__
#define __OPERATOR_E__

#ifndef  __OPERATOR_BASE__
#include "operator_base.h"
#endif

class TOperator_E : public TOperator_base {
 public:
  TOperator_E( int N );
  ~TOperator_E();

  virtual void SetParameter();        // set default parameters
  virtual void SetInit( VectorXd w ); // set information of the goal polygon

  virtual void Cal_u( int k1, int k2, int k3, int k4 ); // optimize tile shapes (Euclidean distance)
  virtual void Cal_u_Procrustes( int k1, int k2, int k3, int k4 ); // optimize tile shapes (Procrustes distance)

  virtual void Update_U_top( double eval, VectorXd& u, int st1 ); // update top solutions
  void Trans_U_center( VectorXd& u ); // Align the center of u with the center of w
  bool Check_Intersection( VectorXd& u ); // self intersection exsist? -> false:not exist, true:exist

  VectorXd fW;     // the coordinates of the goal polygon (w in Eq. 26)
  VectorXd fW_d;   // (w_c in Eq. 26)
  VectorXd *fWj;   // re-indexed fW (w_j)
  VectorXd *fWj_d; // re-indexed fW_d ({w_c}_j)

  double fEval_best; // the evaluation value of the current (fNum_top-th) best solution 
  int fFlagPrune;   // 0: naive EST, 1: efficient EST
  int fFlagDisplay; // display results?  0:no, 1:yes 
  int fFlagCheckIntersection; // allow intersection in the tile shape? 0:allowed, 1:not allowed
  

  // The following functions and variables are used for the decomposition relaxation (Section 4)
  // use: fNp 
  // common: fNv, fBd, fBs, fBv, fBe
  void Search_Comp_SS();
  void Make_Bv_Be_Comp_SS();
  virtual void Make_B_Comp_SS( int k1, int k2 ); // IH
  virtual void Cal_u_Comp_SS( int k1, int k2 );
  void Set_Bd_S_Comp( int ivs, int col_s );

  void Search_Comp_SJS(); // IH
  void Make_Bv_Be_Comp_S(); 
  virtual void Make_B_Comp_S( int k1 ); // IH
  virtual void Cal_u_Comp_Sfb( int k1 );

  int fNp; // the number of points of a partial template
  double ***fEval_Comp_best;
  double **fEval_Comp_best_Sf;
  double **fEval_Comp_best_Sb;

  // Part:
  // use: fNp, fNvp, fBd_part, fBs_part, fBv_part, fBe_part   
  virtual void Cal_u_Part();
  virtual void Cal_u_Procrustes_Part();
  void Set_Be_Jd_Parallel_Part( int ivs1, int ivs2 );
  void Set_Be_S_Part( int ivs ); 
  void Set_Be_Js_Xsymmetry_Part( int ivs1, int ivs2 );
  void Set_Be_J_Part( int ivs );
  void Set_Bd_Jd_Parallel_Part( int ivs1, int ivs2 );
  void Set_Bd_S_Part( int ivs );
  void Set_Bd_Js_Xsymmetry_Part( int ivs1, int ivs2 ); 
  void Set_Bd_J_Part( int ivs );
  void Make_Bv_Be_IH4_Part();
  void Make_B_IH4_Part( int k1, int k2, int k3 );
  void Make_Bv_Be_IH5_Part();
  void Make_B_IH5_Part( int k1, int k2 );
  void Make_Bv_Be_IH6_Part();
  void Make_B_IH6_Part( int k1, int k2 );
  void Orthogonalization_Bd_Part();

  int fNvp; // the number of tiling vertices of a partial template 
  MatrixXd fBv_part; 
  MatrixXd fBe_part; 
  MatrixXd fBd_part;
  SparseMatrix<double> fBs_part; 
  int fMd_part;
  int fMs_part;
  double *fEval_part; 
};


#endif


// Note:
// The equation numbers corresponds to those in [1].
// [1] An Efficient Exhaustive Search Algorithm for the Escherization Problem, Algorithmica
