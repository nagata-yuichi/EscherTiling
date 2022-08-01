/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_I__
#define __OPERATOR_I__

#ifndef  __OPERATOR_BASE__
#include "operator_base.h"
#endif

class TOperator_I : public TOperator_base {
 public:
  TOperator_I( int N, int N_in );
  ~TOperator_I(); 
  virtual void SetParameter(); // set default parameters
  virtual void SetInit( VectorXd w, VectorXd w_in, MatrixXi kki ); // set information of the goal mesh
  virtual void Cal_u( int st1 ); // optimize tile shape (E_I distance)
  virtual void Cal_u_Procrustes( int st1 ); // optimize tile shape (Procrustes distance version)
  virtual void Check_Eval( VectorXd uu,  double eval ); // checking the correctness (for degugging)
  void Set_Bdj_Bsj( int st1 );  // Set fBdj and fBsj 
  void Cal_L_p( VectorXd& p, VectorXd& p_d, MatrixXd& L, int iter ); 
  void Cal_L_gzai( VectorXd& gzai, MatrixXd Lt );
  double Cal_eval_gzai( VectorXd& gzai ); // Solve (Eq. 14) (E_I distance)
  virtual void SetK(); // Set fKK, fG, fH, fG2_inv
  void Check_Orthogonal_K( int st1 );
  bool Check_Intersection( VectorXd& u ); // intersection exsist? -> false:not exist, true:exist
  void Trans_UU_center( VectorXd& uu ); // Align the center of the tile with the center of the goal
  void Trans_UU_wrt_G( VectorXd& uu );  // update top tile shapes
  virtual void Update_UU_top( double eval, VectorXd& uu, int st1 ); // update top tile shapes (mesh)

  SparseMatrix<double> fG;   // (G_I)
  MatrixXd fBdj;             // dens part of (B_{ikj})
  SparseMatrix<double> fBsj; // sparse part of (B_{ikj})
  int fMk;                   // diminsion of the parameters (=fBdj.cols()+fBsj.cols())

  double fAlpha; // A parameter for the weight setting (alpha_i for the inner points)

  MatrixXd fKK;  // Laplacian matrix with weights of the goal mesh (K)
  MatrixXi fKKi; // Laplacian matrix or adjacency matrix 
  MatrixXd fH;   // (- G_2^{-1} * G_1)
  MatrixXd fG2_inv; // (G_2^{-1})

  int fN;         // the number of points of the goal polygon (n)
  int fN_in;      // the number of the inner points (n')
  int fNN;        // the number of all points (n+n')
  VectorXd fW;    // the coordinates of the goal polygon (w)
  VectorXd fW_d;  // (w_c)
  VectorXd fW_in; // the coordinates of the inner points of the goal mesh (w')
  VectorXd fWW;   // the coordinates of all points of the goal mesh (\tilde{w})
  VectorXd fWW_d; // (\tilde{w_c})

  double fEval_best; // the evaluation value of the current (fNum_top-th) best solution 

  MatrixXd fA; // Asymmetric Laplacian matrix with weights of the goal mesh 
  MatrixXd fDin;
  MatrixXd fDout;

  MatrixXi fNei;      // fNei(i,*) is a set of the neighbors of point i (N(i))
  VectorXi fNei_Size; // fNei_Size(i) is the number of the neighbors of point i 
  double fEval_wGw;   // fW.transpose()*fG*fW 


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

