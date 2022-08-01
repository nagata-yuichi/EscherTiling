/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_AD__
#define __OPERATOR_AD__

#ifndef  __OPERATOR_E__
#include "operator_E.h"
#endif

class TOperator_AD : public TOperator_E {
 public:
  TOperator_AD( int N );
  ~TOperator_AD(); 

  void SetInit( VectorXd w ) override;
  
  virtual void Cal_u( int k1, int k2, int k3, int k4 ) override;
  virtual void Cal_u_Procrustes( int k1, int k2, int k3, int k4 ) override;
  void Make_Ba();
  void Make_B_Comp_SS( int k1, int k2 ) override; 
  void Cal_u_Comp_SS( int k1, int k2 ) override;
  void Make_Ba_Comp();
  void Make_B_Comp_S( int k1 ) override; 
  void Cal_u_Comp_Sfb( int k1 ) override;
  void Cal_u_Part() override; 
  void Cal_u_Procrustes_Part() override; 
  void Make_Ba_Part();
  
  MatrixXd fBa; // the matrix B for the AD distance 
  int fMa;      // = fBa.cols()
  VectorXd *fWaj;   // the coordinates of the re-indexed goal polygon (for the AD distance)
  VectorXd *fWaj_d; // w_c versions 
};

#endif

  
