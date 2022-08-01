/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_IR__
#define __OPERATOR_IR__

#ifndef  __OPERATOR_I__
#include "operator_I.h"
#endif

class TOperator_IR : public TOperator_I {
 public:
  TOperator_IR( int N, int N_in );
  ~TOperator_IR(); 
  void SetParameter() override;
  virtual void SetInit( VectorXd w,  VectorXd w_in, MatrixXi kki ) override ;
  virtual void Cal_u( int st1 ) override;
  virtual void Cal_u_Procrustes( int st1 ) override; 
  void Check_Eval( VectorXd uu,  double eval ) override;
  void Set_RotationMatrices( VectorXd uu );  
  double Cal_eval_gzai_Rotation( VectorXd& gzai, int iter ); 
  void SetK() override;

  
  VectorXd fR_cos;    
  VectorXd fR_sin; 
  MatrixXd *fRot; 
  MatrixXd *fRotD;

  MatrixXd fGG;
  double fEval_wGGw;
  SparseMatrix<double> fAA;
  VectorXd fTTWW;
  VectorXd fTTWW_d;
};

#endif

