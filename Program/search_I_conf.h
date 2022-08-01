/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_CONF__
#define __SEARCH_I_CONF__

#define I   // Set I or IR
// I: E_I distance, IR: E_IR distance

#ifdef I
  #ifndef  __OPERATOR_I__
  #include "operator_I.h"
  #define CLASS_NAME TOperator_I
  #endif
#endif
#ifdef IR
  #ifndef  __OPERATOR_IR__
  #include "operator_IR.h"
  #define CLASS_NAME TOperator_IR
  #endif
#endif

#ifndef __DISPLAY__
#include "display.h"
#endif


class TSearch_I_conf : public CLASS_NAME{  
 public:
  TSearch_I_conf( int N, int N_in );
  ~TSearch_I_conf(); 

  void SetParameter();
  void Define();

  void DoIt();
  void IH4();
  void IH5();
  void IH6();
  void IH1();
  void IH2();
  void IH3();
  void IH7();
  void IH21();
  void IH28();
  
  void Update_UU_top( double eval, VectorXd& uu, int st1 ) override; 
  void Set_candi( int num_candi, int* IH_candi, int* Nv_candi, int* st1_candi, int** fn_candi );
  
  int fFlagDisplay;  
  TDisplay *tDisplay; 

  int fNum_top; 
  double *fEval_top; 
  VectorXd *fU_top;  
  VectorXd *fUU_top; 
  int *fIndex_top;   
  int **ffn_top;  
  int *fNv_top;  
  int *fIH_top; 
  int *fSt1_top;

  int fNum_candi;     // the number of template configurations read from the file
  int *fIH_candi;     // [s] -> IH type of the s-th configuration
  int *fSt1_candi;    // [s] -> st1 (=j) of the the s-th configuration
  int **fki_candi;    // [s] -> (k_1, k_2, ...) of the s-th configuration
  int fCurrent_candi; // the number of the template configuration that is currently being optimized
};


#endif
