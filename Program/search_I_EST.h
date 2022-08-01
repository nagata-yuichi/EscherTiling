/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_EST__
#define __SEARCH_I_EST__

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


class TSearch_I_EST : public CLASS_NAME{  
 public:
  TSearch_I_EST( int N, int N_in );
  ~TSearch_I_EST(); 

  void SetParameter();  // set parameter values
  void Define();        // declare arryas 

  void DoIt();          // main process of the EST (exhaustive search of the templates)
  void IH4();           // EST of IH4
  void IH5();
  void IH6();
  void IH1();
  void IH2();
  void IH3();
  void IH7();
  void IH21();
  void IH28();
  
  void Update_UU_top( double eval, VectorXd& uu, int st1 ) override; // updating the top solutions
  void Set_candi( int num_candi, int* IH_candi, int* Nv_candi, int* st1_candi, int** fn_candi ); // set the promising template configurations

  
  int fFlagDisplay;   // display results? 0:no, 1:yes
  TDisplay *tDisplay; // a class used to display results

  int fNum_top;      // the number of top solutions stored in the EST
  double *fEval_top; // the evaluating values of the top solutions 
  VectorXd *fU_top;  // the top solutions (tile shape)
  VectorXd *fUU_top; // the top solutions (mesh representation) 
  int *fIndex_top;   // [i] -> index of the i-th best solution
  int **ffn_top;     // [s] -> fi[] of the tile shape fU_top[s]
  int *fNv_top;      // [s] -> fNv of the tile shape fU_top[s]
  int *fIH_top;      // [s] -> IH type of the tile shape fU_top[s]
  int *fSt1_top;     // [s] -> st1 of the tile shape fU_top[s]
};


#endif
