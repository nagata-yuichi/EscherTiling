/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_E_CONF__
#define __SEARCH_E_CONF__

#define AD   // Set E or AD
// E: Eucliean distance, AD: AD distance

#ifdef E
  #ifndef  __OPERATOR_E__
  #include "operator_E.h"
  #define CLASS_NAME TOperator_E
  #endif
#endif
#ifdef AD
  #ifndef  __OPERATOR_AD__
  #include "operator_AD.h"
  #define CLASS_NAME TOperator_AD
  #endif
#endif

#ifndef __DISPLAY__
#include "display.h"
#endif


class TSearch_E_conf : public CLASS_NAME{ 
 public:
  TSearch_E_conf( int N );
  ~TSearch_E_conf(); 

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
  
  void Update_U_top( double eval, VectorXd& u, int st1 ) override; // modified to address large values of fNum_top     
  void QuickIndex( double* Arg, int* indexOrderd, int begin, int end ); // Returns the sorted index (indexOrdered) of array Arg in ascending order, where the arguments begin and end arethe indecies of the start and end of Arg.
  
  int fFlagDisplay;  
  TDisplay *tDisplay; 

  int fNum_top; 
  double *fEval_top; 
  VectorXd *fU_top;  
  int *fIndex_top;   
  int **ffn_top;  
  int *fNv_top; 
  int *fIH_top; 
  int *fSt1_top;
};


#endif
