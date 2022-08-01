/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_E_EST__
#define __SEARCH_E_EST__

#define E   // Set E or AD
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


class TSearch_E_EST : public CLASS_NAME{ 
 public:
  TSearch_E_EST( int N );
  ~TSearch_E_EST(); 

  void SetParameter(); // set parameter values
  void Define();       // declare arryas 

  void DoIt();         // main process of the EST (exhaustive search of the templates)
  void IH4();          // EST for IH4
  void IH5();
  void IH6();
  void IH1();
  void IH2();
  void IH3();
  void IH7();
  void IH21();
  void IH28();
  
  void Update_U_top( double eval, VectorXd& u, int st1 ) override; // update top solutions
  
  int fFlagDisplay;   // display results? 0:no, 1:yes
  TDisplay *tDisplay; // a class used to display results

  int fNum_top;      // the number of top solutions stored in the EST
  double *fEval_top; // the evaluation values of the top solutions
  VectorXd *fU_top;  // the top solutions (tile shapes)
  int *fIndex_top;   // [i] -> index of the i-th best solution
  int **ffn_top;     // [s] -> fi[] of the tile shape fU_top[s]
  int *fNv_top;      // [s] -> fNv of the tile shape fU_top[s]
  int *fIH_top;      // [s] -> IH type of the tile shape fU_top[s]
  int *fSt1_top;     // [s] -> st1 of the tile shape fU_top[s]
};


#endif
