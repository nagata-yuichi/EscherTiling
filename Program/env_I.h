/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#define __ENVIRONMENT_I__

#define EST   // Set EST or CONF
// EST: exhaustive search of the templates
// CONF: heuristic search 

#ifdef EST
  #ifndef __SEARCH_I_EST__ 
  #include "search_I_EST.h"
  #endif
#endif
#ifdef CONF
  #ifndef __SEARCH_I_CONF__ 
  #include "search_I_conf.h"
  #endif
#endif


class TEnvironment_I {
 public:
  TEnvironment_I(); 
  ~TEnvironment_I(); 
  void SetTarget();          // set the goal figure
  virtual void SetInit();    // initial setting
  virtual void DoIt();       // main procedure
  void Display_top();        // display top tile shapes 
  void Write_Tile();         // write top tilig shapes to a file
  void ReadConfigurations(); // read configuration data 

#ifdef EST
  TSearch_I_EST *tSearch;
#endif
#ifdef CONF
  TSearch_I_conf *tSearch;
#endif
  
  int fN;  // see operation_I.h
  int fN_in; 
  int fNN;
  VectorXd fW; 
  VectorXd fW_in; 
  VectorXd fWW; 
  MatrixXi fKKi;

  clock_t fTimeStart, fTimeEnd; // for execution tme
  char *fInstanceName;          // name of the target figure
  char *fOutputName = NULL;     // output file name
  char *fFileNameConfiguration; // file name of the configuration data

  TDisplay *tDisplay; // a class of drawing the tile shapes
  int fFlagDisplay;   // 0: not display, 1: display tile figures
  double fAlpha;      // parameter alpha for the E_I and E_IR distances
  int fNum_top;        // the number of top solutions stored
};

#endif

  
