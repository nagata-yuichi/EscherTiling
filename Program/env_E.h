/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_E__
#define __ENVIRONMENT_E__

#define EST   // Set EST or CONF
// EST: exhaustive search of the templates
// CONF: get promising configurations for heuristic search

#ifdef EST
  #ifndef __SEARCH_E_EST__ 
  #include "search_E_EST.h"
  #endif
#endif
#ifdef CONF
  #ifndef __SEARCH_E_CONF__ 
  #include "search_E_conf.h"
  #endif
#endif


class TEnvironment_E {
 public:
  TEnvironment_E(); 
  ~TEnvironment_E(); 
  void SetTarget();        // set the target figure
  virtual void SetInit();  // initial setting
  virtual void DoIt();     // main procedure
  void Display_top();      // display top tile shapes 
  void Write_Tile();       // write top tile shapes to a file
  void Write_Conf_data();  // write top configurations to a file

#ifdef EST
  TSearch_E_EST *tSearch;
#endif
#ifdef CONF
  TSearch_E_conf *tSearch;
#endif
  
  int fN;      // the number of vertices of the goal (tile) polygon
  VectorXd fW; // the coordinates of the goal polygon

  clock_t fTimeStart, fTimeEnd; // execution time
  char *fInstanceName;          // name of the goal shape
  char *fOutputName = NULL;     // output file name

  TDisplay *tDisplay;  // a class of displaying tile shapes
  int fFlagDisplay;    // 0: not display, 1: display tile shapes
  int fNum_top;        // the number of top solutions stored
};

#endif

  
