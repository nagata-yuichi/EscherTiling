/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_E__
#include "env_E.h"
#endif

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
  TEnvironment_E* tEnv = NULL;
  tEnv = new TEnvironment_E();

  tEnv->fInstanceName = argv[1]; // file name of the goal figure
  tEnv->fOutputName = argv[2];   // output file name 
  sscanf( argv[3], "%d", &(tEnv->fFlagDisplay) ); // display results? 0:no, 1:yes
  sscanf( argv[4], "%d", &(tEnv->fNum_top) );     // the number of top solutions stored
  tEnv->DoIt();
  
  return 0;
}
