/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env_I.h"
#endif

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
  TEnvironment_I* tEnv = NULL;
  tEnv = new TEnvironment_I();

  tEnv->fInstanceName = argv[1]; // file name of the target shape 
  tEnv->fOutputName = argv[2];   // output file name 
  sscanf( argv[3], "%d", &(tEnv->fFlagDisplay) ); // display results? 0:no, 1:yes
  sscanf( argv[4], "%d", &(tEnv->fNum_top) );     // the number of top solutions stored
  sscanf( argv[5], "%lf", &(tEnv->fAlpha) );      // the value of alpha 
  tEnv->fFileNameConfiguration = argv[6]; // file name of the configuration data
  tEnv->DoIt();
  
  return 0;
}
