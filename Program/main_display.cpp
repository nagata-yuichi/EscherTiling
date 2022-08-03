/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env_display.h"
#endif

#include <stdio.h>
#include <stdlib.h>


int main( int argc, char* argv[] )
{
  TEnvironment* tEnv = NULL;
  tEnv = new TEnvironment();

  tEnv->fResultFileName = argv[1];
  tEnv->DoIt();
  
  return 0;
}
