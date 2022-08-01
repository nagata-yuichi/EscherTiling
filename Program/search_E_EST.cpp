/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_E_EST__
#include "search_E_EST.h"
#endif


TSearch_E_EST::TSearch_E_EST( int N ) : CLASS_NAME ( N ) 
{
  fN = N;
}

				
TSearch_E_EST::~TSearch_E_EST()
{

}

void TSearch_E_EST::SetParameter()
{
  fNum_top = 20;    // this value may be reset in TEnvironment_E::SetInit() 
  fFlagDisplay = 1; // this value may be reset in TEnvironment_E::SetInit() 
  fFlagCheckIntersection = 1; 
  fFlagPrune = 1;
}


void TSearch_E_EST::Define()
{
  fEval_top = new double [ fNum_top ];
  fU_top = new VectorXd [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    fU_top[i] = VectorXd::Zero(2*fN);
  fIndex_top = new int [ fNum_top ];
  fNv_top = new int [ fNum_top ];
  ffn_top = new int* [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    ffn_top[i] = new int [ 6 ];
  fIH_top = new int [ fNum_top ];
  fSt1_top = new int [ fNum_top ];
}


void TSearch_E_EST::DoIt()
{
  fEval_best = 999999999.9;
  for( int i = 0; i < fNum_top; ++i ){
    fEval_top[i] = 999999999.9;
    fIndex_top[i] = i;
  }

  printf( "IH4\n");
  this->IH4();
  printf( "IH5\n");
  this->IH5();  
  printf( "IH6\n");
  this->IH6();
  printf( "IH1\n");
  this->IH1();  
  printf( "IH2\n");
  this->IH2();   
  printf( "IH3\n");
  this->IH3();   
  printf( "IH7\n");
  this->IH7();  
  printf( "IH21\n");
  this->IH21();   
  printf( "IH28\n");
  this->IH28();
}


void TSearch_E_EST::IH4() 
{
  if( fFlagPrune ){ // procedure for pruning the search
    this->Search_Comp_SS();
    this->Make_Bv_Be_IH4_Part();
  }

  this->Make_Bv_Be_IH4(); 

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= fN-fNv-2*k1; k2 = k2+1 ){ 
      for( int k3 = 0; k3 <= fN-fNv-2*k1-k2; k3 = k3+1 ){ 
	if( k2+k3 > fN-fNv - (2*k1+k2+k3) ) // k2+k3 > k4+k5 
	  continue; 

	if( fFlagPrune ){ // procedure for pruning the search
	  this->Make_B_IH4_Part( k1, k2, k3 );
	  this->Cal_u_Part();
	}

     	for( int k4 = 0; k4 <= fN-fNv-2*k1-k2-k3; k4 = k4+1 ){ 
	  if( !fFlagPrune ) 
	    this->Make_B_IH4( k1, k2, k3, k4 ); 
	  this->Cal_u( k1, k2, k3, k4 ); 
	}
      }
    }
  }
}


void TSearch_E_EST::IH5() 
{
  if( fFlagPrune ){ // procedure for pruning the search
    this->Search_Comp_SS();
    this->Make_Bv_Be_IH5_Part();
  }

  this->Make_Bv_Be_IH5();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      if( fFlagPrune ) { // procedure for pruning the search
	this->Make_B_IH5_Part( k1, k2 );
	this->Cal_u_Procrustes_Part();
      }
      for( int k3 = 0; k3 <= fN-fNv-2*k1-2*k2; k3 = k3+1 ){ 
	if( !fFlagPrune ) 
	  this->Make_B_IH5( k1, k2, k3 );
	this->Cal_u_Procrustes( k1, k2, k3, 0 );
      }
    }
  }
}


void TSearch_E_EST::IH6() 
{
  if( fFlagPrune ){ // procedure for pruning the search
    this->Search_Comp_SJS();
    this->Make_Bv_Be_IH6_Part();
  }
  
  this->Make_Bv_Be_IH6();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      if( fFlagPrune ){ // procedure for pruning the search
	this->Make_B_IH6_Part( k1, k2 );
	this->Cal_u_Procrustes_Part();
      }
      for( int k3 = 0; k3 <= fN-fNv-2*k1-2*k2; k3 = k3+1 ){ 
	if( !fFlagPrune )
	  this->Make_B_IH6( k1, k2, k3 );
	this->Cal_u_Procrustes( k1, k2, k3, 0 );
      }
    }
  }
}



void TSearch_E_EST::IH1() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH1();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH1( k1, k2 );
      this->Cal_u( k1, k2, 0, 0 );
    }
  }
}

void TSearch_E_EST::IH2() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH2();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH2( k1, k2 );
      this->Cal_u_Procrustes( k1, k2, 0, 0 );
    }
  }
}

void TSearch_E_EST::IH3() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH3();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH3( k1, k2 );
      this->Cal_u_Procrustes( k1, k2, 0, 0 );
    }
  }
}

void TSearch_E_EST::IH7() 
{
  if( fN % 2 == 1 )
    return;

  int sign;
  sign = -1;

 AAA:
  this->Make_Bv_Be_IH7( sign );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH7( k1, k2, sign );
      this->Cal_u( k1, k2, 0, 0 );
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }
}

void TSearch_E_EST::IH21() 
{
  int sign;
  sign = -1;

 AAA:
  this->Make_Bv_Be_IH21( sign );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH21( k1, k2, sign );
      this->Cal_u( k1, k2, 0, 0 );
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }

  this->Make_Bv_Be_IH21( 1 );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH21( k1, k2, 1 );
      this->Cal_u( k1, k2, 0, 0 );
    }
  }
}

void TSearch_E_EST::IH28() 
{
  int sign;
  sign = -1;

 AAA:
  this->Make_Bv_Be_IH28( sign );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH28( k1, k2, sign );
      this->Cal_u( k1, k2, 0, 0 );
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }
}


void TSearch_E_EST::Update_U_top( double eval, VectorXd& u, int st1 )
{
  static int count_imp = 0;

  for( int s = 0; s < fNum_top; ++s ){
    if( fabs( eval - fEval_top[fIndex_top[s]] ) < 0.0000001 ) 
      return;
  }
  
  int orderIns;
  for( int s = 0; s < fNum_top; ++s ){
    if( eval < fEval_top[fIndex_top[s]] ){ 
      orderIns = s;
      break;
    }
  }

  int index_stock = fIndex_top[fNum_top-1];
  for( int s = fNum_top-1; s >= orderIns+1; --s )
    fIndex_top[s] = fIndex_top[s-1];
  
  fIndex_top[orderIns] = index_stock;
  fEval_top[fIndex_top[orderIns]] = eval;
  fU_top[fIndex_top[orderIns]] = u;
  fNv_top[fIndex_top[orderIns]] = fNv;
  for( int h = 0; h < fNv; ++h )
    ffn_top[fIndex_top[orderIns]][h] = fn[h];
  fEval_best = fEval_top[fIndex_top[fNum_top-1]]; 
  fIH_top[fIndex_top[orderIns]] = fIH;
  fSt1_top[fIndex_top[orderIns]] = st1;
  
  printf( "IH%d(%d): eval=%lf ", fIH, count_imp, eval );
  printf("(");
  for( int i = 0; i < fNv; ++i )
    printf( "%d ", fn[i] );
  printf(")");
  printf( " %d", st1 );
  printf( "\n" );
  ++count_imp;

  if( fFlagDisplay )
    tDisplay->Tile( u );
}

