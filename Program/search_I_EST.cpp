/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_I_EST__
#include "search_I_EST.h"
#endif


TSearch_I_EST::TSearch_I_EST( int N, int N_in ) : CLASS_NAME ( N, N_in )  
{
  fN = N;
  fN_in = N_in;
}

				
TSearch_I_EST::~TSearch_I_EST()
{

}

void TSearch_I_EST::SetParameter()
{
  fNum_top = 20;    // this value may be reset in TEnvironment_E::SetInit() 
  fFlagDisplay = 1; // this value may be reset in TEnvironment_E::SetInit() 
}


void TSearch_I_EST::Define()
{
  fEval_top = new double [ fNum_top ];
  fU_top = new VectorXd [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    fU_top[i] = VectorXd::Zero(2*fN);
  fUU_top = new VectorXd [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    fUU_top[i] = VectorXd::Zero(2*fNN);
  fIndex_top = new int [ fNum_top ];
  fNv_top = new int [ fNum_top ];
  ffn_top = new int* [ fNum_top ];
  for( int i = 0; i < fNum_top; ++i )
    ffn_top[i] = new int [ 6 ];
  fIH_top = new int [ fNum_top ];
  fSt1_top = new int [ fNum_top ];
}


void TSearch_I_EST::DoIt()
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


void TSearch_I_EST::IH4() 
{
  this->Make_Bv_Be_IH4();
  
  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= fN-fNv-2*k1; k2 = k2+1 ){ 
      for( int k3 = 0; k3 <= fN-fNv-2*k1-k2; k3 = k3+1 ){ 
	if( k2+k3 > fN-fNv - (2*k1+k2+k3) ) // fn[1]+fn[2] > fN[4]+fN[5] 
	  continue;

     	for( int k4 = 0; k4 <= fN-fNv-2*k1-k2-k3; k4 = k4+1 ){ 
	  this->Make_B_IH4( k1, k2, k3, k4 );
	  for( int j = 0; j < fN; ++j )
	    this->Cal_u(j); 
	}
      }
    }
  }
}


void TSearch_I_EST::IH5() 
{
  this->Make_Bv_Be_IH5();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      for( int k3 = 0; k3 <= fN-fNv-2*k1-2*k2; k3 = k3+1 ){ 
	this->Make_B_IH5( k1, k2, k3 );
	for( int j = 0; j < fN; ++j )
	  this->Cal_u_Procrustes(j);
      }
    }
  }
}


void TSearch_I_EST::IH6() 
{
  this->Make_Bv_Be_IH6();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      for( int k3 = 0; k3 <= fN-fNv-2*k1-2*k2; k3 = k3+1 ){ 
	this->Make_B_IH6( k1, k2, k3 );
	for( int j = 0; j < fN; ++j )
	  this->Cal_u_Procrustes(j);
      }
    }
  }
}


void TSearch_I_EST::IH1() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH1();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH1( k1, k2 );
      for( int j = 0; j < fN; ++j )
	this->Cal_u(j);
    }
  }
}


void TSearch_I_EST::IH2() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH2();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH2( k1, k2 );
      for( int j = 0; j < fN; ++j )
	this->Cal_u_Procrustes(j);
    }
  }
}


void TSearch_I_EST::IH3() 
{
  if( fN % 2 == 1 )
    return;
  
  this->Make_Bv_Be_IH3();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH3( k1, k2 );
      for( int j = 0; j < fN; ++j )
	this->Cal_u_Procrustes(j);
    }
  }
}

void TSearch_I_EST::IH7() 
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
      for( int j = 0; j < fN; ++j )
	this->Cal_u(j);
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }
}


void TSearch_I_EST::IH21() 
{
  int sign;

  sign = -1;

 AAA:
  this->Make_Bv_Be_IH21( sign );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH21( k1, k2, sign );
      for( int j = 0; j < fN; ++j )
	this->Cal_u(j);
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }
}


void TSearch_I_EST::IH28() 
{
  int sign;
  
  sign = -1;

 AAA:
  this->Make_Bv_Be_IH28( sign );

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      this->Make_B_IH28( k1, k2, sign );
      for( int j = 0; j < fN; ++j )
	this->Cal_u(j);
    }
  }

  if( sign == - 1){
    sign = 1;
    goto AAA;
  }
}


void TSearch_I_EST::Update_UU_top( double eval, VectorXd& uu, int st1 )
{
  VectorXd u(2*fN);
  VectorXd u0(2*fN);

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

  u = uu.head(2*fN);
  for( int i = 0; i < fN; ++i ){ 
    u0(ix(i)) = u(ix(i+st1));
    u0(iy(i)) = u(iy(i+st1));
  }
  fU_top[fIndex_top[orderIns]] = u0;
  fUU_top[fIndex_top[orderIns]] = uu;   

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
    tDisplay->Tile_mesh( uu );
}


void TSearch_I_EST::Set_candi( int num_candi, int* IH_candi, int* Nv_candi, int* st1_candi, int** fn_candi )
{
}
