#ifndef __SEARCH_I_CONF__
#include "search_I_conf.h"
#endif


/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

TSearch_I_conf::TSearch_I_conf( int N, int N_in ) : CLASS_NAME ( N, N_in ) 
{
  fN = N;
  fN_in = N_in;
}

				
TSearch_I_conf::~TSearch_I_conf()
{

}

void TSearch_I_conf::SetParameter()
{
  fNum_top = 20;    // this value may be reset in TEnvironment_E::SetInit() 
  fFlagDisplay = 1; // this value may be reset in TEnvironment_E::SetInit() 
}


void TSearch_I_conf::Define()
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

  fIH_candi = new int [ 100000 ];
  fSt1_candi = new int [ 100000 ];
  fki_candi = new int* [ 100000 ];
  for( int i = 0; i < 100000; ++i )
    fki_candi[i] = new int [ 4 ];
}


void TSearch_I_conf::DoIt()
{
  fEval_best = 999999999.9;
  for( int i = 0; i < fNum_top; ++i ){
    fEval_top[i] = 999999999.9;
    fIndex_top[i] = i;
  }

  fCurrent_candi = 0; 

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


void TSearch_I_conf::IH4() 
{
  int k1, k2, k3, k4, st1;
  this->Make_Bv_Be_IH4();
  
  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 4 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    k3 = fki_candi[fCurrent_candi][2];
    k4 = fki_candi[fCurrent_candi][3];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH4( k1, k2, k3, k4 ); 
    this->Cal_u( st1 ); 
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH5() 
{
  int k1, k2, k3, st1;

  this->Make_Bv_Be_IH5();

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 5 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    k3 = fki_candi[fCurrent_candi][2];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH5( k1, k2, k3 );
    this->Cal_u_Procrustes( st1 );
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH6() 
{
  int k1, k2, k3, st1;

  this->Make_Bv_Be_IH6();

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 6 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    k3 = fki_candi[fCurrent_candi][2];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH6( k1, k2, k3 );
    this->Cal_u_Procrustes( st1 );
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH1() 
{
  if( fN %2 == 1 )
    return;

  int k1, k2, st1;

  this->Make_Bv_Be_IH1();

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 1 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH1( k1, k2 );
    this->Cal_u( st1 );
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH2() 
{
  if( fN %2 == 1 )
    return;

  int k1, k2, st1;

  this->Make_Bv_Be_IH2();

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 2 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH2( k1, k2 );
    this->Cal_u_Procrustes( st1 );
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH3() 
{
  if( fN %2 == 1 )
    return;

  int k1, k2, st1;

  this->Make_Bv_Be_IH3();

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 3 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH3( k1, k2 );
    this->Cal_u_Procrustes( st1 );
    ++fCurrent_candi;
  }
}

void TSearch_I_conf::IH7() 
{
  if( fN %2 == 1 )
    return;

  int k1, k2, st1;
  int sign;
  int stock_current_candi;

  stock_current_candi = fCurrent_candi;

  sign = -1;

 AAA:
  this->Make_Bv_Be_IH7( sign ); 

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 7 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH7( k1, k2, sign );
    this->Cal_u( st1 );
    ++fCurrent_candi;
  }

  if( sign == - 1){
    fCurrent_candi = stock_current_candi;
    sign = 1;
    goto AAA;
  }
}

void TSearch_I_conf::IH21() 
{
  int k1, k2, st1;
  int sign;
  int stock_current_candi;
  
  sign = -1;

 AAA:
  this->Make_Bv_Be_IH21( sign ); 
  
  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 21 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH21( k1, k2, sign );
    this->Cal_u( st1 );
    ++fCurrent_candi;
  }

  if( sign == - 1){
    fCurrent_candi = stock_current_candi;
    sign = 1;
    goto AAA;
  }
}

void TSearch_I_conf::IH28() 
{
  int k1, k2, st1;
  int sign;
  int stock_current_candi;

  sign = -1;

 AAA:
  this->Make_Bv_Be_IH28( sign ); 

  while(1){
    if( fCurrent_candi == fNum_candi )
      break;
    if( fIH_candi[fCurrent_candi] != 28 )
      break;
    k1 = fki_candi[fCurrent_candi][0]; 
    k2 = fki_candi[fCurrent_candi][1];
    st1 = fSt1_candi[fCurrent_candi];
    this->Make_B_IH28( k1, k2, sign );
    this->Cal_u( st1 );
    ++fCurrent_candi;
  }

  if( sign == - 1){
    fCurrent_candi = stock_current_candi;
    sign = 1;
    goto AAA;
  }
}


void TSearch_I_conf::Update_UU_top( double eval, VectorXd& uu, int st1 )
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


void TSearch_I_conf::Set_candi( int num_candi, int* IH_candi, int* Nv_candi, int* st1_candi, int** fn_candi )
{
  int checked;
  int index;

  checked = 0;
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 4 ){
      fIH_candi[checked] = 4;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][2];
      fki_candi[checked][3] = fn_candi[s][4];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 5 ){
      fIH_candi[checked] = 5;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][4];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 6 ){
      fIH_candi[checked] = 6;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fki_candi[checked][2] = fn_candi[s][3];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 1 ){
      fIH_candi[checked] = 1;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 2 ){
      fIH_candi[checked] = 2;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 3 ){
      fIH_candi[checked] = 3;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][1];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 7 ){
      fIH_candi[checked] = 7;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][2];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 21 ){
      fIH_candi[checked] = 21;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][2];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  for( int s = 0; s < num_candi; ++s ){
    if( IH_candi[s] == 28 ){
      fIH_candi[checked] = 28;
      fki_candi[checked][0] = fn_candi[s][0];
      fki_candi[checked][1] = fn_candi[s][2];
      fSt1_candi[checked] = st1_candi[s];
      ++checked;
    }
  }
  assert( checked == num_candi );
  fNum_candi = num_candi;
}
