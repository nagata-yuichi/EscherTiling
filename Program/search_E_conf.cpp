/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __SEARCH_E_CONF__
#include "search_E_conf.h"
#endif


TSearch_E_conf::TSearch_E_conf( int N ) : CLASS_NAME ( N ) 
{
  fN = N;
}

				
TSearch_E_conf::~TSearch_E_conf()
{

}

void TSearch_E_conf::SetParameter()
{
  fNum_top = 50000; // this value may be reset in TEnvironment_E::SetInit() 
  fFlagDisplay = 0; // this value may be reset in TEnvironment_E::SetInit() 
  fFlagCheckIntersection = 0; // intersections must be allowed
  fFlagPrune = 1;
}


void TSearch_E_conf::Define()
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

void TSearch_E_conf::DoIt()
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

  this->Update_U_top( 99999999.9, fU_top[fNum_top-1], -1 ); // needed here
} 


void TSearch_E_conf::IH4() 
{
  if( fFlagPrune ){ 
    this->Search_Comp_SS();
    this->Make_Bv_Be_IH4_Part();
  }

  this->Make_Bv_Be_IH4();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= fN-fNv-2*k1; k2 = k2+1 ){ 
      for( int k3 = 0; k3 <= fN-fNv-2*k1-k2; k3 = k3+1 ){ 
	if( k2+k3 > fN-fNv - (2*k1+k2+k3) ) 
	  continue;

	if( fFlagPrune ){
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


void TSearch_E_conf::IH5() 
{
  if( fFlagPrune ){
    this->Search_Comp_SS();
    this->Make_Bv_Be_IH5_Part();
  }

  this->Make_Bv_Be_IH5();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      if( fFlagPrune ) {
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


void TSearch_E_conf::IH6() 
{
  if( fFlagPrune ){ 
    this->Search_Comp_SJS();
    this->Make_Bv_Be_IH6_Part();
  }
  
  this->Make_Bv_Be_IH6();

  for( int k1 = 0; k1 <= (fN-fNv)/2; k1 = k1+1 ){ 
    for( int k2 = 0; k2 <= (fN-fNv-2*k1)/2; k2 = k2+1 ){ 
      if( fFlagPrune ){ 
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


void TSearch_E_conf::IH1() 
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

void TSearch_E_conf::IH2() 
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

void TSearch_E_conf::IH3() 
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

void TSearch_E_conf::IH7() 
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

void TSearch_E_conf::IH21() 
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

void TSearch_E_conf::IH28() 
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


void TSearch_E_conf::Update_U_top( double eval, VectorXd& u, int st1 )
{
  static int count_tmp = 0;
  static double eval_tmp[100000];
  static int nv_tmp[100000];
  static int st1_tmp[100000];
  static int IH_tmp[100000];
  static int fn_tmp[100000][6];
  // u is not stored 

  assert( fNum_top <= 100000 );

  if( st1 != -1 ){ 
    eval_tmp[count_tmp] = eval;
    nv_tmp[count_tmp] = fNv;
    for( int h = 0; h < fNv; ++h )
      fn_tmp[count_tmp][h] = fn[h];
    st1_tmp[count_tmp] = st1;
    IH_tmp[count_tmp] = fIH;
    ++count_tmp;
  }

  if( count_tmp == 0 ){
    assert( st1 == -1 );
    return;
  }

  if( count_tmp == fNum_top || st1 == -1 ){ 
    int index_tmp[fNum_top]; 
    int index_new[fNum_top];

    for( int i = 0; i < count_tmp; ++i )  
      index_tmp[i] = i;    
    this->QuickIndex( eval_tmp, index_tmp, 0, count_tmp-1 ); 


    double e_top, e_tmp;    
    int i_top, i_tmp;
    int orderIns_tmp = 0;
    int orderIns_top = 0;
    int orderIns_new = 0;
    while( 1 ){

      assert( orderIns_top < fNum_top-orderIns_tmp );
      i_top = fIndex_top[orderIns_top];
      e_top = fEval_top[i_top];
      if( orderIns_tmp < count_tmp ){ 
	i_tmp = index_tmp[orderIns_tmp];
	e_tmp = eval_tmp[i_tmp];
      }
      else{
	e_tmp = e_top + 1.0;
      }


      if( e_top <= e_tmp ){
	index_new[orderIns_new] = i_top;
	++orderIns_top;
      }
      else {
	int i_new = fIndex_top[fNum_top-1-orderIns_tmp];
	index_new[orderIns_new] = i_new;
	fEval_top[i_new] = e_tmp;
	fNv_top[i_new] = nv_tmp[i_tmp];
	for( int h = 0; h < fNv; ++h )
	  ffn_top[i_new][h] = fn_tmp[i_tmp][h];
	fSt1_top[i_new] = st1_tmp[i_tmp];
	fIH_top[i_new] = IH_tmp[i_tmp];
	++orderIns_tmp;
      }

      ++orderIns_new;
      if( orderIns_new == fNum_top )
	break;
    }

    for( int i = 0; i < fNum_top; ++i )
      fIndex_top[i] = index_new[i];

    count_tmp = 0;
    fEval_best = fEval_top[fIndex_top[fNum_top-1]]; 

    for( int i = 0; i < fNum_top-1; ++i )
      assert( fEval_top[fIndex_top[i]] <=fEval_top[fIndex_top[i+1]] );
  }

}

void TSearch_E_conf::QuickIndex( double* Arg, int* indexOrderd, int begin, int end )
{
  int i, j, m;
  double pivot;
  int stock; 
  int flag; 

  flag = 0;

  for( m = (begin + end)/2; m < end; ++m ){
    if( Arg[ indexOrderd[m] ] != Arg[ indexOrderd[m+1] ] ){
      if( Arg[ indexOrderd[m] ] > Arg[ indexOrderd[m+1] ] )
	pivot = Arg[ indexOrderd[m] ];
      else
	pivot = Arg[ indexOrderd[m+1] ];
      flag = 1;
      break;
    }
  }
  if( flag == 0 ){
    for( m = (begin + end)/2; m > begin; --m ){
      if( Arg[ indexOrderd[m] ] != Arg[ indexOrderd[m-1] ] ){
	if( Arg[ indexOrderd[m] ] > Arg[ indexOrderd[m-1] ] )
	  pivot = Arg[ indexOrderd[m] ];
	else
	  pivot = Arg[ indexOrderd[m-1] ];
	flag = 1;
	break;
      }
    }
  }

  if( flag == 0 )
    return;


  i = begin;
  j = end;

  while(1)
  {

    while(1){
      if( Arg[ indexOrderd[i] ] >= pivot )
	break;
      ++i;
    }
    while(1){
      if( Arg[ indexOrderd[j] ] < pivot )
	break;
      --j;
    }

    if( i >= j ) 
      break;

    stock = indexOrderd[i];
    indexOrderd[i] = indexOrderd[j];
    indexOrderd[j] = stock;
  }

  if( begin < i-1 )
    this->QuickIndex( Arg, indexOrderd, begin, i-1 );
  if( j+1 < end )
    this->QuickIndex( Arg, indexOrderd, j+1, end );
}
