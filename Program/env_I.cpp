/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_I__
#include "env_I.h"     
#endif

TEnvironment_I::TEnvironment_I()
{
}


TEnvironment_I::~TEnvironment_I()
{
}


void TEnvironment_I::SetInit()
{
  this->SetTarget();

#ifdef EST
  tSearch = new TSearch_I_EST( fN, fN_in );
#endif
#ifdef CONF
  tSearch = new TSearch_I_conf( fN, fN_in );
#endif
  
  tSearch->SetParameter();
  tSearch->fFlagDisplay = fFlagDisplay;
  tSearch->fNum_top = fNum_top;
  tSearch->fAlpha = fAlpha; 

  tSearch->Define();
  tSearch->SetInit( fW, fW_in, fKKi );

  tDisplay = new TDisplay();
  tDisplay->SetInit( fN, fN_in, fWW, fKKi );
  tSearch->tDisplay = tDisplay; 
}


void TEnvironment_I::DoIt()
{
  this->SetInit();
  if( fFlagDisplay ){
    tDisplay->Goal_mesh();
  }

#ifdef CONF
  this->ReadConfigurations();  
#endif

  fTimeStart = clock();  
  tSearch->DoIt();
  fTimeEnd = clock();  

  printf( "time = %.2f \n", (fTimeEnd-fTimeStart) / (double)CLOCKS_PER_SEC );

  this->Write_Tile();  
  if( fFlagDisplay )
    this->Display_top(); 
}


void TEnvironment_I::SetTarget()
{
  FILE* fp = fopen( fInstanceName, "r" );
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }

  int dumy = fscanf( fp, "%d%d", &fN, &fN_in );
  double xd, yd;

  fNN = fN+fN_in;
  fW = VectorXd(2*fN);
  fW_in = VectorXd(2*fN_in);
  fWW = VectorXd(2*fNN);
  fKKi = MatrixXi(fNN,fNN);
 
  for( int i = 0; i < fNN; ++i ){
    dumy = fscanf( fp, "%lf%lf", &xd, &yd );
    fWW(2*i) = xd;
    fWW(2*i+1) = yd;
  }

  fW = fWW.head(2*fN);
  fW_in = fWW.tail(2*fN_in);

  int d;
  for( int i = 0; i < fNN; ++i ){
    for( int j = 0; j < fNN; ++j ){
      assert( feof( fp ) == 0 );
      dumy = fscanf( fp, "%d", &d );
      fKKi(i,j) = d;
    }
  }
  dumy = fscanf( fp, "%d", &d );
  assert( feof( fp ) != 0 );

  fclose( fp );
}


void TEnvironment_I::Display_top()
{
  int order = 0;
  int index; 
  char ch;
  int mode;
  int fIH;

  mode = 1;
  while( 1 ){
    index = tSearch->fIndex_top[order];
    printf( "order = %2d eval=%lf: IH=%2d", order, tSearch->fEval_top[index], tSearch->fIH_top[index] );
    printf(" (");
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      printf( "%d ", tSearch->ffn_top[index][i] );
    printf(")");
    printf(" %d", tSearch->fSt1_top[index] ) ;
    fIH = tSearch->fIH_top[index];
    tDisplay->Tile_mesh( tSearch->fUU_top[index] );

    if( fIH == 4 )
      tDisplay->Tiling_IH4( tSearch->fU_top[index], tSearch->ffn_top[index], tSearch->fNv_top[index] );
    else if( fIH == 5 )
      tDisplay->Tiling_IH5( tSearch->fU_top[index], tSearch->ffn_top[index], tSearch->fNv_top[index] );
    else if( fIH == 6 )
      tDisplay->Tiling_IH6( tSearch->fU_top[index], tSearch->ffn_top[index], tSearch->fNv_top[index] );


    ch = getc(stdin);
    if( ch == 'f' )
      mode = 1;
    if( ch == 'b' )
      mode = -1;

    order += mode;

    if( order == tSearch->fNum_top ) order = 0;
    if( order == -1 ) order = tSearch->fNum_top-1;
  }
}


void TEnvironment_I::Write_Tile()
{
  FILE *fp;
  char filename[80];
  int index;
  
  for( int order = 0; order < tSearch->fNum_top; ++order ){
    index = tSearch->fIndex_top[order];
    printf( "order = %2d eval=%lf: IH=%2d", order, tSearch->fEval_top[index], tSearch->fIH_top[index] );
    printf(" (");
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      printf( "%d ", tSearch->ffn_top[index][i] );
    printf(") ");
    printf("%d", tSearch->fSt1_top[index] ) ;
    printf("\n");
  }

  if( fOutputName == NULL )
    return;

  sprintf( filename, "%s.tile", fOutputName );
  fp = fopen( filename, "w" );

  // time
  fprintf( fp, "%.2f \n", (fTimeEnd-fTimeStart) / (double)CLOCKS_PER_SEC );

  // instance
  fprintf( fp, "%d %d\n", fN, fN_in ); 
  for( int i = 0; i < fNN; ++i )
    fprintf( fp, "%lf %lf\n", fWW(2*i), fWW(2*i+1) ); 

  // KKi matrix   
  for( int i = 0; i < fNN; ++i ){
    for( int j = 0; j < fNN; ++j ){
      fprintf( fp, "%d ", fKKi(i,j) );
    }
    fprintf( fp, "\n" );
  }

  // tile shapes
  for( int s = 0; s < tSearch->fNum_top; ++s ){
    index = tSearch->fIndex_top[s];
    fprintf( fp, "%d %d %d %d %lf\n", tSearch->fN, tSearch->fN_in, tSearch->fIH_top[index], tSearch->fNv_top[index], tSearch->fEval_top[index] );
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      fprintf( fp, "%d ", tSearch->ffn_top[index][i] );
    fprintf( fp, "%d ", tSearch->fSt1_top[index] );
    fprintf( fp, "\n" );

    for( int i = 0; i < fN; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fU_top[index](2*i), (double)tSearch->fU_top[index](2*i+1) );
    }

    for( int i = 0; i < fN+fN_in; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fUU_top[index](2*i), (double)tSearch->fUU_top[index](2*i+1) );
    }
  } 
  fclose(fp);
}


void TEnvironment_I::ReadConfigurations() 
{
  int ih, nv, ni;
  double eval;
  double dd;
  int n, n_in, nn;
  int dumy; 
  int st1;
  FILE *fp;
  int ffi[7];
  double *eval_candi;
  int *IH_candi;
  int *Nv_candi;
  int *st1_candi;
  int **fn_candi;
  int num_candi;
  double time;
  

  eval_candi = new double [ 100000 ];
  IH_candi = new int [ 100000 ];
  Nv_candi = new int [ 100000 ];
  st1_candi = new int [ 100000 ];
  fn_candi = new int* [ 100000 ];
  for( int i = 0; i < 100000; ++i )
    fn_candi[i] = new int [ 6 ];

  
  fp = fopen( fFileNameConfiguration, "r" );
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }

  dumy = fscanf( fp, "%lf", &time );

  num_candi = 0;
  while( 1 )
  {
    dumy = fscanf( fp, "%lf%d%d", &eval, &ih, &nv );
    if( feof( fp ) != 0 )
      break;
    eval_candi[num_candi] = eval;
    IH_candi[num_candi] = ih;
    Nv_candi[num_candi] = nv;

    for( int i = 0; i < nv; ++i ) {
      dumy = fscanf( fp, "%d", &ni );
      fn_candi[ num_candi ][i] = ni;
    }
    dumy = fscanf( fp, "%d", &st1 );
      st1_candi[ num_candi ] = st1;

    ++num_candi; 
    if( num_candi == 100000 ){
      printf( "input data is too large !!\n" );
      exit(0);
    }
  }
  printf( "num_candi = %d\n", num_candi );
  fclose( fp );

  tSearch->Set_candi( num_candi, IH_candi, Nv_candi, st1_candi, fn_candi );
}
