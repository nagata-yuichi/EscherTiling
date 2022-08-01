/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT_E__
#include "env_E.h"     
#endif

TEnvironment_E::TEnvironment_E()
{
}


TEnvironment_E::~TEnvironment_E()
{
}


void TEnvironment_E::SetInit()
{
  this->SetTarget();

#ifdef EST
  tSearch = new TSearch_E_EST( fN );
#endif
#ifdef CONF
  tSearch = new TSearch_E_conf( fN );
#endif
  
  tSearch->SetParameter();
  tSearch->fFlagDisplay = fFlagDisplay;
  tSearch->fNum_top = fNum_top;

  tSearch->Define();
  tSearch->SetInit( fW );

  tDisplay = new TDisplay();
  tDisplay->SetInit( fW );
  tSearch->tDisplay = tDisplay; 
}


void TEnvironment_E::DoIt()
{
  this->SetInit();
  if( fFlagDisplay ){
    tDisplay->Goal();
  }

  fTimeStart = clock();  
  tSearch->DoIt();
  fTimeEnd = clock();  

  printf( "time = %.2f \n", (fTimeEnd-fTimeStart) / (double)CLOCKS_PER_SEC );

#ifdef EST
  this->Write_Tile();
  if( fFlagDisplay )
    this->Display_top();
#endif
#ifdef CONF
  this->Write_Conf_data();
#endif
}


void TEnvironment_E::SetTarget()
{
  FILE* fp = fopen( fInstanceName, "r" );
  if( fp == NULL ){
    printf( "file not found\n");
    exit(0);
  }

  int dumy = fscanf( fp, "%d", &fN );
  double xd, yd;

  fW = VectorXd(2*fN);
 
  for( int i = 0; i < fN; ++i ){
    dumy = fscanf( fp, "%lf%lf", &xd, &yd );
    fW(2*i) = xd;
    fW(2*i+1) = yd;
  }
  dumy = fscanf( fp, "%lf", &xd );
  assert( feof( fp ) != 0 );

  fclose( fp );
}


void TEnvironment_E::Display_top()
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
    tDisplay->Tile( tSearch->fU_top[index] );

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


void TEnvironment_E::Write_Tile()
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
  fprintf( fp, "%d %d\n", fN, 0 ); 
  for( int i = 0; i < fN; ++i )
    fprintf( fp, "%lf %lf\n", fW(2*i), fW(2*i+1) ); 

  // tile shapes
  for( int s = 0; s < tSearch->fNum_top; ++s ){
    index = tSearch->fIndex_top[s];
    fprintf( fp, "%d %d %d %d %lf\n", tSearch->fN, 0, tSearch->fIH_top[index], tSearch->fNv_top[index], tSearch->fEval_top[index] );
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      fprintf( fp, "%d ", tSearch->ffn_top[index][i] );
    fprintf( fp, "%d ", tSearch->fSt1_top[index] );
    fprintf( fp, "\n" );

    for( int i = 0; i < fN; ++i ){ 
      fprintf( fp, "%lf %lf\n", (double)tSearch->fU_top[index](2*i), (double)tSearch->fU_top[index](2*i+1) );
    }
  } 
  fclose(fp);
}


void TEnvironment_E::Write_Conf_data()
{
  FILE *fp;
  char filename[80];
  int index;

  if( fOutputName == NULL )
    return;

  sprintf( filename, "%s.conf", fOutputName );
  fp = fopen( filename, "w" );

  // time
  fprintf( fp, "%.2f \n", (fTimeEnd-fTimeStart) / (double)CLOCKS_PER_SEC );

  // conf data
  for( int s = 0; s < tSearch->fNum_top; ++s ){
    index = tSearch->fIndex_top[s]; 
    fprintf( fp, "%lf %d %d\n", tSearch->fEval_top[index], tSearch->fIH_top[index], tSearch->fNv_top[index] );
    for( int i = 0; i < tSearch->fNv_top[index]; ++i )
      fprintf( fp, "%d ", tSearch->ffn_top[index][i] );
    fprintf( fp, "%d ", tSearch->fSt1_top[index] );
    fprintf( fp, "\n" );
  } 

  fclose(fp);
}
