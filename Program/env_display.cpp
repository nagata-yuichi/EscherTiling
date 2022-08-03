/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __ENVIRONMENT__
#include "env_display.h"
#endif


TEnvironment::TEnvironment()
{
  
}

TEnvironment::~TEnvironment()
{

}


void TEnvironment::ReadTile()
{
  FILE* fp = fopen( fResultFileName, "r" );
  double x, y;
  int n,n_in, ih, nv;
  double eval;
  int dumy; 
  int st1;

  dumy = fscanf( fp, "%lf", &fTime );
  dumy = fscanf( fp, "%d%d", &fN, &fN_in );
  fNN = fN+fN_in;

  fW = VectorXd(2*fN);
  fW_d = VectorXd(2*fN);
  fWW = VectorXd(2*fNN);  
  fKKi = MatrixXi(fNN,fNN); 

  int num_top = 100;
  fEval_top = new double [ num_top ];
  fU_top = new VectorXd [ num_top ];
  for( int i = 0; i < num_top; ++i )
    fU_top[i] = VectorXd::Zero(2*fN);
  fUU_top = new VectorXd [ num_top ]; 
  for( int i = 0; i < num_top; ++i )
    fUU_top[i] = VectorXd::Zero(2*(fN+fN_in));
  fIndex_top = new int [ num_top ];
  fNv_top = new int [ num_top ];
  ffn_top = new int* [ num_top ];
  for( int i = 0; i < num_top; ++i )
    ffn_top[i] = new int [ 6 ];
  fIH_top = new int [ num_top ];
  fSt1_top = new int [ num_top ];

  int count = 0;
  for( int i = 0; i < fNN; ++i ){
    dumy = fscanf( fp, "%lf%lf", &x, &y );
    fWW(count++) = x;
    fWW(count++) = y;
  }
  fW = fWW.head(2*fN);
  for( int k = 0; k < fN; ++k ){
    fW_d(iy(k))=fW(ix((k)));
    fW_d(ix(k))=-fW(iy((k)));
  }

  if( fN_in != 0 ){
    for( int i = 0; i < fNN; ++i ){
      for( int j = 0; j < fNN; ++j ){
	int a;
	dumy = fscanf( fp, "%d", &a );
	fKKi(i,j) = a;
      }
    }
  }

  fNum_top = 0;
  while( 1 )
  {
    dumy = fscanf( fp, "%d%d%d%d%lf", &n,&n_in,&ih,&nv,&eval );
    if( feof( fp ) != 0 )
      break;
    assert( n == fN );
    assert( n_in == fN_in );
    fIndex_top[ fNum_top ] = fNum_top;
    fEval_top[ fNum_top ] = eval;
    fNv_top[ fNum_top ] = nv;
    fIH_top[ fNum_top ] = ih;
    for( int i = 0; i < nv; ++i )
      dumy = fscanf( fp, "%d", &ffn_top[ fNum_top ][i] );
    dumy = fscanf( fp, "%d", &fSt1_top[ fNum_top ] );

    for( int i = 0; i < fN; ++i ) {
      dumy = fscanf( fp, "%lf%lf", &x, &y );
      fU_top[ fNum_top ](2*i) = x;
      fU_top[ fNum_top ](2*i+1) = y;
    }

    if( fN_in != 0 ){ // mesh 
      for( int i = 0; i < fNN; ++i ) { 
	dumy = fscanf( fp, "%lf%lf", &x, &y );
	fUU_top[ fNum_top ](2*i) = x;
	fUU_top[ fNum_top ](2*i+1) = y;
      }
    }
    ++fNum_top; 
  }
  printf( "num_top = %d\n", fNum_top );
  fclose( fp );

  this->SortIndex( fEval_top, fNum_top, fIndex_top, fNum_top );
}


void TEnvironment::SetInit()
{
  tDisplay = new TDisplay();
  tDisplay->SetInit( fN, fN_in, fWW, fKKi );
}


void TEnvironment::DoIt()
{
  this->ReadTile();
  this->SetInit();

  if( fN_in == 0 )
    tDisplay->Goal();
  else
    tDisplay->Goal_mesh();


  Display_top();     // Display top solutions and tile results one by one
  // tDisplay->All_tiles( fU_top ); // Display top solutions simultaneously

  while(1){
    char ch;
    fprintf(stderr, "Hit Return key to continue.\n");
    ch = getc(stdin);
    if( ch == 'q' )
      break;
  }

}

void TEnvironment::Display_top()
{
  int order = 0;
  int index; 
  char ch;
  int mode;
  int fIH;

  mode = 1;
  while( 1 ){
    index = fIndex_top[order];
    printf( "order = %2d eval=%lf: IH=%2d", order, fEval_top[index], fIH_top[index] );
    printf(" (");
    for( int i = 0; i < fNv_top[index]; ++i )
      printf( "%d ", ffn_top[index][i] );
    printf(")");
    printf(" %d", fSt1_top[index] ) ;
    fIH = fIH_top[index];
    if( fN_in == 0 )
      tDisplay->Tile( fU_top[index] );
    else
      tDisplay->Tile_mesh( fUU_top[index] );

    if( fIH == 4 )
      tDisplay->Tiling_IH4( fU_top[index], ffn_top[index], fNv_top[index] );
    else if( fIH == 5 )
      tDisplay->Tiling_IH5( fU_top[index], ffn_top[index], fNv_top[index] );
    else if( fIH == 6 )
      tDisplay->Tiling_IH6( fU_top[index], ffn_top[index], fNv_top[index] );
    else if( fIH == 1 )
      tDisplay->Tiling_IH1( fU_top[index], ffn_top[index], fNv_top[index] );


    ch = getc(stdin);
    if( ch == 'f' )
      mode = 1;
    if( ch == 'b' )
      mode = -1;

    order += mode;

    if( order == fNum_top ) order = 0;
    if( order == -1 ) order = fNum_top-1;
  }
}

// Arg[]: input array
// numOfArg: the number of elements of Arg[]
// indexOrdered[i] <- the position of the i-th smallest element of Arg[]
// only the top numOfOrd elements are considered
void TEnvironment::SortIndex( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest;
  double valueBest;
  int checked[numOfArg];

  assert( Arg[0] < 99999999999.9 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = 99999999999.9;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] < valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }
}
