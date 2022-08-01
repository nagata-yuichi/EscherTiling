/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __DISPLAY__
#include "display.h"
#endif

int open_flag_tiling = 0;

TDisplay::TDisplay()
{
  fNumDisplay = 1;
  fWidth = 600;
  fHeight = 600;
}

TDisplay::~TDisplay()
{
}

void TDisplay::SetInit( VectorXd w )
{
  fN = w.rows()/2;
  fW = w;
  fN_in = 0;

  fOpen_flag_tiling = 0;
  this->Set_range();
}

void TDisplay::SetInit( int N, int N_in, VectorXd ww, MatrixXi kki )
{
  fN = N;
  fN_in = N_in;
  fNN = fN+fN_in;
  fWW = ww;
  fW = fWW.head(2*fN);
  fKKi = kki;

  fOpen_flag_tiling = 0;
  this->Set_range();
}

void TDisplay::Set_range()
{
  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;

  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9;

  for( int i = 0; i < fN; ++i ){
    if( fW(2*i) < x_min )
      x_min = fW(2*i);
    if( fW(2*i) > x_max )
      x_max = fW(2*i);
    if( fW(2*i+1) < y_min )
      y_min = fW(2*i+1);
    if( fW(2*i+1) > y_max )
      y_max = fW(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  fxs = (x_min+x_max)*0.5 - 0.6*length; 
  fxe = (x_min+x_max)*0.5 + 0.6*length; 
  fys = (y_min+y_max)*0.5 - 0.6*length; 
  fye = (y_min+y_max)*0.5 + 0.6*length;

  fxs_t = (x_min+x_max)*0.5 - 1.01*length; 
  fxe_t = (x_min+x_max)*0.5 + 1.01*length; 
  fys_t = (y_min+y_max)*0.5 - 1.01*length; 
  fye_t = (y_min+y_max)*0.5 + 1.01*length;
}

void TDisplay::Goal() 
{
  static int open_flg = 0;
  static int d_num = 0;

  /* open the window ------------------- */
  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;
    /*--- name the window --------------- */
    g_open((char*)"");
    /*--- set window number and size ----------*/
    g_open_window( d_num, fWidth, fHeight, (char*)"Goal" );
    /*--- define coordinate system in the window ----------*/
    g_window( d_num, fxs, fys, fxe, fye );
    /*--- set the character fonts used in the window */
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    /*--- set the bitmap pattern to be used when filling */
    // g_set_pixpat(d_num);
    usleep(100000);
  }
  g_clearWindow(d_num);

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN; ++i )
    g_line( d_num, fW(ix(i)), fW(iy(i)), fW(ix(i+1)), fW(iy(i+1)), 3 );

  g_flush();
}

void TDisplay::Goal_mesh() 
{
  assert( fN_in != 0 );
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Goal_mesh" );
    g_window( d_num, fxs, fys, fxe, fye );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);    
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN-1; ++i ){
    for( int j = i+1; j < fNN; ++j ){
      if( fKKi(i,j) != 0 ){
	g_line( d_num, fWW(iix(i)), fWW(iiy(i)), fWW(iix(j)), fWW(iiy(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN; ++i )
    g_line( d_num, fWW(ix(i)), fWW(iy(i)), fWW(ix(i+1)), fWW(iy(i+1)), 3 );


  g_flush();
}

void TDisplay::Tile( VectorXd& u ) 
{
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Tile" );
    g_window( d_num, fxs, fys, fxe, fye );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN; ++i )
    g_line( d_num, u(ix(i)), u(iy(i)), u(ix(i+1)), u(iy(i+1)), 3 );

  g_flush();
}

void TDisplay::Tile_mesh( VectorXd& uu ) 
{
  assert( fN_in != 0 );
  static int open_flg = 0;
  static int d_num = 0;

  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, fWidth, fHeight, (char*)"Tile mesh" );
    g_window( d_num, fxs, fys, fxe, fye );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  g_clearWindow(d_num);

  g_setColor(d_num,1); // black
  for( int i = 0; i < fNN-1; ++i ){
    for( int j = i+1; j < fNN; ++j ){
      if( fKKi(i,j) != 0 ){
	g_line( d_num, uu(iix(i)), uu(iiy(i)), uu(iix(j)), uu(iiy(j)), 2 );
      }
    }
  }

  g_setColor(d_num,2); // red
  for( int i = 0; i < fN; ++i )
    g_line( d_num, uu(ix(i)), uu(iy(i)), uu(ix(i+1)), uu(iy(i+1)), 3 );

  g_flush();
}


void TDisplay::Tiling_IH4( VectorXd& u, int* fn, int nv )
{
  int d_num;

  if ( fOpen_flag_tiling == 0 ) 
  {
    fOpen_flag_tiling = fNumDisplay;
    ++fNumDisplay;
    d_num = fOpen_flag_tiling;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, fxs_t, fys_t, fxe_t, fye_t );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  d_num = fOpen_flag_tiling;
  g_clearWindow(d_num);

  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*fN);
  VectorXd u_ori_1(2*fN);
  for( int i = 0; i < fN; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }
  u_ori_1 = - u_ori_0;

  diff(0) = u_ori_1(2*fi[3]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[3]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[fN][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 %2 == 0 ) g_setColor(d_num,10); 
      else g_setColor(d_num,15); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	pointsD[i][0] = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 %2 == 0 ) g_setColor(d_num,19); 
      else g_setColor(d_num,29); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );
    } 
  }

  g_flush();
}



// IH5
void TDisplay::Tiling_IH5( VectorXd& u, int* fn, int nv )
{
  int d_num;

  if ( fOpen_flag_tiling == 0 ) 
  {
    fOpen_flag_tiling = fNumDisplay;
    ++fNumDisplay;
    d_num = fOpen_flag_tiling;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, fxs_t, fys_t, fxe_t, fye_t );
    //g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  d_num = fOpen_flag_tiling;
  g_clearWindow(d_num);

  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;
    
  float cos_th, sin_th;
  VectorXd dir(2); 
  MatrixXd rot(2,2); 
  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  dir(0) = u(2*fi[4]) - u(2*fi[0]); 
  dir(1) = u(2*fi[4]+1) - u(fi[0]+1); 
  dir.normalize();
  cos_th = dir(0);
  sin_th = -dir(1);
  rot(0,0) = cos_th;
  rot(0,1) = -sin_th;
  rot(1,0) = sin_th;
  rot(1,1) = cos_th;

  VectorXd u_ori_0(2*fN);
  VectorXd u_ori_1(2*fN);
  VectorXd u_ori_2(2*fN);
  VectorXd u_ori_3(2*fN);
  for( int i = 0; i < fN; ++i ){
    u_ori_0(2*i) = cos_th*u(2*i) - sin_th*u(2*i+1);
    u_ori_0(2*i+1) = sin_th*u(2*i) + cos_th*u(2*i+1);
  }
  for( int i=0; i < fN; ++i ){ 
    u_ori_1(2*i) = u_ori_0(2*i);
    u_ori_1(2*i+1) = -u_ori_0(2*i+1);
    u_ori_2(2*i) = - u_ori_0(2*i);
    u_ori_2(2*i+1) = u_ori_0(2*i+1);
    u_ori_3(2*i) = - u_ori_0(2*i);
    u_ori_3(2*i+1) = - u_ori_0(2*i+1);
  }


  diff(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[3]); 
  diff(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[3]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_2(2*fi[5]) - u_ori_1(2*fi[0]); 
  diff(1) = u_ori_2(2*fi[5]+1) - u_ori_1(2*fi[0]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_3(2*fi[2]) - u_ori_2(2*fi[1]); 
  diff(1) = u_ori_3(2*fi[2]+1) - u_ori_2(2*fi[1]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_3(2*fi[5]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_3(2*fi[5]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[fN][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN ; ++i ){
	x = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	y = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      if( k1 % 2 == 0 ) g_setColor(d_num,10); 
      else  g_setColor(d_num,19); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	x = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	y = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      if( k1 % 2 == 0 ) g_setColor(d_num,15); 
      else  g_setColor(d_num,29); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	x = u_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	y = u_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      if( k1 % 2 == 0 ) g_setColor(d_num,19); 
      else  g_setColor(d_num,10); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	x = u_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	y = u_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      if( k1 % 2 == 0 ) g_setColor(d_num,29); 
      else  g_setColor(d_num,15); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );
    }
  }

  g_flush();
}


// IH6
void TDisplay::Tiling_IH6( VectorXd& u, int* fn, int nv )
{
  int d_num;

  if ( fOpen_flag_tiling == 0 ) 
  {
    fOpen_flag_tiling = fNumDisplay;
    ++fNumDisplay;
    d_num = fOpen_flag_tiling;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, fxs_t, fys_t, fxe_t, fye_t );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num);
    usleep(100000);   
  }
  d_num = fOpen_flag_tiling;
  g_clearWindow(d_num);

  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;
   
  float cos_th, sin_th;
  VectorXd dir(2); 
  MatrixXd rot(2,2); 
  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  dir(0) = u(2*fi[1]) - u(2*fi[0]) + u(2*fi[3]) - u(2*fi[2]); ; 
  dir(1) = u(2*fi[1]+1) - u(2*fi[0]+1) + u(2*fi[3]+1) - u(2*fi[2]+1); 
  dir.normalize();
  cos_th = dir(0);
  sin_th = -dir(1);
  rot(0,0) = cos_th;
  rot(0,1) = -sin_th;
  rot(1,0) = sin_th;
  rot(1,1) = cos_th;

  VectorXd u_ori_0(2*fN);
  VectorXd u_ori_1(2*fN);
  VectorXd u_ori_2(2*fN);
  VectorXd u_ori_3(2*fN);
  for( int i = 0; i < fN; ++i ){
    u_ori_0(2*i) = cos_th*u(2*i) - sin_th*u(2*i+1);
    u_ori_0(2*i+1) = sin_th*u(2*i) + cos_th*u(2*i+1);
  }
  for( int i=0; i < fN; ++i ){ 
    u_ori_1(2*i) = -u_ori_0(2*i);
    u_ori_1(2*i+1) = u_ori_0(2*i+1);
    u_ori_2(2*i) =  u_ori_0(2*i);
    u_ori_2(2*i+1) = -u_ori_0(2*i+1);
    u_ori_3(2*i) = - u_ori_0(2*i);
    u_ori_3(2*i+1) = - u_ori_0(2*i+1);
  }


  diff(0) = u_ori_1(2*fi[5]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_1(2*fi[5]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_1(2*i) -= diff(0);
    u_ori_1(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_2(2*fi[0]) - u_ori_0(2*fi[2]); 
  diff(1) = u_ori_2(2*fi[0]+1) - u_ori_0(2*fi[2]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_2(2*i) -= diff(0);
    u_ori_2(2*i+1) -= diff(1);
  }
  diff(0) = u_ori_3(2*fi[3]) - u_ori_0(2*fi[4]); 
  diff(1) = u_ori_3(2*fi[3]+1) - u_ori_0(2*fi[4]+1); 
  for( int i = 0; i < fN; ++i ){
    u_ori_3(2*i) -= diff(0);
    u_ori_3(2*i+1) -= diff(1);
  }

  diff_0(0) = u_ori_3(2*fi[0]) - u_ori_0(2*fi[5]);
  diff_0(1) = u_ori_3(2*fi[0]+1) - u_ori_0(2*fi[5]+1);
  diff_1(0) = u_ori_1(2*fi[2]) - u_ori_0(2*fi[5]);
  diff_1(1) = u_ori_1(2*fi[2]+1) - u_ori_0(2*fi[5]+1);

  double pointsD[fN][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN ; ++i ){
	x = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	y = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      g_setColor(d_num,10); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	x = u_ori_1(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	y = u_ori_1(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      g_setColor(d_num,15); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );

      for( int i = 0; i < fN ; ++i ){
	x = u_ori_2(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	y = u_ori_2(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      g_setColor(d_num,19); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );
      
      for( int i = 0; i < fN ; ++i ){
	x = u_ori_3(2*i)+diff_0(0)*k1+diff_1(0)*k2; 
	y = u_ori_3(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
	pointsD[i][0] = cos_th*x + sin_th*y;
	pointsD[i][1] = - sin_th*x + cos_th*y;
      }
      g_setColor(d_num,29); 
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );
    } 
  }

  g_flush();
}

void TDisplay::Tiling_IH1( VectorXd& u, int* fn, int nv )
{
  int d_num;

  if ( fOpen_flag_tiling == 0 ) 
  {
    fOpen_flag_tiling = fNumDisplay;
    ++fNumDisplay;
    d_num = fOpen_flag_tiling;

    g_open((char*)"");
    g_open_window( d_num, 800, 800, (char*)"Tiling" );
    g_window( d_num, fxs_t, fys_t, fxe_t, fye_t );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_set_ColorPixcel();
    // g_set_pixpat(d_num); 
    usleep(100000);   
  }
  d_num = fOpen_flag_tiling;
  g_clearWindow(d_num);

  int fi[7];
  fi[0] = 0;
  for( int k = 0; k < nv; ++k )
    fi[k+1] = fi[k]+fn[k]+1;

  VectorXd diff(2); 
  VectorXd diff_0(2); 
  VectorXd diff_1(2); 
  VectorXd u_ori_0(2*fN);
  for( int i = 0; i < fN; ++i ){
    u_ori_0(2*i) = u(2*i);
    u_ori_0(2*i+1) = u(2*i+1);
  }

  diff_0(0) = u_ori_0(2*fi[2]) - u_ori_0(2*fi[0]);
  diff_0(1) = u_ori_0(2*fi[2]+1) - u_ori_0(2*fi[0]+1);
  diff_1(0) = u_ori_0(2*fi[4]) - u_ori_0(2*fi[0]);
  diff_1(1) = u_ori_0(2*fi[4]+1) - u_ori_0(2*fi[0]+1);

  double pointsD[fN][2];
  double x, y;

  for( int k1 = -5; k1 <= 5; ++k1 ){
    for( int k2 = -5; k2 <= 5; ++k2 ){

      // polygon
      for( int i = 0; i < fN ; ++i ){
	pointsD[i][0] = u_ori_0(2*i)+diff_0(0)*k1+diff_1(0)*k2;
	pointsD[i][1] = u_ori_0(2*i+1)+diff_0(1)*k1+diff_1(1)*k2;
      }
      if( k1 % 2 == 0 ){
	if( k2 % 2 == 0 )
	  g_setColor(d_num,10); 
	else
	  g_setColor(d_num,19);
      }
      else {
	if( k2 % 2 == 0 )
	  g_setColor(d_num,29); 
	else
	  g_setColor(d_num,15);
      }
      g_polygonfill( d_num, pointsD, fN );
      g_setColor(d_num,1); 
      g_polygon( d_num, pointsD, fN, 4 );
    } 
  }

  g_flush();
}

// Size:4*5
void TDisplay::All_tiles( VectorXd* u )
{
  static int open_flg = 0;
  static int d_num = 0;

  double x_min, x_max, y_min, y_max, length;
  double xs, xe, ys, ye;
  x_min = 99999999.9; x_max = -99999999.9; y_min = 99999999.9; y_max = -99999999.9; 

  for( int i = 0; i < fN; ++i ){
    if( fW(2*i) < x_min )
      x_min = fW(2*i);
    if( fW(2*i) > x_max )
      x_max = fW(2*i);
    if( fW(2*i+1) < y_min )
      y_min = fW(2*i+1);
    if( fW(2*i+1) > y_max )
      y_max = fW(2*i+1);
  }

  if( x_max - x_min < y_max - y_min )
    length = y_max - y_min;
  else
    length = x_max - x_min;
  xs = (x_min+x_max)*0.5 - 0.55*length; 
  xe = (x_min+x_max)*0.5 + 0.55*length; 
  ys = (y_min+y_max)*0.5 - 0.55*length; 
  ye = (y_min+y_max)*0.5 + 0.55*length; 
  
  if (open_flg == 0) 
  {
    d_num = fNumDisplay;
    ++fNumDisplay;
    open_flg = 1;

    g_open((char*)"");
    g_open_window( d_num, 1750, 1400, (char*)"Tiling" );
    g_window( d_num, xs, ys, xs+5.0*length*1.1, ys+4.0*length*1.1 );
    // g_setFont(d_num,(char*)"-adobe-helvetica-medium-r-*-*-8-*-*-*-*-*-*-*");
    g_setFont(d_num,(char*)"-*-*-20-*-*-*-*-*-*-*"); // これはエラーでない
    g_set_ColorPixcel();
    /*---塗りつぶしの時使用するビットマップパターンをセット*/
    // g_set_pixpat(d_num); 
    usleep(100000);
  }
  g_clearWindow(d_num);


  
  char str[80];
  double ax = 0.0;
  double  ay = 0.0; 
  for( int s = 0; s < 20; ++s ){
    if( s % 5 == 0 && s != 0 ){
      ax = 0.0;
      ay += length*1.1;
    }

    // surface points
    g_setColor(d_num,2); // red
    for( int i = 0; i < fN; ++i )
      g_line( d_num, u[s](2*i)+ax, u[s](2*i+1)+ay, u[s](2*((i+1)%fN))+ax, u[s](2*((i+1)%fN)+1)+ay, 3 );


    // number
    g_setColor(d_num,1); // black
    sprintf( str, "%d", s+1 );
    g_string( d_num, xs+ax+100, ys+ay+length*1.1-100, str );
    ax += length*1.1;
  }
  g_flush();
}
