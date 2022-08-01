/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_BASE__
#include "operator_base.h"
#endif

TOperator_base::TOperator_base( int N ) 
{
  fN = N;

  fBe = MatrixXd(4*6,2*fN);
  fBd = MatrixXd(2*fN,2*fN);
  fBs = SparseMatrix<double>(2*fN,2*fN);
  fBs.reserve(8*fN);
  fTripletList.reserve(8*fN);

  fi = new int [ 7 ];
  fn = new int [ 6 ];
}


TOperator_base::~TOperator_base()
{

}


MatrixXd TOperator_base::Trans_AtoB( const MatrixXd& A ) 
{
  int m, n;
  m = A.rows();
  n = A.cols();
  
  VectorXd t(n); 
  MatrixXd B;
  B = A.fullPivLu().kernel();

  for( int i = 0; i < B.cols(); ++i ){
    t = B.col(i);
    for( int j = 0; j < i; ++j ){
      t -= t.dot(B.col(j))*B.col(j);
    }
    if( t.norm() < 0.0000001 )
      printf( "Warning: not independent!! \n" );
    B.col(i) = t.normalized();
  } 

  return B;
}


void TOperator_base::Orthogonalization_Bd()
{
  VectorXd t(2*fN); 

  for( int i = 0; i < fMd; ++i ){
    t = fBd.col(i);
    for( int j = 0; j < i; ++j ){
      t -= t.dot(fBd.col(j))*fBd.col(j);
    }
    if( t.norm() < 0.0000001 )
      printf( "Warning: not independent!! \n" );
    fBd.col(i) = t.normalized();
  } 
  // this->Check_Orthogonal();  // for check
}


void TOperator_base::Check_Orthogonal()
{
  // fBd-fBd
  for( int i = 0; i < fMd; ++i ){

    for( int j = 0; j < fMd; ++j ){
    if( i == j )
      assert( fabs( fBd.col(i).dot(fBd.col(j)) - 1.0 ) < 0.000001 );
    else
      assert( fabs( fBd.col(i).dot(fBd.col(j)) ) < 0.000001 );
    }
  }

  // fBs-fBs
  for( int i = 0; i < fMs; ++i ){
    for( int j = 0; j < fMs; ++j ){
    if( i == j )
      assert( fabs( fBs.col(i).dot(fBs.col(j)) - 1.0 ) < 0.000001 );
    else
      assert( fabs( fBs.col(i).dot(fBs.col(j)) ) < 0.000001 );
    }
  }

  // fBs-fBd
  for( int i = 0; i < fMs; ++i ){
    for( int j = 0; j < fMd; ++j ){
      assert( fabs( fBs.col(i).dot(fBd.col(j)) ) < 0.000001 );
    }
  }
}


void TOperator_base::Set_As_Parallel( int ivs, int ive, MatrixXd& As, int& m ) 
{
  assert( 0 <= ivs );  assert( ivs < fNv );
  assert( 0 <= ive );  assert( ive < fNv );

  As(m,ix(ivs+1)) = 1;
  As(m,ix(ivs)) = -1;
  As(m,ix((ive+1)%fNv)) = 1;
  As(m,ix(ive)) = -1;
  ++m;
  As(m,iy(ivs+1)) = 1;
  As(m,iy(ivs)) = -1;
  As(m,iy((ive+1)%fNv)) = 1;
  As(m,iy(ive)) = -1;
  ++m;
}


void TOperator_base::Set_As_Deg( int iv, int deg, MatrixXd& As, int& m ) 
{
  assert( 0 < iv );  assert( iv < fNv );
  double cos_th, sin_th;

  if( deg == 120 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == 90 ){
    sin_th = 1;
    cos_th = 0;
  }
  else if( deg == 60 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = 0.5;
  }
  else if( deg == -120 ){
    sin_th = -0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == -90 ){
    sin_th = -1;
    cos_th = 0;
  }
  else if( deg == -60 ){
    sin_th = -0.5*sqrt(3);  
    cos_th = 0.5;
  }
  else 
    assert( 1 == 2 );
    
  As(m,ix((iv+1)%fNv)) = 1;
  As(m,ix(iv)) = -1+cos_th;
  As(m,ix(iv-1)) = -cos_th;
  As(m,iy(iv)) = -sin_th;
  As(m,iy(iv-1)) = sin_th;
  ++m;

  As(m,iy((iv+1)%fNv)) = 1;
  As(m,iy(iv)) = -1+cos_th;
  As(m,iy(iv-1)) = -cos_th;
  As(m,ix(iv)) = sin_th;
  As(m,ix(iv-1)) = -sin_th;
  ++m;
}


void TOperator_base::Set_Be_Jd_Parallel( int ivs1, int ivs2 ) 
{
  VectorXd t(4*fNv); 
  int ive2;

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  ive2 = ivs2 + 1;
  if( ive2 == fNv )
    ive2 = 0;

  fBe.row(sx_f(ivs1)) = fBv.row(ix(ivs1));  
  fBe.row(sx_b(ive2)) = fBv.row(ix(ive2));  
  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs1)) = 1;
  t(sx_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs1)) = fBv.row(iy(ivs1));  
  fBe.row(sy_b(ive2)) = fBv.row(iy(ive2));  
  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs1)) = 1;
  t(sy_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_S_Parallel( int ivs1, int ivs2 ) 
{
  int ive1, ive2;
  VectorXd t(4*fNv); 

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  ive1 = ivs1 + 1;
  ive2 = ivs2 + 1;
  if( ive2 == fNv )
    ive2 = 0;

  fBe.row(sx_f(ivs1)) = fBv.row(ix(ivs1));  
  fBe.row(sx_b(ive1)) = fBv.row(ix(ive1));  
  fBe.row(sx_f(ivs2)) = fBv.row(ix(ivs2));  
  fBe.row(sx_b(ive2)) = fBv.row(ix(ive2));  

  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs1)) = 1;
  t(sx_b(ive1)) = -1;
  t(sx_f(ivs2)) = -1;
  t(sx_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs1)) = fBv.row(iy(ivs1));  
  fBe.row(sy_b(ive1)) = fBv.row(iy(ive1));  
  fBe.row(sy_f(ivs2)) = fBv.row(iy(ivs2));  
  fBe.row(sy_b(ive2)) = fBv.row(iy(ive2));  

  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs1)) = 1;
  t(sy_b(ive1)) = -1;
  t(sy_f(ivs2)) = -1;
  t(sy_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_J( int ivs ) 
{
  VectorXd t(4*fNv); 

  assert( 0 <= ivs );  assert( ivs < fNv );

  fBe.row(sx_f(ivs)) = fBv.row(ix(ivs));  
  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs)) = fBv.row(iy(ivs));  
  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_S( int ivs )
{
  VectorXd t(4*fNv); 
  int ive;

  assert( 0 <= ivs );  assert( ivs < fNv );

  if( ivs < fNv - 1 ) 
    ive = ivs+1;
  else 
    ive = 0;

  fBe.row(sx_f(ivs)) = fBv.row(ix(ivs));  
  fBe.row(sx_b(ive)) = fBv.row(ix(ive));  
  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs)) = 1;
  t(sx_b(ive)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs)) = fBv.row(iy(ivs));  
  fBe.row(sy_b(ive)) = fBv.row(iy(ive));  
  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs)) = 1;
  t(sy_b(ive)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_Js_Xsymmetry( int ivs1, int ivs2 ) 
{
  VectorXd t(4*fNv); 

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  fBe.row(sx_f(ivs1)) = fBv.row(ix(ivs1));  
  fBe.row(sx_f(ivs2)) = fBv.row(ix(ivs2));  
  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs1)) = 1;
  t(sx_f(ivs2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs1)) = fBv.row(iy(ivs1));  
  fBe.row(sy_f(ivs2)) = fBv.row(iy(ivs2));  
  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs1)) = 1;
  t(sy_f(ivs2)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_Js_Ysymmetry( int ivs1, int ivs2 ) 
{
  VectorXd t(4*fNv); 

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  fBe.row(sx_f(ivs1)) = fBv.row(ix(ivs1));  
  fBe.row(sx_f(ivs2)) = fBv.row(ix(ivs2));  
  t = VectorXd::Zero(4*fNv);
  t(sx_f(ivs1)) = 1;
  t(sx_f(ivs2)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  fBe.row(sy_f(ivs1)) = fBv.row(iy(ivs1));  
  fBe.row(sy_f(ivs2)) = fBv.row(iy(ivs2));  
  t = VectorXd::Zero(4*fNv);
  t(sy_f(ivs1)) = 1;
  t(sy_f(ivs2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Be_Jd_Deg( int iv, int deg ) 
{
  VectorXd t(4*fNv); 
  double cos_th, sin_th;

  assert( 0 < iv ); assert( iv < fNv ); 

  if( deg == 120 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == 90 ){
    sin_th = 1;
    cos_th = 0;
  }
  else if( deg == 60 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = 0.5;
  }
  else if( deg == -120 ){
    sin_th = -0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == -90 ){
    sin_th = -1;
    cos_th = 0;
  }
  else if( deg == -60 ){
    sin_th = -0.5*sqrt(3);  
    cos_th = 0.5;
  }

  else 
    assert( 1 == 2 );

  fBe.row(sx_b(iv)) = fBv.row(ix(iv));  
  fBe.row(sy_b(iv)) = fBv.row(iy(iv));  
  fBe.row(sx_f(iv)) = fBv.row(ix(iv));  
  fBe.row(sy_f(iv)) = fBv.row(iy(iv));  

  t = VectorXd::Zero(4*fNv);
  t(sx_b(iv)) = 1;
  t(sy_b(iv)) = 0;
  t(sx_f(iv)) = cos_th;
  t(sy_f(iv)) = sin_th;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;

  t = VectorXd::Zero(4*fNv);
  t(sx_b(iv)) = 0;
  t(sy_b(iv)) = 1;
  t(sx_f(iv)) = -sin_th;
  t(sy_f(iv)) = cos_th;
  t = t.normalized();
  for( int i = 0; i < fMd; ++i )
    fBe.col(i) -= fBe.col(i).dot(t)*t;
}


void TOperator_base::Set_Bd_Jd_Parallel( int ivs1, int ivs2 ) 
{
  int n, is1, ie2;
  int ive2;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );
  
  ive2 = ivs2 + 1;
  if( ive2 == fNv ) 
    ive2 = 0;
  
  assert( fn[ivs1] == fn[ivs2] );
  n = fn[ivs1];
  is1 = fi[ivs1];
  ie2 = fi[ivs2+1]; 

  for( int i = 1; i <= n; ++i ){ 
    fBd.row(ix(is1+i)) = fBe.block(sx_f(ivs1),0,1,fMd); 
    fTripletList.push_back(T(ix(is1+i),fMs,inv_sq2));
    fBd.row(ix(ie2-i)) = fBe.block(sx_b(ive2),0,1,fMd); 
    fTripletList.push_back(T(ix(ie2-i),fMs,inv_sq2)); 
    ++fMs;

    fBd.row(iy(is1+i)) = fBe.block(sy_f(ivs1),0,1,fMd);
    fTripletList.push_back(T(iy(is1+i),fMs,inv_sq2)); 
    fBd.row(iy(ie2-i)) = fBe.block(sy_b(ive2),0,1,fMd);
    fTripletList.push_back(T(iy(ie2-i),fMs,inv_sq2));
    ++fMs;
  }
}


void TOperator_base::Set_Bd_S_Parallel( int ivs1, int ivs2 )
{
  int n;
  int ive1, ive2;
  int is1, ie1, is2, ie2;
  const double inv_sq4 = 0.5;

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  ive1 = ivs1 + 1;
  ive2 = ivs2 + 1;
  if( ive2 == fNv )
    ive2 = 0;
  
  n = fn[ivs1];
  is1 = fi[ivs1];
  ie1 = fi[ivs1+1];
  is2 = fi[ivs2];
  ie2 = fi[ivs2+1];

  for( int i = 1; i <= n/2; ++i ){
    fBd.row(ix(is1+i)) = fBe.block(sx_f(ivs1),0,1,fMd); 
    fTripletList.push_back(T(ix(is1+i),fMs,inv_sq4));
    fBd.row(ix(ie1-i)) = fBe.block(sx_b(ive1),0,1,fMd); 
    fTripletList.push_back(T(ix(ie1-i),fMs,-inv_sq4));
    fBd.row(ix(is2+i)) = fBe.block(sx_f(ivs2),0,1,fMd); 
    fTripletList.push_back(T(ix(is2+i),fMs,-inv_sq4));
    fBd.row(ix(ie2-i)) = fBe.block(sx_b(ive2),0,1,fMd); 
    fTripletList.push_back(T(ix(ie2-i),fMs,inv_sq4));
    ++fMs;

    fBd.row(iy(is1+i)) = fBe.block(sy_f(ivs1),0,1,fMd); 
    fTripletList.push_back(T(iy(is1+i),fMs,inv_sq4));
    fBd.row(iy(ie1-i)) = fBe.block(sy_b(ive1),0,1,fMd); 
    fTripletList.push_back(T(iy(ie1-i),fMs,-inv_sq4));
    fBd.row(iy(is2+i)) = fBe.block(sy_f(ivs2),0,1,fMd); 
    fTripletList.push_back(T(iy(is2+i),fMs,-inv_sq4));
    fBd.row(iy(ie2-i)) = fBe.block(sy_b(ive2),0,1,fMd); 
    fTripletList.push_back(T(iy(ie2-i),fMs,inv_sq4));
    ++fMs;
  }

  if( n % 2 == 1 ){
    fBd.row(ix(is1+(n+1)/2)) = 0.5*(fBv.block(ix(ivs1),0,1,fMd)+fBv.block(ix(ive1),0,1,fMd));
    fBd.row(iy(is1+(n+1)/2)) = 0.5*(fBv.block(iy(ivs1),0,1,fMd)+fBv.block(iy(ive1),0,1,fMd));
    fBd.row(ix(is2+(n+1)/2)) = 0.5*(fBv.block(ix(ivs2),0,1,fMd)+fBv.block(ix(ive2),0,1,fMd));
    fBd.row(iy(is2+(n+1)/2)) = 0.5*(fBv.block(iy(ivs2),0,1,fMd)+fBv.block(iy(ive2),0,1,fMd));
  }
}


void TOperator_base::Set_Bd_S( int ivs ) 
{
  int n, ive, is, ie;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs ); assert( ivs < fNv ); 
  
  if( ivs < fNv - 1 ) 
    ive = ivs+1;
  else 
    ive = 0;

  n = fn[ivs];
  is = fi[ivs];
  ie = fi[ivs+1];

  for( int i = 1; i <= n/2; ++i ){
    fBd.row(ix(is+i)) = fBe.block(sx_f(ivs),0,1,fMd); 
    fTripletList.push_back(T(ix(is+i),fMs,inv_sq2));
    fBd.row(ix(ie-i)) = fBe.block(sx_b(ive),0,1,fMd); 
    fTripletList.push_back(T(ix(ie-i),fMs,-inv_sq2));
    ++fMs;

    fBd.row(iy(is+i)) = fBe.block(sy_f(ivs),0,1,fMd); 
    fTripletList.push_back(T(iy(is+i),fMs,inv_sq2));
    fBd.row(iy(ie-i)) = fBe.block(sy_b(ive),0,1,fMd); 
    fTripletList.push_back(T(iy(ie-i),fMs,-inv_sq2));
    ++fMs;
  }
  if( n % 2 == 1 ){
    fBd.row(ix(is+(n+1)/2)) = 0.5*(fBv.block(ix(ivs),0,1,fMd)+fBv.block(ix(ive),0,1,fMd));
    fBd.row(iy(is+(n+1)/2)) = 0.5*(fBv.block(iy(ivs),0,1,fMd)+fBv.block(iy(ive),0,1,fMd));
  }
}


void TOperator_base::Set_Bd_Js_Xsymmetry( int ivs1, int ivs2 ) 
{
  int n, is1, is2;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  n = fn[ivs1];
  is1 = fi[ivs1];
  is2 = fi[ivs2];

  for( int i = 1; i <= n; ++i ){ 
    fBd.row(ix(is1+i)) = fBe.block(sx_f(ivs1),0,1,fMd); 
    fTripletList.push_back(T(ix(is1+i),fMs,inv_sq2));
    fBd.row(ix(is2+i)) = fBe.block(sx_f(ivs2),0,1,fMd); 
    fTripletList.push_back(T(ix(is2+i),fMs,inv_sq2)); 
    ++fMs;

    fBd.row(iy(is1+i)) = fBe.block(sy_f(ivs1),0,1,fMd);
    fTripletList.push_back(T(iy(is1+i),fMs,inv_sq2)); 
    fBd.row(iy(is2+i)) = fBe.block(sy_f(ivs2),0,1,fMd);
    fTripletList.push_back(T(iy(is2+i),fMs,-inv_sq2));
    ++fMs;
  }
}


void TOperator_base::Set_Bd_Js_Ysymmetry( int ivs1, int ivs2 ) 
{
  int n, is1, is2;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs1 );  assert( ivs1 < fNv );
  assert( 0 <= ivs2 );  assert( ivs2 < fNv );
  assert( ivs1 < ivs2 );

  n = fn[ivs1];
  is1 = fi[ivs1];
  is2 = fi[ivs2];

  for( int i = 1; i <= n; ++i ){ 
    fBd.row(ix(is1+i)) = fBe.block(sx_f(ivs1),0,1,fMd); 
    fTripletList.push_back(T(ix(is1+i),fMs,inv_sq2));
    fBd.row(ix(is2+i)) = fBe.block(sx_f(ivs2),0,1,fMd); 
    fTripletList.push_back(T(ix(is2+i),fMs,-inv_sq2)); 
    ++fMs;

    fBd.row(iy(is1+i)) = fBe.block(sy_f(ivs1),0,1,fMd);
    fTripletList.push_back(T(iy(is1+i),fMs,inv_sq2)); 
    fBd.row(iy(is2+i)) = fBe.block(sy_f(ivs2),0,1,fMd);
    fTripletList.push_back(T(iy(is2+i),fMs,inv_sq2));
    ++fMs;
  }
}


void TOperator_base::Set_Bd_Jd_Deg( int iv, int deg ) 
{
  int n, is;
  const double sq3 = sqrt(3.0);
  const double inv_sq2 = 1.0/sqrt(2.0);
  double cos_th, sin_th;

  assert( 0 < iv );  assert( iv < fNv );

  if( deg == 120 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == 90 ){
    sin_th = 1;
    cos_th = 0;
  }
  else if( deg == 60 ){
    sin_th = 0.5*sqrt(3);  
    cos_th = 0.5;
  }
  else if( deg == -120 ){
    sin_th = - 0.5*sqrt(3);  
    cos_th = -0.5;
  }
  else if( deg == -90 ){
    sin_th = -1;
    cos_th = 0;
  }
  else if( deg == -60 ){
    sin_th = -0.5*sqrt(3);  
    cos_th = 0.5;
  }
  else 
    assert( 1 == 2 );

  n = fn[iv];
  is = fi[iv];

  for( int i = 1; i <= n; ++i ){ 
    fBd.row(ix(is-i)) = fBe.block(sx_b(iv),0,1,fMd); 
    fBd.row(iy(is-i)) = fBe.block(sy_b(iv),0,1,fMd); 
    fBd.row(ix(is+i)) = fBe.block(sx_f(iv),0,1,fMd); 
    fBd.row(iy(is+i)) = fBe.block(sy_f(iv),0,1,fMd); 

    fTripletList.push_back(T(ix(is-i),fMs,inv_sq2));
    fTripletList.push_back(T(ix(is+i),fMs,inv_sq2*cos_th));
    fTripletList.push_back(T(iy(is+i),fMs,inv_sq2*sin_th));
    ++fMs;

    fTripletList.push_back(T(iy(is-i),fMs,inv_sq2));
    fTripletList.push_back(T(ix(is+i),fMs,-inv_sq2*sin_th));
    fTripletList.push_back(T(iy(is+i),fMs,inv_sq2*cos_th));
    ++fMs;
  }
}

///////////////////////////////////////////////////////////////////////////
//// IH4 ////
void TOperator_base::Make_Bv_Be_IH4() 
{
  int m;

  fIH = 4;
  fNv = 6; 
  int m_max = 2; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Parallel(0,3,As,m); 
  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Parallel(0,3); 
  Set_Be_S(1); 
  Set_Be_S(2); 
  Set_Be_S(4); 
  Set_Be_S(5); 
}

void TOperator_base::Make_B_IH4( int k1, int k2, int k3, int k4 )
{
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k3; 
  fn[3] = k1; 
  fn[4] = k4; 
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Parallel(0,3); // IH
  Set_Bd_S(1); 
  Set_Bd_S(2); 
  Set_Bd_S(4); 
  Set_Bd_S(5); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH5 ////
void TOperator_base::Make_Bv_Be_IH5() 
{
  int m;

  fIH = 5;
  fNv = 6; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Parallel(0,3,As,m);  

  As(m,iy(1)) = 1; 
  As(m,iy(3)) = -1;
  ++m;
  As(m,ix(2)) = 1;
  As(m,ix(1)) = -0.5;
  As(m,ix(3)) = -0.5;
  ++m;

  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Parallel(0,3); 
  Set_Be_Js_Xsymmetry(1,2); 
  Set_Be_S(4); 
  Set_Be_S(5); 
}

void TOperator_base::Make_B_IH5( int k1, int k2, int k3 ) 
{
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k2; 
  fn[3] = k1; 
  fn[4] = k3; 
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Parallel(0,3); 
  Set_Bd_Js_Xsymmetry(1,2); 
  Set_Bd_S(4); 
  Set_Bd_S(5); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH6 ////
void TOperator_base::Make_Bv_Be_IH6() 
{
  int m;

  fIH = 6;
  fNv = 6; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  As(m,ix(1)) = 1;
  As(m,ix(0)) = -1;
  As(m,ix(3)) = -1;
  As(m,ix(2)) = 1;
  ++m;
  As(m,iy(1)) = 1;
  As(m,iy(0)) = -1;
  As(m,iy(3)) = 1;
  As(m,iy(2)) = -1;
  ++m;
  As(m,ix(2)) = 1;
  As(m,ix(1)) = -1;
  As(m,ix(5)) = 1;
  As(m,ix(4)) = -1;
  ++m;
  As(m,iy(2)) = 1;
  As(m,iy(1)) = -1;
  As(m,iy(5)) = -1;
  As(m,iy(4)) = 1;
  ++m;

  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Js_Xsymmetry(0,2); 
  Set_Be_Js_Ysymmetry(1,4); 
  Set_Be_S(3); 
  Set_Be_S(5); 
}

void TOperator_base::Make_B_IH6( int k1, int k2, int k3 ) 
{
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k1; 
  fn[3] = k3; 
  fn[4] = k2; 
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Js_Xsymmetry(0,2); 
  Set_Bd_Js_Ysymmetry(1,4); 
  Set_Bd_S(3); 
  Set_Bd_S(5); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH1 ////
void TOperator_base::Make_Bv_Be_IH1() 
{
  int m;

  fIH = 1;
  fNv = 6; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Parallel(0,3,As,m);  
  Set_As_Parallel(1,4,As,m);  
  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Parallel(0,3); 
  Set_Be_Jd_Parallel(1,4); 
  Set_Be_Jd_Parallel(2,5); 
}

void TOperator_base::Make_B_IH1( int k1, int k2 ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = (fN-fNv-2*k1-2*k2)/2;  
  fn[3] = k1; 
  fn[4] = k2;
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 
  assert( fn[2] == fn[5] );

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Parallel(0,3); 
  Set_Bd_Jd_Parallel(1,4); 
  Set_Bd_Jd_Parallel(2,5); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH2 ////
void TOperator_base::Make_Bv_Be_IH2() 
{
  int m;

  fIH = 2;
  fNv = 6; 
  int m_max = 5; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Parallel(0,3,As,m);  

  As(m,iy(1)) = 1; 
  As(m,iy(3)) = -1;
  ++m;
  As(m,ix(2)) = 1;
  As(m,ix(1)) = -0.5;
  As(m,ix(3)) = -0.5;
  ++m;
  As(m,ix(5)) = 1;
  As(m,ix(4)) = -0.5;
  As(m,ix(0)) = -0.5;
  ++m;

  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Parallel(0,3); 
  Set_Be_Js_Xsymmetry(1,2); 
  Set_Be_Js_Xsymmetry(4,5); 
}

void TOperator_base::Make_B_IH2( int k1, int k2 ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k2; 
  fn[3] = k1; 
  fn[4] = (fN-fNv-2*k1-2*k2)/2; 
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 
  assert( fn[4] == fn[5] );

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Parallel(0,3); 
  Set_Bd_Js_Xsymmetry(1,2); 
  Set_Bd_Js_Xsymmetry(4,5); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH3 ////
void TOperator_base::Make_Bv_Be_IH3() 
{
  int m;

  fIH = 3;
  fNv = 6; 
  int m_max = 5; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Parallel(0,3,As,m);  

  As(m,iy(1)) = 1; 
  As(m,iy(3)) = -1;
  ++m;
  As(m,ix(2)) = 1;
  As(m,ix(1)) = -1;
  As(m,ix(5)) = -1;
  As(m,ix(0)) = 1;
  ++m;
  As(m,iy(2)) = 1;
  As(m,iy(1)) = -1;
  As(m,iy(5)) = 1;
  As(m,iy(0)) = -1;
  ++m;

  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Parallel(0,3); 
  Set_Be_Js_Ysymmetry(1,5); 
  Set_Be_Js_Ysymmetry(2,4); 
}

void TOperator_base::Make_B_IH3( int k1, int k2 ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = (fN-fNv-2*k1-2*k2)/2; 
  fn[3] = k1; 
  fn[4] = (fN-fNv-2*k1-2*k2)/2; 
  fn[5] = k2;
  assert( fn[2] == fn[4] );

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Parallel(0,3); 
  Set_Bd_Js_Ysymmetry(1,5); 
  Set_Bd_Js_Ysymmetry(2,4); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH7 ////
void TOperator_base::Make_Bv_Be_IH7( int sign ) 
{
  int m;

  fIH = 7;
  fNv = 6; 
  int m_max = 6; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Deg(1,sign*120,As,m);  
  Set_As_Deg(3,sign*120,As,m);  
  Set_As_Deg(5,sign*120,As,m);  
  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Deg(1,sign*120); 
  Set_Be_Jd_Deg(3,sign*120); 
  Set_Be_Jd_Deg(5,sign*120); 
}

void TOperator_base::Make_B_IH7( int k1, int k2, int sign ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k1; 
  fn[2] = k2; 
  fn[3] = k2; 
  fn[4] = (fN-fNv-fn[0]-fn[1]-fn[2]-fn[3])/2; 
  fn[5] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]-fn[4]; 
  assert( fn[4] == fn[5] );

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Deg(1,sign*120); 
  Set_Bd_Jd_Deg(3,sign*120); 
  Set_Bd_Jd_Deg(5,sign*120); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH21 ////
void TOperator_base::Make_Bv_Be_IH21( int sign ) 
{
  int m;

  fIH = 21;
  fNv = 5; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Deg(1,sign*120,As,m);  
  Set_As_Deg(3,sign*60,As,m);  
  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Deg(1,sign*120); 
  Set_Be_Jd_Deg(3,sign*60); 
  Set_Be_S(4); 
}

void TOperator_base::Make_B_IH21( int k1, int k2, int sign ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k1; 
  fn[2] = k2; 
  fn[3] = k2; 
  fn[4] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]; 

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Deg(1,sign*120); 
  Set_Bd_Jd_Deg(3,sign*60); 
  Set_Bd_S(4); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

//// IH28 ////
void TOperator_base::Make_Bv_Be_IH28( int sign ) 
{
  int m;

  fIH = 28;
  fNv = 5; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNv);
  m = 0;
  Set_As_Deg(1,sign*90,As,m);  
  Set_As_Deg(3,sign*90,As,m);  
  assert( m_max == m );

  fBv = this->Trans_AtoB( As );
  fMd = fBv.cols();
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBd.resize(2*fN,fMd); 
  fBs.resize(2*fN,2*fN);

  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_Jd_Deg(1,sign*90); 
  Set_Be_Jd_Deg(3,sign*90); 
  Set_Be_S(4); 
}

void TOperator_base::Make_B_IH28( int k1, int k2, int sign ) 
{
  assert( fN %2 == 0 );
  fn[0] = k1; 
  fn[1] = k1; 
  fn[2] = k2; 
  fn[3] = k2; 
  fn[4] = fN-fNv-fn[0]-fn[1]-fn[2]-fn[3]; 

  fi[0] = 0;
  for( int s = 1; s <= fNv; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  assert( fi[fNv] == fN );

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_Jd_Deg(1,sign*90); 
  Set_Bd_Jd_Deg(3,sign*90); 
  Set_Bd_S(4); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

