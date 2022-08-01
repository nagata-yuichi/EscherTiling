/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_E__
#include "operator_E.h"
#endif

TOperator_E::TOperator_E( int N ) : TOperator_base( N )
{
  fN = N;
  fW = VectorXd::Zero(2*fN);
  fW_d = VectorXd::Zero(2*fN);
  fWj = new VectorXd [fN];
  fWj_d = new VectorXd [fN];
  for( int i = 0; i < fN; ++i ){
    fWj[i] = VectorXd::Zero(2*fN);
    fWj_d[i] = VectorXd::Zero(2*fN);
  }

  fEval_Comp_best = new double** [ fN ];
  for( int i = 0; i < fN; ++i )
    fEval_Comp_best[ i ] = new double* [ fN ];
  for( int i = 0; i < fN; ++i )
    for( int j = 0; j < fN; ++j )
      fEval_Comp_best[ i ][ j ] = new double [ fN ];

  fEval_Comp_best_Sf = new double* [ fN ];
  for( int i = 0; i < fN; ++i )
    fEval_Comp_best_Sf[ i ] = new double [ fN ];
  fEval_Comp_best_Sb = new double* [ fN ];
  for( int i = 0; i < fN; ++i )
    fEval_Comp_best_Sb[ i ] = new double [ fN ];
  fEval_part = new double [ fN ]; 

  fBd_part = MatrixXd(2*fN,2*fN);
  fBs_part = SparseMatrix<double>(2*fN,2*fN);
  fBs_part.reserve(8*fN);
  fBe_part = MatrixXd(4*6,2*fN);
  fBv_part = MatrixXd(2*fN,2*fN);
}


TOperator_E::~TOperator_E()
{

}


void TOperator_E::SetParameter() 
{
  fFlagPrune = 1;             
  fFlagCheckIntersection = 1; // reset in search_E_conf.cpp
}


void TOperator_E::SetInit( VectorXd w ) 
{
  fW = w;
  for( int i = 0; i < fN; ++i ){
    fW_d(ix(i)) = fW(iy(i)); 
    fW_d(iy(i)) = -fW(ix(i)); 
  }
  
  for( int st1 = 0; st1 < fN; ++st1 ){
    for( int k = 0; k < fN; ++k ){
      fWj[st1](ix(k))=fW(ix(k+st1));
      fWj[st1](iy(k))=fW(iy(k+st1));
    }
    
    for( int k = 0; k < fN; ++k ){
      fWj_d[st1](ix(k))= fW(iy(k+st1)); 
      fWj_d[st1](iy(k))= -fW(ix(k+st1)); 
    }
  }
}


void TOperator_E::Cal_u( int k1, int k2, int k3, int k4 ) 
{
  int st1, st2;
  int numOfCheckS;
  int checkS[fN];

  if( fIH == 4 && fFlagPrune ){ // prune the search with the decomposition relaxation (IH4)
    int k5 = fN - fNv - (2*k1 + k2 + k3 + k4);
    numOfCheckS = 0;
    for( int s = 0; s < fN; ++s ){
      st1 = s;
      st2 = (st1 + 2*k1+k2+k3+4)%fN; 
      if( fEval_part[st1] + fEval_Comp_best[ st2 ][ k4 ][ k5 ] < fEval_best ){
	checkS[numOfCheckS++] = st1;
      }
    }
    if( numOfCheckS == 0 )
      return;

    this->Make_B_IH4( k1, k2, k3, k4 ); 
  }
  else{ 
    numOfCheckS = fN;
    for( int s = 0; s < fN; ++s )
      checkS[s] = s;
  }

  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fN);
  double eval;
  static double eval_ww_full = fW.dot(fW);

  // compute the tile shapes for only the selected st1 values
  // st1: the start point for the numbering of the goal polygon (index j in the paper)
  for( int s = 0; s < numOfCheckS; ++s ){
    st1 = checkS[s]; 

    gzai.head(fMd) = fBd.transpose() * fWj[st1]; 
    gzai.tail(fMs) = fBs.leftCols(fMs).transpose() * fWj[st1];
    eval = -gzai.dot(gzai) + eval_ww_full;

    if( eval < fEval_best ){
      u = fBd*gzai.head(fMd) + fBs.leftCols(fMs)*gzai.tail(fMs);
      if( fFlagCheckIntersection == 0 || Check_Intersection( u ) == false ){
	this->Update_U_top( eval, u, st1 );
      }
    }
  }
}

void TOperator_E::Cal_u_Procrustes( int k1, int k2, int k3, int k4 ) 
{
  int st1, st2;
  int numOfCheckS;
  int checkS[fN];

  if( fIH == 5 && fFlagPrune ){ // prune the search with the decomposition relaxation (IH5)
    k4 =  fN - fNv - (2*k1 + 2*k2 + k3);
    numOfCheckS = 0;
    for( int s = 0; s < fN; ++s ){
      st1 = s;
      st2 = (st1 + 2*k1+2*k2+4)%fN; 
      if( fEval_part[st1] + fEval_Comp_best[ st2 ][ k3 ][ k4 ] < fEval_best ){
	checkS[numOfCheckS++] = st1;
      }
    }
    if( numOfCheckS == 0 )
      return;

    this->Make_B_IH5( k1, k2, k3 );  
  }
  else if( fIH == 6 && fFlagPrune ){ // prune the search with the decomposition relaxation (IH6)
    k4 = fN - fNv - (2*k1 + 2*k2 + k3);
    numOfCheckS = 0;
    for( int s = 0; s < fN; ++s ){
      st1 = s;
      st2 = (st1 + 2*k1+k2+3)%fN; 
      if( fEval_part[st1] + fEval_Comp_best_Sf[ st2 ][ k3 ] + fEval_Comp_best_Sb[ st1 ][ k4 ] < fEval_best ){
	checkS[numOfCheckS++] = st1;
      }
    }

    if( numOfCheckS == 0 )
      return;

    this->Make_B_IH6( k1, k2, k3 );  
  }
  else{
    numOfCheckS = fN;
    for( int s = 0; s < fN; ++s )
      checkS[s] = s;
  }

  VectorXd gzai(fMd+fMs);
  VectorXd ur(2*fN);
  VectorXd u(2*fN);
  VectorXd p(fMd+fMs);
  VectorXd p_d(fMd+fMs);
  double eval;
  static double eval_ww_full = fW.dot(fW);
  double tr, rt, cos_th, sin_th;
  double lambda;
  int ne;
  VectorXd x_unit(fMd+fMs); 
  VectorXd y_unit(fMd+fMs); 
  VectorXd p2(2);
  VectorXd p_d2(2); 
  MatrixXd A2(2,2);
  double a_tr;

  // compute the tile shapes for only the selected st1 values
  // st1: the start point for the numbering of the goal polygon (index j in the paper)
  for( int s = 0; s < numOfCheckS; ++s ){
    st1 = checkS[s];

    p.head(fMd) = fBd.transpose()*fWj[st1];
    p.tail(fMs) = fBs.leftCols(fMs).transpose()*fWj[st1];
    p_d.head(fMd) = fBd.transpose()*fWj_d[st1];
    p_d.tail(fMs) = fBs.leftCols(fMs).transpose()*fWj_d[st1];
    x_unit = p.normalized();
    y_unit = p_d - p_d.dot(x_unit)*x_unit;
    y_unit = y_unit.normalized();
    p2(0) = p.norm(); 
    p2(1) = 0;        
    p_d2(0) = p_d.dot(x_unit);
    p_d2(1) = p_d.dot(y_unit);
    A2 = p2*p2.transpose() + p_d2*p_d2.transpose();
    a_tr = A2(0,0)+A2(1,1);
    lambda = 0.5*(a_tr+sqrt(a_tr*a_tr - 4.0*(A2(0,0)*A2(1,1)-A2(0,1)*A2(1,0))));

    eval = eval_ww_full - lambda;

    if( eval < fEval_best ){
      gzai = A2(0,1)*x_unit + (lambda-A2(0,0))*y_unit;
      gzai = gzai.normalized()*sqrt(lambda);

      u = fBd*gzai.head(fMd) + fBs.leftCols(fMs)*gzai.tail(fMs);

      if( fFlagCheckIntersection == 0 || Check_Intersection( u ) == false ){
	tr = u.transpose()*fWj[st1]; 
	rt = u.transpose()*fWj_d[st1]; 
	cos_th = tr/sqrt(tr*tr + rt*rt);   
	sin_th = rt/sqrt(tr*tr + rt*rt);
      
	for( int k = 0; k < fN; ++k ){
	  ur(ix(k)) = cos_th*u(ix(k)) - sin_th*u(iy(k)); 
	  ur(iy(k)) = sin_th*u(ix(k)) + cos_th*u(iy(k)); 
	}
	
	this->Update_U_top( eval, ur, st1 );
      }
    }
  }
}


void TOperator_E::Update_U_top( double eval, VectorXd& u, int st1 ) 
{

}


void TOperator_E::Trans_U_center( VectorXd& u )
{
  VectorXd centerU(2);
  VectorXd centerW(2);

  centerU = VectorXd::Zero(2); 
  centerW = VectorXd::Zero(2); 
  
  for( int i = 0; i < fN; ++i ){  
    centerW += fW.segment(2*i,2);
    centerU += u.segment(2*i,2);
  }
  centerW /= (double)fN;
  centerU /= (double)fN;

  for( int i = 0; i < fN; ++i )  
    u.segment(2*i,2) += (-centerU + centerW);
}


bool TOperator_E::Check_Intersection( VectorXd& u )
{
  double x1, x2, x3, x4, y1, y2, y3, y4, tc, td, value1, value2;
  for( int i = 0; i < fN; ++ i ){
    for( int j = 0; j < i; ++ j ){
      x1 = u(ix(i));
      y1 = u(iy(i));
      x2 = u(ix((i+1)%fN));
      y2 = u(iy((i+1)%fN));
      x3 = u(ix(j));
      y3 = u(iy(j));
      x4 = u(ix((j+1)%fN));
      y4 = u(iy((j+1)%fN));

      tc = (x1-x2)*(y3-y1)+(y1-y2)*(x1-x3);
      td = (x1-x2)*(y4-y1)+(y1-y2)*(x1-x4);
      value1 = tc*td;

      tc = (x3-x4)*(y1-y3)+(y3-y4)*(x3-x1);
      td = (x3-x4)*(y2-y3)+(y3-y4)*(x3-x2);
      value2 = tc*td;

      if( value1 < 0 && value2 < 0 )
	return true;
    }
  }
  return false;
}


/////////////////// Comp //////////////////////
void TOperator_E::Search_Comp_SS() 
{
  this->Make_Bv_Be_Comp_SS();

  for( int i = 0; i < fN; ++i ) 
    for( int j = 0; j < fN; ++j ){
      for( int k = 0; k < fN; ++k ) 
	fEval_Comp_best[ i ][ j ][ k ] = 99999999.9;
    }

  for( int k1 = 0; k1 < fN-fNv; ++k1 ){ 
    for( int k2 = 0; k2 < fN-fNv-k1; ++k2 ){ 
      this->Make_B_Comp_SS( k1, k2 );
      this->Cal_u_Comp_SS( k1, k2 );
    }
  }
}

void TOperator_E::Make_Bv_Be_Comp_SS() 
{
  int m;

  fNv = 3; 

  int m_max = 0;
  fBv = MatrixXd::Identity(2*fNv,2*fNv); 
  fMd = fBv.cols();
  fMv = fMd; // comp
  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_S(0); 
  Set_Be_S(1); 
}

void TOperator_E::Make_B_Comp_SS( int k1, int k2 ) 
{
  const double inv_sq2 = 1.0/sqrt(2.0);
  int col_s; 
  
  fn[0] = k1; 
  fn[1] = k2; 

  fi[0] = 0;
  for( int s = 1; s <= fNv-1; ++s )  
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNv-1]-1; 

  fMd = fMv; 
  if( k1 != 0 && k2 != 0 ){ // fMd is not changed
    fBd.resize(2*fNp,fMd); 
    col_s = 0;
  }
  else if( k1 == 0 && k2 != 0 ){ // fMd is changed
    fMd -= 2;
    fBd.resize(2*fNp,fMd); 
    col_s = 2;
  }
  else if( k1 != 0 && k2 == 0 ){ // fMd is changed
    fMd -= 2;
    fBd.resize(2*fNp,fMd); 
    col_s = 0;
  }
  else{
    fMd -= 4;
    fBd.resize(2*fNp,fMd); 
    col_s = 2;
  }

  fBd.row(ix(fi[1])-2) = fBv.block(ix(1),col_s,1,fMd);  
  fBd.row(iy(fi[1])-2) = fBv.block(iy(1),col_s,1,fMd); 

  fBs.resize(2*fNp,2*fNp);
  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_S_Comp(0, col_s); 
  Set_Bd_S_Comp(1, col_s); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

void TOperator_E::Cal_u_Comp_SS( int k1, int k2 ) 
{
  VectorXd w(2*fNp); 
  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fNp); 
  static double eval;
  static double eval_ww;
  double d;
  int st;

  for( int i = 0; i < fN; ++i ){
    st = i;
    for( int k = 0; k < fNp; ++k ){
      w(ix(k))=fW(ix((k+st+1)%fN));
      w(iy(k))=fW(iy((k+st+1)%fN));
    }
    eval_ww = w.dot(w); 

    gzai.head(fMd) = fBd.transpose() * w; 
    gzai.tail(fMs) = fBs.leftCols(fMs).transpose() * w;

    assert( fNp == k1+k2+1 );
    eval = -gzai.dot(gzai) + eval_ww;
 
    fEval_Comp_best[ st ][ k1 ][ k2 ] = eval;
  }
}

void TOperator_E::Set_Bd_S_Comp( int ivs, int col_s ) 
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
    fBd.row(ix(is+i)-2) = fBe.block(sx_f(ivs),col_s,1,fMd); 
    fTripletList.push_back(T(ix(is+i)-2,fMs,inv_sq2));
    fBd.row(ix(ie-i)-2) = fBe.block(sx_b(ive),col_s,1,fMd); 
    fTripletList.push_back(T(ix(ie-i)-2,fMs,-inv_sq2));
    ++fMs;

    fBd.row(iy(is+i)-2) = fBe.block(sy_f(ivs),col_s,1,fMd); 
    fTripletList.push_back(T(iy(is+i)-2,fMs,inv_sq2));
    fBd.row(iy(ie-i)-2) = fBe.block(sy_b(ive),col_s,1,fMd); 
    fTripletList.push_back(T(iy(ie-i)-2,fMs,-inv_sq2));
    ++fMs;
  }
  if( n % 2 == 1 ){
    fBd.row(ix(is+(n+1)/2)-2) = 0.5*(fBv.block(ix(ivs),col_s,1,fMd)+fBv.block(ix(ive),col_s,1,fMd));
    fBd.row(iy(is+(n+1)/2)-2) = 0.5*(fBv.block(iy(ivs),col_s,1,fMd)+fBv.block(iy(ive),col_s,1,fMd));
  }
}

void TOperator_E::Search_Comp_SJS() 
{
  this->Make_Bv_Be_Comp_S();

  for( int i = 0; i < fN; ++i ) {
    for( int j = 0; j < fN; ++j ) {
      fEval_Comp_best_Sf[ i ][ j ] = 99999999.9;
      fEval_Comp_best_Sb[ i ][ j ] = 99999999.9;
    }
  }

  for( int k1 = 0; k1 < fN-fNv; ++k1 ){ 
    this->Make_B_Comp_S( k1 );
    this->Cal_u_Comp_Sfb( k1 );
  }
}

void TOperator_E::Make_Bv_Be_Comp_S() 
{
  int m;

  fNv = 2; 

  int m_max = 0;
  fBv = MatrixXd::Identity(2*fNv,2*fNv);
  fMd = fBv.cols(); 
  fMv = fMd; // comp

  assert( fBv.rows() == 2*fNv );
  assert( fBv.cols() == 2*fNv-m_max );

  fBe.resize(4*fNv,fMd);
  fBe = MatrixXd::Zero(4*fNv,fMd);

  Set_Be_S(0); 
}

void TOperator_E::Make_B_Comp_S( int k1 ) 
{
  const double inv_sq2 = 1.0/sqrt(2.0);
  int col_s; 
  
  fn[0] = k1; 

  fi[0] = 0;
  for( int s = 1; s <= fNv-1; ++s )  
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNv-1];

  fMd = fMv; 
  if( k1 != 0 ){ // fMd is not changed
    fBd.resize(2*fNp,fMd); 
    col_s = 0;
  }
  else if( k1 == 0 ){ // fMd is changed
    fMd -= 2;
    fBd.resize(2*fNp,fMd); 
    col_s = 2;
  }

  fBd.row(ix(fi[1])-2) = fBv.block(ix(1),col_s,1,fMd);  
  fBd.row(iy(fi[1])-2) = fBv.block(iy(1),col_s,1,fMd); 

  fBs.resize(2*fNp,2*fNp);  
  fBs.setZero();
  fTripletList.clear();
  fMs = 0;

  Set_Bd_S_Comp(0, col_s); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

void TOperator_E::Cal_u_Comp_Sfb( int k1 ) 
{
  VectorXd w(2*fNp); 
  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fNp); 
  static double eval;
  static double eval_ww;
  double d;
  int st;

  // forward
  for( int i = 0; i < fN; ++i ){
    st = i;
    for( int k = 0; k < fNp; ++k ){
      w(ix(k))=fW(ix((k+st+1)%fN));  
      w(iy(k))=fW(iy((k+st+1)%fN));  
    }
    eval_ww = w.dot(w); 

    gzai.head(fMd) = fBd.transpose() * w; 
    gzai.tail(fMs) = fBs.leftCols(fMs).transpose() * w;

    assert( fNp == k1+1 );
    eval = -gzai.dot(gzai) + eval_ww;
    if( eval < fEval_Comp_best_Sf[ st ][ k1 ] )
      fEval_Comp_best_Sf[ st ][ k1 ] = eval;
  }

  // backward
  int index; 
  for( int i = 0; i < fN; ++i ){
    st = i;
    for( int k = 0; k < fNp; ++k ){
      index = (-k+st-1);
      if( index < 0 )
	index += fN;
      w(ix(k))=fW(ix(index));  
      w(iy(k))=fW(iy(index));  
    }
    eval_ww = w.dot(w); 

    gzai.head(fMd) = fBd.transpose() * w; 
    gzai.tail(fMs) = fBs.leftCols(fMs).transpose() * w;

    assert( fNp == k1+1 );
    eval = -gzai.dot(gzai) + eval_ww;
    if( eval < fEval_Comp_best_Sb[ st ][ k1 ] )
      fEval_Comp_best_Sb[ st ][ k1 ] = eval;
  }

}


///////////////// Partial ////////////////////
void TOperator_E::Cal_u_Part() 
{
  VectorXd gzai(fMd_part+fMs_part); 
  static double eval;
  static double eval_ww;
  int st1;

  for( int s = 0; s < fN; ++s ){
    st1 = s;
    eval_ww = fWj[st1].head(2*fNp).dot(fWj[st1].head(2*fNp)); 

    gzai.head(fMd_part) = fBd_part.transpose() * fWj[st1].head(2*fNp); 
    gzai.tail(fMs_part) = fBs_part.leftCols(fMs_part).transpose() * fWj[st1].head(2*fNp);
    eval = -gzai.dot(gzai) + eval_ww;

    fEval_part[st1] = eval; 
  }
}

void TOperator_E::Cal_u_Procrustes_Part() 
{
  VectorXd gzai(fMd_part+fMs_part);
  VectorXd p(fMd_part+fMs_part);
  VectorXd p_d(fMd_part+fMs_part);
  static double eval;
  static double eval_ww;
  int st1;
  double lambda;
  int ne;
  static int count_imp = 0;
  VectorXd x_unit(fMd_part+fMs_part); 
  VectorXd y_unit(fMd_part+fMs_part); 
  static VectorXd p2(2);
  static VectorXd p_d2(2); 
  static MatrixXd A2(2,2);
  double a_tr;

  for( int s = 0; s < fN; ++s ){
    st1 = s;
    eval_ww = fWj[st1].head(2*fNp).dot(fWj[st1].head(2*fNp)); 

    p.head(fMd_part) = fBd_part.transpose()*fWj[st1].head(2*fNp);
    p.tail(fMs_part) = fBs_part.leftCols(fMs_part).transpose()*fWj[st1].head(2*fNp);
    p_d.head(fMd_part) = fBd_part.transpose()*fWj_d[st1].head(2*fNp);
    p_d.tail(fMs_part) = fBs_part.leftCols(fMs_part).transpose()*fWj_d[st1].head(2*fNp);
    x_unit = p.normalized();
    y_unit = p_d - p_d.dot(x_unit)*x_unit;
    y_unit = y_unit.normalized();
    p2(0) = p.norm(); 
    p2(1) = 0;        
    p_d2(0) = p_d.dot(x_unit);
    p_d2(1) = p_d.dot(y_unit);
    A2 = p2*p2.transpose() + p_d2*p_d2.transpose();
    a_tr = A2(0,0)+A2(1,1);
    lambda = 0.5*(a_tr+sqrt(a_tr*a_tr - 4.0*(A2(0,0)*A2(1,1)-A2(0,1)*A2(1,0))));

    eval = eval_ww - lambda;
    fEval_part[st1] = eval; 
  }
}

void TOperator_E::Orthogonalization_Bd_Part() 
{
  VectorXd t(2*fNp); 

  for( int i = 0; i < fMd_part; ++i ){
    t = fBd_part.col(i);
    for( int j = 0; j < i; ++j ){
      t -= t.dot(fBd_part.col(j))*fBd_part.col(j);
    }
    if( t.norm() < 0.0000001 )
      printf( "Warning: not independent!! \n" );
    fBd_part.col(i) = t.normalized();
  } 
}

void TOperator_E::Set_Be_Jd_Parallel_Part( int ivs1, int ivs2 ) 
{
  VectorXd t(4*fNvp); 
  int ive2;

  assert( 0 <= ivs1 );  assert( ivs1 < fNvp );
  assert( 0 <= ivs2 );  assert( ivs2 < fNvp );
  assert( ivs1 < ivs2 );

  ive2 = ivs2 + 1;
  assert( ive2 < fNvp );

  fBe_part.row(sx_f(ivs1)) = fBv_part.row(ix(ivs1));  
  fBe_part.row(sx_b(ive2)) = fBv_part.row(ix(ive2));  
  t = VectorXd::Zero(4*fNvp);
  t(sx_f(ivs1)) = 1;
  t(sx_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;

  fBe_part.row(sy_f(ivs1)) = fBv_part.row(iy(ivs1));  
  fBe_part.row(sy_b(ive2)) = fBv_part.row(iy(ive2));  
  t = VectorXd::Zero(4*fNvp);
  t(sy_f(ivs1)) = 1;
  t(sy_b(ive2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;
}

void TOperator_E::Set_Be_S_Part( int ivs )
{
  VectorXd t(4*fNvp); 
  int ive;

  assert( 0 <= ivs );  assert( ivs < fNvp );

  assert( ivs < fNvp - 1 ); 
  ive = ivs+1;

  fBe_part.row(sx_f(ivs)) = fBv_part.row(ix(ivs));  
  fBe_part.row(sx_b(ive)) = fBv_part.row(ix(ive));  
  t = VectorXd::Zero(4*fNvp);
  t(sx_f(ivs)) = 1;
  t(sx_b(ive)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;

  fBe_part.row(sy_f(ivs)) = fBv_part.row(iy(ivs));  
  fBe_part.row(sy_b(ive)) = fBv_part.row(iy(ive));  
  t = VectorXd::Zero(4*fNvp);
  t(sy_f(ivs)) = 1;
  t(sy_b(ive)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;
}

void TOperator_E::Set_Be_Js_Xsymmetry_Part( int ivs1, int ivs2 ) 
{
  VectorXd t(4*fNvp); 

  assert( 0 <= ivs1 );  assert( ivs1 < fNvp );
  assert( 0 <= ivs2 );  assert( ivs2 < fNvp );
  assert( ivs1 < ivs2 );

  fBe_part.row(sx_f(ivs1)) = fBv_part.row(ix(ivs1));  
  fBe_part.row(sx_f(ivs2)) = fBv_part.row(ix(ivs2));  
  t = VectorXd::Zero(4*fNvp);
  t(sx_f(ivs1)) = 1;
  t(sx_f(ivs2)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;

  fBe_part.row(sy_f(ivs1)) = fBv_part.row(iy(ivs1));  
  fBe_part.row(sy_f(ivs2)) = fBv_part.row(iy(ivs2));  
  t = VectorXd::Zero(4*fNvp);
  t(sy_f(ivs1)) = 1;
  t(sy_f(ivs2)) = -1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;
}

void TOperator_E::Set_Be_J_Part( int ivs ) 
{
  VectorXd t(4*fNvp); 

  assert( 0 <= ivs );  assert( ivs < fNv );

  fBe_part.row(sx_f(ivs)) = fBv_part.row(ix(ivs));  
  t = VectorXd::Zero(4*fNvp);
  t(sx_f(ivs)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;

  fBe_part.row(sy_f(ivs)) = fBv_part.row(iy(ivs));  
  t = VectorXd::Zero(4*fNvp);
  t(sy_f(ivs)) = 1;
  t = t.normalized();
  for( int i = 0; i < fMd_part; ++i )
    fBe_part.col(i) -= fBe_part.col(i).dot(t)*t;
}

void TOperator_E::Set_Bd_Jd_Parallel_Part( int ivs1, int ivs2 ) 
{
  int n, is1, ie2;
  int ive2;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs1 );  assert( ivs1 < fNvp );
  assert( 0 <= ivs2 );  assert( ivs2 < fNvp );
  assert( ivs1 < ivs2 );
  
  ive2 = ivs2 + 1;
  assert( ive2 < fNvp );
  
  assert( fn[ivs1] == fn[ivs2] );
  n = fn[ivs1];
  is1 = fi[ivs1];
  ie2 = fi[ivs2+1]; 

  for( int i = 1; i <= n; ++i ){ 
    fBd_part.row(ix(is1+i)) = fBe_part.block(sx_f(ivs1),0,1,fMd_part); 
    fTripletList.push_back(T(ix(is1+i),fMs_part,inv_sq2));
    fBd_part.row(ix(ie2-i)) = fBe_part.block(sx_b(ive2),0,1,fMd_part); 
    fTripletList.push_back(T(ix(ie2-i),fMs_part,inv_sq2)); 
    ++fMs_part;

    fBd_part.row(iy(is1+i)) = fBe_part.block(sy_f(ivs1),0,1,fMd_part);
    fTripletList.push_back(T(iy(is1+i),fMs_part,inv_sq2)); 
    fBd_part.row(iy(ie2-i)) = fBe_part.block(sy_b(ive2),0,1,fMd_part);
    fTripletList.push_back(T(iy(ie2-i),fMs_part,inv_sq2));
    ++fMs_part;
  }
}

void TOperator_E::Set_Bd_S_Part( int ivs ) 
{
  int n, ive, is, ie;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs ); assert( ivs < fNvp ); 
  
  ive = ivs+1;
  assert( ive < fNvp );

  n = fn[ivs];
  is = fi[ivs];
  ie = fi[ivs+1];

  for( int i = 1; i <= n/2; ++i ){
    fBd_part.row(ix(is+i)) = fBe_part.block(sx_f(ivs),0,1,fMd_part); 
    fTripletList.push_back(T(ix(is+i),fMs_part,inv_sq2));
    fBd_part.row(ix(ie-i)) = fBe_part.block(sx_b(ive),0,1,fMd_part); 
    fTripletList.push_back(T(ix(ie-i),fMs_part,-inv_sq2));
    ++fMs_part;

    fBd_part.row(iy(is+i)) = fBe_part.block(sy_f(ivs),0,1,fMd_part); 
    fTripletList.push_back(T(iy(is+i),fMs_part,inv_sq2));
    fBd_part.row(iy(ie-i)) = fBe_part.block(sy_b(ive),0,1,fMd_part); 
    fTripletList.push_back(T(iy(ie-i),fMs_part,-inv_sq2));
    ++fMs_part;
  }
  if( n % 2 == 1 ){
    fBd_part.row(ix(is+(n+1)/2)) = 0.5*(fBv_part.block(ix(ivs),0,1,fMd_part)+fBv_part.block(ix(ive),0,1,fMd_part));
    fBd_part.row(iy(is+(n+1)/2)) = 0.5*(fBv_part.block(iy(ivs),0,1,fMd_part)+fBv_part.block(iy(ive),0,1,fMd_part));
  }
}

void TOperator_E::Set_Bd_Js_Xsymmetry_Part( int ivs1, int ivs2 ) 
{
  int n, is1, is2;
  const double inv_sq2 = 1.0/sqrt(2.0);

  assert( 0 <= ivs1 );  assert( ivs1 < fNvp );
  assert( 0 <= ivs2 );  assert( ivs2 < fNvp );
  assert( ivs1 < ivs2 );

  n = fn[ivs1];
  is1 = fi[ivs1];
  is2 = fi[ivs2];

  for( int i = 1; i <= n; ++i ){ 
    fBd_part.row(ix(is1+i)) = fBe_part.block(sx_f(ivs1),0,1,fMd_part); 
    fTripletList.push_back(T(ix(is1+i),fMs_part,inv_sq2));
    fBd_part.row(ix(is2+i)) = fBe_part.block(sx_f(ivs2),0,1,fMd_part); 
    fTripletList.push_back(T(ix(is2+i),fMs_part,inv_sq2)); 
    ++fMs_part;

    fBd_part.row(iy(is1+i)) = fBe_part.block(sy_f(ivs1),0,1,fMd_part);
    fTripletList.push_back(T(iy(is1+i),fMs_part,inv_sq2)); 
    fBd_part.row(iy(is2+i)) = fBe_part.block(sy_f(ivs2),0,1,fMd_part);
    fTripletList.push_back(T(iy(is2+i),fMs_part,-inv_sq2));
    ++fMs_part;
  }
}

void TOperator_E::Set_Bd_J_Part( int ivs ) 
{
  int is;
  int n;
  n = fn[ivs];
  is = fi[ivs];

  for( int i = 1; i <= n; ++i ){ 
    fBd_part.row(ix(is+i)) = fBe_part.block(sx_f(ivs),0,1,fMd_part); 
    fTripletList.push_back(T(ix(is+i),fMs_part,1));
    ++fMs_part;
    fBd_part.row(iy(is+i)) = fBe_part.block(sy_f(ivs),0,1,fMd_part);
    fTripletList.push_back(T(iy(is+i),fMs_part,1)); 
    ++fMs_part;
  }
}

void TOperator_E::Make_Bv_Be_IH4_Part() 
{
  int m;

  fNvp = 5; 
  int m_max = 2;
  MatrixXd As = MatrixXd::Zero(m_max,2*fNvp);
  m = 0;
  int ivs, ive;
  ivs = 0; ive =3;
  // Set_As_Parallel(0,3,As,m); 
  As(m,ix(ivs+1)) = 1;
  As(m,ix(ivs)) = -1;
  As(m,ix(ive+1)) = 1;
  As(m,ix(ive)) = -1;
  ++m;
  As(m,iy(ivs+1)) = 1;
  As(m,iy(ivs)) = -1;
  As(m,iy(ive+1)) = 1;
  As(m,iy(ive)) = -1;
  ++m;
  assert( m_max == m );

  fBv_part = this->Trans_AtoB( As );
  fMd_part = fBv_part.cols();
  assert( fBv_part.rows() == 2*fNvp );
  assert( fBv_part.cols() == 2*fNvp-m_max );

  fBe_part.resize(4*fNvp,fMd_part);
  fBe_part = MatrixXd::Zero(4*fNvp,fMd_part);

  Set_Be_Jd_Parallel_Part(0,3); 
  Set_Be_S_Part(1); 
  Set_Be_S_Part(2); 
}

void TOperator_E::Make_B_IH4_Part( int k1, int k2, int k3 ) 
{
  fn[0] = k1; // IH
  fn[1] = k2; 
  fn[2] = k3; 
  fn[3] = k1; 

  fi[0] = 0;
  for( int s = 1; s <= fNvp; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNvp-1]+1;

  fBd_part.resize(2*fNp,fMd_part); 

  for( int s = 0; s < fNvp; ++s ){
    fBd_part.row(ix(fi[s])) = fBv_part.block(ix(s),0,1,fMd_part); 
    fBd_part.row(iy(fi[s])) = fBv_part.block(iy(s),0,1,fMd_part); 
  }

  fBs_part.resize(2*fNp,2*fNp);
  fBs_part.setZero();
  fTripletList.clear();
  fMs_part = 0;
  
  Set_Bd_Jd_Parallel_Part(0,3); 
  Set_Bd_S_Part(1); 
  Set_Bd_S_Part(2); 

  fBs_part.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd_Part();
}

void TOperator_E::Make_Bv_Be_IH5_Part() 
{
  int m;

  fNvp = 5; 
  int m_max = 4; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNvp);
  m = 0;
  int ivs, ive;
  // Set_As_Parallel(0,3,As,m);  
  ivs = 0; ive = 3; 
  As(m,ix(ivs+1)) = 1;
  As(m,ix(ivs)) = -1;
  As(m,ix(ive+1)) = 1;
  As(m,ix(ive)) = -1;
  ++m;
  As(m,iy(ivs+1)) = 1;
  As(m,iy(ivs)) = -1;
  As(m,iy(ive+1)) = 1;
  As(m,iy(ive)) = -1;
  ++m;

  As(m,iy(1)) = 1; 
  As(m,iy(3)) = -1;
  ++m;
  As(m,ix(2)) = 1;
  As(m,ix(1)) = -0.5;
  As(m,ix(3)) = -0.5;
  ++m;

  assert( m_max == m );

  fBv_part = this->Trans_AtoB( As );
  fMd_part = fBv_part.cols();
  assert( fBv_part.rows() == 2*fNvp );
  assert( fBv_part.cols() == 2*fNvp-m_max );

  fBe_part.resize(4*fNvp,fMd_part);
  fBe_part = MatrixXd::Zero(4*fNvp,fMd_part);

  Set_Be_Jd_Parallel_Part(0,3); 
  Set_Be_Js_Xsymmetry_Part(1,2); 
}

void TOperator_E::Make_B_IH5_Part( int k1, int k2 ) 
{
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k2; 
  fn[3] = k1; 

  fi[0] = 0;
  for( int s = 1; s <= fNvp; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;
  
  fNp = fi[fNvp-1]+1;

  fBd_part.resize(2*fNp,fMd_part); 
  for( int s = 0; s < fNvp; ++s ){
    fBd_part.row(ix(fi[s])) = fBv_part.block(ix(s),0,1,fMd_part); 
    fBd_part.row(iy(fi[s])) = fBv_part.block(iy(s),0,1,fMd_part); 
  }
  
  fBs_part.resize(2*fNp,2*fNp);
  fBs_part.setZero();
  fTripletList.clear();
  fMs_part = 0;
  
  Set_Bd_Jd_Parallel_Part(0,3); 
  Set_Bd_Js_Xsymmetry_Part(1,2); 

  fBs_part.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd_Part();
}


void TOperator_E::Make_Bv_Be_IH6_Part() 
{
  int m;

  fNvp = 4; 
  int m_max = 2; 
  MatrixXd As = MatrixXd::Zero(m_max,2*fNvp);
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

  assert( m_max == m );

  fBv_part = this->Trans_AtoB( As );
  fMd_part = fBv_part.cols();
  assert( fBv_part.rows() == 2*fNvp );
  assert( fBv_part.cols() == 2*fNvp-m_max );

  fBe_part.resize(4*fNvp,fMd_part);
  fBe_part = MatrixXd::Zero(4*fNvp,fMd_part);

  Set_Be_Js_Xsymmetry_Part(0,2);
  Set_Be_J_Part(1); 
}


void TOperator_E::Make_B_IH6_Part( int k1, int k2 ) 
{
  fn[0] = k1; 
  fn[1] = k2; 
  fn[2] = k1; 

  fi[0] = 0;
  for( int s = 1; s <= fNvp; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNvp-1]+1;

  fBd_part.resize(2*fNp,fMd_part); 

  for( int s = 0; s < fNvp; ++s ){
    fBd_part.row(ix(fi[s])) = fBv_part.block(ix(s),0,1,fMd_part); 
    fBd_part.row(iy(fi[s])) = fBv_part.block(iy(s),0,1,fMd_part); 
  }

  fBs_part.resize(2*fNp,2*fNp);
  fBs_part.setZero();
  fTripletList.clear();
  fMs_part = 0;

  Set_Bd_Js_Xsymmetry_Part(0,2); 
  Set_Bd_J_Part(1); 

  fBs_part.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd_Part();
}
