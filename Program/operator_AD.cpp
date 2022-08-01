/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_AD__
#include "operator_AD.h"
#endif

TOperator_AD::TOperator_AD( int N ) : TOperator_E( N )  
{
  fWaj = new VectorXd [fN]; 
  fWaj_d = new VectorXd [fN]; 
  for( int i = 0; i < fN; ++i ){
    fWaj[i] = VectorXd::Zero(2*fN); 
    fWaj_d[i] = VectorXd::Zero(2*fN); 
  }

  fBa = MatrixXd(2*fN,2*fN); 
} 

TOperator_AD::~TOperator_AD() 
{

} 


// SetParameter(): defined in TOperator_E 
// Define(): defined in TOperator_base


void TOperator_AD::SetInit( VectorXd w )
{
  TOperator_E::SetInit( w );

  for( int st1 = 0; st1 < fN; ++st1 ){ 
    for( int k = 0; k < fN; ++k ){
      fWaj[st1](ix(k)) = fWj[st1](ix(k+1)) - fWj[st1](ix(k));
      fWaj[st1](iy(k)) = fWj[st1](iy(k+1)) - fWj[st1](iy(k));
      fWaj_d[st1](ix(k)) = fWj[st1](iy(k+1)) - fWj[st1](iy(k));
      fWaj_d[st1](iy(k)) = -(fWj[st1](ix(k+1)) - fWj[st1](ix(k)));
    }
  }
}


void TOperator_AD::Cal_u( int k1, int k2, int k3, int k4 ) 
{
  int st1, st2;
  int numOfCheckS;
  int checkS[fN];

  if( fIH == 4 && fFlagPrune ){ 
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

    this->Make_B_IH4( k1, k2, k3, k4 );  // (IH4)
    this->Make_Ba(); // adj
  }
  else{
    numOfCheckS = fN;
    for( int s = 0; s < fN; ++s )
      checkS[s] = s;
    this->Make_Ba();
  }

  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fN);
  double eval;
  VectorXd tmp(2*fN);  

  for( int s = 0; s < numOfCheckS; ++s ){
    st1 = checkS[s];

    gzai.head(fMa) = fBa.transpose() * fWaj[st1]; 
    eval = -gzai.head(fMa).dot(gzai.head(fMa)) + fWaj[st1].transpose()*fWaj[st1];

    if( eval < fEval_best ){
      tmp = fBa*gzai.head(fMa); // adj
      u(ix(0)) = 0.0;  // adj 
      u(iy(0))= 0.0;   // adj
      for( int i = 1; i < fN; ++i ){ // adj
	u(ix(i)) = u(ix(i-1)) + tmp(ix(i-1));
	u(iy(i)) = u(iy(i-1)) + tmp(iy(i-1));
      }

      if( fFlagCheckIntersection == 0 || Check_Intersection( u ) == false ){
	this->Trans_U_center( u );
	this->Update_U_top( eval, u, st1 );
      }
    }
  }
}

void TOperator_AD::Cal_u_Procrustes( int k1, int k2, int k3, int k4 ) 
{
  int st1, st2;
  int numOfCheckS;
  int checkS[fN];

  if( fIH == 5 && fFlagPrune ){ // (IH5)
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
    this->Make_Ba(); // adj
  }
  else if( fIH == 6 && fFlagPrune ){ // (IH6)
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
    this->Make_Ba(); // adj
  }
  else{
    numOfCheckS = fN;
    for( int s = 0; s < fN; ++s )
      checkS[s] = s;
    this->Make_Ba(); // adj
  }

  VectorXd gzai(fMd+fMs);
  VectorXd ur(2*fN);
  VectorXd u(2*fN);
  VectorXd p(fMd+fMs);
  VectorXd p_d(fMd+fMs);
  double eval;
  double tr, rt, cos_th, sin_th;
  double lambda;
  int ne;
  VectorXd x_unit(fMd+fMs); 
  VectorXd y_unit(fMd+fMs); 
  VectorXd p2(2);
  VectorXd p_d2(2); 
  MatrixXd A2(2,2);
  double a_tr;
  VectorXd tmp(2*fN);  

  for( int s = 0; s < numOfCheckS; ++s ){
    st1 = checkS[s];

    p.head(fMa) = fBa.transpose()*fWaj[st1]; // adj
    p_d.head(fMa) = fBa.transpose()*fWaj_d[st1]; // adj 
    x_unit.head(fMa) = p.head(fMa).normalized();
    y_unit.head(fMa) = p_d.head(fMa) - p_d.head(fMa).dot(x_unit.head(fMa))*x_unit.head(fMa);
    y_unit.head(fMa) = y_unit.head(fMa).normalized();
    p2(0) = p.head(fMa).norm(); 
    p2(1) = 0;        
    p_d2(0) = p_d.head(fMa).dot(x_unit.head(fMa));
    p_d2(1) = p_d.head(fMa).dot(y_unit.head(fMa));
    A2 = p2*p2.transpose() + p_d2*p_d2.transpose();
    a_tr = A2(0,0)+A2(1,1);
    lambda = 0.5*(a_tr+sqrt(a_tr*a_tr - 4.0*(A2(0,0)*A2(1,1)-A2(0,1)*A2(1,0))));
    eval = fWaj[st1].transpose()*fWaj[st1] - lambda; // adj

    if( eval < fEval_best){
      gzai.head(fMa) = A2(0,1)*x_unit.head(fMa) + (lambda-A2(0,0))*y_unit.head(fMa);
      gzai.head(fMa) = gzai.head(fMa).normalized()*sqrt(lambda);
      
      tmp = fBa*gzai.head(fMa); // adj
      u(ix(0)) = 0.0;  // adj
      u(iy(0)) = 0.0;  // adj
      for( int i = 1; i < fN; ++i ){  // adj
	u(ix(i)) = u(ix(i-1)) + tmp(ix(i-1));
	u(iy(i)) = u(iy(i-1)) + tmp(iy(i-1));
      }

      if( fFlagCheckIntersection == 0 || Check_Intersection( u ) == false ){
	tr = tmp.transpose()*fWaj[st1]; 
	rt = tmp.transpose()*fWaj_d[st1]; 
	cos_th = tr/sqrt(tr*tr + rt*rt);
	sin_th = rt/sqrt(tr*tr + rt*rt);
	for( int k = 0; k < fN; ++k ){
	  ur(ix(k)) = cos_th*u(ix(k)) - sin_th*u(iy(k)); 
	  ur(iy(k)) = sin_th*u(ix(k)) + cos_th*u(iy(k)); 
	}

	this->Trans_U_center( ur );
	this->Update_U_top( eval, ur, st1 );
      }
    }
  }
}


void TOperator_AD::Make_Ba() 
{
  VectorXd t(2*fN); 
  MatrixXd B_tmp(2*fN,fMd+fMs);  
  MatrixXd B_tmp2(2*fN,fMd+fMs);  
  double value;
  VectorXi ZeroIndex = VectorXi::Zero(fMd+fMs);
  VectorXd BGB(fMd+fMs);  

  fMa = fMd+fMs; 
  B_tmp.leftCols(fMd) = fBd.leftCols(fMd); 
  B_tmp.rightCols(fMs) = fBs.leftCols(fMs); 

  for( int i = 0; i < fN; ++i ){
    B_tmp2.row(ix(i)) = B_tmp.row(ix(i+1)) - B_tmp.row(ix(i));
    B_tmp2.row(iy(i)) = B_tmp.row(iy(i+1)) - B_tmp.row(iy(i));
  }

  fMa = 0;
  for( int i = 0; i < fMd+fMs; ++i ){
    if( i != 0 )
      BGB.head(i) = B_tmp2.col(i).transpose()*B_tmp2.leftCols(i);
    for( int j = 0; j < i; ++j ){
      if( fabs(BGB(j)) > 0.00000001 ) 
	B_tmp2.col(i) -= BGB(j) * B_tmp2.col(j);
    }

    value = double(B_tmp2.col(i).transpose()*B_tmp2.col(i));
    if( fabs(value) > 0.000001 ){ 
      B_tmp2.col(i) *= 1.0/sqrt(value);
      ZeroIndex(i) = 1;
      ++fMa;
    }
  } 

  fBa.resize(2*fN,fMa); 
  int k = 0;
  for( int i = 0; i < fMd+fMs; ++i ){
    if( ZeroIndex(i) == 1 )
      fBa.col(k++) = B_tmp2.col(i);    
  }
  assert( k == fMa );
}

void TOperator_AD::Make_B_Comp_SS( int k1, int k2 ) 
{
  const double inv_sq2 = 1.0/sqrt(2.0);
  
  fn[0] = k1; 
  fn[1] = k2; 

  fi[0] = 0;
  for( int s = 1; s <= fNv-1; ++s )  
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNv-1]+1; 

  fMd = fMv; // comp
  fBd.resize(2*fNp,fMd); 

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd); 
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.resize(2*fNp,2*fNp);
  fBs.setZero();
  fTripletList.clear();
  fMs = 0;
  
  Set_Bd_S(0); 
  Set_Bd_S(1); 

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}

void TOperator_AD::Cal_u_Comp_SS( int k1, int k2 ) 
{
  VectorXd w(2*fNp); 
  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fNp); 
  static double eval;
  static double eval_ww;
  double d;
  int st;

  this->Make_Ba_Comp();
  for( int i = 0; i < fN; ++i ){
    st = i;

    eval_ww = fWaj[st].head(2*(fNp-1)).dot(fWaj[st].head(2*(fNp-1))); 
    gzai.head(fMa) = fBa.transpose() * fWaj[st].head(2*(fNp-1)); 
    eval = -gzai.head(fMa).dot(gzai.head(fMa)) + eval_ww;

    assert( fNp == k1+k2+3 );
 
    fEval_Comp_best[ st ][ k1 ][ k2 ] = eval;
  }
}

void TOperator_AD::Make_Ba_Comp() 
{
  VectorXd t(2*(fNp-1)); 
  MatrixXd B_tmp(2*fNp,fMd+fMs);  
  MatrixXd B_tmp2(2*(fNp-1),fMd+fMs);  
  double value;

  fMa = fMd+fMs; 
  B_tmp.leftCols(fMd) = fBd.leftCols(fMd); 
  B_tmp.rightCols(fMs) = fBs.leftCols(fMs); 

  for( int i = 0; i < fNp-1; ++i ){
    assert( i < fN-1 );
    B_tmp2.row(ix(i)) = B_tmp.row(ix(i+1)) - B_tmp.row(ix(i));
    B_tmp2.row(iy(i)) = B_tmp.row(iy(i+1)) - B_tmp.row(iy(i));
  }

  for( int i = 0; i < fMa; ++i )
    B_tmp2.col(i).normalize();

  for( int i = 0; i < fMa; ++i ){
    t = B_tmp2.col(i);
    for( int j = 0; j < i; ++j ){
      t -= (t.transpose() * B_tmp2.col(j)) * B_tmp2.col(j); 
    }
    if( fabs(t.transpose()*t) < 0.000001 ){ 
      // printf( "Warning: not independent!!\n" );
      for( int h = i; h < fMa-1; ++h ){
	B_tmp2.col(h) = B_tmp2.col(h+1);
      }
      --fMa;
      --i;
    }
    else {
      value = double(t.transpose()*t);
      assert( value > 0.000001 ); 
      B_tmp2.col(i) = 1.0/sqrt(value) * t; 
    }
  } 

  fBa.resize(2*(fNp-1),fMa); 
  fBa = B_tmp2.leftCols(fMa);
}

void TOperator_AD::Make_B_Comp_S( int k1 ) 
{
  const double inv_sq2 = 1.0/sqrt(2.0);
  int col_s;
  
  fn[0] = k1; 

  fi[0] = 0;
  for( int s = 1; s <= fNv-1; ++s )
    fi[s] = fi[s-1] + fn[s-1] + 1;

  fNp = fi[fNv-1]+1; // adj

  fMd = fMv; 
  fBd.resize(2*fNp,fMd); 

  for( int s = 0; s < fNv; ++s ){
    fBd.row(ix(fi[s])) = fBv.block(ix(s),0,1,fMd);  // adj
    fBd.row(iy(fi[s])) = fBv.block(iy(s),0,1,fMd); 
  }

  fBs.resize(2*fNp,2*fNp);  
  fBs.setZero();
  fTripletList.clear();
  fMs = 0;

  Set_Bd_S(0);  // adj

  fBs.setFromTriplets(fTripletList.begin(), fTripletList.end());

  this->Orthogonalization_Bd();
}


void TOperator_AD::Cal_u_Comp_Sfb( int k1 ) 
{
  VectorXd w(2*(fNp-1));  // adj
  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fNp); 
  static double eval;
  static double eval_ww;
  double d;
  int st;

  this->Make_Ba_Comp(); // adj
  // forward
  for( int i = 0; i < fN; ++i ){
    st = i;
    eval_ww = fWaj[st].head(2*(fNp-1)).dot(fWaj[st].head(2*(fNp-1))); // adj
    gzai.head(fMa) = fBa.transpose() * fWaj[st].head(2*(fNp-1));  // adj
    eval = -gzai.head(fMa).dot(gzai.head(fMa)) + eval_ww; // adj

    assert( fNp == k1+2 );
    if( eval < fEval_Comp_best_Sf[ st ][ k1 ] )
      fEval_Comp_best_Sf[ st ][ k1 ] = eval;
  }

  // backward
  int index; 
  for( int i = 0; i < fN; ++i ){
    st = i;
    for( int k = 0; k < fNp-1; ++k ){ // adj
      index = fN-1-k; // adj
      w(ix(k))=fWaj[st](ix(index));  // adj
      w(iy(k))=fWaj[st](iy(index));  // adj
    }
    eval_ww = w.dot(w); 
    gzai.head(fMa) = fBa.transpose() * w;  // adj
    eval = -gzai.head(fMa).dot(gzai.head(fMa)) + eval_ww; // adj

    assert( fNp == k1+2 );
    if( eval < fEval_Comp_best_Sb[ st ][ k1 ] )
      fEval_Comp_best_Sb[ st ][ k1 ] = eval;
  }
}


void TOperator_AD::Cal_u_Part() 
{
  this->Make_Ba_Part(); // adj

  VectorXd gzai(fMa); 
  static double eval;
  static double eval_ww;
  int st1;

  for( int s = 0; s < fN; ++s ){
    st1 = s;
    eval_ww = fWaj[st1].head(2*(fNp-1)).dot(fWaj[st1].head(2*(fNp-1))); // adj

    gzai = fBa.transpose() * fWaj[st1].head(2*(fNp-1)); // adj
    eval = -gzai.dot(gzai) + eval_ww; // adj
    assert( -9999999.9 < eval || eval < 99999999.9 );

    fEval_part[st1] = eval; 
  }
}


void TOperator_AD::Cal_u_Procrustes_Part() 
{
  this->Make_Ba_Part(); // adj

  VectorXd gzai(fMa);
  VectorXd p(fMa);
  VectorXd p_d(fMa);
  static double eval;
  static double eval_ww;
  int st1;
  double lambda;
  int ne;
  static int count_imp = 0;
  VectorXd x_unit(fMa); 
  VectorXd y_unit(fMa); 
  static VectorXd p2(2);
  static VectorXd p_d2(2); 
  static MatrixXd A2(2,2);
  double a_tr;
  
  for( int s = 0; s < fN; ++s ){
    st1 = s;
    eval_ww = fWaj[st1].head(2*(fNp-1)).dot(fWaj[st1].head(2*(fNp-1))); // adj

    p = fBa.transpose()*fWaj[st1].head(2*(fNp-1)); // adj
    p_d = fBa.transpose()*fWaj_d[st1].head(2*(fNp-1)); // adj
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
    assert( -9999999.9 < eval || eval < 99999999.9 );
  }
}

void TOperator_AD::Make_Ba_Part() 
{
  VectorXd t(2*(fNp-1)); 
  MatrixXd B_tmp(2*fNp,fMd_part+fMs_part);  
  MatrixXd B_tmp2(2*(fNp-1),fMd_part+fMs_part);  
  double value;

  fMa = fMd_part+fMs_part; 
  B_tmp.leftCols(fMd_part) = fBd_part.leftCols(fMd_part); 
  B_tmp.rightCols(fMs_part) = fBs_part.leftCols(fMs_part); 

  for( int i = 0; i < fNp-1; ++i ){
    assert( i < fN-1 );
    B_tmp2.row(ix(i)) = B_tmp.row(ix(i+1)) - B_tmp.row(ix(i));
    B_tmp2.row(iy(i)) = B_tmp.row(iy(i+1)) - B_tmp.row(iy(i));
  }

  for( int i = 0; i < fMa; ++i )
    B_tmp2.col(i).normalize();

  for( int i = 0; i < fMa; ++i ){
    t = B_tmp2.col(i);
    for( int j = 0; j < i; ++j ){
      t -= (t.transpose() * B_tmp2.col(j)) * B_tmp2.col(j); 
    }
    if( fabs(t.transpose()*t) < 0.000001 ){ 
      // printf( "Warning: not independent!!\n" );
      for( int h = i; h < fMa-1; ++h ){
	B_tmp2.col(h) = B_tmp2.col(h+1);
      }
      --fMa;
      --i;
    }
    else {
      value = double(t.transpose()*t);
      assert( value > 0.000001 ); 
      B_tmp2.col(i) = 1.0/sqrt(value) * t; 
    }
  } 

  fBa.resize(2*(fNp-1),fMa); 
  fBa = B_tmp2.leftCols(fMa);
}
