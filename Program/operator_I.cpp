/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_I__
#include "operator_I.h"
#endif

TOperator_I::TOperator_I( int N, int N_in ) : TOperator_base( N )
{
  fN = N;
  fN_in = N_in; 
  fNN = fN + fN_in;
  fW = VectorXd::Zero(2*fN);
  fW_d = VectorXd::Zero(2*fN);
  fW_in = VectorXd::Zero(2*fN_in); 
  fWW = VectorXd::Zero(2*fNN); 
  fWW_d = VectorXd::Zero(2*fNN); 
  fKKi = MatrixXi::Zero(fNN,fNN);
  fKK = MatrixXd::Zero(fNN,fNN); 
  fH = MatrixXd::Zero(2*fN_in, 2*fN);
  fA = MatrixXd::Zero(fNN, fNN); 
  fG2_inv = MatrixXd::Zero(2*fN_in, 2*fN_in);

  fG = SparseMatrix<double>(2*fN,2*fN);
  fG.reserve(8*fN);
  fTripletList.reserve(8*fN);
  fBdj = MatrixXd(2*fN,2*fN);
  fBsj = SparseMatrix<double>(2*fN,2*fN);
  fBsj.reserve(8*fN);
  
  fNei = MatrixXi::Zero(fNN,20);
  fNei_Size = VectorXi::Zero(fNN);
} 

TOperator_I::~TOperator_I() 
{

} 

void TOperator_I::SetParameter() 
{
  fAlpha = 1.0; // this value may be reset in TEnvironment_I::SetInit()
}

void TOperator_I::SetInit( VectorXd w, VectorXd w_in, MatrixXi kki ) 
{
  fW = w;
  for( int i = 0; i < fN; ++i ){
    fW_d(ix(i))= fW(iy(i)); 
    fW_d(iy(i))= -fW(ix(i)); 
  }

  fW_in = w_in;  
  for( int i = 0; i < fN; ++i ) {
    fWW(2*i) = fW(2*i);
    fWW(2*i+1) = fW(2*i+1);
  }
  for( int i = 0; i < fN_in; ++i ) {
    fWW(2*(fN+i)) = fW_in(2*i);
    fWW(2*(fN+i)+1) = fW_in(2*i+1);
  }

  for( int i = 0; i < fNN; ++i ) { 
    fWW_d(iix(i)) = fWW(iiy(i)); 
    fWW_d(iiy(i)) = -fWW(iix(i));
  }

  fKKi = kki;
  this->SetK(); 
}

void TOperator_I::Cal_u( int st1 ) 
{
  VectorXd gzai(fMd+fMs); 
  VectorXd u(2*fN);
  VectorXd u_in(2*fN_in);
  VectorXd uu(2*fNN); 
  double eval;

  ///////////////////////////////////////////////////////
  fMk = fMd+fMs;
  fBdj.resize(2*fN,fMd);
  fBsj.resize(2*fN,fMs);

  this->Set_Bdj_Bsj( st1 ); // Set fBdj and fBsj

  eval = this->Cal_eval_gzai( gzai ); 
  assert( eval > 0.0 );

  if( eval < fEval_best ){
    u = fBdj*gzai.head(fMd); 
    u += fBsj*gzai.tail(fMs);  

    if( Check_Intersection( u ) == false ){
      u_in = fH*(u-fW) + fW_in; 
      uu.head(2*fN) = u;        
      uu.tail(2*fN_in) = u_in;
      // this->Check_Eval( uu, eval ); // check
	
      this->Trans_UU_center( uu );     
      this->Update_UU_top( eval, uu, st1 ); 
    }
  }
}


void TOperator_I::Cal_u_Procrustes( int st1 ) 
{
  VectorXd gzai(fMd+fMs);
  VectorXd u_best(2*fN);
  VectorXd u(2*fN);
  VectorXd p(fMd+fMs);
  VectorXd p_d(fMd+fMs);
  double eval;
  double lambda;
  int ne;
  VectorXd x_unit(fMd+fMs); 
  VectorXd y_unit(fMd+fMs); 
  VectorXd p2(2);
  VectorXd p_d2(2); 
  MatrixXd A2(2,2);
  double a_tr;
  MatrixXd L(fMd+fMs,fMd+fMs);
  VectorXd u_in(2*fN_in); // Inner
  VectorXd uu(2*fNN); // Inner
  
  ///////////////////////////////////////////////////////
  fMk = fMd+fMs;
  fBdj.resize(2*fN,fMd);
  fBsj.resize(2*fN,fMs);
  static VectorXd fGW = fG*fW;
  static VectorXd fGW_d = fG*fW_d;
  double tr, rt, cos_th, sin_th; 
  VectorXd ur(2*fN);

  this->Set_Bdj_Bsj( st1 ); // Set fBdj and fBsj

  p.head(fMd) = fBdj.transpose()*fGW;
  p.tail(fMs) = fBsj.transpose()*fGW;
  p_d.head(fMd) = fBdj.transpose()*fGW_d;
  p_d.tail(fMs) = fBsj.transpose()*fGW_d;
  this->Cal_L_p( p, p_d, L, 0 );

  x_unit.head(fMk) = p.head(fMk).normalized();
  y_unit.head(fMk) = p_d.head(fMk) - p_d.head(fMk).dot(x_unit.head(fMk))*x_unit.head(fMk);
  y_unit.head(fMk) = y_unit.head(fMk).normalized();

  p2(0) = p.head(fMk).norm(); 
  p2(1) = 0;        
  p_d2(0) = p_d.head(fMk).dot(x_unit.head(fMk));
  p_d2(1) = p_d.head(fMk).dot(y_unit.head(fMk));
  A2 = p2*p2.transpose() + p_d2*p_d2.transpose();
  a_tr = A2(0,0)+A2(1,1);
  lambda = 0.5*(a_tr+sqrt(a_tr*a_tr - 4.0*(A2(0,0)*A2(1,1)-A2(0,1)*A2(1,0))));
  eval = fEval_wGw - lambda;
  assert( eval > 0.0 );
    
  gzai.head(fMk) = A2(0,1)*x_unit.head(fMk) + (lambda-A2(0,0))*y_unit.head(fMk);
  gzai.head(fMk) = gzai.head(fMk).normalized()*sqrt(lambda);
  this->Cal_L_gzai( gzai, L.transpose() );

  if( eval < fEval_best ){
    u = fBdj*gzai.head(fMd);
    u += fBsj*gzai.tail(fMs);

    if( Check_Intersection( u ) == false ){
      tr = u.transpose()*fGW;
      rt = u.transpose()*fGW_d;
      cos_th = tr/sqrt(tr*tr + rt*rt);
      sin_th = rt/sqrt(tr*tr + rt*rt);
      for( int k = 0; k < fN; ++k ){
	ur(ix(k)) = cos_th*u(ix(k)) - sin_th*u(iy(k)); 
	ur(iy(k)) = sin_th*u(ix(k)) + cos_th*u(iy(k)); 
      }
      u = ur; 
    
      u_in = fH*(u-fW) + fW_in; 
      uu.head(2*fN) = u;        
      uu.tail(2*fN_in) = u_in;
      // this->Check_Eval( uu, eval ); // check

      this->Trans_UU_center( uu );     
      this->Update_UU_top( eval, uu, st1 ); 
    }
  }
}


void TOperator_I::Check_Eval( VectorXd uu,  double eval ) 
{
  VectorXd v(2); 
  double sum = 0.0;
  for( int i = 0; i < fNN; ++i ){
    int size = fNei_Size(i);

    for( int s = 0; s < size; ++s ){
      int j = fNei(i,s); 
      double aij = fA(i,j);
      assert( i != j ); 
      v(0) = uu(iix(j)) - uu(iix(i)) - (fWW(iix(j)) - fWW(iix(i)));
      v(1) = uu(iiy(j)) - uu(iiy(i)) - (fWW(iiy(j)) - fWW(iiy(i)));
      sum += aij*(v(0)*v(0)+v(1)*v(1));
    }
  }
	
  v(0) = uu(iix(0)) - fWW(iix(0));
  v(1) = uu(iiy(0)) - fWW(iiy(0));
  sum += v(0)*v(0)+v(1)*v(1);   

  printf( "%lf %lf -> %lf\n", sum , eval, eval-sum ); 
  assert( fabs(eval-sum) < 0.000001 ); 
}


void TOperator_I::Set_Bdj_Bsj( int st1 )
{
  int i, j;
  double val;

  for( int k = 0; k < fN; ++k ){
    fBdj.row(ix(k)) = fBd.row(ix(k-st1));
    fBdj.row(iy(k)) = fBd.row(iy(k-st1));
  } 

  fTripletList.clear();
  for ( int k = 0 ; k < fBs.outerSize(); ++k ) { 
    for ( SparseMatrix<double>::InnerIterator it(fBs,k); it; ++it ){
      val = it.value();
      i = (int)it.row();
      j = (int)it.col();
      fTripletList.push_back(T((i+2*st1)%(2*fN),j,val));
    }
  }
  fBsj.setFromTriplets(fTripletList.begin(), fTripletList.end());
}


void TOperator_I::Cal_L_p( VectorXd& p, VectorXd& p_d, MatrixXd& L, int iter )
{
  MatrixXd C(fMd+fMs,fMd+fMs);
  VectorXd x(fMd+fMs);
  static LLT<MatrixXd> llt;

  if( iter == 0 ){
    C.topLeftCorner(fMd,fMd) = fBdj.transpose()*fG*fBdj; 
    C.topRightCorner(fMd,fMs) = fBdj.transpose()*fG*fBsj;
    C.bottomLeftCorner(fMs,fMd) = fBsj.transpose()*fG*fBdj;
    C.bottomRightCorner(fMs,fMs) = fBsj.transpose()*fG*fBsj;
    llt.compute(C);
  }

  L = llt.matrixL();
  MatrixXd sum(1,1);

  x(0) = p(0)/L(0,0); 
  for( int i = 1; i < fMk; ++i ){
    sum = L.block(i,0,1,i)*x.head(i);
    assert( fabs(L(i,i)) > 0.000001 );
    x(i) = (p(i)-sum(0,0)) / L(i,i); 
  }
  p = x;

  x(0) = p_d(0)/L(0,0); 
  for( int i = 1; i < fMk; ++i ){
    sum = L.block(i,0,1,i)*x.head(i);
    assert( fabs(L(i,i)) > 0.000001 );
    x(i) = (p_d(i)-sum(0,0)) / L(i,i); 
  }
  p_d = x;
}

void TOperator_I::Cal_L_gzai( VectorXd& gzai, MatrixXd Lt )
{
  VectorXd x(fMk);
  MatrixXd sum(1,1);

  x(fMk-1) = gzai(fMk-1)/Lt(fMk-1,fMk-1); 
  for( int i = 1; i < fMk; ++i ){
    sum = Lt.block(fMk-1-i,fMk-i,1,i)*x.tail(i);
    assert( fabs(Lt(fMk-1-i,fMk-1-i)) > 0.000001 );
    x(fMk-1-i) = (gzai(fMk-1-i)-sum(0,0)) / Lt(fMk-1-i,fMk-1-i); 
  }

  gzai = x;
}

double TOperator_I::Cal_eval_gzai( VectorXd& gzai ) 
{
  MatrixXd C(fMd+fMs,fMd+fMs);
  VectorXd b(fMd+fMs);
  VectorXd tmp(2*fN);

  C.topLeftCorner(fMd,fMd) = fBdj.transpose()*fG*fBdj; 
  C.topRightCorner(fMd,fMs) = fBdj.transpose()*fG*fBsj;
  C.bottomLeftCorner(fMs,fMd) = fBsj.transpose()*fG*fBdj;
  C.bottomRightCorner(fMs,fMs) = fBsj.transpose()*fG*fBsj;

  tmp = fG*fW;
  b.topRows(fMd) = fBdj.transpose()*tmp;   
  b.bottomRows(fMs) = fBsj.transpose()*tmp;

  LLT<MatrixXd> llt(C);
  gzai = llt.solve(b) ;

  double eval;
  eval = -gzai.dot(b) + fEval_wGw;
  
  return eval;
}


void TOperator_I::SetK() 
{
  /////// set fKK, etc  //////
  MatrixXd K0 = MatrixXd::Zero(fN,fN); 
  MatrixXd K1 = MatrixXd::Zero(fN,fN_in); 
  MatrixXd K2 = MatrixXd::Zero(fN_in,fN_in); 
  MatrixXd K2_inv = MatrixXd::Zero(fN_in,fN_in);
  MatrixXd H = MatrixXd::Zero(fN_in,fN); 

  for( int i = 0; i < fNN; ++i ){
    int num = 0;
    for( int j = 0; j < fNN; ++j ){
      if( fKKi(i,j) != 0 && i != j )
	fNei(i,num++) = j;
      assert( num < 20 );
    }
    fNei_Size(i) = num;
  }

  double alpha;
  fA = MatrixXd::Zero(fNN,fNN); 
  for( int i = 0; i < fNN; ++i ){
    for( int s = 0; s < fNei_Size(i); ++ s ){
      int j = fNei(i,s);
      if( i < fN ) 
	alpha = 1.0;
      else 
	alpha = fAlpha;

      fA(i,j) = alpha;
    }
  }

  fDout = MatrixXd::Zero(fNN,fNN); 
  fDin = MatrixXd::Zero(fNN,fNN); 
  for( int i = 0; i < fNN; ++i ){
    for( int s = 0; s < fNei_Size(i); ++ s ){
      int j = fNei(i,s);
      fDout(i,i) += fA(i,j);
      fDin(i,i) += fA(j,i);
    }
  }

  fKK = fDout + fDin - 2.0*fA;
  for( int i = 0; i < fNN; ++i ){
    for( int j = i+1; j < fNN; ++j ){
      double tmp = 0.5*(fKK(i,j)+fKK(j,i));
      fKK(i,j) = tmp;
      fKK(j,i) = tmp;
    }
  }
  fKK(0,0) += 1.0;

  
  K0 = fKK.topLeftCorner(fN,fN);
  K1 = fKK.topRightCorner(fN,fN_in);
  K2 = fKK.bottomRightCorner(fN_in,fN_in);
  K2_inv = K2.inverse();

  // check
  MatrixXd TT = MatrixXd::Zero(fN_in,fN_in); 
  TT = K2*K2_inv; 
  for( int i = 0; i < fN_in; ++i ){
    for( int j = 0; j < fN_in; ++j ){
      if( i == j ) assert( fabs(TT(i,j)-1.0) < 0.000001 );
      else assert( fabs(TT(i,j)) < 0.000001 );
    }
  }

  MatrixXd K = MatrixXd::Zero(fN,fN);
  K = K0 - K1*K2_inv*K1.transpose();
  // check
  for( int i = 0; i < fN; ++i )
    for( int j = 0; j < fN; ++j )
      assert( fabs(K(i,j)-K(j,i)) < 0.000001 );

  /////// set fG and fH //////
  MatrixXd fg;
  fg = MatrixXd::Zero(2*fN,2*fN);
  for( int i = 0; i < fN; ++i ){
    for( int j = 0; j < fN; ++j ){
      fg(ix(i),ix(j)) = K(i,j);
      fg(iy(i),iy(j)) = K(i,j);
    }
  }

  double val;
  fTripletList.clear();
  for( int i = 0; i < 2*fN; ++i ){
    for( int j = 0; j < 2*fN; ++j ){
      val = fg(i,j); 
      if( fabs(val) > 0.000000001 ) 
	fTripletList.push_back(T(i,j,val));
    }
  }
  fG.setFromTriplets(fTripletList.begin(), fTripletList.end());

  fEval_wGw = fW.transpose()*fG*fW; 

  // check
  SelfAdjointEigenSolver<MatrixXd> es(fg,false);
  int  ne = es.eigenvalues().size();
  assert( es.eigenvalues()(0) >= -0.0000001 );

  // H, G2_inv
  H = - K2_inv*K1.transpose();
  for( int i = 0; i < fN_in; ++i ){
    for( int j = 0; j < fN; ++j ){
      fH(2*i,ix(j)) = H(i,j%fN);
      fH(2*i+1,iy(j)) = H(i,j%fN);
    }
  }

  fG2_inv = MatrixXd::Zero(2*fN_in, 2*fN_in);
  for( int i = 0; i < fN_in; ++i ){
    for( int j = 0; j < fN_in; ++j ){
      fG2_inv(2*i,2*j) = K2_inv(i,j);
      fG2_inv(2*i+1,2*j+1) = K2_inv(i,j);
    }
  }
}

bool TOperator_I::Check_Intersection( VectorXd& u ){
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


void TOperator_I::Trans_UU_center( VectorXd& uu )
{
  VectorXd centerU(2);
  VectorXd centerW(2);

  centerU = VectorXd::Zero(2); 
  centerW = VectorXd::Zero(2); 
  
  for( int i = 0; i < fN; ++i ){  
    centerW += fW.segment(2*i,2);
    centerU += uu.segment(2*i,2);
  }
  centerW /= (double)fN;
  centerU /= (double)fN;

  for( int i = 0; i < fNN; ++i )  
    uu.segment(2*i,2) += (-centerU + centerW);
}

void TOperator_I::Trans_UU_wrt_G( VectorXd& uu )
{
  VectorXd uu1(2*fNN);
  double tr, rt, cos_th, sin_th;

 
  tr = uu.head(2*fN).transpose()*fG*fW;
  rt = uu.head(2*fN).transpose()*fG*fW_d;
  cos_th = tr/sqrt(tr*tr + rt*rt);
  sin_th = rt/sqrt(tr*tr + rt*rt);

  for( int k = 0; k < fNN; ++k ){
    uu1(2*k) = cos_th*uu(2*k) + sin_th*uu(2*k+1); 
    uu1(2*k+1) = -sin_th*uu(2*k) + cos_th*uu(2*k+1); 
  }

  uu = uu1;

  this->Trans_UU_center( uu );
}


void TOperator_I::Update_UU_top( double eval, VectorXd& uu, int st1 )
{
 
}

