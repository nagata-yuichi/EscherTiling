/*
Author: Yuichi Nagata
Copyright (c) 2022, Yuichi Nagata
This code is released under the MIT License.
*/

#ifndef __OPERATOR_IR__
#include "operator_IR.h"
#endif

TOperator_IR::TOperator_IR( int N, int N_in ) : TOperator_I( N, N_in )
{
  fR_cos = VectorXd::Zero(fNN); 
  fR_sin = VectorXd::Zero(fNN); 
  fRot = new MatrixXd [fNN];
  for( int i = 0; i < fNN; ++i )
    fRot[i] = MatrixXd(2,2);
  fRotD = new MatrixXd [fNN];
  for( int i = 0; i < fNN; ++i )
    fRotD[i] = MatrixXd(2,2);

  fGG = MatrixXd::Zero(2*fNN, 2*fNN);
  fAA = SparseMatrix<double>(2*fNN,2*fNN);
  fAA.reserve(8*fNN);
  fTTWW = VectorXd::Zero(2*fNN); 
  fTTWW_d = VectorXd::Zero(2*fNN);
} 

TOperator_IR::~TOperator_IR() 
{

} 

void TOperator_IR::SetParameter() 
{
  TOperator_I::SetParameter();
}

void TOperator_IR::SetInit( VectorXd w,  VectorXd w_in, MatrixXi kki ) 
{
  TOperator_I::SetInit( w, w_in, kki );
}


void TOperator_IR::Cal_u( int st1 ) 
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

  fTTWW = fGG*fWW;
  for( int iter = 0; iter < 3; ++iter ){
    eval = this->Cal_eval_gzai_Rotation( gzai, iter ); 
    assert( eval > 0.0 );
    u = fBdj*gzai.head(fMd); 
    u += fBsj*gzai.tail(fMs); 
    u_in = fG2_inv*fTTWW.bottomRows(2*fN_in) + fH*u;
    uu.head(2*fN) = u;        
    uu.tail(2*fN_in) = u_in;  

    this->Set_RotationMatrices( uu );
  }

  eval = (uu.transpose()*fGG*uu)(0,0) - 2.0*(uu.transpose()*fTTWW)(0,0) + fEval_wGGw;
  assert( eval > 0.0 );
    
  // this->Check_Eval( uu, eval ); // check

  if( eval < fEval_best ){
    if( Check_Intersection( u ) == false ){ 
      this->Trans_UU_center( uu );     
      this->Update_UU_top( eval, uu, st1 );
    }
  }
}


void TOperator_IR::Cal_u_Procrustes( int st1 ) 
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
  VectorXd tmp(2*fN);
  VectorXd tmp_d(2*fN);
  double tr, rt, cos_th, sin_th; 
  VectorXd ur(2*fN);

  this->Set_Bdj_Bsj( st1 ); // Set fBdj and fBsj

  fTTWW = fGG*fWW;
  fTTWW_d = fGG*fWW_d;
  for( int iter = 0; iter < 3; ++iter ){ 
    tmp = fTTWW.topRows(2*fN) + fH.transpose()*fTTWW.bottomRows(2*fN_in); 
    p.head(fMd) = fBdj.transpose()*tmp;
    p.tail(fMs) = fBsj.transpose()*tmp;

    tmp_d = fTTWW_d.topRows(2*fN) + fH.transpose()*fTTWW_d.bottomRows(2*fN_in); 
    p_d.head(fMd) = fBdj.transpose()*tmp_d;
    p_d.tail(fMs) = fBsj.transpose()*tmp_d;
    this->Cal_L_p( p, p_d, L, iter );

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
    eval = fEval_wGGw - fTTWW.bottomRows(2*fN_in).transpose()*fG2_inv*fTTWW.bottomRows(2*fN_in) - lambda;
    assert( eval > 0.0 );

    gzai.head(fMk) = A2(0,1)*x_unit.head(fMk) + (lambda-A2(0,0))*y_unit.head(fMk);
    gzai.head(fMk) = gzai.head(fMk).normalized()*sqrt(lambda);
    this->Cal_L_gzai( gzai, L.transpose() );
    u = fBdj*gzai.head(fMd);
    u += fBsj*gzai.tail(fMs); 

    tr = u.transpose()*tmp;
    rt = u.transpose()*tmp_d;
    cos_th = tr/sqrt(tr*tr + rt*rt);
    sin_th = rt/sqrt(tr*tr + rt*rt);

    for( int k = 0; k < fN; ++k ){
      ur(ix(k)) = cos_th*u(ix(k)) - sin_th*u(iy(k)); 
      ur(iy(k)) = sin_th*u(ix(k)) + cos_th*u(iy(k)); 
    }
    u = ur; 
    u_in = fG2_inv*fTTWW.bottomRows(2*fN_in) + fH*u;  

    uu.head(2*fN) = u;        
    uu.tail(2*fN_in) = u_in;  

    this->Set_RotationMatrices( uu );  
  }

  eval = (uu.transpose()*fGG*uu)(0,0) - 2.0*(uu.transpose()*fTTWW)(0,0) + fEval_wGGw;
  assert( eval > 0.0 );

  // this->Check_Eval( uu, eval ); // check

  if( eval < fEval_best ){
    if( Check_Intersection( u ) == false ){
      this->Trans_UU_center( uu );     
      this->Update_UU_top( eval, uu, st1 );
    }
  }
}

void TOperator_IR::Check_Eval( VectorXd uu,  double eval )
{
  VectorXd v(2); 
  double sum = 0.0;
  for( int i = 0; i < fNN; ++i ){
    int size = fNei_Size(i);

    for( int s = 0; s < size; ++s ){
      int j = fNei(i,s); 
      double aij = fA(i,j);
      assert( i != j ); 
      v(0) = uu(iix(j)) - uu(iix(i)) - fR_cos(i)*(fWW(iix(j)) - fWW(iix(i))) + fR_sin(i)*(fWW(iiy(j)) - fWW(iiy(i)));
      v(1) = uu(iiy(j)) - uu(iiy(i)) - fR_sin(i)*(fWW(iix(j)) - fWW(iix(i))) - fR_cos(i)*(fWW(iiy(j)) - fWW(iiy(i)));
      sum += aij*(v(0)*v(0)+v(1)*v(1));
    }
  }
	
  v(0) = uu(iix(0)) - fWW(iix(0));
  v(1) = uu(iiy(0)) - fWW(iiy(0));
  sum += v(0)*v(0)+v(1)*v(1);   
  
  printf( "%lf %lf -> %lf\n", sum , eval, eval-sum ); 
  assert( fabs(eval-sum) < 0.000001 ); 
}


void TOperator_IR::Set_RotationMatrices( VectorXd uu )
{
  int j,size;
  VectorXd ud(2*fNN); 
  VectorXd wd(2*fNN);
  VectorXd wd_c(2*fNN);
  double tr, rt, sqr;
  double di, aij;
  VectorXd RRWW(2*fNN);

  for( int i = 0; i < fNN; ++i ){
    size = fNei_Size(i);
    for( int s = 0; s < size; ++s ){
      j = fNei(i,s); 
      aij = sqrt(fA(i,j)); 
      ud(2*s) = aij*(uu(iix(j)) - uu(iix(i)));  
      ud(2*s+1) = aij*(uu(iiy(j)) - uu(iiy(i))); 
      wd(2*s) = aij*(fWW(iix(j)) - fWW(iix(i)));
      wd(2*s+1) = aij*(fWW(iiy(j)) - fWW(iiy(i)));
      wd_c(2*s) = wd(2*s+1); 
      wd_c(2*s+1) = -wd(2*s);
    }

    tr =ud.head(2*size).dot(wd.head(2*size));
    rt =ud.head(2*size).dot(wd_c.head(2*size));
    sqr = sqrt( tr*tr + rt*rt );
    fR_cos(i) = tr/sqr;  
    fR_sin(i) = -rt/sqr; 
    fRot[i](0,0) = fR_cos(i);  
    fRot[i](0,1) = -fR_sin(i);
    fRot[i](1,0) = fR_sin(i);
    fRot[i](1,1) = fR_cos(i);
  }

  /*
  //// naive calculation
  TT = MatrixXd::Zero(2*fNN,2*fNN);
  for( int i = 0; i < fNN; ++i ){
    size = fNei_Size(i);
    for( int s = 0; s < size; ++s ){
      j = fNei(i,s); 
      TT.block(2*i,2*i,2,2) += fA(i,j)*fRot[i];
      TT.block(2*j,2*j,2,2) += fA(i,j)*fRot[i];
      TT.block(2*i,2*j,2,2) -= fA(i,j)*fRot[i];
      TT.block(2*j,2*i,2,2) -= fA(i,j)*fRot[i];
    }
  }
  TT(0,0) += 1.0;
  TT(1,1) += 1.0;

  fTTWW = TT*fWW;
  */
  ////

  for( int i = 0; i < fNN; ++i ){
    size = fNei_Size(i);
    di = fDout(i,i);
    fRotD[i] = di*fRot[i]; 

    for( int s = 0; s < size; ++s ){
      j = fNei(i,s); 
      fRotD[i] += fA(j,i)*fRot[j]; 
    }
  }
  fRotD[0](0,0) += 1.0; 
  fRotD[0](1,1) += 1.0;

  // fTT = (D1*R -A^t*R + bar_R - R*A); 

  fTTWW = VectorXd::Zero(2*fNN);
  RRWW = VectorXd::Zero(2*fNN);
  for( int i = 0; i < fNN; ++i ){
    fTTWW.segment(2*i,2) += fRotD[i]*fWW.segment(2*i,2); 
    RRWW.segment(2*i,2) += fRot[i]*fWW.segment(2*i,2); 
  }
  fTTWW += -fAA.transpose()*RRWW; 

  VectorXd tmp(2*fNN);
  tmp = fAA*fWW;
  for( int i = 0; i < fNN; ++i )
    fTTWW.segment(2*i,2) -= fRot[i]*tmp.segment(2*i,2); 

  for( int i = 0; i < fNN; ++i ){
    fTTWW_d(iix(i)) = fTTWW(iiy(i));
    fTTWW_d(iiy(i)) = -fTTWW(iix(i));
  }
}


double TOperator_IR::Cal_eval_gzai_Rotation( VectorXd& gzai, int iter )  
{
  MatrixXd C(fMd+fMs,fMd+fMs);
  VectorXd b(fMd+fMs);
  static LLT<MatrixXd> llt;

  if( iter == 0 ){
    C.topLeftCorner(fMd,fMd) = fBdj.transpose()*fG*fBdj; // Inner
    C.topRightCorner(fMd,fMs) = fBdj.transpose()*fG*fBsj;
    C.bottomLeftCorner(fMs,fMd) = fBsj.transpose()*fG*fBdj;
    C.bottomRightCorner(fMs,fMs) = fBsj.transpose()*fG*fBsj;
    llt.compute(C);
  }

  VectorXd tmp(2*fN);
  tmp = fTTWW.topRows(2*fN) + fH.transpose()*fTTWW.bottomRows(2*fN_in);

  b.topRows(fMd) = fBdj.transpose()*tmp;
  b.bottomRows(fMs) = fBsj.transpose()*tmp;

  gzai = llt.solve(b) ;
  // assert(fabs((C*gzai - b).norm()) < 0.000001 );

  double eval;
  eval = -gzai.dot(b) + fEval_wGGw - fTTWW.bottomRows(2*fN_in).transpose()*fG2_inv*fTTWW.bottomRows(2*fN_in);

  return eval;
}


void TOperator_IR::SetK() 
{
  TOperator_I::SetK();

  fGG = MatrixXd::Zero(2*fNN,2*fNN);
  for( int i = 0; i < fNN; ++i ){
    for( int j = 0; j < fNN; ++j ){
      fGG(iix(i),iix(j)) = fKK(i,j);
      fGG(iiy(i),iiy(j)) = fKK(i,j);
    }
  }

  fEval_wGGw = fWW.transpose()*fGG*fWW;

  MatrixXd AA =  MatrixXd::Zero(2*fNN,2*fNN); 
  for( int i = 0; i < fNN; ++i ){
    for( int j = 0; j < fNN; ++j ){
      AA(iix(i),iix(j)) = fA(i,j);
      AA(iiy(i),iiy(j)) = fA(i,j);
    }
  }

  double val;
  fTripletList.clear();
  for( int j = 0; j < 2*fNN; ++j ){
    for( int i = 0; i < 2*fNN; ++i ){
      val = AA(i,j); 
      if( fabs(val) > 0.000000001 ) 
	fTripletList.push_back(T(i,j,val)); 
    }
  }
  fAA.setFromTriplets(fTripletList.begin(), fTripletList.end());
}


