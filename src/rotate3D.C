#include "TVector3.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "TRandom3.h"

//do
//.L rotate3D.C++
//rotate3D()
//

TRandom3 grand(1234);

double getPhi(const TVector3 vec, const TVector3 xx, const TVector3 yy)
{
  /*

.L ../rotate3D.C 
TVector3 xx(30, 0, 0), yy(0,90, 0)
TVector3 vec;
vec.SetXYZ(1,0,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(1,1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(0,1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(-1,1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(-1,0,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(-1,-1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(0,-1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(1,-1,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
vec.SetXYZ(1,-1E-10,6); getPhi(vec,xx,yy)*TMath::RadToDeg()
   */
  
  const TVector3 ux = xx.Unit();
  const TVector3 uy = yy.Unit();

  const double vecX = vec.Dot(ux);
  const double vecY = vec.Dot(uy);

  const double EPS = 1E-40;
  double phi = -999;
  if(TMath::Abs(vecX)<EPS){
    phi = (vecY>0?1:-1) * TMath::PiOver2();
  }
  else{
    phi = TMath::ATan(vecY/vecX);
    if(vecX<0){
      phi += TMath::Pi();
    }
  }

  while(phi>TMath::Pi()) phi-= TMath::TwoPi();

  return phi;
}

TVector3 DirectRotate(const TVector3 v0, const double tmpTheta, const double tmpPhi)
{
  TVector3 vec(v0);
  vec.RotateY(tmpTheta);
  vec.RotateZ(tmpPhi);
  return vec;
}

void getRMatrix(TMatrixD &aa, const double tmpTheta, const double tmpPhi, const bool kprint = false)
{
  //e0: initial frame with initial neutrino direction in [2]
  TVector3 e0[3];
  for(int ii=0; ii<3; ii++){
    e0[ii].SetXYZ(ii==0,ii==1,ii==2);

    if(kprint) e0[ii].Print();
  }

  //e1: final frame with final neutrino direction in [2]
  TVector3 e1[3];
  for(int ii=0; ii<3; ii++){
    e1[ii] = DirectRotate(e0[ii], tmpTheta, tmpPhi);
    
    if(kprint) e1[ii].Print();
  }

  //use e0 and e1 define roration matrix so that e1 = aa*e0, https://en.wikipedia.org/wiki/Three-dimensional_rotation_operator
  for(int ii=0; ii<3; ii++){
    for(int jj=0; jj<3; jj++){
      aa[ii][jj]=e0[ii].Dot(e1[jj]);
    }
  }
  if(kprint) aa.Print();

  //This is a VERY dumb method!! Try to find the analytical formula for aa[ii][jj] directly as a function of tmpPhi and tmpTheta
}

TVector3 MatrixRotate(const TMatrixD aa, const TVector3 v0, const bool kprint = kFALSE)
{
  //for any vector V0, final V1 = aa* V0

  TMatrixD b0(3,1);
  for(int ii=0; ii<3; ii++){
    b0[ii][0]=v0[ii];
  }
  if(kprint) b0.Print();

  const TMatrixD b1=aa*b0;
  TVector3 v1;
  v1.SetXYZ(b1[0][0],b1[1][0],b1[2][0]);

  return v1;
}

void rotate3D()
{

  TMatrixD aa(3,3);

  //test1: check with initial neutrino direction b0=e0[2], b1=aa*b0, b1 should be identical to e1[2]
  getRMatrix(aa, 30*TMath::DegToRad(), 45*TMath::DegToRad(), 1);
  const TVector3 v0(0,0,1);
  const TVector3 v1 = MatrixRotate(aa, v0, 1);
  v1.Print();

  //test2: check particle p and neutrino dot product before and after rotation
  TVector3 nu0;
  nu0.SetMagThetaPhi(4.26, 0, 0);
  const TVector3 p0(1,2,3);
  
  printf("\n\nbefore %f:\n", nu0.Dot(p0));
  nu0.Print();
  p0.Print();

  const double testTheta = 15*TMath::DegToRad();
  const double testPhi = 162*TMath::DegToRad();
  getRMatrix(aa, testTheta, testPhi);  
  
  const TVector3 nuTest1=MatrixRotate(aa, nu0);
  const TVector3 pTest1 =MatrixRotate(aa, p0);
  printf("\n\nafter %f:\n", nuTest1.Dot(pTest1));
  nuTest1.Print();
  pTest1.Print();

  //test3: show that aa (MatrixRotate) is the same as RotateY+RotateZ (DirectRotate)
  const TVector3 pTest2 = DirectRotate(p0, testTheta, testPhi);
  pTest2.Print();

  //test4: flat phi distribution before, show after rotation
  TH1D * h0=new TH1D("h0","before",30, -180, 180);
  TH1D * h1=new TH1D("h1","after", 30, -180, 180);

  TVector3 postNu, postX(1,0,0), postY(0,1,0);
  postNu.SetMagThetaPhi(40, testTheta, testPhi);
  postX=DirectRotate(postX, testTheta, testPhi);
  postY=DirectRotate(postY, testTheta, testPhi);
  
  for(int ii=0; ii<10000; ii++){
    TVector3 veci;
    veci.SetMagThetaPhi(5, 45*TMath::DegToRad(), TMath::TwoPi()*grand.Rndm());
    h0->Fill(veci.Phi()*TMath::RadToDeg());

    veci = DirectRotate(veci, testTheta, testPhi);
    const double phi1 = getPhi(veci, postX, postY)*TMath::RadToDeg();
    h1->Fill(phi1);
  }

  gStyle->SetHistMinimumZero();
  gStyle->SetOptStat(0);//"enoumr");
  h0->SetLineStyle(kDashed);
  h0->SetLineWidth(3);
  h0->Draw();
  gPad->SetGrid();
  h1->SetLineColor(kRed);
  h1->Draw("same");
}
