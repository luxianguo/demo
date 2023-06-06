#include "TVector3.h"
#include "TMatrixD.h"

//do
//.L rotate3D.C++
//rotate3D()
//

void getRMatrix(TMatrixD &aa, const double gphi, const double gtheta, const bool kprint = false)
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
    e1[ii]=e0[ii];

    e1[ii].RotateY(gtheta);
    e1[ii].RotateZ(gphi);
    
    if(kprint) e1[ii].Print();
  }

  //use e0 and e1 define roration matrix so that e1 = aa*e0, https://en.wikipedia.org/wiki/Three-dimensional_rotation_operator
  for(int ii=0; ii<3; ii++){
    for(int jj=0; jj<3; jj++){
      aa[ii][jj]=e0[ii].Dot(e1[jj]);
    }
  }
  if(kprint) aa.Print();

  //This is a VERY dumb method!! Try to find the analytical formula for aa[ii][jj] directly as a function of gphi and gtheta
}

TVector3 getRotated(const TMatrixD aa, const TVector3 v0, const bool kprint = kFALSE)
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
  getRMatrix(aa, 45*TMath::DegToRad(), 30*TMath::DegToRad(), 1);
  const TVector3 v0(0,0,1);
  const TVector3 v1 = getRotated(aa, v0, 1);
  v1.Print();

  //test2: check particle p and neutrino dot product before and after rotation
  TVector3 nu0;
  nu0.SetMagThetaPhi(4.26, 0, 0);
  const TVector3 p0(1,2,3);
  
  printf("\n\nbefore %f:\n", nu0.Dot(p0));
  nu0.Print();
  p0.Print();

  getRMatrix(aa, 162*TMath::DegToRad(), 15*TMath::DegToRad());  
  
  const TVector3 nuTest=getRotated(aa, nu0);
  const TVector3 pTest =getRotated(aa, p0);
  printf("\n\nafter %f:\n", nuTest.Dot(pTest));
  nuTest.Print();
  pTest.Print();
}
