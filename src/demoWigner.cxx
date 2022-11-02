#define EPS 1E-5

double GetTwoBoostAdlerPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull, TLorentzVector iniNfull)
{
  const bool kprint = true;

  /*
first to boost everything to the initial nucleon rest frame, record the nu-mu, or equivalently the delta, direction (called v0);
then boost everything (except v0) to the delta rest frame and use the initial-nucleon-rest-frame v0 as z-axis.
Then y should still be nu-mu in the delta-rest-frame, it should be equal to v0 cross mu.
  */
  //lab frame
  TLorentzVector delta = pifull + nucleonfull;

  const TLorentzVector nuold = nufull;
  const TLorentzVector muonold = muonfull;
  const TLorentzVector piold = pifull;
  const TLorentzVector nucleonold = nucleonfull;
  const TLorentzVector iniNold = iniNfull;
  const TLorentzVector delold = delta;

  if(kprint){
    printf("\n\n\n==============================================================================================================================\n");
    printf("******************************* AnaFunctions::GetTwoBoostAdlerPhi in the Lab Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
  }

  //check 4-momentum conserved at vertex
  const TLorentzVector p4balance = nufull + iniNfull - muonfull - pifull - nucleonfull;
  if(p4balance.P()>EPS){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi 4 momentum not conserved at vertes!\n\n\n");
    cout<<"iniNfull "<<endl; iniNfull.Print();
    cout<<"nufull "<<endl; nufull.Print();
    cout<<"muonfull "<<endl; muonfull.Print();
    cout<<"pifull "<<endl; pifull.Print();
    cout<<"nucleonfull "<<endl; nucleonfull.Print();
    p4balance.Print();
    exit(1);
  }

  //go to initial nucleon rest frame first
  const TVector3 boostToIniN = -iniNfull.BoostVector();

  nufull.Boost(boostToIniN);
  muonfull.Boost(boostToIniN);
  pifull.Boost(boostToIniN);
  nucleonfull.Boost(boostToIniN);
  delta.Boost(boostToIniN);
  iniNfull.Boost(boostToIniN);

  if(iniNfull.P()>EPS){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerAngle something wrong with iniN boosting!\n\n\n");
    iniNfull.Print();
    exit(1);
  }

  const TVector3 Zaxis = delta.Vect().Unit();//should be equal to nu-mu

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetTwoBoostAdlerPhi in the Initial Nucleon Rest Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
  }

  //from iniN rest frame to delta rest frame
  const TVector3 boostToDelta = -delta.BoostVector();

  nufull.Boost(boostToDelta);
  muonfull.Boost(boostToDelta);
  pifull.Boost(boostToDelta);
  nucleonfull.Boost(boostToDelta);
  delta.Boost(boostToDelta);
  iniNfull.Boost(boostToDelta);

  //boost to delta rest frame, check delta at rest after boost
  if(delta.P()>EPS){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerAngle something wrong with boosting!\n\n\n"); 
    delta.Print(); 
    exit(1);
  }

  const TVector3 Yaxis = (nufull.Vect().Cross(muonfull.Vect())).Unit();
  const double yzdot = Yaxis.Dot(Zaxis);
  if(fabs(yzdot)>EPS){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi Y and Z not perpendicular! %f\n\n\n", yzdot); 
    cout<<"iniNfull "<<endl; iniNfull.Print();
    cout<<"nuold "<<endl; nuold.Print();
    cout<<"muonold "<<endl; muonold.Print();
    cout<<"piold "<<endl; piold.Print();
    cout<<"nucleonold "<<endl; nucleonold.Print();
    cout<<"nufull "<<endl; nufull.Print();
    cout<<"muonfull "<<endl; muonfull.Print();
    cout<<"pifull "<<endl; pifull.Print();
    cout<<"nucleonfull "<<endl; nucleonfull.Print();
    cout<<"Yaxis "<<endl; Yaxis.Print();
    cout<<"Zaxis "<<endl; Zaxis.Print();
    exit(1);
  }

  const TVector3 Xaxis = Yaxis.Cross(Zaxis);

  const double nucleonX = nucleonfull.Vect().Dot(Xaxis);
  const double nucleonY = nucleonfull.Vect().Dot(Yaxis);
  const double nucleonR = TMath::Sqrt(nucleonX*nucleonX + nucleonY*nucleonY);

  //in 0-180
  double phi = TMath::ACos(nucleonX/nucleonR)*TMath::RadToDeg();
  if(phi<0 || phi>180){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi wrong domain in ACos %f %f %f\n\n\n", nucleonX, nucleonY, phi); exit(1);
  }
  if(nucleonY<0){
    phi = 360-phi;
  }

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetTwoBoostAdlerPhi in Delta Rest Frame Adler Phi %f\n", phi);
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
    cout<<"Xaxis "<<endl; Xaxis.Print(); cout<<endl;
    cout<<"Yaxis "<<endl; Yaxis.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
    cout<<"nu-mu"<<endl; (nufull-muonfull).Vect().Print(); cout<<endl;
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Check this final nucleon 3 component in The Two-Boost Frame\n");
    const double nucleonZ = nucleonfull.Vect().Dot(Zaxis);
    const TVector3 tmp(nucleonX, nucleonY, nucleonZ);
    tmp.Print();
    printf("==============================================================================================================================\n");
  }

  return phi;
}

double GetOneBoostAdlerPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull, TLorentzVector iniNfull)
{
  const bool kprint = true;

  //same as two-boost, just differ a Wigner Rotation which doesn't change the relative position of particles and therefore Adler angles
  //lab frame
  const TLorentzVector delta = pifull + nucleonfull;

  if(kprint){
    printf("\n\n\n==============================================================================================================================\n");
    printf("******************************* AnaFunctions::GetOneBoostAdlerPhi in the Lab Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
  }

  const TVector3 boost = -delta.BoostVector();
  
  nufull.Boost(boost);
  muonfull.Boost(boost);
  nucleonfull.Boost(boost);

  const TVector3 Zaxis = (nufull-muonfull).Vect().Unit();//only after boost! has to be q direction, otherwise won't be perpendicular to nu cross mu when there if Fermi motion
  const TVector3 Yaxis = (nufull.Vect().Cross(muonfull.Vect())).Unit();
  const TVector3 Xaxis = Yaxis.Cross(Zaxis);

  const double nucleonX = nucleonfull.Vect().Dot(Xaxis);
  const double nucleonY = nucleonfull.Vect().Dot(Yaxis);


  //in 0-180
  double phi = TMath::ATan2(nucleonY, nucleonX)*TMath::RadToDeg();
  if(phi<0){
    phi+=360;
  }

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetOneBoostAdlerPhi in Delta Rest Frame Adler Phi %f\n", phi);
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    TLorentzVector tmpdelta = delta;
    tmpdelta.Boost(boost);
    cout<<"delta "<<endl; tmpdelta.Print(); cout<<endl;
    cout<<"Xaxis "<<endl; Xaxis.Print(); cout<<endl;
    cout<<"Yaxis "<<endl; Yaxis.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
    cout<<"nu-mu"<<endl; (nufull-muonfull).Vect().Print(); cout<<endl;
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Check this final nucleon 3 component in The One-Boost Frame\n");
    const double nucleonZ = nucleonfull.Vect().Dot(Zaxis);
    const TVector3 tmp(nucleonX, nucleonY, nucleonZ);
    tmp.Print();
    printf("==============================================================================================================================\n");
  }

  return phi;
}

void demoWigner()
{
  const TLorentzVector iniNfull(-0.154022,0.076381,-0.046338,0.919298);
  const TLorentzVector nufull(0.000000,0.000000,3.533616,3.533616);
  const TLorentzVector muonfull(-0.304323,0.272725,3.243601,3.270948);
  const TLorentzVector pifull(0.045749,-0.011444,0.162628,0.216542);
  const TLorentzVector nucleonfull(0.104551,-0.184900,0.081049,0.965423);

  GetTwoBoostAdlerPhi(nufull, muonfull, pifull, nucleonfull, iniNfull);
  GetOneBoostAdlerPhi(nufull, muonfull, pifull, nucleonfull, iniNfull);

}
