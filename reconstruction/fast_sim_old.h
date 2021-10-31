TLorentzVector smearMomDM(TLorentzVector const& mom4, const int Bfield_type = 0)
{ // smear magnitude only
  // B field type 0 -- Barbar, 1 -- Beast
  assert(Bfield_type<2); // safe guard

  float const pt = mom4.Perp();
  float const p = mom4.P();
  float sPt = -9999;
  float sP = -9999;
  float eta = mom4.PseudoRapidity();

  const int p_reso_etabin = 3;
  const double p_reso_etalo[p_reso_etabin] = {0, 1, 2.5};
  const double p_reso_etahi[p_reso_etabin] = {1, 2.5, 3.5};

  int ietabin = -9999;
  for (int ieta = 0; ieta < p_reso_etabin; ++ieta)
  {
    if (fabs(eta)>=p_reso_etalo[ieta] && fabs(eta)<p_reso_etahi[ieta]) ietabin = ieta;
  }

  TLorentzVector sMom4;
  sMom4.SetXYZM(-9999,-9999,-9999,0);
  if (ietabin<0) return sMom4; // if outside eta range

  // low pt limit (DM values coming from Hybrid simulation)
  float PT_LO_Barbar = -9999, PT_LO_Beast = -9999;
  if (fabs(eta)<1) {PT_LO_Barbar = 0.2; PT_LO_Beast = 0.4;}
  if (fabs(eta)>=1 && fabs(eta)<1.5) {PT_LO_Barbar = 0.15; PT_LO_Beast = 0.3;}
  if (fabs(eta)>=1.5 && fabs(eta)<2.0) {PT_LO_Barbar = 0.07; PT_LO_Beast = 0.16;}
  if (fabs(eta)>=2.0 && fabs(eta)<2.5) {PT_LO_Barbar = 0.13; PT_LO_Beast = 0.22;}
  if (fabs(eta)>=2.5 && fabs(eta)<3.0) {PT_LO_Barbar = 0.1; PT_LO_Beast = 0.15;}

  if (Bfield_type==0 && pt<PT_LO_Barbar) return sMom4;
  if (Bfield_type==1 && pt<PT_LO_Beast) return sMom4;

  // sigma_p/p = A*p + B (unit %)
  const double A_Barbar[p_reso_etabin] = {0.04, 0.04, 0.2};
  const double B_Barbar[p_reso_etabin] = {1, 2, 5};

  const double A_Beast[p_reso_etabin] = {0.02, 0.02, 0.1};
  const double B_Beast[p_reso_etabin] = {0.5, 1, 2};

  if (Bfield_type==0) sP = gRandom->Gaus(p, p*sqrt(A_Barbar[ietabin]*A_Barbar[ietabin]*p*p+B_Barbar[ietabin]*B_Barbar[ietabin])*1E-2);
  if (Bfield_type==1) sP = gRandom->Gaus(p, p*sqrt(A_Beast[ietabin]*A_Beast[ietabin]*p*p+B_Beast[ietabin]*B_Beast[ietabin])*1E-2);

  sPt = sP*TMath::Sin(mom4.Theta());
  sMom4.SetXYZM(sPt * cos(mom4.Phi()), sPt * sin(mom4.Phi()), sPt * sinh(eta), mom4.M());
  return sMom4; 
}

TVector3 smearPosDM(TVector3 const& mom, TVector3 const& pos)
{
  float const pt = mom.Perp();
  float sigmaPosXY = 0, sigmaPosZ;
  double eta =  TMath::ATanH(mom.Pz()/mom.Mag());//Doing it this way to suppress TVector3 Warnings; depends on any pre-cuts    

  const int pos_reso_etabin = 3;
  const double pos_reso_etalo[pos_reso_etabin] = {0, 1, 2};
  const double pos_reso_etahi[pos_reso_etabin] = {1, 2, 3};

  int ietabin = -9999;
  for (int ieta = 0; ieta < pos_reso_etabin; ++ieta)
  {
    if (fabs(eta)>=pos_reso_etalo[ieta] && fabs(eta)<pos_reso_etahi[ieta]) ietabin = ieta;
  }

  TVector3 newPos(-9999,-9999,-9999);
  if (ietabin<0) return newPos; // if outside eta range

  // sigma_posxy = A/pT + B (unit um)
  const double Axy[pos_reso_etabin] = {30, 40, 60};
  const double Bxy[pos_reso_etabin] = {5, 10, 15};

  // sigma_posz = A/pT + B (unit um) 
  // const double Az[pos_reso_etabin] = {30, 100, 350};
  // const double Bz[pos_reso_etabin] = {5, 20, 40};

  sigmaPosXY=gRandom->Gaus(0, sqrt(Axy[ietabin]*Axy[ietabin]/pt/pt+Bxy[ietabin]*Bxy[ietabin])*1E-3);

  // FIX ME: double check the relationship between the pointing resolution along XY and Z
  if (fabs(eta)<=1) sigmaPosZ = gRandom->Gaus(0, sqrt(30.*30./pt/pt+5.*5.)*1E-3);
  else sigmaPosZ = sigmaPosXY/sin(mom.Theta());

  newPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  newPos -= momPerp.Unit() * sigmaPosXY;

  return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TLorentzVector smearMomLBL(TLorentzVector const& mom4, const int Bfield_type = 0)
{ // smear magnitude only
  // B field type 0 -- Barbar, 1 -- Beast
  assert(Bfield_type<2); // safe guard

  float const pt = mom4.Perp();
  float const p = mom4.P();
  float sPt = -9999;
  float sP = -9999;
  float eta = mom4.PseudoRapidity();

  // absolute value binning
  const int p_reso_etabin = 8;
  const double p_reso_etalo[p_reso_etabin] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5};
  const double p_reso_etahi[p_reso_etabin] = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};

  int ietabin = -9999;
  for (int ieta = 0; ieta < p_reso_etabin; ++ieta)
  {
    if (fabs(eta)>=p_reso_etalo[ieta] && fabs(eta)<p_reso_etahi[ieta]) ietabin = ieta;
  }

  TLorentzVector sMom4;
  sMom4.SetXYZM(-9999,-9999,-9999,0);
  if (ietabin<0) return sMom4; // if outside eta range

  // low pt limit (FIX ME: just copy pasted values from hybrid, cross-check with Rey about what would be the R threshold)
  float PT_LO_Barbar = -9999, PT_LO_Beast = -9999;
  if (fabs(eta)<1) {PT_LO_Barbar = 0.09; PT_LO_Beast = 0.195;} // this is to hit the outer layer
  if (fabs(eta)>=1 && fabs(eta)<1.5) {PT_LO_Barbar = 0.15; PT_LO_Beast = 0.3;}
  if (fabs(eta)>=1.5 && fabs(eta)<2.0) {PT_LO_Barbar = 0.07; PT_LO_Beast = 0.16;}
  if (fabs(eta)>=2.0 && fabs(eta)<2.5) {PT_LO_Barbar = 0.13; PT_LO_Beast = 0.22;}
  if (fabs(eta)>=2.5 && fabs(eta)<3.0) {PT_LO_Barbar = 0.1; PT_LO_Beast = 0.15;}

  if (Bfield_type==0 && pt<PT_LO_Barbar) return sMom4;
  if (Bfield_type==1 && pt<PT_LO_Beast) return sMom4;

  // sigma_p/p = A*p + B (unit %)
  const double A_Barbar[p_reso_etabin] = {0.041, 0.034, 0.034, 0.026, 0.041, 0.085, 0.215, 0.642};
  const double B_Barbar[p_reso_etabin] = {0.773, 0.906, 0.922, 1.000, 1.551, 2.853, 5.254, 9.657};

  const double A_Beast[p_reso_etabin] = {0.018, 0.016, 0.016, 0.012, 0.018, 0.039, 0.103, 0.281};
  const double B_Beast[p_reso_etabin] = {0.382, 0.431, 0.424, 0.462, 0.721, 1.331, 2.441, 4.716};

  if (Bfield_type==0) sP = gRandom->Gaus(p, p*sqrt(A_Barbar[ietabin]*A_Barbar[ietabin]*p*p+B_Barbar[ietabin]*B_Barbar[ietabin])*1E-2);
  if (Bfield_type==1) sP = gRandom->Gaus(p, p*sqrt(A_Beast[ietabin]*A_Beast[ietabin]*p*p+B_Beast[ietabin]*B_Beast[ietabin])*1E-2);

  sPt = sP*TMath::Sin(mom4.Theta());
  sMom4.SetXYZM(sPt * cos(mom4.Phi()), sPt * sin(mom4.Phi()), sPt * sinh(eta), mom4.M());
  return sMom4; 
}

TVector3 smearPosLBL(TVector3 const& mom, TVector3 const& pos)
{
  float const pt = mom.Perp();
  float sigmaPosXY = 0, sigmaPosZ;
  double eta =  TMath::ATanH(mom.Pz()/mom.Mag()); // Doing it this way to suppress TVector3 Warnings; depends on any pre-cuts   

  // absolute value binning
  const int pos_reso_etabin = 8; 
  const double pos_reso_etalo[pos_reso_etabin] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5};
  const double pos_reso_etahi[pos_reso_etabin] = {0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};

  int ietabin = -9999;
  for (int ieta = 0; ieta < pos_reso_etabin; ++ieta)
  {
    if (fabs(eta)>=pos_reso_etalo[ieta] && fabs(eta)<pos_reso_etahi[ieta]) ietabin = ieta;
  }

  TVector3 newPos(-9999,-9999,-9999);
  if (ietabin<0) return newPos; // if outside eta range

  // sigma_posxy = A/pT + B (unit um)
  const double Axy[pos_reso_etabin] = {26, 31, 35, 41, 48, 59, 66, 69};
  const double Bxy[pos_reso_etabin] = {3.9, 4.0, 5.1, 4.9, 7.7, 11.2, 25.3, 72.1};

  // sigma_posz = A/pT + B (unit um)
  // const double Az[pos_reso_etabin] = {27, 35, 56, 112, 212, 373, 732, 1057};
  // const double Bz[pos_reso_etabin] = {3.3, 3.8, 5.4, 7.1, 16.0, 37.9, 87.7, 221};

  sigmaPosXY=gRandom->Gaus(0, sqrt(Axy[ietabin]*Axy[ietabin]/pt/pt+Bxy[ietabin]*Bxy[ietabin])*1E-3);
  sigmaPosZ = sigmaPosXY/sin(mom.Theta()); // consistent with the parameters for posz

  newPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  newPos -= momPerp.Unit() * sigmaPosXY;

  newPos.SetXYZ(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);

  return newPos;
}

TLorentzVector smearMomHybrid(TLorentzVector const& mom4, const int Bfield_type = 0)
{ // smear magnitude only
  // B field type 0 -- Barbar, 1 -- Beast
  assert(Bfield_type<2); // safe guard

  float const pt = mom4.Perp();
  float const p = mom4.P();
  float sPt = -9999;
  float sP = -9999;
  float eta = mom4.PseudoRapidity();

  TLorentzVector sMom4;
  sMom4.SetXYZM(-9999,-9999,-9999,0);
  float A_Barbar = -9999, B_Barbar = -9999;
  float A_Beast = -9999, B_Beast = -9999;

  // low pt limit
  float PT_LO_Barbar = -9999, PT_LO_Beast = -9999;
  if (fabs(eta)<1) {PT_LO_Barbar = 0.2; PT_LO_Beast = 0.4;}
  if (fabs(eta)>=1 && fabs(eta)<1.5) {PT_LO_Barbar = 0.15; PT_LO_Beast = 0.3;}
  if (fabs(eta)>=1.5 && fabs(eta)<2.0) {PT_LO_Barbar = 0.07; PT_LO_Beast = 0.16;}
  if (fabs(eta)>=2.0 && fabs(eta)<2.5) {PT_LO_Barbar = 0.13; PT_LO_Beast = 0.22;}
  if (fabs(eta)>=2.5 && fabs(eta)<3.0) {PT_LO_Barbar = 0.1; PT_LO_Beast = 0.15;}

  if (Bfield_type==0 && pt<PT_LO_Barbar) return sMom4;
  if (Bfield_type==1 && pt<PT_LO_Beast) return sMom4;

  if (fabs(eta)<1)
  {
    if (p<5)
    {
      A_Barbar = 0.21; B_Barbar = 0.34;
      A_Beast = 0.11; B_Beast = 0.18;
    }
    else if (p<30)
    {
      A_Barbar = 0.06; B_Barbar = 1.09;
      A_Beast = 0.03; B_Beast = 0.54;
    }
    else return sMom4;
  }
  else if (fabs(eta)<2.5)
  {
    if (p<8)
    {
      A_Barbar = 0.24; B_Barbar = 0.67;
      A_Beast = 0.11; B_Beast = 0.33;
    }
    else if (p<30)
    {
      A_Barbar = 0.07; B_Barbar = 1.81;
      A_Beast = 0.04; B_Beast = 0.88;
    }
    else return sMom4;
  }
  else if (fabs(eta)<3.5)
  {
    A_Barbar = 0.09; B_Barbar = 3.71;
    A_Beast = 0.05; B_Beast = 1.90;
  }
  else return sMom4;

  if (Bfield_type==0) sP = gRandom->Gaus(p, p*sqrt(A_Barbar*A_Barbar*p*p+B_Barbar*B_Barbar)*1E-2);
  if (Bfield_type==1) sP = gRandom->Gaus(p, p*sqrt(A_Beast*A_Beast*p*p+B_Beast*B_Beast)*1E-2);

  sPt = sP*TMath::Sin(mom4.Theta());
  sMom4.SetXYZM(sPt * cos(mom4.Phi()), sPt * sin(mom4.Phi()), sPt * sinh(eta), mom4.M());
  return sMom4; 
}

TVector3 smearPosHybrid(TVector3 const& mom, TVector3 const& pos)
{
  float const pt = mom.Perp();
  float sigmaPosXY = 0, sigmaPosZ;
  double eta =  TMath::ATanH(mom.Pz()/mom.Mag()); // Doing it this way to suppress TVector3 Warnings; depends on any pre-cuts   

  // absolute value binning
  const int pos_reso_etabin = 3; 
  const double pos_reso_etalo[pos_reso_etabin] = {0, 1, 2.5};
  const double pos_reso_etahi[pos_reso_etabin] = {1, 2.5, 3.5};

  int ietabin = -9999;
  for (int ieta = 0; ieta < pos_reso_etabin; ++ieta)
  {
    if (fabs(eta)>=pos_reso_etalo[ieta] && fabs(eta)<pos_reso_etahi[ieta]) ietabin = ieta;
  }

  TVector3 newPos(-9999,-9999,-9999);
  if (ietabin<0) return newPos; // if outside eta range

  // sigma_posxy = A/pT + B (unit um)
  const double Axy[pos_reso_etabin] = {14.1, 23.3, 49.3};
  const double Bxy[pos_reso_etabin] = {2.11, 3.32, 9.64};

  // // sigma_posz = A/pT + B (unit um)
  // const double Az[pos_reso_etabin] = {23.2, 78.3, 596.9};
  // const double Bz[pos_reso_etabin] = {2.64, 3.11, 41.05};

  sigmaPosXY=gRandom->Gaus(0, sqrt(Axy[ietabin]*Axy[ietabin]/pt/pt+Bxy[ietabin]*Bxy[ietabin])*1E-3);
  sigmaPosZ = sigmaPosXY/sin(mom.Theta()); // consistent with the parameters for posz

  newPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  newPos -= momPerp.Unit() * sigmaPosXY;

  newPos.SetXYZ(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);

  return newPos;
}

void passing_DIRC(TLorentzVector const& mom4, bitset<4>& binary_id, const float Bfield_type = 0, const float R_in = 0.9)
{ // FIX ME: only consider firing or not, need to consider sigma separation later 
  if (Bfield_type>1) return;

  // const float R_in = 0.5; // (move closest possible)
  float Bfield = 1.4;
  if (Bfield_type==1) Bfield = 3.0;
  float thr_kin = 0.3*Bfield*R_in/2; 

  float pt = mom4.Pt();
  float eta = mom4.PseudoRapidity();
  if (pt<=thr_kin || fabs(eta)>1) return; // no new info added from DIRC since particle does not hit it

  float p = mom4.P();
  if (p>0.00048) binary_id[0] = 1; // e fired
  if (p>0.13) binary_id[1] = 1; // pi fired
  if (p>0.47) binary_id[2] = 1; // K fired
  if (p>0.88) binary_id[3] = 1; // p fired
  if (p>6) for (int i = 0; i < 4; ++i) binary_id[i] = 0; // FIX ME: brute-force hard cut
}

void passing_hside_dRICH(TLorentzVector const& mom4, bitset<4>& binary_id)
{
  float eta = mom4.PseudoRapidity();
  if (eta<=1 || eta>4) return; // no new info added from dRICH since particle does not hit it

  float p = mom4.P();
  if (p>0.00245) binary_id[0] = 1; // e fired
  if (p>0.69) binary_id[1] = 1; // pi fired
  if (p>2.46) binary_id[2] = 1; // K fired
  if (p>4.67) binary_id[3] = 1; // p fired
  if (p>60) for (int i = 0; i < 4; ++i) binary_id[i] = 0; // FIX ME: brute-force hard cut
}

void passing_eside_dRICH(TLorentzVector const& mom4, bitset<4>& binary_id)
{
  float eta = mom4.PseudoRapidity();
  if (eta>=-1 || eta<-4) return; // no new info added from dRICH since particle does not hit it

  float p = mom4.P();
  if (p>0.00245) binary_id[0] = 1; // e fired
  if (p>0.69) binary_id[1] = 1; // pi fired
  if (p>2.46) binary_id[2] = 1; // K fired
  if (p>4.67) binary_id[3] = 1; // p fired
  if (p>60) for (int i = 0; i < 4; ++i) binary_id[i] = 0; // FIX ME: brute-force hard cut
}

void passing_DM_PID(TLorentzVector const& mom4, bitset<4>& binary_id)
{ // FIX ME: not considering kinematic reach for now
  // NB: ALL-Si stop at 43cm, Hybrid TPC 20-78cm, maybe r @ 50cm can be a reasonble PID reference
  bool flag_PID = false;

  float eta = mom4.PseudoRapidity();
  float p = mom4.P();
  if (fabs(eta)<1 && p<6) flag_PID = true;
  if (eta<-1 && p<10) flag_PID = true;
  if (eta>1 && p<50) flag_PID = true;

  if (flag_PID) for (int i = 0; i < 4; ++i) binary_id[i] = 1;

  return;
}

void identify_charged_hadrons(const int truth_id, bitset<4> const& binary_id, float& _prob_e, float& _prob_pi, float& _prob_K, float& _prob_p)
{
  if (binary_id.to_ulong()>=5)
  {
    if (fabs(truth_id)==11)
    {
      _prob_e = 1.0;
      _prob_pi = 0;
      _prob_K = 0;
      _prob_p = 0;
    } 
    else if (fabs(truth_id)==211 || fabs(truth_id)==13)
    {
      _prob_e = 0;
      _prob_pi = 1.0;
      _prob_K = 0;
      _prob_p = 0;
    } 
    else if (fabs(truth_id)==321)
    {
      _prob_e = 0;
      _prob_pi = 0;
      _prob_K = 1.0;
      _prob_p = 0;
    } 
    else
    {
      _prob_e = 0;
      _prob_pi = 0;
      _prob_K = 0;
      _prob_p = 1.0;
    }
  }
  else if (binary_id.to_ulong()>=3)
  {
    if (fabs(truth_id)==11)
    {
      _prob_e = 1.0;
      _prob_pi = 0;
      _prob_K = 0;
      _prob_p = 0;
    } 
    else if (fabs(truth_id)==211 || fabs(truth_id)==13)
    {
      _prob_e = 0;
      _prob_pi = 1.0;
      _prob_K = 0;
      _prob_p = 0;
    } 
    else
    {
      _prob_e = 0;
      _prob_pi = 0;
      _prob_K = 0.6;
      _prob_p = 0.4;
    } 
  }
  else if (binary_id.to_ulong()>=1)
  {
    if (fabs(truth_id)==11)
    {
      _prob_e = 1.0;
      _prob_pi = 0;
      _prob_K = 0;
      _prob_p = 0;
    } 
    else
    {
      _prob_e = 0;
      _prob_pi = 0.7;
      _prob_K = 0.2;
      _prob_p = 0.1;
    } 
  }
}















