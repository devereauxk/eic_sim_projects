// ATHENA single track smearing
const int NETA_MAX = 500;
TH1F* Res_Handler = NULL;
TGraph *gmom_res[NETA_MAX];
TGraph *gdca_rphi_res[NETA_MAX];
TGraph *gdca_z_res[NETA_MAX];
// ATHENA primary vertex smearing
TH1F* VertexRes_X = NULL;
TH1F* VertexRes_Y = NULL;
TH1F* VertexRes_Z = NULL;

void setup_ATHENA_single_track_smearing(TFile* fin)
{
  cout << "Current max eta bin is " << NETA_MAX << ", adjust it inside fastsim.h when neccesary" << endl;

  if (!fin)
  {
    cout << "Cannot find the ATHENA setup file, abort..." << endl;
    exit(0);
  }

  Res_Handler = (TH1F*)fin->Get("Res_Handler");
  if (!Res_Handler)
  {
    cout << "No Res_Handler found, abort..." << endl;
    exit(0);
  }

  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    gmom_res[ibin] = (TGraph*)fin->Get(Form("gmom_res_%i",ibin));
    gdca_rphi_res[ibin] = (TGraph*)fin->Get(Form("gdca_rphi_res_%i",ibin));
    gdca_z_res[ibin] = (TGraph*)fin->Get(Form("gdca_z_res_%i",ibin));
  }
}

void setup_ATHENA_PV_smearing(TFile* fin)
{
  if (!fin)
  {
    cout << "Cannot find the ATHENA setup file, abort..." << endl;
    exit(0);
  }

  VertexRes_X = (TH1F*)fin->Get("VertexRes_X");
  VertexRes_Y = (TH1F*)fin->Get("VertexRes_Y");
  VertexRes_Z = (TH1F*)fin->Get("VertexRes_Z");
}

TVector3 smearPVTATHENA(TVector3 const& vtx, const float multi_35 = 3.6)
{ // smear primary vertex in transverse plane
  // multi_35 is the charged track multiplicity in |eta| < 3.5
  if (!VertexRes_X || !VertexRes_Y)
  {
    cout << "No PV smearing input found, abort..." << endl;
    exit(0);
  }

  TVector3 sVtx(-9999,-9999,vtx.Z());
  if (multi_35<2)
  {
    cout << "Event with too low multiplicty (<2), no primary vertex reconstructed" << endl;
    return sVtx;
  } 

  // NB: parameterization only goes up to multi ~ 20. Since it's relatively flat at higher multiplicty, for now we just use the multi==20 resolution for multi>20
  if (multi_35>20) multi_35 = 20;

  float rand_x = VertexRes_X->GetBinContent( VertexRes_X->FindBin(multi_35) ); // in unit of um
  float rand_y = VertexRes_Y->GetBinContent( VertexRes_Y->FindBin(multi_35) );
  rand_x /= 1E3; // convert from um to mm
  rand_y /= 1E3;

  sVtx.SetX(vtx.X() + rand_x);
  sVtx.SetY(vtx.Y() + rand_y);

  return sVtx;
}

TVector3 smearPVZATHENA(TVector3 const& vtx, const float multi_10 = 1.8)
{ // smear primary vertex in longitudinal direction
  // multi_10 is the charged track multiplicity in |eta| < 1.0
  if (!VertexRes_Z)
  {
    cout << "No PV smearing input found, abort..." << endl;
    exit(0);
  }

  TVector3 sVtx(vtx.X(),vtx.Y(),-9999);
  if (multi_10<2)
  {
    cout << "Event with too low multiplicty (<2), no primary vertex reconstructed" << endl;
    return sVtx;
  }

  // NB: parameterization only goes up to multi ~ 20. Since it's relatively flat at higher multiplicty, for now we just use the multi==20 resolution for multi>20
  if (multi_10>20) multi_10 = 20;

  float rand_z = VertexRes_Z->GetBinContent( VertexRes_Z->FindBin(multi_10) ); // in unit of um
  rand_z /= 1E3; // convert um to mm

  sVtx.SetZ(vtx.Z() + rand_z);

  return sVtx;
}

bool passTrackingATHENA(double p, double eta)
{
  if (fabs(eta)>3.5) return false;
  if (p<0.5) return false;
  if (p>100) return false; // NB: just because the current parameterization goes up to 100
  return true;
}

TLorentzVector smearMomATHENA(TLorentzVector const& mom4)
{ // takes true 4-momentum mom4 as input, return smeared 4-momentum sMom4
  TLorentzVector sMom4;
  sMom4.SetXYZM(-9999,-9999,-9999,0);

  int smear_graph_bin = Res_Handler->FindBin(mom4->PseudoRapidity());
  if (smear_graph_bin<1) return sMom4; // if not valid bin
  if (!gmom_res[smear_graph_bin-1]) return sMom4; // if no smearing parameter found

  float const pt = mom4.Perp();
  float const p = mom4.P();
  float const eta = mom4.PseudoRapidity();

  if (!passTrackingATHENA(p,eta)) return sMom4; // if track not in acceptance

  float rel_p_reso = gmom_res[smear_graph_bin-1]->Eval(p); // in unit of 1
  float p_reso = p*rel_p_reso;

  sP = gRandom->Gaus(p,p_reso);
  sPt = sP*TMath::Sin(mom4.Theta());

  sMom->SetXYZM(sPt * cos(mom4.Phi()), sPt * sin(mom4.Phi()), sPt * sinh(eta), mom4.M());
  return sMom; 
}

TVector3 smearPosATHENA(TVector3 const& mom, TVector3 const& pos)
{ // takes true mom and vertex position vector as input, return smeared vertex position sPos
  TVector3 sPos(-9999,-9999,-9999);

  int smear_graph_bin = Res_Handler->FindBin(mom4->PseudoRapidity());
  if (smear_graph_bin<1) return sPos; // if not valid bin
  if (!gdca_rphi_res[smear_graph_bin-1]) return sPos; // if no smearing parameter found
  if (!gdca_z_res[smear_graph_bin-1]) return sPos; // if no smearing parameter found

  float const pt = mom.Perp();
  float const p = mom.Mag();
  float const eta =  TMath::ATanH(mom.Pz()/mom.Mag());//Doing it this way to suppress TVector3 warnings; depends on any pre-cuts

  if (!passTrackingATHENA(p,eta)) return sPos;

  // resolutions are in microns, need in mm
  float dca_rphi_reso = gdca_rphi_res[smear_graph_bin-1]->Eval(p); // in unit of um
  float dca_z_reso = gdca_z_res[smear_graph_bin-1]->Eval(p); // in unit of um
  dca_rphi_reso /= 1E3; // convert to mm
  dca_z_reso /= 1E3;

  float rand_xy = gRandom->Gaus(0,dca_rphi_reso);
  float rand_z = gRandom->Gaus(0,dca_z_reso);

  // calculate new vertex position in transverse plane
  sPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  sPos -= momPerp.Unit() * rand_xy;

  // add back the z dimension
  sPos.SetZ(pos.Z() + rand_z);
  return sPos;
}

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

  TVector3 sPos(-9999,-9999,-9999);
  if (ietabin<0) return sPos; // if outside eta range

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

  sPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  sPos -= momPerp.Unit() * sigmaPosXY;

  return TVector3(sPos.X(), sPos.Y(), pos.Z() + sigmaPosZ);
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

  TVector3 sPos(-9999,-9999,-9999);
  if (ietabin<0) return sPos; // if outside eta range

  // sigma_posxy = A/pT + B (unit um)
  const double Axy[pos_reso_etabin] = {26, 31, 35, 41, 48, 59, 66, 69};
  const double Bxy[pos_reso_etabin] = {3.9, 4.0, 5.1, 4.9, 7.7, 11.2, 25.3, 72.1};

  // sigma_posz = A/pT + B (unit um)
  // const double Az[pos_reso_etabin] = {27, 35, 56, 112, 212, 373, 732, 1057};
  // const double Bz[pos_reso_etabin] = {3.3, 3.8, 5.4, 7.1, 16.0, 37.9, 87.7, 221};

  sigmaPosXY=gRandom->Gaus(0, sqrt(Axy[ietabin]*Axy[ietabin]/pt/pt+Bxy[ietabin]*Bxy[ietabin])*1E-3);
  sigmaPosZ = sigmaPosXY/sin(mom.Theta()); // consistent with the parameters for posz

  sPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  sPos -= momPerp.Unit() * sigmaPosXY;

  sPos.SetXYZ(sPos.X(), sPos.Y(), pos.Z() + sigmaPosZ);

  return sPos;
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

  TVector3 sPos(-9999,-9999,-9999);
  if (ietabin<0) return sPos; // if outside eta range

  // sigma_posxy = A/pT + B (unit um)
  const double Axy[pos_reso_etabin] = {14.1, 23.3, 49.3};
  const double Bxy[pos_reso_etabin] = {2.11, 3.32, 9.64};

  // // sigma_posz = A/pT + B (unit um)
  // const double Az[pos_reso_etabin] = {23.2, 78.3, 596.9};
  // const double Bz[pos_reso_etabin] = {2.64, 3.11, 41.05};

  sigmaPosXY=gRandom->Gaus(0, sqrt(Axy[ietabin]*Axy[ietabin]/pt/pt+Bxy[ietabin]*Bxy[ietabin])*1E-3);
  sigmaPosZ = sigmaPosXY/sin(mom.Theta()); // consistent with the parameters for posz

  sPos.SetXYZ(pos.X(),pos.Y(),0);
  TVector3 momPerp(-mom.Y(), mom.X(), 0.0);
  sPos -= momPerp.Unit() * sigmaPosXY;

  sPos.SetXYZ(sPos.X(), sPos.Y(), pos.Z() + sigmaPosZ);

  return sPos;
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















