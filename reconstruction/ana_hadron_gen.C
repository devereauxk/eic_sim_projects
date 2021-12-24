R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

const double degree = 180./TMath::Pi();

const int verbosity = 1; // use this to control debugging messages

using namespace std;

class hadron_gen
{
  private:
    // evt info
    int hadron_id;

    double x_true;
    double Q2_true;
    double nu_true;

    TLorentzVector hadron_beam;
    TLorentzVector struck_quark;

    TH1D* h1d_parent_id;

    int iQ2bin;
    int ixbin;
    TH2D* h2d_hadron_pt_vs_eta_gen[Q2bin][xbin];
    TH2D* h2d_hadron_z_vs_eta_gen[Q2bin][xbin];

  public:
    hadron_gen(int _hadron_id)
    {
      cout << "Constructing analyzing module to study hadronization of particle with ID " << _hadron_id << endl;
      hadron_id = abs(_hadron_id);
      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      iQ2bin = -9999;
      ixbin = -9999;

      h1d_parent_id = new TH1D(Form("h1d_hadron_%d_parent_id",hadron_id),"parent ID",10000,-0.5,9999.5);
      h1d_parent_id->Sumw2();

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_%d_%d",hadron_id,iQ2,ix),"D0 pt vs eta",100,0,10,160,-4,4);
          h2d_hadron_pt_vs_eta_gen[iQ2][ix]->Sumw2();

          h2d_hadron_z_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_%d_%d",hadron_id,iQ2,ix),"D0 z vs eta",100,0,1,160,-4,4);
          h2d_hadron_z_vs_eta_gen[iQ2][ix]->Sumw2();
        }
      }
    }
    virtual ~hadron_gen()
    {
      delete h1d_parent_id;
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          delete h2d_hadron_pt_vs_eta_gen[iQ2][ix];
          delete h2d_hadron_z_vs_eta_gen[iQ2][ix];
        }
      }
    };

    void Reset()
    {
      x_true = 1E1;
      Q2_true = 1E-5;
      nu_true = -9999;

      iQ2bin = -9999;
      ixbin = -9999;

      h1d_parent_id->Reset("ICESM");
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen[iQ2][ix]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen[iQ2][ix]->Reset("ICESM");
        }
      }
    }

    void SetQ2True(double _Q2_true)
    { 
      Q2_true = _Q2_true; 

      iQ2bin = -9999;
      for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
      { // NB: do not loop the last bin which is inclusive
        if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
      }      
    }
    void SetXTrue(double _x_true)
    { 
      x_true = _x_true; 

      ixbin = -9999;
      for (int ix = 0; ix < xbin-1; ++ix)
      { // NB: do not loop the last bin which is inclusive
        if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
      }
    }

    void SetNuTrue(double _nu_true) { nu_true = _nu_true; };

    void FillGenKin(erhic::EventPythia* py_evt)
    {
      erhic::ParticleMC* proton = py_evt->GetTrack(1);
      if (proton!=NULL)
      {
        assert(abs(proton->Id())!=2212);
        hadron_beam = proton->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (abs(part->Id())!=hadron_id) continue;

        h1d_parent_id->Fill(abs(part->GetParentId()));

        TLorentzVector hadron_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(hadron_mom4_gen)/(nu_true*hadron_beam.M());

        if (verbosity>1) std::cout << "Particle id " << hadron_id << " with pt " << hadron_mom4_gen.Pt() << " z " << frag_z << " eta " << hadron_mom4_gen.PseudoRapidity() << std::endl;

        if (verbosity>1) std::cout << "Q2, x " << Q2_true << ", " << x_true << " bin " << iQ2bin << ", " << ixbin << std::endl;

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_hadron_pt_vs_eta_gen[iQ2bin][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen[iQ2bin][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen[iQ2bin][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen[iQ2bin][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen[Q2bin-1][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen[Q2bin-1][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen[Q2bin-1][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen[Q2bin-1][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
        } 
      }
    }

    void Write()
    {
      h1d_parent_id->Write();
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen[iQ2][ix]->Write();
          h2d_hadron_z_vs_eta_gen[iQ2][ix]->Write();
        }
      }
    }
};

void ana_hadron_gen(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hist.root", int nevt = 0)
{ 
  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event);

  hadron_gen ana_D0(421); // which will also analyze -421
  hadron_gen ana_Lc(4122);
  hadron_gen ana_charged_pion(211);
  hadron_gen ana_neutral_pion(111);
  hadron_gen ana_charged_kaon(321);
  hadron_gen ana_proton(2212);

  TH1D* h1d_nevt[Q2bin][xbin] = {0};
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt[iQ2][ix] = new TH1D(Form("h1d_nevt_%d_%d",iQ2,ix),"event counter",1,0.5,1.5);
      h1d_nevt[iQ2][ix]->Sumw2(); 
    }
  }
  TH2D* h2d_logx_logQ2 = new TH2D("h2d_logx_logQ2","inclusive event counter",100,-3,0,100,0,3);
  h2d_logx_logQ2->Sumw2();
  
  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {    
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    //Write Out Q2
    double Q2 = event->GetQ2(); // Can also do event->QSquared
    double x = event->GetX();

    h2d_logx_logQ2->Fill(log10(x),log10(Q2));
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);
    // if (Q2>10) continue; // process low Q2

    if (event->GetTrueY()<0.05 || event->GetTrueY()>0.8) continue;

    int iQ2bin = -9999, ixbin = -9999;
    for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
    { // NB: do not loop the last bin which is inclusive
      if (Q2>=Q2_lo[iQ2] && Q2<Q2_hi[iQ2]) iQ2bin = iQ2;
    } 
    for (int ix = 0; ix < xbin-1; ++ix)
    { // NB: do not loop the last bin which is inclusive
      if (x>=x_lo[ix] && x<x_hi[ix]) ixbin = ix;
    } 

    if (iQ2bin>=0 && ixbin>=0)
    {
      h1d_nevt[iQ2bin][ixbin]->Fill(1);
      h1d_nevt[iQ2bin][xbin-1]->Fill(1);
      h1d_nevt[Q2bin-1][ixbin]->Fill(1);
      h1d_nevt[Q2bin-1][xbin-1]->Fill(1);
    } 
    
    ana_D0.SetQ2True(event->GetQ2());
    ana_D0.SetXTrue(event->GetX());
    ana_D0.SetNuTrue(event->GetNu());
    ana_D0.FillGenKin(event);

    ana_Lc.SetQ2True(event->GetQ2());
    ana_Lc.SetXTrue(event->GetX());
    ana_Lc.SetNuTrue(event->GetNu());
    ana_Lc.FillGenKin(event);

    ana_charged_pion.SetQ2True(event->GetQ2());
    ana_charged_pion.SetXTrue(event->GetX());
    ana_charged_pion.SetNuTrue(event->GetNu());
    ana_charged_pion.FillGenKin(event);

    ana_neutral_pion.SetQ2True(event->GetQ2());
    ana_neutral_pion.SetXTrue(event->GetX());
    ana_neutral_pion.SetNuTrue(event->GetNu());
    ana_neutral_pion.FillGenKin(event);

    ana_charged_kaon.SetQ2True(event->GetQ2());
    ana_charged_kaon.SetXTrue(event->GetX());
    ana_charged_kaon.SetNuTrue(event->GetNu());
    ana_charged_kaon.FillGenKin(event);

    ana_proton.SetQ2True(event->GetQ2());
    ana_proton.SetXTrue(event->GetX());
    ana_proton.SetNuTrue(event->GetNu());
    ana_proton.FillGenKin(event);
  }

  TFile* fout = new TFile(outFile,"recreate");
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt[iQ2][ix]->Write();
    }
  }
  h2d_logx_logQ2->Write();
  ana_D0.Write();
  ana_Lc.Write();
  ana_charged_pion.Write();
  ana_neutral_pion.Write();
  ana_charged_kaon.Write();
  ana_proton.Write();

  fout->Write();
}
