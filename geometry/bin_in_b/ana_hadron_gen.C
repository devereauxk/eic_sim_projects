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
    int inubin;

    TH2D* h2d_hadron_pt_vs_eta_gen_in_x[Q2bin][xbin];
    TH2D* h2d_hadron_z_vs_eta_gen_in_x[Q2bin][xbin];

    TH2D* h2d_hadron_pt_vs_eta_gen_in_nu[Q2bin][nubin];
    TH2D* h2d_hadron_z_vs_eta_gen_in_nu[Q2bin][nubin];

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
      inubin = -9999;

      h1d_parent_id = new TH1D(Form("h1d_hadron_%d_parent_id",hadron_id),"parent ID",10000,-0.5,9999.5);
      h1d_parent_id->Sumw2();

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_in_x_%d_%d",hadron_id,iQ2,ix),"D0 pt vs eta",100,0,10,40,-10,10);
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2][ix]->Sumw2();

          h2d_hadron_z_vs_eta_gen_in_x[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_in_x_%d_%d",hadron_id,iQ2,ix),"D0 z vs eta",100,0,1,40,-10,10);
          h2d_hadron_z_vs_eta_gen_in_x[iQ2][ix]->Sumw2();
        }

        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2][inu] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_in_nu_%d_%d",hadron_id,iQ2,inu),"D0 pt vs eta",100,0,10,40,-10,10);
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2][inu]->Sumw2();

          h2d_hadron_z_vs_eta_gen_in_nu[iQ2][inu] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_in_nu_%d_%d",hadron_id,iQ2,inu),"D0 z vs eta",100,0,1,40,-10,10);
          h2d_hadron_z_vs_eta_gen_in_nu[iQ2][inu]->Sumw2();
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
          delete h2d_hadron_pt_vs_eta_gen_in_x[iQ2][ix];
          delete h2d_hadron_z_vs_eta_gen_in_x[iQ2][ix];
        }

        for (int inu = 0; inu < nubin; ++inu)
        {
          delete h2d_hadron_pt_vs_eta_gen_in_nu[iQ2][inu];
          delete h2d_hadron_z_vs_eta_gen_in_nu[iQ2][inu];
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
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2][ix]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen_in_x[iQ2][ix]->Reset("ICESM");
        }

        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2][inu]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen_in_nu[iQ2][inu]->Reset("ICESM");
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

    void SetNuTrue(double _nu_true)
    {
      nu_true = _nu_true;

      inubin = -9999;
      for (int inu = 0; inu < nubin-1; ++inu)
      { // NB: do not loop the last bin which is inclusive
        if (nu_true>=nu_lo[inu] && nu_true<nu_hi[inu]) inubin = inu;
      }
    }

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
        if (part->GetStatus()!=1 && part->GetStatus()!=11 && part->GetStatus()!=2) continue; // 1 -- stable particles, 11 -- decay particles in Pythia, 2 -- decay particles in BeAGLE

        h1d_parent_id->Fill(abs(part->GetParentId()));

        TLorentzVector hadron_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(hadron_mom4_gen)/(nu_true*hadron_beam.M());

        // if (hadron_mom4_gen.Pt()<0.1) continue;

        if (verbosity>1) std::cout << "Particle id " << hadron_id << " with pt " << hadron_mom4_gen.Pt() << " z " << frag_z << " eta " << hadron_mom4_gen.PseudoRapidity() << std::endl;

        if (verbosity>1) std::cout << "Q2, x " << Q2_true << ", " << x_true << " bin " << iQ2bin << ", " << ixbin << std::endl;

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2bin][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x[iQ2bin][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2bin][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x[iQ2bin][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_x[Q2bin-1][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x[Q2bin-1][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_x[Q2bin-1][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x[Q2bin-1][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
        }

        if (iQ2bin>=0 && inubin>=0)
        {
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2bin][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_nu[iQ2bin][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2bin][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_nu[iQ2bin][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_nu[Q2bin-1][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_nu[Q2bin-1][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_nu[Q2bin-1][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_nu[Q2bin-1][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
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
          h2d_hadron_pt_vs_eta_gen_in_x[iQ2][ix]->Write();
          h2d_hadron_z_vs_eta_gen_in_x[iQ2][ix]->Write();
        }

        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_nu[iQ2][inu]->Write();
          h2d_hadron_z_vs_eta_gen_in_nu[iQ2][inu]->Write();
        }
      }
    }
};

bool event_w_charm(erhic::EventPythia* event, int gen_type)
{
  if (gen_type==0)
  { // Pythia 6
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=21 && !flag_search_group ) continue;
      if ( part->GetStatus()==21 )
      { // entering into the group of particles with KS=21
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=21 && flag_search_group ) break;
    }
  }
  else
  { // BeAGLE
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = event->GetNTracks()-1; ipart >=0; ipart--)
    { // faster looping backwards
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=3 && !flag_search_group ) continue;
      if ( part->GetStatus()==3 )
      { // entering into the group of particles with KS=3
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=3 && flag_search_group ) break;
    }
  }

  return false;
}

void ana_hadron_gen(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hist.root", int nevt = 0, int data_type = 0, int gen_type = 0, int b_lo = 0, int b_hi = 15)
{
  cout << "Data Type: ";
  if (data_type==0) cout << "EIC" << endl;
  else if (data_type==1) cout << "HERMES" << endl;
  else cout << "CLAS" << endl;

  cout << "Generator Type: ";
  if (gen_type==0) cout << "Pythia6" << endl;
  else cout << "BeAGLE" << endl;

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  // Event Class
  erhic::EventBeagle *event(NULL);
  if (gen_type==0) erhic::EventPythia *event(NULL); //Note that I use Pointer

  // Access event Branch
  tree->SetBranchAddress("event",&event);

  // Inclusive event counter before any event selections
  TH2D* h2d_logx_logQ2 = new TH2D("h2d_logx_logQ2","inclusive event counter",100,-3,0,100,0,3);
  h2d_logx_logQ2->Sumw2();

  // Event counter after event selections
  TH1D* h1d_nevt_in_x[Q2bin][xbin] = {0};
  TH1D* h1d_nevt_w_charm_in_x[Q2bin][xbin] = {0};
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_in_x[iQ2][ix] = new TH1D(Form("h1d_nevt_in_x_%d_%d",iQ2,ix),"event counter",1,0.5,1.5);
      h1d_nevt_in_x[iQ2][ix]->Sumw2();

      h1d_nevt_w_charm_in_x[iQ2][ix] = new TH1D(Form("h1d_nevt_w_charm_in_x_%d_%d",iQ2,ix),"event counter",1,0.5,1.5);
      h1d_nevt_w_charm_in_x[iQ2][ix]->Sumw2();
    }
  }

  TH1D* h1d_nevt_in_nu[Q2bin][nubin] = {0};
  TH1D* h1d_nevt_w_charm_in_nu[Q2bin][nubin] = {0};
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_nu[iQ2][inu] = new TH1D(Form("h1d_nevt_in_nu_%d_%d",iQ2,inu),"event counter",1,0.5,1.5);
      h1d_nevt_in_nu[iQ2][inu]->Sumw2();

      h1d_nevt_w_charm_in_nu[iQ2][inu] = new TH1D(Form("h1d_nevt_w_charm_in_nu_%d_%d",iQ2,inu),"event counter",1,0.5,1.5);
      h1d_nevt_w_charm_in_nu[iQ2][inu]->Sumw2();
    }
  }

  // Hadron analyzers
  hadron_gen ana_D0(421); // which will also analyze -421
  hadron_gen ana_Lc(4122);
  hadron_gen ana_charged_pion(211);
  hadron_gen ana_neutral_pion(111);
  hadron_gen ana_charged_kaon(321);
  hadron_gen ana_proton(2212);

  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    //Write Out Q2
    double Q2 = event->GetTrueQ2(); // Can also do event->QSquared
    double x = event->GetTrueX();
    double nu = event->GetTrueNu();

    h2d_logx_logQ2->Fill(log10(x),log10(Q2));
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);

    bool flag_event_select = true;
    // thickness cut
    // if (gen_type != 0 && (event->Thickness < thickness_lo || event->Thickness > thickness_hi)) flag_event_select = false;
    // impact parameter (b) cut
    if (gen_type != 0 && (event->b < b_lo || event->b > b_hi)) flag_event_select = false;
    // kinematics cuts
    if (data_type==0)
    { // EIC event selection (from Yuxiang)
      if (event->GetTrueY()<0.05 || event->GetTrueY()>0.8) flag_event_select = false;
    }
    else if (data_type==1)
    { // HERMES event selections (from HERMES paper)
      if (event->GetTrueY()>0.85) flag_event_select = false;
      if (event->GetTrueW2()<4) flag_event_select = false;
      if (event->GetTrueNu()<6) flag_event_select = false;
    }
    else
    { // CLAS event selection (to be implemented)

    }

    if (!flag_event_select) continue;

    int iQ2bin = -9999, ixbin = -9999, inubin = -9999;
    for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
    { // NB: do not loop the last bin which is inclusive
      if (Q2>=Q2_lo[iQ2] && Q2<Q2_hi[iQ2]) iQ2bin = iQ2;
    }
    for (int ix = 0; ix < xbin-1; ++ix)
    { // NB: do not loop the last bin which is inclusive
      if (x>=x_lo[ix] && x<x_hi[ix]) ixbin = ix;
    }
    for (int inu = 0; inu < nubin-1; ++inu)
    { // NB: do not loop the last bin which is inclusive
      if (nu>=nu_lo[inu] && nu<nu_hi[inu]) inubin = inu;
    }

    if (iQ2bin>=0 && ixbin>=0)
    {
      h1d_nevt_in_x[iQ2bin][ixbin]->Fill(1);
      h1d_nevt_in_x[iQ2bin][xbin-1]->Fill(1);
      h1d_nevt_in_x[Q2bin-1][ixbin]->Fill(1);
      h1d_nevt_in_x[Q2bin-1][xbin-1]->Fill(1);

      if ( event_w_charm(event,gen_type) )
      {
        h1d_nevt_w_charm_in_x[iQ2bin][ixbin]->Fill(1);
        h1d_nevt_w_charm_in_x[iQ2bin][xbin-1]->Fill(1);
        h1d_nevt_w_charm_in_x[Q2bin-1][ixbin]->Fill(1);
        h1d_nevt_w_charm_in_x[Q2bin-1][xbin-1]->Fill(1);
      }
    }

    if (iQ2bin>=0 && inubin>=0)
    {
      h1d_nevt_in_nu[iQ2bin][inubin]->Fill(1);
      h1d_nevt_in_nu[iQ2bin][nubin-1]->Fill(1);
      h1d_nevt_in_nu[Q2bin-1][inubin]->Fill(1);
      h1d_nevt_in_nu[Q2bin-1][nubin-1]->Fill(1);

      if ( event_w_charm(event,gen_type) )
      {
        h1d_nevt_w_charm_in_nu[iQ2bin][inubin]->Fill(1);
        h1d_nevt_w_charm_in_nu[iQ2bin][nubin-1]->Fill(1);
        h1d_nevt_w_charm_in_nu[Q2bin-1][inubin]->Fill(1);
        h1d_nevt_w_charm_in_nu[Q2bin-1][nubin-1]->Fill(1);
      }
    }

    ana_D0.SetQ2True(event->GetTrueQ2());
    ana_D0.SetXTrue(event->GetTrueX());
    ana_D0.SetNuTrue(event->GetTrueNu());
    ana_D0.FillGenKin(event);

    ana_Lc.SetQ2True(event->GetTrueQ2());
    ana_Lc.SetXTrue(event->GetTrueX());
    ana_Lc.SetNuTrue(event->GetTrueNu());
    ana_Lc.FillGenKin(event);

    ana_charged_pion.SetQ2True(event->GetTrueQ2());
    ana_charged_pion.SetXTrue(event->GetTrueX());
    ana_charged_pion.SetNuTrue(event->GetTrueNu());
    ana_charged_pion.FillGenKin(event);

    ana_neutral_pion.SetQ2True(event->GetTrueQ2());
    ana_neutral_pion.SetXTrue(event->GetTrueX());
    ana_neutral_pion.SetNuTrue(event->GetTrueNu());
    ana_neutral_pion.FillGenKin(event);

    ana_charged_kaon.SetQ2True(event->GetTrueQ2());
    ana_charged_kaon.SetXTrue(event->GetTrueX());
    ana_charged_kaon.SetNuTrue(event->GetTrueNu());
    ana_charged_kaon.FillGenKin(event);

    ana_proton.SetQ2True(event->GetTrueQ2());
    ana_proton.SetXTrue(event->GetTrueX());
    ana_proton.SetNuTrue(event->GetTrueNu());
    ana_proton.FillGenKin(event);
  }

  TFile* fout = new TFile(outFile,"recreate");
  h2d_logx_logQ2->Write();
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_in_x[iQ2][ix]->Write();
      h1d_nevt_w_charm_in_x[iQ2][ix]->Write();
    }

    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_nu[iQ2][inu]->Write();
      h1d_nevt_w_charm_in_nu[iQ2][inu]->Write();
    }
  }

  ana_D0.Write();
  ana_Lc.Write();
  ana_charged_pion.Write();
  ana_neutral_pion.Write();
  ana_charged_kaon.Write();
  ana_proton.Write();

  fout->Write();
}
