R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);

#include "fastjet/ClusterSequence.hh"
#include <iostream>
using namespace fastjet;
using namespace std;

const int verbosity = 0;

TH1D* h1d_jet_eec = NULL;
TH1D* h1d_jet_pt = NULL;

double calculate_distance(PseudoJet p0, PseudoJet p1)
{
  float dphiabs = fabs(p0.phi() - p1.phi());

  if (dphiabs > TMath::Pi()) dphiabs = 2*TMath::Pi() - dphiabs;

  float dy = p0.eta() - p1.eta();
  return sqrt(dy*dy + dphiabs*dphiabs);
}

class Correlator_Builder
{
  private:
    vector<PseudoJet> particle_list;
    int mult;
    vector<vector<double>> pair_list;
    double scale;

  public:
    Correlator_Builder(vector<PseudoJet> _particle_list, float _scale)
    {
      particle_list = _particle_list;
      mult = particle_list.size();
      scale = _scale;
    }

    void make_pairs()
    {
      // produces mult X mult upper triangular matrix with elements being the
      // distance between partices i and j
      for (int i = 0; i < mult; i++)
      {
        vector<double> inner_list(mult, 0.0);
        for (int j = i; j < mult; j++)
        {
          inner_list[j] = calculate_distance(particle_list[i], particle_list[j]);
        }
        pair_list.push_back(inner_list);
      }
    }

    void construct_EEC(TH1D* h1d_jet_eec)
    {
      int overlap = 0;
      for (int i = 0; i < mult; i++)
      {
        for (int j = i; j < mult; j++)
        {
          double dist12 = pair_list[i][j];
          double eec_weight = (particle_list[i].pt() * particle_list[j].pt()) / pow(scale, 2);

          if (i == j) overlap++;
          if (overlap == 0) eec_weight = eec_weight*2;
          if (overlap > 0) eec_weight = eec_weight*1;

          h1d_jet_eec->Fill(dist12, eec_weight);
        }
      }
    }
};

void eec_hists(const char* inFile = "merged.root", const char* outFile = "hists_eec.root")
{
  //Event Class
  erhic::EventPythia *event(NULL);

  TFile *f = new TFile(inFile);

  // compute log bins for eec histogram
  double xmin = 1E-4;
  double xmax = 1;
  int nbins = 50;
  Double_t lbins[nbins+1];
  double binwidth = (log10(xmax) - log10(xmin)) / nbins;
  for (int i = 0; i < nbins+1; i++)
  {
    lbins[i] = TMath::Power(10, log10(xmin) + binwidth * i);
  }

  // histogram definitions
  h1d_jet_pt = new TH1D("h1d_jet_pt","jet pt",800,0,800);
  h1d_jet_pt->Sumw2();

  h1d_jet_eec = new TH1D("h1d_jet_eec","jet eec",50,lbins);
  h1d_jet_eec->Sumw2();

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nevt = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nevt<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  //Loop Over Events
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    tree->GetEntry(ievt);

    if (ievt%1000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    /*
    //Write Out Q2
    double Q2 = event->GetQ2(); //Can also do event->QSquared
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);

    // process cuts
    bool flag_direct = false;
    if (event->GetProcess()==99) flag_direct = true;
    if (event->GetProcess()==131 || event->GetProcess()==132) flag_direct = true;
    if (event->GetProcess()==135 || event->GetProcess()==136) flag_direct = true;
    if (!flag_direct) continue; // only process direct processes
    */


    // particle enumeration, addition to jet reco setup, and total pt calculation
    erhic::ParticleMC* particle;
    vector<PseudoJet> jet_constits;
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      particle = event->GetTrack(ipart);

      // use all fsp particles w/ < 3.5 eta, not including scattered electron, for jet reconstruction
      if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11)
      {
        PseudoJet constit = new PseudoJet(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
        constit.set_user_index(ipart);
        jet_constits.push_back(constit);
      }

    }

    // jet reconstruction
    JetDefinition R1jetdef(antikt_algorithm, 1.0);
    ClusterSequence cs(jet_constits, R1jetdef);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    // jet processing
    for (unsigned ijet = 0; ijet < jets.size(); ijet++)
    {
      // cuts on jet kinematics, require jet_pt >= 5GeV, |jet_eta| <= 2.5
      if (jets[ijet].pt() < 5 || fabs(jets[ijet].eta()) > 2.5) continue;
      h1d_jet_pt->Fill(jets[ijet].pt());

      // cuts on jet constituent kinematics, require consitituents_pt >= 0.5GeV, |consitituents_eta| <= 3.5
      // take only charged constituents for eec calculation
      vector<PseudoJet> constituents = jets[ijet].constituents();
      vector<PseudoJet> charged_constituents;
      for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
      {
        if (constituents[iconstit].pt() < 0.5 || fabs(constituents[iconstit].eta()) > 3.5) continue;

        //int ip = constituents[iconstit].user_index();
        //float charge = event->GetTrack(ip)->Id().Info()->Charge();
        //cout<<"constituent pt:"<<constituents[iconstit].pt()<<" track pt:"<<event->GetTrack(ip)->GetPt()<<" charge:"<<charge<<endl;
        //if (charge != 0) charged_constituents.push_back(constituents[iconstit]);
      }

      if (charged_constituents.size() < 1) continue;

      // eec calculation
      Correlator_Builder cb(charged_constituents, jets[ijet].pt());
      cb.make_pairs();
      cb.construct_EEC(h1d_jet_eec);

    }

  }

  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  h1d_jet_pt->Write();
  h1d_jet_eec->Write();

}
