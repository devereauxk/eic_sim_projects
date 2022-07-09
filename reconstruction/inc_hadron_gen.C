R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

const double degree = 180./TMath::Pi();

const int verbosity = 1; // use this to control debugging messages

using namespace std;

void inc_hadron_gen(const char* inFile = "merged.root", const char* outFile = "inc_merged.root", int nevt = 0, int data_type = 0, int gen_type = 0)
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

  // Event Class note EventBeagle used instead of EventPythia
  erhic::EventBeagle *event(NULL); //Note that I use Pointer

  // Access event Branch
  tree->SetBranchAddress("event",&event);

  // Initialize histograms
  TH3D* h3d_event_Ninc_vs_thickness_vs_b = new TH3D("h3d_event_Ninc_vs_thickness_vs_b","Ninc vs thickness vs b",20,0,4,20,0,4,20,0,4);
  h3d_event_Ninc_vs_thickness_vs_b->Sumw2();
  TH3D* h3d_event_Nincch_vs_thickness_vs_b = new TH3D("h3d_event_Nincch_vs_thickness_vs_b","Nincch vs thickness vs b",20,0,4,20,0,4,20,0,4);
  h3d_event_Nincch_vs_thickness_vs_b->Sumw2();

  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    //event level variables
    double Q2 = event->GetTrueQ2();
    double x = event->GetTrueX();
    double nu = event->GetTrueNu();
    double Ninc = event->NINC;
    double Nincch = event->NINCch;
    double thickness = event->Thickness;
    double b = event->b;

    //kinematic cuts
    bool flag_event_select = true;
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

    // fill histograms
    h3d_event_Ninc_vs_thickness_vs_b->Fill(Ninc, thickness, b);
    h3d_event_Nincch_vs_thickness_vs_b->Fill(Nincch, thickness, b);

  } // end of event by event processing

  //write histograms to outfile
  TFile* fout = new TFile(outFile,"recreate");

  h3d_event_Ninc_vs_thickness_vs_b->Write();
  h3d_event_Nincch_vs_thickness_vs_b->Write();

  fout->Write();
}
