R__LOAD_LIBRARY(libeicsmear);
#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int sys_bins = 4;
const char* fin_dirs[sys_bins] = {"./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eAu_10_100_tauforOff_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eC_10_100_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eC_10_100_tauforOff_qhat0_nlo/outForPythiaMode/"};
const char* sys_name[sys_bins] = {"e+Au INC on", "e+Au INC off", "e+C INC on", "e+C INC off"};
const char* sys_abbr[sys_bins] = {"eAuINCon", "eAuINCff", "eCINCon", "eCINCoff"};

TH1D* D0_p[sys_bins] = {0};
TH1D* D0_pt[sys_bins] = {0};

void D0_mom_INC_diff_gen(const char* fout_name = "D0_mom_INC_diff.root")
{

  for (int isys = 0; isys < sys_bins; isys++)
  {
    TFile* f = new TFile(Form("%smerged.root", fin_dirs[isys]), "READ");

    TTree *tree = (TTree*)f->Get("EICTree");
    erhic::EventBeagle *event(NULL);
    tree->SetBranchAddress("event",&event);

    // initialize histograms
    D0_p[isys] = new TH1D(Form("h1d_D0_p_%s", sys_abbr[isys]), Form("h1d_D0_p_%s", sys_abbr[isys]), 100, 0, 50);
    D0_pt[isys] = new TH1D(Form("h1d_D0_pt_%s", sys_abbr[isys]), Form("h1d_D0_pt_%s", sys_abbr[isys]), 100, 0, 10);

    // fill histograms
    for(int ievt = 0; ievt < tree->GetEntries(); ievt++)
    {
      tree->GetEntry(ievt);
      for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = event->GetTrack(ipart);
        if (abs(part->Id()) == 421)
        {
          D0_p[isys]->Fill(part->GetP());
          D0_pt[isys]->Fill(part->GetPt());
        }
      }
    }
    cout<<sys_name[isys]<<" hists filled"<<endl;
  }

  // save histograms
  TFile* fout = new TFile(fout_name, "recreate");
  for(int isys = 0; isys < sys_bins; isys++)
  {
    D0_p[isys]->Write();
    D0_pt[isys]->Write();
  }
  fout->Write();

}
