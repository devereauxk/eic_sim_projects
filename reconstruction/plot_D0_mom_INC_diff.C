R__LOAD_LIBRARY(libeicsmear);
#include "bins.h"
using namespace std;

const int sys_bins = 2
const char* fin_dirs[sys_bins] = {"./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/"};

TH1D* D0_p[sys_bins] = {0};
TH1D* D0_pt[sys_bins] = {0};

static int cno = 0;
mcs(-1);

void standardLatex()
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.50,0.84,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8");
}

TFile* fins[sys_bins] = {0};
for (int isys = 0; isys < sys_bins; isys++)
{
  fins[isys] = new TFile(Form("%smerged.root", fin_dirs[isys]), "READ");

  TTree *tree = (TTree*)f->Get("EICTree");
  erhic::EventBeagle *event(NULL);
  tree->SetBranchAddress("event",&event);

  // initialize histograms
  D0_p[isys] = new TH1D(Form("D0_p_%d", isys), Form("D0_p_%d", isys), 100, 0, 60);
  D0_pt[isys] = new TH1D(Form("D0_pt_%d", isys), Form("D0_pt_%d", isys), 100, 0, 10);

  // fill histograms
  for(int ievt = 0; ievt < tree->GetEntries(); ievt++)
  {
    for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
    {
      erhic::ParticleMC* part = py_evt->GetTrack(ipart);
      if (abs(part->Id()) == 421)
      {
        D0_p[isys]->Fill(sqrt(pow(part->GetPx(),2) + pow(part->GetPy(),2) + pow(part->GetPz(),2)));
        D0_pt[isys]->Fill(part->GetPt());
      }
    }
  }

  // print histograms
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 60;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 24000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}}");
    htemp.GetYaxis()->SetTitle("normalized counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    for (int isys = 0; isys < sys_bins; isys++)
    {
      D0_p[isys]->SetMarkerColor(sys_color[isys]);
      D0_p[isys]->SetLineColor(sys_color[isys]);
      D0_p[isys]->Draw("hsame");
      leg.AddEntry(D0_p[isys],Form("%s", sys_name[isys]), "l");
    }

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%stemp.pdf\")", cno-1, outDir) );
  }


}
