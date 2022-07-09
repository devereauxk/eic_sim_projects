#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 11;
const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)", "e+Pb"};
const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia", "ePb"};

const int energy_bins = 6;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

static int cno = 0;

const int MaxVarbin = 10;


void plot_inc_hadron(const char* inFile = "inc_merged.root", const int sys_option = 0, const int energy_option = 0, const char* outDir = "figs/")
{
  mcs(-1);

  TFile* fin = new TFile(inFile,"READ");

  // load histograms
  TH3D* h3d_event_Ninc_vs_thickness_vs_b = (TH3D*)fin->Get("h3d_event_Ninc_vs_thickness_vs_b");
  TH3D* h3d_event_Nincch_vs_thickness_vs_b = (TH3D*)fin->Get("h3d_event_Nincch_vs_thickness_vs_b");

  // 2d plots
  TH2D* h2d_event_thickness_vs_b = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("xy");
  TH2D* h2d_event_Ninc_vs_thickness = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("yx");
  TH2D* h2d_event_Ninc_vs_b = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("zx");

  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 10;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 10;

    TH2F* htemp = new TH2F("htemp","",20,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("thickness [fm]");
    htemp->GetYaxis()->SetTitle("b [fm]");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h2d_event_thickness_vs_b->Draw("colz");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sthickness_vs_b.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }

}
