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


void plot_inc_hadron(const char* inFile = "inc_merged.root", const int sys_option = 0, const int energy_option = 0, const char* outDir = "figs/")
{
  mcs(-1);

  TFile* fin = new TFile(inFile,"READ");

  // load histograms
  TH3D* h3d_event_Ninc_vs_thickness_vs_b = (TH3D*)fin->Get("h3d_event_Ninc_vs_thickness_vs_b");
  TH3D* h3d_event_Nincch_vs_thickness_vs_b = (TH3D*)fin->Get("h3d_event_Nincch_vs_thickness_vs_b");

  // 1d plots
  TH1D* h1d_event_thickness = (TH1D*)h3d_event_Ninc_vs_thickness_vs_b->ProjectionY("thickness");
  TH1D* h1d_event_b = (TH1D*)h3d_event_Ninc_vs_thickness_vs_b->ProjectionZ("b");
  TH1D* h1d_event_Ninc = (TH1D*)h3d_event_Ninc_vs_thickness_vs_b->ProjectionX("Ninc");
  TH1D* h1d_event_Nincch = (TH1D*)h3d_event_Nincch_vs_thickness_vs_b->ProjectionX("Nincc");

  // 2d plots
  TH2D* h2d_event_thickness_vs_b = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("yz");
  TH2D* h2d_event_Ninc_vs_thickness = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("xy");
  TH2D* h2d_event_Ninc_vs_b = (TH2D*)h3d_event_Ninc_vs_thickness_vs_b->Project3D("xz");
  TH2D* h2d_event_Nincch_vs_thickness = (TH2D*)h3d_event_Nincch_vs_thickness_vs_b->Project3D("xy");
  TH2D* h2d_event_Nincch_vs_b = (TH2D*)h3d_event_Nincch_vs_thickness_vs_b->Project3D("xz");

  // profiles of 2d plots
  TProfile* prof_thickness_in_b = (TProfile*)h2d_event_thickness_vs_b->ProfileX("prof_thickness_in_b");
  TProfile* prof_Ninc_in_thickness = (TProfile*)h2d_event_Ninc_vs_thickness->ProfileX("prof_Ninc_in_thickness");
  TProfile* prof_Ninc_in_b = (TProfile*)h2d_event_Ninc_vs_b->ProfileX("prof_Ninc_in_b");
  TProfile* prof_Nincch_in_thickness = (TProfile*)h2d_event_Nincch_vs_thickness->ProfileX("prof_Nincch_in_thickness");
  TProfile* prof_Nincch_in_b = (TProfile*)h2d_event_Nincch_vs_b->ProfileX("prof_Nincch_in_b");

  // setting histogram bounds
  float thickness_lo = 0;
  float thickness_hi = 6;
  float b_lo = 0;
  float b_hi = 6;
  float Ninc_lo = 0;
  float Ninc_hi = 6;
  float Nincch_lo = 0;
  float Nincch_hi = 6;

  // writing profiles to outfile
  TFile* fout = new TFile(Form("%sinc_hists_gen.root", outDir),"recreate");
  prof_thickness_in_b->Write();
  prof_Ninc_in_thickness->Write();
  prof_Ninc_in_b->Write();
  prof_Nincch_in_thickness->Write();
  prof_Nincch_in_b->Write();
  fout->Write();


  mcs(cno++);
  {
    float plot_xrange_lo = thickness_lo;
    float plot_xrange_hi = thickness_hi;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.2*h1d_event_thickness->GetMaximum();

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("thickness [fm]");
    htemp->GetYaxis()->SetTitle("counts");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h1d_event_thickness->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sthickness.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }
  mcs(cno++);
  {
    float plot_xrange_lo = b_lo;
    float plot_xrange_hi = b_hi;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.2*h1d_event_b->GetMaximum();

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("b [fm]");
    htemp->GetYaxis()->SetTitle("counts");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h1d_event_b->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sb.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }
  mcs(cno++);
  {
    float plot_xrange_lo = Ninc_lo;
    float plot_xrange_hi = Ninc_hi;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.2*h1d_event_Ninc->GetMaximum();

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("N_{inc}");
    htemp->GetYaxis()->SetTitle("counts");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h1d_event_Ninc->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNinc.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }
  mcs(cno++);
  {
    float plot_xrange_lo = Nincch_lo;
    float plot_xrange_hi = Nincch_hi;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.2*h1d_event_Nincch->GetMaximum();

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("N_{inc}");
    htemp->GetYaxis()->SetTitle("counts");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h1d_event_Nincch->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNincch.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }

  // ############################################################################

  mcs(cno++);
  {
    float plot_xrange_lo = b_lo;
    float plot_xrange_hi = b_hi;

    float plot_yrange_lo = thickness_lo;
    float plot_yrange_hi = thickness_hi;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("b [fm]");
    htemp->GetYaxis()->SetTitle("thickness [fm]");
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
  mcs(cno++);
  {
    float plot_xrange_lo = b_lo;
    float plot_xrange_hi = b_hi;

    float plot_yrange_lo = thickness_lo;
    float plot_yrange_hi = thickness_hi;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("b [fm]");
    htemp->GetYaxis()->SetTitle("thickness [fm]");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    prof_thickness_in_b->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sthickness_vs_b_profile.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }

  mcs(cno++);
  {
    float plot_xrange_lo = thickness_lo;
    float plot_xrange_hi = thickness_hi;

    float plot_yrange_lo = Ninc_lo;
    float plot_yrange_hi = Ninc_hi;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("thickness [fm]");
    htemp->GetYaxis()->SetTitle("N_{inc}");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h2d_event_Ninc_vs_thickness->Draw("colz");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNinc_vs_thickness.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }
  mcs(cno++);
  {
    float plot_xrange_lo = thickness_lo;
    float plot_xrange_hi = thickness_hi;

    float plot_yrange_lo = Ninc_lo;
    float plot_yrange_hi = 2;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("thickness [fm]");
    htemp->GetYaxis()->SetTitle("N_{inc}");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    prof_Ninc_in_thickness->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNinc_vs_thickness_profile.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }

  mcs(cno++);
  {
    float plot_xrange_lo = b_lo;
    float plot_xrange_hi = b_hi;

    float plot_yrange_lo = Ninc_lo;
    float plot_yrange_hi = Ninc_hi;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("b [fm]");
    htemp->GetYaxis()->SetTitle("N_{inc}");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    h2d_event_Ninc_vs_b->Draw("colz");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNinc_vs_b.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }
  mcs(cno++);
  {
    float plot_xrange_lo = b_lo;
    float plot_xrange_hi = b_hi;

    float plot_yrange_lo = Ninc_lo;
    float plot_yrange_hi = 2;

    TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp->Draw();
    htemp->GetXaxis()->SetTitle("b [fm]");
    htemp->GetYaxis()->SetTitle("N_{inc}");
    myhset(htemp, 1.2, 1.6, 0.05, 0.045);

    prof_Ninc_in_b->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_option],energy_name[energy_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sNinc_vs_b_profile.pdf\")", cno-1, outDir) );

    delete htemp;
    delete tl;
  }

}
