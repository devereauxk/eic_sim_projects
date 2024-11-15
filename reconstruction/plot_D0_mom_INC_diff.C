R__LOAD_LIBRARY(libeicsmear);
#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int sys_bins = 2;
const char* fin_dirs[sys_bins] = {"./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eAu_10_100_tauforOff_qhat0_nlo/outForPythiaMode/"};
const char* sys_name[sys_bins] = {"e+Au INC on", "e+Au INC off"};
const int sys_color[sys_bins] = {kBlack, kRed};
// e+Au
const char* sys_abbr[sys_bins] = {"eAuINCon", "eAuINCff"};
const double sys_cross_section[sys_bins] = {3.7930063E-02, 3.5348639E-02};
// e+C
//const char* sys_abbr[sys_bins] = {"eCINCon", "eCINCoff"};
//const double sys_cross_section[sys_bins] = {3.4899661E-02, 3.7974608E-02};

TH1D* D0_p[sys_bins] = {0};
TH1D* D0_pt[sys_bins] = {0};

TH1D* D0_p_diff;
TH1D* D0_pt_diff;
TH1D* D0_p_ratio;
TH1D* D0_pt_ratio;

static int cno = 0;

void standardLatex()
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.50,0.84,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8");
}

void individual_hists()
{
  for(int isys = 0; isys < sys_bins; isys++)
  {
    mclogy(cno++); // p
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 30;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1400000;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
      htemp.GetYaxis()->SetTitle("counts");
      myhset(&htemp,1.2,1.6,0.05,0.05);

      TLegend leg(0.55,0.66,0.84,0.82);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      D0_p[isys]->SetMarkerColor(sys_color[isys]);
      D0_p[isys]->SetLineColor(sys_color[isys]);
      D0_p[isys]->Draw("hsame");
      leg.AddEntry(D0_p[isys],Form("%s", sys_name[isys]), "l");

      leg.Draw("hsame");

      standardLatex();

      gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_p_logy.pdf\")", cno-1, fin_dirs[isys]) );
    }
    mclogy(cno++); // pt
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 6;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1400000;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
      htemp.GetYaxis()->SetTitle("counts");
      myhset(&htemp,1.2,1.6,0.05,0.05);

      TLegend leg(0.55,0.66,0.84,0.82);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      D0_pt[isys]->SetMarkerColor(sys_color[isys]);
      D0_pt[isys]->SetLineColor(sys_color[isys]);
      D0_pt[isys]->Draw("hsame");
      leg.AddEntry(D0_pt[isys],Form("%s", sys_name[isys]), "l");

      leg.Draw("hsame");

      standardLatex();

      gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_pt_logy.pdf\")", cno-1, fin_dirs[isys]) );
    }
  }
}

void overlay_hists(const char* out_dir = "./", const char* label = "e+Au, 1E4 events")
{
  mclogy(cno++); // p
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 30;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 14000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}");
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

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_p_diff_logy_noscale.pdf\")", cno-1, out_dir) );
  }
  mclogy(cno++); // pt
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 6;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 14000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    for (int isys = 0; isys < sys_bins; isys++)
    {
      D0_pt[isys]->SetMarkerColor(sys_color[isys]);
      D0_pt[isys]->SetLineColor(sys_color[isys]);
      D0_pt[isys]->Draw("hsame");
      leg.AddEntry(D0_pt[isys],Form("%s", sys_name[isys]), "l");
    }

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_pt_diff_logy_noscale.pdf\")", cno-1, out_dir) );
  }

  // difference hists
  mcs(cno++); // p
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 30;

    float plot_yrange_lo = -4000;
    float plot_yrange_hi = 2000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}_{INC on} - N^{D^{0}}_{INC off}");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    D0_p_diff->SetMarkerColor(kRed);
    D0_p_diff->SetLineColor(kRed);
    D0_p_diff->Draw("hsame");
    leg.AddEntry(D0_p_diff, label, "l");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("hsame");

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_p_diff.pdf\")", cno-1, out_dir) );
  }
  mcs(cno++); // pt
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 6;

    float plot_yrange_lo = -4000;
    float plot_yrange_hi = 2000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}_{INC on} - N^{D^{0}}_{INC off}");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    D0_pt_diff->SetMarkerColor(kRed);
    D0_pt_diff->SetLineColor(kRed);
    D0_pt_diff->Draw("hsame");
    leg.AddEntry(D0_pt_diff, label, "l");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("hsame");

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_pt_diff.pdf\")", cno-1, out_dir) );
  }

  // ratio hists
  mcs(cno++); // p
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 30;

    float plot_yrange_lo = 0.5;
    float plot_yrange_hi = 1.5;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}_{INC on} / N^{D^{0}}_{INC off}");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    D0_p_ratio->SetMarkerColor(kRed);
    D0_p_ratio->SetLineColor(kRed);
    D0_p_ratio->Draw("hsame");
    leg.AddEntry(D0_p_diff, label, "l");

    TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("hsame");

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_p_ratio.pdf\")", cno-1, out_dir) );
  }
  mcs(cno++); // pt
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 6;

    float plot_yrange_lo = 0.25;
    float plot_yrange_hi = 1.25;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("N^{D^{0}}_{INC on} / N^{D^{0}}_{INC off}");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    D0_pt_ratio->SetMarkerColor(kRed);
    D0_pt_ratio->SetLineColor(kRed);
    D0_pt_ratio->Draw("hsame");
    leg.AddEntry(D0_pt_diff, label, "l");

    TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("hsame");

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_pt_ratio.pdf\")", cno-1, out_dir) );
  }

}

void plot_D0_mom_INC_diff(const char* fin_name = "D0_mom_INC_diff.root", const char* overlay_out_dir = "./", const char* label = "e+Au, 1E4 events")
{
  mclogy(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  for(int isys = 0; isys < sys_bins; isys++)
  {
    D0_p[isys] = (TH1D*) fin->Get(Form("h1d_D0_p_%s", sys_abbr[isys]));
    D0_p[isys]->SetName(Form("h1d_D0_p_%s", sys_abbr[isys]));

    D0_pt[isys] = (TH1D*) fin->Get(Form("h1d_D0_pt_%s", sys_abbr[isys]));
    D0_pt[isys]->SetName(Form("h1d_D0_pt_%s", sys_abbr[isys]));
  }

  // print histograms for individual systems
  individual_hists();


  // normalize histograms st [1] has same integral as [0]
  /*
    D0_p[1]->Scale(D0_p[0]->Integral() / D0_p[1]->Integral());
    D0_pt[1]->Scale(D0_pt[0]->Integral() / D0_pt[1]->Integral());
  */

  // normalize histograms st [1] is scaled to have the same cross section as [0]
  D0_p[1]->Scale(sys_cross_section[0] / sys_cross_section[1]);
  D0_pt[1]->Scale(sys_cross_section[0] / sys_cross_section[1]);


  //generate diff and ratio hists
  // p
  D0_p_diff = (TH1D*) D0_p[0]->Clone("D0_p_diff");
  D0_p_diff->Add(D0_p[1], -1);
  D0_p_ratio = (TH1D*) D0_p[0]->Clone("D0_p_ratio");
  D0_p_ratio->Divide(D0_p[1]);

  // pt
  D0_pt_diff = (TH1D*) D0_pt[0]->Clone("D0_pt_diff");
  D0_pt_diff->Add(D0_pt[1], -1);
  D0_pt_ratio = (TH1D*) D0_pt[0]->Clone("D0_pt_ratio");
  D0_pt_ratio->Divide(D0_pt[1]);


  // print overlayed histograms for all systems
  overlay_hists(overlay_out_dir, label);


}
