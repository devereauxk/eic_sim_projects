#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 4;
const char* sys_name[sys_bins] = {"e+C", "e+Cu", "e+Au", "e+Pb"};
const char* sys_abbr[sys_bins] = {"eC", "eCu", "eAu", "ePb"};
const int sys_color[sys_bins] = {kBlack, kRed, kOrange+1, kGreen+1};

const int energy_bins = 6;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

TProfile* prof_thickness_in_b[sys_bins] = {0};
TProfile* prof_Ninc_in_thickness[sys_bins] = {0};
TProfile* prof_Ninc_in_b[sys_bins] = {0};
TProfile* prof_Nincch_in_thickness[sys_bins] = {0};
TProfile* prof_Nincch_in_b[sys_bins] = {0};

static int cno = 0;

void standardLatex()
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.19,0.17,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8");
}

void plot_comparison(const int energy_option = 1)
{
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 1;

    float plot_yrange_lo = 0.;
    float plot_yrange_hi = 2;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("b [fm]");
    htemp.GetYaxis()->SetTitle("thickness [fm]");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.71,0.84,0.87);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    for (int isys = 0; isys < sys_bins; isys++)
    {
      prof_thickness_in_b[isys]->SetMarkerStyle(21);
      prof_thickness_in_b[isys]->SetMarkerSize(0.7);
      prof_thickness_in_b[isys]->SetMarkerColor(sys_color[isys]);
      prof_thickness_in_b[isys]->SetLineColor(sys_color[isys]);
      prof_thickness_in_b[isys]->Draw("same");
      leg.AddEntry(prof_thickness_in_b[isys],Form("%s @ %s", sys_name[isys], energy_name[isys], "p"));
    }

    leg.Draw("same");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"%scompare_thickness_vs_b.pdf\")", cno-1, outDir) );
  }


}

void compare_plot_inc_hadron_diff_sys(const int energy_option = 1, const char* outDir = "figs/")
{
  mcs(-1);

  TFile* fin[sys_bins] = {0};
  for (int isys = 0; isys < sys_bins; ++isys)
  {
    fin[isys] = new TFile(Form("inc_hists_gen_%s.root",sys_abbr[isys]),"READ");

    prof_thickness_in_b[isys] = (TProfile*)fin[isys]->Get("prof_thickness_in_b");
    prof_thickness_in_b[isys]->SetName(Form("prof_thickness_in_b_%s", sys_abbr[isys]));

    prof_Ninc_in_thickness[isys] = (TProfile*)fin[isys]->Get("prof_Ninc_in_thickness");
    prof_Ninc_in_thickness[isys]->SetName(Form("prof_Ninc_in_thickness_%s", sys_abbr[isys]));

    prof_Ninc_in_b[isys] = (TProfile*)fin[isys]->Get("prof_Ninc_in_b");
    prof_Ninc_in_b[isys]->SetName(Form("prof_Ninc_in_b_%s", sys_abbr[isys]));

    prof_Nincch_in_thickness[isys] = (TProfile*)fin[isys]->Get("prof_Nincch_in_thickness");
    prof_Nincch_in_thickness[isys]->SetName(Form("prof_Nincch_in_thickness_%s", sys_abbr[isys]));

    prof_Nincch_in_b[isys] = (TProfile*)fin[isys]->Get("prof_Nincch_in_b");
    prof_Nincch_in_b[isys]->SetName(Form("prof_Nincch_in_b_%s", sys_abbr[isys]));


    cout<<"system #"<<isys<<" works"<<endl;
  }

  plot_comparison(energy_option);

}
