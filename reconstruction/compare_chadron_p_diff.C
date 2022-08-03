#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 2;
const char* sys_name[sys_bins] = {"e+Au INC on", "e+Au INC off"};
const char* sys_abbr[sys_bins] = {"eAuINCon", "eAuINCoff"};
const char* sys_fin[sys_bins] = {"./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/ana_merged.root", "./BeAGLE_v102/eAu_10_100_tauforOff_qhat0_nlo/outForPythiaMode/ana_merged.root"};
const int sys_color[sys_bins] = {kBlack, kRed};

TH2D* h2d_D0_p_vs_eta_gen_in_x[sys_bins][Q2bin][xbin];
TH2D* h2d_D0_pt_vs_eta_gen_in_x[sys_bins][Q2bin][xbin];

TH1D* h1d_D0_p[sys_bins] = {0};
TH1D* h1d_D0_pt[sys_bins] = {0};

static int cno = 0;

void make_histograms()
{
  for(int isys = 0; isys < sys_bins; isys++)
  {
    // initializing / projecting inclusive th2ds
    h1d_D0_p[isys] = h2d_D0_p_vs_eta_gen_in_x[isys][Q2bin-1][xbin-1]->ProjectionY("h1d_D0_p");
    h1d_D0_pt[isys] = h2d_D0_pt_vs_eta_gen_in_x[isys][Q2bin-1][xbin-1]->ProjectionY("h1d_D0_pt");

    // if isys!=0, scalling histograms to same # entries of sys=0 histograms
    if(isys!=0)
    {
      h1d_D0_p[isys]->Scale(h1d_D0_p[0]->Integral() / h1d_D0_p[isys]->Integral());
      h1d_D0_pt[isys]->Scale(h1d_D0_pt[0]->Integral() / h1d_D0_pt[isys]->Integral());
    }
  }
}

void plot_comparison()
{
    mclogy(cno++);
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

      TLegend leg(0.55,0.71,0.84,0.87);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      for (int isys = 0; isys < sys_bins; ++isys)
      {
        h1d_D0_p[isys]->SetMarkerColor(sys_color[isys]);
        h1d_D0_p[isys]->SetLineColor(sys_color[isys]);
        h1d_D0_p[isys]->Draw("same b");
        leg.AddEntry(h1d_D0_p[isys],Form("%s",sys_name[isys]),"p");
      }

      leg.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.03);
      tl->DrawLatexNDC(0.19,0.17,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

      gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_p_compare.pdf\")", cno-1) );
    }

}

void compare_chadron_p_diff()
{
  mcs(-1);

  TFile* fin[sys_bins] = {0};
  for (int isys = 0; isys < sys_bins; ++isys)
  { // skipping e+p system
    fin[isys] = new TFile(sys_fin[isys],"READ");
    for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
    {
      for (int ix = 0; ix < xbin; ++ix)
      {
        h2d_D0_p_vs_eta_gen_in_x[isys][iQ2][ix] = (TH2D*)fin[isys]->Get( Form("h2d_hadron_421_p_vs_eta_gen_in_x_%d_%d",iQ2,ix) );
        h2d_D0_p_vs_eta_gen_in_x[isys][iQ2][ix]->SetName( Form("h2d_hadron_421_p_vs_eta_gen_in_x_%s_Q2%d_x%d",sys_abbr[isys],iQ2,ix) );

        h2d_D0_pt_vs_eta_gen_in_x[isys][iQ2][ix] = (TH2D*)fin[isys]->Get( Form("h2d_hadron_421_pt_vs_eta_gen_in_x_%d_%d",iQ2,ix) );
        h2d_D0_pt_vs_eta_gen_in_x[isys][iQ2][ix]->SetName( Form("h2d_hadron_421_pt_vs_eta_gen_in_x_%s_Q2%d_x%d",sys_abbr[isys],iQ2,ix) );
      }
    }

    cout<<"system #"<<isys<<" works"<<endl;
  }

  make_histograms();

  plot_comparison();
}
