R__LOAD_LIBRARY(libeicsmear);
#include "bins_fine.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

TH2D* h2d_ztheo_vs_zjet[Q2bin][etabin][processbin] = {0};

static int cno = 0;

void individual_hists(const char* out_dir)
{
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      for (int iprocess = 0; iprocess < processbin; ++iprocess)
      {
        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 0;
          float plot_yrange_hi = 1;

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.SetTitle(Form("Q2: [%f,%f], eta: [%f,%f]", Q2_lo[iQ2], Q2_hi[iQ2], eta_lo[ieta], eta_hi[ieta]));
          htemp.GetXaxis()->SetTitle("z (true)");
          htemp.GetYaxis()->SetTitle("z (calculated)");
          myhset(&htemp,1.2,1.6,0.05,0.05);

          h2d_ztheo_vs_zjet[iQ2][ieta][iprocess]->Draw("hsame");

          gROOT->ProcessLine( Form("cc%d->Print(\"%sz_def_%d_%d_%d.pdf\")", cno-1, out_dir, iQ2, ieta, iprocess) );
        }
      }
    }
  }
}


void plot_z_def_hists(const char* fin_name = "hists.root", const char* out_dir = "./")
{
  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      for (int iprocess = 0; iprocess < processbin; ++iprocess)
      {
        h2d_ztheo_vs_zjet[iQ2][ieta][iprocess] = (TH2D*) fin->Get(Form("h2d_z_frag_%d_%d_%d", iQ2, ieta, iprocess));
        h2d_ztheo_vs_zjet[iQ2][ieta][iprocess]->SetName(Form("h2d_z_frag_%d_%d_%d", iQ2, ieta, iprocess));
      }
    }
  }

  // print individual 2D histograms
  individual_hists(out_dir);

}
