R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

TH1D* h1d_jet_pt = NULL;
TH1D* h1d_jet_eec = NULL;

static int cno = 0;

void individual_hists(const char* out_dir)
{
  // 1d jet pt histogram
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 25;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 500;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("jet p_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    h1d_jet_pt->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt.pdf\")", cno-1, out_dir) );
  }

  // 1d jet eec histogram, log bins

  TCanvas* c = new TCanvas("c","c", 800, 800);
  c->SetLogx();
  //c->SetLogy();
  c->Range(0,0,1,1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);

  h1d_jet_eec->Draw("hist e");

  c->SaveAs( Form("%sh1d_jet_eec.pdf", out_dir));

  /*
  mclogx(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 1;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 70;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("R_{L}");
    htemp.GetYaxis()->SetTitle("normalized EEC");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    h1d_jet_eec->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec.pdf\")", cno-1, out_dir) );
  }
  */
}


void plot_eec_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./")
{
  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_jet_pt = (TH1D*) fin->Get("h1d_jet_pt");
  h1d_jet_pt->SetName("h1d_jet_pt");

  h1d_jet_eec = (TH1D*) fin->Get("h1d_jet_eec");
  h1d_jet_eec->SetName("h1d_jet_eec");


  // print individual 2D histograms
  individual_hists(out_dir);

}
