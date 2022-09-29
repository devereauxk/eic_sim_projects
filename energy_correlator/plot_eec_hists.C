R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int ptbin = 3; // inclusive on last bin, inclusive on lower limit, explusive on upper
static double pt_lo[ptbin] = {5, 10, 5};
static double pt_hi[ptbin] = {10, 20, 20};
const int pt_color[ptbin] = {kGreen+1, kBlue};

TH1D* h1d_jet_eec[ptbin] = {};
TH1D* h1d_jet_pt = NULL;

static int cno = 0;

void individual_hists(const char* out_dir)
{
  // 1d jet pt histogram
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 30;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = h1d_jet_pt->GetMaximum()*1.15;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("jet p_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    h1d_jet_pt->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt.pdf\")", cno-1, out_dir) );
  }

  // 1d jet eec histogram, log bins

  /*
  TCanvas* c = new TCanvas("c","c", 800, 800);
  c->SetLogx();
  c->SetLogy();
  c->Range(0,0,1,1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);

  h1d_jet_eec->Draw("hist e");

  c->SaveAs( Form("%sh1d_jet_eec.pdf", out_dir));
  */
  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    mclogxy(cno++);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 1;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = h1d_jet_eec[ipt]->GetMaximum()*1.50;

      /*
      TH2F htemp("htemp","",50,plot_xrange_lo,plot_xrange_hi,50,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw("hsame");
      htemp.GetXaxis()->SetTitle("R_{L}");
      htemp.GetYaxis()->SetTitle("normalized EEC");
      myhset(&htemp,1.2,1.6,0.05,0.05);
      */

      h1d_jet_eec[ipt]->GetXaxis()->SetTitle("R_{L}");
      h1d_jet_eec[ipt]->GetYaxis()->SetTitle("normalized EEC");
      h1d_jet_eec[ipt]->GetXaxis()->SetTitleOffset(1.3);
      h1d_jet_eec[ipt]->GetYaxis()->SetTitleOffset(1.5);
      h1d_jet_eec[ipt]->Draw("same hist e");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_%d.pdf\")", cno-1, out_dir, ipt) );
    }
  }

}

void overlay_hists(const char* out_dir)
{
  mclogxy(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 1;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = h1d_jet_eec[ptbin-1]->GetMaximum()*1.50;

    TH2F htemp("htemp","",50,plot_xrange_lo,plot_xrange_hi,50,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("R_{L}");
    htemp.GetYaxis()->SetTitle("normalized EEC");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend* leg = new TLegend(0.21,0.17,0.51,0.29);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-1; ipt++)
    {
      h1d_jet_eec[ipt]->SetMarkerColor(pt_color[ipt]);
      h1d_jet_eec[ipt]->SetLineColor(pt_color[ipt]);
      h1d_jet_eec[ipt]->SetMarkerSize(0.5);
      h1d_jet_eec[ipt]->SetMarkerStyle(21);
      h1d_jet_eec[ipt]->Draw("same hist e");
      leg->AddEntry(h1d_jet_eec[ipt],Form("%.1f < p_{T} < %.1f",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay.pdf\")", cno-1, out_dir) );
  }
}


void plot_eec_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./")
{
  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_jet_pt = (TH1D*) fin->Get("h1d_jet_pt");
  h1d_jet_pt->SetName("h1d_jet_pt");

  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    h1d_jet_eec[ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d", ipt));
    h1d_jet_eec[ipt]->SetName(Form("h1d_jet_eec_%d", ipt));

    h1d_jet_eec[ipt]->Scale(1/h1d_jet_eec[ipt]->Integral()); // normalization
  }

  // print individual 2D histograms
  individual_hists(out_dir);

  // print overlay hists for: pt
  overlay_hists(out_dir);

}
