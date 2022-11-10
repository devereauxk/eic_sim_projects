R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 40, 5};
static double pt_hi[ptbin] = {10, 20, 40, 60, 60};
const int pt_color[ptbin] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed};

const int etabin = 4; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5};
const int eta_color[etabin] = {kGreen+1, kBlue, kViolet, kOrange+1};

TH1D* h1d_jet_pt = NULL;
TH1D* h1d_jet_eta = NULL;

TH1D* h1d_jet_eec[etabin][ptbin] = {};
TH1D* h1d_jet_eec_norm[etabin][ptbin] = {};
TH1D* h1d_jet_eec_baseline[etabin][ptbin] = {};
TH1D* h1d_jet_eec_baseline_norm[etabin][ptbin] = {};

TH1D* h1d_jet_eec_rlsqrtpt[ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt_baseline[ptbin] = {};

static int cno = 0;

void individual_hists(const char* out_dir)
{
  // 1d jet pt histogram
  mclogy(cno++);
  {
    h1d_jet_pt->Draw("same");

    h1d_jet_pt->GetXaxis()->SetRangeUser(0,70);
    h1d_jet_pt->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    h1d_jet_pt->GetYaxis()->SetTitle("counts");
    h1d_jet_pt->GetXaxis()->SetTitleOffset(1.3);
    h1d_jet_pt->GetYaxis()->SetTitleOffset(1.5);
    h1d_jet_pt->Draw("same hist e");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt.pdf\")", cno-1, out_dir) );
  }

  // 1d jet eta histogram
  mcs(cno++);
  {
    float plot_xrange_lo = -5;
    float plot_xrange_hi = 5;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = h1d_jet_eta->GetMaximum()*1.15;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("jet eta");
    htemp.GetYaxis()->SetTitle("counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    h1d_jet_eta->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eta.pdf\")", cno-1, out_dir) );
  }

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      mclogxy(cno++);
      {
        h1d_jet_eec_norm[ieta][ipt]->GetXaxis()->SetTitle("R_{L}");
        h1d_jet_eec_norm[ieta][ipt]->GetYaxis()->SetTitle("normalized EEC");
        h1d_jet_eec_norm[ieta][ipt]->GetXaxis()->SetTitleOffset(1.3);
        h1d_jet_eec_norm[ieta][ipt]->GetYaxis()->SetTitleOffset(1.5);
        h1d_jet_eec_norm[ieta][ipt]->Draw("same hist e");

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_%d_%d.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

void overlay_hists(const char* out_dir)
{
  // overlay h1d_jet_pt with pt binnings, one plot per eta binning
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mclogxy(cno++);
    {
      float plot_xrange_lo = 1E-2;
      float plot_xrange_hi = 1;
      float plot_yrange_lo = 1E-5;
      float plot_yrange_hi = 5E-1;

      TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.025);
      leg->SetFillStyle(0);
      leg->SetMargin(0.1);

      for (int ipt = 0; ipt < ptbin-1; ipt++)
      {
        h1d_jet_eec_norm[ieta][ipt]->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        h1d_jet_eec_norm[ieta][ipt]->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        h1d_jet_eec_norm[ieta][ipt]->SetMarkerColor(pt_color[ipt]);
        h1d_jet_eec_norm[ieta][ipt]->SetLineColor(pt_color[ipt]);
        h1d_jet_eec_norm[ieta][ipt]->SetMarkerSize(0.5);
        h1d_jet_eec_norm[ieta][ipt]->SetMarkerStyle(21);
        h1d_jet_eec_norm[ieta][ipt]->Draw("same hist e");
        leg->AddEntry(h1d_jet_eec_norm[ieta][ipt],Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
      }
      leg->Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.025);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay_%d.pdf\")", cno-1, out_dir, ieta) );
    }
  }

  // overlay h1d_jet_eec_rlsqrtpt with pt binnings
  mclogxy(cno++);
  {
    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-1; ipt++)
    {
      h1d_jet_eec_rlsqrtpt[ipt]->SetMarkerColor(pt_color[ipt]);
      h1d_jet_eec_rlsqrtpt[ipt]->SetLineColor(pt_color[ipt]);
      h1d_jet_eec_rlsqrtpt[ipt]->SetMarkerSize(0.5);
      h1d_jet_eec_rlsqrtpt[ipt]->SetMarkerStyle(21);
      h1d_jet_eec_rlsqrtpt[ipt]->Draw("hist same");
      h1d_jet_eec_rlsqrtpt[ipt]->GetXaxis()->SetTitle("R_{L} #sqrt{p_{T,jet}}");
      h1d_jet_eec_rlsqrtpt[ipt]->GetYaxis()->SetTitle("EEC (no normalization)");
      leg->AddEntry(h1d_jet_eec_rlsqrtpt[ipt],Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_overlay.pdf\")", cno-1, out_dir) );
  }
}

void ratio_hists(const char* out_dir)
{
  // ratio hists for h1d_jet_eec, one plot per eta binning
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mclogxy(cno++);
    {
      float plot_xrange_lo = 1E-2;
      float plot_xrange_hi = 1;
      float plot_yrange_lo = 0.4;
      float plot_yrange_hi = 2.5;

      TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.025);
      leg->SetFillStyle(0);
      leg->SetMargin(0.1);

      for (int ipt = 0; ipt < ptbin-2; ipt++)
      {
        // calculate ratio
        TH1D* ratio = (TH1D*) h1d_jet_eec_norm[ieta][ipt]->Clone("ratio");
        ratio->Divide(h1d_jet_eec_baseline_norm[ieta][ipt]);

        // plot
        ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        ratio->GetYaxis()->SetTitle("normalized EEC; energy loss on / off");
        ratio->SetMarkerColor(pt_color[ipt]);
        ratio->SetLineColor(pt_color[ipt]);
        ratio->SetMarkerSize(0.5);
        ratio->SetMarkerStyle(21);
        ratio->Draw("same hist e");
        leg->AddEntry(ratio,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
      }
      leg->Draw("same");

      TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
      l1.SetLineStyle(7);
      l1.SetLineColor(kGray+2);
      l1.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.025);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_ratio_%d.pdf\")", cno-1, out_dir, ieta) );
    }
  }

  // ratio hists for h1d_jet_eec_rlsqrtpt, on / off, on and off normalized
  mclogxy(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = 0.5;
    float plot_yrange_hi = 5;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[ipt]->Clone("ratio");
      ratio->Divide(h1d_jet_eec_rlsqrtpt_baseline[ipt]);
      ratio->Scale(h1d_jet_eec_rlsqrtpt_baseline[ipt]->Integral()/h1d_jet_eec_rlsqrtpt[ipt]->Integral());

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      ratio->GetYaxis()->SetTitle("normalized EEC; energy loss on / off");
      ratio->SetMarkerColor(pt_color[ipt]);
      ratio->SetLineColor(pt_color[ipt]);
      ratio->SetMarkerSize(0.5);
      ratio->SetMarkerStyle(21);
      ratio->Draw("same hist");
      leg->AddEntry(ratio,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_ratio.pdf\")", cno-1, out_dir) );
  }

  // ratio hists for h1d_jet_eec_rlsqrtpt, (on - off) / int dR_L off
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.06;
    float plot_yrange_hi = 0.02;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[ipt]->Clone("ratio");
      ratio->Add(h1d_jet_eec_rlsqrtpt_baseline[ipt], -1);
      ratio->Scale(1/h1d_jet_eec_rlsqrtpt_baseline[ipt]->Integral());

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      ratio->GetYaxis()->SetTitle("normalized EEC (on - off)");
      ratio->SetMarkerColor(pt_color[ipt]);
      ratio->SetLineColor(pt_color[ipt]);
      ratio->SetMarkerSize(0.5);
      ratio->SetMarkerStyle(21);
      ratio->Draw("same hist");
      leg->AddEntry(ratio,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_ratio_shifted.pdf\")", cno-1, out_dir) );
  }

  // ratio hists for h1d_jet_eec_rlsqrtpt, (relative normalization * on - off) / int dR_L off
  // Determine “relative normalization” by making sure that K=10 and K=0 are on top of each other in the region where we know there is no modification (small angle region).
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.02;
    float plot_yrange_hi = 0.07;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[ipt]->Clone("ratio");
      float norm_range_lo = 1E-2;
      float norm_range_hi = 1E-1;
      double relative_normalization =  h1d_jet_eec_rlsqrtpt_baseline[ipt]->Integral(norm_range_lo,norm_range_hi) / h1d_jet_eec_rlsqrtpt[ipt]->Integral(norm_range_lo,norm_range_hi);
      cout<<"relative normalization: "<<relative_normalization<<endl;
      ratio->Scale(relative_normalization);
      ratio->Add(h1d_jet_eec_rlsqrtpt_baseline[ipt], -1);
      ratio->Scale(1/h1d_jet_eec_rlsqrtpt_baseline[ipt]->Integral());

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      ratio->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      ratio->SetMarkerColor(pt_color[ipt]);
      ratio->SetLineColor(pt_color[ipt]);
      ratio->SetMarkerSize(0.5);
      ratio->SetMarkerStyle(21);
      ratio->Draw("same hist");
      leg->AddEntry(ratio,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_ratio_shifted_relnorm.pdf\")", cno-1, out_dir) );
  }


}


void plot_eec_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./", int make_ratios = 0, const char* fin_name_baseline = "")
{
  // if make_ratios == 1, uses fin_name_baseline to calculate ratios as baseline, else no ratios calculated

  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_jet_pt = (TH1D*) fin->Get("h1d_jet_pt");
  h1d_jet_pt->SetName("h1d_jet_pt");
  h1d_jet_eta = (TH1D*) fin->Get("h1d_jet_eta");
  h1d_jet_eta->SetName("h1d_jet_eta");

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      // raw data histogram
      h1d_jet_eec[ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
      h1d_jet_eec[ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d", ieta, ipt));

      // normalized histogram
      h1d_jet_eec_norm[ieta][ipt] = (TH1D*) h1d_jet_eec[ieta][ipt]->Clone(Form("h1d_jet_eec_%d_%d_norm", ieta, ipt));
      h1d_jet_eec_norm[ieta][ipt]->Scale(1/h1d_jet_eec_norm[ieta][ipt]->Integral()); // normalization
    }
  }

  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    // raw data histogram
    h1d_jet_eec_rlsqrtpt[ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d", ipt));
    h1d_jet_eec_rlsqrtpt[ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d", ipt));
  }


  // print individual 2D histograms
  individual_hists(out_dir);

  // print overlay hists for: pt
  overlay_hists(out_dir);

  // print ratio hists for EEC over py, needs baseline file
  if (make_ratios == 1)
  {
    TFile* fin_baseline = new TFile(fin_name_baseline, "READ");

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        // raw data histogram
        h1d_jet_eec_baseline[ieta][ipt] = (TH1D*) fin_baseline->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
        h1d_jet_eec_baseline[ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d_baseline", ieta, ipt));

        // normalized histogram
        h1d_jet_eec_baseline_norm[ieta][ipt] = (TH1D*) h1d_jet_eec_baseline[ieta][ipt]->Clone(Form("h1d_jet_eec_%d_%d_norm", ieta, ipt));
        h1d_jet_eec_baseline_norm[ieta][ipt]->Scale(1/h1d_jet_eec_baseline_norm[ieta][ipt]->Integral()); // normalization
      }
    }

    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h1d_jet_eec_rlsqrtpt_baseline[ipt] = (TH1D*) fin_baseline->Get(Form("h1d_jet_eec_rlsqrtpt_%d", ipt));
      h1d_jet_eec_rlsqrtpt_baseline[ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d_baseline", ipt));
    }

    ratio_hists(out_dir);
  }

}
