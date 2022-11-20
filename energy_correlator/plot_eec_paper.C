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

const int speciesnum = 3;
static char* species[speciesnum] = {(char*)"e+C", (char*)"e+Cu", (char*)"e+Au"};

const int knum = 4;
static int k[knum] = {0,2,4,10};

const int energynum = 3;
static char* energy[energynum] = {(char*)"5 on 41 GeV", (char*)"10 on 100 GeV", (char*)"18 on 110 GeV"};

static char* fname_eC_by_K[knum] = {(char*)"./eHIJING/eC_1E8_K0/merged.root", (char*)"", (char*)"./eHIJING/eC_1E8_K4/merged.root", (char*)""};
static char* fname_eCu_by_K[knum] = {(char*)"./eHIJING/eCu_1E8_K0/merged.root", (char*)"", (char*)"./eHIJING/eCu_1E8_K4/merged.root", (char*)""};
static char* fname_eAu_by_K[knum] = {(char*)"./eHIJING/eAu_1E8_K0_condor_v2/merged.root", (char*)"./eHIJING/eAu_1E8_K2/merged.root", (char*)"./eHIJING/eAu_1E8_K4/merged.root", (char*)"./eHIJING/eAu_1E8_condor_v2/merged.root"};
static char** fname_eA_by_K[speciesnum] = {fname_eC_by_K, fname_eCu_by_K, fname_eAu_by_K};

const char* out_dir = "./";

TH1D* h1d_jet_eec[speciesnum][knum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt[speciesnum][knum][ptbin] = {};

static int cno = 0;

void pt_eta_3by3_hists()
{
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      mclogxy(cno++);
      {
        float plot_xrange_lo = 1E-2;
        float plot_xrange_hi = 1;
        float plot_yrange_lo = 1E-5;
        float plot_yrange_hi = 5E-1;
        float legend_x = 0.5;
        float legend_y = 0.2;

        TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.2);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.025);
        leg->SetFillStyle(0);
        leg->SetMargin(0.1);

        TH1D* temp;
        TH1D* temp_baseline = (TH1D*) h1d_jet_eec[2][0][ieta][ipt]->Clone();
        for (int ik = 0; ik < knum; ik++)
        {
          temp = (TH1D*) h1d_jet_eec[2][ik][ieta][ipt]->Clone();

          // calculate relative normalization ratio
          int norm_binrange_lo = temp->FindBin(1E-2);
          int norm_binrange_hi = temp->FindBin(0.2);
          double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
          temp->Scale(relative_normalization);
          temp->Add(temp_baseline, -1);
          temp->Scale(1/temp_baseline->Integral());

          // plot histogram
          temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
          temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
          temp->GetXaxis()->SetTitle("R_{L}");
          temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
          temp->SetMarkerColor(pt_color[ik]);
          temp->SetLineColor(pt_color[ik]);
          temp->SetMarkerSize(0.5);
          temp->SetMarkerStyle(21);
          temp->Draw("same hist e");
          leg->AddEntry(temp,Form("K = %i",k[ik]));
        }

        leg->Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.025);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.22,0.86,"eHIJING, e+Au @ 10+100 GeV, 10^{8} events");
        tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.82,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay_%d_%d.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

/*
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

  // ratio hists for h1d_jet_eec, (on - off) / int dR_L off, only calculated for inclusive eta bin (etabin-1)
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.05;
    float plot_yrange_hi = 0.03;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec[etabin-1][ipt]->Clone("ratio");
      ratio->Add(h1d_jet_eec_baseline[etabin-1][ipt], -1);
      ratio->Scale(1/h1d_jet_eec_baseline[etabin-1][ipt]->Integral());

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}");
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

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.025);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin-1],eta_hi[etabin-1]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_ratio_shifted.pdf\")", cno-1, out_dir) );
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
  // small angle region set manually to R_L\sqrt(jet pt) \in [1E-2, 1E-1]
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.04;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[ipt]->Clone("ratio");
      int norm_binrange_lo = h1d_jet_eec_rlsqrtpt[ipt]->FindBin(1E-2);
      int norm_binrange_hi = h1d_jet_eec_rlsqrtpt[ipt]->FindBin(1);
      cout<<"hi bin high edge "<<h1d_jet_eec_rlsqrtpt[ipt]->GetBinCenter(norm_binrange_hi) + h1d_jet_eec_rlsqrtpt[ipt]->GetBinWidth(norm_binrange_hi)<<endl;
      double relative_normalization =  h1d_jet_eec_rlsqrtpt_baseline[ipt]->Integral(norm_binrange_lo,norm_binrange_hi) / h1d_jet_eec_rlsqrtpt[ipt]->Integral(norm_binrange_lo,norm_binrange_hi);
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

  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.04;

    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      // calculate ratio
      TH1D* ratio = (TH1D*) h1d_jet_eec[etabin-1][ipt]->Clone("ratio");
      int norm_binrange_lo = h1d_jet_eec[etabin-1][ipt]->FindBin(1E-2);
      int norm_binrange_hi = h1d_jet_eec[etabin-1][ipt]->FindBin(0.2);
      cout<<"hi bin high edge "<<h1d_jet_eec[etabin-1][ipt]->GetBinCenter(norm_binrange_hi) + h1d_jet_eec[etabin-1][ipt]->GetBinWidth(norm_binrange_hi)<<endl;
      double relative_normalization =  h1d_jet_eec_baseline[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi) / h1d_jet_eec[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi);
      ratio->Scale(relative_normalization);
      ratio->Add(h1d_jet_eec_baseline[etabin-1][ipt], -1);
      ratio->Scale(1/h1d_jet_eec_baseline[etabin-1][ipt]->Integral());

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}");
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

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_ratio_shifted_relnorm.pdf\")", cno-1, out_dir) );
  }


}
*/

void plot_eec_paper()
{

  mcs(-1);

  // load histograms
  TFile* fin;
  char* fin_name;
  for (int ispecies = 0; ispecies < speciesnum; ispecies++)
  {
    for (int ik = 0; ik < knum; ik++)
    {
      fin_name = fname_eA_by_K[ispecies][ik];

      if (strcmp(fin_name,"") != 0)
      {
        fin = new TFile(fin_name, "READ");

        for (int ieta = 0; ieta < etabin; ieta++)
        {
          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            // raw data histograms
            h1d_jet_eec[ispecies][ik][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
            h1d_jet_eec[ispecies][ik][ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d_%d_%d", ieta, ipt, ispecies, ik));
          }
        }
        for (int ipt = 0; ipt < ptbin; ipt++)
        {
          // raw data histogram
          h1d_jet_eec_rlsqrtpt[ispecies][ik][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d", ipt));
          h1d_jet_eec_rlsqrtpt[ispecies][ik][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d_%d_%d", ipt, ispecies, ik));
        }
      }

    }
  }

  // plot individual panels

  pt_eta_3by3_hists();

  //nuclei_hists();

  //energy_hists();

  //power_hists();

}
