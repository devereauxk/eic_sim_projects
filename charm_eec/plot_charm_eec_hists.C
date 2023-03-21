R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 40, 5};
static double pt_hi[ptbin] = {10, 20, 40, 60, 60};
const int pt_color[9] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed, kCyan+1, kAzure+7, kViolet+7, kBlack};

const int etabin = 6; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5, -1, 0};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5, 0, 1};
const int eta_color[9] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed, kCyan+1, kAzure+7, kViolet+7, kBlack};

TH1D* h1d_jet_pt[etabin] = {};
TH1D* h1d_jet_eta = NULL;

TH1D* h1d_jet_eec[etabin][ptbin] = {};
TH1D* h1d_jet_eec_norm[etabin][ptbin] = {};
TH1D* h1d_jet_eec_baseline[etabin][ptbin] = {};
TH1D* h1d_jet_eec_baseline_norm[etabin][ptbin] = {};

TH1D* h1d_jet_eec_rlsqrtpt[etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt_baseline[etabin][ptbin] = {};

TH2D* h2d_jet_Q2_x[etabin][ptbin] = {};

TH1D* h1d_part_pt[etabin] = {};
TH1D* h1d_part_eta[ptbin] = {};
TH1D* h1d_part_mult = NULL;
TH1D* h1d_fixed_event_mult = NULL;
TH1D* h1d_fixed_jet_mult = NULL;

static int cno = 0;

void hists_to_csv(const char* outfile_name, vector<TH1*> hists)
{
   ofstream outfile;
   outfile.open(Form("%s", outfile_name));

   int n = hists[0]->GetNbinsX();

   for (int i = 1; i <= n; i++) // loop over lines
   {
     outfile << hists[0]->GetBinCenter(i) << "," << hists[0]->GetBinWidth(i);

     for (int ihist = 0; ihist < hists.size(); ihist++)
     {
       outfile << "," << hists[ihist]->GetBinContent(i) << "," << hists[ihist]->GetBinError(i);
     }
     outfile << "\n";

   }
   outfile.close();
}

void individual_hists(const char* out_dir)
{
  TLegend* leg = new TLegend(0.61,0.7,0.91,0.82);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetFillStyle(0);
  leg->SetMargin(0.1);

  // 1d jet pt histogram
  mclogy(cno++);
  {
    for (int ieta = 0; ieta < 3; ieta++)
    {
      h1d_jet_pt[ieta]->Draw("same");

      h1d_jet_pt[ieta]->GetXaxis()->SetRangeUser(0,40);
      h1d_jet_pt[ieta]->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      h1d_jet_pt[ieta]->GetYaxis()->SetRangeUser(1,1.3*h1d_jet_pt[3]->GetMaximum());
      h1d_jet_pt[ieta]->GetYaxis()->SetTitle("counts");
      h1d_jet_pt[ieta]->GetXaxis()->SetTitleOffset(1.3);
      h1d_jet_pt[ieta]->GetYaxis()->SetTitleOffset(1.5);
      h1d_jet_pt[ieta]->SetMarkerColor(eta_color[ieta]);
      h1d_jet_pt[ieta]->SetLineColor(eta_color[ieta]);
      h1d_jet_pt[ieta]->Draw("same hist e");
      leg->AddEntry(h1d_jet_pt[ieta],Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt.pdf\")", cno-1, out_dir) );
  }

  // 1d jet pt histogram normalized by number of jets, inclusive on eta
  mclogy(cno++);
  {
    TLegend* leg = new TLegend(0.61,0.7,0.91,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      TH1D* temp = (TH1D*) h1d_jet_pt[ieta]->Clone("temp");
      temp->Scale(1/temp->GetEntries());

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(0,40);
      temp->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->SetMarkerColor(eta_color[ieta]);
      temp->SetLineColor(eta_color[ieta]);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt_selfnorm.pdf\")", cno-1, out_dir) );
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

  // single eec plots
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      mclogxy(cno++);
      {
        TH1D* temp = (TH1D*) h1d_jet_eec[ieta][ipt]->Clone("temp");
        temp->Scale(1/temp->GetEntries());

        temp->GetXaxis()->SetTitle("R_{L}");
        temp->GetYaxis()->SetTitle("normalized EEC");
        temp->GetXaxis()->SetTitleOffset(1.3);
        temp->GetYaxis()->SetTitleOffset(1.5);
        temp->Draw("same hist e");

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_%d_%d.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

void overlay_hists(const char* out_dir)
{
  // overlay h1d_jet_eec with pt binnings, one plot per eta binning
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mclogxy(cno++);
    {
      vector<TH1*> hists;

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
        TH1D* temp = (TH1D*) h1d_jet_eec[ieta][ipt]->Clone("temp");
        hists.push_back(temp);

        temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        //temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->SetMarkerColor(pt_color[ipt]);
        temp->SetLineColor(pt_color[ipt]);
        temp->SetMarkerSize(0.5);
        temp->SetMarkerStyle(21);
        temp->Draw("same hist");

        leg->AddEntry(temp,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
      }
      leg->Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.025);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay_%d.pdf\")", cno-1, out_dir, ieta) );
      hists_to_csv(Form("%seec_overlay_%d.csv", out_dir, ieta), hists);
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
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->SetMarkerColor(pt_color[ipt]);
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->SetLineColor(pt_color[ipt]);
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->SetMarkerSize(0.5);
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->SetMarkerStyle(21);
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Draw("hist same");
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->GetXaxis()->SetTitle("R_{L} #sqrt{p_{T,jet}}");
      h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->GetYaxis()->SetTitle("EEC (no normalization)");
      leg->AddEntry(h1d_jet_eec_rlsqrtpt[etabin-1][ipt],Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
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
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Clone("ratio");
      ratio->Divide(h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt]);
      ratio->Scale(h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt]->Integral()/h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Integral());

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
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Clone("ratio");
      ratio->Add(h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt], -1);
      ratio->Scale(1/h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt]->Integral());

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
      TH1D* ratio = (TH1D*) h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Clone("ratio");
      int norm_binrange_lo = h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->FindBin(1E-2);
      int norm_binrange_hi = h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->FindBin(1);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > ratio->GetNbinsX())
      {
        norm_binrange_lo = ratio->GetNbinsX();
        cout<<"bin range hi too high; set to "<<ratio->GetNbinsX()<<endl;
      }
      cout<<"hi bin high edge "<<h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->GetBinCenter(norm_binrange_hi) + h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->GetBinWidth(norm_binrange_hi)<<endl;
      double relative_normalization =  h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi) / h1d_jet_eec_rlsqrtpt[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi);
      ratio->Scale(relative_normalization);
      ratio->Add(h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt], -1);
      ratio->Scale(1/h1d_jet_eec_rlsqrtpt_baseline[etabin-1][ipt]->Integral());

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
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > ratio->GetNbinsX())
      {
        norm_binrange_lo = ratio->GetNbinsX();
        cout<<"bin range hi too high; set to "<<ratio->GetNbinsX()<<endl;
      }
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

  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.1;
    float plot_yrange_hi = 0.2;

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
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > ratio->GetNbinsX())
      {
        norm_binrange_lo = ratio->GetNbinsX();
        cout<<"bin range hi too high; set to "<<ratio->GetNbinsX()<<endl;
      }
      cout<<"hi bin high edge "<<h1d_jet_eec[etabin-1][ipt]->GetBinCenter(norm_binrange_hi) + h1d_jet_eec[etabin-1][ipt]->GetBinWidth(norm_binrange_hi)<<endl;
      double relative_normalization =  h1d_jet_eec_baseline[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi) / h1d_jet_eec[etabin-1][ipt]->Integral(norm_binrange_lo,norm_binrange_hi);
      ratio->Scale(relative_normalization);
      ratio->Add(h1d_jet_eec_baseline[etabin-1][ipt], -1);
      ratio->Divide(h1d_jet_eec_baseline[etabin-1][ipt]);

      // plot
      ratio->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      ratio->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      ratio->GetXaxis()->SetTitle("R_{L}");
      ratio->GetYaxis()->SetTitle("% deviation from baseline");
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

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_error_from_baseline.pdf\")", cno-1, out_dir) );
  }

}

void Q2_x_panel(const char* out_dir)
{

  float plot_xrange_lo = 0.8;
  float plot_xrange_hi = 2.7;
  float plot_yrange_lo = -0.005;
  float plot_yrange_hi = 0.04;
  float plot_zrange_lo = 0;
  float plot_zrange_hi = 750E3;

  for (int ieta = 0; ieta < 3; ieta++)
  {
    for (int ipt = 0; ipt < 3; ipt++)
    {
      mclogxyz(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
      {
        TH2D* temp = (TH2D*) h2d_jet_Q2_x[ieta][ipt]->Clone();

        // plot
        //temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        //temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);
        temp->GetXaxis()->SetTitle("x_{B}");
        temp->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
        temp->Draw("colz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.028);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+p, 4*10^{8} events");
        tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_Q2_x_%i_%i.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

void particle_hists(const char* out_dir)
{
  // 1d particle pt histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TLegend* leg = new TLegend(0.51,0.7,0.81,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;

    for (int ieta = 0; ieta < 3; ieta++)
    {
      temp = (TH1D*) h1d_part_pt[ieta]->Clone("temp");
      hists.push_back(temp);

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(0,25);
      temp->GetXaxis()->SetTitle("p_{T} [GeV]");
      temp->GetYaxis()->SetRangeUser(1,1.3*h1d_part_pt[3]->GetMaximum());
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->SetMarkerColor(eta_color[ieta]);
      temp->SetLineColor(eta_color[ieta]);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_pt.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%spt.csv", out_dir), hists); hists.push_back(temp);
  }

  // 1d particle eta histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TLegend* leg = new TLegend(0.51,0.7,0.81,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;

    for (int ipt = 0; ipt < 3; ipt++)
    {
      temp = (TH1D*) h1d_part_eta[ipt]->Clone("temp");
      hists.push_back(temp);

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(-5,5);
      temp->GetXaxis()->SetTitle("#eta");
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->SetMarkerColor(pt_color[ipt]);
      temp->SetLineColor(pt_color[ipt]);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < pt < %1.1f",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_eta.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%seta.csv", out_dir), hists);
  }

  // 1d event multiplicity histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TH1D* temp = (TH1D*) h1d_part_mult->Clone("temp");
    hists.push_back(temp);

    temp->Draw("same");

    temp->GetXaxis()->SetRangeUser(0,120);
    temp->GetXaxis()->SetTitle("event multiplicity");
    temp->GetYaxis()->SetTitle("counts");
    temp->GetXaxis()->SetTitleOffset(1.3);
    temp->GetYaxis()->SetTitleOffset(1.5);
    temp->Draw("same hist e");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_mult.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%smult.csv", out_dir), hists);
  }

  // 1d D0 multiplicity per event histogram
  mclogy(cno++);
  {

    TH1D* temp = (TH1D*) h1d_fixed_event_mult->Clone("temp");

    temp->Draw("same");

    temp->GetXaxis()->SetRangeUser(0,10);
    temp->GetXaxis()->SetTitle("D0 event multiplicity");
    temp->GetYaxis()->SetTitle("counts");
    temp->GetXaxis()->SetTitleOffset(1.3);
    temp->GetYaxis()->SetTitleOffset(1.5);
    temp->Draw("same hist e");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_D0_inevent_mult.pdf\")", cno-1, out_dir) );
  }

  // 1d D0 multiplicity per event histogram
  mclogy(cno++);
  {

    TH1D* temp = (TH1D*) h1d_fixed_jet_mult->Clone("temp");

    temp->Draw("same");

    temp->GetXaxis()->SetRangeUser(0,10);
    temp->GetXaxis()->SetTitle("D0 jet multiplicity");
    temp->GetYaxis()->SetTitle("counts");
    temp->GetXaxis()->SetTitleOffset(1.3);
    temp->GetYaxis()->SetTitleOffset(1.5);
    temp->Draw("same hist e");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_D0_injet_mult.pdf\")", cno-1, out_dir) );
  }

}


void plot_charm_eec_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./", int make_ratios = 0, const char* fin_name_baseline = "")
{
  // if make_ratios == 1, uses fin_name_baseline to calculate ratios as baseline, else no ratios calculated

  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_jet_eta = (TH1D*) fin->Get("h1d_jet_eta");
  h1d_jet_eta->SetName("h1d_jet_eta");

  h1d_part_mult = (TH1D*) fin->Get("h1d_part_mult");
  h1d_fixed_event_mult = (TH1D*) fin->Get("h1d_fixed_event_mult");
  h1d_fixed_jet_mult = (TH1D*) fin->Get("h1d_fixed_jet_mult");

  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    h1d_part_eta[ipt] = (TH1D*) fin->Get(Form("h1d_part_eta_%d", ipt));
    h1d_part_eta[ipt]->SetName(Form("h1d_part_eta_%d", ipt));
  }

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    h1d_jet_pt[ieta] = (TH1D*) fin->Get(Form("h1d_jet_pt_%d", ieta));
    h1d_jet_pt[ieta]->SetName(Form("h1d_jet_pt_%d", ieta));

    h1d_part_pt[ieta] = (TH1D*) fin->Get(Form("h1d_part_pt_%d", ieta));
    h1d_part_pt[ieta]->SetName(Form("h1d_part_pt_%d", ieta));

    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      // raw data histogram
      h1d_jet_eec[ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
      h1d_jet_eec[ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d", ieta, ipt));

      h1d_jet_eec_rlsqrtpt[ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));
      h1d_jet_eec_rlsqrtpt[ieta][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));

      h2d_jet_Q2_x[ieta][ipt] = (TH2D*) fin->Get(Form("h2d_Q2_x_%d_%d", ieta, ipt));
      h2d_jet_Q2_x[ieta][ipt]->SetName(Form("h2d_Q2_x_%d_%d", ieta, ipt));
    }
  }

  // print particle level histograms
  particle_hists(out_dir);

  // print individual 2D histograms
  individual_hists(out_dir);

  // print overlay hists for: pt
  overlay_hists(out_dir);

  // print 2D Q2-x distribution
  Q2_x_panel(out_dir);

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

        h1d_jet_eec_rlsqrtpt_baseline[ieta][ipt] = (TH1D*) fin_baseline->Get(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));
        h1d_jet_eec_rlsqrtpt_baseline[ieta][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d_%d_baseline", ieta, ipt));

        // normalized histogram
        h1d_jet_eec_baseline_norm[ieta][ipt] = (TH1D*) h1d_jet_eec_baseline[ieta][ipt]->Clone(Form("h1d_jet_eec_%d_%d_norm", ieta, ipt));
        h1d_jet_eec_baseline_norm[ieta][ipt]->Scale(1/h1d_jet_eec_baseline_norm[ieta][ipt]->Integral()); // normalization
      }
    }

    ratio_hists(out_dir);
  }

}
