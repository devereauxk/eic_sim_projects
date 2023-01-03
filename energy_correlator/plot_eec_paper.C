R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include <string.h>
using namespace std;

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 30, 5};
static double pt_hi[ptbin] = {10, 20, 30, 60, 60};
const int pt_color[7] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed, kCyan+1, kAzure+7};

const int etabin = 4; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5};
const int eta_color[etabin] = {kGreen+1, kBlue, kViolet, kOrange+1};

const int speciesnum = 9;
static char* species[speciesnum] = {(char*)"e+p", (char*)"e+D", (char*)"e+He-3", (char*)"e+He-4", (char*)"e+C", (char*)"e+Ca", (char*)"e+Cu", (char*)"e+Au", (char*)"e+U"};
static int species_A[speciesnum] = {1, 2, 3, 4, 12, 40, 64, 197, 238};
//static double species_A16[speciesnum] = {1, 1.12246, 1.51309, 1.84931, 2, 2.41219, 2.48941};

const int knum = 4;
static int k[knum] = {0,2,4,10};

const int energynum = 3;
static char* energy[energynum] = {(char*)"5 on 41 GeV", (char*)"10 on 100 GeV", (char*)"18 on 110 GeV"};

const int pownum = 4;
static double power[pownum] = {0.5, 1, 1.5, 2};

// power=0.5
static char* fname_ep_by_K[knum] = {(char*)"./eHIJING/ep_10_100_K0_pow05/merged.root", (char*)"", (char*)"", (char*)""};
static char* fname_eD_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eD_10_100_K4_pow05/merged.root", (char*)""};
static char* fname_eHe3_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eHe3_10_100_pdf0/merged.root", (char*)""};
static char* fname_eHe4_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eHe4_10_100_pdf0/merged.root", (char*)""};
static char* fname_eC_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eC_1E8_K4_pow05/merged.root", (char*)""};
static char* fname_eCa_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eCa_10_100_K4_pow05/merged.root", (char*)""};
static char* fname_eCu_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eCu_1E8_K4_pow05/merged.root", (char*)""};
static char* fname_eAu_by_K[knum] = {(char*)"", (char*)"./eHIJING/eAu_1E8_K2_pow05/merged.root", (char*)"./eHIJING/eAu_1E8_K4_pow05/merged.root", (char*)"./eHIJING/eAu_1E8_condor_pow05/merged.root"};
static char* fname_eU_by_K[knum] = {(char*)"", (char*)"", (char*)"./eHIJING/eU_10_100_K4_pow05/merged.root", (char*)""};
static char** fname_eA_by_K[speciesnum] = {fname_ep_by_K, fname_eD_by_K, fname_eHe3_by_K, fname_eHe4_by_K, fname_eC_by_K, fname_eCa_by_K, fname_eCu_by_K, fname_eAu_by_K, fname_eU_by_K}; // K=4 cases are 2E8 events, 1E8 events otherwise

static char* fname_eA_isospin[speciesnum] {(char*)"./eHIJING/ep_10_100_K0_pow05/merged.root", (char*)"./eHIJING/eD_10_100_K4_pow05/merged.root", (char*)"./eHIJING/eHe3_10_100_pdf0/merged.root",\
  (char*)"./eHIJING/eHe4_10_100_pdf0/merged.root", (char*)"./eHIJING/eC_10_100_pdf0/merged.root", (char*)"./eHIJING/eCa_10_100_pdf0/merged.root",\
  (char*)"./eHIJING/eCu_10_100_pdf0/merged.root", (char*)"./eHIJING/eAu_10_100_pdf0/merged.root", (char*)"./eHIJING/eU_10_100_pdf0/merged.root"};

static char* fname_ep_by_power[pownum] = {(char*)"./eHIJING/ep_10_100_K0_pow05/merged.root", (char*)"./eHIJING/ep_10_100_K0/merged.root", (char*)"./eHIJING/ep_10_100_K0_pow15/merged.root", (char*)"./eHIJING/ep_10_100_K0_pow2/merged.root"}; // all cases are K=0, 2E8 events
static char* fname_eAu_by_power[pownum] = {(char*)"./eHIJING/eAu_10_100_K4_pow05/merged.root", (char*)"./eHIJING/eAu_10_100_K4/merged.root", (char*)"./eHIJING/eAu_10_100_K4_pow15/merged.root", (char*)"./eHIJING/eAu_10_100_K4_pow2/merged.root"}; // all cases are K=4, 2E8 events

static char* fname_ep_Q2_x = "./eHIJING/ep_10_100_K0_Q2x/merged.root";

const char* out_dir = "./paperplots/";

TH1D* h1d_jet_pt[speciesnum] = {}; // K=4, 10x100, 2E8 events

TH1D* h1d_jet_eec[speciesnum][knum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt[speciesnum][knum][etabin][ptbin] = {};

TH1D* h1d_jet_eec_ep_by_power[pownum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_eAu_by_power[pownum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt_ep_by_power[pownum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt_eAu_by_power[pownum][etabin][ptbin] = {};

TH1D* h1d_jet_eec_isospin[speciesnum][etabin][ptbin] = {};

TH2D* h2d_jet_Q2_x[etabin][ptbin] = {};

const float rl_norm_hi = 0.05; //0.08;
const float rl_norm_lo = 1E-3;

const float rlsqrtpt_norm_hi = rl_norm_hi * sqrt(20); //1;
const float rlsqrtpt_norm_lo = rl_norm_lo * sqrt(20); //rl_norm_lo * sqrt(20);

static int cno = 0;

// output file initialization; one file per graph


// convert an array of TH1Ds (of same binsize and range) to csv file
void hists_to_csv(const char* outfile_name, vector<TH1*> hists)
{
   ofstream outfile;
   outfile.open(Form("%s%s", out_dir, outfile_name));

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

// convert a TGraph to csv file
void graph_to_csv(const char* outfile_name, TGraph* graph)
{
  ofstream outfile;
  outfile.open(Form("%s%s", out_dir, outfile_name));

  int n = graph->GetN();
  Double_t* x = graph->GetX();
  Double_t* y = graph->GetY();

  for (int i = 0; i < n; i++) // loop over lines
  {
    outfile << x[i] << "," << y[i] << "\n";
  }
  outfile.close();

}


void pt_eta_3by3_hists()
{
  int species_pick = 7;
  int etabin_pick = 2;
  int ptbin_pick = 2;

  // using ep K=0 w/ 2E8 events the baseline

  // with R_L on the x-axis, plotting (alpha_i * K=i) / (int dR_L K=0)
  mclogxy(cno++);
  {
    vector<TH1*> hists;

    float plot_xrange_lo = 0.05;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = 1E-3;
    float plot_yrange_hi = 5E-1;
    float legend_x = 0.7;
    float legend_y = 0.2;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    // ep, K=0
    temp = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();
    temp_baseline = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();

    // calculate relative normalization ratio
    int norm_binrange_lo = temp->FindBin(rl_norm_lo);
    int norm_binrange_hi = temp->FindBin(rl_norm_hi);
    if (norm_binrange_lo == 0)
    {
      norm_binrange_lo = 1;
      cout<<"bin range lo too low; set to 1"<<endl;
    }
    if (norm_binrange_hi > temp->GetNbinsX())
    {
      norm_binrange_lo = temp->GetNbinsX();
      cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
    }
    double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
    temp->Scale(relative_normalization);
    temp->Scale(1/temp_baseline->Integral());

    // plot histogram
    temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
    temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
    temp->GetXaxis()->SetTitle("R_{L}");
    temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on)");
    temp->SetMarkerColor(pt_color[0]);
    temp->SetLineColor(pt_color[0]);
    temp->SetMarkerSize(0.5);
    temp->SetMarkerStyle(21);
    temp->Draw("same hist e");
    leg->AddEntry(temp,"e+p, K = 0");

    hists.push_back(temp);

    // eAu, K=ik
    for (int ik = 1; ik < knum; ik++)
    {
      temp = (TH1D*) h1d_jet_eec[species_pick][ik][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Scale(1/temp_baseline->Integral());

      hists.push_back(temp);

      // plot histogram
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on)");
      temp->SetMarkerColor(pt_color[ik]);
      temp->SetLineColor(pt_color[ik]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("e+Au, K = %i",k[ik]));
    }

    leg->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+Au @ 10+100 GeV, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay.pdf\")", cno-1, out_dir) );

    hists_to_csv("fig2.csv", hists);
  }

}

void pt_bin_side_by_side()
{
  int species_pick = 7;
  int etabin_pick = 2;
  int k_pick = 2;

  // using ep K=0 w/ 2E8 events the baseline

  // with R_L on the x-axis, plotting (alpha_i * K=i - K=0) / (int dR_L K=0), pt binnings for eAu forward eta selection
  /*
  mclogx(cno++);
  {
    float plot_xrange_lo = 5E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.08; //0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      temp = (TH1D*) h1d_jet_eec[species_pick][k_pick][etabin_pick][ipt]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ipt]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      for (int i = 0; i < temp->GetNbinsX(); i++)
      {
        cout<<temp->GetBinError(i)<<" ";
      }
      cout<<endl;

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[ipt]);
      temp->SetLineColor(pt_color[ipt]);
      //temp->SetMarkerSize(0.5);
      //temp->SetMarkerStyle(21);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("p_{T} #in [%.1f, %.1f)",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,Form("eHIJING, %s @ 10+100 GeV, 4*10^{8} events", species[species_pick]));
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_eAu_pt_rl.pdf\")", cno-1, out_dir) );
  }
  */

  // with R_L*sqrt(pt) on the x-axis, plotting (alpha_i * K=i - K=0) / (int dR_L K=0), pt binnings for eAu forward eta selection
  mclogx(cno++);
  {
    vector<TH1*> hists;

    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.08; //0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ipt = 0; ipt < ptbin-2; ipt++)
    {
      temp = (TH1D*) h1d_jet_eec_rlsqrtpt[species_pick][k_pick][etabin_pick][ipt]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_rlsqrtpt[0][0][etabin_pick][ipt]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rlsqrtpt_norm_lo);
      int norm_binrange_hi = temp->FindBin(rlsqrtpt_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      hists.push_back(temp);

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[ipt]);
      temp->SetLineColor(pt_color[ipt]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("p_{T} #in [%.1f, %.1f)",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,Form("eHIJING, %s @ 10+100 GeV, 4*10^{8} events", species[species_pick]));
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_eAu_pt_rlsqrtpt.pdf\")", cno-1, out_dir) );

  }

}

void baseline_comparison()
{
  int etabin_pick = 3;
  int k_pick = 0;
  const int nspecies_picks = 3;
  static int species_picks[nspecies_picks] = {0, 1, 7};

  // with R_L on the x-axis, plotting (alpha_i * K=0) / (int dR_L K=0) for eD K=0 with eAu K=0 as baseline
  // three plots - one for each pt bin, all have inclusive eta
  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    mclogxy(cno++);
    {
      vector<TH1*> hists;

      float plot_xrange_lo = 0.05;
      float plot_xrange_hi = 1;
      float plot_yrange_lo = 5E-3;
      float plot_yrange_hi = 5E-1;
      float legend_x = 0.7;
      float legend_y = 0.2;

      TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.028);
      leg->SetFillStyle(0);
      leg->SetMargin(0.1);

      TH1D* temp;

      for (int ispecies = 0; ispecies < nspecies_picks; ispecies++)
      {
        temp = (TH1D*) h1d_jet_eec[species_picks[ispecies]][k_pick][etabin_pick][ipt]->Clone();

        // calculate relative normalization ratio
        temp->Scale(1/temp->Integral());

        hists.push_back(temp);

        // plot histogram
        temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->GetXaxis()->SetTitle("R_{L}");
        temp->GetYaxis()->SetTitle("normalized EEC");
        temp->SetMarkerColor(pt_color[ispecies]);
        temp->SetLineColor(pt_color[ispecies]);
        temp->SetMarkerSize(0.5);
        temp->SetMarkerStyle(21);
        temp->Draw("same hist e");
        leg->AddEntry(temp,Form("%s, K = %i",species[species_picks[ispecies]], k[k_pick]));
      }

      leg->Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.028);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+Au @ 10+100 GeV, 4*10^{8} events");
      tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
      tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_baseline_comparison_%d.pdf\")", cno-1, out_dir, ipt) );

      hists_to_csv(Form("baseline_%i.csv", ipt), hists);
    }
  }

}

void nuclei_hists()
{
  int k_pick = 2;
  int etabin_pick = 2;
  int ptbin_pick = 2;
  const int nspecies_picks = 6;
  static int species_picks[nspecies_picks] = {1, 4, 5, 6, 7, 8};

  // using ep K=0 w/ 2E8 events the baseline

  // with R_L on the x-axis, plotting (alpha_i * K=i - K=0) / (int R_L K=0)
  /*
  mclogx(cno++);
  {
    float plot_xrange_lo = 5E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.08; //0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    // proton
    temp = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();
    temp->Add(temp, -1);
    temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
    temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
    temp->GetXaxis()->SetTitle("R_{L}");
    temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
    temp->SetMarkerColor(pt_color[0]);
    temp->SetLineColor(pt_color[0]);
    temp->SetMarkerSize(0.5);
    temp->SetMarkerStyle(21);
    temp->Draw("same hist");
    leg->AddEntry(temp,"e+p, K = 0");

    for (int i = 0; i < nspecies_picks; i++)
    {
      int ispecies = species_picks[i];
      temp = (TH1D*) h1d_jet_eec[ispecies][k_pick][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[i+1]);
      temp->SetLineColor(pt_color[i+1]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("%s, K = %i",species[ispecies], k[k_pick]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+A @ 10+100 GeV, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_by_nuclei_ratio.pdf\")", cno-1, out_dir) );

  }
  */

  // with R_L*sqrt(pt) on the x-axis, plotting (alpha_i * K=i - K=0) / (int R_L K=0)
  mclogx(cno++);
  {
    vector<TH1*> hists;

    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.08; //0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    // proton
    temp = (TH1D*) h1d_jet_eec_rlsqrtpt[0][0][etabin_pick][ptbin_pick]->Clone();
    temp_baseline = (TH1D*) h1d_jet_eec_rlsqrtpt[0][0][etabin_pick][ptbin_pick]->Clone();
    temp->Add(temp, -1);
    temp->Scale(1/temp_baseline->Integral());
    temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
    temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
    temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
    temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
    temp->SetMarkerColor(pt_color[0]);
    temp->SetLineColor(pt_color[0]);
    temp->SetMarkerSize(0.5);
    temp->SetMarkerStyle(21);
    temp->Draw("same hist");
    leg->AddEntry(temp,"e+p, K = 0");
    hists.push_back(temp);

    for (int i = 0; i < nspecies_picks; i++)
    {
      int ispecies = species_picks[i];
      temp = (TH1D*) h1d_jet_eec_rlsqrtpt[ispecies][k_pick][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_rlsqrtpt[0][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rlsqrtpt_norm_lo);
      int norm_binrange_hi = temp->FindBin(rlsqrtpt_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      hists.push_back(temp);

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[i+1]);
      temp->SetLineColor(pt_color[i+1]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("%s, K = %i",species[ispecies], k[k_pick]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+A @ 10+100 GeV, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_by_nuclei_ratio.pdf\")", cno-1, out_dir) );

    hists_to_csv("fig5.csv", hists);
  }

  /*
  // with R_L*sqrt(pt) on the x-axis, plotting (alpha_i * K=i) / (int R_L K=0)
  mclogxy(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = 1E-3;
    float plot_yrange_hi = 5E-1;
    float legend_x = 0.7;
    float legend_y = 0.2;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ispecies = 1; ispecies < speciesnum; ispecies++)
    {
      temp = (TH1D*) h1d_jet_eec_rlsqrtpt[ispecies][k_pick][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_rlsqrtpt[0][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rlsqrtpt_norm_lo);
      int norm_binrange_hi = temp->FindBin(rlsqrtpt_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Scale(1/temp_baseline->Integral());

      // plot histogram
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on)");
      temp->SetMarkerColor(pt_color[ispecies]);
      temp->SetLineColor(pt_color[ispecies]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("%s, K = %i",species[ispecies], k[k_pick]));
    }

    leg->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+A @ 10+100 GeV, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_by_nuclei.pdf\")", cno-1, out_dir) );
  }
  */

}


void power_hists()
{
  int ptbin_pick = 2;
  int etabin_pick = 2;

  // with R_L on the x-axis, plotting (alpha_i * K=i - K=0) / (int R_L K=0)
  mclogx(cno++);
  {
    vector<TH1*> hists;

    float plot_xrange_lo = 5E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.025; //-0.015;
    float plot_yrange_hi = 0.08; // 0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ipower = 0; ipower < pownum; ipower++)
    {
      temp = (TH1D*) h1d_jet_eec_eAu_by_power[ipower][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_ep_by_power[ipower][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      hists.push_back(temp);

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[ipower]);
      temp->SetLineColor(pt_color[ipower]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("power = %.1f", power[ipower]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+Au, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_by_pow_ratio.pdf\")", cno-1, out_dir) );

    // hists_to_csv("fig4.csv", hists);
  }

  // with R_L*sqrt(pt) on the x-axis, plotting (alpha_i * K=i - K=0) / (int R_L K=0)
  /*
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-1;
    float plot_xrange_hi = 5;
    float plot_yrange_lo = -0.025; //-0.015;
    float plot_yrange_hi = 0.08; // 0.04;
    float legend_x = 0.22;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ipower = 0; ipower < pownum; ipower++)
    {
      temp = (TH1D*) h1d_jet_eec_rlsqrtpt_eAu_by_power[ipower][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_rlsqrtpt_ep_by_power[ipower][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rlsqrtpt_norm_lo);
      int norm_binrange_hi = temp->FindBin(rlsqrtpt_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T,jet}}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[ipower]);
      temp->SetLineColor(pt_color[ipower]);
      temp->SetMarkerSize(0.5);
      temp->SetMarkerStyle(21);
      temp->Draw("same hist");
      leg->AddEntry(temp,Form("power = %.1f", power[ipower]));
    }
    leg->Draw("same");

    TLine l1(plot_xrange_lo,0,plot_xrange_hi,0);
    l1.SetLineStyle(7);
    l1.SetLineColor(kGray+2);
    l1.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.028);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+Au, 4*10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_rlsqrtpt_by_pow_ratio.pdf\")", cno-1, out_dir) );

  }
  */
}

void pt_spectra()
{
  // overlay of jet pt distribution for each species, K=4 10x100
  mclogy(cno++);
  {
    float legend_x = 0.6;
    float legend_y = 0.6;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.20);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;

    for (int ispecies = 1; ispecies < speciesnum; ispecies++)
    {
      temp = (TH1D*) h1d_jet_pt[ispecies]->Clone();
      temp->GetXaxis()->SetRangeUser(0,50);
      temp->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      temp->GetYaxis()->SetTitle("counts");
      temp->SetMarkerColor(pt_color[ispecies-1]);
      temp->SetLineColor(pt_color[ispecies-1]);
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->Draw("same hist e");
      leg->AddEntry(temp,species[ispecies]);
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt_spectra.pdf\")", cno-1, out_dir) );
  }
}

void peak_height_vs_A()
{
  int k_pick = 2;
  int ptbin_pick = 2;
  int etabin_pick = 2;

  // y-value of (alpha_i * K=i - K=0) / (int R_L K=0) at R_L = 1 vs A^{1/6} of nucleus. Each point is an eA nuclei species
  mcs(cno++);
  {
    float plot_xrange_lo = 0.8;
    float plot_xrange_hi = 2.7;
    float plot_yrange_lo = -0.005;
    float plot_yrange_hi = 0.04;

    // get x values
    static double species_logA[speciesnum] = {};
    for (int ispecies = 0; ispecies < speciesnum; ispecies++)
    {
      species_logA[ispecies] = TMath::Log10(species_A[ispecies]);
    }

    // get y values
    double peak_height_by_A[speciesnum] = {};
    peak_height_by_A[0] = 0;

    TH1D* temp;
    TH1D* temp_baseline;
    for (int ispecies = 1; ispecies < speciesnum; ispecies++)
    {
      temp = (TH1D*) h1d_jet_eec[ispecies][k_pick][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec[0][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      peak_height_by_A[ispecies] = temp->GetBinContent(temp->FindBin(0.999));
    }

    auto g = new TGraph(speciesnum, species_logA, peak_height_by_A);

    //g->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
    //g->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
    g->GetXaxis()->SetTitle("log(A)");
    g->GetYaxis()->SetTitle("normalized EEC value @ R_{L}=1");
    g->SetMarkerColor(pt_color[0]);
    g->SetLineColor(pt_color[0]);
    g->SetMarkerSize(0.5);
    g->SetMarkerStyle(21);
    g->Draw("");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_height_by_logA.pdf\")", cno-1, out_dir) );

    graph_to_csv("fig5_insert.csv", g);
  }

}

void peak_height_vs_A_isospin()
{
  int ptbin_pick = 2;
  int etabin_pick = 2;

  // nPDFset = 0 for all these data sets
  // y-value of (alpha_i * K=i - K=0) / (int R_L K=0) at R_L = 1 vs A^{1/6} of nucleus. Each point is an eA nuclei species
  mcs(cno++);
  {
    float plot_xrange_lo = 0.8;
    float plot_xrange_hi = 2.7;
    float plot_yrange_lo = -0.005;
    float plot_yrange_hi = 0.04;

    // get x values
    static double species_logA[speciesnum] = {};
    for (int ispecies = 0; ispecies < speciesnum; ispecies++)
    {
      species_logA[ispecies] = TMath::Log10(species_A[ispecies]);
    }

    // get y values
    double peak_height_by_A[speciesnum] = {};
    peak_height_by_A[0] = 0;

    TH1D* temp;
    TH1D* temp_baseline;
    for (int ispecies = 1; ispecies < speciesnum; ispecies++)
    {
      temp = (TH1D*) h1d_jet_eec_isospin[ispecies][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec_isospin[0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(rl_norm_lo);
      int norm_binrange_hi = temp->FindBin(rl_norm_hi);
      if (norm_binrange_lo == 0)
      {
        norm_binrange_lo = 1;
        cout<<"bin range lo too low; set to 1"<<endl;
      }
      if (norm_binrange_hi > temp->GetNbinsX())
      {
        norm_binrange_lo = temp->GetNbinsX();
        cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
      }
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      peak_height_by_A[ispecies] = temp->GetBinContent(temp->FindBin(0.999));
    }

    auto g = new TGraph(speciesnum, species_logA, peak_height_by_A);

    //g->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
    //g->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
    g->GetXaxis()->SetTitle("log(A)");
    g->GetYaxis()->SetTitle("normalized EEC value @ R_{L}=1");
    g->SetMarkerColor(pt_color[0]);
    g->SetLineColor(pt_color[0]);
    g->SetMarkerSize(0.5);
    g->SetMarkerStyle(21);
    g->Draw("");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_height_by_logA_isospin.pdf\")", cno-1, out_dir) );

    graph_to_csv("fig5_insert_isospin.csv", g);
  }

}

void Q2_x_panel()
{
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      mclogxy(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
      {
        float plot_xrange_lo = 0.8;
        float plot_xrange_hi = 2.7;
        float plot_yrange_lo = -0.005;
        float plot_yrange_hi = 0.04;
        float plot_zrange_lo = 0;
        float plot_zrange_hi = 750E3;

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

        if (k[ik] == 4)
        {
          h1d_jet_pt[ispecies] = (TH1D*) fin->Get("h1d_jet_pt");
          h1d_jet_pt[ispecies]->SetName(Form("h1d_jet_pt_%d", ispecies));
        }

        for (int ieta = 0; ieta < etabin; ieta++)
        {
          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            // raw data histograms
            h1d_jet_eec[ispecies][ik][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
            h1d_jet_eec[ispecies][ik][ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d_%d_%d", ieta, ipt, ispecies, ik));

            h1d_jet_eec_rlsqrtpt[ispecies][ik][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));
            h1d_jet_eec_rlsqrtpt[ispecies][ik][ieta][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_%d_%d_%d_%d", ieta, ipt, ispecies, ik));
          }
        }
        cout<<fin_name<<" loaded!"<<endl;
      }
      else
      {
        cout<<"couldn't find file for "<<species[ispecies]<<" with K = "<<k[ik]<<endl;
      }

    }
  }

  for (int ipower = 0; ipower < pownum; ipower++)
  {
    // e+p
    fin_name = fname_ep_by_power[ipower];
    fin = new TFile(fin_name, "READ");

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        // raw data histograms
        h1d_jet_eec_ep_by_power[ipower][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
        h1d_jet_eec_ep_by_power[ipower][ieta][ipt]->SetName(Form("h1d_jet_eec_eAubypow_ep_%d_%d_%d", ieta, ipt, ipower));

        h1d_jet_eec_rlsqrtpt_ep_by_power[ipower][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));
        h1d_jet_eec_rlsqrtpt_ep_by_power[ipower][ieta][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_eAubypow_ep_%d_%d_%d", ieta, ipt, ipower));
      }
    }

    cout<<fin_name<<" loaded!"<<endl;

    // e+Au
    fin_name = fname_eAu_by_power[ipower];
    fin = new TFile(fin_name, "READ");

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        // raw data histograms
        h1d_jet_eec_eAu_by_power[ipower][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
        h1d_jet_eec_eAu_by_power[ipower][ieta][ipt]->SetName(Form("h1d_jet_eec_eAubypow_eAu_%d_%d_%d", ieta, ipt, ipower));

        h1d_jet_eec_rlsqrtpt_eAu_by_power[ipower][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt));
        h1d_jet_eec_rlsqrtpt_eAu_by_power[ipower][ieta][ipt]->SetName(Form("h1d_jet_eec_rlsqrtpt_eAubypow_eAu_%d_%d_%d", ieta, ipt, ipower));
      }
    }

    cout<<fin_name<<" loaded!"<<endl;
  }

  for (int ispecies = 0; ispecies < speciesnum; ispecies++)
  {
    fin_name = fname_eA_isospin[ispecies];
    fin = new TFile(fin_name, "READ");

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        // raw data histograms
        h1d_jet_eec_isospin[ispecies][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
        h1d_jet_eec_isospin[ispecies][ieta][ipt]->SetName(Form("h1d_jet_eec_isospin_%d_%d_%d", ieta, ipt, ispecies));
      }
    }

    cout<<fin_name<<" loaded!"<<endl;
  }

  fin_name = fname_ep_Q2_x;
  fin = new TFile(fin_name, "READ");

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      // raw data histograms
      h2d_jet_Q2_x[ieta][ipt] = (TH2D*) fin->Get(Form("h2d_Q2_x_%d_%d", ieta, ipt));
      h2d_jet_Q2_x[ieta][ipt]->SetName(Form("h2d_Q2_x_%d_%d", ieta, ipt));
    }
  }

  cout<<fin_name<<" loaded!"<<endl;

  // plot individual panels

  pt_eta_3by3_hists();

  pt_bin_side_by_side();

  baseline_comparison();

  nuclei_hists();

  power_hists();

  //pt_spectra();

  peak_height_vs_A();

  peak_height_vs_A_isospin();

  Q2_x_panel();

}
