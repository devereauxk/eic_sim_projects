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

const char* out_dir = "./paperplots/";

TH1D* h1d_jet_eec[speciesnum][knum][etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt[speciesnum][knum][ptbin] = {};

static int cno = 0;

void pt_eta_3by3_hists()
{
  // with R_L on the x-axis, plotting (alpha_i * K=i) / (int dR_L K=0)
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      mclogxy(cno++);
      {
        float plot_xrange_lo = 0.02;
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

        for (int ik = 0; ik < knum; ik++)
        {
          temp = (TH1D*) h1d_jet_eec[2][ik][ieta][ipt]->Clone();
          temp_baseline = (TH1D*) h1d_jet_eec[2][0][ieta][ipt]->Clone();

          // calculate relative normalization ratio
          int norm_binrange_lo = temp->FindBin(1E-3);
          int norm_binrange_hi = temp->FindBin(1E-1);
          double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
          temp->Scale(relative_normalization);
          temp->Scale(1/temp_baseline->Integral());

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
          leg->AddEntry(temp,Form("K = %i",k[ik]));
        }

        leg->Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.028);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+Au @ 10+100 GeV, 10^{8} events");
        tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay_%d_%d.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

void nuclei_hists()
{
  // with R_L on the x-axis, plotting (alpha_i * K=i - K=0) / (int R_L K=0)
  int k_pick = 2;
  int ptbin_pick = 2;
  int etabin_pick = 2;

  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-2;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = -0.015;
    float plot_yrange_hi = 0.04;
    float legend_x = 0.7;
    float legend_y = 0.2;

    TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.028);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;
    TH1D* temp_baseline;

    for (int ispecies = 0; ispecies < speciesnum; ispecies++)
    {
      temp = (TH1D*) h1d_jet_eec[ispecies][k_pick][etabin_pick][ptbin_pick]->Clone();
      temp_baseline = (TH1D*) h1d_jet_eec[ispecies][0][etabin_pick][ptbin_pick]->Clone();

      // calculate relative normalization ratio
      int norm_binrange_lo = temp->FindBin(1E-3);
      int norm_binrange_hi = temp->FindBin(1E-1);
      double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
      temp->Scale(relative_normalization);
      temp->Add(temp_baseline, -1);
      temp->Scale(1/temp_baseline->Integral());

      // plot
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("R_{L}");
      temp->GetYaxis()->SetTitle("normalized EEC (rel. norm. * on - off)");
      temp->SetMarkerColor(pt_color[ispecies]);
      temp->SetLineColor(pt_color[ispecies]);
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
    tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+A @ 10+100 GeV, 10^{8} events");
    tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick]));
    tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ptbin_pick],pt_hi[ptbin_pick]));

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_by_nuclei.pdf\")", cno-1, out_dir) );

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
        cout<<fin_name<<" loaded!"<<endl;
      }
      else
      {
        cout<<"couldn't find file for "<<species[ispecies]<<" with K = "<<k[ik]<<endl;
      }

    }
  }

  // plot individual panels

  pt_eta_3by3_hists();

  nuclei_hists();

  //energy_hists();

  //power_hists();

}
