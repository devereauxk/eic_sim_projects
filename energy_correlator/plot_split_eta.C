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
const int pt_color[9] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed, kCyan+1, kAzure+7, kViolet+7, kBlack};

const int etabin = 6; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5, -1, 0};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5, 0, 1};
const int eta_color[9] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed, kCyan+1, kAzure+7, kViolet+7, kBlack};

const int speciesnum = 2;
static char* species[speciesnum] = {(char*)"e+p", (char*)"e+Au"};
static int species_A[speciesnum] = {1, 197};
//static double species_A16[speciesnum] = {1, 1.12246, 1.51309, 1.84931, 2, 2.41219, 2.48941};

const int knum = 4;
static int k[knum] = {0,2,4,10};

// power=0.5
static char* fname_ep_by_K[knum] = {(char*)"./eHIJING/ep_10_100_K0_pow05/merged.root", (char*)"", (char*)"", (char*)""};
static char* fname_eAu_by_K[knum] = {(char*)"", (char*)"./eHIJING/eAu_1E8_K2_pow05/merged.root", (char*)"./eHIJING/eAu_1E8_K4_pow05/merged.root", (char*)"./eHIJING/eAu_1E8_condor_pow05/merged.root"};
static char** fname_eA_by_K[speciesnum] = {fname_ep_by_K, fname_eAu_by_K}; // K=4 cases are 2E8 events, 1E8 events otherwise

const char* out_dir = "./paperplots/test/";

TH1D* h1d_jet_eec[speciesnum][knum][etabin][ptbin] = {};

const float rl_norm_hi = 0.2; // 0.3 used for 3by3 top left and second col //0.08;
const float rl_norm_lo = 1E-3;

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

void pt_eta_3by3_hists()
{
  int species_pick = 1;

  // using ep K=0 w/ 2E8 events the baseline

  // 3x3 panel
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < 3; ipt++)
    {
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
        temp = (TH1D*) h1d_jet_eec[0][0][ieta][ipt]->Clone();
        temp_baseline = (TH1D*) h1d_jet_eec[0][0][ieta][ipt]->Clone();

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
          temp = (TH1D*) h1d_jet_eec[species_pick][ik][ieta][ipt]->Clone();
          temp_baseline = (TH1D*) h1d_jet_eec[0][0][ieta][ipt]->Clone();

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
        tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eec_overlay_%i_%i.pdf\")", cno-1, out_dir, ieta, ipt) );

        hists_to_csv(Form("raw_eec_%i_%i.csv", ieta, ipt), hists);
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

        for (int ieta = 0; ieta < etabin; ieta++)
        {
          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            // raw data histograms
            h1d_jet_eec[ispecies][ik][ieta][ipt] = (TH1D*) fin->Get(Form("h1d_jet_eec_%d_%d", ieta, ipt));
            h1d_jet_eec[ispecies][ik][ieta][ipt]->SetName(Form("h1d_jet_eec_%d_%d_%d_%d", ieta, ipt, ispecies, ik));
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

  // plot individual panels

  pt_eta_3by3_hists();


}
