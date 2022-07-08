#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

/*
const int sys_bins = 5;
const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C"};
const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC"};
const int sys_color[sys_bins] = {kBlack, kRed, kOrange+1, kGreen+1, kBlue};
*/

const int sys_bins = 4;
const char* sys_name[sys_bins] = {"e+C", "e+Cu", "e+Au", "e+Pb"};
const char* sys_abbr[sys_bins] = {"eC", "eCu", "eAu", "ePb"};
const int sys_color[sys_bins] = {kBlack, kRed, kOrange+1, kGreen+1};

const int energy_bins = 5;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275"};

TH1D* h1d_D0_z_in_eta_gen_ratio[sys_bins][Q2bin][xbin][etabin] = {0};
TH1D* h1d_Lc_z_in_eta_gen_ratio[sys_bins][Q2bin][xbin][etabin] = {0};
TH1D* h1d_pion_z_in_eta_gen_ratio[sys_bins][Q2bin][xbin][etabin] = {0};

static int cno = 0;

void plot_comparison(const int energy_eA_option = 3)
{
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    if (iQ2<Q2bin-1) continue;

    for (int ix = 0; ix < xbin; ++ix)
    {
      if (ix<xbin-1) continue;

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (ieta<etabin-1) continue;

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 0.;
          float plot_yrange_hi = 2;

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{D^{0}}");
          htemp.GetYaxis()->SetTitle("R^{D^{0}}_{eA}");
          myhset(&htemp,1.2,1.6,0.05,0.05);

          TLegend leg(0.55,0.71,0.84,0.87);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          for (int isys = 1; isys < sys_bins; ++isys)
          {
            h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetMarkerColor(sys_color[isys]);
            h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetLineColor(sys_color[isys]);
            h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->Draw("same");
            leg.AddEntry(h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta],Form("%s @ %s",sys_name[isys],energy_name[energy_eA_option]),"p");
          }

          leg.Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.19,0.17,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 0.;
          float plot_yrange_hi = 2;

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
          htemp.GetYaxis()->SetTitle("R^{#Lambda_{c}}_{eA}");
          myhset(&htemp,1.2,1.6,0.05,0.05);

          TLegend leg(0.55,0.71,0.84,0.87);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          for (int isys = 1; isys < sys_bins; ++isys)
          {
            h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetMarkerColor(sys_color[isys]);
            h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetLineColor(sys_color[isys]);
            h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->Draw("same");
            leg.AddEntry(h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta],Form("%s @ %s",sys_name[isys],energy_name[energy_eA_option]),"p");
          }

          leg.Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.19,0.17,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 0.5;
          float plot_yrange_hi = 1.5;

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
          htemp.GetYaxis()->SetTitle("R^{#pi^{#pm}}_{eA}");
          myhset(&htemp,1.2,1.6,0.05,0.05);

          TLegend leg(0.55,0.71,0.84,0.87);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          for (int isys = 1; isys < sys_bins; ++isys)
          {
            h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetMarkerColor(sys_color[isys]);
            h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetLineColor(sys_color[isys]);
            h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->Draw("same");
            leg.AddEntry(h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta],Form("%s @ %s",sys_name[isys],energy_name[energy_eA_option]),"p");
          }

          leg.Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.19,0.17,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }
      }
    }
  }
}

void compare_chadron_gen_diff_sys(const int energy_eA_option = 3)
{
  mcs(-1);

  TFile* fin[sys_bins] = {0};
  for (int isys = 1; isys < sys_bins; ++isys)
  { // skipping e+p system
    fin[isys] = new TFile(Form("hists_gen_%s.root",sys_abbr[isys]),"READ");
    for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
    {
      for (int ix = 0; ix < xbin; ++ix)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta] = (TH1D*)fin[isys]->Get( Form("h1d_hadron_421_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
          h1d_D0_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetName( Form("h1d_hadron_421_z_in_eta_gen_ratio_%s_Q2%d_x%d_eta%d",sys_abbr[isys],iQ2,ix,ieta) );

          h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta] = (TH1D*)fin[isys]->Get( Form("h1d_hadron_4122_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
          h1d_Lc_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetName( Form("h1d_hadron_4122_z_in_eta_gen_ratio_%s_Q2%d_x%d_eta%d",sys_abbr[isys],iQ2,ix,ieta) );

          h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta] = (TH1D*)fin[isys]->Get( Form("h1d_hadron_211_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
          h1d_pion_z_in_eta_gen_ratio[isys][iQ2][ix][ieta]->SetName( Form("h1d_hadron_211_z_in_eta_gen_ratio_%s_Q2%d_x%d_eta%d",sys_abbr[isys],iQ2,ix,ieta) );
        }
      }
    }

    cout<<"system #"<<isys<<" works"<<endl;
  }

  plot_comparison(energy_eA_option);
}
