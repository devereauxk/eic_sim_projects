#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 7;
const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+d"};
const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "ed"};

const int energy_bins = 5;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275"};

TH2D* h2d_D0_pt_vs_eta_gen_ep[Q2bin][xbin] = {0};
TH2D* h2d_D0_z_vs_eta_gen_ep[Q2bin][xbin] = {0};
TH2D* h2d_Lc_pt_vs_eta_gen_ep[Q2bin][xbin] = {0};
TH2D* h2d_Lc_z_vs_eta_gen_ep[Q2bin][xbin] = {0};
TH2D* h2d_pion_pt_vs_eta_gen_ep[Q2bin][xbin] = {0};
TH2D* h2d_pion_z_vs_eta_gen_ep[Q2bin][xbin] = {0};

TH1D* h1d_nevt_eA[Q2bin][xbin] = {0};
TH2D* h2d_D0_pt_vs_eta_gen_eA[Q2bin][xbin] = {0};
TH2D* h2d_D0_z_vs_eta_gen_eA[Q2bin][xbin] = {0};
TH2D* h2d_Lc_pt_vs_eta_gen_eA[Q2bin][xbin] = {0};
TH2D* h2d_Lc_z_vs_eta_gen_eA[Q2bin][xbin] = {0};
TH2D* h2d_pion_pt_vs_eta_gen_eA[Q2bin][xbin] = {0};
TH2D* h2d_pion_z_vs_eta_gen_eA[Q2bin][xbin] = {0};

TH1D* h1d_nevt_ep[Q2bin][xbin] = {0};
TH1D* h1d_D0_z_in_eta_gen_ep[Q2bin][xbin][etabin] = {0};
TH1D* h1d_Lc_z_in_eta_gen_ep[Q2bin][xbin][etabin] = {0};
TH1D* h1d_pion_z_in_eta_gen_ep[Q2bin][xbin][etabin] = {0};
TH1D* h1d_D0_z_in_eta_gen_eA[Q2bin][xbin][etabin] = {0};
TH1D* h1d_Lc_z_in_eta_gen_eA[Q2bin][xbin][etabin] = {0};
TH1D* h1d_pion_z_in_eta_gen_eA[Q2bin][xbin][etabin] = {0};

TH1D* h1d_D0_z_in_eta_gen_ratio[Q2bin][xbin][etabin] = {0};
TH1D* h1d_Lc_z_in_eta_gen_ratio[Q2bin][xbin][etabin] = {0};
TH1D* h1d_pion_z_in_eta_gen_ratio[Q2bin][xbin][etabin] = {0};

static int cno = 0;

void slide_in_eta(const int sys_ep_option, const int energy_ep_option, const int sys_eA_option, const int energy_eA_option)
{ // make sure slice_lo/hi are double before slicing
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta] = (TH1D*)h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->ProjectionX( Form("h1d_D0_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix,ieta) );

        h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta] = (TH1D*)h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->ProjectionX( Form("h1d_Lc_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix,ieta) );

        h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta] = (TH1D*)h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->ProjectionX( Form("h1d_pion_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix,ieta) );

        h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta] = (TH1D*)h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->ProjectionX( Form("h1d_D0_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix,ieta) );

        h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta] = (TH1D*)h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->ProjectionX( Form("h1d_Lc_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix,ieta) );

        h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
        h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta] = (TH1D*)h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->ProjectionX( Form("h1d_pion_z_in_eta_gen_%s_%s_Q2%d_x%d_eta%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix,ieta) );

        h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->Rebin(4);
        h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->Rebin(4);
        h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->Rebin(4);
        h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->Rebin(4);
        h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->Rebin(4);
        h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->Rebin(4);
      }

      // set back to the original range
      h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
      h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
      h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
      h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
      h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
      h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->GetYaxis()->SetRangeUser(-4,4);
    }
  }
}

void set_ratio_hists()
{ // NB: may add names later
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta] = (TH1D*)h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->Clone();
        h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetName( Form("h1d_D0_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
        h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Divide(h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]);

        h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta] = (TH1D*)h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->Clone();
        h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetName( Form("h1d_Lc_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
        h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Divide(h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]);

        h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta] = (TH1D*)h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->Clone();
        h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetName( Form("h1d_pion_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",iQ2,ix,ieta) );
        h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Divide(h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]);
      }
    }
  }
}

void plot_1D(const int sys_ep_option, const int energy_ep_option, const int sys_eA_option, const int energy_eA_option)
{
  TGaxis::SetMaxDigits(3);

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    // if (Q2_hi[iQ2]<=1E1) continue;
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

          float plot_yrange_lo = 1E-6;
          float plot_yrange_hi = 1.5*h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->GetMaximum();

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{D^{0}}");
          htemp.GetYaxis()->SetTitle("N^{D^{0}}/N^{e}");
          myhset(&htemp,1.2,1.6,0.05,0.045);

          TLegend leg(0.2,0.72,0.83,0.82);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->Draw("same");
          // leg.AddEntry(h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta],"Pythia e+p 10+100GeV","lp");
          leg.AddEntry(h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]),"lp");

          h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->SetLineColor(kRed);
          h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerColor(kRed);
          h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->Draw("same");
          // leg.AddEntry(h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta],"BeAGLE e+Au 10+110GeV","lp");
          leg.AddEntry(h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]),"lp");

          leg.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.85,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 1E-6;
          float plot_yrange_hi = 1.5*h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->GetMaximum();

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
          htemp.GetYaxis()->SetTitle("N^{#Lambda_{c}}/N^{e}");
          myhset(&htemp,1.2,1.6,0.05,0.045);

          TLegend leg(0.2,0.72,0.83,0.82);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]),"lp");

          h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->SetLineColor(kRed);
          h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerColor(kRed);
          h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]),"lp");

          leg.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.85,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 1E-6;
          float plot_yrange_hi = 1.5*h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->GetMaximum();

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
          htemp.GetYaxis()->SetTitle("N^{#pi^{#pm}}/N^{e}");
          myhset(&htemp,1.2,1.6,0.05,0.045);

          TLegend leg(0.2,0.72,0.83,0.82);
          leg.SetBorderSize(0);
          leg.SetTextSize(0.03);
          leg.SetFillStyle(0);
          leg.SetMargin(0.1);

          h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]),"lp");

          h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->SetLineColor(kRed);
          h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->SetMarkerColor(kRed);
          h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta],Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]),"lp");

          leg.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.85,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }
      }
    }
  }
}

void plot_ratio(const int sys_ep_option, const int energy_ep_option, const int sys_eA_option, const int energy_eA_option)
{
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    // if (Q2_hi[iQ2]<=1E1) continue;
    if (iQ2!=(Q2bin-1)) continue;

    for (int ix = 0; ix < xbin; ++ix)
    {
      if (ix!=(xbin-1)) continue;

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
          // htemp.GetYaxis()->SetTitle("(N^{D^{0}}/N^{e} in e+p)/(N^{D^{0}}/N^{e} in e+Au)");
          myhset(&htemp,1.2,1.6,0.05,0.05);

          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.84,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
          tl->DrawLatexNDC(0.23,0.80,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
          tl->DrawLatexNDC(0.23,0.75,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

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
          // htemp.GetYaxis()->SetTitle("(N^{#Lambda_{c}}/N^{e} in e+p)/(N^{#Lambda_{c}}/N^{e} in e+Au)");
          myhset(&htemp,1.2,1.6,0.05,0.045);

          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.84,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
          tl->DrawLatexNDC(0.23,0.80,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
          tl->DrawLatexNDC(0.23,0.75,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }

        mcs(cno++);
        {
          float plot_xrange_lo = 0;
          float plot_xrange_hi = 1;

          float plot_yrange_lo = 0.;
          float plot_yrange_hi = 2;

          TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp.Draw();
          htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
          htemp.GetYaxis()->SetTitle("R^{#pi^{#pm}}_{eA}");
          // htemp.GetYaxis()->SetTitle("(N^{#pi^{#pm}}/N^{e} in e+p)/(N^{#pi^{#pm}}/N^{e} in e+Au)");
          myhset(&htemp,1.2,1.6,0.05,0.045);

          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.7);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");

          TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
          l1.SetLineStyle(7);
          l1.SetLineColor(kGray+2);
          l1.Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.23,0.84,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
          tl->DrawLatexNDC(0.23,0.80,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
          tl->DrawLatexNDC(0.23,0.75,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, iQ2, ix, ieta) );
        }
      }
    }
  }
}

void plot_ratio_comparison(const int sys_ep_option, const int energy_ep_option, const int sys_eA_option, const int energy_eA_option)
{
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
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
        // htemp.GetYaxis()->SetTitle("(N^{D^{0}}/N^{e} in e+p)/(N^{D^{0}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ix = 0; ix < xbin-1; ++ix)
        { // excluding inclusive bin
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(x_color[ix]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(x_color[ix]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < x < %.1e",x_lo[ix],x_hi[ix]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_ratio_z_in_eta_allx_%d_%d.pdf\")", cno-1, iQ2, ieta) );
      }
    }
  }

  for (int ix = 0; ix < xbin; ++ix)
  {
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
        // htemp.GetYaxis()->SetTitle("(N^{D^{0}}/N^{e} in e+p)/(N^{D^{0}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        { // excluding inclusive bin
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(Q2_color[iQ2]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(Q2_color[iQ2]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < Q^{2} < %.1e",Q2_lo[iQ2],Q2_hi[iQ2]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_ratio_z_in_eta_allQ2_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }
/*
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    if (iQ2<Q2bin-1) continue;

    for (int ix = 0; ix < xbin; ++ix)
    {
      if (ix<xbin-1) continue;

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
        // htemp.GetYaxis()->SetTitle("(N^{D^{0}}/N^{e} in e+p)/(N^{D^{0}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ieta = 0; ieta < etabin-1; ++ieta)
        { // excluding inclusive bin
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(eta_color[ieta]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(eta_color[ieta]);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/D0_ratio_z_in_eta_alleta_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }
*/

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
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
        htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
        htemp.GetYaxis()->SetTitle("R^{#Lambda_{c}}_{eA}");
        // htemp.GetYaxis()->SetTitle("(N^{#Lambda_{c}}/N^{e} in e+p)/(N^{#Lambda_{c}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ix = 0; ix < xbin-1; ++ix)
        { // excluding inclusive bin
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(x_color[ix]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(x_color[ix]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < x < %.1e",x_lo[ix],x_hi[ix]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_ratio_z_in_eta_allx_%d_%d.pdf\")", cno-1, iQ2, ieta) );
      }
    }
  }

  for (int ix = 0; ix < xbin; ++ix)
  {
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
        htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
        htemp.GetYaxis()->SetTitle("R^{#Lambda_{c}}_{eA}");
        // htemp.GetYaxis()->SetTitle("(N^{#Lambda_{c}}/N^{e} in e+p)/(N^{#Lambda_{c}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        { // excluding inclusive bin
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(Q2_color[iQ2]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(Q2_color[iQ2]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < Q^{2} < %.1e",Q2_lo[iQ2],Q2_hi[iQ2]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_ratio_z_in_eta_allQ2_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }

/*
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    if (iQ2<Q2bin-1) continue;

    for (int ix = 0; ix < xbin; ++ix)
    {
      if (ix<xbin-1) continue;

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
        // htemp.GetYaxis()->SetTitle("(N^{#Lambda_{c}}/N^{e} in e+p)/(N^{#Lambda_{c}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ieta = 0; ieta < etabin-1; ++ieta)
        { // excluding inclusive bin
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(eta_color[ieta]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(eta_color[ieta]);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/Lc_ratio_z_in_eta_alleta_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }*/

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
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
        htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
        htemp.GetYaxis()->SetTitle("R^{#pi^{#pm}}_{eA}");
        // htemp.GetYaxis()->SetTitle("(N^{#pi^{#pm}}/N^{e} in e+p)/(N^{#pi^{#pm}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ix = 0; ix < xbin-1; ++ix)
        { // excluding inclusive bin
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(x_color[ix]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(x_color[ix]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < x < %.1e",x_lo[ix],x_hi[ix]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_ratio_z_in_eta_allx_%d_%d.pdf\")", cno-1, iQ2, ieta) );
      }
    }
  }

  for (int ix = 0; ix < xbin; ++ix)
  {
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
        htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
        htemp.GetYaxis()->SetTitle("R^{#pi^{#pm}}_{eA}");
        // htemp.GetYaxis()->SetTitle("(N^{#pi^{#pm}}/N^{e} in e+p)/(N^{#pi^{#pm}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        { // excluding inclusive bin
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(Q2_color[iQ2]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(Q2_color[iQ2]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1e < Q^{2} < %.1e",Q2_lo[iQ2],Q2_hi[iQ2]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_ratio_z_in_eta_allQ2_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }
/*
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    if (iQ2<Q2bin-1) continue;

    for (int ix = 0; ix < xbin; ++ix)
    {
      if (ix<xbin-1) continue;

      mcs(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = 0.;
        float plot_yrange_hi = 2;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
        htemp.GetYaxis()->SetTitle("R^{#pi^{#pm}}_{eA}");
        // htemp.GetYaxis()->SetTitle("(N^{#pi^{#pm}}/N^{e} in e+p)/(N^{#pi^{#pm}}/N^{e} in e+Au)");
        myhset(&htemp,1.2,1.6,0.05,0.05);

        TLegend leg(0.55,0.73,0.84,0.85);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        for (int ieta = 0; ieta < etabin-1; ++ieta)
        { // excluding inclusive bin
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerColor(eta_color[ieta]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetLineColor(eta_color[ieta]);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerSize(0.5);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->SetMarkerStyle(21);
          h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Draw("same");
          leg.AddEntry(h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta],Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]),"p");
        }
        leg.Draw("same");

        TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
        l1.SetLineStyle(7);
        l1.SetLineColor(kGray+2);
        l1.Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.19,0.24,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
        tl->DrawLatexNDC(0.19,0.20,Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]));
        tl->DrawLatexNDC(0.19,0.15,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8, |#eta|<3.5");

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/pion_ratio_z_in_eta_alleta_%d_%d.pdf\")", cno-1, ix, ieta) );
      }
    }
  }
  */
}

void plot_2D(const int sys_ep_option, const int energy_ep_option, const int sys_eA_option, const int energy_eA_option)
{
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      // if (Q2_hi[iQ2]<=1E1) continue;
      // if (iQ2<Q2bin-1) continue;
      // if (ix<xbin-1) continue;

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{D^{0}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_ep_option],energy_name[energy_ep_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/ep_D0_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{D^{0}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_eA_option],energy_name[energy_eA_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/eA_D0_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_ep_option],energy_name[energy_ep_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/ep_Lc_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{#Lambda_{c}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_eA_option],energy_name[energy_eA_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/eA_Lc_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_ep_option],energy_name[energy_ep_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/ep_pion_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }

      mclogz(cno++);
      {
        float plot_xrange_lo = 0;
        float plot_xrange_hi = 1;

        float plot_yrange_lo = -4;
        float plot_yrange_hi = 4;

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("z^{#pi^{#pm}}");
        htemp.GetYaxis()->SetTitle("#eta");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.8);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->Draw("samecolz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.2,0.92,Form("%s @ %s, %.0f < Q^{2} < %.0f GeV^{2}, %.0e < x < %.0e ",sys_name[sys_eA_option],energy_name[energy_eA_option],Q2_lo[iQ2],Q2_hi[iQ2],x_lo[ix],x_hi[ix]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/eA_pion_z_vs_eta_%d_%d.pdf\")", cno-1, iQ2, ix) );
      }
    }
  }
}

void plot_chadron_gen(const char* inFile_ep = "hists_gen_ep.root", const int sys_ep_option = 0, const int energy_ep_option = 0, const char* inFile_eA = "hists_gen_eA.root", const int sys_eA_option = 0, const int energy_eA_option = 0, const char* outFile = "hists_gen.root")
{
  mcs(-1);

  TFile* fin_ep = new TFile(inFile_ep,"READ");
  // const int nevt_ep = h1d_nevt_ep->Integral();
  // const int nevt_ep = 1.1E8;
  // fin_ep->ls();
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_ep[iQ2][ix] = (TH1D*)fin_ep->Get(Form("h1d_nevt_%d_%d",iQ2,ix));
      h1d_nevt_ep[iQ2][ix]->SetName(Form("h1d_nevt_%s_%s_Q2%d_x%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix));
      int nevt_ep = h1d_nevt_ep[iQ2][ix]->Integral();
      cout << "ep # of events " << nevt_ep << " (Q2, x) bin (" << iQ2 << ", " << ix << ")" <<endl;

      // h2d_D0_z_vs_eta_gen_ep[iQ2][ix] = (TH2D*)fin_ep->Get(Form("h2d_D0_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_D0_z_vs_eta_gen_ep[iQ2][ix] = (TH2D*)fin_ep->Get(Form("h2d_hadron_421_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->SetName(Form("h2d_D0_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix));
      // h2d_Lc_z_vs_eta_gen_ep[iQ2][ix] = (TH2D*)fin_ep->Get(Form("h2d_Lc_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_Lc_z_vs_eta_gen_ep[iQ2][ix] = (TH2D*)fin_ep->Get(Form("h2d_hadron_4122_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->SetName(Form("h2d_Lc_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix));

      h2d_pion_z_vs_eta_gen_ep[iQ2][ix] = (TH2D*)fin_ep->Get(Form("h2d_hadron_211_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->SetName(Form("h2d_pion_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],iQ2,ix));

      if (nevt_ep>0)
      {
        h2d_D0_z_vs_eta_gen_ep[iQ2][ix]->Scale(1./nevt_ep);
        h2d_Lc_z_vs_eta_gen_ep[iQ2][ix]->Scale(1./nevt_ep);
        h2d_pion_z_vs_eta_gen_ep[iQ2][ix]->Scale(1./nevt_ep);
      }
    }
  }

  TFile* fin_eA = new TFile(inFile_eA,"READ");
  // h1d_nevt_eA = (TH1D*)fin_eA->Get("h1d_nevt");
  // const int nevt_eA = h1d_nevt_eA->Integral();
  // const int nevt_eA = 4.99E7 /5.; // 5x41 = 499*2E4
  // const int nevt_eA = 1E6; // BeAGLE ep 10x100GeV
  // const int nevt_eA = 4.96E7; // BeAGLE eAu 10x110 = 496*1E5

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_eA[iQ2][ix] = (TH1D*)fin_eA->Get(Form("h1d_nevt_%d_%d",iQ2,ix));
      h1d_nevt_eA[iQ2][ix]->SetName(Form("h1d_nevt_%s_%s_Q2%d_x%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix));
      int nevt_eA = h1d_nevt_eA[iQ2][ix]->Integral();
      cout << "eA # of events " << nevt_eA << " (Q2, x) bin (" << iQ2 << ", " << ix << ")" <<endl;

      // h2d_D0_z_vs_eta_gen_eA[iQ2][ix] = (TH2D*)fin_eA->Get(Form("h2d_D0_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_D0_z_vs_eta_gen_eA[iQ2][ix] = (TH2D*)fin_eA->Get(Form("h2d_hadron_421_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->SetName(Form("h2d_D0_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix));
      // h2d_Lc_z_vs_eta_gen_eA[iQ2][ix] = (TH2D*)fin_eA->Get(Form("h2d_Lc_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_Lc_z_vs_eta_gen_eA[iQ2][ix] = (TH2D*)fin_eA->Get(Form("h2d_hadron_4122_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->SetName(Form("h2d_Lc_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix));

      h2d_pion_z_vs_eta_gen_eA[iQ2][ix] = (TH2D*)fin_eA->Get(Form("h2d_hadron_211_z_vs_eta_gen_%d_%d",iQ2,ix));
      h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->SetName(Form("h2d_pion_z_vs_eta_gen_%s_%s_Q2%d_x%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],iQ2,ix));

      if (nevt_eA>0)
      {
        h2d_D0_z_vs_eta_gen_eA[iQ2][ix]->Scale(1./nevt_eA);
        h2d_Lc_z_vs_eta_gen_eA[iQ2][ix]->Scale(1./nevt_eA);
        h2d_pion_z_vs_eta_gen_eA[iQ2][ix]->Scale(1./nevt_eA);
      }
    }
  }

  slide_in_eta(sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option);

  plot_1D(sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option);

  plot_2D(sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option);

  set_ratio_hists();

  plot_ratio(sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option);
  plot_ratio_comparison(sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option);

  TFile* fout = new TFile(outFile,"recreate");
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        h1d_D0_z_in_eta_gen_ep[iQ2][ix][ieta]->Write();
        h1d_Lc_z_in_eta_gen_ep[iQ2][ix][ieta]->Write();
        h1d_pion_z_in_eta_gen_ep[iQ2][ix][ieta]->Write();
        h1d_D0_z_in_eta_gen_eA[iQ2][ix][ieta]->Write();
        h1d_Lc_z_in_eta_gen_eA[iQ2][ix][ieta]->Write();
        h1d_pion_z_in_eta_gen_eA[iQ2][ix][ieta]->Write();

        h1d_D0_z_in_eta_gen_ratio[iQ2][ix][ieta]->Write();
        h1d_Lc_z_in_eta_gen_ratio[iQ2][ix][ieta]->Write();
        h1d_pion_z_in_eta_gen_ratio[iQ2][ix][ieta]->Write();
      }
    }
  }
  fout->Write();

  // TFile* fout_proposal = new TFile("hists_gen_D0_pion.root","recreate");
  // h1d_D0_z_in_eta_gen_ep[Q2bin-1][xbin-1][0]->Write("h1d_D0_z_gen_ep");
  // h1d_D0_z_in_eta_gen_eA[Q2bin-1][xbin-1][0]->Write("h1d_D0_z_gen_eA");
  // h1d_pion_z_in_eta_gen_ep[Q2bin-1][xbin-1][0]->Write("h1d_pion_z_gen_ep");
  // h1d_pion_z_in_eta_gen_eA[Q2bin-1][xbin-1][0]->Write("h1d_pion_z_gen_eA");
  // h1d_D0_z_in_eta_gen_ratio[Q2bin-1][xbin-1][0]->Write("h1d_D0_z_gen_raio");
  // h1d_pion_z_in_eta_gen_ratio[Q2bin-1][xbin-1][0]->Write("h1d_pion_z_gen_raio");
  // fout_proposal->Write();
}
