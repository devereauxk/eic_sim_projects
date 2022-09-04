//#include "charm_H1.h"
#include "plot_events.h"
using namespace std;
unsigned int verbosity = 0;

PlotXsec* plot_Pythia6_ep = NULL;
//PlotXsec* plot_Pythia6_eAu = NULL;
//PlotXsec* plot_BeAGLE_ep = NULL;
//PlotXsec* plot_BeAGLE_eAu = NULL;

static int cno = 0;

/*
void set_charm_H1()
{
  double display_x_err[5] = {0}; // no more than 4 x bins, so 5 should be enough

  // Q2 = 12
  for (int ixbin = 0; ixbin < xbins_H1_12; ++ixbin)
  {
    rcs_staterr_cc_H1_12[ixbin] *= rcs_cc_H1_12[ixbin]*0.01;
    rcs_syserr_cc_H1_12[ixbin] *= rcs_cc_H1_12[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_12[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_12 = new TGraphErrors(xbins_H1_12, x_H1_12, rcs_cc_H1_12, 0, rcs_staterr_cc_H1_12);
  g_rcs_stat_vs_x_H1_12->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_12->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_12->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_12->SetLineColor(2);

  g_rcs_sys_vs_x_H1_12 = new TGraphErrors(xbins_H1_12, x_H1_12, rcs_cc_H1_12, display_x_err, rcs_syserr_cc_H1_12);
  g_rcs_sys_vs_x_H1_12->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_12->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_12->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_12->SetLineColor(2);
  g_rcs_sys_vs_x_H1_12->SetFillColorAlpha(2,0.3);

  // Q2 = 20
  for (int ixbin = 0; ixbin < xbins_H1_20; ++ixbin)
  {
    rcs_staterr_cc_H1_20[ixbin] *= rcs_cc_H1_20[ixbin]*0.01;
    rcs_syserr_cc_H1_20[ixbin] *= rcs_cc_H1_20[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_20[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_20 = new TGraphErrors(xbins_H1_20, x_H1_20, rcs_cc_H1_20, 0, rcs_staterr_cc_H1_20);
  g_rcs_stat_vs_x_H1_20->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_20->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_20->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_20->SetLineColor(2);

  g_rcs_sys_vs_x_H1_20 = new TGraphErrors(xbins_H1_20, x_H1_20, rcs_cc_H1_20, display_x_err, rcs_syserr_cc_H1_20);
  g_rcs_sys_vs_x_H1_20->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_20->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_20->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_20->SetLineColor(2);
  g_rcs_sys_vs_x_H1_20->SetFillColorAlpha(2,0.3);

  // Q2 = 35
  for (int ixbin = 0; ixbin < xbins_H1_35; ++ixbin)
  {
    rcs_staterr_cc_H1_35[ixbin] *= rcs_cc_H1_35[ixbin]*0.01;
    rcs_syserr_cc_H1_35[ixbin] *= rcs_cc_H1_35[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_35[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_35 = new TGraphErrors(xbins_H1_35, x_H1_35, rcs_cc_H1_35, 0, rcs_staterr_cc_H1_35);
  g_rcs_stat_vs_x_H1_35->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_35->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_35->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_35->SetLineColor(2);

  g_rcs_sys_vs_x_H1_35 = new TGraphErrors(xbins_H1_35, x_H1_35, rcs_cc_H1_35, display_x_err, rcs_syserr_cc_H1_35);
  g_rcs_sys_vs_x_H1_35->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_35->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_35->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_35->SetLineColor(2);
  g_rcs_sys_vs_x_H1_35->SetFillColorAlpha(2,0.3);

  // Q2 = 60
  for (int ixbin = 0; ixbin < xbins_H1_60; ++ixbin)
  {
    rcs_staterr_cc_H1_60[ixbin] *= rcs_cc_H1_60[ixbin]*0.01;
    rcs_syserr_cc_H1_60[ixbin] *= rcs_cc_H1_60[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_60[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_60 = new TGraphErrors(xbins_H1_60, x_H1_60, rcs_cc_H1_60, 0, rcs_staterr_cc_H1_60);
  g_rcs_stat_vs_x_H1_60->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_60->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_60->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_60->SetLineColor(2);

  g_rcs_sys_vs_x_H1_60 = new TGraphErrors(xbins_H1_60, x_H1_60, rcs_cc_H1_60, display_x_err, rcs_syserr_cc_H1_60);
  g_rcs_sys_vs_x_H1_60->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_60->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_60->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_60->SetLineColor(2);
  g_rcs_sys_vs_x_H1_60->SetFillColorAlpha(2,0.3);

  // Q2 = 120
  for (int ixbin = 0; ixbin < xbins_H1_120; ++ixbin)
  {
    rcs_staterr_cc_H1_120[ixbin] *= rcs_cc_H1_120[ixbin]*0.01;
    rcs_syserr_cc_H1_120[ixbin] *= rcs_cc_H1_120[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_120[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_120 = new TGraphErrors(xbins_H1_120, x_H1_120, rcs_cc_H1_120, 0, rcs_staterr_cc_H1_120);
  g_rcs_stat_vs_x_H1_120->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_120->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_120->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_120->SetLineColor(2);

  g_rcs_sys_vs_x_H1_120 = new TGraphErrors(xbins_H1_120, x_H1_120, rcs_cc_H1_120, display_x_err, rcs_syserr_cc_H1_120);
  g_rcs_sys_vs_x_H1_120->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_120->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_120->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_120->SetLineColor(2);
  g_rcs_sys_vs_x_H1_120->SetFillColorAlpha(2,0.3);

  // Q2 = 200
  for (int ixbin = 0; ixbin < xbins_H1_200; ++ixbin)
  {
    rcs_staterr_cc_H1_200[ixbin] *= rcs_cc_H1_200[ixbin]*0.01;
    rcs_syserr_cc_H1_200[ixbin] *= rcs_cc_H1_200[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_200[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_200 = new TGraphErrors(xbins_H1_200, x_H1_200, rcs_cc_H1_200, 0, rcs_staterr_cc_H1_200);
  g_rcs_stat_vs_x_H1_200->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_200->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_200->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_200->SetLineColor(2);

  g_rcs_sys_vs_x_H1_200 = new TGraphErrors(xbins_H1_200, x_H1_200, rcs_cc_H1_200, display_x_err, rcs_syserr_cc_H1_200);
  g_rcs_sys_vs_x_H1_200->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_200->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_200->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_200->SetLineColor(2);
  g_rcs_sys_vs_x_H1_200->SetFillColorAlpha(2,0.3);

  // Q2 = 300
  for (int ixbin = 0; ixbin < xbins_H1_300; ++ixbin)
  {
    rcs_staterr_cc_H1_300[ixbin] *= rcs_cc_H1_300[ixbin]*0.01;
    rcs_syserr_cc_H1_300[ixbin] *= rcs_cc_H1_300[ixbin]*0.01;
    display_x_err[ixbin] = x_H1_300[ixbin]*0.1;
  }
  g_rcs_stat_vs_x_H1_300 = new TGraphErrors(xbins_H1_300, x_H1_300, rcs_cc_H1_300, 0, rcs_staterr_cc_H1_300);
  g_rcs_stat_vs_x_H1_300->SetMarkerStyle(20);
  g_rcs_stat_vs_x_H1_300->SetMarkerSize(0.7);
  g_rcs_stat_vs_x_H1_300->SetMarkerColor(2);
  g_rcs_stat_vs_x_H1_300->SetLineColor(2);

  g_rcs_sys_vs_x_H1_300 = new TGraphErrors(xbins_H1_300, x_H1_300, rcs_cc_H1_300, display_x_err, rcs_syserr_cc_H1_300);
  g_rcs_sys_vs_x_H1_300->SetMarkerStyle(20);
  g_rcs_sys_vs_x_H1_300->SetMarkerSize(0.7);
  g_rcs_sys_vs_x_H1_300->SetMarkerColor(2);
  g_rcs_sys_vs_x_H1_300->SetLineColor(2);
  g_rcs_sys_vs_x_H1_300->SetFillColorAlpha(2,0.3);
}
*/

void set_Q2_binning()
{
  float size = 0.1;
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    Q2_lo[iQ2] = (1-size)*Q2_mid[iQ2];
    Q2_hi[iQ2] = (1+size)*Q2_mid[iQ2];
  }
}

void plot_graph(const char* outDir)
{
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-3;
    float plot_xrange_hi = 1E0;
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 0.5;
    for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
    {
      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("x");
      htemp.GetYaxis()->SetTitle("#tilde{#sigma}");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      TLegend leg(0.48,0.65,0.9,0.87);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.035);
      leg.SetFillStyle(0);

      plot_Pythia6_ep->g_cs_vs_x_c[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_Pythia6_ep->g_cs_vs_x_c[iQ2],Form("%s, #tilde{#sigma}_{c#bar{c}}",plot_Pythia6_ep->sys_latex),"p");

      plot_Pythia6_ep->g_cs_vs_x_s[iQ2]->SetLineColor(kBlue);
      plot_Pythia6_ep->g_cs_vs_x_s[iQ2]->SetMarkerColor(kBlue);
      plot_Pythia6_ep->g_cs_vs_x_s[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_Pythia6_ep->g_cs_vs_x_s[iQ2],Form("%s, #tilde{#sigma}_{s#bar{s}}",plot_Pythia6_ep->sys_latex),"p");

      /*
      plot_Pythia6_eAu->g_cs_vs_x_c[iQ2]->SetMarkerStyle(24);
      plot_Pythia6_eAu->g_cs_vs_x_c[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_Pythia6_eAu->g_cs_vs_x_c[iQ2],plot_Pythia6_eAu->sys_latex,"p");

      plot_BeAGLE_ep->g_cs_vs_x_c[iQ2]->SetLineColor(kBlue);
      plot_BeAGLE_ep->g_cs_vs_x_c[iQ2]->SetMarkerColor(kBlue);
      plot_BeAGLE_ep->g_cs_vs_x_c[iQ2]->SetMarkerStyle(21);
      plot_BeAGLE_ep->g_cs_vs_x_c[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_BeAGLE_ep->g_cs_vs_x_c[iQ2],plot_BeAGLE_ep->sys_latex,"p");

      plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2]->SetLineColor(kBlue);
      plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2]->SetMarkerColor(kBlue);
      plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2]->SetMarkerStyle(25);
      plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2],plot_BeAGLE_eAu->sys_latex,"p");
      */

      /*
      if (iQ2==0)
      {
        g_rcs_sys_vs_x_H1_12->Draw("samee2");
        g_rcs_stat_vs_x_H1_12->Draw("psame");
      }
      if (iQ2==1)
      {
        g_rcs_sys_vs_x_H1_20->Draw("samee2");
        g_rcs_stat_vs_x_H1_20->Draw("psame");
      }
      if (iQ2==2)
      {
        g_rcs_sys_vs_x_H1_35->Draw("samee2");
        g_rcs_stat_vs_x_H1_35->Draw("psame");
      }
      if (iQ2==3)
      {
        g_rcs_sys_vs_x_H1_60->Draw("samee2");
        g_rcs_stat_vs_x_H1_60->Draw("psame");
      }
      if (iQ2==4)
      {
        g_rcs_sys_vs_x_H1_120->Draw("samee2");
        g_rcs_stat_vs_x_H1_120->Draw("psame");
      }
      if (iQ2==5)
      {
        g_rcs_sys_vs_x_H1_200->Draw("samee2");
        g_rcs_stat_vs_x_H1_200->Draw("psame");
      }
      if (iQ2==6)
      {
        g_rcs_sys_vs_x_H1_300->Draw("samee2");
        g_rcs_stat_vs_x_H1_300->Draw("psame");
      }

      leg.AddEntry(g_rcs_stat_vs_x_H1_12,"H1","p");
      */

      leg.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.04);

      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.20,0.85,Form("Q^{2} = %.0f GeV^{2}",Q2_mid[iQ2]));
      // tl->DrawLatexNDC(0.20,0.80,"0.1 < y < 0.9");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sEventCounts_strange_charm_%d.pdf\")", cno-1, outDir, iQ2) );
    }
  }

  /*
  mclogx(cno++);
  {
    float plot_xrange_lo = 1E-3;
    float plot_xrange_hi = 1E0;
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1;
    for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
    {
      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("x");
      htemp.GetYaxis()->SetTitle("#tilde{#sigma}_{incl}");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      TLegend leg(0.48,0.7,0.9,0.87);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.035);
      leg.SetFillStyle(0);

      plot_Pythia6_ep->g_cs_vs_x_e[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_Pythia6_ep->g_cs_vs_x_e[iQ2],plot_Pythia6_ep->sys_latex,"p");

      plot_Pythia6_eAu->g_cs_vs_x_e[iQ2]->SetMarkerStyle(24);
      plot_Pythia6_eAu->g_cs_vs_x_e[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_Pythia6_eAu->g_cs_vs_x_e[iQ2],plot_Pythia6_eAu->sys_latex,"p");

      plot_BeAGLE_ep->g_cs_vs_x_e[iQ2]->SetLineColor(kBlue);
      plot_BeAGLE_ep->g_cs_vs_x_e[iQ2]->SetMarkerColor(kBlue);
      plot_BeAGLE_ep->g_cs_vs_x_e[iQ2]->SetMarkerStyle(21);
      plot_BeAGLE_ep->g_cs_vs_x_e[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_BeAGLE_ep->g_cs_vs_x_e[iQ2],plot_BeAGLE_ep->sys_latex,"p");

      plot_BeAGLE_eAu->g_cs_vs_x_e[iQ2]->SetLineColor(kBlue);
      plot_BeAGLE_eAu->g_cs_vs_x_e[iQ2]->SetMarkerColor(kBlue);
      plot_BeAGLE_eAu->g_cs_vs_x_e[iQ2]->SetMarkerStyle(25);
      plot_BeAGLE_eAu->g_cs_vs_x_e[iQ2]->Draw("cpsame");
      leg.AddEntry(plot_BeAGLE_eAu->g_cs_vs_x_c[iQ2],plot_BeAGLE_eAu->sys_latex,"p");


      leg.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.04);

      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.20,0.85,Form("Q^{2} = %.0f GeV^{2}",Q2_mid[iQ2]));
      // tl->DrawLatexNDC(0.20,0.80,"0.1 < y < 0.9");

      gROOT->ProcessLine( Form("cc%d->Print(\"EventCounts_incl_%d.pdf\")", cno-1, iQ2) );
    }
  }
  */


}

void plot_events(const char* inFile = "hists_eventcounts_ep.root", const char* outDir = "./")
{
  //set_charm_H1();

  set_Q2_binning();

  mcs(-1);

  // Pythia 6 e+p input
  TFile* fin1 = new TFile(inFile, "read");
  plot_Pythia6_ep = new PlotXsec(1);
  plot_Pythia6_ep->ReadEvtHists(fin1);
  plot_Pythia6_ep->CalculateReducedXsec();

  // Pythia 6 e+Au input
  // TFile* fin2 = new TFile("hists_nevt_Pythia6_eAu.root","read");
  // plot_Pythia6_eAu = new PlotXsec(2);
  // plot_Pythia6_eAu->ReadEvtHists(fin2);
  // plot_Pythia6_eAu->CalculateReducedXsec();

  // TFile* fin2 = new TFile("hists_nevt_Pythia6_eAu_10042.root","read");
  // plot_Pythia6_eAu = new PlotXsec(5);
  // plot_Pythia6_eAu->ReadEvtHists(fin2);
  // plot_Pythia6_eAu->CalculateReducedXsec();

  // TFile* fin2 = new TFile("hists_nevt_Pythia6_eAu_10042_EPSLO.root","read");
  // plot_Pythia6_eAu = new PlotXsec(7);
  // plot_Pythia6_eAu->ReadEvtHists(fin2);
  // plot_Pythia6_eAu->CalculateReducedXsec();

  // BeAGLE e+p input
  // TFile* fin3 = new TFile("hists_nevt_BeAGLE_ep.root","read");
  // plot_BeAGLE_ep = new PlotXsec(3);
  // plot_BeAGLE_ep->ReadEvtHists(fin3);
  // plot_BeAGLE_ep->CalculateReducedXsec();

  // BeAGLE e+Au input
  // TFile* fin4 = new TFile("hists_nevt_BeAGLE_INCoff_NLO_eAu.root","read");
  // plot_BeAGLE_eAu = new PlotXsec(6);
  // plot_BeAGLE_eAu->ReadEvtHists(fin4);
  // plot_BeAGLE_eAu->CalculateReducedXsec();

  // TFile* fin4 = new TFile("hists_nevt_BeAGLE_genShd1_eAu.root","read");
  // plot_BeAGLE_eAu = new PlotXsec(8);
  // plot_BeAGLE_eAu->ReadEvtHists(fin4);
  // plot_BeAGLE_eAu->CalculateReducedXsec();

  // TFile* fin4 = new TFile("hists_nevt_BeAGLE_eAu.root","read");
  // plot_BeAGLE_eAu = new PlotXsec(4);
  // plot_BeAGLE_eAu->ReadEvtHists(fin4);
  // plot_BeAGLE_eAu->CalculateReducedXsec();

  plot_graph(outDir);
}
