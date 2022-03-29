R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

static int cno = 0;

void plot_histogram(const char* inFile, const char* outDir, const char* title = "", const double events = 1000000)
{

  TFile* f = new TFile(inFile,"read");

  TTree* t = (TTree*) f->Get("EICTree");

  TGaxis::SetMaxDigits(3);
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mcs(cno++);
    {
      TH1D* fg1d_Kpimass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta)))->ProjectionX("fg1d");
      TH1D* bg1d_Kpimass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("bg2d_Kpimass_vs_p_2_%d", ieta)))->ProjectionX("bg1d");
      TH1D* sg1d_Kpimass_vs_p = (TH1D*) fg1d_Kpimass_vs_p->Clone();
      sg1d_Kpimass_vs_p->SetName("sg1d");
      sg1d_Kpimass_vs_p->Add(bg1d_Kpimass_vs_p, -1);

      float plot_xrange_lo = 1.7;
      float plot_xrange_hi = 2.05;
      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1.7*fg1d_Kpimass_vs_p->GetMaximum();

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("m_{K^{#pm}#pi^{#pm}} [GeV/c^{2}]");
      htemp.GetYaxis()->SetTitle("Counts");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      TLegend leg(0.2,0.65,0.83,0.80);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      fg1d_Kpimass_vs_p->SetLineColor(kBlue);
      bg1d_Kpimass_vs_p->SetLineColor(kRed);
      sg1d_Kpimass_vs_p->SetLineColor(kGreen+1);
      fg1d_Kpimass_vs_p->Draw("hsame");
      bg1d_Kpimass_vs_p->Draw("hsame");
      sg1d_Kpimass_vs_p->Draw("hsame");

      float temp_mean = -9999;
      float temp_sigma = 0;
      sg1d_Kpimass_vs_p->Fit("gaus","0R","",1.8,1.95);
      TF1* gaus = sg1d_Kpimass_vs_p->GetFunction("gaus");

      if (gaus!=NULL)
      {
        temp_mean = gaus->GetParameter(1);
        temp_sigma = gaus->GetParameter(2);

        TF1* func_peak = new TF1("func_peak","[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*pow(x-[4],2)/pow([5],2))",temp_mean-3*temp_sigma,temp_mean+3*temp_sigma);
        //func_peak->SetLineColor(pid_color[ipid]);
        func_peak->FixParameter(3,gaus->GetParameter(0));
        func_peak->FixParameter(4,gaus->GetParameter(1));
        func_peak->FixParameter(5,gaus->GetParameter(2));
        fg1d_Kpimass_vs_p->Fit(func_peak,"R","",temp_mean-8*temp_sigma,temp_mean+8*temp_sigma);
      }

      float int_range_lo = sg1d_Kpimass_vs_p->FindBin(temp_mean-3*temp_sigma);
      float int_range_hi = sg1d_Kpimass_vs_p->FindBin(temp_mean+3*temp_sigma);
      float N_SG = sg1d_Kpimass_vs_p->Integral(int_range_lo, int_range_hi);
      float N_BG = bg1d_Kpimass_vs_p->Integral(int_range_lo, int_range_hi);

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.2,0.85,title);
      tl->DrawLatexNDC(0.2,0.80,Form("%.1e events",events));
      tl->DrawLatexNDC(0.2,0.75,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
      tl->DrawLatexNDC(0.2,0.70,Form("signal = %.1e, background = %.1e", N_SG, N_BG));
      gROOT->ProcessLine( Form("cc%d->Print(\"%skpimass_vs_p_%d.pdf\")", cno-1, outDir, ieta) );

      gROOT->ProcessLine( Form("cc%d->Clear()", cno-1) );

      //Lc
      TH1D* fg1d_Kpipmass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("fg2d_Kpipmass_vs_p_2_%d", ieta)))->ProjectionX("fg1d");
      TH1D* bg1d_Kpipmass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("bg2d_Kpipmass_vs_p_2_%d", ieta)))->ProjectionX("bg1d");
      TH1D* sg1d_Kpipmass_vs_p = (TH1D*) fg1d_Kpimass_vs_p->Clone();
      sg1d_Kpipmass_vs_p->SetName("sg1d");
      sg1d_Kpipmass_vs_p->Add(bg1d_Kpipmass_vs_p, -1);

      plot_xrange_lo = 2.1;
      plot_xrange_hi = 2.45;
      plot_yrange_lo = 0;
      plot_yrange_hi = 1.7*fg1d_Kpipmass_vs_p->GetMaximum();

      TH2F htemp2("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp2.Draw();
      htemp2.GetXaxis()->SetTitle("m_{K^{#pm}#pi^{#mp}p^{#mp}} [GeV/c^{2}]");
      htemp2.GetYaxis()->SetTitle("Counts");
      myhset(&htemp2,1.2,1.6,0.05,0.045);

      TLegend leg2(0.2,0.65,0.83,0.80);
      leg2.SetBorderSize(0);
      leg2.SetTextSize(0.03);
      leg2.SetFillStyle(0);
      leg2.SetMargin(0.1);

      fg1d_Kpipmass_vs_p->SetLineColor(kBlue);
      bg1d_Kpipmass_vs_p->SetLineColor(kRed);
      sg1d_Kpipmass_vs_p->SetLineColor(kGreen+1);
      fg1d_Kpipmass_vs_p->Draw("hsame");
      bg1d_Kpipmass_vs_p->Draw("hsame");
      sg1d_Kpipmass_vs_p->Draw("hsame");

      temp_mean = -9999;
      temp_sigma = 0;
      sg1d_Kpipmass_vs_p->Fit("gaus","0R","",2.1,2.4);
      gaus = sg1d_Kpipmass_vs_p->GetFunction("gaus");

      if (gaus!=NULL)
      {
        temp_mean = gaus->GetParameter(1);
        temp_sigma = gaus->GetParameter(2);

        TF1* func_peak = new TF1("func_peak","[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*pow(x-[4],2)/pow([5],2))",temp_mean-3*temp_sigma,temp_mean+3*temp_sigma);
        //func_peak->SetLineColor(pid_color[ipid]);
        func_peak->FixParameter(3,gaus->GetParameter(0));
        func_peak->FixParameter(4,gaus->GetParameter(1));
        func_peak->FixParameter(5,gaus->GetParameter(2));
        fg1d_Kpipmass_vs_p->Fit(func_peak,"R","",temp_mean-8*temp_sigma,temp_mean+8*temp_sigma);
      }

      int_range_lo = sg1d_Kpipmass_vs_p->FindBin(temp_mean-3*temp_sigma);
      int_range_hi = sg1d_Kpipmass_vs_p->FindBin(temp_mean+3*temp_sigma);
      N_SG = sg1d_Kpipmass_vs_p->Integral(int_range_lo, int_range_hi);
      N_BG = bg1d_Kpipmass_vs_p->Integral(int_range_lo, int_range_hi);

      tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.2,0.85,title);
      tl->DrawLatexNDC(0.2,0.80,Form("%.1e events",events));
      tl->DrawLatexNDC(0.2,0.75,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
      tl->DrawLatexNDC(0.2,0.70,Form("signal = %.1e, background = %.1e", N_SG, N_BG));
      gROOT->ProcessLine( Form("cc%d->Print(\"%skpipmass_vs_p_%d.pdf\")", cno-1, outDir, ieta) );


    }

  }

}
