#include "bins.h"

void plot_histogram(const char* inFile, const char* outDir, const char* title = "", const double events = 1000000)
{

  TFile* f = new TFile(inFile,"read");

  TTree* t = (TTree*) f->Get("EICTree");

  TCanvas* c = new TCanvas("c3","c3",800,800);
  c->Range(0,0,1,1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);

  for (int ieta = 0; ieta < etabin; ieta++) {

    // D0
    c = new TCanvas("c3","c3",800,800);
    TH2D* fg2d_Kpimass_vs_p = (TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta));
    fg2d_Kpimass_vs_p->Draw("colz");
    fg2d_Kpimass_vs_p->SetStats(0);
    c->SetLogz();
    c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_%d.pdf", outDir, ieta));

    TH1D* fg1d_Kpimass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta)))->ProjectionX("x");
    TH1D* bg1d_Kpimass_vs_p = (TH1D*) ((TH2D*) f->Get(Form("bg2d_Kpimass_vs_p_2_%d", ieta)))->ProjectionX("x");
    TH1D* sg1d_Kpimass_vs_p = (TH1D*) fg1d_Kpimass_vs_p->Clone();
    sg1d_Kpimass_vs_p->SetName(Form("sg2d_Kpimass_vs_p_2_%d", ieta));
    sg1d_Kpimass_vs_p->Add(bg1d_Kpimass_vs_p, -1);

    c = new TCanvas("c3","c3",800,800);
    fg1d_Kpimass_vs_p->Draw("hsame");
    fg1d_Kpimass_vs_p->GetXaxis()->SetRangeUser(1.75,1.95);
    fg1d_Kpimass_vs_p->SetTitle("");
    fg1d_Kpimass_vs_p->GetXaxis()->SetTitle("M(K^{#pm}#pi^{#mp}) [GeV/c^{2}]");
    fg1d_Kpimass_vs_p->GetYaxis()->SetTitle("counts");
    fg1d_Kpimass_vs_p->GetXaxis()->SetTitleOffset(1.3);
    fg1d_Kpimass_vs_p->GetYaxis()->SetTitleOffset(1.5);
    fg1d_Kpimass_vs_p->SetStats(0);

    /*
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
    */

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.035);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.15,0.85,title);
    tl->DrawLatexNDC(0.15,0.80,Form("%.1e events",events));
    tl->DrawLatexNDC(0.15,0.75,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
    //tl->DrawLatexNDC(0.15,0.70,Form("signal = %.1e, background = %.1e", N_SG, N_BG));
    c->SaveAs(Form("%sfg1d_Kpimass_2_%d.pdf", outDir, ieta));
    //c->SaveAs(Form("%sfg1d_Kpimass_2_%d.pdf", outDir, ieta));



    /*
    //Lc
    c = new TCanvas("c3","c3",800,800);
    x = (TH1D*) ((TH2D*) f->Get(Form("fg2d_Kpipmass_vs_p_2_%d", ieta)))->ProjectionX("x");
    x->Draw("hsame");
    x->GetXaxis()->SetRangeUser(2.2,2.4);
    x->SetTitle("");
    x->GetXaxis()->SetTitle("M(K^{#pm}#pi^{#mp}p^{#mp}) [GeV/c^{2}]");
    x->GetYaxis()->SetTitle("counts");
    x->GetXaxis()->SetTitleOffset(1.3);
    x->GetYaxis()->SetTitleOffset(1.5);
    x->SetStats(0);
    tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.035);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.15,0.85,title);
    tl->DrawLatexNDC(0.15,0.80,Form("%.1e events",events));
    tl->DrawLatexNDC(0.15,0.75,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
    c->SaveAs(Form("%sfg2d_Kpipmass_2_%d.pdf", outDir, ieta));

    x2 = (TH2D*) f->Get(Form("fg2d_Kpipmass_vs_p_2_%d", ieta));
    x2->Draw("colz");
    x->SetStats(0);
    c->SetLogz();
    c->SaveAs(Form("%sfg2d_Kpipmass_vs_p_2_%d.pdf", outDir, ieta));*/


  }

}
