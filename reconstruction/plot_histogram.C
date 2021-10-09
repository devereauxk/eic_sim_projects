const int etabin = 3;
static double eta_lo[etabin] = {-3,-1,1};
static double eta_hi[etabin] = {-1,1,3};

void plot_histogram(const char* inFile, const char* outDir, const char* title = "", const double events = 1000000)
{

  TFile* f = new TFile(inFile,"read");

  TTree* t = (TTree*) f->Get("EICTree");

  TCanvas* c = new TCanvas("c3","c3",800,800);
  c->Range(0,0,1,1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);

  for (int ieta = 0; ieta < 3; ieta++) {

    c = new TCanvas("c3","c3",800,800);
    TH1D* x = (TH1D*) ((TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta)))->ProjectionX("x");
    x->Draw("hsame");
    x->GetXaxis()->SetRangeUser(1.7,2);
    x->SetTitle("");
    x->GetXaxis()->SetTitle("M(K^{#pm}#pi^{#mp}) [GeV/c^{2}]");
    x->GetYaxis()->SetTitle("counts");
    x->GetXaxis()->SetTitleOffset(1.3);
    x->GetYaxis()->SetTitleOffset(1.5);
    x->SetStats(0);
    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.035);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.15,0.85,title);
    tl->DrawLatexNDC(0.15,0.80,Form("%.1e events",events));
    tl->DrawLatexNDC(0.15,0.75,Form("%.0f < #eta < %.0f",eta_lo[ieta],eta_hi[ieta]));
    c->SaveAs(Form("%sfg2d_Kpimass_2_%d.pdf", outDir, ieta));

    TH2D* x2 = (TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta));
    x2->Draw("colz");
    c->SetLogz();
    x->SetTitle("");
    x->GetXaxis()->SetTitle("M(K^{#pm}#pi^{#mp}) [GeV/c^{2}]");
    x->GetYaxis()->SetTitle("<p_{T}>(K^{#pm}#pi^{#mp}) [GeV/c]");
    x->GetXaxis()->SetTitleOffset(1.3);
    x->GetYaxis()->SetTitleOffset(1.5);
    x->SetStats(0);
    c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_%d.pdf", outDir, ieta));

  }

}
