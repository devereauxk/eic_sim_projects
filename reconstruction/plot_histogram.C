void plot_histogram(const char* inFile, const char* outDir)
{

  TFile* f = new TFile(inFile,"read");

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
    c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_%d_proj.pdf", outDir, ieta));

    TH2D* x2 = (TH2D*) f->Get(Form("fg2d_Kpimass_vs_p_2_%d", ieta));
    x2->Draw("colz");
    c->SetLogz();
    c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_%d.pdf", outDir, ieta));

  }

}
