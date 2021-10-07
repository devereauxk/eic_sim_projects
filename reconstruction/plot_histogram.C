void plot_histogram(const char* inFile, const char* outDir)
{

  TFile* f = new TFile(inFile,"read");

  TCanvas* c = new TCanvas("c3","c3",800,800);
  c->Range(0,0,1,1);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.1);

  TH1D* x = (TH1D*) ((TH2D*) f->Get("fg2d_Kpimass_vs_p_2_0"))->ProjectionX("x")
  x->Draw("hsame");
  x->GetXaxis()->SetRangeUser(1.7,2);
  c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_0_proj.pdf", outDir));

  TH2D* x2 = (TH2D*) f->Get("fg2d_Kpimass_vs_p_2_0");
  x2->Draw("colz");
  c->SetLogz();
  c->SaveAs(Form("%sfg2d_Kpimass_vs_p_2_0.pdf", outDir));

}
