void plot_histogram(const char* inFile)
{
  // run with root -l 'plot_histogram("histogram_dir/access_tree_output.root", "histogram_dir/")'
  // p_T distribution

  TFile* fin = new TFile(inFile,"read");
  cout<<"Print input file content"<<endl;
  fin->ls();

  TH2D* h2d_kaon = (TH2D*)fin->Get("h2d_kaon");
  TH1D* h1d_kaon_pos = (TH1D*) h2d_kaon->ProjectionX("h1d_kaon_pos");
  TH1D* h1d_kaon_neg = (TH1D*) h2d_kaon->ProjectionY("h1d_kaon_neg");
  TH1D* h1d_kaon_total = new TH1D("h1d_kaon_total", "charged kaon multiplicity", 12, -0.5, 11.5);
  h1d_kaon_total->Add(h1d_kaon_pos, h1d_kaon_neg);

  TH2D* h2d_pion = (TH2D*)fin->Get("h2d_pion");
  TH1D* h1d_pion_pos = (TH1D*) h2d_pion->ProjectionX("h1d_pion_pos");
  TH1D* h1d_pion_neg = (TH1D*) h2d_pion->ProjectionY("h1d_pion_neg");
  TH1D* h1d_pion_total = new TH1D("h1d_pion_total", "charged pion multiplicity", 10000, -0.5, 9999.5);
  h1d_pion_total->Add(h1d_pion_pos, h1d_pion_neg);

  TH1D* h1d_proton = (TH1D*)fin->Get("h1d_proton");


  // charged kaon multiplicity

  TCanvas * c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_total->GetXaxis()->SetTitle("K^{#pm} multiplicity");
  h1d_kaon_total->GetYaxis()->SetTitle("events");
  h1d_kaon_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_total->Draw("hsame");
  c1->SaveAs("total_kaon_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_pos->GetXaxis()->SetTitle("K^{+} multiplicity");
  h1d_kaon_pos->GetYaxis()->SetTitle("events");
  h1d_kaon_pos->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_pos->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_pos->Draw("hsame");
  c1->SaveAs("pos_kaon_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_neg->GetXaxis()->SetTitle("K^{-} multiplicity");
  h1d_kaon_neg->GetYaxis()->SetTitle("events");
  h1d_kaon_neg->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_neg->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_neg->Draw("hsame");
  c1->SaveAs("neg_kaon_mul.pdf");

  // charged pion multiplicity

  TCanvas* c2 = new TCanvas("c2","c2",800,800); // create new canvas
  c2->Range(0,0,1,1);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.1);
  h1d_pion_total->GetXaxis()->SetTitle("#pi ^ {#pm} multiplicity");
  h1d_pion_total->GetYaxis()->SetTitle("events");
  h1d_pion_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_total->Draw("hsame");
  c2->SaveAs("total_pion_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_pion_pos->GetXaxis()->SetTitle("#pi ^{+} multiplicity");
  h1d_pion_pos->GetYaxis()->SetTitle("events");
  h1d_pion_pos->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_pos->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_pos->Draw("hsame");
  c1->SaveAs("pos_pion_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_pion_neg->GetXaxis()->SetTitle("#pi ^{-} multiplicity");
  h1d_pion_neg->GetYaxis()->SetTitle("events");
  h1d_pion_neg->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_neg->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_neg->Draw("hsame");
  c1->SaveAs("neg_pion_mul.pdf");


  // proton multiplicity

  TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
  c3->Range(0,0,1,1);
  c3->SetLeftMargin(0.15);
  c3->SetBottomMargin(0.1);
  h1d_proton->GetXaxis()->SetTitle("proton multiplicity");
  h1d_proton->GetYaxis()->SetTitle("events");
  h1d_proton->GetXaxis()->SetTitleOffset(1.3);
  h1d_proton->GetYaxis()->SetTitleOffset(1.5);
  h1d_proton->Draw("hsame");
  c3->SaveAs("proton_mul.pdf");


}
/*
THStack *hs = new THStack("hs", "");
hs->Add(h1d_kaon_pos);
hs->Add(h1d_kaon_neg);
TCanvas *cs = new TCanvas("cs", "cs", 1400, 700);
cs->Divide(2);
cs->cd(1); h2d_kaon->Draw("colz");
cs->cd(2); hs->Draw("nostackb");
cs->SaveAs("comp_kaon_mul.pdf");
*/
