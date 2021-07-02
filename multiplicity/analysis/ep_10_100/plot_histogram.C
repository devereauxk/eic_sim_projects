void insert_text(TH1D * h1d)
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.035);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.2,0.85,"e+p @ 10GeV, 100GeV");
  tl->SetTextColor(kBlue);
  tl->DrawLatexNDC(0.2,0.80,Form("Total counts on plot is %.3e",h1d->GetEntries()));
}

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
  TH1D* h1d_kaon_total = (TH1D*)fin->Get("h1d_kaon_total");

  TH2D* h2d_pion = (TH2D*)fin->Get("h2d_pion");
  TH1D* h1d_pion_pos = (TH1D*) h2d_pion->ProjectionX("h1d_pion_pos");
  TH1D* h1d_pion_neg = (TH1D*) h2d_pion->ProjectionY("h1d_pion_neg");
  TH1D* h1d_pion_total = (TH1D*)fin->Get("h1d_pion_total");

  TH2D* h2d_proton = (TH2D*)fin->Get("h2d_proton");
  TH1D* h1d_proton = (TH1D*) h2d_proton->ProjectionX("h1d_proton");
  TH1D* h1d_anti_proton = (TH1D*) h2d_proton->ProjectionY("h1d_anti_proton");
  TH1D* h1d_proton_total = (TH1D*)fin->Get("h1d_proton_total");


  // charged kaon multiplicity

  TCanvas * c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_total->GetXaxis()->SetTitle("K^{#pm} multiplicity [counts]");
  h1d_kaon_total->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_kaon_total->Scale(1 / h1d_kaon_total->GetEntries());
  h1d_kaon_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_total->Draw("hsame");
  insert_text(h1d_kaon_total);
  c1->SaveAs("total_kaon_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_pos->GetXaxis()->SetTitle("K^{+} multiplicity [counts]");
  h1d_kaon_pos->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_kaon_pos->Scale(1 / h1d_kaon_pos->GetEntries());
  h1d_kaon_pos->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_pos->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_pos->Draw("hsame");
  insert_text(h1d_kaon_pos);
  c1->SaveAs("pos_kaon_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_kaon_neg->GetXaxis()->SetTitle("K^{-} multiplicity [counts]");
  h1d_kaon_neg->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_kaon_neg->Scale(1 / h1d_kaon_neg->GetEntries());
  h1d_kaon_neg->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_neg->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_neg->Draw("hsame");
  insert_text(h1d_kaon_neg);
  c1->SaveAs("neg_kaon_mul.pdf");

  // charged pion multiplicity [counts]

  TCanvas* c2 = new TCanvas("c2","c2",800,800); // create new canvas
  c2->Range(0,0,1,1);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.1);
  h1d_pion_total->GetXaxis()->SetTitle("#pi^{#pm} multiplicity [counts]");
  h1d_pion_total->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_pion_total->Scale(1 / h1d_pion_total->GetEntries());
  h1d_pion_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_total->Draw("hsame");
  insert_text(h1d_pion_total);
  c2->SaveAs("total_pion_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_pion_pos->GetXaxis()->SetTitle("#pi^{+} multiplicity [counts]");
  h1d_pion_pos->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_pion_pos->Scale(1 / h1d_pion_pos->GetEntries());
  h1d_pion_pos->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_pos->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_pos->Draw("hsame");
  insert_text(h1d_pion_pos);
  c1->SaveAs("pos_pion_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_pion_neg->GetXaxis()->SetTitle("#pi^{-} multiplicity [counts]");
  h1d_pion_neg->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_pion_neg->Scale(1 / h1d_pion_neg->GetEntries());
  h1d_pion_neg->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_neg->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_neg->Draw("hsame");
  insert_text(h1d_pion_neg);
  c1->SaveAs("neg_pion_mul.pdf");


  // proton multiplicity

  TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
  c3->Range(0,0,1,1);
  c3->SetLeftMargin(0.15);
  c3->SetBottomMargin(0.1);
  h1d_proton->GetXaxis()->SetTitle("p+#overline{p} multiplicity [counts]");
  h1d_proton->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_proton->Scale(1 / h1d_proton->GetEntries());
  h1d_proton->GetXaxis()->SetTitleOffset(1.3);
  h1d_proton->GetYaxis()->SetTitleOffset(1.5);
  h1d_proton->Draw("hsame");
  insert_text(h1d_proton_total);
  c3->SaveAs("total_proton_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_proton->GetXaxis()->SetTitle("p [counts]");
  h1d_proton->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_proton->Scale(1 / h1d_proton->GetEntries());
  h1d_proton->GetXaxis()->SetTitleOffset(1.3);
  h1d_proton->GetYaxis()->SetTitleOffset(1.5);
  h1d_proton->Draw("hsame");
  insert_text(h1d_proton);
  c1->SaveAs("proton_mul.pdf");

  c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  h1d_anti_proton->GetXaxis()->SetTitle("#overline{p} [counts]");
  h1d_anti_proton->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_anti_proton->Scale(1 / h1d_anti_proton->GetEntries());
  h1d_anti_proton->GetXaxis()->SetTitleOffset(1.3);
  h1d_anti_proton->GetYaxis()->SetTitleOffset(1.5);
  h1d_anti_proton->Draw("hsame");
  insert_text(h1d_anti_proton);
  c1->SaveAs("anti_proton_mul.pdf");


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
