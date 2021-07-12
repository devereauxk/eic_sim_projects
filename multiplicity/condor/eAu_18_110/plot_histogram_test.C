void insert_text(TH1D * h1d)
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.035);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.2,0.85,"e+Au @ 18GeV, 110GeV");
  tl->SetTextColor(kBlue);
  tl->DrawLatexNDC(0.2,0.80,Form("Total counts on plot is %.3e",h1d->GetEntries()));
}

void plot_histogram_test(const char* inFile)
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

  TCanvas * c_all = new TCanvas("c_all", "c_all", 2400, 1600);
  c_all->Divide(3, 2);

  // proton multiplicity

  c_all->cd(1);
  h1d_proton_total->GetXaxis()->SetTitle("proton+antiproton multiplicity [counts]");
  h1d_proton_total->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_proton_total->Scale(1 / h1d_proton_total->GetEntries());
  h1d_proton_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_proton_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_proton_total->SetStats(0);
  h1d_proton_total->Draw("hsame");

  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.035);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.2,0.85,"e + Au @ 18 + 110 GeV");
  tl->DrawLatexNDC(0.2,1.2,Form("$.0f events", h1d_proton_total->GetEntries()));
  tl->DrawLatexNDC(0.2,1.4,Form("avg: ", h1d_proton_total->GetMean());

  c_all->cd(4);
  TH2F htemp("htemp","",12,-0.5,11.5,10,0,1.2*fmaxf(h1d_proton->GetMaximum(), h1d_anti_proton->GetMaximum()));
  htemp.SetStats(0);
  htemp.Draw();
  htemp.GetXaxis()->SetTitle("multiplicity [counts]");
  htemp.GetYaxis()->SetTitle("fraction of events [%]");
  htemp.GetXaxis()->SetTitleOffset(1.3);
  htemp.GetYaxis()->SetTitleOffset(1.5);

  TLegend* leg = new TLegend(0.60,0.70,0.80,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);

  h1d_proton->SetLineColor(kRed);
  h1d_proton->SetStats(0);
  h1d_proton->Scale(1 / h1d_proton->GetEntries());
  h1d_proton->Draw("hsame");
  leg->AddEntry(h1d_proton,"proton","l");

  h1d_anti_proton->SetLineColor(kBlue);
  h1d_anti_proton->SetStats(0);
  h1d_anti_proton->Scale(1 / h1d_anti_proton->GetEntries());
  h1d_anti_proton->Draw("hsame");
  leg->AddEntry(h1d_proton,"antiproton","l");

  c_all->SaveAs("temp.pdf");


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
