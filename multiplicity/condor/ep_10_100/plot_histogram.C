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

  // proton multiplicity ------------------------------------------------------------------------
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
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.4,0.80,"e + p @ 10 + 100 GeV");
  tl->DrawLatexNDC(0.4,0.75,Form("%.0f events", h1d_proton_total->GetEntries()));
  tl->DrawLatexNDC(0.4,0.70,Form("avg: %1.4f", h1d_proton_total->GetMean()));

  c_all->cd(4);
  TH2F htemp("htemp","",12,-0.5,11.5,10,0,1.2*fmaxf(h1d_proton->GetMaximum() / h1d_proton->GetEntries(), h1d_anti_proton->GetMaximum() / h1d_anti_proton->GetEntries()));
  htemp.SetStats(0);
  htemp.Draw();
  htemp.GetXaxis()->SetTitle("multiplicity [counts]");
  htemp.GetYaxis()->SetTitle("fraction of events [%]");
  htemp.GetXaxis()->SetTitleOffset(1.3);
  htemp.GetYaxis()->SetTitleOffset(1.5);

  TLegend* leg = new TLegend(0.35,0.70,0.80,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);
  tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.35,0.60,Form("proton avg: %1.4f", h1d_proton->GetMean()));
  tl->DrawLatexNDC(0.35,0.55,Form("antiproton avg: %1.4f", h1d_anti_proton->GetMean()));

  h1d_proton->SetLineColor(kRed);
  h1d_proton->SetStats(0);
  h1d_proton->Scale(1 / h1d_proton->GetEntries());
  h1d_proton->Draw("hsame");
  leg->AddEntry(h1d_proton,"proton","l");

  h1d_anti_proton->SetLineColor(kBlue);
  h1d_anti_proton->SetStats(0);
  h1d_anti_proton->Scale(1 / h1d_anti_proton->GetEntries());
  h1d_anti_proton->Draw("hsame");
  leg->AddEntry(h1d_anti_proton,"antiproton","l");

  leg->Draw("same");

  // kaon multiplicity ----------------------------------------------------------------------------
  c_all->cd(2);
  h1d_kaon_total->GetXaxis()->SetTitle("K^{#pm} multiplicity [counts]");
  h1d_kaon_total->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_kaon_total->Scale(1 / h1d_kaon_total->GetEntries());
  h1d_kaon_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_kaon_total->SetStats(0);
  h1d_kaon_total->Draw("hsame");

  tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.4,0.80,"e + p @ 10 + 100 GeV");
  tl->DrawLatexNDC(0.4,0.75,Form("%.0f events", h1d_kaon_total->GetEntries()));
  tl->DrawLatexNDC(0.4,0.70,Form("avg: %1.4f", h1d_kaon_total->GetMean()));

  c_all->cd(5);
  TH2F htemp5("htemp5","",12,-0.5,11.5,10,0,1.2*fmaxf(h1d_kaon_pos->GetMaximum() / h1d_kaon_pos->GetEntries(), h1d_kaon_neg->GetMaximum() / h1d_kaon_neg->GetEntries()));
  htemp5.SetStats(0);
  htemp5.Draw();
  htemp5.GetXaxis()->SetTitle("multiplicity [counts]");
  htemp5.GetYaxis()->SetTitle("fraction of events [%]");
  htemp5.GetXaxis()->SetTitleOffset(1.3);
  htemp5.GetYaxis()->SetTitleOffset(1.5);

  leg = new TLegend(0.35,0.70,0.80,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);
  tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.35,0.60,Form("K^{+} avg: %1.4f", h1d_kaon_pos->GetMean()));
  tl->DrawLatexNDC(0.35,0.55,Form("K^{-} avg: %1.4f", h1d_kaon_neg->GetMean()));

  h1d_kaon_pos->SetLineColor(kRed);
  h1d_kaon_pos->SetStats(0);
  h1d_kaon_pos->Scale(1 / h1d_kaon_pos->GetEntries());
  h1d_kaon_pos->Draw("hsame");
  leg->AddEntry(h1d_kaon_pos,"K^{+}","l");

  h1d_kaon_neg->SetLineColor(kBlue);
  h1d_kaon_neg->SetStats(0);
  h1d_kaon_neg->Scale(1 / h1d_kaon_neg->GetEntries());
  h1d_kaon_neg->Draw("hsame");
  leg->AddEntry(h1d_kaon_neg,"K^{-}","l");

  leg->Draw("same");




  //pion multiplicity ---------------------------------------------------------------------------
  c_all->cd(3);
  h1d_pion_total->GetXaxis()->SetTitle("#pi^{#pm} multiplicity [counts]");
  h1d_pion_total->GetYaxis()->SetTitle("fraction of events [%]");
  h1d_pion_total->Scale(1 / h1d_pion_total->GetEntries());
  h1d_pion_total->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion_total->GetYaxis()->SetTitleOffset(1.5);
  h1d_pion_total->SetStats(0);
  h1d_pion_total->Draw("hsame");

  tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.4,0.80,"e + p @ 10 + 100 GeV");
  tl->DrawLatexNDC(0.4,0.75,Form("%.0f events", h1d_pion_total->GetEntries()));
  tl->DrawLatexNDC(0.4,0.70,Form("avg: %1.4f", h1d_pion_total->GetMean()));

  c_all->cd(6);
  TH2F htemp6("htemp6","", 40, -0.5, 39.5,10,0,1.2*fmaxf(h1d_pion_pos->GetMaximum() / h1d_pion_pos->GetEntries(), h1d_pion_neg->GetMaximum() / h1d_pion_neg->GetEntries()));
  htemp6.SetStats(0);
  htemp6.Draw();
  htemp6.GetXaxis()->SetTitle("multiplicity [counts]");
  htemp6.GetYaxis()->SetTitle("fraction of events [%]");
  htemp6.GetXaxis()->SetTitleOffset(1.3);
  htemp6.GetYaxis()->SetTitleOffset(1.5);

  leg = new TLegend(0.35,0.70,0.80,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);
  tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.045);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.35,0.60,Form("#pi^{+} avg: %1.4f", h1d_pion_pos->GetMean()));
  tl->DrawLatexNDC(0.35,0.55,Form("#pi^{-} avg: %1.4f", h1d_pion_neg->GetMean()));

  h1d_pion_pos->SetLineColor(kRed);
  h1d_pion_pos->SetStats(0);
  h1d_pion_pos->Scale(1 / h1d_pion_pos->GetEntries());
  h1d_pion_pos->Draw("hsame");
  leg->AddEntry(h1d_pion_pos,"#pi^{+}","l");

  h1d_pion_neg->SetLineColor(kBlue);
  h1d_pion_neg->SetStats(0);
  h1d_pion_neg->Scale(1 / h1d_pion_neg->GetEntries());
  h1d_pion_neg->Draw("hsame");
  leg->AddEntry(h1d_pion_neg,"#pi^{-}","l");

  leg->Draw("same");


  // save plot ------------------------------------------------------------------------
  c_all->SaveAs("multiplicities.pdf");


}
