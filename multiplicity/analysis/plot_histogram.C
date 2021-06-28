void plot_histogram(const char* inFile, const char* outDir)
{
  // p_T distribution

  TFile* fin = new TFile(inFile,"read");
  cout<<"Print input file content"<<endl;
  fin->ls();

  TH1D* h1d_kaon = (TH1D*)fin->Get("h1d_kaon");
  TH1D* h1d_pion = (TH1D*)fin->Get("h1d_pion");
  TH1D* h1d_proton = (TH1D*)fin->Get("h1d_proton");


  // charged kaon multiplicity

  TCanvas* c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  //c1->SetLogy(); // set y axis to log scale

  h1d_kaon->GetXaxis()->SetTitle("charge kaon multiplicity (counts)");
  h1d_kaon->GetYaxis()->SetTitle("events");
  h1d_kaon->GetXaxis()->SetTitleOffset(1.3);
  h1d_kaon->GetYaxis()->SetTitleOffset(1.5);

  h1d_kaon->Draw("hsame");

  c1->SaveAs(outDir+"charged_kaon_mul.pdf");

  // charged pion multiplicity

  TCanvas* c2 = new TCanvas("c2","c2",800,800); // create new canvas
  c2->Range(0,0,1,1);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.1);
  //c2->SetLogy(); // set y axis to log scale

  h1d_pion->GetXaxis()->SetTitle("charge pion multiplicity (counts)");
  h1d_pion->GetYaxis()->SetTitle("events");
  h1d_pion->GetXaxis()->SetTitleOffset(1.3);
  h1d_pion->GetYaxis()->SetTitleOffset(1.5);

  h1d_pion->Draw("hsame");

  c2->SaveAs(outDir+"charged_pion_mul.pdf");

  // charged kaon multiplicity

  TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
  c3->Range(0,0,1,1);
  c3->SetLeftMargin(0.15);
  c3->SetBottomMargin(0.1);
  //c3->SetLogy(); // set y axis to log scale

  h1d_proton->GetXaxis()->SetTitle("proton multiplicity (counts)");
  h1d_proton->GetYaxis()->SetTitle("events");
  h1d_proton->GetXaxis()->SetTitleOffset(1.3);
  h1d_proton->GetYaxis()->SetTitleOffset(1.5);

  h1d_proton->Draw("hsame");

  c3->SaveAs(outDir+"proton_mul.pdf");


}
