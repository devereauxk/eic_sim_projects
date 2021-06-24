void plot_histogram(const char* inFile = "output.root")
{
  // p_T distribution

  TFile* fin = new TFile(inFile,"read");
  cout<<"Print input file content"<<endl;
  fin->ls();

  TH1D* h1d_part_pt = (TH1D*)fin->Get("h1d_part_pt");

  TCanvas* c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  c1->SetLogy(); // set y axis to log scale

  h1d_part_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h1d_part_pt->GetYaxis()->SetTitle("Counts");
  h1d_part_pt->GetXaxis()->SetTitleOffset(1.3);
  h1d_part_pt->GetYaxis()->SetTitleOffset(1.5);

  h1d_part_pt->Draw("hsame");

  c1->SaveAs("part_pt.pdf");


  //p_T vs eta

  TH2D* h2d_pt_vs_eta = (TH2D*)fin->Get("h2d_pt_vs_eta");

  TCanvas* c2 = new TCanvas("c2","c2",800,800); // create new canvas
  c2->Range(0,0,1,1);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.1);

  h2d_pt_vs_eta->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h2d_pt_vs_eta->GetYaxis()->SetTitle("eta");
  h2d_pt_vs_eta->GetXaxis()->SetTitleOffset(1.3);
  h2d_pt_vs_eta->GetYaxis()->SetTitleOffset(1.5);

  h2d_pt_vs_eta->Draw("colz");

  c2->SaveAs("pt_vs_eta.pdf");

}
