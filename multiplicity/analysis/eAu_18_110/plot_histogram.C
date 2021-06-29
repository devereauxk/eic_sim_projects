void plot_h1d(TH1D * h1d, const char * x, const char * y, const char * filename)
{
  TCanvas* c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  //c1->SetLogy(); // set y axis to log scale

  h1d->GetXaxis()->SetTitle(x);
  h1d->GetYaxis()->SetTitle(y);
  h1d->GetXaxis()->SetTitleOffset(1.3);
  h1d->GetYaxis()->SetTitleOffset(1.5);

  h1d->Draw("hsame");

  c1->SaveAs(filename);
}

void plot_h2d(TH1D * h2d, const char * x, const char * y, const char * filename)
{
  TCanvas* c1 = new TCanvas("c1","c1",800,800); // create new canvas
  c1->Range(0,0,1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.1);
  //c1->SetLogy(); // set y axis to log scale

  h2d->GetXaxis()->SetTitle(x);
  h2d->GetYaxis()->SetTitle(y);
  h2d->GetXaxis()->SetTitleOffset(1.3);
  h2d->GetYaxis()->SetTitleOffset(1.5);

  h2d->Draw("colz");

  c1->SaveAs(filename);
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
  TH1D* h1d_kaon_total = h1d_kaon_pos + h1d_kaon_neg;
  TH2D* h2d_pion = (TH2D*)fin->Get("h2d_pion");
  TH1D* h1d_pion_pos = (TH1D*) h2d_pion->ProjectionX("h1d_pion_pos");
  TH1D* h1d_pion_neg = (TH1D*) h2d_pion->ProjectionY("h1d_pion_neg");
  TH1D* h1d_pion_total = h1d_pion_pos + h1d_pion_neg;
  TH1D* h1d_proton = (TH1D*)fin->Get("h1d_proton");


  // charged kaon multiplicity

  plot_h2d(h2d_kaon, "positive kaon multiplicity (counts)", "negative kaon multiplicity (counts)", "kaon_mul_2d.pdf");
  plot_h1d(h1d_kaon_total, "charged kaon multiplicity (counts)", "events", "kaon_mul_2d.pdf");


  // charged pion multiplicity

  plot_h2d(h2d_pion, "positive pion multiplicity (counts)", "negative pion multiplicity (counts)", "pion_mul_2d.pdf");
  plot_h1d(h1d_pion_total, "charged pion multiplicity (counts)", "events", "pion_mul_2d.pdf");


  // proton multiplicity

  plot_h1d(h1d_proton, "proton multiplicity (counts)", "events", "proton_mul.pdf");


}
