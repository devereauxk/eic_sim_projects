const int etabin = 3;
const double eta_lo[etabin] = {-4,-1,1};
const double eta_hi[etabin] = {-1,1,4};
const int eta_color[etabin] = {kRed, kBlue, kGreen+1};
TH1D* h1d_part_pt = NULL;
TH2D* h2d_pt_vs_eta = NULL;
TH1D* h1d_pt_in_eta[etabin] = {0};
void slice_2D_hist()
{
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    // select the range of slice that you want to project
    h2d_pt_vs_eta->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
    // project to X axis and set name of the new 1D histograms (Form() function return a string), Make sure all the projected hisotgrams get a different name.
    h1d_pt_in_eta[ieta] = (TH1D*)h2d_pt_vs_eta->ProjectionX( Form("h1d_pt_in_eta_%d",ieta) );
    h1d_pt_in_eta[ieta]->SetName( Form("h1d_pt_in_eta_%d",ieta) );
    h1d_pt_in_eta[ieta]->SetTitle( "particle p_{T} ditributiton" );
  }
  // set back to original range
  h2d_pt_vs_eta->GetYaxis()->SetRangeUser(-4,4);
}
void plot_histogram_projection(const char* inFile = "output.root")
{
  TFile* fin = new TFile(inFile,"read");
  cout<<"Print input file content"<<endl;
  fin->ls();
  h1d_part_pt = (TH1D*)fin->Get("h1d_part_pt");
  h2d_pt_vs_eta = (TH2D*)fin->Get("h2d_pt_vs_eta");
  slice_2D_hist();
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
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    TCanvas* c2 = new TCanvas("c2","c2",800,800); // create new canvas
    c2->Range(0,0,1,1);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.1);
    c2->SetLogy(); // set y axis to log scale
    h1d_pt_in_eta[ieta]->GetXaxis()->SetRangeUser(0,5); // zoomin to pt range of 0 to 5GeV
    h1d_pt_in_eta[ieta]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    h1d_pt_in_eta[ieta]->GetYaxis()->SetTitle("Counts");
    h1d_pt_in_eta[ieta]->GetXaxis()->SetTitleOffset(1.3);
    h1d_pt_in_eta[ieta]->GetYaxis()->SetTitleOffset(1.5);
    h1d_pt_in_eta[ieta]->Draw("hsame");
    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.035);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.2,0.85,Form("#eta range %.1f to %.1f",eta_lo[ieta],eta_hi[ieta]));
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.2,0.80,Form("Total counts on plot is %.3e",h1d_pt_in_eta[ieta]->Integral()));
    c2->SaveAs( Form("part_pt_in_eta_%d.pdf",ieta) );
  }
  // plot different eta bin slices on the same canvas
  TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
  c3->Range(0,0,1,1);
  c3->SetLeftMargin(0.15);
  c3->SetBottomMargin(0.1);
  c3->SetLogy(); // set y axis to log scale
  TLegend* leg = new TLegend(0.60,0.70,0.80,0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);
  float plot_xrange_lo = 0, plot_xrange_hi = 5;
  float plot_yrange_lo = 1E0, plot_yrange_hi = 1.5*h1d_pt_in_eta[2]->GetMaximum(); // when using log axis, cannot use 0 as start as plot range
  // use the empty 2D histogram htemp as a frame
  TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
  htemp.SetStats(0); // not showing the box on the top right corner
  htemp.Draw();
  htemp.GetXaxis()->SetTitle("p_{T} [GeV/c]");
  htemp.GetYaxis()->SetTitle("Counts");
  htemp.GetXaxis()->SetTitleOffset(1.3);
  htemp.GetYaxis()->SetTitleOffset(1.5);
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    h1d_pt_in_eta[ieta]->SetStats(0); // not showing the box on the top right corner
    h1d_pt_in_eta[ieta]->SetLineColor(eta_color[ieta]); // eta_color array defined at the beginning
    h1d_pt_in_eta[ieta]->SetMarkerColor(eta_color[ieta]);
    h1d_pt_in_eta[ieta]->Draw("hsame");
    leg->AddEntry(h1d_pt_in_eta[ieta],Form("%.0f < #eta < %.0f",eta_lo[ieta],eta_hi[ieta]),"l"); // "l" means line
  }
  leg->Draw("same");
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.035);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.2,0.85,"e + p @ 10 + 110 GeV");
  c3->SaveAs( Form("part_pt_in_eta_all.pdf") );

  // all in one pdf
  TCanvas * c_all = new TCanvas("canvas_with_all_of_them", "canvas_with_all_of_them", 1600, 1600);
  c_all->Divide(2, 2);
  for (int ieta = 0; ieta < 3; ++ieta)
  {
    h1d_pt_in_eta[ieta]->SetLineColor(kBlack);
    c_all->cd(ieta+1); h1d_pt_in_eta[ieta]->Draw("hsame");
  }
  c_all->cd(4); htemp.Draw();
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    h1d_pt_in_eta[ieta]->SetStats(0); // not showing the box on the top right corner
    h1d_pt_in_eta[ieta]->SetLineColor(eta_color[ieta]); // eta_color array defined at the beginning
    h1d_pt_in_eta[ieta]->SetMarkerColor(eta_color[ieta]);
    h1d_pt_in_eta[ieta]->Draw("hsame");
  }
  c_all->SaveAs(Form("all_the_plots_on_one_pdf.pdf"));


}
