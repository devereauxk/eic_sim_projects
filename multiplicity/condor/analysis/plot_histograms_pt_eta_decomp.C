const int etabin = 7;
const double eta_lo[etabin] = {-3,-1,1,-3,-4,-7,-10};
const double eta_hi[etabin] = {-1,1,3,3,4,7,10};
const int eta_color[etabin] = {kRed, kBlue, kGreen+1, kBlack, kOrange, kYellow, kCyan};
TH2D* h2d_pt_vs_eta = NULL;
TH1D* h1d_pt_in_eta[etabin] = {0};
TGraphErrors* g_pt_vs_eta = NULL;

const int num_species = 4;
std::string dirs[num_species] = {"../ep_10_100/outfiles/", "../eD_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/"};
std::string names[num_species] = {"e + p @ 10 + 100 GeV", "e + D @ 18 + 110 GeV", "e + C @ 18 + 110 GeV", "e + Au @ 18 + 110 GeV"};
std::string outdirs[num_species] = {"../ep_10_100/", "../eD_18_110/", "../eC_18_110/", "../eAu_18_110/"};
std::string inFileName = "pt_eta_decomp.root";

const int num_histograms = 6;
std::string histogram_titles[num_histograms] = {"kaon_pos", "kaon_neg", "pion_pos", "pion_neg", "proton_pos", "proton_neg"};
std::string histogram_symbols[num_histograms] = {"K^{+}", "K^{-}", "#pi^{+}", "#pi^{-}", "proton", "antiproton"};


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
  h2d_pt_vs_eta->GetYaxis()->SetRangeUser(-7,7);
}

void plot_histograms_pt_eta_decomp()
{
  for(int i_titles; i_titles < num_histograms; i_titles++) {

    TFile * fin;
    std::string inFile;
    for (int i = 0; i < num_species; i++) {
      inFile = dirs[i] + inFileName;
      fin = new TFile(inFile.c_str(), "read");
      h2d_pt_vs_eta = (TH2D*)fin->Get(histogram_titles[i_titles].c_str());
      Double_t nEntries = h2d_pt_vs_eta->GetEntries();    // make sure the eta range on the h2d_pt_vs_eta is large enough so that all particles are included
      slice_2D_hist();

      // plots pt and eta 2d histogram
      TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
      c3->Range(0,0,1,1);
      c3->SetLeftMargin(0.15);
      c3->SetRightMargin(0.15);
      c3->SetBottomMargin(0.1);
      h2d_pt_vs_eta = (TH2D*) h2d_pt_vs_eta->Rebin2D(2, 2);
      h2d_pt_vs_eta->SetTitle(Form("%s multiplicity, %s", histogram_symbols[i_titles].c_str(), names[i].c_str()));
      h2d_pt_vs_eta->GetXaxis()->SetTitle("<p_{T}> [GeV/c]");
      h2d_pt_vs_eta->GetYaxis()->SetTitle("#eta");
      h2d_pt_vs_eta->GetXaxis()->SetTitleOffset(1.3);
      h2d_pt_vs_eta->GetYaxis()->SetTitleOffset(1.5);
      c3->SetLogz();
      h2d_pt_vs_eta->Draw("colz");
      h2d_pt_vs_eta->SetStats(0);
      c3->SaveAs(Form("%spt_vs_eta_2D_histo_%s.pdf", outdirs[i].c_str(), histogram_titles[i_titles].c_str()));

      // plot different eta bin slices on the same canvas
      c3 = new TCanvas("c3","c3",800,800); // create new canvas
      c3->Range(0,0,1,1);
      c3->SetLeftMargin(0.15);
      c3->SetBottomMargin(0.1);
      //c3->SetLogy(); // set y axis to log scale
      TLegend* leg = new TLegend(0.65,0.60,0.75,0.80);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.035);
      leg->SetFillStyle(0);
      leg->SetMargin(0.3);
      float plot_xrange_lo = 0, plot_xrange_hi = 7;
      float plot_yrange_lo = 0, plot_yrange_hi = 0.7; // when using log axis, cannot use 0 as start as plot range
      // use the empty 2D histogram htemp as a frame
      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.SetStats(0); // not showing the box on the top right corner
      htemp.Draw();
      htemp.SetTitle(Form("%s multiplicity", histogram_symbols[i_titles].c_str()));
      htemp.GetXaxis()->SetTitle("<p_{T}> [GeV/c]");
      htemp.GetYaxis()->SetTitle("fraction of counts [%]");
      htemp.GetXaxis()->SetTitleOffset(1.3);
      htemp.GetYaxis()->SetTitleOffset(1.5);
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        h1d_pt_in_eta[ieta]->SetStats(0); // not showing the box on the top right corner
        h1d_pt_in_eta[ieta]->SetLineColor(eta_color[ieta]); // eta_color array defined at the beginning
        h1d_pt_in_eta[ieta]->SetMarkerColor(eta_color[ieta]);
        h1d_pt_in_eta[ieta]->Scale(1 / nEntries);
        h1d_pt_in_eta[ieta]->Rebin(2);
        //cout<<h1d_pt_in_eta[ieta]->Integral(0,10)<<endl; // integral should be 1 if all particles included in eta range!
        h1d_pt_in_eta[ieta]->Draw("hsame");
        leg->AddEntry(h1d_pt_in_eta[ieta],Form("%.0f < #eta < %.0f",eta_lo[ieta],eta_hi[ieta]),"l"); // "l" means line
      }
      leg->Draw("same");
      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.5,0.85,names[i].c_str());
      c3->SaveAs( Form("%spt_eta_binned_single_%s.pdf", outdirs[i].c_str(), histogram_titles[i_titles].c_str()) );

      /*
      // all in one pdf
      TCanvas * c_all = new TCanvas("canvas_with_all_of_them", "canvas_with_all_of_them", 2400, 800);
      c_all->Divide(3);
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        c_all->cd(ieta+1); h1d_pt_in_eta[ieta]->Draw("hsame");
        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.035);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.2,0.85,names[i].c_str());
        tl->DrawLatexNDC(0.2,0.80,Form("%.0f < #eta < %.0f",eta_lo[ieta],eta_hi[ieta]));
        h1d_pt_in_eta[ieta]->GetXaxis()->SetTitle("<p_{T}> [GeV/c]");
        h1d_pt_in_eta[ieta]->GetYaxis()->SetTitle("counts");
        h1d_pt_in_eta[ieta]->GetXaxis()->SetTitleOffset(1.3);
        h1d_pt_in_eta[ieta]->GetYaxis()->SetTitleOffset(1.5);
        h1d_pt_in_eta[ieta]->GetXaxis()->SetRangeUser(0,10);
        h1d_pt_in_eta[ieta]->SetTitle("");
        h1d_pt_in_eta[ieta]->SetStats(0); // not showing the box on the top right corner
        h1d_pt_in_eta[ieta]->SetLineColor(eta_color[ieta]); // eta_color array defined at the beginning
        h1d_pt_in_eta[ieta]->SetMarkerColor(eta_color[ieta]);
      }
      c_all->SaveAs(Form("%spt_eta_binned_decomp_%s.pdf", outdirs[i].c_str(), histogram_titles[i_titles].c_str()));

      // plot average pt as function of eta using TGraph
      double temp_pt[etabin] = {0}, temp_pt_err[etabin] = {0};
      double temp_eta[etabin] = {0}, temp_eta_err[etabin] = {0};
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        temp_pt[ieta] = h1d_pt_in_eta[ieta]->GetMean();
        temp_eta[ieta] = 0.5*(eta_lo[ieta]+eta_hi[ieta]);
      }
      g_pt_vs_eta = new TGraphErrors(etabin,temp_eta,temp_pt,temp_eta_err,temp_pt_err);
      g_pt_vs_eta->SetMarkerColor(kRed);
      g_pt_vs_eta->SetMarkerStyle(20);
      g_pt_vs_eta->SetMarkerSize(1);
      TCanvas* c4 = new TCanvas("c4","c4",800,800); // create new canvas
      c4->Range(0,0,1,1);
      c4->SetLeftMargin(0.15);
      c4->SetBottomMargin(0.1);
      leg = new TLegend(0.60,0.70,0.80,0.80);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.035);
      leg->SetFillStyle(0);
      leg->SetMargin(0.3);
      plot_xrange_lo = eta_lo[0]-0.5, plot_xrange_hi = eta_hi[etabin-1]+0.5;
      plot_yrange_lo = 0, plot_yrange_hi = 4;
      // use the empty 2D histogram htemp as a frame
      TH2F htemp1("htemp1","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp1.SetStats(0); // not showing the box on the top right corner
      htemp1.Draw();
      htemp1.GetYaxis()->SetTitle("average <p_{T}> [GeV/c]");
      htemp1.GetXaxis()->SetTitle("#eta");
      htemp1.GetXaxis()->SetTitleOffset(1.3);
      htemp1.GetYaxis()->SetTitleOffset(1.5);
      g_pt_vs_eta->Draw("psame");
      tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.035);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.2,0.85,names[i].c_str());
      c4->SaveAs( Form("%savg_pt_vs_eta_%s.pdf", outdirs[i].c_str(), histogram_titles[i_titles].c_str()) );
      */
    }

  }

}
