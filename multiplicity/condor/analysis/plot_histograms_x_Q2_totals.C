TH2D* h2d_x_vs_q2 = NULL;

const int num_species = 7;
std::string dirs[num_species] = {"../ep_10_100/outfiles/", "../eD_18_110/outForPythiaMode/",  "../eHe4_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eCa_18_110/outForPythiaMode/", "../eCu_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/"};
std::string names[num_species] = {"e + p @ 10 + 100 GeV", "e + D @ 18 + 110 GeV", "e + He4 @ 18 + 110 GeV", "e + C @ 18 + 110 GeV", "e + Ca @ 18 + 110 GeV", "e + Cu @ 18 + 110 GeV", "e + Au @ 18 + 110 GeV"};
std::string outdirs[num_species] = {"../ep_10_100/", "../eD_18_110/", "../eHe4_18_110/", "../eC_18_110/", "../eCa_18_110/", "../eCu_18_110/", "../eAu_18_110/"};
std::string inFileName = "x_Q2_total_TH2D.root";


void plot_histograms_x_Q2_totals()
{
  TFile * fin;
  std::string inFile;
  for (int i = 0; i < num_species; i++) {
    inFile = dirs[i] + inFileName;
    fin = new TFile(inFile.c_str(), "read");
    fin->ls();
    h2d_x_vs_q2 = (TH2D*)fin->Get("h2d_x_vs_q2");

    // plots pt and eta 2d histogram
    TCanvas* c3 = new TCanvas("c3","c3",800,800); // create new canvas
    c3->Range(0,0,1,1);
    c3->SetLeftMargin(0.15);
    c3->SetRightMargin(0.15);
    c3->SetBottomMargin(0.1);
    h2d_x_vs_q2 = (TH2D*) h2d_x_vs_q2->Rebin2D(2, 2);
    h2d_x_vs_q2->SetTitle(names[i].c_str());
    h2d_x_vs_q2->GetXaxis()->SetTitle("log(x)");
    h2d_x_vs_q2->GetYaxis()->SetTitle("Q^{2}");
    h2d_x_vs_q2->GetXaxis()->SetTitleOffset(1.3);
    h2d_x_vs_q2->GetYaxis()->SetTitleOffset(1.5);
    c3->SetLogy();
    c3->SetLogz();
    h2d_x_vs_q2->Draw("colz");
    h2d_x_vs_q2->SetStats(0);
    c3->SaveAs(Form("%sx_vs_Q2_2D_histo.pdf", outdirs[i].c_str()));

  }

}
