const int xbin = 5;
const double x_lo[xbin] = {1e-5, 1e-4 1e-3, 1e-2, 1e-1};
const double x_hi[xbin] = {1e-4, 1e-3, 1e-2, 1e-1, 1};
const int x_color[xbin] = {kRed, kBlue, kGreen+1, kBlack, kOrange};
TH1D * kaon[xbin];
TH1D * pion[xbin];
TH1D * proton[xbin];

const int num_species = 7;
std::string dirs[num_species] = {"../ep_10_100/outfiles/", "../eD_18_110/outForPythiaMode/",  "../eHe4_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eCa_18_110/outForPythiaMode/", "../eCu_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/"};
std::string names[num_species] = {"e + p @ 10 + 100 GeV", "e + D @ 18 + 110 GeV", "e + He4 @ 18 + 110 GeV", "e + C @ 18 + 110 GeV", "e + Ca @ 18 + 110 GeV", "e + Cu @ 18 + 110 GeV", "e + Au @ 18 + 110 GeV"};
std::string outdirs[num_species] = {"../ep_10_100/", "../eD_18_110/", "../eHe4_18_110/", "../eC_18_110/", "../eCa_18_110/", "../eCu_18_110/", "../eAu_18_110/"};
std::string inFileName = "mult_x_binned.root";

void plot_histograms_mult_x_binned()
{
  TFile * fin;
  std::string inFile;
  for (int i = 0; i < num_species; i++) {
    for (int ix = 0; ix < xbin; ix++) {

      inFile = dirs[i] + inFileName;
      fin = new TFile(inFile.c_str(),"read");

      TH2D* h2d_kaon = (TH2D*)fin->Get(Form("h2d_kaon_%d", ix));
      TH1D* h1d_kaon_pos = (TH1D*) h2d_kaon->ProjectionX("h1d_kaon_pos");
      TH1D* h1d_kaon_neg = (TH1D*) h2d_kaon->ProjectionY("h1d_kaon_neg");
      TH1D* h1d_kaon_total = (TH1D*) h1d_kaon_pos->Clone();
      h1d_kaon_total->Add(h1d_kaon_neg);

      TH2D* h2d_pion = (TH2D*)fin->Get(Form("h2d_pion_%d", ix));
      TH1D* h1d_pion_pos = (TH1D*) h2d_pion->ProjectionX("h1d_pion_pos");
      TH1D* h1d_pion_neg = (TH1D*) h2d_pion->ProjectionY("h1d_pion_neg");
      TH1D* h1d_pion_total = (TH1D*) h1d_pion_pos->Clone();
      h1d_pion_total->Add(h1d_pion_neg);

      TH2D* h2d_proton = (TH2D*)fin->Get(Form("h2d_proton_%d", ix));
      TH1D* h1d_proton = (TH1D*) h2d_proton->ProjectionX("h1d_proton");
      TH1D* h1d_anti_proton = (TH1D*) h2d_proton->ProjectionY("h1d_anti_proton");
      TH1D* h1d_proton_total = (TH1D*) h1d_proton->Clone();
      h1d_proton_total->Add(h1d_anti_proton);

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
      tl->DrawLatexNDC(0.4,0.85, Form("%e < x < %e",x_lo[ix],x_hi[ix]));
      tl->DrawLatexNDC(0.4,0.80, names[i].c_str());
      tl->DrawLatexNDC(0.4,0.75,Form("%e events", h1d_proton_total->GetEntries() / 2));
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
      tl->DrawLatexNDC(0.4,0.85, Form("%e < x < %e",x_lo[ix],x_hi[ix]));
      tl->DrawLatexNDC(0.4,0.80, names[i].c_str());
      tl->DrawLatexNDC(0.4,0.75,Form("%e events", h1d_kaon_total->GetEntries() / 2));
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
      tl->DrawLatexNDC(0.4,0.85, Form("%e < x < %e",x_lo[ix],x_hi[ix]));
      tl->DrawLatexNDC(0.4,0.80, names[i].c_str());
      tl->DrawLatexNDC(0.4,0.75,Form("%e events", h1d_pion_total->GetEntries() / 2));
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
      c_all->SaveAs(Form("%smultiplicities_x_binned_%d.pdf", outdirs[i].c_str(), ix));

    }
  }

}
