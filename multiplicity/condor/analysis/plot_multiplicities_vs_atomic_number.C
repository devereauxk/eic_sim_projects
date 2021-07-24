void plot_multiplicities_vs_atomic_number()
{

  const int num_species = 7;
  std::string dirs[num_species] = {"../ep_10_100/output.root", "../eD_18_110/output.root", "../eHe4_18_110/output.root", "../eC_18_110/output.root", "../eCa_18_110/output.root", "../eCu_18_110/output.root", "../eAu_18_110/output.root"};
  double atomic_weight[num_species] = {1, 2, 4, 12, 40, 63, 197};

  double kaon_pos_mean[num_species] = {0}, kaon_pos_mean_err[num_species] = {0};
  double kaon_neg_mean[num_species] = {0}, kaon_neg_mean_err[num_species] = {0};
  double kaon_pos_rms[num_species] = {0}, kaon_pos_rms_err[num_species] = {0};
  double kaon_neg_rms[num_species] = {0}, kaon_neg_rms_err[num_species] = {0};

  double pion_pos_mean[num_species] = {0}, pion_pos_mean_err[num_species] = {0};
  double pion_neg_mean[num_species] = {0}, pion_neg_mean_err[num_species] = {0};
  double pion_pos_rms[num_species] = {0}, pion_pos_rms_err[num_species] = {0};
  double pion_neg_rms[num_species] = {0}, pion_neg_rms_err[num_species] = {0};

  double proton_pos_mean[num_species] = {0}, proton_pos_mean_err[num_species] = {0};
  double proton_neg_mean[num_species] = {0}, proton_neg_mean_err[num_species] = {0};
  double proton_pos_rms[num_species] = {0}, proton_pos_rms_err[num_species] = {0};
  double proton_neg_rms[num_species] = {0}, proton_neg_rms_err[num_species] = {0};

  const char * dir_char;
  TFile * fin;
  TH2D * h2d_kaon, *h2d_pion, *h2d_proton;
  TH1D * h1d_kaon_pos, *h1d_kaon_neg, *h1d_pion_pos, *h1d_pion_neg, *h1d_proton_pos, *h1d_proton_neg;
  for (int i = 0; i < num_species; i++) {
    cout<<"processing " + dirs[i] + " ..."<<endl;

    dir_char = dirs[i].c_str();
    fin = new TFile(dir_char, "read");

    h2d_kaon = (TH2D*)fin->Get("h2d_kaon");
    h1d_kaon_pos = (TH1D*) h2d_kaon->ProjectionX("h1d_kaon_pos");
    h1d_kaon_neg = (TH1D*) h2d_kaon->ProjectionY("h1d_kaon_neg");

    h2d_pion = (TH2D*)fin->Get("h2d_pion");
    h1d_pion_pos = (TH1D*) h2d_pion->ProjectionX("h1d_pion_pos");
    h1d_pion_neg = (TH1D*) h2d_pion->ProjectionY("h1d_pion_neg");

    h2d_proton = (TH2D*)fin->Get("h2d_proton");
    h1d_proton_pos = (TH1D*) h2d_proton->ProjectionX("h1d_proton");
    h1d_proton_neg = (TH1D*) h2d_proton->ProjectionY("h1d_anti_proton");

    kaon_pos_mean[i] = h1d_kaon_pos->GetMean(); kaon_pos_mean_err[i] = h1d_kaon_pos->GetMeanError();
    kaon_neg_mean[i] = h1d_kaon_neg->GetMean(); kaon_neg_mean_err[i] = h1d_kaon_neg->GetMeanError();
    kaon_pos_rms[i] = h1d_kaon_pos->GetRMS(); kaon_pos_rms_err[i] = h1d_kaon_pos->GetRMSError();
    kaon_neg_rms[i] = h1d_kaon_neg->GetRMS(); kaon_neg_rms_err[i] = h1d_kaon_neg->GetRMSError();

    pion_pos_mean[i] = h1d_pion_pos->GetMean(); pion_pos_mean_err[i] = h1d_pion_pos->GetMeanError();
    pion_neg_mean[i] = h1d_pion_neg->GetMean(); pion_neg_mean_err[i] = h1d_pion_neg->GetMeanError();
    pion_pos_rms[i] = h1d_pion_pos->GetRMS(); pion_pos_rms_err[i] = h1d_pion_pos->GetRMSError();
    pion_neg_rms[i] = h1d_pion_neg->GetRMS(); pion_neg_rms_err[i] = h1d_pion_neg->GetRMSError();

    proton_pos_mean[i] = h1d_proton_pos->GetMean(); proton_pos_mean_err[i] = h1d_proton_pos->GetMeanError();
    proton_neg_mean[i] = h1d_proton_neg->GetMean(); proton_neg_mean_err[i] = h1d_proton_neg->GetMeanError();
    proton_pos_rms[i] = h1d_proton_pos->GetRMS(); proton_pos_rms_err[i] = h1d_proton_pos->GetRMSError();
    proton_neg_rms[i] = h1d_proton_neg->GetRMS(); proton_neg_rms_err[i] = h1d_proton_neg->GetRMSError();

    fin->Close();

  }

  TGraphErrors* kaon_pos_mean_vs_an = new TGraphErrors(num_species, atomic_weight, kaon_pos_mean, NULL, kaon_pos_mean_err);
  TGraphErrors* kaon_neg_mean_vs_an = new TGraphErrors(num_species, atomic_weight, kaon_neg_mean, NULL, kaon_neg_mean_err);
  TGraphErrors* kaon_pos_rms_vs_an = new TGraphErrors(num_species, atomic_weight, kaon_pos_rms, NULL, kaon_pos_rms_err);
  TGraphErrors* kaon_neg_rms_vs_an = new TGraphErrors(num_species, atomic_weight, kaon_neg_rms, NULL, kaon_neg_rms_err);
  TGraphErrors* pion_pos_mean_vs_an = new TGraphErrors(num_species, atomic_weight, pion_pos_mean, NULL, pion_pos_mean_err);
  TGraphErrors* pion_neg_mean_vs_an = new TGraphErrors(num_species, atomic_weight, pion_neg_mean, NULL, pion_neg_mean_err);
  TGraphErrors* pion_pos_rms_vs_an = new TGraphErrors(num_species, atomic_weight, pion_pos_rms, NULL, pion_pos_rms_err);
  TGraphErrors* pion_neg_rms_vs_an = new TGraphErrors(num_species, atomic_weight, pion_neg_rms, NULL, pion_neg_rms_err);
  TGraphErrors* proton_pos_mean_vs_an = new TGraphErrors(num_species, atomic_weight, proton_pos_mean, NULL, proton_pos_mean_err);
  TGraphErrors* proton_neg_mean_vs_an = new TGraphErrors(num_species, atomic_weight, proton_neg_mean, NULL, proton_neg_mean_err);
  TGraphErrors* proton_pos_rms_vs_an = new TGraphErrors(num_species, atomic_weight, proton_pos_rms, NULL, proton_pos_rms_err);
  TGraphErrors* proton_neg_rms_vs_an = new TGraphErrors(num_species, atomic_weight, proton_neg_rms, NULL, proton_neg_rms_err);

  TGraphErrors* graph_arr[12] = {kaon_pos_mean_vs_an, kaon_neg_mean_vs_an, kaon_pos_rms_vs_an, kaon_neg_rms_vs_an, pion_pos_mean_vs_an, pion_neg_mean_vs_an, pion_pos_rms_vs_an, pion_neg_rms_vs_an, proton_pos_mean_vs_an, proton_neg_mean_vs_an, proton_pos_rms_vs_an, proton_neg_rms_vs_an};

  TCanvas* c_main;
  TGraphErrors * pos_graph;
  TGraphErrors * neg_graph;
  TLegend* leg;
  TLatex* t1;
  float plot_xrange_lo = 0, plot_xrange_hi = 10;
  float plot_yrange_lo = 0, plot_yrange_hi = 10;

  // kaon --------------------------------------------------------------------------------------------------
  c_main = new TCanvas("c_main", "c_main", 800, 800);
  c_main->Range(0,0,1,1);
  c_main->SetLeftMargin(0.15);
  c_main->SetBottomMargin(0.1);
  c_main->Divide(1,2);

  c_main->cd(1);
  pos_graph = graph_arr[0];
  pos_graph->SetMarkerColor(kRed);
  pos_graph->SetMarkerStyle(20);
  pos_graph->SetMarkerSize(1);

  neg_graph = graph_arr[1];
  neg_graph->SetMarkerColor(kBlue);
  neg_graph->SetMarkerStyle(20);
  neg_graph->SetMarkerSize(1);

  plot_xrange_lo = 0;
  plot_xrange_hi = 220;
  plot_yrange_lo = min(TMath::MinElement(pos_graph->GetN(),pos_graph->GetY()), TMath::MinElement(neg_graph->GetN(),neg_graph->GetY())) * 0.92;
  plot_yrange_hi = max(TMath::MaxElement(pos_graph->GetN(),pos_graph->GetY()), TMath::MaxElement(neg_graph->GetN(),neg_graph->GetY())) * 1.08;

  TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
  htemp.SetStats(0);
  htemp.Draw();
  htemp.GetXaxis()->SetTitle("mass number of A");
  htemp.GetYaxis()->SetTitle("average particle multiplicity [counts]");
  htemp.GetXaxis()->SetTitleOffset(1.3);
  htemp.GetYaxis()->SetTitleOffset(1.5);

  pos_graph->Draw("psame");
  neg_graph->Draw("psame");

  leg = new TLegend(0.60, 0.70, 0.80, 0.80);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillStyle(0);
  leg->SetMargin(0.3);
  leg->AddEntry(pos_graph, "K^+", "l");
  leg->AddEntry(neg_graph, "K^-", "l");

  t1 = new TLatex();
  t1->SetTextAlign(11);
  t1->SetTextSize(0.045);
  t1->SetTextColor(kBlack);
  t1->DrawLatexNDC(0.4,0.85,"e + p @ 10 + 100 GeV");
  t1->DrawLatexNDC(0.4,0.80,"e + A @ 18 + 110 Gev");

  c_main->cd(2);
  pos_graph = graph_arr[2];
  pos_graph->SetMarkerColor(kRed);
  pos_graph->SetMarkerStyle(20);
  pos_graph->SetMarkerSize(1);

  neg_graph = graph_arr[3];
  neg_graph->SetMarkerColor(kBlue);
  neg_graph->SetMarkerStyle(20);
  neg_graph->SetMarkerSize(1);

  plot_xrange_lo = 0;
  plot_xrange_hi = 220;
  plot_yrange_lo = min(TMath::MinElement(pos_graph->GetN(),pos_graph->GetY()), TMath::MinElement(neg_graph->GetN(),neg_graph->GetY())) * 0.92;
  plot_yrange_hi = max(TMath::MaxElement(pos_graph->GetN(),pos_graph->GetY()), TMath::MaxElement(neg_graph->GetN(),neg_graph->GetY())) * 1.08;

  TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
  htemp.SetStats(0);
  htemp.Draw();
  htemp.GetXaxis()->SetTitle("mass number of A");
  htemp.GetYaxis()->SetTitle("particle multiplicity rms [counts]");
  htemp.GetXaxis()->SetTitleOffset(1.3);
  htemp.GetYaxis()->SetTitleOffset(1.5);

  pos_graph->Draw("psame");
  neg_graph->Draw("psame");

  c_main->SaveAs( Form("kaon_vs_atomic_number.pdf") );

}
