R__LOAD_LIBRARY(libeicsmear);
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 40, 5};
static double pt_hi[ptbin] = {10, 20, 40, 60, 60};
const int pt_color[ptbin] = {kGreen+1, kBlue, kViolet, kOrange+1, kRed};

const int etabin = 6; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5, -1, 0};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5, 0, 1};
const int eta_color[etabin] = {kGreen+1, kBlue, kViolet, kOrange+1};

// jet level histograms
TH3D* h3d_jet_eec_rl_xi_phi[etabin][ptbin] = {};
TH3D* h3d_jet_eec_rlsqrtpt_xi_phi[etabin][ptbin] = {};
TH1D* h1d_jet_pt[etabin] = {};
TH1D* h1d_jet_eta = NULL;
TH1D* h1d_jet_multiplicity[etabin][ptbin] = {};
TH1D* h1d_jet_multiplicity_charged[etabin][ptbin] = {};
TH2D* h2d_Q2_x[etabin][ptbin] = {};

// particle level histograms
TH1D* h1d_part_pt[etabin] = {};
TH1D* h1d_part_eta[ptbin] = {};
TH1D* h1d_part_mult = NULL;

static int cno = 0;

void hists_to_csv(const char* outfile_name, vector<TH1*> hists)
{
   ofstream outfile;
   outfile.open(Form("%s", outfile_name));

   int n = hists[0]->GetNbinsX();

   for (int i = 1; i <= n; i++) // loop over lines
   {
     outfile << hists[0]->GetBinCenter(i) << "," << hists[0]->GetBinWidth(i);

     for (int ihist = 0; ihist < hists.size(); ihist++)
     {
       outfile << "," << hists[ihist]->GetBinContent(i) << "," << hists[ihist]->GetBinError(i);
     }
     outfile << "\n";

   }
   outfile.close();
}

void hists_to_csv_2d(const char* outfile_name, vector<TH2*> hists)
{
  ofstream outfile;
  outfile.open(Form("%s", outfile_name));

  int nx = hists[0]->GetNbinsX();
  int ny = hists[0]->GetNbinsY();

  for (int i = 1; i <= nx; i++) // loop over lines, x axis loop
  {
    for (int j = 1; j <= ny; j++) // y axis loop
    {
      outfile << hists[0]->GetXaxis()->GetBinCenter(i) << "," << hists[0]->GetXaxis()->GetBinWidth(i) << ",";
      outfile << hists[0]->GetYaxis()->GetBinCenter(j) << "," << hists[0]->GetYaxis()->GetBinWidth(j);

      for (int ihist = 0; ihist < hists.size(); ihist++)
      {
        outfile << "," << hists[ihist]->GetBinContent(i,j) << "," << hists[ihist]->GetBinError(i,j);
      }
      outfile << "\n";
    }

  }
  outfile.close();
}

void hists_to_csv_3d(const char* outfile_name, vector<TH3*> hists)
{
  ofstream outfile;
  outfile.open(Form("%s", outfile_name));

  int nx = hists[0]->GetNbinsX();
  int ny = hists[0]->GetNbinsY();
  int nz = hists[0]->GetNbinsZ();

  for (int i = 1; i <= nx; i++) // loop over lines, x axis loop
  {
    for (int j = 1; j <= ny; j++) // y axis loop
    {
      for (int k = 1; k <= nz; k++)
      {
        outfile << hists[0]->GetXaxis()->GetBinCenter(i) << "," << hists[0]->GetXaxis()->GetBinWidth(i) << ",";
        outfile << hists[0]->GetYaxis()->GetBinCenter(j) << "," << hists[0]->GetYaxis()->GetBinWidth(j) << ",";
        outfile << hists[0]->GetZaxis()->GetBinCenter(k) << "," << hists[0]->GetZaxis()->GetBinWidth(k);

        for (int ihist = 0; ihist < hists.size(); ihist++)
        {
          outfile << "," << hists[ihist]->GetBinContent(i,j,k) << "," << hists[ihist]->GetBinError(i,j,k);
        }
        outfile << "\n";
      }
    }
  }
  
  outfile.close();
}

void individual_hists(const char* out_dir)
{
  TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);
  leg->SetFillStyle(0);
  leg->SetMargin(0.1);

  // 1d jet pt histogram
  mclogy(cno++);
  {
    for (int ieta = 0; ieta < etabin; ieta++)
    {
      h1d_jet_pt[ieta]->Draw("same");

      h1d_jet_pt[ieta]->GetXaxis()->SetRangeUser(0,70);
      h1d_jet_pt[ieta]->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      h1d_jet_pt[ieta]->GetYaxis()->SetTitle("counts");
      h1d_jet_pt[ieta]->GetXaxis()->SetTitleOffset(1.3);
      h1d_jet_pt[ieta]->GetYaxis()->SetTitleOffset(1.5);
      h1d_jet_pt[ieta]->Draw("same hist e");
      leg->AddEntry(h1d_jet_pt[ieta],Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt.pdf\")", cno-1, out_dir) );
  }

  // 1d jet pt histogram normalized by number of jets, inclusive on eta
  mclogy(cno++);
  {
    TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      TH1D* temp = (TH1D*) h1d_jet_pt[ieta]->Clone("temp");
      temp->Scale(1/temp->GetEntries());

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(0,70);
      temp->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_pt_selfnorm.pdf\")", cno-1, out_dir) );
  }

  // 1d jet eta histogram
  mcs(cno++);
  {
    float plot_xrange_lo = -5;
    float plot_xrange_hi = 5;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = h1d_jet_eta->GetMaximum()*1.15;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw("hsame");
    htemp.GetXaxis()->SetTitle("jet eta");
    htemp.GetYaxis()->SetTitle("counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    h1d_jet_eta->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_eta.pdf\")", cno-1, out_dir) );
  }

}

void overlay_hists(const char* out_dir)
{
  // overlay h1d_jet_eec_rl_xi_phi projected to 1D h1d_jet_rl, for each eta bin curves for pt
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mclogxy(cno++);
    {
      float plot_xrange_lo = 1E-2;
      float plot_xrange_hi = 1;
      float plot_yrange_lo = 1E-5;
      float plot_yrange_hi = 5E-1;

      TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.025);
      leg->SetFillStyle(0);
      leg->SetMargin(0.1);

      TH1D* temp;

      for (int ipt = 0; ipt < ptbin-1; ipt++)
      {
        temp = (TH1D*) h3d_jet_eec_rl_xi_phi[ieta][ipt]->ProjectionX();
        temp->Scale(1/temp->Integral());

        temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        //temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->SetMarkerColor(pt_color[ipt]);
        temp->SetLineColor(pt_color[ipt]);
        temp->GetXaxis()->SetTitle("R_{L}");
        temp->GetYaxis()->SetTitle("normalized E3C");
        temp->SetMarkerSize(0.5);
        temp->SetMarkerStyle(21);
        temp->Draw("same hist");
        leg->AddEntry(temp,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
      }
      leg->Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.025);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_e3c_rl_%d.pdf\")", cno-1, out_dir, ieta) );
    }
  }

  // overlay h1d_jet_eec_rlsqrtpt_xi_phi projected to 1D h1d_jet_rlsqrtpt, for each eta bin curves for pt
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    mclogxy(cno++);
    {
      float plot_xrange_lo = 1E-2*sqrt(20);
      float plot_xrange_hi = 1*sqrt(20);
      float plot_yrange_lo = 1E-5;
      float plot_yrange_hi = 5E-1;

      TLegend* leg = new TLegend(0.21,0.7,0.51,0.82);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.025);
      leg->SetFillStyle(0);
      leg->SetMargin(0.1);

      TH1D* temp;

      for (int ipt = 0; ipt < ptbin-2; ipt++)
      {
        temp = (TH1D*) h3d_jet_eec_rlsqrtpt_xi_phi[ieta][ipt]->ProjectionX();
        temp->Scale(1/temp->Integral());

        temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        //temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->SetMarkerColor(pt_color[ipt]);
        temp->SetLineColor(pt_color[ipt]);
        temp->GetXaxis()->SetTitle("R_{L}#sqrt{p_{T}}");
        temp->GetYaxis()->SetTitle("normalized E3C");
        temp->SetMarkerSize(0.5);
        temp->SetMarkerStyle(21);
        temp->Draw("hist same");
        leg->AddEntry(temp,Form("%.1f GeV < p_{T} < %.1f GeV",pt_lo[ipt],pt_hi[ipt]));
      }
      leg->Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.025);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_e3c_rlsqrtpt_%d.pdf\")", cno-1, out_dir, ieta) );
    }
  }

}

void raw_nice_3d_plots(const char* out_dir)
{
  // xi-phi plots for rl bins within range, forward eta, 20-30GeV pt bin
  int etabin_pick = 2;
  int ptbin_pick = 1; // changed to 3 for eU, use 2 for all other cases
  float rl_range_lo = 1E-2;
  float rl_range_hi = 1;
  TH3D* picked = h3d_jet_eec_rl_xi_phi[etabin_pick][ptbin_pick];
  int norm_binrange_lo = picked->GetXaxis()->FindBin(rl_range_lo);
  cout<<"low bin: "<<norm_binrange_lo<<endl;
  int norm_binrange_hi = picked->GetXaxis()->FindBin(rl_range_hi);
  cout<<"high bin: "<<norm_binrange_hi<<endl;
  if (norm_binrange_lo == 0)
  {
    norm_binrange_lo = 1;
    cout<<"bin range lo too low; set to 1"<<endl;
  }
  if (norm_binrange_hi > picked->GetNbinsX()+1)
  {
    norm_binrange_hi = picked->GetNbinsX()+1;
    cout<<"bin range hi too high; set to "<<picked->GetNbinsX()<<endl;
  }

  vector<TH3*> hists;
  hists.push_back(picked);
  hists_to_csv_3d(Form("%s3d_hists.csv", out_dir), hists);

  // determine zrange for plotting, looking at RL=40th bin and RL=50th bin...
  float bin_center = picked->GetXaxis()->GetBinCenter(40);
  sliced = (TH3D*) picked->Clone("temp3d");
  sliced->GetXaxis()->SetRange(ibin,ibin);
  temp = (TH2D*) sliced->Project3D("zy");
  float plot_zrange_lo = temp->GetMinimum();

  float bin_center = picked->GetXaxis()->GetBinCenter(50);
  sliced = (TH3D*) picked->Clone("temp3d");
  sliced->GetXaxis()->SetRange(ibin,ibin);
  temp = (TH2D*) sliced->Project3D("zy");
  float plot_zrange_hi = temp->GetMaximum();

  TH3D* sliced;
  TH2D* temp;
  for (int ibin = norm_binrange_lo; ibin <= norm_binrange_hi; ibin++)
  {
    // raw distribution 2d projection
    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 1;
      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1.5; // TMath::Pi() / 2.0;

      float bin_center = picked->GetXaxis()->GetBinCenter(ibin);
      sliced = (TH3D*) picked->Clone("temp3d");
      sliced->GetXaxis()->SetRange(ibin,ibin);
      temp = (TH2D*) sliced->Project3D("zy");

      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetXaxis()->SetTitle("#xi");
      temp->GetYaxis()->SetTitle("#phi");
      temp->Draw("colz");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.03);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %.1f), p_{T} #in [%.1f, %.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick],pt_lo[ptbin_pick],pt_hi[ptbin_pick]));
      tl->DrawLatexNDC(0.22,0.81,Form("R_{L} ~ %.3f", bin_center));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_eec_xi_phi_%d.pdf\")", cno-1, out_dir, ibin) );
    }

    // raw distribution nice 3d plot
    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13, -30);
    {
      temp->Draw("SURF2Z");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d_surface.pdf\")", cno-1, out_dir, ibin) );
    }

    // raw distribution nice 3d plot, fixed z-axis
    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13, -30);
    {
      temp->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);

      temp->Draw("SURF2Z");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d_surface_zfixed.pdf\")", cno-1, out_dir, ibin) );
    }


    // "modification of uniform distro wrt raw distro"
    // calculating essentially (1 - raw) / raw 
    """
    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13, -30);
    {
      TH2D* temp_for_calc = (TH2D*) temp->Clone("temp temp");
      Int_t nBinsX = temp_for_calc->GetNbinsX();
      Int_t nBinsY = temp_for_calc->GetNbinsY();
      for (Int_t iBinX = 1; iBinX <= nBinsX; ++iBinX) {
        for (Int_t iBinY = 1; iBinY <= nBinsY; ++iBinY) {
          double original = temp->GetBinContent(iBinX, iBinY);
          temp_for_calc->SetBinContent(iBinX, iBinY, 1 / original);
        }
      }

      temp_for_calc->Draw("SURF2Z");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3cModUniform_xi_phi_%d_surface.pdf\")", cno-1, out_dir, ibin) );
    }
    """

  }
}

void Q2_x_panel(const char* out_dir)
{

  float plot_xrange_lo = 0.8;
  float plot_xrange_hi = 2.7;
  float plot_yrange_lo = -0.005;
  float plot_yrange_hi = 0.04;
  float plot_zrange_lo = 0;
  float plot_zrange_hi = 750E3;

  for (int ieta = 0; ieta < 3; ieta++)
  {
    for (int ipt = 0; ipt < 3; ipt++)
    {
      mclogxyz(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
      {
        TH2D* temp = (TH2D*) h2d_Q2_x[ieta][ipt]->Clone();

        // plot
        //temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        //temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
        temp->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);
        temp->GetXaxis()->SetTitle("x_{B}");
        temp->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
        temp->Draw("colz");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.028);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.22,0.84,"eHIJING, e+p, 4*10^{8} events");
        tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_Q2_x_%i_%i.pdf\")", cno-1, out_dir, ieta, ipt) );
      }
    }
  }

}

void particle_hists(const char* out_dir)
{
  // 1d particle pt histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TLegend* leg = new TLegend(0.51,0.7,0.81,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;

    for (int ieta = 0; ieta < 3; ieta++)
    {
      temp = (TH1D*) h1d_part_pt[ieta]->Clone("temp");
      hists.push_back(temp);

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(0,25);
      temp->GetXaxis()->SetTitle("jet p_{T} [GeV]");
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < eta < %1.1f",eta_lo[ieta],eta_hi[ieta]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_pt.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%spt.csv", out_dir), hists);
  }

  // 1d particle eta histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TLegend* leg = new TLegend(0.51,0.7,0.81,0.82);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.025);
    leg->SetFillStyle(0);
    leg->SetMargin(0.1);

    TH1D* temp;

    for (int ipt = 0; ipt < 3; ipt++)
    {
      temp = (TH1D*) h1d_part_eta[ipt]->Clone("temp");
      hists.push_back(temp);

      temp->Draw("same");

      temp->GetXaxis()->SetRangeUser(-5,5);
      temp->GetXaxis()->SetTitle("#eta");
      temp->GetYaxis()->SetTitle("counts");
      temp->GetXaxis()->SetTitleOffset(1.3);
      temp->GetYaxis()->SetTitleOffset(1.5);
      temp->Draw("same hist e");
      leg->AddEntry(temp,Form("%1.1f < pt < %1.1f",pt_lo[ipt],pt_hi[ipt]));
    }
    leg->Draw("same");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_eta.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%seta.csv", out_dir), hists);
  }

  // 1d event multiplicity histogram
  mclogy(cno++);
  {
    vector<TH1*> hists;

    TH1D* temp = (TH1D*) h1d_part_mult->Clone("temp");
    hists.push_back(temp);

    temp->Draw("same");

    temp->GetXaxis()->SetRangeUser(0,200);
    temp->GetXaxis()->SetTitle("event multiplicity");
    temp->GetYaxis()->SetTitle("counts");
    temp->GetXaxis()->SetTitleOffset(1.3);
    temp->GetYaxis()->SetTitleOffset(1.5);
    temp->Draw("same hist e");

    gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_part_mult.pdf\")", cno-1, out_dir) );

    hists_to_csv(Form("%smult.csv", out_dir), hists);
  }

}


void plot_e3c_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./")
{
  // if make_ratios == 1, uses fin_name_baseline to calculate ratios as baseline, else no ratios calculated

  mcs(-1);

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_jet_eta = (TH1D*) fin->Get("h1d_jet_eta");
  h1d_jet_eta->SetName("h1d_jet_eta");

  h1d_part_mult = (TH1D*) fin->Get("h1d_part_mult");

  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    h1d_part_eta[ipt] = (TH1D*) fin->Get(Form("h1d_part_eta_%d", ipt));
    h1d_part_eta[ipt]->SetName(Form("h1d_part_eta_%d", ipt));
  }

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    h1d_jet_pt[ieta] = (TH1D*) fin->Get(Form("h1d_jet_pt_%d", ieta));
    h1d_jet_pt[ieta]->SetName(Form("h1d_jet_pt_%d", ieta));

    h1d_part_pt[ieta] = (TH1D*) fin->Get(Form("h1d_part_pt_%d", ieta));
    h1d_part_pt[ieta]->SetName(Form("h1d_part_pt_%d", ieta));

    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h3d_jet_eec_rl_xi_phi[ieta][ipt] = (TH3D*) fin->Get(Form("h3d_jet_eec_rl_xi_phi_%d_%d", ieta, ipt));
      h3d_jet_eec_rl_xi_phi[ieta][ipt]->SetName(Form("h3d_jet_eec_rl_xi_phi_%d_%d", ieta, ipt));

      h3d_jet_eec_rlsqrtpt_xi_phi[ieta][ipt] = (TH3D*) fin->Get(Form("h3d_jet_eec_rlsqrtpt_xi_phi_%d_%d", ieta, ipt));
      h3d_jet_eec_rlsqrtpt_xi_phi[ieta][ipt]->SetName(Form("h3d_jet_eec_rlsqrtpt_xi_phi_%d_%d", ieta, ipt));

      h2d_Q2_x[ieta][ipt] = (TH2D*) fin->Get(Form("h2d_Q2_x_%d_%d", ieta, ipt));
      h2d_Q2_x[ieta][ipt]->SetName(Form("h2d_Q2_x_%d_%d", ieta, ipt));
    }
  }

  particle_hists(out_dir);

  individual_hists(out_dir);

  overlay_hists(out_dir);

  raw_nice_3d_plots(out_dir);

  Q2_x_panel(out_dir);

}
