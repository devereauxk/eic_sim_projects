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

const int speciesnum = 9;
static char* species[speciesnum] = {(char*)"e+p", (char*)"e+D", (char*)"e+He-3", (char*)"e+He-4", (char*)"e+C", (char*)"e+Ca", (char*)"e+Cu", (char*)"e+Au", (char*)"e+U"};
static int species_A[speciesnum] = {1, 2, 3, 4, 12, 40, 64, 197, 238};

static char* fname_eA[speciesnum] = {(char*)"./analysis/ep_10_100_K0_density_pow05/merged.root", (char*)"", (char*)"", (char*)"", (char*)"./analysis/eC_10_100_K4_density_pow05/merged.root", (char*)"", (char*)"", \
  (char*)"./analysis/eAu_10_100_K4_density_pow05/merged.root", (char*)""};

const char* out_dir = "./paperplots/";

// jet level histograms
TH3D* h3d_jet_eec_rl_xi_phi[speciesnum][etabin][ptbin] = {};
/*
TH3D* h3d_jet_eec_rlsqrtpt_xi_phi[etabin][ptbin] = {};
TH1D* h1d_jet_pt[etabin] = {};
TH1D* h1d_jet_eta = NULL;
TH1D* h1d_jet_multiplicity[etabin][ptbin] = {};
TH1D* h1d_jet_multiplicity_charged[etabin][ptbin] = {};
TH2D* h2d_Q2_x[etabin][ptbin] = {};
*/

// particle level histograms
/*
TH1D* h1d_part_pt[etabin] = {};
TH1D* h1d_part_eta[ptbin] = {};
TH1D* h1d_part_mult = NULL;
*/

const float rl_norm_hi = 0.08; // 0.3 used for 3by3 top left and second col //0.08;
const float rl_norm_lo = 1E-3;

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

void e3c_projected_hists()
{
  const int nspecies_picks = 3;
  static int species_picks[nspecies_picks] = {0, 4, 7}; // {0, 7, 9};

  // 3x3 panel
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < 3; ipt++)
    {
      mclogxy(cno++);
      {
        float plot_xrange_lo = 1E-1;
        float plot_xrange_hi = 1;
        float plot_yrange_lo = 1E-2;
        float plot_yrange_hi = 3E-1;
        float legend_x = 0.7;
        float legend_y = 0.2;

        TLegend* leg = new TLegend(legend_x,legend_y,legend_x+0.3,legend_y+0.15);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.028);
        leg->SetFillStyle(0);
        leg->SetMargin(0.1);

        int species_pick;
        TH1D* temp;
        TH1D* temp_baseline = (TH1D*) h3d_jet_eec_rl_xi_phi[0][ieta][ipt]->ProjectionX();
        for (int ispecies = 0; ispecies < nspecies_picks; ispecies++)
        {
          species_pick = species_picks[ispecies];
          temp = (TH1D*) h3d_jet_eec_rl_xi_phi[species_pick][ieta][ipt]->ProjectionX();

          int norm_binrange_lo = temp->FindBin(rl_norm_lo);
          int norm_binrange_hi = temp->FindBin(rl_norm_hi);
          if (norm_binrange_lo == 0)
          {
            norm_binrange_lo = 1;
            cout<<"bin range lo too low; set to 1"<<endl;
          }
          if (norm_binrange_hi > temp->GetNbinsX())
          {
            norm_binrange_hi = temp->GetNbinsX();
            cout<<"bin range hi too high; set to "<<temp->GetNbinsX()<<endl;
          }
          double relative_normalization =  temp_baseline->Integral(norm_binrange_lo,norm_binrange_hi) / temp->Integral(norm_binrange_lo,norm_binrange_hi);
          temp->Scale(relative_normalization);
          temp->Scale(1/temp_baseline->Integral());

          // plot histogram
          temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
          temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
          temp->GetXaxis()->SetTitle("R_{L}");
          temp->GetYaxis()->SetTitle("Projected E3C (normalized at low R_{L})");
          temp->SetMarkerColor(pt_color[ispecies]);
          temp->SetLineColor(pt_color[ispecies]);
          temp->SetLineWidth(2);
          temp->SetMarkerSize(1);
          temp->SetMarkerStyle(21);
          temp->Draw("same hist e");
          leg->AddEntry(temp,Form("%s", species[species_pick]));
        }
        leg->Draw("same");

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.028);
        tl->SetTextColor(kBlack);
        tl->DrawLatexNDC(0.22,0.84,"10+100 GeV, 4*10^{8} events");
        tl->DrawLatexNDC(0.22,0.81,Form("#eta #in [%.1f, %0.1f)",eta_lo[ieta],eta_hi[ieta]));
        tl->DrawLatexNDC(0.22,0.78,Form("p_{T,jet} #in [%.1f, %0.1f)",pt_lo[ipt],pt_hi[ipt]));

        gROOT->ProcessLine( Form("cc%d->Print(\"%sh1d_jet_e3c_overlay_%i_%i.pdf\")", cno-1, out_dir, ieta, ipt) );

      }
    }
  }

}

void xi_phi_difference_hists(int num_species_pick, int denom_species_pick)
{
  // calculates (alpha*num_species_pick - denom_species_pick) * bin-by-bin normalization
  int etabin_pick = 2;
  int ptbin_pick = 1;
  float rl_range_lo = 1E-2;
  float rl_range_hi = 1;

  TH3D* picked = (TH3D*) h3d_jet_eec_rl_xi_phi[num_species_pick][etabin_pick][ptbin_pick]->Clone();
  TH3D* baseline = (TH3D*) h3d_jet_eec_rl_xi_phi[denom_species_pick][etabin_pick][ptbin_pick]->Clone();

  TH1D* picked_x = h3d_jet_eec_rl_xi_phi[num_species_pick][etabin_pick][ptbin_pick]->ProjectionX();
  TH1D* baseline_x = h3d_jet_eec_rl_xi_phi[denom_species_pick][etabin_pick][ptbin_pick]->ProjectionX();

  int norm_binrange_lo = picked_x->FindBin(rl_range_lo);
  cout<<"low bin: "<<norm_binrange_lo<<endl;
  int norm_binrange_hi = picked_x->FindBin(rl_range_hi);
  cout<<"high bin: "<<norm_binrange_hi<<endl;
  if (norm_binrange_lo == 0)
  {
    norm_binrange_lo = 1;
    cout<<"bin range lo too low; set to 1"<<endl;
  }
  if (norm_binrange_hi > picked_x->GetNbinsX()+1)
  {
    norm_binrange_hi = picked_x->GetNbinsX()+1;
    cout<<"bin range hi too high; set to "<<picked_x->GetNbinsX()<<endl;
  }
  double relative_normalization =  baseline_x->Integral(norm_binrange_lo,norm_binrange_hi) / picked_x->Integral(norm_binrange_lo,norm_binrange_hi);
  picked->Scale(relative_normalization);
  picked->Add(baseline, -1);
  picked->Scale(1/baseline_x->Integral());

  TH3D* sliced;
  TH2D* temp;
  for (int ibin = norm_binrange_lo; ibin <= norm_binrange_hi; ibin++)
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.5; // TMath::Pi() / 2.0;
    float plot_zrange_lo = -0.1E-3;
    float plot_zrange_hi = 0.3E-3;

    float bin_center = picked->GetXaxis()->GetBinCenter(ibin);
    sliced = (TH3D*) picked->Clone("temp3d");
    sliced->GetXaxis()->SetRange(ibin,ibin);
    temp = (TH2D*) sliced->Project3D("zy");

    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
    {
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);
      temp->GetXaxis()->SetTitle("#xi");
      temp->GetYaxis()->SetTitle("#phi");
      temp->Draw("colz");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.03);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %.1f), p_{T} #in [%.1f, %.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick],pt_lo[ptbin_pick],pt_hi[ptbin_pick]));
      tl->DrawLatexNDC(0.22,0.81,Form("R_{L} ~ %.3f", bin_center));
      tl->DrawLatexNDC(0.22,0.78,Form("%s / %s",species[num_species_pick], species[denom_species_pick]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d-diff-%d_%d.pdf\")", cno-1, out_dir, num_species_pick, denom_species_pick, ibin) );
    }

    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
    {
      temp->Draw("SURF2Z");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d-diff-%d_%d_surface.pdf\")", cno-1, out_dir, num_species_pick, denom_species_pick, ibin) );
    }

  }

}

void xi_phi_ratio_hists(int num_species_pick, int denom_species_pick)
{
  // calculates (alpha*num_species_pick - denom_species_pick) / denom_species_pick
  int etabin_pick = 2;
  int ptbin_pick = 1;
  float rl_range_lo = 1E-2;
  float rl_range_hi = 1;

  TH3D* picked = (TH3D*) h3d_jet_eec_rl_xi_phi[num_species_pick][etabin_pick][ptbin_pick]->Clone();
  TH3D* baseline = (TH3D*) h3d_jet_eec_rl_xi_phi[denom_species_pick][etabin_pick][ptbin_pick]->Clone();

  TH1D* picked_x = h3d_jet_eec_rl_xi_phi[num_species_pick][etabin_pick][ptbin_pick]->ProjectionX();
  TH1D* baseline_x = h3d_jet_eec_rl_xi_phi[denom_species_pick][etabin_pick][ptbin_pick]->ProjectionX();

  int norm_binrange_lo = picked_x->FindBin(rl_range_lo);
  cout<<"low bin: "<<norm_binrange_lo<<endl;
  int norm_binrange_hi = picked_x->FindBin(rl_range_hi);
  cout<<"high bin: "<<norm_binrange_hi<<endl;
  if (norm_binrange_lo == 0)
  {
    norm_binrange_lo = 1;
    cout<<"bin range lo too low; set to 1"<<endl;
  }
  if (norm_binrange_hi > picked_x->GetNbinsX()+1)
  {
    norm_binrange_hi = picked_x->GetNbinsX()+1;
    cout<<"bin range hi too high; set to "<<picked_x->GetNbinsX()<<endl;
  }
  double relative_normalization =  baseline_x->Integral(norm_binrange_lo,norm_binrange_hi) / picked_x->Integral(norm_binrange_lo,norm_binrange_hi);
  picked->Scale(relative_normalization);
  picked->Add(baseline, -1);
  picked->Divide(baseline);

  TH3D* sliced;
  TH2D* temp;
  for (int ibin = norm_binrange_lo; ibin <= norm_binrange_hi; ibin++)
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 1;
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.5; // TMath::Pi() / 2.0;
    float plot_zrange_lo = -0.5;
    float plot_zrange_hi = 1.5;

    float bin_center = picked->GetXaxis()->GetBinCenter(ibin);
    sliced = (TH3D*) picked->Clone("temp3d");
    sliced->GetXaxis()->SetRange(ibin,ibin);
    temp = (TH2D*) sliced->Project3D("zy");

    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
    {
      temp->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
      temp->GetYaxis()->SetRangeUser(plot_yrange_lo,plot_yrange_hi);
      temp->GetZaxis()->SetRangeUser(plot_zrange_lo, plot_zrange_hi);
      temp->GetXaxis()->SetTitle("#xi");
      temp->GetYaxis()->SetTitle("#phi");
      temp->Draw("colz");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.03);
      tl->SetTextColor(kBlack);
      tl->DrawLatexNDC(0.22,0.84,Form("#eta #in [%.1f, %.1f), p_{T} #in [%.1f, %.1f)",eta_lo[etabin_pick],eta_hi[etabin_pick],pt_lo[ptbin_pick],pt_hi[ptbin_pick]));
      tl->DrawLatexNDC(0.22,0.81,Form("R_{L} ~ %.3f", bin_center));
      tl->DrawLatexNDC(0.22,0.78,Form("%s / %s",species[num_species_pick], species[denom_species_pick]));

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d-ratio-%d_%d.pdf\")", cno-1, out_dir, num_species_pick, denom_species_pick, ibin) );
    }

    mcs(cno++, 0, 0, 400, 400, 0.12, 0.15, 0.1, 0.13);
    {
      temp->Draw("SURF2Z");

      gROOT->ProcessLine( Form("cc%d->Print(\"%sh2d_jet_e3c_xi_phi_%d-ratio-%d_%d_surface.pdf\")", cno-1, out_dir, num_species_pick, denom_species_pick, ibin) );
    }

  }

}

void plot_e3c_ratios()
{
  // if make_ratios == 1, uses fin_name_baseline to calculate ratios as baseline, else no ratios calculated

  mcs(-1);

  // load histograms
  TFile* fin;
  char* fin_name;
  for (int ispecies = 0; ispecies < speciesnum; ispecies++)
  {
    fin_name = fname_eA[ispecies];

    if (strcmp(fin_name,"") == 0) continue;

    fin = new TFile(fin_name, "READ");

    for (int ieta = 0; ieta < etabin; ieta++)
    {
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        h3d_jet_eec_rl_xi_phi[ispecies][ieta][ipt] = (TH3D*) fin->Get(Form("h3d_jet_eec_rl_xi_phi_%d_%d", ieta, ipt));
        h3d_jet_eec_rl_xi_phi[ispecies][ieta][ipt]->SetName(Form("h3d_jet_eec_rl_xi_phi_%d_%d_%d", ispecies, ieta, ipt));
      }
    }
    cout<<fin_name<<" loaded!"<<endl;

  }

  // make the plots
  e3c_projected_hists();

  xi_phi_difference_hists(7, 0); // eAu - ep, where eAu is without density implementation
  xi_phi_difference_hists(9, 0); // eAu - ep, where eAu is with density implementation

  xi_phi_ratio_hists(7, 0); // (eAu - ep) / ep, where eAu is without density implementation
  xi_phi_ratio_hists(9, 0); // (eAu - ep) / ep, where eAu is with density implementation

  xi_phi_ratio_hists(9, 7); // (eAu[with density] - eAu[without denisty]) / eAu[without density]

}
