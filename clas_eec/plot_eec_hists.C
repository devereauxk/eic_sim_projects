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

TH1D* h1d_event_eec = NULL;
TH1D* h1d_event_eec_rlQ = NULL;
TH2D* h2d_Q2_x = NULL;

// debugging histograms
TH1D* h1d_part_pt = NULL;
TH1D* h1d_part_eta = NULL;
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

void push_to_csvs(const char* out_dir)
{

  // eec binned in rl dist
  vector<TH1*> hists;
  TH1D* temp = (TH1D*) h1d_event_eec->Clone("temp");
  hists.push_back(temp);
  hists_to_csv(Form("%seec.csv", out_dir), hists);
  hists.clear()

  // eec binned in rlQ dist
  TH1D* temp = (TH1D*) h1d_event_eec_rlQ->Clone("temp");
  hists.push_back(temp);
  hists_to_csv(Form("%seec_rlQ.csv", out_dir), hists);
  hists.clear()

  // 1d event particle pt histogram
  TH1D* temp = (TH1D*) h1d_part_pt->Clone("temp");
  hists.push_back(temp);
  hists_to_csv(Form("%spt.csv", out_dir), hists);
  hists.clear()

  // 1d event particle eta histogram
  TH1D* temp = (TH1D*) h1d_part_eta->Clone("temp");
  hists.push_back(temp);
  hists_to_csv(Form("%seta.csv", out_dir), hists);
  hists.clear()

  // 1d event multiplicity histogram
  TH1D* temp = (TH1D*) h1d_part_mult->Clone("temp");
  hists.push_back(temp);
  hists_to_csv(Form("%smult.csv", out_dir), hists);
  hists.clear()

}


void plot_eec_hists(const char* fin_name = "hists_eec.root", const char* out_dir = "./")
{

  // load histograms
  TFile* fin = new TFile(fin_name, "READ");

  h1d_event_eec = (TH1D*) fin->Get("h1d_event_eec");
  h1d_event_eec->SetName("h1d_event_eec");

  h1d_event_eec_rlQ = (TH1D*) fin->Get("h1d_event_eec_rlQ");
  h1d_event_eec_rlQ->SetName("h1d_event_eec_rlQ");

  h1d_part_pt = (TH1D*) fin->Get("h1d_part_pt");
  h1d_part_pt->SetName("h1d_part_pt");

  h1d_part_eta = (TH1D*) fin->Get("h1d_part_eta");
  h1d_part_eta->SetName("h1d_part_eta");

  h1d_part_mult = (TH1D*) fin->Get("h1d_part_mult");
  h1d_part_mult->SetName("h1d_part_mult");

  push_to_csvs(outdir);


}
