R__LOAD_LIBRARY(libeicsmear);
using namespace std;

const double ETA_LO = -3.525;
const double ETA_HI = 3.525;
const int N_BINS = 141;

// using exact same setup for Res_Handler as ATHENA_resolutions_r.root. However,
// However, for_Wenqing.csv has data only for -3 < eta < 3 instead of
// -3.5 < eta < 3.5 as is the case with ATHEA_resolutions_r.root
// Note: TH1D is inclusive on lower bin edge and exclusive on upper

TH1F* Res_Handler = NULL;
TGraph *gmom_res[N_BINS];
TGraph *gdca_rphi_res[N_BINS];
TGraph *gdca_z_res[N_BINS];

void convert_res_csv_to_root(const char* inFile = "for_Wenqing.csv", const char* outFile = "resolutions.root")
{

  // setup Res_Handler
  Res_Handler = new TH1F("Res_Handler", "Res_Handler", N_BINS, ETA_LO, ETA_HI);

  //setup and fill TGraphs
  TTree* tree = new TTree("tree from csv", "tree from csv");
  tree->ReadFile(inFile, "Momentum/D:Theta:Eta:DeltaP:DCArPhi:DCAz");

  int n;
  double bin_lo;
  double bin_hi;
  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    bin_lo = Res_Handler->GetBinLowEdge(ibin);
    bin_hi = bin_lo + Res_Handler->GetBinWidth(ibin);

    // gmom_res
    n = tree->Draw("Momentum:DeltaP", Form("Eta >= %f && Eta < %f", bin_lo, bin_hi), "goff");
    gmom_res[ibin] = new TGraph(n, tree->GetV1(), tree->GetV2());
    gmom_res[ibin]->SetName(Form("gmom_res_%i", ibin));

    // gdca_rphi_res
    n = tree->Draw("Momentum:DCArPhi", Form("Eta >= %f && Eta < %f", bin_lo, bin_hi), "goff");
    gdca_rphi_res[ibin] = new TGraph(n, tree->GetV1(), tree->GetV2());
    gdca_rphi_res[ibin]->SetName(Form("gdca_rphi_res_%i", ibin));

    //gdca_z_res
    n = tree->Draw("Momentum:DCAz", Form("Eta >= %f && Eta < %f", bin_lo, bin_hi), "goff");
    gdca_z_res[ibin] = new TGraph(n, tree->GetV1(), tree->GetV2());
    gdca_z_res[ibin]->SetName(Form("gdca_z_res_%i", ibin));
  }


  // add objects to outfile
  TFile* fout = new TFile(outFile,"recreate");
  Res_Handler->Write();

  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    gmom_res[ibin]->Write();
    gdca_rphi_res[ibin]->Write();
    gdca_z_res[ibin]->Write();
  }

  fout->Write();

}
