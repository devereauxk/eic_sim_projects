R__LOAD_LIBRARY(libeicsmear);
using namespace std;

const double ETA_LO = -3.5;
const double ETA_HI = 3.5;
const int N_BINS = 141;

// TH1D is inclusive on lower bin edge and exclusive on upper

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
  tree->ReadFile(inFile, "Momentum/D:Theta:Eta:Deltap-p:DCA-rPhi:DCA-z");

  int n;
  double bin_lo;
  double bin_hi;
  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    bin_lo = Res_Handler->GetBinLowEdge(ibin);
    bin_hi = bin_lo + Res_Handler->GetBinWidth(ibin);

    // gmom_res
    n = tree->Draw(Form("Momentum:Deltap-p", "Eta >= %d && Eta < %d ", "goff"), bin_lo, bin_hi);
    gmom_res[ibin] = new TGraph(n, tree->GetV1(), tree->GetV2());
    gmom_res[ibin]->SetName(Form("gmom_res_%i", ibin));

    /*
    gdca_rphi_res[ibin] = TGraph();
    gdca_rphi_res[ibin]->SetName(Form("gdca_rphi_res_%i", ibin));
    gdca_z_res[ibin] = TGraph();
    gdca_z_res[ibin]->SetName(Form("gdca_z_res_%i", ibin));
    */
  }


  // add objects to outfile
  TFile* fout = new TFile(outFile,"recreate");
  Res_Handler->Write();

  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    gmom_res[ibin]->Write();
  }

  fout->Write();

}
