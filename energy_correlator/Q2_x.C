R__LOAD_LIBRARY(libeicsmear);

#include <iostream>
using namespace std;
const int verbosity = 0;

const Double_t Mp(0.9383); // in GeV/c^2

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 40, 5};
static double pt_hi[ptbin] = {10, 20, 40, 60, 60};

const int etabin = 4; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5};

TH2D* h2d_Q2_x[etabin][ptbin] = {};

void read_csv(const char* inFile = "merged.csv")
{
  // csv must be in the following format - eHIJING standard
  // each particle has the line
  // F << evtn << "," << Q2 << "," << xB << std::endl;

  // set up file and content vector
  fstream fin;
  fin.open(inFile, ios::in);
  vector<vector<string>> content;
  string line_str, element_str;

  // read in file, write to content 2D vector
  // note all elements of 2D vector are strings
  cout<<"reading in file..."<<endl;
  while (getline(fin, line_str)) // iterate over lines
  {
    stringstream str(line_str);

    vector<string> line;
    while (getline(str, element_str, ',')) // iterate over elements in line
    {
      line.push_back(element_str);
    }
    content.push_back(line);
  }
  cout<<"file read"<<endl;

  // initialize particle level variables
  Double_t Q2, xB;

  // number of lines
  int iline = 0;
  int nlines = content.size();

  cout << "nlines: " << nlines<<endl;

  // loop over lines
  while (iline < nlines)
  {

    if (iline%10000==0) cout<<"Processing event = "<<iline<<endl;

    // read content for this line, make type conversions
    vector<string> line;
    line = content[iline];
    Q2 = stod(line[1]);
    xB = stod(line[2]);

    cout<<iline<<" "<<Q2<<" "<<xB<<endl;

    iline++;

  }

  fin.close();

}


void Q2_x(const char* inFile = "merged.root", const char* outFile = "hists_Q2_x.root")
{

  // compute log bins for Q2-x
  double xmin = 1E-6;
  double xmax = 1;
  int xnbins = 50;
  Double_t xlbins[xnbins+1];
  double xbinwidth = (log10(xmax) - log10(xmin)) / xnbins;
  for (int i = 0; i < xnbins+1; i++)
  {
    xlbins[i] = TMath::Power(10, log10(xmin) + xbinwidth * i);
  }

  double ymin = 1;
  double ymax = 1E4;
  int ynbins = 50;
  Double_t ylbins[ynbins+1];
  double ybinwidth = (log10(ymax) - log10(ymin)) / ynbins;
  for (int i = 0; i < ynbins+1; i++)
  {
    ylbins[i] = TMath::Power(10, log10(ymin) + ybinwidth * i);
  }

  // histogram definition
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h2d_Q2_x[ieta][ipt] = new TH2D(Form("h2d_Q2_x_%d_%d", ieta, ipt),"Q2_x",50,xlbins,50,ylbins);
      h2d_Q2_x[ieta][ipt]->Sumw2();
    }
  }


  // reads file
  read_csv(inFile);
  cout<<"@kdebug last"<<endl;

  // write out histograms
  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h2d_Q2_x[ieta][ipt]->Write();
      cout<<"h1d_Q2_x_"<<ieta<<"_"<<ipt<<" entries:"<<h2d_Q2_x[ieta][ipt]->GetEntries()<<endl;
    }
  }

}
