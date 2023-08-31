R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);

#include "fastjet/ClusterSequence.hh"
#include <iostream>
using namespace fastjet;
using namespace std;
const int verbosity = 0;

const Double_t Mp(0.9383); // in GeV/c^2
const Double_t Me(0.511E-3); // in GeV/c^2
const Double_t MAu(183.4343); // in GeV/c^2
const Double_t MCarbon(11.26703); // in GeV/c^2
const Double_t MCu(60.09468); // in GeV/c^2
const Double_t Md(1.87783); // in GeV/c^2
const Double_t Mu(223.49758); // in GeV/c^2
const Double_t MCa(37.55675); // in GeV/c^2
const Double_t MHe3(2.8161095); // in GeV/c^2
const Double_t MHe4(3.755675); // in GeV/c^2

const int nspecies = 9;
static double targ_A[nspecies] = {1, 197, 12, 64, 2, 40, 238, 3, 4};
static double targ_m[nspecies] = {Mp, MAu, MCarbon, MCu, Md, MCa, Mu, MHe3, MHe4};

TH1D* h1d_event_eec = NULL;
TH1D* h1d_event_eec_rlQ = NULL;
TH2D* h2d_Q2_x = NULL;

// debugging histograms
TH1D* h1d_part_pt = NULL;
TH1D* h1d_part_eta = NULL;
TH1D* h1d_part_mult = NULL;

double calculate_distance(PseudoJet p0, PseudoJet p1)
{
  float dphiabs = fabs(p0.phi() - p1.phi());

  if (dphiabs > TMath::Pi()) dphiabs = 2*TMath::Pi() - dphiabs;

  float dy = p0.eta() - p1.eta();
  return sqrt(dy*dy + dphiabs*dphiabs);
}

class Correlator_Builder
{
  private:
    vector<PseudoJet> particle_list;
    int mult;
    vector<vector<double>> pair_list;
    double Q2;
    double weight_pow;

  public:
    Correlator_Builder(vector<PseudoJet> _particle_list, float _Q2, double _weight_pow)
    {
      particle_list = _particle_list;
      mult = particle_list.size();
      Q2 = _Q2;
      weight_pow = _weight_pow; // power of the pt1 * pt2 / pt^2jet factor used in calculation
    }

    void make_pairs()
    {
      // produces mult X mult upper triangular matrix with elements being the
      // distance between partices i and j
      for (int i = 0; i < mult; i++)
      {
        vector<double> inner_list(mult, 0.0);
        for (int j = 0; j < mult; j++)
        {
          inner_list[j] = calculate_distance(particle_list[i], particle_list[j]);
        }
        pair_list.push_back(inner_list);
      }
    }

    void construct_EEC()
    {
      for (int i = 0; i < mult; i++)
      {
        for (int j = 0; j < mult; j++)
        {
          double dist12 = pair_list[i][j]; // no change in boost
          double eec_weight = pow((particle_list[i].E() * particle_list[j].E()) / Q2, weight_pow);

          // filling histograms
          h1d_event_eec->Fill(dist12, eec_weight);

          float rlQ = dist12 * sqrt(Q2);
          h1d_event_eec_rlQ->Fill(rlQ, eec_weight);

        }
      }
    }
};


void read_csv(const char* inFile = "merged.csv", int boost = 1, double proj_rest_e = 10, double targ_lab_e = 100, int targ_species = 0,
    double eec_weight_power = 1)
{
  // csv must be in the following format - eHIJING standard
  // each particle has the line
  // F << evtn << "," << p.id() << "," << p.charge() << ","
  // << p.px() << "," << p.py() << "," << p.pz() << "," << p.m() << std::endl;

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

  // boost calculation
  // calculation forces target to be 100 Gev proton, electron projectile has whatever energy neccesary to satisfy this
  TLorentzVector part;
  TLorentzVector Ei, Ef, Pf;
  Ei.SetXYZM(0, 0, -proj_rest_e, Me);
  Pf.SetXYZM(0, 0, targ_lab_e * targ_A[targ_species], targ_m[targ_species]);
  TVector3 boost_vec = Pf.BoostVector();
  Ef = Ei; Ef.Boost(boost_vec); // electron 4-vector after boost (in lab frame)
  cout<<"projectile in lab frame: (should be what you expect)"<<endl;
  Ef.Print();
  cout<<"target in lab frame: (should be what you expect)"<<endl;
  Pf.Print();

  // initialize particle level variables
  Int_t Id;
  Double_t Charge, Px, Py, Pz, Mass, Q2, xB, Pt, Eta, Mult;

  // number of lines
  int iline = 0;
  int nlines = content.size();
  int ievt;

  // loop over lines
  while (iline < nlines)
  {
    try {
      ievt = stoi(content[iline][0]); // get event number for this new event
    } catch (invalid_argument& e) {
      cout<<"@kdebug -1.5"<<endl;
      break;
    }
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<endl;

    vector<PseudoJet> charged_constituents;

    // loop over particles with this event number
    Mult = 0;
    while (iline < nlines && stoi(content[iline][0]) == ievt)
    {
      Mult++;

      // read content for this line, make type conversions
      vector<string> line;
      line = content[iline];
      Id = stoi(line[1]);
      Charge = stod(line[2]);
      Px = stod(line[3]);
      Py = stod(line[4]);
      Pz = stod(line[5]);
      Mass = stod(line[6]);
      Q2 = stod(line[7]);
      xB = stod(line[8]);

      //cout<<iline<<" "<<Id<<" "<<Charge<<" "<<Px<<" "<<Py<<" "<<Pz<<" "<<Mass<<endl;

      // apply boost to particle (boost it into lab frame)
      part.SetXYZM(Px, Py, Pz, Mass);
      if (boost == 1) part.Boost(boost_vec);

      // debug histograms
      Pt = part.Pt();
      Eta = part.Eta();
      h1d_part_pt->Fill(Pt);
      h1d_part_eta->Fill(Eta);

      // use all charged fsp particles w/ in theta [10,70] region with momentum >= 0.2 GeV/c
      if (part.Theta() > 10 && part.Theta() < 70 && Id!=11 && part.E() >= 0.2 && Charge != 0)
      {
        PseudoJet constit = PseudoJet(part.Px(),part.Py(),part.Pz(),part.E());
        charged_constituents.push_back(constit);
      }

      iline++;
    }
    h1d_part_mult->Fill(Mult);
    h2d_Q2_x->Fill(Q2, xB);

    if (charged_constituents.size() < 1) continue;

    // eec calculation
    Correlator_Builder cb(charged_constituents, Q2, eec_weight_power);
    cb.make_pairs();
    cb.construct_EEC();

  }

  fin.close();

}


void eec_hists(const char* inFile = "merged.root", const char* outFile = "hists_eec.root", const int gen_type = 1,
    double proj_rest_e = 2131.56, double targ_lab_e = 100, int targ_species = 0, double eec_weight_power = 1,
    int boost = 0)
{
  // proj_rest_e = energy of projectile beam in target rest frame, leave blank if pythia
  // targ_lab_e = energy of target beam in lab frame, leave blank if pythia
  // all energies positive and in GeV units
  // only for e+A collsions, specify A with targ_species, =0 for p, =1 for Au
  // gen_type = 0 for pythia6, =-1 for pyhtia8, =1 or anything for eHIJING (DIFFERENT FROM Q2_x.C settings)
  // boost = 0, do not apply boost =1 apply boost according to targ_lab_e and targ_species

  cout << "Generator Type: ";
  if (gen_type==0) cout << "Pythia6" << endl;
  else if (gen_type==-1) cout << "Pythia8" << endl;
  else cout << "eHIJING" << endl;

  // compute log bins for eec histogram
  double xmin = 1E-4;
  double xmax = 1;
  int nbins = 50;
  Double_t lbins[nbins+1];
  double binwidth = (log10(xmax) - log10(xmin)) / nbins;
  for (int i = 0; i < nbins+1; i++)
  {
    lbins[i] = TMath::Power(10, log10(xmin) + binwidth * i);
  }

  xmin = 5E-3;
  xmax = 50;
  nbins = 50;
  Double_t lbins_rlQ[nbins+1];
  binwidth = (log10(xmax) - log10(xmin)) / nbins;
  for (int i = 0; i < nbins+1; i++)
  {
    lbins_rlQ[i] = TMath::Power(10, log10(xmin) + binwidth * i);
  }

  // compute log bins for Q2-x
  // xbins correspond to xB
  xmin = 1E-8;
  xmax = 1;
  int xnbins = 50;
  Double_t xlbins[xnbins+1];
  double xbinwidth = (log10(xmax) - log10(xmin)) / xnbins;
  for (int i = 0; i < xnbins+1; i++)
  {
    xlbins[i] = TMath::Power(10, log10(xmin) + xbinwidth * i);
  }

  // ybins correspond to Q2
  double ymin = 1;
  double ymax = 1E2;
  int ynbins = 50;
  Double_t ylbins[ynbins+1];
  double ybinwidth = (log10(ymax) - log10(ymin)) / ynbins;
  for (int i = 0; i < ynbins+1; i++)
  {
    ylbins[i] = TMath::Power(10, log10(ymin) + ybinwidth * i);
  }

  h1d_part_pt = new TH1D("h1d_part_pt","particle pt",400,0,100);
  h1d_part_pt->Sumw2();

  h1d_part_eta = new TH1D("h1d_part_eta","particle eta",1600,-15,5);
  h1d_part_eta->Sumw2();
  
  h1d_part_mult = new TH1D("h1d_part_mult", "event particle multiplicity",200,0,200);
  h1d_part_mult->Sumw2();

  h2d_Q2_x = new TH2D("h2d_Q2_x","Q2_x",50,xlbins,50,ylbins);
  h2d_Q2_x->Sumw2();

  // reads file and fills in jet_constits
  read_csv(inFile, boost, proj_rest_e, targ_lab_e, targ_species, eec_weight_power); // eHIJING, assumes target frame and boosts to EIC
  cout<<"@kdebug last"<<endl;

  // write out histograms
  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();

  h1d_part_pt->Write();
  h1d_part_eta->Write();
  h1d_part_mult->Write();
  h1d_Q2_x->Write();

  h1d_event_eec->Write();
  h1d_event_eec_rlQ->Write();

}
