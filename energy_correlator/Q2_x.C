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

const int ptbin = 5; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double pt_lo[ptbin] = {5, 10, 20, 40, 5};
static double pt_hi[ptbin] = {10, 20, 40, 60, 60};

const int etabin = 4; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5};

TH2D* h2d_Q2_x[etabin][ptbin] = {};

void read_root(const char* inFile = "merged.root", int gen_type = 0)
{
  //Event Class
  erhic::EventHepMC *event(NULL);
  if (gen_type==0) erhic::EventPythia *event(NULL); // pythia6

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nevt = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nevt<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  int total_jets = 0;

  //Loop Over Events
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    tree->GetEntry(ievt);

    if (ievt%1000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    double xB = event->GetX();
    double Q2 = event->GetQ2();

    /*
    //Write Out Q2
    double Q2 = event->GetQ2(); //Can also do event->QSquared
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);

    // process cuts
    bool flag_direct = false;
    if (event->GetProcess()==99) flag_direct = true;
    if (event->GetProcess()==131 || event->GetProcess()==132) flag_direct = true;
    if (event->GetProcess()==135 || event->GetProcess()==136) flag_direct = true;
    if (!flag_direct) continue; // only process direct processes
    */

    // particle enumeration, addition to jet reco setup, and total pt calculation
    erhic::ParticleMC* particle;
    vector<PseudoJet> jet_constits;
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      particle = event->GetTrack(ipart);

      // use all fsp particles w/ < 3.5 eta, not including scattered electron, for jet reconstruction
      if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11)
      {
        PseudoJet constit = PseudoJet(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
        constit.set_user_index(ipart);
        jet_constits.push_back(constit);
      }

    }

    // jet reconstruction
    JetDefinition R1jetdef(antikt_algorithm, 1.0);
    ClusterSequence cs(jet_constits, R1jetdef);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    // jet processing
    for (unsigned ijet = 0; ijet < jets.size(); ijet++)
    {
      // cuts on jet kinematics, require jet_pt >= 5GeV, |jet_eta| <= 2.5
      if (jets[ijet].pt() < 5 || fabs(jets[ijet].eta()) > 2.5) continue;

      total_jets++;

      // fill Q2-x histograms
      for (int ieta = 0; ieta < etabin; ieta++)
      {
        if (jets[ijet].eta() >= eta_lo[ieta] && jets[ijet].eta() < eta_hi[ieta])
        {
          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            if (jets[ijet].pt() >= pt_lo[ipt] && jets[ijet].pt() < pt_hi[ipt])
            {
              h2d_Q2_x[ieta][ipt]->Fill(xB, Q2);
            }
          }
        }
      }

    }

  }

  cout<<"total num jets = "<<total_jets<<endl;

}

void read_csv(const char* inFile = "merged.csv", double proj_rest_e = 10, double targ_lab_e = 100, int targ_species = 0)
{
  // csv must be in the following format - eHIJING standard
  // each particle has the line
  // F << evtn << "," << p.id() << "," << p.charge() << ","
  // << p.px() << "," << p.py() << "," << p.pz() << "," << p.m() << std::endl;
  // << Q2 << xB

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
  TLorentzVector part_rest, part_lab;
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
  Double_t Charge, Px, Py, Pz, Mass, Q2, xB;

  // number of lines
  int iline = 0;
  int nlines = content.size();
  int ievt;

  int total_jets = 0;

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

    vector<PseudoJet> jet_constits;

    // loop over particles with this event number
    while (iline < nlines && stoi(content[iline][0]) == ievt)
    {
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
      part_rest.SetXYZM(Px, Py, Pz, Mass);
      //part_rest.Print();
      part_lab = part_rest; part_lab.Boost(boost_vec);
      //part_lab.Print();

      // use all fsp particles w/ < 3.5 eta, not including scattered electron, for jet reconstruction
      //cout<<"part_lab eta:"<<part_lab.Eta()<<endl;
      if (fabs(part_lab.Eta())<3.5 && Id!=11)
      {
        PseudoJet constit = PseudoJet(part_lab.Px(),part_lab.Py(),part_lab.Pz(),part_lab.E());
        constit.set_user_index(iline);
        jet_constits.push_back(constit);
      }

      iline++;
    }

    // jet reconstruction
    JetDefinition R1jetdef(antikt_algorithm, 1.0);
    ClusterSequence cs(jet_constits, R1jetdef);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    //cout<<"n jets:"<<jets.size()<<endl;

    // jet processing
    for (unsigned ijet = 0; ijet < jets.size(); ijet++)
    {
      // cuts on jet kinematics, require jet_pt >= 5GeV, |jet_eta| <= 2.5
      if (jets[ijet].pt() < 5 || fabs(jets[ijet].eta()) > 2.5) continue;

      total_jets++;

      // fill Q2-x histograms
      for (int ieta = 0; ieta < etabin; ieta++)
      {
        if (jets[ijet].eta() >= eta_lo[ieta] && jets[ijet].eta() < eta_hi[ieta])
        {
          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            if (jets[ijet].pt() >= pt_lo[ipt] && jets[ijet].pt() < pt_hi[ipt])
            {
              h2d_Q2_x[ieta][ipt]->Fill(xB, Q2);
            }
          }
        }
      }

    }

  }

  cout<<"total num jets = "<<total_jets<<endl;

  fin.close();

}


void Q2_x(const char* inFile = "merged.root", const char* outFile = "hists_eec.root", int gen_type = 2,
    double proj_rest_e = 2131.56, double targ_lab_e = 100, int targ_species = 0)
{

  cout << "Generator Type: ";
  if (gen_type==0) cout << "Pythia6" << endl;
  else if (gen_type==1) cout << "Pythia8" << endl;
  else cout << "eHIJING" << endl;

  // compute log bins for Q2-x
  // xbins correspond to xB
  double xmin = 1E-6;
  double xmax = 1;
  int xnbins = 50;
  Double_t xlbins[xnbins+1];
  double xbinwidth = (log10(xmax) - log10(xmin)) / xnbins;
  for (int i = 0; i < xnbins+1; i++)
  {
    xlbins[i] = TMath::Power(10, log10(xmin) + xbinwidth * i);
  }

  // ybins correspond to Q2
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
  if (gen_type <= 1) read_root(inFile, gen_type); // assumes lab frame, pythia6 (EventPythia) or pythia8 (EventHepMC)
  else read_csv(inFile, proj_rest_e, targ_lab_e, targ_species); // assumes target frame, eHIJING
  cout<<"@kdebug last"<<endl;


  // write out histograms
  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h2d_Q2_x[ieta][ipt]->Write();
      cout<<"h2d_Q2_x_"<<ieta<<"_"<<ipt<<" entries:"<<h2d_Q2_x[ieta][ipt]->GetEntries()<<endl;
    }
  }

}
