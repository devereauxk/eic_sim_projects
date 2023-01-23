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

TH1D* h1d_jet_eec[etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt[etabin][ptbin] = {};
TH1D* h1d_jet_pt = NULL;
TH1D* h1d_jet_eta = NULL;

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
    double jet_pt;
    double jet_eta;
    double weight_pow;

  public:
    Correlator_Builder(vector<PseudoJet> _particle_list, float _jet_pt, float _jet_eta, double _weight_pow)
    {
      particle_list = _particle_list;
      mult = particle_list.size();
      jet_pt = _jet_pt;
      jet_eta = _jet_eta;
      weight_pow = _weight_pow; // power of the pt1 * pt2 / pt^2jet factor used in calculation
    }

    void make_pairs()
    {
      // produces mult X mult upper triangular matrix with elements being the
      // distance between partices i and j
      for (int i = 0; i < mult; i++)
      {
        vector<double> inner_list(mult, 0.0);
        for (int j = i; j < mult; j++)
        {
          inner_list[j] = calculate_distance(particle_list[i], particle_list[j]);
        }
        pair_list.push_back(inner_list);
      }
    }

    void construct_EEC()
    {
      int overlap = 0;
      for (int i = 0; i < mult; i++)
      {
        for (int j = i; j < mult; j++)
        {
          double dist12 = pair_list[i][j]; // only thing to change in boost
          double eec_weight = pow((particle_list[i].pt() * particle_list[j].pt()) / pow(jet_pt, 2), weight_pow);

          if (i == j) overlap++;
          if (overlap == 0) eec_weight = eec_weight*2;
          if (overlap > 0) eec_weight = eec_weight*1;

          // filling h1d_jet_eec[etabin][ptbin] histograms
          // filing h1d_jet_eec_rlsqrtpt[etabin][ptbin] histograms
          for (int ieta = 0; ieta < etabin; ieta++)
          {
            if (jet_eta >= eta_lo[ieta] && jet_eta < eta_hi[ieta])
            {
              for (int ipt = 0; ipt < ptbin; ipt++)
              {
                if (jet_pt >= pt_lo[ipt] && jet_pt < pt_hi[ipt])
                {
                  h1d_jet_eec[ieta][ipt]->Fill(dist12, eec_weight);

                  float rlsqrtpt = dist12 * sqrt(jet_pt);
                  h1d_jet_eec_rlsqrtpt[ieta][ipt]->Fill(rlsqrtpt, eec_weight);
                }
              }
            }
          }

        }
      }
    }
};

void read_root(const char* inFile = "merged.root", double eec_weight_power = 1)
{
  //Event Class
  erhic::EventPythia *event(NULL);

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nevt = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nevt<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  //Loop Over Events
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    tree->GetEntry(ievt);

    if (ievt%1000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

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
      // jet histograms filled on inclusive jet information
      h1d_jet_pt->Fill(jets[ijet].pt());
      h1d_jet_eta->Fill(jets[ijet].eta());

      // cuts on jet kinematics, require jet_pt >= 5GeV, |jet_eta| <= 2.5
      if (jets[ijet].pt() < 5 || fabs(jets[ijet].eta()) > 2.5) continue;

      // cuts on jet constituent kinematics, require consitituents_pt >= 0.5GeV, |consitituents_eta| <= 3.5
      // take only charged constituents for eec calculation
      vector<PseudoJet> constituents = jets[ijet].constituents();
      vector<PseudoJet> charged_constituents;
      for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
      {
        if (constituents[iconstit].pt() < 0.5 || fabs(constituents[iconstit].eta()) > 3.5) continue;

        int ip = constituents[iconstit].user_index();
        float charge = event->GetTrack(ip)->Id().Info()->Charge();
        //cout<<"constituent pt:"<<constituents[iconstit].pt()<<" track pt:"<<event->GetTrack(ip)->GetPt()<<" charge:"<<charge<<endl;
        if (charge != 0) charged_constituents.push_back(constituents[iconstit]);
      }

      if (charged_constituents.size() < 1) continue;

      // eec calculation
      Correlator_Builder cb(charged_constituents, jets[ijet].pt(), jets[ijet].eta(), eec_weight_power);
      cb.make_pairs();
      cb.construct_EEC();
    }

  }
}

void read_csv(const char* inFile = "merged.csv", double proj_rest_e = 10, double targ_lab_e = 100, int targ_species = 0,
    double eec_weight_power = 1, int in_lab_frame = 0)
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
  Double_t Charge, Px, Py, Pz, Mass;

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
      // jet histograms filled on inclusive jet information
      h1d_jet_pt->Fill(jets[ijet].pt());
      h1d_jet_eta->Fill(jets[ijet].eta());
      //cout<<"jet pt:"<<jets[ijet].pt()<<" jet eta:"<<jets[ijet].eta()<<endl;

      // cuts on jet kinematics, require jet_pt >= 5GeV, |jet_eta| <= 2.5
      if (jets[ijet].pt() < 5 || fabs(jets[ijet].eta()) > 2.5) continue;

      // TODO for calculation in nuclear rest frame, boost constituents here

      // cuts on jet constituent kinematics, require consitituents_pt >= 0.5GeV, |consitituents_eta| <= 3.5
      // take only charged constituents for eec calculation
      vector<PseudoJet> constituents = jets[ijet].constituents();
      vector<PseudoJet> charged_constituents;
      for (unsigned iconstit = 0; iconstit < constituents.size(); iconstit++)
      {
        if (constituents[iconstit].pt() < 0.5 || fabs(constituents[iconstit].eta()) > 3.5) continue;

        int il = constituents[iconstit].user_index();
        vector<string> line = content[il];
        Charge = stod(line[2]);
        //cout<<"constituent pt:"<<constituents[iconstit].pt()<<" charge:"<<Charge<<endl;
        if (Charge != 0)
        {
          if (in_lab_frame == 0)
          {
            charged_constituents.push_back(constituents[iconstit]); // keep in lab frame
          }
          else if (in_lab_frame == 1) // boost back charged constituent to nuclear rest frane
          {
            Px = stod(line[3]);
            Py = stod(line[4]);
            Pz = stod(line[5]);
            Mass = stod(line[6]);

            cout<<"this part runs OMGGGGG "<<endl;

            part_rest.SetXYZM(Px, Py, Pz, Mass);
            PseudoJet constituent_rest_frame = PseudoJet(part_rest.Px(),part_rest.Py(),part_rest.Pz(),part_rest.E()); // in original nuclear rest frame
            charged_constituents.push_back(constituent_rest_frame);
          }
        }
      }

      if (charged_constituents.size() < 1) continue;

      // eec calculation
      Correlator_Builder cb(charged_constituents, jets[ijet].pt(), jets[ijet].eta(), eec_weight_power);
      cb.make_pairs();
      cb.construct_EEC();

    }

  }

  fin.close();

}


void eec_hists(const char* inFile = "merged.root", const char* outFile = "hists_eec.root", const int gen_type = 1,
    double proj_rest_e = 2131.56, double targ_lab_e = 100, int targ_species = 0, double eec_weight_power = 1,
    int in_lab_frame = 0)
{
  // proj_rest_e = energy of projectile beam in target rest frame, leave blank if pythia
  // targ_lab_e = energy of target beam in lab frame, leave blank if pythia
  // all energies positive and in GeV units
  // only for e+A collsions, specify A with targ_species, =0 for p, =1 for Au
  // in_lab_frame = 0 for dist12 calculation in lab frame, =1 for dist12 calculation in nuclear rest frame

  cout << "Generator Type: ";
  if (gen_type==0) cout << "Pythia6" << endl;
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
  Double_t lbins_rlsqrtpt[nbins+1];
  binwidth = (log10(xmax) - log10(xmin)) / nbins;
  for (int i = 0; i < nbins+1; i++)
  {
    lbins_rlsqrtpt[i] = TMath::Power(10, log10(xmin) + binwidth * i);
  }

  // histogram definitions
  h1d_jet_pt = new TH1D("h1d_jet_pt","jet pt",800,0,800);
  h1d_jet_pt->Sumw2();
  h1d_jet_eta = new TH1D("h1d_jet_eta", "jet eta",800,-5,5);
  h1d_jet_eta->Sumw2();

  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h1d_jet_eec[ieta][ipt] = new TH1D(Form("h1d_jet_eec_%d_%d", ieta, ipt),"jet eec",50,lbins);
      h1d_jet_eec[ieta][ipt]->Sumw2();

      h1d_jet_eec_rlsqrtpt[ieta][ipt] = new TH1D(Form("h1d_jet_eec_rlsqrtpt_%d_%d", ieta, ipt),"jet eec rlsqrtpt",50,lbins_rlsqrtpt);
      h1d_jet_eec_rlsqrtpt[ieta][ipt]->Sumw2();
    }
  }


  // reads file and fills in jet_constits
  if (gen_type == 0) read_root(inFile, eec_weight_power);
  else read_csv(inFile, proj_rest_e, targ_lab_e, targ_species, eec_weight_power, in_lab_frame);
  cout<<"@kdebug last"<<endl;

  // write out histograms
  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  h1d_jet_pt->Write();
  cout<<"h1d_jet_pt entries:"<<h1d_jet_pt->GetEntries()<<endl;
  h1d_jet_eta->Write();
  cout<<"h1d_jet_eta entries:"<<h1d_jet_eta->GetEntries()<<endl;
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h1d_jet_eec[ieta][ipt]->Write();
      cout<<"h1d_jet_eec_"<<ieta<<"_"<<ipt<<" entries:"<<h1d_jet_eec[ieta][ipt]->GetEntries()<<endl;

      h1d_jet_eec_rlsqrtpt[ieta][ipt]->Write();
      cout<<"h1d_jet_eec_rlsqrtpt_"<<ieta<<"_"<<ipt<<" entries:"<<h1d_jet_eec_rlsqrtpt[ieta][ipt]->GetEntries()<<endl;
    }
  }

}
