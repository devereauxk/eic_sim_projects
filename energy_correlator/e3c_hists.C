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

const int etabin = 6; // inclusive on last bin, inclusive on lower limit, exclusive on upper
static double eta_lo[etabin] = {-3.5, -1, 1, -3.5, -1, 0};
static double eta_hi[etabin] = {-1, 1, 3.5, 3.5, 0, 1};

TH1D* h1d_jet_eec[etabin][ptbin] = {};
TH1D* h1d_jet_eec_rlsqrtpt[etabin][ptbin] = {};
TH1D* h1d_jet_pt[etabin] = {};
TH1D* h1d_jet_eta = NULL;
TH1D* h1d_jet_multiplicity[etabin][ptbin] = {};
TH1D* h1d_jet_multiplicity_charged[etabin][ptbin] = {};
TH2D* h2d_Q2_x[etabin][ptbin] = {};

// debugging histograms
TH1D* h1d_part_pt[etabin] = {};
TH1D* h1d_part_eta[ptbin] = {};
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
      // produces mult X mult matrix with elements being the
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
          for (int k = 0; k < mult; k++)
          {

          }
          double distij = pair_list[i][j];
          double distjk = pair_list[j][k];
          double distik = pair_list[i][k];

          double Rl = Max(distij, Max(distjk, distik));
          double Rs = Min(distij, Min(distjk, distik));
          double Rm = Max(Min(distij, distjk) + Max(Min(distij, distik) + Min(distjk, distik)));

          double Xi = Rs / Rm;
          double phi = ASin(Sqrt)

          double eec_weight = pow((particle_list[i].pt() * particle_list[j].pt() * particle_list[k].pt()) / pow(jet_pt, 3), weight_pow);

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

void read_root(const char* inFile = "merged.root", double eec_weight_power = 1, int gen_type = 0,
    int boost = 0, double proj_rest_e = 10, double targ_lab_e = 100, int targ_species = 0)
{
  //Event Class
  erhic::EventHepMC *event(NULL); //pythia8
  if (gen_type==0) erhic::EventPythia *event(NULL); // pythia6

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nevt = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nevt<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  // boost calculation (used iff boost == 1)
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

  Double_t Px, Py, Pz, Mass, Q2, xB, charge, Pt, Eta, Mult;

  int total_jets = 0;

  //Loop Over Events
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    tree->GetEntry(ievt);

    if (ievt%1000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    // Q2-cut
    if (event->GetQ2() < 10) continue;

    Q2 = event->GetQ2();
    xB = event->GetX();

    Mult = 0;

    // particle enumeration, addition to jet reco setup, and total pt calculation
    erhic::ParticleMC* particle;
    vector<PseudoJet> jet_constits;
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      particle = event->GetTrack(ipart);

      if (particle->GetStatus()==1) Mult++;

      // get particle kinematics
      Px = particle->GetPx();
      Py = particle->GetPy();
      Pz = particle->GetPz();
      Mass = particle->GetM();
      part.SetXYZM(Px, Py, Pz, Mass);

      // apply boost, if required. boosted according to proj_rest_e, targ_lab_e, targ_species
      if (boost == 1) part.Boost(boost_vec);

      // debug histograms
      if (particle->GetStatus()==1)
      {
        Pt = part.Pt();
        Eta = part.Eta();
        for (int ipt = 0; ipt < ptbin; ipt++)
        {
          if (Pt >= pt_lo[ipt] && Pt < pt_hi[ipt]) h1d_part_eta[ipt]->Fill(Eta);
        }
        for (int ieta = 0; ieta < etabin; ieta++)
        {
          if (Eta >= eta_lo[ieta] && Eta < eta_hi[ieta]) h1d_part_pt[ieta]->Fill(Pt);
        }
      }

      // use all fsp particles w/ < 3.5 eta, not including scattered electron, for jet reconstruction
      if (particle->GetStatus()==1 && fabs(particle->GetEta())<3.5 && particle->Id()!=11)
      {
        PseudoJet constit = PseudoJet(part.Px(),part.Py(),part.Pz(),part.E());
        constit.set_user_index(ipart);
        jet_constits.push_back(constit);
      }

    }
    h1d_part_mult->Fill(Mult);

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
        //if (charge_flag != 0) charged_constituents.push_back(constituents[iconstit]);
      }

      // jet histograms filled on inclusive jet information
      for (int ieta = 0; ieta < etabin; ieta++)
      {
        if (jets[ijet].eta() >= eta_lo[ieta] && jets[ijet].eta() < eta_hi[ieta])
        {
          h1d_jet_pt[ieta]->Fill(jets[ijet].pt());

          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            if (jets[ijet].pt() >= pt_lo[ipt] && jets[ijet].pt() < pt_hi[ipt])
            {
              h1d_jet_multiplicity[ieta][ipt]->Fill(constituents.size());
              h1d_jet_multiplicity_charged[ieta][ipt]->Fill(charged_constituents.size());

              h2d_Q2_x[ieta][ipt]->Fill(xB, Q2);
            }
          }
        }

      }
      h1d_jet_eta->Fill(jets[ijet].eta());

      if (charged_constituents.size() < 1) continue;

      // eec calculation
      Correlator_Builder cb(charged_constituents, jets[ijet].pt(), jets[ijet].eta(), eec_weight_power);
      cb.make_pairs();
      cb.construct_EEC();
    }

  }

  cout<<"total num jets = "<<total_jets<<endl;

}

void read_csv(const char* inFile = "merged.csv", int boost = 1, double proj_rest_e = 10, double targ_lab_e = 100, int targ_species = 0,
    double eec_weight_power = 1, int calc_Q2x = 0)
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

  int total_jets = 0;

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

    vector<PseudoJet> jet_constits;

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
      if (calc_Q2x == 1)
      {
        Q2 = stod(line[7]);
        xB = stod(line[8]);
      }

      //cout<<iline<<" "<<Id<<" "<<Charge<<" "<<Px<<" "<<Py<<" "<<Pz<<" "<<Mass<<endl;

      // apply boost to particle (boost it into lab frame)
      part.SetXYZM(Px, Py, Pz, Mass);
      if (boost == 1) part.Boost(boost_vec);

      // debug histograms
      Pt = part.Pt();
      Eta = part.Eta();
      for (int ipt = 0; ipt < ptbin; ipt++)
      {
        if (Pt >= pt_lo[ipt] && Pt < pt_hi[ipt]) h1d_part_eta[ipt]->Fill(Eta);
      }
      for (int ieta = 0; ieta < etabin; ieta++)
      {
        if (Eta >= eta_lo[ieta] && Eta < eta_hi[ieta]) h1d_part_pt[ieta]->Fill(Pt);
      }

      // use all fsp particles w/ < 3.5 eta, not including scattered electron, for jet reconstruction
      //cout<<"part_lab eta:"<<part_lab.Eta()<<endl;
      if (fabs(part.Eta())<3.5 && Id!=11)
      {
        PseudoJet constit = PseudoJet(part.Px(),part.Py(),part.Pz(),part.E());
        constit.set_user_index(iline);
        jet_constits.push_back(constit);
      }

      iline++;
    }
    h1d_part_mult->Fill(Mult);

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
        if (Charge != 0) charged_constituents.push_back(constituents[iconstit]);
      }

      // jet histograms filled on inclusive jet information
      for (int ieta = 0; ieta < etabin; ieta++)
      {
        if (jets[ijet].eta() >= eta_lo[ieta] && jets[ijet].eta() < eta_hi[ieta])
        {
          h1d_jet_pt[ieta]->Fill(jets[ijet].pt());

          for (int ipt = 0; ipt < ptbin; ipt++)
          {
            if (jets[ijet].pt() >= pt_lo[ipt] && jets[ijet].pt() < pt_hi[ipt])
            {
              h1d_jet_multiplicity[ieta][ipt]->Fill(constituents.size());
              h1d_jet_multiplicity_charged[ieta][ipt]->Fill(charged_constituents.size());

              if (calc_Q2x == 1) h2d_Q2_x[ieta][ipt]->Fill(xB, Q2);
            }
          }
        }

      }
      h1d_jet_eta->Fill(jets[ijet].eta());

      if (charged_constituents.size() < 1) continue;

      // eec calculation
      Correlator_Builder cb(charged_constituents, jets[ijet].pt(), jets[ijet].eta(), eec_weight_power);
      cb.make_pairs();
      cb.construct_EEC();

    }

  }

  cout<<"total num jets = "<<total_jets<<endl;

  fin.close();

}


void eec_hists(const char* inFile = "merged.root", const char* outFile = "hists_eec.root", const int gen_type = 1,
    double proj_rest_e = 2131.56, double targ_lab_e = 100, int targ_species = 0, double eec_weight_power = 1,
    int boost = 0, int calc_Q2x = 0)
{
  // proj_rest_e = energy of projectile beam in target rest frame, leave blank if pythia
  // targ_lab_e = energy of target beam in lab frame, leave blank if pythia
  // all energies positive and in GeV units
  // only for e+A collsions, specify A with targ_species, =0 for p, =1 for Au
  // gen_type = 0 for pythia6, =-1 for pyhtia8, =1 or anything for eHIJING (DIFFERENT FROM Q2_x.C settings)
  // boost = 0, do not apply boost =1 apply boost according to targ_lab_e and targ_species
  // calc_Q2x = 1 fills Q2x histogram (eHIJING must have Q2 and x printout), =0 doesn't filled

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
  Double_t lbins_rlsqrtpt[nbins+1];
  binwidth = (log10(xmax) - log10(xmin)) / nbins;
  for (int i = 0; i < nbins+1; i++)
  {
    lbins_rlsqrtpt[i] = TMath::Power(10, log10(xmin) + binwidth * i);
  }

  // compute log bins for Q2-x
  // xbins correspond to xB
  xmin = 1E-6;
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
  double ymax = 1E4;
  int ynbins = 50;
  Double_t ylbins[ynbins+1];
  double ybinwidth = (log10(ymax) - log10(ymin)) / ynbins;
  for (int i = 0; i < ynbins+1; i++)
  {
    ylbins[i] = TMath::Power(10, log10(ymin) + ybinwidth * i);
  }

  // histogram definitions
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    h1d_jet_pt[ieta] = new TH1D(Form("h1d_jet_pt_%d", ieta),"jet pt",800,0,800);
    h1d_jet_pt[ieta]->Sumw2();

    h1d_part_pt[ieta] = new TH1D(Form("h1d_part_pt_%d", ieta),"particle pt",400,0,100);
    h1d_part_pt[ieta]->Sumw2();
  }
  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    h1d_part_eta[ipt] = new TH1D(Form("h1d_part_eta_%d", ipt),"particle eta",1600,-15,5);
    h1d_part_eta[ipt]->Sumw2();
  }
  h1d_part_mult = new TH1D("h1d_part_mult", "event particle multiplicity",200,0,200);
  h1d_part_mult->Sumw2();

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

      h1d_jet_multiplicity[ieta][ipt] = new TH1D(Form("h1d_jet_multiplicity_%d_%d", ieta, ipt), "jet multiplicity", 50, 0, 50);
      h1d_jet_multiplicity[ieta][ipt]->Sumw2();

      h1d_jet_multiplicity_charged[ieta][ipt] = new TH1D(Form("h1d_jet_multiplicity_charged_%d_%d", ieta, ipt), "jet charged multiplicity", 50, 0, 50);
      h1d_jet_multiplicity_charged[ieta][ipt]->Sumw2();
    }
  }
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    for (int ipt = 0; ipt < ptbin; ipt++)
    {
      h2d_Q2_x[ieta][ipt] = new TH2D(Form("h2d_Q2_x_%d_%d", ieta, ipt),"Q2_x",50,xlbins,50,ylbins);
      h2d_Q2_x[ieta][ipt]->Sumw2();
    }
  }


  // reads file and fills in jet_constits
  if (gen_type == 0 || gen_type == -1) read_root(inFile, eec_weight_power, gen_type, boost, proj_rest_e, targ_lab_e, targ_species); // pythia6 (EventPythia) or pythia8 (EventHepMC), boosts acording to boost variable
  else read_csv(inFile, boost, proj_rest_e, targ_lab_e, targ_species, eec_weight_power, calc_Q2x); // eHIJING, assumes target frame and boosts to EIC
  cout<<"@kdebug last"<<endl;

  // write out histograms
  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  for (int ieta = 0; ieta < etabin; ieta++)
  {
    h1d_jet_pt[ieta]->Write();
    cout<<"h1d_jet_pt_"<<ieta<<" entries:"<<h1d_jet_pt[ieta]->GetEntries()<<endl;

    h1d_part_pt[ieta]->Write();
  }
  for (int ipt = 0; ipt < ptbin; ipt++)
  {
    h1d_part_eta[ipt]->Write();
  }
  h1d_part_mult->Write();
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

      h1d_jet_multiplicity[ieta][ipt]->Write();
      cout<<"h1d_jet_multiplicity_"<<ieta<<"_"<<ipt<<" entries:"<<h1d_jet_multiplicity[ieta][ipt]->GetEntries()<<endl;
      h1d_jet_multiplicity_charged[ieta][ipt]->Write();
      cout<<"h1d_jet_multiplicity_charged_"<<ieta<<"_"<<ipt<<" entries:"<<h1d_jet_multiplicity_charged[ieta][ipt]->GetEntries()<<endl;

      h2d_Q2_x[ieta][ipt]->Write();
      cout<<"h2d_Q2_x_"<<ieta<<"_"<<ipt<<" entries:"<<h2d_Q2_x[ieta][ipt]->GetEntries()<<endl;
    }
  }

}
