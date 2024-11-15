R__LOAD_LIBRARY(libeicsmear);
void access_tree_pt_eta_TH2Ds()
{
  // run with root -l 'access_tree("sim_dir/sim_output.root", "histogram_dir/output.root")'
  // arg: inFile: directory + filename of .root file outputed from sim

  // makes num_species amount of 2D histograms counting multiplicity wrt pT and eta

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  const int num_species = 7;
  std::string dirs[num_species] = {"../ep_10_100/outfiles/", "../eD_18_110/outForPythiaMode/", "../eHe4_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eCa_18_110/outForPythiaMode/", "../eCu_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/"};

  TH2D * pt_vs_eta_kaon;
  TH2D * pt_vs_eta_pion;
  TH2D * pt_vs_eta_proton;

  //loop over files
  TFile *f;
  TFile *fout;
  std::string inFile;
  std::string outFile;
  TTree *tree;
  Int_t nEntries;
  Int_t nParticles(0);
  Int_t id;
  Int_t status;
  Double_t pT;
  Double_t eta;

  //loop over each merged.root file
  for (int i = 0; i < num_species; i++) {
    //Load ROOT File for pythia
    inFile = dirs[i] + "merged.root";
    f = new TFile(inFile.c_str(), "read");

    //Get EICTree Tree
    tree = (TTree*)f->Get("EICTree");
    nEntries = tree->GetEntries();

    //Access event Branch
    tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

    //load new out file
    outFile = dirs[i] + "pt_eta_TH2Ds.root";
    fout = new TFile(outFile.c_str(), "recreate");

    //Initialize new histogram
    pt_vs_eta_kaon = new TH2D("pt_vs_eta_kaon", "particle multiplicity, pt vs eta, kaon", 100, 0, 20, 3, -3, 3);
    pt_vs_eta_kaon->Sumw2();
    pt_vs_eta_pion = new TH2D("pt_vs_eta_pion", "particle multiplicity, pt vs eta, pion", 100, 0, 20, 3, -3, 3);
    pt_vs_eta_pion->Sumw2();
    pt_vs_eta_proton = new TH2D("pt_vs_eta_proton", "particle multiplicity, pt vs eta, proton", 100, 0, 20, 3, -3, 3);
    pt_vs_eta_proton->Sumw2();

    //loop over events
    for (int j = 0; j < nEntries; j++) {

      tree->GetEntry(j);

      //Get Total Number of Particles
      nParticles = event->GetNTracks();

      //Loop Over Each Particle
      for(int k = 0; k < nParticles; k++) {

        particle = event->GetTrack(k);
        status = (Int_t) particle->GetStatus(); //Can also do particle->KS
        id = (Int_t) particle->Id();
        pT = (Double_t) particle->GetPt();
        eta = (Double_t) particle->GetEta();

        if (status == 1) {
          if (id == 321 || id == -321) {
            pt_vs_eta_kaon->Fill(pT, eta);
          } else if (id == 211 || id == -211) {
            pt_vs_eta_pion->Fill(pT, eta);
          } else if (id == 2212 || id == -2212) {
            pt_vs_eta_proton->Fill(pT, eta);
          }
        }

      }
    }

    pt_vs_eta_kaon->Write();
    pt_vs_eta_pion->Write();
    pt_vs_eta_proton->Write();
    fout->Write();
    fout->Close();

  }

}
