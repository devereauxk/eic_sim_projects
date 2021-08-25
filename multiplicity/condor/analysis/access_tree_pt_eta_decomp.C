R__LOAD_LIBRARY(libeicsmear);
void access_tree_pt_eta_decomp()
{
  // run with root -l 'access_tree("sim_dir/sim_output.root", "histogram_dir/output.root")'
  // arg: inFile: directory + filename of .root file outputed from sim

  // makes num_species amount of 2D histograms counting multiplicity wrt pT and eta

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  const int num_species = 4;
  std::string dirs[num_species] = {"../ep_10_100/outfiles/", "../eD_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/"};
  std::string inFileNames[num_species] = {"merged.root", "merged.root", "merged.root", "merged.root"};

  TH2D * kaon_pos;
  TH2D * kaon_neg;
  TH2D * pion_pos;
  TH2D * pion_neg;
  TH2D * proton_pos;
  TH2D * proton_neg;

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
    inFile = dirs[i] + inFileNames[i];
    f = new TFile(inFile.c_str(), "read");

    //Get EICTree Tree
    tree = (TTree*)f->Get("EICTree");
    nEntries = tree->GetEntries();

    //Access event Branch
    tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

    //load new out file
    outFile = dirs[i] + "pt_eta_decomp.root";
    fout = new TFile(outFile.c_str(), "recreate");

    //Initialize new histogram
    kaon_pos = new TH2D("kaon_pos","particle pt vs pseudo-rapidity, kaon pos",100,0,20,80,-10,10);
    kaon_pos->Sumw2();
    kaon_neg = new TH2D("kaon_neg","particle pt vs pseudo-rapidity, kaon neg",100,0,20,80,-10,10);
    kaon_neg->Sumw2();
    pion_pos = new TH2D("pion_pos","particle pt vs pseudo-rapidity, pion pos",100,0,20,80,-10,10);
    pion_pos->Sumw2();
    pion_neg = new TH2D("pion_neg","particle pt vs pseudo-rapidity, pion neg",100,0,20,80,-10,10);
    pion_neg->Sumw2();
    proton_pos = new TH2D("proton_pos","particle pt vs pseudo-rapidity, proton pos",100,0,20,80,-10,10);
    proton_pos->Sumw2();
    proton_neg = new TH2D("proton_neg","particle pt vs pseudo-rapidity, proton neg",100,0,20,80,-10,10);
    proton_neg->Sumw2();

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
          switch(id) {
            case 321:
              kaon_pos->Fill(pT, eta);
              break;
            case -321:
              kaon_neg->Fill(pT, eta);
              break;
            case 211:
              pion_pos->Fill(pT, eta);
              break;
            case -211:
              pion_neg->Fill(pT, eta);
              break;
            case 2212:
              proton_pos->Fill(pT, eta);
              break;
            case -2212:
              proton_neg->Fill(pT, eta);
              break;
            default: break;
          }
        }
      }
    }

    kaon_pos->Write();
    kaon_neg->Write();
    pion_pos->Write();
    pion_neg->Write();
    proton_pos->Write();
    proton_neg->Write();
    fout->Write();
    fout->Close();

  }

}
