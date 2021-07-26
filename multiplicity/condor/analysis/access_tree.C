R__LOAD_LIBRARY(libeicsmear);
void access_tree(const char* outFile)
{
  // run with root -l 'access_tree("sim_dir/sim_output.root", "histogram_dir/output.root")'
  // arg: inFile: directory + filename of .root file outputed from sim

  // makes num_species amount of 2D histograms counting multiplicity wrt pT and eta

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  TH3D* kaon_pos = new TH2D("kaon_pos","kaon multiplicity vs pT and eta and A",18,-4.5,4.5,30,-0.5,14.5,220,-0.5,219.5);
  kaon_pos->Sumw2(); // to handle error propagation correctly later



  int num_species = 7;
  std::string dirs[num_species] = {"../ep_10_100/outfiles/merged.root", "../eD_18_110/outForPythiaMode/merged.root", "../eHe4_18_110/outForPythiaMode/merged.root", "../eC_18_110/outForPythiaMode/merged.root", "../eCa_18_110/outForPythiaMode/merged.root", "../eCu_18_110/outForPythiaMode/merged.root", "../eAu_18_110/outForPythiaMode/merged.root"};
  Double_t atomic_weight[num_species] = {1, 2, 4, 12, 40, 63, 197};

  //loop over files
  TFile *f;
  TTree *tree;
  Int_t nEntries;
  Int_t nParticles(0);
  Int_t id;
  Int_t status;
  Double_t pT;
  Double_t eta;
  Double_t a_number;
  //loop over each merged.root file
  for (int i = 0; i < num_species; i++) {
    //Load ROOT File for pythia
    f = new TFile(dirs[i].c_str(), "read");

    //Get EICTree Tree
    tree = (TTree*)f->Get("EICTree");
    nEntries = tree->GetEntries();

    //Access event Branch
    tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

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
        a_number = (Double_t) atomic_weight[i];

        if (status == 1) {
          switch(id) {
            case 321:
              kaon_pos->Fill(pT, eta, a_number);
              break;
            case -321:
              //nNegKaons++;
              break;
            case 211:
              //nPosPions++;
              break;
            case -211:
              //nNegPions++;
              break;
            case 2212:
              //nProtons++;
              break;
            case -2212:
              //nAntiProtons++;
              break;
            default: break;
          }
        }

      }
    }

  }

  cout<<"-------------------------------"<<endl;
  cout<<"Write output to root file"<<endl;
  TFile* fout = new TFile(outFile,"recreate");
  kaon_pos->Write();
  fout->Write();
  fout->Close();

}
