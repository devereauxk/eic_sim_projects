R__LOAD_LIBRARY(libeicsmear);
void access_tree_mult_eta_binned()
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

  const int etabin = 3;
  const double eta_lo[etabin] = {-3, -1, 1};
  const double eta_hi[etabin] = {-1, 1, 3};
  TH2D * kaon[etabin];
  TH2D * pion[etabin];
  TH2D * proton[etabin];

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
  Double_t eta;
  Int_t nPosKaons[etabin] = {0};
  Int_t nNegKaons[etabin] = {0};
  Int_t nPosPions[etabin] = {0};
  Int_t nNegPions[etabin] = {0};
  Int_t nProtons[etabin] = {0};
  Int_t nAntiProtons[etabin] = {0};

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
    outFile = dirs[i] + "mult_eta_binned.root";
    fout = new TFile(outFile.c_str(), "recreate");

    //Initialize new histograms
    for (int ieta = 0; ieta < etabin; ieta++) {
      kaon[ieta] = new TH2D(Form("h2d_kaon_%d", ieta),"charged kaon multiplicity",12,-0.5,11.5,12,-0.5,11.5);
      kaon[ieta]->Sumw2();
      pion[ieta] = new TH2D(Form("h2d_pion_%d", ieta),"charged pion multiplicity",25,-0.5,24.5,25,-0.5,24.5);
      pion[ieta]->Sumw2();
      proton[ieta] = new TH2D(Form("h2d_proton_%d", ieta),"proton/antiproton multiplicity",15, -0.5, 14.5, 15, -0.5, 14.5);
      proton[ieta]->Sumw2();
    }

    //loop over events
    for (int j = 0; j < nEntries; j++) {

      Int_t nPosKaons[etabin] = {0};
      Int_t nNegKaons[etabin] = {0};
      Int_t nPosPions[etabin] = {0};
      Int_t nNegPions[etabin] = {0};
      Int_t nProtons[etabin] = {0};
      Int_t nAntiProtons[etabin] = {0};

      tree->GetEntry(j);

      //Get Total Number of Particles
      nParticles = event->GetNTracks();

      //Loop Over Each Particle
      for(int k = 0; k < nParticles; k++) {

        particle = event->GetTrack(k);
        status = (Int_t) particle->GetStatus(); //Can also do particle->KS
        id = (Int_t) particle->Id();
        eta = (Double_t) particle->GetEta();

        if (status == 1) {
          for (int ieta = 0; ieta < etabin; ieta++) {
            if (eta_lo[ieta] <= eta && eta <= eta_hi[ieta]) {
              switch(id) {
                case 321:
                  nPosKaons[ieta]++;
                  break;
                case -321:
                  nNegKaons[ieta]++;
                  break;
                case 211:
                  nPosPions[ieta]++;
                  break;
                case -211:
                  nNegPions[ieta]++;
                  break;
                case 2212:
                  nProtons[ieta]++;
                  break;
                case -2212:
                  nAntiProtons[ieta]++;
                  break;
                default: break;
              }
            }
          }
        }
      }

      for (int ieta = 0; ieta < etabin; ieta++) {
        kaon[ieta]->Fill(nPosKaons[ieta], nNegKaons[ieta]);
        pion[ieta]->Fill(nPosPions[ieta], nNegKaons[ieta]);
        proton[ieta]->Fill(nProtons[ieta], nAntiProtons[ieta]);
      }

    }

    for (int ieta = 0; ieta < etabin; ieta++) {
      kaon[ieta]->Write();
      pion[ieta]->Write();
      proton[ieta]->Write();
    }

    fout->Write();
    fout->Close();

  }

}
