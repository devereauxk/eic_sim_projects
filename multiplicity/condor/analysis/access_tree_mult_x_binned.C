R__LOAD_LIBRARY(libeicsmear);
void access_tree_mult_x_binned()
{
  // run with root -l 'access_tree("sim_dir/sim_output.root", "histogram_dir/output.root")'
  // arg: inFile: directory + filename of .root file outputed from sim

  // makes num_species amount of 2D histograms counting multiplicity wrt pT and eta

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  const int num_species = 7;
  std::string dirs[num_species] = {"../eD_18_110/outForPythiaMode/",  "../eHe4_18_110/outForPythiaMode/", "../eC_18_110/outForPythiaMode/", "../eCa_18_110/outForPythiaMode/", "../eCu_18_110/outForPythiaMode/", "../eAu_18_110/outForPythiaMode/", "../ep_10_100/outfiles/"};
  std::string inFileNames[num_species] = {"merged.root", "merged.root", "merged.root", "merged.root", "merged.root", "merged.root", "merged.root"};

  const int xbin = 5;
  const double x_lo[xbin] = {1/100000, 1/10000, 1/1000, 1/100, 1/10};
  const double x_hi[xbin] = {1/10000, 1/1000, 1/100, 1/10, 1};
  TH2D * kaon[xbin];
  TH2D * pion[xbin];
  TH2D * proton[xbin];

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
  Double_t x;
  Int_t nPosKaons[xbin] = {0};
  Int_t nNegKaons[xbin] = {0};
  Int_t nPosPions[xbin] = {0};
  Int_t nNegPions[xbin] = {0};
  Int_t nProtons[xbin] = {0};
  Int_t nAntiProtons[xbin] = {0};

  //loop over each merged.root file
  for (Int_t i = 0; i < num_species; i++) {
    cout<<dirs[i]<<endl;
    //Load ROOT File for pythia
    inFile = dirs[i] + inFileNames[i];
    f = new TFile(inFile.c_str(), "read");

    //Get EICTree Tree
    tree = (TTree*)f->Get("EICTree");
    nEntries = tree->GetEntries();

    //Access event Branch
    tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

    //load new out file
    outFile = dirs[i] + "mult_x_binned.root";
    fout = new TFile(outFile.c_str(), "recreate");

    //Initialize new histograms
    for (Int_t ix = 0; ix < xbin; ix++) {
      kaon[ix] = new TH2D(Form("h2d_kaon_%d", ix),"charged kaon multiplicity",12,-0.5,11.5,12,-0.5,11.5);
      kaon[ix]->Sumw2();
      pion[ix] = new TH2D(Form("h2d_pion_%d", ix),"charged pion multiplicity",25,-0.5,24.5,25,-0.5,24.5);
      pion[ix]->Sumw2();
      proton[ix] = new TH2D(Form("h2d_proton_%d", ix),"proton/antiproton multiplicity",15, -0.5, 14.5, 15, -0.5, 14.5);
      proton[ix]->Sumw2();
    }

    //loop over events
    for (Int_t j = 0; j < nEntries; j++) {

      Int_t nPosKaons[xbin] = {0};
      Int_t nNegKaons[xbin] = {0};
      Int_t nPosPions[xbin] = {0};
      Int_t nNegPions[xbin] = {0};
      Int_t nProtons[xbin] = {0};
      Int_t nAntiProtons[xbin] = {0};

      tree->GetEntry(j);

      //Get Total Number of Particles
      nParticles = event->GetNTracks();

      x = event->GetX();

      if (event->GetProcess() == 99) {
        for (Int_t ix = 0; ix < xbin; ix++) {
          if (x_lo[ix] <= x && x <= x_hi[ix]) {

            //Loop Over Each Particle
            for(Int_t k = 0; k < nParticles; k++) {

              particle = event->GetTrack(k);
              status = (Int_t) particle->GetStatus(); //Can also do particle->KS
              id = (Int_t) particle->Id();

              if (status == 1) {
                switch(id) {
                  case 321:
                    nPosKaons[ix]++;
                    break;
                  case -321:
                    nNegKaons[ix]++;
                    break;
                  case 211:
                    nPosPions[ix]++;
                    break;
                  case -211:
                    nNegPions[ix]++;
                    break;
                  case 2212:
                    nProtons[ix]++;
                    break;
                  case -2212:
                    nAntiProtons[ix]++;
                    break;
                  default: break;
                }
              }
            }
          }
          break;
          
        }
      }

      for (Int_t ix = 0; ix < xbin; ix++) {
        kaon[ix]->Fill(nPosKaons[ix], nNegKaons[ix]);
        pion[ix]->Fill(nPosPions[ix], nNegPions[ix]);
        proton[ix]->Fill(nProtons[ix], nAntiProtons[ix]);
      }

    }

    for (Int_t ix = 0; ix < xbin; ix++) {
      kaon[ix]->Write();
      pion[ix]->Write();
      proton[ix]->Write();
    }

    fout->Write();
    fout->Close();

  }

}
