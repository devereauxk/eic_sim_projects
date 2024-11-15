R__LOAD_LIBRARY(libeicsmear);
void access_tree(const char* inFile, const char* outFile)
{
  // run with root -l 'access_tree("sim_dir/sim_output.root", "histogram_dir/output.root")'
  // arg: inFile: directory + filename of .root file outputed from sim

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  //Load ROOT File for pythia
  TFile *f = new TFile(inFile, "read");

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");
  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  TH2D* h2d_kaon = new TH2D("h2d_kaon","charged kaon multiplicity",12,-0.5,11.5,12,-0.5,11.5);
  h2d_kaon->Sumw2(); // to handle error propagation correctly later

  TH2D* h2d_pion = new TH2D("h2d_pion","charged pion multiplicity",25,-0.5,24.5,25,-0.5,24.5);
  h2d_pion->Sumw2(); // to handle error propagation correctly later

  TH2D* h2d_proton = new TH2D("h2d_proton","proton/antiproton multiplicity",10, -0.5, 9.5, 10, -0.5, 9.5);
  h2d_proton->Sumw2(); // to handle error propagation correctly later

  TH1D* h1d_kaon = new TH1D("h1d_kaon_total", "charged kaon multiplicity", 12, -0.5, 11.5);
  h2d_pion->Sumw2();

  TH1D* h1d_pion = new TH1D("h1d_pion_total", "charged pion multiplicity", 40, -0.5, 39.5);
  h2d_pion->Sumw2();

  TH1D* h1d_proton = new TH1D("h1d_proton_total","proton/antiproton multiplicity",10,-0.5,9.5);
  h1d_proton->Sumw2(); // to handle error propagation correctly later

  //Define Some Variables
  Int_t nParticles(0);
  Int_t id;
  Int_t status;
  Int_t nPosKaons;
  Int_t nNegKaons;
  Int_t nPosPions;
  Int_t nNegPions;
  Int_t nProtons;
  Int_t nAntiProtons;

  //Loop Over Events
  for(Int_t i = 0; i < nEntries; i++) {

    nPosKaons = 0;
    nNegKaons = 0;
    nPosPions = 0;
    nNegPions = 0;
    nProtons = 0;
    nAntiProtons = 0;

    tree->GetEntry(i);

    //Get Total Number of Particles
    nParticles = event->GetNTracks();

    //Loop Over Each Particle
    for(Int_t j = 0; j < nParticles; j++) {

      particle = event->GetTrack(j);
      status = (Int_t) particle->GetStatus(); //Can also do particle->KS
      id = (Int_t) particle->Id();

      if (status == 1) {
        switch(id) {
          case 321:
            nPosKaons++;
            break;
          case -321:
            nNegKaons++;
            break;
          case 211:
            nPosPions++;
            break;
          case -211:
            nNegPions++;
            break;
          case 2212:
            nProtons++;
            break;
          case -2212:
            nAntiProtons++;
            break;
          default: break;
        }
      }

    }

    h2d_kaon->Fill(nPosKaons, nNegKaons);
    h1d_kaon->Fill(nPosKaons+nNegKaons);
    h2d_pion->Fill(nPosPions, nNegPions);
    h1d_pion->Fill(nPosPions+nNegPions);
    h2d_proton->Fill(nProtons, nAntiProtons);
    h1d_proton->Fill(nProtons+nAntiProtons);

  }

  cout<<"-------------------------------"<<endl;
  cout<<"Write output to root file"<<endl;
  TFile* fout = new TFile(outFile,"recreate");
  h2d_kaon->Write();
  h1d_kaon->Write();
  h2d_pion->Write();
  h1d_pion->Write();
  h2d_proton->Write();
  h1d_proton->Write();
  fout->Write();
  fout->Close();

}
