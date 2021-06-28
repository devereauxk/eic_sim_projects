R__LOAD_LIBRARY(libeicsmear);
void access_tree(const char* inFile, const char* outDir)
{
  // arg: inFile: directory + filename of .root file outputed from sim
  //arg: outDir: directory to output output.root file containing histogram data

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

  TH1D* h1d_kaon = new TH1D("h1d_kaon","charged kaon multiplicity",100,0,20,0,100);
  h1d_kaon->Sumw2(); // to handle error propagation correctly later

  TH1D* h1d_pion = new TH1D("h1d_pion","charged pion multiplicity",100,0,20,0,100);
  h1d_pion->Sumw2(); // to handle error propagation correctly later

  TH1D* h1d_proton = new TH1D("h1d_proton","proton multiplicity",100,0,20,0,100);
  h1d_proton->Sumw2(); // to handle error propagation correctly later

  //Define Some Variables
  Int_t nParticles(0);
  Int_t status;
  Int_t id;
  Int_t nChargedKaons;
  Int_t nChargedPions;
  Int_t nProtons;

  //Loop Over Events
  for(Int_t i = 0; i < 10; i++) {

    nChargedKaons = 0;
    nChargedKaons = 0;
    nProtons = 0;

    tree->GetEntry(i);

    //Get Total Number of Particles
    nParticles = event->GetNTracks();
    printf("For Event %d, we have %d particles!\n",i,nParticles);

    //Loop Over Each Particle
    for(Int_t j = 0; j < nParticles; j++) {

      particle = event->GetTrack(j);
      status = (Int_t) particle->GetStatus(); //Can also do particle->KS
      id = (Int_t) particle->Id();

      if (id == 321 || id == -321) { // if charged kaon (+ or -, resp.)
        nChargedKaons++;
      } else if (id == 211 || id == -211) { // if charged pion (+ or -, resp.)
        nChargedPions++;
      } else if (id == 2212) { // if proton
        nProtons++;
      }

      if(Status[j]==1){
      	printf("For Event %d, particle %d Eta = %.3f!\n",i,j,Eta[j]);
      	printf("For Event %d, particle %d Px = %.3f GeV/c!\n",i,j,Px[j]);
      	printf("For Event %d, particle %d Status = %d!\n",i,j,Status[j]);
      	printf("For Event %d, particle %d Id = %d!\n",i,j,id[j]);
        h2d_pt_vs_eta->Fill(particle->GetPt(),particle->GetEta());
        h1d_part_pt->Fill(particle->GetPt());
      }

      h1d_kaon->Fill(nChargedKaons);
      h1d_pion->Fill(nChargedPions);
      h1d_proton->Fill(nProtons);
    }
  }

  cout<<"-------------------------------"<<endl;
  cout<<"Write output to root file"<<endl;
  TFile* fout = new TFile(outDir+"output.root","recreate");
  h1d_kaon->Write();
  h1d_pion->Write();
  h1d_proton->Write();
  fout->Write();
  fout->Close();

}
