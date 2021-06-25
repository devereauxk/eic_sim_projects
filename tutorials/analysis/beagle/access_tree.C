R__LOAD_LIBRARY(libeicsmear);
void access_tree(){
  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer
  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer
  //Load ROOT File for pythia
  //TFile *f = new TFile("../pythia/outfiles/ep_10_100_norad.root");
  //Load ROOT File for beagle
  TFile *f = new TFile("../../beagle/outForPythiaMode/eAu.root");
  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");
  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;
  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer
  //Initialize a 2D histogram to store particle momentum vs eta
  TH2D* h2d_pt_vs_eta = new TH2D("h2d_pt_vs_eta","particle pt vs pseudo-rapidity",100,0,20,80,-4,4);
  h2d_pt_vs_eta->Sumw2(); // to handle error propagation correctly later
  //Initialize a 1D histogram to store particle transeverse momentum distribution
  TH1D* h1d_part_pt = new TH1D("h1d_part_pt","particle pt distribution",100,0,20);
  h1d_part_pt->Sumw2(); // to handle error propagation correctly later
  //Define Some Variables
  Float_t Q2(0);
  Int_t nParticles(0);
  Float_t Eta[500];
  Float_t Px[500];
  Int_t Status[500];
  Int_t id[500];
  //Loop Over Events
  for(Int_t i=0;i<10;i++){
    tree->GetEntry(i);
    //Write Out Q2
    Q2 = (Float_t) event->GetQ2(); //Can also do event->QSquared
    printf("For Event %d, Q^2 = %.3f GeV^2!\n",i,Q2);
    //Get Total Number of Particles
    nParticles = event->GetNTracks();
    printf("For Event %d, we have %d particles!\n",i,nParticles);
    //Loop Over Each Particle
    for(Int_t j=0;j<nParticles;j++){
      particle = event->GetTrack(j);
      Eta[j] =  (Float_t) particle->GetEta();
      Px[j] = (Float_t) particle->GetPx();
      Status[j] = (Int_t) particle->GetStatus(); //Can also do particle->KS
      id[j] = (Int_t) particle->Id();
      if(Status[j]==1){
      	printf("For Event %d, particle %d Eta = %.3f!\n",i,j,Eta[j]);
      	printf("For Event %d, particle %d Px = %.3f GeV/c!\n",i,j,Px[j]);
      	printf("For Event %d, particle %d Status = %d!\n",i,j,Status[j]);
      	printf("For Event %d, particle %d Id = %d!\n",i,j,id[j]); 
        h2d_pt_vs_eta->Fill(particle->GetPt(),particle->GetEta());
        h1d_part_pt->Fill(particle->GetPt());
      }
    }
  }
  cout<<"-------------------------------"<<endl;
  cout<<"Write output to root file"<<endl;
  TFile* fout = new TFile("output.root","recreate");
  h1d_part_pt->Write();
  h2d_pt_vs_eta->Write();
  fout->Write();
  fout->Close();
}
