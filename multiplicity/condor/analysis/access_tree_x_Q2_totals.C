R__LOAD_LIBRARY(libeicsmear);
void access_tree_x_Q2_totals()
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

  TH2D * h2d_x_vs_q2;

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
  Double_t q2;

  //loop over each merged.root file
  for (int i = 0; i < num_species; i++) {
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
    outFile = dirs[i] + "x_Q2_total_TH2D.root";
    fout = new TFile(outFile.c_str(), "recreate");

    //Initialize new histogram
    h2d_x_vs_q2 = new TH2D("h2d_x_vs_q2","event x vs Q^2",100,-5,1,100,10,100);
    h2d_x_vs_q2->Sumw2();

    //loop over events
    for (int j = 0; j < nEntries; j++) {

      tree->GetEntry(j);

      //Get Total Number of Particles
      nParticles = event->GetNTracks();

      x = event->GetX();
      q2 = event->GetQ2();

      if (event->GetProcess() == 99) {
        h2d_x_vs_q2->Fill(log10(x), q2); // use actual log scalling
      }

    }

    h2d_x_vs_q2->Write();
    fout->Write();
    fout->Close();

  }

}
