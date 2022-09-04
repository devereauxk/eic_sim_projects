R__LOAD_LIBRARY(libeicsmear);

using namespace std;

TH1D* h1d_nevt = NULL;

// Cross Section Bins
const int Q2bin = 7;
const double Q2_mid[Q2bin] = {12, 20, 35, 60, 120, 200, 300};
double Q2_lo[Q2bin] = {0};
double Q2_hi[Q2bin] = {0};

const int xbin = 15;
const double x_min = 1E-3;
const double x_max = 1E0;
double x_lo[xbin] = {0};
double x_hi[xbin] = {0};
double edgesx[xbin+1] = {0};

TH1D* h1d_x_e[Q2bin] = {0}; // inclusive events
TH1D* h1d_x_c[Q2bin] = {0}; // events with charm
TH1D* h1d_x_s[Q2bin] = {0}; // events with strange

const int verbosity = 0;

void set_Q2_binning()
{
  float size = 0.1;
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    Q2_lo[iQ2] = (1-size)*Q2_mid[iQ2];
    Q2_hi[iQ2] = (1+size)*Q2_mid[iQ2];
  }
}

void set_x_binning()
{
  double log_bw = (log10(x_max) - log10(x_min))/xbin;
  double log_lo,log_hi;

  for(int ibin = 0; ibin < xbin; ibin++)
  {
    log_lo = log10(x_min) + (ibin*log_bw);
    log_hi = log10(x_min) + (ibin+1)*log_bw;
    x_lo[ibin] = pow(10,log_lo);
    x_hi[ibin] = pow(10,log_hi);

    edgesx[ibin] = x_lo[ibin];
    if (ibin==xbin-1) edgesx[ibin+1] = x_hi[ibin];
  }
}

bool event_w_charm(erhic::EventPythia* event, int gen_type)
{
  if (gen_type==0)
  { // Pythia 6
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=21 && !flag_search_group ) continue;
      if ( part->GetStatus()==21 )
      { // entering into the group of particles with KS=21
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=21 && flag_search_group ) break;
    }
  }
  else
  { // BeAGLE
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = event->GetNTracks()-1; ipart >=0; ipart--)
    { // faster looping backwards
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=3 && !flag_search_group ) continue;
      if ( part->GetStatus()==3 )
      { // entering into the group of particles with KS=3
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=3 && flag_search_group ) break;
    }
  }

  return false;
}

bool event_w_strange(erhic::EventPythia* event, int gen_type)
{
  if (gen_type==0)
  { // Pythia 6
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=21 && !flag_search_group ) continue;
      if ( part->GetStatus()==21 )
      { // entering into the group of particles with KS=21
        flag_search_group = true;
        if ( part->Id()==3 || part->Id()==-3 ) return true;
      }
      if ( part->GetStatus()!=21 && flag_search_group ) break;
    }
  }
  else
  { // BeAGLE
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = event->GetNTracks()-1; ipart >=0; ipart--)
    { // faster looping backwards
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=3 && !flag_search_group ) continue;
      if ( part->GetStatus()==3 )
      { // entering into the group of particles with KS=3
        flag_search_group = true;
        if ( part->Id()==3 || part->Id()==-3 ) return true;
      }
      if ( part->GetStatus()!=3 && flag_search_group ) break;
    }
  }

  return false;
}

void EventCounts_EIC(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hists_eventcounts_ep.root", int nevt = 0, int gen_type = 0)
{
  cout << "Generator Type: ";
  if (gen_type==0) cout << "Pythia6" << endl;
  else cout << "BeAGLE" << endl;

  //Event Class
  erhic::EventPythia *event(NULL);

  TFile *f = new TFile(inFile);

  h1d_nevt = new TH1D("h1d_nevt","# of events",1,0.5,1.5); // 1: inclusive events
  h1d_nevt->Sumw2();

  set_Q2_binning();
  set_x_binning();

  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    h1d_x_e[iQ2] = new TH1D(Form("h1d_x_e_Q2_%d",iQ2),"N^{inclusive}_{evt} vs x_{true}",xbin,edgesx);
    h1d_x_e[iQ2]->Sumw2();

    h1d_x_c[iQ2] = new TH1D(Form("h1d_x_c_Q2_%d",iQ2),"N^{charm}_{evt} vs x_{true}",xbin,edgesx);
    h1d_x_c[iQ2]->Sumw2();

    h1d_x_s[iQ2] = new TH1D(Form("h1d_x_s_Q2_%d",iQ2),"N^{strange}_{evt} vs x_{true}",xbin,edgesx);
    h1d_x_s[iQ2]->Sumw2();
  }

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    h1d_nevt->Fill(1);

    tree->GetEntry(ievt);

    if (ievt%1000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    //Write Out Q2
    double Q2 = event->GetQ2(); //Can also do event->QSquared
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);

    int iQ2bin = -9999;
    for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
    {
      if (Q2>=Q2_lo[iQ2] && Q2<Q2_hi[iQ2]) iQ2bin = iQ2;
    }

    if (iQ2bin<0) continue;

    bool flag_direct = false;
    if (event->GetProcess()==99) flag_direct = true;
    if (event->GetProcess()==131 || event->GetProcess()==132) flag_direct = true;
    if (event->GetProcess()==135 || event->GetProcess()==136) flag_direct = true;

    if (!flag_direct) continue; // only process direct processes

    h1d_x_e[iQ2bin]->Fill( event->GetX() );
    if ( event_w_charm(event,gen_type) )
    {
      h1d_x_c[iQ2bin]->Fill( event->GetX() );
    }
    if ( event_w_strange(event,gen_type) )
    {
      h1d_x_s[iQ2bin]->Fill( event->GetX() );
    }
  }

  TFile* fout = new TFile(outFile,"recreate");
  fout->Write();
  h1d_nevt->Write();
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    h1d_x_e[iQ2]->Write();
    h1d_x_c[iQ2]->Write();
    h1d_x_s[iQ2]->Write();
  }
}
