R__LOAD_LIBRARY(libeicsmear);

#include "fast_sim.h"

TDatabasePDG* pdg = NULL;

const int chargebin = 3; // 0: -, 1: +, 2:+/-

const int etabin = 3;
static double eta_lo[etabin] = {-3,-1,1}; 
static double eta_hi[etabin] = {-1,1,3}; 

// parent pt bin
const int pptbin = 5;
static double ppt_lo[pptbin] = {0.2, 0.5, 1, 2, 4}; 
static double ppt_hi[pptbin] = {0.5, 1, 2, 4, 10}; 

const int Q2bin = 5;
static double Q2_lo[Q2bin] = {1E0, 5E0, 1E1, 5E1, 1E2}; 
static double Q2_hi[Q2bin] = {5E0, 1E1, 5E1, 1E2, 5E2}; 

const int xbin = 4;
static double x_lo[xbin] = {1E-4, 1E-3, 1E-2, 1E-1}; 
static double x_hi[xbin] = {1E-3, 1E-2, 1E-1, 1E0}; 

const double degree = 180./TMath::Pi();

const int verbosity = 1;

const double KMASS = 0.493677; // charged K mass, unit GeV
const double PIMASS = 0.139570; // charged pi mass, unit GeV
const double PMASS = 0.938272; // proton, unit GeV

using namespace std;

double dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{ // calculate the signed DCA value
  TVector3 posDiff = pos - vertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
  return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

int find_quark_origin(erhic::EventPythia* evt, erhic::ParticleMC* part)
{ // find the quark origin of a particle
  // is current particle quark or anti-quark?
  if (part->Id()>=1 && part->Id()<=8) return part->GetIndex(); // current particle is quark
  if (part->Id()<=-1 && part->Id()>=-8) return part->GetIndex(); // current particle is anti-quark

  // if not, does this curent particle has parent?
  int parent_indx = part->GetParentIndex();
  if (parent_indx==0) return 0;

  erhic::ParticleMC* parent = evt->GetTrack(parent_indx-1);
  if (verbosity>2)
  {
    if (pdg->GetParticle(parent->Id())) cout<< "parent id is " << pdg->GetParticle(parent->Id())->GetName() << endl;
    else cout<< "parent id is " << parent->Id() << endl;
  }
  
  return find_quark_origin(evt, parent);
}

class D0_reco
{
  private:
    // evt info
    TVector3 evt_vtx_true;
    TVector3 evt_vtx_reco;

    double x_true;
    double Q2_true;
    double nu_true;

    TLorentzVector hadron_beam;
    TLorentzVector struck_quark;

    // negative charged track true info
    vector<TLorentzVector> negl_p_true;
    vector<TVector3> negl_vtx_true;
    vector<int> negl_id_true;
    vector<int> negl_parent_index_true;
    vector<int> negl_parent_id_true;
    vector<TLorentzVector> negl_quark_p;

    // positive charged track true info
    vector<TLorentzVector> posl_p_true;
    vector<TVector3> posl_vtx_true;
    vector<int> posl_id_true;
    vector<int> posl_parent_index_true;
    vector<int> posl_parent_id_true;
    vector<TLorentzVector> posl_quark_p;

    // negative charged track reco info
    vector<TLorentzVector> negl_p_reco;
    vector<TVector3> negl_vtx_reco;
    vector<float> negl_prob_e;
    vector<float> negl_prob_pi;
    vector<float> negl_prob_K;
    vector<float> negl_prob_p;

    // positive charged track reco info
    vector<TLorentzVector> posl_p_reco;
    vector<TVector3> posl_vtx_reco;
    vector<float> posl_prob_e;
    vector<float> posl_prob_pi;
    vector<float> posl_prob_K;
    vector<float> posl_prob_p;

    // single track cut info
    float TRK_P_LO;
    float TRK_DCA;
    int ID_OPTION; // no PID, PID w/o assumption, PID w/ assumption
    int SMEAR_OPTION;
    int BFIELD_TYPE;

    // pair cut info
    float PAIR_DCA;
    float DECAY_L;
    float D0_DCA;
    float D0_COSTHETA;

    // mass pair vs pt
    TH2D* fg2d_Kpimass_vs_p[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpimass_vs_p[chargebin][etabin]; // 0: K-pi+

    // mass pair vs z
    TH2D* fg2d_Kpimass_vs_z[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpimass_vs_z[chargebin][etabin]; // 0: K-pi+

    TH2D* h2d_K_D0_p_vs_eta[etabin][pptbin];
    TH2D* h2d_pi_D0_p_vs_eta[etabin][pptbin];

    TH2D* h2d_D0_pt_vs_eta[Q2bin][xbin];
    TH2D* h2d_D0_z_vs_eta[Q2bin][xbin];

    TH2D* h2d_ztheo_vs_zjet;

  public:
    D0_reco()
    {
      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0; // 0--no smearing, 1--DM smearing, 1--LBL smearing, 2--Hybrid smearing
      BFIELD_TYPE = 0; // 0--Barbar, 1--Beast

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      D0_DCA = -9999;
      D0_COSTHETA = -9999;

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta] = new TH2D(Form("fg2d_Kpimass_vs_p_%d_%d",icharge,ieta),"Kpi mass vs Kpi p",480,0.6,3,100,0,10);
          fg2d_Kpimass_vs_p[icharge][ieta]->Sumw2();

          bg2d_Kpimass_vs_p[icharge][ieta] = new TH2D(Form("bg2d_Kpimass_vs_p_%d_%d",icharge,ieta),"Kpi mass vs Kpi p",480,0.6,3,100,0,10);
          bg2d_Kpimass_vs_p[icharge][ieta]->Sumw2();

          fg2d_Kpimass_vs_z[icharge][ieta] = new TH2D(Form("fg2d_Kpimass_vs_z_%d_%d",icharge,ieta),"Kpi mass vs Kpi z",480,0.6,3,100,0,1);
          fg2d_Kpimass_vs_z[icharge][ieta]->Sumw2();

          bg2d_Kpimass_vs_z[icharge][ieta] = new TH2D(Form("bg2d_Kpimass_vs_z_%d_%d",icharge,ieta),"Kpi mass vs Kpi z",480,0.6,3,100,0,1);
          bg2d_Kpimass_vs_z[icharge][ieta]->Sumw2();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_K_D0_p_vs_eta_%d_%d",ieta,ipt),"K p vs eta",100,0,10,160,-4,4);
          h2d_K_D0_p_vs_eta[ieta][ipt]->Sumw2();
          h2d_pi_D0_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_pi_D0_p_vs_eta_%d_%d",ieta,ipt),"pi p vs eta",100,0,10,160,-4,4);
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Sumw2();
        }
      }
      

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix] = new TH2D(Form("h2d_D0_pt_vs_eta_%d_%d",iQ2,ix),"D0 pt vs eta",100,0,10,160,-4,4);
          h2d_D0_pt_vs_eta[iQ2][ix]->Sumw2();

          h2d_D0_z_vs_eta[iQ2][ix] = new TH2D(Form("h2d_D0_z_vs_eta_%d_%d",iQ2,ix),"D0 z vs eta",100,0,1,160,-4,4);
          h2d_D0_z_vs_eta[iQ2][ix]->Sumw2();
        }
      }

      h2d_ztheo_vs_zjet = new TH2D("h2d_z_frag","z_{theo} vs z_{jet}",100,0,1,100,0,1);
      h2d_ztheo_vs_zjet->Sumw2();
    }
    virtual ~D0_reco()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          delete fg2d_Kpimass_vs_p[icharge][ieta];
          delete bg2d_Kpimass_vs_p[icharge][ieta];

          delete fg2d_Kpimass_vs_z[icharge][ieta];
          delete bg2d_Kpimass_vs_z[icharge][ieta];
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          delete h2d_K_D0_p_vs_eta[ieta][ipt];
          delete h2d_pi_D0_p_vs_eta[ieta][ipt];
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          delete h2d_D0_pt_vs_eta[iQ2][ix];
          delete h2d_D0_z_vs_eta[iQ2][ix];
        }
      }

      delete h2d_ztheo_vs_zjet;
    };

    void Reset()
    {
      x_true = 1E1;
      Q2_true = 1E-5;
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0;
      BFIELD_TYPE = 0;

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      D0_DCA = -9999;
      D0_COSTHETA = -9999;

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta]->Reset("ICESM");
          bg2d_Kpimass_vs_p[icharge][ieta]->Reset("ICESM");

          fg2d_Kpimass_vs_z[icharge][ieta]->Reset("ICESM");
          bg2d_Kpimass_vs_z[icharge][ieta]->Reset("ICESM");
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt]->Reset("ICESM");
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Reset("ICESM"); 
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix]->Reset("ICESM");
          h2d_D0_z_vs_eta[iQ2][ix]->Reset("ICESM");
        }
      }

      h2d_ztheo_vs_zjet->Reset("ICESM");
    }

    void SetLowP(const float mom_thr) { TRK_P_LO = mom_thr; }

    void SetDCACuts()
    {
      // single
      TRK_DCA = 0.03;

      // pair
      PAIR_DCA = 0.12; // 120um in unit of mm
      DECAY_L = 0.04;
      D0_DCA = 0.01;
      D0_COSTHETA = 0.98;
    }

    void SetIDCuts(const int id_opt) { ID_OPTION = id_opt; }

    void SetSmearType(const int smear_opt) { SMEAR_OPTION = smear_opt; }

    void SetBFieldType(const int b_opt) { BFIELD_TYPE = b_opt; }

    void SetVectTrue(TVector3& _evt_vtx_true) { evt_vtx_true = _evt_vtx_true; }
    void SetVectReco(TVector3& _evt_vtx_reco) { evt_vtx_reco = _evt_vtx_reco; }

    void SetQ2True(double _Q2_true) { Q2_true = _Q2_true; }
    void SetXTrue(double _x_true) { x_true = _x_true; }

    void SetNuTrue(double _nu_true) { nu_true = _nu_true; };

    void ClearTracks()
    {
      negl_p_true.clear();
      negl_vtx_true.clear();
      negl_id_true.clear();
      negl_parent_index_true.clear();
      negl_parent_id_true.clear();
      negl_quark_p.clear();
      
      negl_p_reco.clear();
      negl_vtx_reco.clear();
      negl_prob_e.clear();
      negl_prob_pi.clear();
      negl_prob_K.clear();
      negl_prob_p.clear();

      posl_p_true.clear();
      posl_vtx_true.clear();
      posl_id_true.clear();
      posl_parent_index_true.clear();
      posl_parent_id_true.clear();
      posl_quark_p.clear();

      posl_p_reco.clear();
      posl_vtx_reco.clear();
      posl_prob_e.clear();
      posl_prob_pi.clear();
      posl_prob_K.clear();
      posl_prob_p.clear();
    }

    void FillSingleTracks(erhic::EventPythia* py_evt)
    {
      ClearTracks();

      erhic::ParticleMC* proton = py_evt->GetTrack(1);
      if (proton!=NULL)
      {
        assert(abs(proton->Id())!=2212);
        hadron_beam = proton->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      erhic::ParticleMC* quark = py_evt->GetTrack(9);
      if (quark!=NULL && py_evt->GetProcess()==99)
      {
        struck_quark = quark->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (part->GetStatus()!=1) continue; // only loop through final stable particles

        if (abs(part->Id())!=2212 && abs(part->Id())!=211 && abs(part->Id())!=321 && abs(part->Id())!=11 && abs(part->Id())!=13) continue; // proton, pion, kaon, electron, muon

        if (TRK_P_LO>-99 && part->GetPt()<TRK_P_LO) continue; 

        TLorentzVector track_mom4_true = part->Get4Vector();
        TLorentzVector track_mom4_reco = track_mom4_true;

        TVector3 track_vtx_true = part->GetVertex();
        TVector3 track_vtx_reco = track_vtx_true;

        if (SMEAR_OPTION==0);
        if (SMEAR_OPTION==1)
        {
          track_mom4_reco = smearMomDM(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosDM(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==2)
        {
          track_mom4_reco = smearMomLBL(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosLBL(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==3)
        {
          track_mom4_reco = smearMomHybrid(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosHybrid(track_mom4_true.Vect(), track_vtx_true);
        }
        if (track_mom4_reco.E()>1000 || track_vtx_reco.Mag()>1000) continue; // outside eta or momentum range 

        // single track DCA cut
        double track_dca = dcaSigned(track_mom4_reco.Vect(),track_vtx_reco,evt_vtx_reco);
        if (TRK_DCA>-99 && fabs(track_dca)<TRK_DCA) continue;

        //==========================
        //      PID selection
        //==========================
        bitset<4> track_binary_id;         
        if (ID_OPTION==-1);
        else if (ID_OPTION==0) passing_DM_PID(track_mom4_reco,track_binary_id);
        else if (ID_OPTION==1 || ID_OPTION==2)
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        else
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE,0.5);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        if (verbosity>2) cout << "track_binary_id " <<track_binary_id.to_ulong() << endl;

        //==============================================================================
        //    Assumption: Cherekov detectors (mass ordering, no mu/pi separation)
        // assign probability according to the multiplicity (NB: qualitative for now)
        //==============================================================================
        float prob_pi = -999;
        float prob_K = -999;
        float prob_e = -999;
        float prob_p = -999;
        identify_charged_hadrons(part->Id(), track_binary_id, prob_e, prob_pi, prob_K, prob_p);

        int quark_index = find_quark_origin(py_evt,part);
        erhic::ParticleMC* quark_orig = py_evt->GetTrack(quark_index-1);

        //============================
        //  opposite charge selection
        //============================
        bool flag_pos = false, flag_neg = false;
        // no charge info in the input tree, use id to separate charge
        if (abs(part->Id())==11 && part->Id()<0) flag_pos = true;
        else if (abs(part->Id())==11 && part->Id()>0) flag_neg = true;
        else if (part->Id()>0) flag_pos = true;
        else flag_neg = true;

        if (flag_neg)
        {
          negl_p_true.push_back(track_mom4_true);
          negl_vtx_true.push_back(track_vtx_true);
          negl_id_true.push_back(part->Id());
          negl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL) negl_parent_id_true.push_back(parent->Id());
          else negl_parent_id_true.push_back(-9999);

          negl_p_reco.push_back(track_mom4_reco);
          negl_vtx_reco.push_back(track_vtx_reco);
          negl_prob_e.push_back(prob_e);
          negl_prob_pi.push_back(prob_pi);
          negl_prob_K.push_back(prob_K);
          negl_prob_p.push_back(prob_p);

          if (quark_orig!=NULL) negl_quark_p.push_back(quark_orig->Get4Vector());
          else negl_quark_p.push_back(TLorentzVector(0,0,0,0));
        }

        if (flag_pos)
        {
          posl_p_true.push_back(track_mom4_true);
          posl_vtx_true.push_back(track_vtx_true);
          posl_id_true.push_back(part->Id());
          posl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL) posl_parent_id_true.push_back(parent->Id());
          else posl_parent_id_true.push_back(-9999);

          posl_p_reco.push_back(track_mom4_reco);
          posl_vtx_reco.push_back(track_vtx_reco);
          posl_prob_e.push_back(prob_e);
          posl_prob_pi.push_back(prob_pi);
          posl_prob_K.push_back(prob_K);
          posl_prob_p.push_back(prob_p);

          if (quark_orig!=NULL) posl_quark_p.push_back(quark_orig->Get4Vector());
          else posl_quark_p.push_back(TLorentzVector(0,0,0,0));
        }
      }
    }

    void fill_Kpi_mass(const int charge_type, const int ipart1, const int ipart2)
    {
      if (charge_type>1) return;
      if (ipart1<0 || ipart2<0) return;

      bool is_SG = true;
      TLorentzVector kaon_p, pion_p;
      TLorentzVector quark_p;
      if (charge_type==0)
      { // K-pi+
        if (fabs(negl_parent_id_true[ipart1])!=421) is_SG = false;
        if (fabs(posl_parent_id_true[ipart2])!=421) is_SG = false;
        if (negl_parent_index_true[ipart1]!=posl_parent_index_true[ipart2]) is_SG = false;

        kaon_p = negl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(negl_id_true[ipart1])!=321) is_SG = false;

        pion_p = posl_p_reco[ipart2];
        pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
        if (fabs(posl_id_true[ipart2])!=211) is_SG = false;

        quark_p = negl_quark_p[ipart1];
      }
      else
      { // K+pi-
        if (fabs(posl_parent_id_true[ipart1])!=421) is_SG = false;
        if (fabs(negl_parent_id_true[ipart2])!=421) is_SG = false;
        if (posl_parent_index_true[ipart1]!=negl_parent_index_true[ipart2]) is_SG = false;

        kaon_p = posl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(posl_id_true[ipart1])!=321) is_SG = false;

        pion_p = negl_p_reco[ipart2];
        pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
        if (fabs(negl_id_true[ipart2])!=211) is_SG = false;

        quark_p = posl_quark_p[ipart1];
      }

      TLorentzVector pair(kaon_p+pion_p);
      int ietabin = -9999;
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (pair.PseudoRapidity()>=eta_lo[ieta] && pair.PseudoRapidity()<eta_hi[ieta]) ietabin = ieta;
      }
      if (ietabin<0) return;

      double frag_z = hadron_beam.Dot(pair)/(nu_true*hadron_beam.M());  // z= Pp/Pq where Pq=nuM

      fg2d_Kpimass_vs_p[charge_type][ietabin]->Fill(pair.M(),pair.Pt());
      fg2d_Kpimass_vs_p[2][ietabin]->Fill(pair.M(),pair.Pt());
      fg2d_Kpimass_vs_z[charge_type][ietabin]->Fill(pair.M(),frag_z);
      fg2d_Kpimass_vs_z[2][ietabin]->Fill(pair.M(),frag_z);
      if (!is_SG)
      {
        bg2d_Kpimass_vs_p[charge_type][ietabin]->Fill(pair.M(),pair.Pt());
        bg2d_Kpimass_vs_p[2][ietabin]->Fill(pair.M(),pair.Pt());
        bg2d_Kpimass_vs_z[charge_type][ietabin]->Fill(pair.M(),frag_z);
        bg2d_Kpimass_vs_z[2][ietabin]->Fill(pair.M(),frag_z);
      }
      else
      {
        // cout<<"D0 pt = "<<pair.Pt()<<" eta = "<<pair.PseudoRapidity()<<" z = "<<frag_z<<endl;
        // if (struck_quark.E()>0.001)
        // {
        //   cout<<"D0 energy = "<<pair.E()<<" struck quark energy = "<<struck_quark.E()<<" frac = "<<pair.E()/struck_quark.E()<<endl;
        //   // cout<<"D0 Pt = "<<pair.Pt()<<" struck quark Pt = "<<struck_quark.Pt()<<" frac = "<<
        //   cout<<"frac = "<<(pair.Vect()).Dot(struck_quark.Vect())/(struck_quark.Vect()).Dot(struck_quark.Vect())<<endl;
        //   cout<<"**********************"<<endl;
        // } 

        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          if (pair.Pt()>=ppt_lo[ipt] && pair.Pt()<ppt_hi[ipt])
          {
            h2d_K_D0_p_vs_eta[ietabin][ipt]->Fill(kaon_p.P(),kaon_p.PseudoRapidity());
            h2d_pi_D0_p_vs_eta[ietabin][ipt]->Fill(pion_p.P(),pion_p.PseudoRapidity());
          }
        }

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_D0_pt_vs_eta[iQ2bin][ixbin]->Fill(pair.Pt(),pair.PseudoRapidity());
          h2d_D0_z_vs_eta[iQ2bin][ixbin]->Fill(frag_z,pair.PseudoRapidity());
        } 

        h2d_ztheo_vs_zjet->Fill((pair.Vect()).Dot(quark_p.Vect())/(quark_p.Vect()).Dot(quark_p.Vect()),frag_z);
      }
    }

    void FillD0Pairs()
    {
      for (int ineg = 0; ineg < negl_id_true.size(); ++ineg)
      {
        for (int ipos = 0; ipos < posl_id_true.size(); ++ipos)
        {
          //==========================
          //    decay topology cut
          //==========================
          // pair DCA < cut value
          TVector3 dca_pair = negl_vtx_reco[ineg]-posl_vtx_reco[ipos];
          if (PAIR_DCA>-99 && dca_pair.Mag()>PAIR_DCA) continue;

          // Decay length > cut value
          TVector3 decay_l = (negl_vtx_reco[ineg]+posl_vtx_reco[ipos])*0.5-evt_vtx_reco;
          if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue;

          // D0 DCA > cut value
          TVector3 D0_vec = negl_p_reco[ineg].Vect()+posl_p_reco[ipos].Vect();
          double D0_dca = dcaSigned(D0_vec, decay_l, evt_vtx_reco);
          if (D0_DCA>-99 && fabs(D0_dca)<D0_DCA) continue;

          // D0 cos(theta) > cut value
          double D0_costheta = TMath::Cos(D0_vec.Angle(decay_l));
          if (D0_COSTHETA>-99 && D0_costheta<D0_costheta) continue;
          
          if (ID_OPTION==-1)
          { // no hID (but with eID)
            if (fabs(negl_id_true[ineg])!=11 && fabs(posl_id_true[ipos])!=11)
            {
              fill_Kpi_mass(0,ineg,ipos);
              fill_Kpi_mass(1,ipos,ineg);
            }
          } 

          if (ID_OPTION==0)
          { // PID with no low momentum cutoff
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          } 

          if (ID_OPTION==1)
          { // PID with low momentum cutoff & some mis-identified pi, K
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]>0.5)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]>0.5)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }   

          if (ID_OPTION==2)
          { // PID with low momentum cutoff & all identified pi, K
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }

          if (ID_OPTION==3)
          { // PID with low momentum cutoff & all identified pi, K (DIRC at 50cm)
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }          
        }
      }
    }

    void Write()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta]->Write();
          bg2d_Kpimass_vs_p[icharge][ieta]->Write();
          fg2d_Kpimass_vs_z[icharge][ieta]->Write();
          bg2d_Kpimass_vs_z[icharge][ieta]->Write();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt]->Write();
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Write();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix]->Write();
          h2d_D0_z_vs_eta[iQ2][ix]->Write();
        }
      }

      h2d_ztheo_vs_zjet->Write();
    }
};

class Lc_reco
{
  private:
    // evt info
    TVector3 evt_vtx_true;
    TVector3 evt_vtx_reco;

    double x_true;
    double Q2_true;
    double nu_true;

    TLorentzVector hadron_beam;

    // negative charged track true info
    vector<TLorentzVector> negl_p_true;
    vector<TVector3> negl_vtx_true;
    vector<int> negl_id_true;
    vector<int> negl_parent_index_true;
    vector<int> negl_parent_id_true;
    vector<int> negl_gparent_index_true;
    vector<int> negl_gparent_id_true;

    // positive charged track true info
    vector<TLorentzVector> posl_p_true;
    vector<TVector3> posl_vtx_true;
    vector<int> posl_id_true;
    vector<int> posl_parent_index_true;
    vector<int> posl_parent_id_true;
    vector<int> posl_gparent_index_true;
    vector<int> posl_gparent_id_true;

    // negative charged track reco info
    vector<TLorentzVector> negl_p_reco;
    vector<TVector3> negl_vtx_reco;
    vector<float> negl_prob_e;
    vector<float> negl_prob_pi;
    vector<float> negl_prob_K;
    vector<float> negl_prob_p;

    // positive charged track reco info
    vector<TLorentzVector> posl_p_reco;
    vector<TVector3> posl_vtx_reco;
    vector<float> posl_prob_e;
    vector<float> posl_prob_pi;
    vector<float> posl_prob_K;
    vector<float> posl_prob_p;

    // single track cut info
    float TRK_P_LO;
    float TRK_DCA;
    int ID_OPTION; // no PID, PID w/o assumption, PID w/ assumption
    int SMEAR_OPTION;
    int BFIELD_TYPE;

    // pair cut info
    float PAIR_DCA;
    float DECAY_L;
    float Lc_DCA;
    float Lc_COSTHETA;

    // mass pair vs pt
    TH2D* fg2d_Kpipmass_vs_p[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpipmass_vs_p[chargebin][etabin]; // 0: K-pi+

    // mass pair vs z
    TH2D* fg2d_Kpipmass_vs_z[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpipmass_vs_z[chargebin][etabin]; // 0: K-pi+

    TH2D* h2d_K_Lc_p_vs_eta[etabin][pptbin];
    TH2D* h2d_pi_Lc_p_vs_eta[etabin][pptbin];
    TH2D* h2d_p_Lc_p_vs_eta[etabin][pptbin];

    TH2D* h2d_Lc_pt_vs_eta[Q2bin][xbin];
    TH2D* h2d_Lc_z_vs_eta[Q2bin][xbin];

  public:
    Lc_reco()
    {
      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0;  // 0--no smearing, 1--DM smearing, 1--LBL smearing, 2--Hybrid smearing
      BFIELD_TYPE = 0; // 0--Barbar, 1--Beast

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      Lc_DCA = -9999;
      Lc_COSTHETA = -9999;

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpipmass_vs_p[icharge][ieta] = new TH2D(Form("fg2d_Kpipmass_vs_p_%d_%d",icharge,ieta),"Kpip mass vs Kpip p",600,1.5,3,100,0,10);
          fg2d_Kpipmass_vs_p[icharge][ieta]->Sumw2();

          bg2d_Kpipmass_vs_p[icharge][ieta] = new TH2D(Form("bg2d_Kpipmass_vs_p_%d_%d",icharge,ieta),"Kpip mass vs Kpip p",600,1.5,3,100,0,10);
          bg2d_Kpipmass_vs_p[icharge][ieta]->Sumw2();

          fg2d_Kpipmass_vs_z[icharge][ieta] = new TH2D(Form("fg2d_Kpipmass_vs_z_%d_%d",icharge,ieta),"Kpip mass vs Kpip z",600,1.5,3,100,0,1);
          fg2d_Kpipmass_vs_z[icharge][ieta]->Sumw2();

          bg2d_Kpipmass_vs_z[icharge][ieta] = new TH2D(Form("bg2d_Kpipmass_vs_z_%d_%d",icharge,ieta),"Kpip mass vs Kpip z",600,1.5,3,100,0,1);
          bg2d_Kpipmass_vs_z[icharge][ieta]->Sumw2();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_Lc_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_K_Lc_p_vs_eta_%d_%d",ieta,ipt),"K p vs eta",100,0,10,160,-4,4);
          h2d_K_Lc_p_vs_eta[ieta][ipt]->Sumw2();
          h2d_pi_Lc_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_pi_Lc_p_vs_eta_%d_%d",ieta,ipt),"pi p vs eta",100,0,10,160,-4,4);
          h2d_pi_Lc_p_vs_eta[ieta][ipt]->Sumw2();
          h2d_p_Lc_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_p_Lc_p_vs_eta_%d_%d",ieta,ipt),"p p vs eta",100,0,10,160,-4,4);
          h2d_p_Lc_p_vs_eta[ieta][ipt]->Sumw2();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_Lc_pt_vs_eta[iQ2][ix] = new TH2D(Form("h2d_Lc_pt_vs_eta_%d_%d",iQ2,ix),"Lc pt vs eta",100,0,10,160,-4,4);
          h2d_Lc_pt_vs_eta[iQ2][ix]->Sumw2();

          h2d_Lc_z_vs_eta[iQ2][ix] = new TH2D(Form("h2d_Lc_z_vs_eta_%d_%d",iQ2,ix),"Lc pt vs eta",100,0,1,160,-4,4);
          h2d_Lc_z_vs_eta[iQ2][ix]->Sumw2();
        }
      }
    }
    virtual ~Lc_reco()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          delete fg2d_Kpipmass_vs_p[icharge][ieta];
          delete bg2d_Kpipmass_vs_p[icharge][ieta];

          delete fg2d_Kpipmass_vs_z[icharge][ieta];
          delete bg2d_Kpipmass_vs_z[icharge][ieta];
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          delete h2d_K_Lc_p_vs_eta[ieta][ipt];
          delete h2d_pi_Lc_p_vs_eta[ieta][ipt];
          delete h2d_p_Lc_p_vs_eta[ieta][ipt];
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          delete h2d_Lc_pt_vs_eta[iQ2][ix];
          delete h2d_Lc_z_vs_eta[iQ2][ix];
        }
      }
    };

    void Reset()
    {
      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0;
      BFIELD_TYPE = 0;

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      Lc_DCA = -9999;
      Lc_COSTHETA = -9999;

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpipmass_vs_p[icharge][ieta]->Reset("ICESM");
          bg2d_Kpipmass_vs_p[icharge][ieta]->Reset("ICESM");
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_Lc_p_vs_eta[ieta][ipt]->Reset("ICESM");
          h2d_pi_Lc_p_vs_eta[ieta][ipt]->Reset("ICESM");
          h2d_p_Lc_p_vs_eta[ieta][ipt]->Reset("ICESM");
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_Lc_pt_vs_eta[iQ2][ix]->Reset("ICESM");
          h2d_Lc_z_vs_eta[iQ2][ix]->Reset("ICESM");
        }
      }
    }

    void SetLowP(const float mom_thr) { TRK_P_LO = mom_thr; }

    void SetDCACuts()
    {
      // single
      TRK_DCA = -9999;

      // pair
      PAIR_DCA = 0.3; // 300um in unit of mm
      DECAY_L = 0.01;
      Lc_DCA = 0.15;
      Lc_COSTHETA = -9999;
    }

    void SetIDCuts(const int id_opt) { ID_OPTION = id_opt; }

    void SetSmearType(const int smear_opt) { SMEAR_OPTION = smear_opt; }

    void SetBFieldType(const int b_opt) { BFIELD_TYPE = b_opt; }

    void SetVectTrue(TVector3& _evt_vtx_true) { evt_vtx_true = _evt_vtx_true; }
    void SetVectReco(TVector3& _evt_vtx_reco) { evt_vtx_reco = _evt_vtx_reco; }

    void SetQ2True(double _Q2_true) { Q2_true = _Q2_true; }
    void SetXTrue(double _x_true) { x_true = _x_true; }

    void SetNuTrue(double _nu_true) { nu_true = _nu_true; };

    void ClearTracks()
    {
      negl_p_true.clear();
      negl_vtx_true.clear();
      negl_id_true.clear();
      negl_parent_index_true.clear();
      negl_parent_id_true.clear();
      negl_gparent_index_true.clear();
      negl_gparent_id_true.clear();
      
      negl_p_reco.clear();
      negl_vtx_reco.clear();
      negl_prob_e.clear();
      negl_prob_pi.clear();
      negl_prob_K.clear();
      negl_prob_p.clear();

      posl_p_true.clear();
      posl_vtx_true.clear();
      posl_id_true.clear();
      posl_parent_index_true.clear();
      posl_parent_id_true.clear();
      posl_gparent_index_true.clear();
      posl_gparent_id_true.clear();

      posl_p_reco.clear();
      posl_vtx_reco.clear();
      posl_prob_e.clear();
      posl_prob_pi.clear();
      posl_prob_K.clear();
      posl_prob_p.clear();
    }

    void FillSingleTracks(erhic::EventPythia* py_evt)
    {
      ClearTracks();

      erhic::ParticleMC* proton = py_evt->GetTrack(1);
      if (proton!=NULL)
      {
        assert(abs(proton->Id())!=2212);
        hadron_beam = proton->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (part->GetStatus()!=1) continue; // only loop through final stable particles

        if (abs(part->Id())!=2212 && abs(part->Id())!=211 && abs(part->Id())!=321 && abs(part->Id())!=11 && abs(part->Id())!=13) continue; // proton, pion, kaon, electron, muon

        if (TRK_P_LO>-99 && part->GetPt()<TRK_P_LO) continue; 

        TLorentzVector track_mom4_true = part->Get4Vector();
        TLorentzVector track_mom4_reco = track_mom4_true;

        TVector3 track_vtx_true = part->GetVertex();
        TVector3 track_vtx_reco = track_vtx_true;

        if (SMEAR_OPTION==0);
        if (SMEAR_OPTION==1)
        {
          track_mom4_reco = smearMomDM(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosDM(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==2)
        {
          track_mom4_reco = smearMomLBL(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosLBL(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==3)
        {
          track_mom4_reco = smearMomHybrid(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosHybrid(track_mom4_true.Vect(), track_vtx_true);
        }
        if (track_mom4_reco.E()>1000 || track_vtx_reco.Mag()>1000) continue; // outside eta or momentum range 

        // single track DCA cut
        double track_dca = dcaSigned(track_mom4_reco.Vect(),track_vtx_reco,evt_vtx_reco);
        if (TRK_DCA>-99 && fabs(track_dca)<TRK_DCA) continue;

        //==========================
        //      PID selection
        //==========================
        bitset<4> track_binary_id;         
        if (ID_OPTION==-1);
        else if (ID_OPTION==0) passing_DM_PID(track_mom4_reco,track_binary_id);
        else if (ID_OPTION==1 || ID_OPTION==2)
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        else
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE,0.5);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        if (verbosity>2) cout << "track_binary_id " <<track_binary_id.to_ulong() << endl;

        //==============================================================================
        //    Assumption: Cherekov detectors (mass ordering, no mu/pi separation)
        // assign probability according to the multiplicity (NB: qualitative for now)
        //==============================================================================
        float prob_pi = -999;
        float prob_K = -999;
        float prob_e = -999;
        float prob_p = -999;
        identify_charged_hadrons(part->Id(), track_binary_id, prob_e, prob_pi, prob_K, prob_p);

        //============================
        //  opposite charge selection
        //============================
        bool flag_pos = false, flag_neg = false;
        // no charge info in the input tree, use id to separate charge
        if (abs(part->Id())==11 && part->Id()<0) flag_pos = true;
        else if (abs(part->Id())==11 && part->Id()>0) flag_neg = true;
        else if (part->Id()>0) flag_pos = true;
        else flag_neg = true;

        if (flag_neg)
        {
          negl_p_true.push_back(track_mom4_true);
          negl_vtx_true.push_back(track_vtx_true);
          negl_id_true.push_back(part->Id());
          negl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL)
          {
            negl_parent_id_true.push_back(parent->Id());
            negl_gparent_index_true.push_back(parent->GetParentIndex());

            erhic::ParticleMC* gparent = py_evt->GetTrack(parent->GetParentIndex()-1);
            if (gparent!=NULL) negl_gparent_id_true.push_back(gparent->Id());
            else negl_gparent_id_true.push_back(-9999);
          } 
          else
          {
            negl_parent_id_true.push_back(-9999);
            negl_gparent_index_true.push_back(-9999);
            negl_gparent_id_true.push_back(-9999);
          } 

          negl_p_reco.push_back(track_mom4_reco);
          negl_vtx_reco.push_back(track_vtx_reco);
          negl_prob_e.push_back(prob_e);
          negl_prob_pi.push_back(prob_pi);
          negl_prob_K.push_back(prob_K);
          negl_prob_p.push_back(prob_p);
        }

        if (flag_pos)
        {
          posl_p_true.push_back(track_mom4_true);
          posl_vtx_true.push_back(track_vtx_true);
          posl_id_true.push_back(part->Id());
          posl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL)
          {
            posl_parent_id_true.push_back(parent->Id());
            posl_gparent_index_true.push_back(parent->GetParentIndex());

            erhic::ParticleMC* gparent = py_evt->GetTrack(parent->GetParentIndex()-1);
            if (gparent!=NULL) posl_gparent_id_true.push_back(gparent->Id());
            else posl_gparent_id_true.push_back(-9999);
          } 
          else
          {
            posl_parent_id_true.push_back(-9999);
            posl_gparent_index_true.push_back(-9999);
            posl_gparent_id_true.push_back(-9999);
          } 

          posl_p_reco.push_back(track_mom4_reco);
          posl_vtx_reco.push_back(track_vtx_reco);
          posl_prob_e.push_back(prob_e);
          posl_prob_pi.push_back(prob_pi);
          posl_prob_K.push_back(prob_K);
          posl_prob_p.push_back(prob_p);
        }
      }
    }

    void fill_Kpip_mass(const int charge_type, const int mass_order, const int ipart1, const int ipart2, const int ipart3)
    {
      if (charge_type>1 || mass_order>1) return;
      if (ipart1<0 || ipart2<0 || ipart3<0) return;

      bool is_SG = true;
      TLorentzVector kaon_p, pion_p, proton_p;
      if (charge_type==0)
      { // K-
        int Lc_index_part1 = -9999;
        if (fabs(negl_parent_id_true[ipart1])==4122) Lc_index_part1 = negl_parent_index_true[ipart1];
        if (fabs(negl_gparent_id_true[ipart1])==4122) Lc_index_part1 = negl_gparent_index_true[ipart1];

        int Lc_index_part2 = -9999;
        if (fabs(posl_parent_id_true[ipart2])==4122) Lc_index_part2 = posl_parent_index_true[ipart2];
        if (fabs(posl_gparent_id_true[ipart2])==4122) Lc_index_part2 = posl_gparent_index_true[ipart2];

        int Lc_index_part3 = -9999;
        if (fabs(posl_parent_id_true[ipart3])==4122) Lc_index_part3 = posl_parent_index_true[ipart3];
        if (fabs(posl_gparent_id_true[ipart3])==4122) Lc_index_part3 = posl_gparent_index_true[ipart3];

        if (Lc_index_part1<0 || Lc_index_part2<0 || Lc_index_part3<0) is_SG = false;
        if (Lc_index_part1!=Lc_index_part2 || Lc_index_part2!=Lc_index_part3) is_SG = false;

        kaon_p = negl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(negl_id_true[ipart1])!=321) is_SG = false;
        
        if (mass_order==0)
        { // pi+p+
          pion_p = posl_p_reco[ipart2];
          pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
          if (fabs(posl_id_true[ipart2])!=211) is_SG = false;

          proton_p = posl_p_reco[ipart3];
          proton_p.SetXYZM(proton_p.X(),proton_p.Y(),proton_p.Z(),PMASS);
          if (fabs(posl_id_true[ipart3])!=2212) is_SG = false;
        }
        else
        { // p+pi+
          pion_p = posl_p_reco[ipart3];
          pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
          if (fabs(posl_id_true[ipart3])!=211) is_SG = false;

          proton_p = posl_p_reco[ipart2];
          proton_p.SetXYZM(proton_p.X(),proton_p.Y(),proton_p.Z(),PMASS);
          if (fabs(posl_id_true[ipart2])!=2212) is_SG = false;
        }
      }
      else
      { // K+
        int Lc_index_part1 = -9999;
        if (fabs(posl_parent_id_true[ipart1])==4122) Lc_index_part1 = posl_parent_index_true[ipart1];
        if (fabs(posl_gparent_id_true[ipart1])==4122) Lc_index_part1 = posl_gparent_index_true[ipart1];

        int Lc_index_part2 = -9999;
        if (fabs(negl_parent_id_true[ipart2])==4122) Lc_index_part2 = negl_parent_index_true[ipart2];
        if (fabs(negl_gparent_id_true[ipart2])==4122) Lc_index_part2 = negl_gparent_index_true[ipart2];

        int Lc_index_part3 = -9999;
        if (fabs(negl_parent_id_true[ipart3])==4122) Lc_index_part3 = negl_parent_index_true[ipart3];
        if (fabs(negl_gparent_id_true[ipart3])==4122) Lc_index_part3 = negl_gparent_index_true[ipart3];

        if (Lc_index_part1<0 || Lc_index_part2<0 || Lc_index_part3<0) return;
        if (Lc_index_part1!=Lc_index_part2 || Lc_index_part2!=Lc_index_part3) return;

        kaon_p = posl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(posl_id_true[ipart1])!=321) is_SG = false;
        
        if (mass_order==0)
        { // pi-p-
          pion_p = negl_p_reco[ipart2];
          pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
          if (fabs(negl_id_true[ipart2])!=211) is_SG = false;

          proton_p = negl_p_reco[ipart3];
          proton_p.SetXYZM(proton_p.X(),proton_p.Y(),proton_p.Z(),PMASS);
          if (fabs(negl_id_true[ipart3])!=2212) is_SG = false;
        }
        else
        { // p-pi-
          pion_p = negl_p_reco[ipart3];
          pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
          if (fabs(negl_id_true[ipart3])!=211) is_SG = false;

          proton_p = negl_p_reco[ipart2];
          proton_p.SetXYZM(proton_p.X(),proton_p.Y(),proton_p.Z(),PMASS);
          if (fabs(negl_id_true[ipart2])!=2212) is_SG = false;
        }
      }

      TLorentzVector trip(kaon_p+pion_p+proton_p);
      // if (is_SG) cout<<"mass "<<trip.M()<<endl;
      int ietabin = -9999;
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (trip.PseudoRapidity()>=eta_lo[ieta] && trip.PseudoRapidity()<eta_hi[ieta]) ietabin = ieta;
      }
      if (ietabin<0) return;

      double frag_z = hadron_beam.Dot(trip)/(nu_true*hadron_beam.M());  // z= Pp/Pq where Pq=nuM
        
      fg2d_Kpipmass_vs_p[charge_type][ietabin]->Fill(trip.M(),trip.Pt());
      fg2d_Kpipmass_vs_p[2][ietabin]->Fill(trip.M(),trip.Pt());
      fg2d_Kpipmass_vs_z[charge_type][ietabin]->Fill(trip.M(),frag_z);
      fg2d_Kpipmass_vs_z[2][ietabin]->Fill(trip.M(),frag_z);
      if (!is_SG)
      {
        bg2d_Kpipmass_vs_p[charge_type][ietabin]->Fill(trip.M(),trip.Pt());
        bg2d_Kpipmass_vs_p[2][ietabin]->Fill(trip.M(),trip.Pt());
        bg2d_Kpipmass_vs_z[charge_type][ietabin]->Fill(trip.M(),frag_z);
        bg2d_Kpipmass_vs_z[2][ietabin]->Fill(trip.M(),frag_z);
      }
      else
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          if (trip.Pt()>=ppt_lo[ipt] && trip.Pt()<ppt_hi[ipt])
          {
            cout<<"Lc pt "<<trip.Pt()<<" eta "<<trip.PseudoRapidity()<<endl;
            h2d_K_Lc_p_vs_eta[ietabin][ipt]->Fill(kaon_p.P(),kaon_p.PseudoRapidity());
            h2d_pi_Lc_p_vs_eta[ietabin][ipt]->Fill(pion_p.P(),pion_p.PseudoRapidity());
            h2d_p_Lc_p_vs_eta[ietabin][ipt]->Fill(proton_p.P(),proton_p.PseudoRapidity());
          }
        }

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_Lc_pt_vs_eta[iQ2bin][ixbin]->Fill(trip.Pt(),trip.PseudoRapidity());
          h2d_Lc_z_vs_eta[iQ2bin][ixbin]->Fill(frag_z,trip.PseudoRapidity());
        } 
      }
    }

    void FillLcTriplets()
    {
      for (int ineg = 0; ineg < negl_id_true.size(); ++ineg)
      {
        for (int ipos1 = 0; ipos1 < posl_id_true.size(); ++ipos1)
        {
          for (int ipos2 = ipos1+1; ipos2 < posl_id_true.size(); ++ipos2)
          {
            //==========================
            //    decay topology cut
            //==========================
            // pair DCA < cut value
            TVector3 dca_pair1 = negl_vtx_reco[ineg]-posl_vtx_reco[ipos1];
            TVector3 dca_pair2 = negl_vtx_reco[ineg]-posl_vtx_reco[ipos2];
            TVector3 dca_pair3 = posl_vtx_reco[ipos1]-posl_vtx_reco[ipos2];
            double dca_pair[3] = {dca_pair1.Mag(),dca_pair2.Mag(),dca_pair3.Mag()};
            if (PAIR_DCA>-99 && TMath::MaxElement(3,dca_pair)>PAIR_DCA) continue;

            // Decay length > cut value
            TVector3 decay_l = (negl_vtx_reco[ineg]+posl_vtx_reco[ipos1]+posl_vtx_reco[ipos2])*(1./3)-evt_vtx_reco;
            if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue;

            // Lc DCA > cut value
            TVector3 Lc_vec = negl_p_reco[ineg].Vect()+posl_p_reco[ipos1].Vect()+posl_p_reco[ipos2].Vect();
            double Lc_dca = dcaSigned(Lc_vec, decay_l, evt_vtx_reco);
            if (Lc_DCA>-99 && fabs(Lc_dca)>Lc_DCA) continue;

            // Lc cos(theta) > cut value
            double Lc_costheta = TMath::Cos(Lc_vec.Angle(decay_l));
            if (Lc_COSTHETA>-99 && Lc_costheta<Lc_costheta) continue;

            if (ID_OPTION==-1)
            { // no hID (but with eID)
              if (fabs(negl_id_true[ineg])!=11 && fabs(posl_id_true[ipos1])!=11 && fabs(posl_id_true[ipos2])!=11)
              {
                // K-pi+p+
                fill_Kpip_mass(0,0,ineg,ipos1,ipos2);

                // K-p+pi+
                fill_Kpip_mass(0,1,ineg,ipos1,ipos2);
              }
            } 

            if (ID_OPTION==0)
            { // PID with no low momentum cutoff
              if (negl_prob_K[ineg]==1 && posl_prob_pi[ipos1]==1 && posl_prob_p[ipos2]==1)
              { // K-pi+p+
                fill_Kpip_mass(0,0,ineg,ipos1,ipos2);
                // if (fabs(negl_parent_id_true[ineg])==4122) cout<<"Lambda c"<<endl;
              }
              if (negl_prob_K[ineg]==1 && posl_prob_p[ipos1]==1 && posl_prob_pi[ipos2]==1)
              { // K-p+pi+
                fill_Kpip_mass(0,1,ineg,ipos1,ipos2);
              }
            } 

            if (ID_OPTION==1)
            { // PID with low momentum cutoff & some mis-identified K, p
              if (negl_prob_K[ineg]>0.5 && posl_prob_pi[ipos1]==1 && posl_prob_p[ipos2]>0.1)
              { // K-pi+p+
                fill_Kpip_mass(0,0,ineg,ipos1,ipos2);
              }
              if (negl_prob_K[ineg]>0.5 && posl_prob_p[ipos1]>0.1 && posl_prob_pi[ipos2]==1)
              { // K-p+pi+
                fill_Kpip_mass(0,1,ineg,ipos1,ipos2);
              }
            }   

            if (ID_OPTION==2)
            { // PID with low momentum cutoff & all identified pi, K, p
              if (negl_prob_K[ineg]==1 && posl_prob_pi[ipos1]==1 && posl_prob_p[ipos2]==1)
              { // K-pi+p+
                fill_Kpip_mass(0,0,ineg,ipos1,ipos2);
              }
              if (negl_prob_K[ineg]==1 && posl_prob_p[ipos1]==1 && posl_prob_pi[ipos2]==1)
              { // K-p+pi+
                fill_Kpip_mass(0,1,ineg,ipos1,ipos2);
              }
            }

            if (ID_OPTION==3)
            { // PID with low momentum cutoff & all identified pi, K, p (DIRC @ 50cm)
              if (negl_prob_K[ineg]==1 && posl_prob_pi[ipos1]==1 && posl_prob_p[ipos2]==1)
              { // K-pi+p+
                fill_Kpip_mass(0,0,ineg,ipos1,ipos2);
              }
              if (negl_prob_K[ineg]==1 && posl_prob_p[ipos1]==1 && posl_prob_pi[ipos2]==1)
              { // K-p+pi+
                fill_Kpip_mass(0,1,ineg,ipos1,ipos2);
              }
            }   
          }
        }
      }

      for (int ipos = 0; ipos < posl_id_true.size(); ++ipos)
      {
        for (int ineg1 = 0; ineg1 < negl_id_true.size(); ++ineg1)
        {
          for (int ineg2 = ineg1+1; ineg2 < negl_id_true.size(); ++ineg2)
          {
            //==========================
            //    decay topology cut
            //==========================
            // pair DCA < cut value
            TVector3 dca_pair1 = posl_vtx_reco[ipos]-negl_vtx_reco[ineg1];
            TVector3 dca_pair2 = posl_vtx_reco[ipos]-negl_vtx_reco[ineg2];
            TVector3 dca_pair3 = negl_vtx_reco[ineg1]-negl_vtx_reco[ineg2];
            double dca_pair[3] = {dca_pair1.Mag(),dca_pair2.Mag(),dca_pair3.Mag()};
            if (PAIR_DCA>-99 && TMath::MaxElement(3,dca_pair)>PAIR_DCA) continue;

            // Decay length > cut value
            TVector3 decay_l = (posl_vtx_reco[ipos]+negl_vtx_reco[ineg1]+negl_vtx_reco[ineg2])*(1./3)-evt_vtx_reco;
            if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue;

            // Lc DCA > cut value
            TVector3 Lc_vec = posl_p_reco[ipos].Vect()+negl_p_reco[ineg1].Vect()+negl_p_reco[ineg2].Vect();
            double Lc_dca = dcaSigned(Lc_vec, decay_l, evt_vtx_reco);
            if (Lc_DCA>-99 && fabs(Lc_dca)>Lc_DCA) continue;

            // Lc cos(theta) > cut value
            double Lc_costheta = TMath::Cos(Lc_vec.Angle(decay_l));
            if (Lc_COSTHETA>-99 && Lc_costheta<Lc_costheta) continue;

            if (ID_OPTION==-1)
            { // no hID (but with eID)
              if (fabs(posl_id_true[ipos])!=11 && fabs(negl_id_true[ineg1])!=11 && fabs(negl_id_true[ineg2])!=11)
              {
                // K+pi-p-
                fill_Kpip_mass(1,0,ipos,ineg1,ineg2);

                // K+p-pi-
                fill_Kpip_mass(1,1,ipos,ineg1,ineg2);
              }
            } 

            if (ID_OPTION==0)
            { // PID with no low momentum cutoff
              if (posl_prob_K[ipos]==1 && negl_prob_pi[ineg1]==1 && negl_prob_p[ineg2]==1)
              { // K+pi-p-
                fill_Kpip_mass(1,0,ipos,ineg1,ineg2);
              }
              if (posl_prob_K[ipos]==1 && negl_prob_p[ineg1]==1 && negl_prob_pi[ineg2]==1)
              { // K+p-pi-
                fill_Kpip_mass(1,1,ipos,ineg1,ineg2);
              }
            } 

            if (ID_OPTION==1)
            { // PID with low momentum cutoff & some mis-identified K, p
              if (posl_prob_K[ipos]>0.5 && negl_prob_pi[ineg1]==1 && negl_prob_p[ineg2]>0.1)
              { // K+pi-p-
                fill_Kpip_mass(1,0,ipos,ineg1,ineg2);
              }
              if (posl_prob_K[ipos]>0.5 && negl_prob_p[ineg1]>0.1 && negl_prob_pi[ineg2]==1)
              { // K+p-pi-
                fill_Kpip_mass(1,1,ipos,ineg1,ineg2);
              }
            }   

            if (ID_OPTION==2)
            { // PID with low momentum cutoff & all identified pi, K, p
              if (posl_prob_K[ipos]==1 && negl_prob_pi[ineg1]==1 && negl_prob_p[ineg2]==1)
              { // K-pi+p+
                fill_Kpip_mass(1,0,ipos,ineg1,ineg2);
              }
              if (posl_prob_K[ipos]==1 && negl_prob_p[ineg1]==1 && negl_prob_pi[ineg2]==1)
              { // K-p+pi+
                fill_Kpip_mass(1,1,ipos,ineg1,ineg2);
              }
            }

            if (ID_OPTION==3)
            { // PID with low momentum cutoff & all identified pi, K, p (DIRC @ 50cm)
              if (posl_prob_K[ipos]==1 && negl_prob_pi[ineg1]==1 && negl_prob_p[ineg2]==1)
              { // K-pi+p+
                fill_Kpip_mass(1,0,ipos,ineg1,ineg2);
              }
              if (posl_prob_K[ipos]==1 && negl_prob_p[ineg1]==1 && negl_prob_pi[ineg2]==1)
              { // K-p+pi+
                fill_Kpip_mass(1,1,ipos,ineg1,ineg2);
              }
            } 
          }
        }
      }
    }

    void Write()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpipmass_vs_p[icharge][ieta]->Write();
          bg2d_Kpipmass_vs_p[icharge][ieta]->Write();

          fg2d_Kpipmass_vs_z[icharge][ieta]->Write();
          bg2d_Kpipmass_vs_z[icharge][ieta]->Write();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_Lc_p_vs_eta[ieta][ipt]->Write();
          h2d_pi_Lc_p_vs_eta[ieta][ipt]->Write();
          h2d_p_Lc_p_vs_eta[ieta][ipt]->Write();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_Lc_pt_vs_eta[iQ2][ix]->Write();
          h2d_Lc_z_vs_eta[iQ2][ix]->Write();
        }
      }
    }
};

void D0_tree(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hist.root", int nevt = 0, const int smear_option = 0, const int Bfield_type = 0, const int PID_option = 0)
{ // smear_2nd_vtx & momentum: 0--no smearing, 1--DM smearing, 2--LBL smearing, 3--Hybrid smearing
  // Bfield_type: 0--Barbar, 1--Beast
  // 0--no hID (but with eID), 1--PID with no low momentum cutoff, 2--PID with low momentum cutoff & some mis-identified pi, K, 3--PID with low momentum cutoff & all identified pi, K

  // PDG data table
  pdg = new TDatabasePDG();

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  // erhic::ParticleMC *child(NULL);
  
  //Load ROOT File
  //TFile *f = new TFile("../pythia/outfiles/ep_10_100_norad_def.root"); //Not created by eic or jlab version
  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  //Define Some Variables
  Float_t Q2(0);
  Int_t nParticles(0);

  D0_reco ana_D0;
  ana_D0.SetLowP(0.1);
  ana_D0.SetDCACuts();
  ana_D0.SetIDCuts(PID_option);
  ana_D0.SetSmearType(smear_option);
  ana_D0.SetBFieldType(Bfield_type);

  Lc_reco ana_Lc;
  ana_Lc.SetDCACuts();
  ana_Lc.SetIDCuts(PID_option);
  ana_Lc.SetSmearType(smear_option);
  ana_Lc.SetBFieldType(Bfield_type);
  
  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {    
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    //Write Out Q2
    Q2 = (Float_t) event->GetQ2(); //Can also do event->QSquared
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);
    // if (Q2>10) continue; // process low Q2
    
    //Get Total Number of Particles
    nParticles = event->GetNTracks();
    // printf("For Event %d, we have %d particles!\n",ievt,nParticles);

    TVector3 evt_vtx(0,0,0);
    ana_D0.SetVectTrue(evt_vtx);
    ana_D0.SetVectReco(evt_vtx); // NB: no smear on primary vertex yet

    ana_D0.SetQ2True(event->GetQ2());
    ana_D0.SetXTrue(event->GetX());

    ana_D0.SetNuTrue(event->GetNu());

    ana_D0.FillSingleTracks(event);
    ana_D0.FillD0Pairs();

    ana_Lc.SetVectTrue(evt_vtx);
    ana_Lc.SetVectReco(evt_vtx); // NB: no smear on primary vertex yet

    ana_Lc.SetQ2True(event->GetQ2());
    ana_Lc.SetXTrue(event->GetX());

    ana_Lc.SetNuTrue(event->GetNu());

    ana_Lc.FillSingleTracks(event);
    ana_Lc.FillLcTriplets();
  }

  TFile* fout = new TFile(outFile,"recreate");

  ana_D0.Write();
  ana_Lc.Write();

  fout->Write();
}
