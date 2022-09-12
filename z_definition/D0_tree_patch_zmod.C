// changes made:
// remove the condition before applying ATHENA smearing
// Lc DCA cut change to < cutting value (better to confirm with Yuanjing)
// remove struck quark requirement in D0 FillSingleTracks
// change all cut to transverse cuts
// remove single track DCA cut and D0 DCA cut

// modded D0_tree_patch to support h2d_ztheo_vs_zjet bins and cuts on event types
// of the [processid] bin of this histogram,
// first bin is process=99, second bin is process=135 or 136, third bin is inclusive
// only for D0reco
// index of the bin is encodded in PROCESS_INDEX

R__LOAD_LIBRARY(libeicsmear);

#include "fast_sim.h"
#include "bins_fine.h"

TDatabasePDG* pdg = NULL;

const double degree = 180./TMath::Pi();

const int verbosity = 1;

const double KMASS = 0.493677; // charged K mass, unit GeV
const double PIMASS = 0.139570; // charged pi mass, unit GeV
const double PMASS = 0.938272; // proton, unit GeV
const double LIGHT_SPEED = 299792458; // unit m/s
const double D0_MEAN_LIFE = 410.1E-15; // \pm 1.5 10^{-15} (seconds)
const double LC_MEAN_LIFE = 202.4E-15; // \pm 3.1 10^{-15} (seconds)

using namespace std;

double dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{ // calculate the signed DCA value
  TVector3 posDiff = pos - vertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
  return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

double dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
 TVector3 newPos(pos);
 newPos.SetZ(0);

 TVector3 newP(p);
 newP.SetZ(0);

 TVector3 newVertex(vertex);
 newVertex.SetZ(0);

 TVector3 posDiff = newPos - newVertex;
 double sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
 return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
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

void correct_D0_verticies(erhic::EventPythia* py_evt)
{ // patches the D0 vertex issue where with BeAGLE, D0 verticies are much smaller than they should be. Manually generates random and apporpraite verticies. Modifies verticies for child particles as well

  TF1* func_D0_decay_time = new TF1("func_D0_decay_length", "exp(-x/([0]*[1]))", 0, 1E-10); // output unitless, x units of 10^{-15} s
  // [0] = gamma [unitless], [1] = MEAN_LIFE [seconds]

  //changing D0 vertecies to patch zero vertex issue
  for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
  {

    erhic::ParticleMC* part = py_evt->GetTrack(ipart);

    if(abs(part->Id()) == 421)
    {
      if (verbosity>1)
      {
        cout<<"This is a D0 ================================================================"<<endl;
        (part->GetVertex()).Print();
        cout<<(part->GetVertex()).Mag()<<endl;
      }

      TLorentzVector track_mom4_true = part->Get4Vector();

      //calculate new vertex coords
      double_t velocity_mag = LIGHT_SPEED * (sqrt(1 - (1 / pow(track_mom4_true.Gamma(), 2)))); // units of m/s
      //cout<<"velocity: "<<velocity_mag<<endl;
      func_D0_decay_time->SetParameters(track_mom4_true.Gamma(), D0_MEAN_LIFE);  //gamma, D0_MEAN_LIFE
      double_t decay_length = (func_D0_decay_time->GetRandom()) * velocity_mag * 1E3; // units of mm
      //cout<<"new decay length(mm): "<<decay_length<<endl;
      double_t decay_dir_phi = track_mom4_true.Phi();
      double_t decay_dir_theta = track_mom4_true.Theta();

      //make new vertex
      TVector3 new_vtx_true;
      new_vtx_true.SetMagThetaPhi(decay_length, decay_dir_theta, decay_dir_phi);

      //set new vertex for all D0 children
      erhic::ParticleMC* child_part;
      for (int ichild = 0; ichild < part->GetNChildren(); ichild++)
      {
        if (part->GetChild1Index() == 0) break; // breaks if part has no child particles
        child_part = py_evt->GetTrack(part->GetChild1Index() + ichild - 1);
        //cout<<"supposed child particle. parent id: "<<child_part->GetParentId()<<" child id: "<<child_part->Id()<<" track index: "<<(part->GetChild1Index() + ichild)<<endl;
        child_part->SetVertex(new_vtx_true);
        //cout<<"child particle vertex corrected:"<<endl;
      }
      // TODO confirm setting
    }
  }
}

void correct_Lc_verticies(erhic::EventPythia* py_evt)
{ // patches the D0 vertex issue where with BeAGLE, D0 verticies are much smaller than they should be. Manually generates random and apporpraite verticies. Modifies verticies for child particles as well

  TF1* func_Lc_decay_time = new TF1("func_Lc_decay_length", "exp(-x/([0]*[1]))", 0, 1E-10); // output unitless, x units of 10^{-15} s
  // [0] = gamma [unitless], [1] = MEAN_LIFE [seconds]

  //changing Lc vertecies to patch zero vertex issue
  for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
  {

    erhic::ParticleMC* part = py_evt->GetTrack(ipart);

    if(abs(part->Id()) == 4122)
    {
      if (verbosity>1)
      {
        cout<<"This is a Lc ================================================================"<<endl;
        (part->GetVertex()).Print();
        cout<<(part->GetVertex()).Mag()<<endl;
      }

      TLorentzVector track_mom4_true = part->Get4Vector();

      //calculate new vertex coords
      double_t velocity_mag = LIGHT_SPEED * (sqrt(1 - (1 / pow(track_mom4_true.Gamma(), 2)))); // units of m/s
      //cout<<"velocity: "<<velocity_mag<<endl;
      func_Lc_decay_time->SetParameters(track_mom4_true.Gamma(), LC_MEAN_LIFE);  //gamma, LC_MEAN_LIFE
      double_t decay_length = (func_Lc_decay_time->GetRandom()) * velocity_mag * 1E3; // units of mm
      //cout<<"new decay length(mm): "<<decay_length<<endl;
      double_t decay_dir_phi = track_mom4_true.Phi();
      double_t decay_dir_theta = track_mom4_true.Theta();

      //make new vertex
      TVector3 new_vtx_true;
      new_vtx_true.SetMagThetaPhi(decay_length, decay_dir_theta, decay_dir_phi);

      //set new vertex for all Lc children
      erhic::ParticleMC* child_part;
      for (int ichild = 0; ichild < part->GetNChildren(); ichild++)
      {
        if (part->GetChild1Index() == 0) break; // breaks if part has no child particles
        child_part = py_evt->GetTrack(part->GetChild1Index() + ichild - 1);
        //cout<<"supposed child particle. parent id: "<<child_part->GetParentId()<<" child id: "<<child_part->Id()<<" track index: "<<(part->GetChild1Index() + ichild)<<endl;
        child_part->SetVertex(new_vtx_true);
        //cout<<"child particle vertex corrected:"<<endl;
      }
      // TODO confirm setting
    }
  }
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

    // if D0 vertex issue needs to be corrected, 1 if yes, 0 if no
    int do_correct_vertex;

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

    TH2D* h2d_D0_pt_vs_eta_gen[Q2bin][xbin];
    TH2D* h2d_D0_z_vs_eta_gen[Q2bin][xbin];

    int PROCESS_INDEX;
    TH2D* h2d_ztheo_vs_zjet[Q2bin][etabin][processbin];

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

      PROCESS_INDEX = -9999;

      do_correct_vertex = 0; // corrects D0 verticies

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

          h2d_D0_pt_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_D0_pt_vs_eta_gen_%d_%d",iQ2,ix),"D0 pt vs eta",100,0,10,160,-4,4);
          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Sumw2();

          h2d_D0_z_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_D0_z_vs_eta_gen_%d_%d",iQ2,ix),"D0 z vs eta",100,0,1,160,-4,4);
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Sumw2();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          for (int iprocess = 0; iprocess < processbin; ++iprocess)
          {
            h2d_ztheo_vs_zjet[iQ2][ieta][iprocess] = new TH2D(Form("h2d_z_frag_%d_%d_%d", iQ2, ieta, iprocess),"z_{theo} vs z_{jet}",100,0,1,100,0,1);
            h2d_ztheo_vs_zjet[iQ2][ieta][iprocess]->Sumw2();
          }
        }
      }

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

          delete h2d_D0_pt_vs_eta_gen[iQ2][ix];
          delete h2d_D0_z_vs_eta_gen[iQ2][ix];
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          for (int iprocess = 0; iprocess < processbin; ++iprocess)
          {
            delete h2d_ztheo_vs_zjet[iQ2][ieta][iprocess];
          }
        }
      }

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

      PROCESS_INDEX = -9999;

      do_correct_vertex = 0;

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

          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Reset("ICESM");
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Reset("ICESM");
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          for (int iprocess = 0; iprocess < processbin; ++iprocess)
          {
            h2d_ztheo_vs_zjet[iQ2][ieta][iprocess]->Reset("ICESM");
          }
        }
      }

    }

    void SetProcessID(const int id)
    {
      cout<<"process id = "<<id<<endl;
      if (id == 99)
      {
        PROCESS_INDEX = 0;
      } else if (id == 135 || id == 136)
      {
        PROCESS_INDEX = 1;
      } else
      {
        PROCESS_INDEX = -9999;
      }
    }

    void SetLowP(const float mom_thr) { TRK_P_LO = mom_thr; }

    void SetDCACuts()
    {
      // single
      TRK_DCA = -9999; //0.03; // in units of mm

      // pair
      PAIR_DCA = 0.13; //0.12; // 120um in unit of mm
      DECAY_L = 0.04;
      D0_DCA = -9999; //0.01;
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

    void SetDoCorrectVertex(int _do_correct) { do_correct_vertex = _do_correct; };

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

    void FillGenKin(erhic::EventPythia* py_evt)
    {
      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (abs(part->Id())!=421) continue; // FIX ME: for now fill all the D0 (maybe fill the D0->Kpi later)
        // NB: however at generation level there shouldn't be kin difference for different decay channels (it can be different once including the detector acceptance)

        TLorentzVector D0_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(D0_mom4_gen)/(nu_true*hadron_beam.M());

        if (verbosity>1) std::cout << "D0 with pt " << D0_mom4_gen.Pt() << " z " << frag_z << " eta " << D0_mom4_gen.PseudoRapidity() << std::endl;

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (verbosity>1) std::cout << "Q2, x " << Q2_true << ", " << x_true << " bin " << iQ2bin << ", " << ixbin << std::endl;

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_D0_pt_vs_eta_gen[iQ2bin][ixbin]->Fill(D0_mom4_gen.Pt(),D0_mom4_gen.PseudoRapidity());
          h2d_D0_z_vs_eta_gen[iQ2bin][ixbin]->Fill(frag_z,D0_mom4_gen.PseudoRapidity());

          h2d_D0_pt_vs_eta_gen[Q2bin-1][xbin-1]->Fill(D0_mom4_gen.Pt(),D0_mom4_gen.PseudoRapidity());
          h2d_D0_z_vs_eta_gen[Q2bin-1][xbin-1]->Fill(frag_z,D0_mom4_gen.PseudoRapidity());
        }
      }
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

      if (do_correct_vertex == 1) correct_D0_verticies(py_evt);

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
        if (SMEAR_OPTION==4)
        {
          track_mom4_reco = smearMomATHENA(track_mom4_true);
          track_vtx_reco = smearPosATHENA(track_mom4_true.Vect(), track_vtx_true);
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

      cout<<"D0 pt "<<pair.Pt()<<" eta "<<pair.PseudoRapidity()<<endl;

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
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_D0_pt_vs_eta[iQ2bin][ixbin]->Fill(pair.Pt(),pair.PseudoRapidity());
          h2d_D0_z_vs_eta[iQ2bin][ixbin]->Fill(frag_z,pair.PseudoRapidity());

          h2d_D0_pt_vs_eta[Q2bin-1][xbin-1]->Fill(pair.Pt(),pair.PseudoRapidity());
          h2d_D0_z_vs_eta[Q2bin-1][xbin-1]->Fill(frag_z,pair.PseudoRapidity());
        }

        cout<<"PROCESS_INDEX = "<<PROCESS_INDEX<<endl;
        if (iQ2bin>=0)
        {
          h2d_ztheo_vs_zjet[Q2bin-1][ietabin-1][processbin-1]->Fill((pair.Vect()).Dot(quark_p.Vect())/(quark_p.Vect()).Dot(quark_p.Vect()),frag_z);

          if (PROCESS_INDEX >= 0)
          {
            h2d_ztheo_vs_zjet[iQ2bin][ietabin][PROCESS_INDEX]->Fill((pair.Vect()).Dot(quark_p.Vect())/(quark_p.Vect()).Dot(quark_p.Vect()),frag_z);
          }
        }

      }
    }

    void FillD0Pairs()
    {
      for (int ineg = 0; ineg < negl_id_true.size(); ++ineg)
      {
        for (int ipos = 0; ipos < posl_id_true.size(); ++ipos)
        {
          // ==========================
          //    decay topology cut
          // ==========================
          // pair DCA < cut value
          TVector3 dca_pair = negl_vtx_reco[ineg]-posl_vtx_reco[ipos];
          dca_pair.SetZ(0);
          // cout << "dca.pair " << dca_pair.Mag() << endl;
          if (PAIR_DCA>-99 && dca_pair.Mag()>PAIR_DCA) continue;
          // cout << "double check dca.pair " << dca_pair.Mag() << endl;

          // Decay length > cut value
          TVector3 decay_vtx_reco = (negl_vtx_reco[ineg]+posl_vtx_reco[ipos])*0.5;
          TVector3 decay_l = decay_vtx_reco-evt_vtx_reco;
          decay_l.SetZ(0);
          if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue;

          // D0 DCA > cut value
          TVector3 D0_vec = negl_p_reco[ineg].Vect()+posl_p_reco[ipos].Vect();
          double D0_dca = dcaXY(D0_vec, decay_vtx_reco, evt_vtx_reco);
          if (D0_DCA>-99 && fabs(D0_dca)<D0_DCA) continue;

          // D0 cos(theta) > cut value
          TVector3 D0_vecT = D0_vec;
          D0_vecT.SetZ(0);
          double D0_costheta = TMath::Cos(D0_vecT.Angle(decay_l));
          if (D0_COSTHETA>-99 && D0_costheta<D0_COSTHETA) continue;

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

          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Write();
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Write();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          for (int iprocess = 0; iprocess < processbin; ++iprocess)
          {
            h2d_ztheo_vs_zjet[iQ2][ieta][iprocess]->Write();
          }
        }
      }

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

    // if Lc vertex issue needs to be corrected, 1 if yes, 0 if no
    int do_correct_vertex;

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

    TH2D* h2d_Lc_pt_vs_eta_gen[Q2bin][xbin];
    TH2D* h2d_Lc_z_vs_eta_gen[Q2bin][xbin];

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

      do_correct_vertex = 0;

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

          h2d_Lc_pt_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_Lc_pt_vs_eta_gen_%d_%d",iQ2,ix),"Lc pt vs eta",100,0,10,160,-4,4);
          h2d_Lc_pt_vs_eta_gen[iQ2][ix]->Sumw2();

          h2d_Lc_z_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_Lc_z_vs_eta_gen_%d_%d",iQ2,ix),"Lc pt vs eta",100,0,1,160,-4,4);
          h2d_Lc_z_vs_eta_gen[iQ2][ix]->Sumw2();
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

          delete h2d_Lc_pt_vs_eta_gen[iQ2][ix];
          delete h2d_Lc_z_vs_eta_gen[iQ2][ix];
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

      do_correct_vertex = 0;

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

          h2d_Lc_pt_vs_eta_gen[iQ2][ix]->Reset("ICESM");
          h2d_Lc_z_vs_eta_gen[iQ2][ix]->Reset("ICESM");
        }
      }
    }

    void SetLowP(const float mom_thr) { TRK_P_LO = mom_thr; }

    void SetDCACuts()
    {
      TRK_DCA = -9999; // by default no cut

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

    void SetDoCorrectVertex(int _do_correct) { do_correct_vertex = _do_correct; };

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

    void FillGenKin(erhic::EventPythia* py_evt)
    {
      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (abs(part->Id())!=4122) continue; // FIX ME: for now fill all the D0 (maybe fill the Lc->pKpi later)

        TLorentzVector Lc_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(Lc_mom4_gen)/(nu_true*hadron_beam.M());

        if (verbosity>1) std::cout << "Lc with pt " << Lc_mom4_gen.Pt() << " z " << frag_z << " eta " << Lc_mom4_gen.PseudoRapidity() << std::endl;

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_Lc_pt_vs_eta_gen[iQ2bin][ixbin]->Fill(Lc_mom4_gen.Pt(),Lc_mom4_gen.PseudoRapidity());
          h2d_Lc_z_vs_eta_gen[iQ2bin][ixbin]->Fill(frag_z,Lc_mom4_gen.PseudoRapidity());

          h2d_Lc_pt_vs_eta_gen[Q2bin-1][xbin-1]->Fill(Lc_mom4_gen.Pt(),Lc_mom4_gen.PseudoRapidity());
          h2d_Lc_z_vs_eta_gen[Q2bin-1][xbin-1]->Fill(frag_z,Lc_mom4_gen.PseudoRapidity());
        }
      }
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

      if (do_correct_vertex == 1) correct_Lc_verticies(py_evt);

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
        if (SMEAR_OPTION==4)
        {
          track_mom4_reco = smearMomATHENA(track_mom4_true);
          track_vtx_reco = smearPosATHENA(track_mom4_true.Vect(), track_vtx_true);
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
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_Lc_pt_vs_eta[iQ2bin][ixbin]->Fill(trip.Pt(),trip.PseudoRapidity());
          h2d_Lc_z_vs_eta[iQ2bin][ixbin]->Fill(frag_z,trip.PseudoRapidity());

          h2d_Lc_pt_vs_eta_gen[iQ2bin][ixbin]->Fill(trip.Pt(),trip.PseudoRapidity());
          h2d_Lc_z_vs_eta_gen[iQ2bin][ixbin]->Fill(frag_z,trip.PseudoRapidity());

          h2d_Lc_pt_vs_eta[Q2bin-1][xbin-1]->Fill(trip.Pt(),trip.PseudoRapidity());
          h2d_Lc_z_vs_eta[Q2bin-1][xbin-1]->Fill(frag_z,trip.PseudoRapidity());

          h2d_Lc_pt_vs_eta_gen[Q2bin-1][xbin-1]->Fill(trip.Pt(),trip.PseudoRapidity());
          h2d_Lc_z_vs_eta_gen[Q2bin-1][xbin-1]->Fill(frag_z,trip.PseudoRapidity());
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
            dca_pair1.SetZ(0);
            dca_pair2.SetZ(0);
            dca_pair3.SetZ(0);
            double dca_pair[3] = {dca_pair1.Mag(),dca_pair2.Mag(),dca_pair3.Mag()};
            if (PAIR_DCA>-99 && TMath::MaxElement(3,dca_pair)>PAIR_DCA) continue;

            // Decay length > cut value
            TVector3 decay_vtx_reco = (negl_vtx_reco[ineg]+posl_vtx_reco[ipos1]+posl_vtx_reco[ipos2])*(1./3);
            TVector3 decay_l = decay_vtx_reco-evt_vtx_reco;
            decay_l.SetZ(0);
            if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue; //commented for eA

            // Lc DCA < cut value
            TVector3 Lc_vec = negl_p_reco[ineg].Vect()+posl_p_reco[ipos1].Vect()+posl_p_reco[ipos2].Vect();
            double Lc_dca = dcaXY(Lc_vec, decay_vtx_reco, evt_vtx_reco);
            if (Lc_DCA>-99 && fabs(Lc_dca)>Lc_DCA) continue;

            // Lc cos(theta) > cut value
            TVector3 Lc_vecT = Lc_vec;
            Lc_vecT.SetZ(0);
            double Lc_costheta = TMath::Cos(Lc_vecT.Angle(decay_l));
            if (Lc_COSTHETA>-99 && Lc_costheta<Lc_COSTHETA) continue;

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

          h2d_Lc_pt_vs_eta_gen[iQ2][ix]->Write();
          h2d_Lc_z_vs_eta_gen[iQ2][ix]->Write();
        }
      }
    }
};

void D0_tree_patch_zmod(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hist.root", int nevt = 0, const int smear_option = 0, const int Bfield_type = 0, const int PID_option = 0)
{ // smear_2nd_vtx & momentum: 0--no smearing, 1--DM smearing, 2--LBL smearing, 3--Hybrid smearing, 4--ATHENA smearing
  // Bfield_type: 0--Barbar, 1--Beast
  // 0--no hID (but with eID), 1--PID with no low momentum cutoff, 2--PID with low momentum cutoff & some mis-identified pi, K, 3--PID with low momentum cutoff & all identified pi, K
  // DCA_cut: 0--no cut, 1--cut on DCA
  // do_correct_vertex: 0--no correction applied, 1--D0 and Lc vertex decay length recalculated, modifications applied to input root files and data used in above analysis

  // PDG data table
  pdg = new TDatabasePDG();

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  // erhic::ParticleMC *child(NULL);

  //Load ROOT File
  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  // ATHENA smeaing
  TFile* f_athena_tracks = new TFile("ATHENA_Resolutions_r.root","READ");
  TFile* f_athena_vertex = new TFile("VertexRes_ATHENA.root","READ");
  if (smear_option==4)
  {
    cout << "Setup ATHENA smearing parameters" << endl;
    setup_ATHENA_single_track_smearing(f_athena_tracks);
    setup_ATHENA_PV_smearing(f_athena_vertex);
  }

  //Define Some Variables
  Float_t Q2(0);
  Int_t nParticles(0);

  D0_reco ana_D0;
  ana_D0.SetLowP(0.1);
  ana_D0.SetDCACuts();
  ana_D0.SetIDCuts(PID_option);
  ana_D0.SetSmearType(smear_option);
  ana_D0.SetBFieldType(Bfield_type);
  // ana_D0.SetDoCorrectVertex(do_correct_vertex);

  Lc_reco ana_Lc;
  ana_Lc.SetDCACuts();
  ana_Lc.SetIDCuts(PID_option);
  ana_Lc.SetSmearType(smear_option);
  ana_Lc.SetBFieldType(Bfield_type);
  // ana_Lc.SetDoCorrectVertex(do_correct_vertex);

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

    TVector3 evt_vtx_true(0,0,0); // not sure how to get event vertex
    TVector3 evt_vtx_reco = evt_vtx_true;
    if (smear_option==4)
    {
      float multi_charged_hadron_35 = 0;
      float multi_charged_hadron_10 = 0;
      for (int ipart = 0; ipart < nParticles; ++ipart)
      {
        particle = event->GetTrack(ipart);
        if (particle->GetStatus()!=1) continue; // only loop through final stable particles

        if (abs(particle->Id())!=2212 && abs(particle->Id())!=211 && abs(particle->Id())!=321) continue; // proton, pion, kaon
        if(fabs(particle->GetEta())<3.5 && particle->GetP()>0.5) multi_charged_hadron_35++;
        if(fabs(particle->GetEta())<1.0 && particle->GetP()>0.5) multi_charged_hadron_10++;
      }
      evt_vtx_reco = smearPVTATHENA(evt_vtx_true, multi_charged_hadron_35); // smear transverse direction
      // comment out the longgitunal smearing for now
      // evt_vtx_reco = smearPVZATHENA(evt_vtx_true, multi_charged_hadron_10); // smear longitudinal direction
    }
    if (evt_vtx_reco.X()<-999 || evt_vtx_reco.Y()<-999)
    {
      // cout << evt_vtx_reco.X() << ", " << evt_vtx_reco.Y() << ", " << evt_vtx_reco.Z() << " mm" << endl;
      continue; // No valid PV smearing in transverse direction, jump to next event
    }

    // used for z defintion study, cutting on process id
    ana_D0.SetProcessID(event->GetProcess());

    ana_D0.SetVectTrue(evt_vtx_true);
    ana_D0.SetVectReco(evt_vtx_reco); // NB: no smear on primary vertex yet

    ana_D0.SetQ2True(event->GetQ2());
    ana_D0.SetXTrue(event->GetX());

    ana_D0.SetNuTrue(event->GetNu());

    ana_D0.FillGenKin(event);

    ana_D0.FillSingleTracks(event);
    ana_D0.FillD0Pairs();

    ana_Lc.SetVectTrue(evt_vtx_true);
    ana_Lc.SetVectReco(evt_vtx_reco); // NB: no smear on primary vertex yet

    ana_Lc.SetQ2True(event->GetQ2());
    ana_Lc.SetXTrue(event->GetX());

    ana_Lc.SetNuTrue(event->GetNu());

    ana_Lc.FillGenKin(event);

    ana_Lc.FillSingleTracks(event);
    ana_Lc.FillLcTriplets();
  }

  TFile* fout = new TFile(outFile,"recreate");

  ana_D0.Write();
  ana_Lc.Write();

  fout->Write();
}
