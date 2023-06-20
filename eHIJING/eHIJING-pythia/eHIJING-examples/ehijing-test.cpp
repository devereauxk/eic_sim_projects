#include "Pythia8/Pythia.h"
using namespace Pythia8;
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <unistd.h>

// A separate Pythia instance that only handles hadronization
class hadronizer{
public:
    // Constructor, set Pythia, random generator
    hadronizer():pythia(),rd(),gen(rd()),dist(0.,1.){
      pythia.readString("ProcessLevel:all = off");
      pythia.readString("Print:quiet = on");
      pythia.readString("Next:numberShowInfo = 0");
      pythia.readString("Next:numberShowProcess = 0");
      pythia.readString("Next:numberShowEvent = 0");
      // prton tune and PDF set, please check
      pythia.readString("Tune:pp = 19");
      pythia.readString("PDF:pSet = 12");
      // These two parameters have something to do with the
      // Lund String interative breaking conditions
      // ** We have changed stopMass=0.0 GeV, different from Pythia8 Default
      // too match HERMES FF measurements of pion and Kaon.
      pythia.readString("StringFragmentation:stopMass = 0.0");
      pythia.readString("HadronLevel:mStringMin = 0.5");
      // Set some hadronic & decay specific channls
      pythia.readString("HadronLevel:Decay = on");
      pythia.readString("111:mayDecay=off");
      pythia.readString("211:mayDecay=off");
      pythia.readString("321:mayDecay=off");
      pythia.init();
    }

    // This function takes the shower PythiaIn (ep-shower with recoil particles)
    // and assume it fragments in an evneronment of nucleus with proton Z and mass A
    std::vector<Particle> hadronize(Pythia & pythiaIn, int Z, int A){
       // proton fraction;
       double ZoverA = Z*1./A;
       std::vector<Particle> FinalParticles; FinalParticles.clear();
       // Get the initial hard parton ID
       int hardid = pythiaIn.event[5].id();
       // Reset the hadronizer Pythia
       pythia.event.reset();
       // Loop over the partons in the shower PythiaIn, and put them in the hadronizer
       for (int i=0; i<pythiaIn.event.size(); i++){
         auto & p = pythiaIn.event[i];
	 // Find the final-state parton stuff (this drops the deflected lepton in the event)
         if (! ( p.isFinal() && p.isParton() ) ) continue;
         // These are di-quark remnants of the proton
         if (p.status()==63 && 1000<p.idAbs() && p.idAbs()<3000) {
                 // valence stuff, the remnants will contain the rest flavor compoennt.
                 // note that the hard quark has already been sampled accorrding to the
                 // the isospin content of the nuclear PDF;
		 // *** However, the remanent is generated assuming the rest stuff comes
		 // from a proton. Therefore, we need to resample it according to the Z/A
		 // ratio this nuclei
                 // 1) decide wither it is from a neutron or proton
                 if (dist(gen) < ZoverA) { // From a proton 2212
                     if (hardid==1) { // produce 2203
                       p.id(2203);
                     }
                     if (hardid==2) { // produce 2101 and 2103 with ratio 3:1
                       if (dist(gen) < 0.75) p.id(2101);
                       else p.id(2103);
                     }
                 } else { // From a neutron 2112
                     if (hardid==1) { // produce 2101 and 2103 with ratio 3:1
                       if (dist(gen) < 0.75) p.id(2101);
                       else p.id(2103);
                    }
                     if (hardid==2) { // produce 1103
                       p.id(1103);
                     }
                 }
         }
	 // For other partons, just put it in the shower
         pythia.event.append(p.id(), 23, p.col(), p.acol(),
                       p.px(), p.py(), p.pz(), p.e(), p.m());
       }
       // Do Hadronization
       pythia.next();
       // Return only final-state particles.
       for (int i=0; i<pythia.event.size(); i++) {
         auto & p = pythia.event[i];
         if (p.isFinal()) FinalParticles.push_back(p);
      }
      return FinalParticles;
    }

private:
    Pythia pythia;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dist;
};

// Helper function to set Pythia settings
template <class T>
void add_arg(Pythia & pythia, std::string name, T value){
  std::stringstream ss;
  ss << name << " = " << value;
  std::cout << ss.str() << std::endl;
  pythia.readString(ss.str());
}

// You can put whatever triggers of Pythia events here for e-p.
bool trigger(Pythia & pythia) {
    Vec4 pProton = pythia.event[1].p(); // four-momentum of proton
    Vec4 peIn    = pythia.event[4].p(); // incoming electron
    Vec4 peOut   = pythia.event[6].p(); // outgoing electron
    Vec4 pGamma = peIn - peOut; // virtual boson photon/Z^0/W^+-

    // Q2, W2, Bjorken x, y.
    double Q2 = - pGamma.m2Calc(); // hard scale square
    double W2 = (pProton + pGamma).m2Calc();
    double x  = Q2 / (2. * pProton * pGamma); // Bjorken x
    double nu = pGamma.e();
    double y     = (pProton * pGamma) / (pProton * peIn);

    return (0.01<y) && (y<0.95) && (1.0<Q2);
}

// old output function
/*
void Output(Pythia & pythia, std::vector<Particle> & plist, std::ofstream & F){
    // Compute four-momenta of proton, electron, virtual
    Vec4 pProton = pythia.event[1].p(); // four-momentum of proton
    Vec4 peIn    = pythia.event[4].p(); // incoming electron
    Vec4 peOut   = pythia.event[6].p(); // outgoing electron
    Vec4 pGamma = peIn - peOut; // virtual boson photon/Z^0/W^+-
    double Q2 = - pGamma.m2Calc(); // hard scale square
    double xB  = Q2 / (2. * pProton * pGamma); // Bjorken x
    Vec4 pCoM = pGamma + pProton;
    double nu = pGamma.e();
    double theta = - pGamma.theta();
    double phi = - pGamma.phi();
    F << "# " << Q2 << " " << xB << std::endl;
    for (auto & p : plist) {
        if (p.isFinal() && p.isHadron()){
            // Fixed target frame
            double z = p.e()/nu;
            auto prot = p.p();
            prot.rot(0, phi);
            prot.rot(theta, 0);
            double pT = p.pT();
	    double kT = prot.pT();
	    F << " " << p.id() << " "
		 // photon-z frame
		 << z << " " << kT << " " << prot.phi() << " "
		 // Fixed target frame
		 << p.pT() << " " << p.y() << " " << p.phi() << std::endl;
	 }
    }
}
*/

// modified output function
int evtn = 0;
void Output(Pythia & pythia, std::vector<Particle> & plist, std::ofstream & F){
    // Compute four-momenta of proton, electron, virtual
    Vec4 pProton = pythia.event[1].p(); // four-momentum of proton
    Vec4 peIn    = pythia.event[4].p(); // incoming electron
    Vec4 peOut   = pythia.event[6].p(); // outgoing electron
    Vec4 pGamma = peIn - peOut; // virtual boson photon/Z^0/W^+-
    double Q2 = - pGamma.m2Calc(); // hard scale square
    double xB  = Q2 / (2. * pProton * pGamma); // Bjorken x
    Vec4 pCoM = pGamma + pProton;
    double nu = pGamma.e();
    double theta = - pGamma.theta();
    double phi = - pGamma.phi();

    for (auto & p : plist) {
        if (p.isFinal())
        {
            int id = p.id();
            double charge = p.charge();

            // kinematics in target rest frame (to be boosted in analysis script)
      	    F << evtn << "," << p.id() << "," << p.charge() << ","
            << p.px() << "," << p.py() << "," << p.pz() << "," << p.m() << "," << Q2 << "," << xB << std::endl;

	      }
    }
    evtn++;
}

// low-Q2 medium correction (a Monte Carlo version of the the modified FF model)
class Modified_FF{
public:
  Modified_FF(int mode_, int Z_, int A_, double K, double n, double lambda, std::string TablePath)
    : mode(mode_), Z(Z_), A(A_), ZoverA(Z*1./A),
      Coll(K, n, lambda),
      eHIJING_Geometry(A, Z),
      rd(), gen(rd()), dist(0.,1.)	{
    Coll.Tabulate(TablePath);
  };
  void sample_FF_partons(Event & event);

private:
    int mode, Z, A;
    double ZoverA;
    EHIJING::MultipleCollision Coll;
    EHIJING::NuclearGeometry eHIJING_Geometry;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dist;
};

// The main routine:
int main(int argc, char *argv[]) {

    // Read commandline arguments
    // number of events, atomic Z and A
    int nEvent = atoi(argv[1]);
    int Z = atoi(argv[2]);
    int A = atoi(argv[3]);
    // build the nucleus ID used by Pythia & the PDF (the isospin effect)
    int inuclei = 100000000
                +   Z*10000
                   +   A*10;
    // Shadowing effect:
    //   nPDFset=0: only isospin
    //   nPDFset = 1: EPS09 LO
    //   nPDFset = 2: EPS09 NLO
    //   nPDFset = 3: EPPS16 NLO
    // We will use only isospin for deuteron,
    // and EPPS16 NLO for heavier nucleus
    int nPDFset = (A>2)?3:0;
    if (A == 238)
    {
      nPDFset = 2;
    }
    //nPDFset = 0; // override nPDFset = 0 for isospin study REMEMBER TO CHANGE BACK

    // mode=0: higher-twist, in the soft-gluon-emission limit
    // mode=1: generalized higher-twist, in the soft-gluon-emission limit
    int mode = atof(argv[4]);

    // K-factor of the gluon distribution.
    double K = atof(argv[5]);

    // eHIJING table path
    std::string TablePath(argv[6]);

    // Pythia8+other eHijing params configurations
    auto configfile = std::string(argv[8]);

    // header/folder of the output
    // use CPU process id used to name the output file
    auto header = std::string(argv[7]);
    int process_id = getpid();
    std::stringstream sout;
    //sout << header << "/" << process_id << ".dat";
    sout << header << "/" << "events" << ".dat";
    std::ofstream fout(sout.str());

    // initialize the hadronizer instance:
    hadronizer HZ;

    // Initialize the eHIJING-pythia for high-Q parton shower in medium
    Pythia pythia;        // Generator
    Event& event = pythia.event; // Event record
    pythia.readFile(configfile); // read settings
    add_arg<int>(pythia, "eHIJING:Mode", mode);
    add_arg<int>(pythia, "PDF:nPDFSetA", nPDFset);
    add_arg<int>(pythia, "PDF:nPDFBeamA", inuclei);
    add_arg<int>(pythia, "eHIJING:AtomicNumber", A);
    add_arg<int>(pythia, "eHIJING:ChargeNumber", Z);
    add_arg<double>(pythia, "eHIJING:Kfactor", K);
    add_arg<std::string>(pythia, "eHIJING:TablePath", TablePath);
    pythia.init();

    // initialize the modified FF class
    Modified_FF MFF(mode, Z, A, K, pythia.settings.parm("eHIJING:xG-n"),
                               pythia.settings.parm("eHIJING:xG-lambda"),
			       TablePath);

    // Begin event loop.
    int Ntriggered = 0;
    int Ntotal = 0, Nfailed=0;
    while(Ntriggered<nEvent){
        Ntotal ++;
        if (!pythia.next()) {
            Nfailed ++;
            continue;
            // count failed events
        };

        if (!trigger(pythia)) continue; // only study triggered events
        Ntriggered ++;
        if (Ntriggered%1000==0)
        std::cout << "# of trigged events: " << Ntriggered << std::endl;

        // Modify the final shower with low-Q2 medium corrections
        MFF.sample_FF_partons(event);

        // put the parton level event into a separate hadronizer
        auto event2 = HZ.hadronize(pythia, Z, A);

        // output
        Output(pythia, event2, fout);
    }
    // Check the trigger rate
    std::cout << "Trigger Rate = " << Ntriggered*100./Ntotal << "%" << std::endl;
    // Check the rate of failed events
    std::cout << "Failed Rate = " << Nfailed*100./Ntotal << "%" << std::endl;

    // Done.
    return 0;
}


// The realization of the modified FF
void Modified_FF::sample_FF_partons(Event & event){
    Vec4 pProton = event[1].p(); // four-momentum of proton
    Vec4 peIn    = event[4].p(); // incoming electron
    Vec4 peOut   = event[6].p(); // outgoing electron
    Vec4 pGamma = peIn - peOut; // virtual boson photon/Z^0/W^+-
    double Q20 = - pGamma.m2Calc(); // hard scale square
    double Q0 = std::sqrt(Q20);
    double xB  = Q20 / (2. * pProton * pGamma); // Bjorken x
    double nu = pGamma.e();
    double W2 = (pProton + pGamma).m2Calc();
    auto & hardP = event[5];
    
    // Use fixed coupling at Qs of this event for anything medium-induced below Qs
    double kt2max_now = event.SeparationScale();
    double alpha_fix = EHIJING::alphas(kt2max_now);
    double alphabar = alpha_fix * EHIJING::CA/M_PI;
    double emin = .2;

    // This will hold the gluons emitted from the medium-contrbuted FF in this stage
    // as well as recoil partons that keeps the entire system color neutral.
    std::vector<Particle> new_particles; new_particles.clear();

    // For each parton from the high-Q shower, sample gluons and recoils from low-Q medium
    // corrections according to either HT or Generlized HT in the soft gluon limit
    for (int i=0; i<event.size(); i++) {

        // only find final-state quarks and gluons with energy above
        // 2*emin in the nuclear rest frame
        auto & p = event[i];
        auto abspid = p.idAbs();
        if (!p.isFinal() ) continue;
        if (p.e()<2*emin) continue;

        // we are skipping modificationsof heavy quarks right now!!!
	bool isLightParton = (abspid==1) || (abspid==2) || (abspid==3) || (abspid==21);
        if (!isLightParton) continue;

        // Get its collision history
        std::vector<double> qt2s = p.coll_qt2s(), ts = p.coll_ts(), phis = p.coll_phis();
        // if no collision, nothing to do
	int Ncolls = ts.size();
        if (Ncolls==0) continue;

     
	double vx = p.px()/p.e(), vy = p.py()/p.e(), vz = p.pz()/p.e();
        double L = eHIJING_Geometry.compute_L(event.Rx(), event.Ry(), event.Rz(), vx, vy, vz);
        double TA = eHIJING_Geometry.compute_TA(event.Rx(), event.Ry(), event.Rz(), vx, vy, vz);

	double sumq2 = 0.; // useful quantity for H-T approach
        for (int i=0; i<Ncolls; i++) sumq2 += qt2s[i];
        if (sumq2<1e-9) continue; // negelect too soft momentum kicks
	   

          // tauf ordered fragmentation gluon
	  // A very large cut off, since the LPM effect will effective regulate the tauf divergence
          double taufmax = p.e()/EHIJING::mu2;
	  // Start from the minimum tauf ~ 1/Qs, which is the separation scale;
          double taufmin = 1.0/std::sqrt(event.SeparationScale());
	  double tauf = taufmin;

	  double acceptance = 0;
          // z, kt2, lt2, phi_k, and dphiqk = phi_q - phik
	  double z, kt2, lt2, phik, dphiqk;

	  // holds fragmented gluons and recoiled beam remnants
          std::vector<Particle> frag_gluons, recoil_remnants;
          frag_gluons.clear();
          recoil_remnants.clear();

	  // Formation time loop
          while(tauf < taufmax && p.e()>2*emin){
              double zmin = std::min(emin / p.e(), .4);
              double zmax = 1.-zmin;
              if (zmax<zmin) break;
              double maxlogz =  std::log(zmax/zmin);
              double maxdiffz = 1./zmin - 1./zmax + 2.*maxlogz;
              // step1: find the next tauf
              double r = dist(gen);
              if (mode==1){
                  double invrpower = alphabar * maxlogz * 4. * Ncolls;
                  double step_factor = std::pow(1./r, 1./invrpower);
                  tauf = tauf * step_factor;
              }
              else {
                  double coeff = alphabar * maxdiffz * 4. * sumq2 / 2. / p.e();
                  tauf = tauf + std::log(1./r)/coeff;
              }
              if (tauf > taufmax || tauf<taufmin) break;
              acceptance = 0.;
              if (mode==1) {
                  for (int j=0; j<Ncolls; j++){
                      double phase = (1.-std::cos(ts[j]/tauf));
                      double z1mz = tauf * qt2s[j] / 2. / p.e();
                      if (z1mz>.25) acceptance += phase * maxlogz;
                      else {
                          double dz = std::sqrt(.25 - z1mz);
                          double z1 = .5 - dz;
                          double z2 = .5 + dz;
                          if (z1>zmin) acceptance += phase * std::log(z1/zmin);
                          if (z2<zmax) acceptance += phase * std::log(zmax/z2);
                      }
                  }
                  acceptance /= (maxlogz * 2. * Ncolls);
              }
              else{
                  for (int j=0; j<Ncolls; j++) acceptance += qt2s[j]*(1.-std::cos(ts[j]/tauf));
                  acceptance /= (2.*sumq2);
              }
              if (acceptance < dist(gen)) continue;

              // step 2: sample z, which also determines kt2
              acceptance = 0.;
              if (mode==0) {
                  double N1 = 2*(1./zmin-2.);
                  double N2 = -4*std::log(2.*(1.-zmax));
                  double Ntot = N1+N2;
                  double r0 = N1/Ntot;

                  double acceptance = 0.;
                  while(acceptance<dist(gen)){
                      double r = dist(gen);
                      if (r<r0){
                          z = zmin/(1. - zmin*r*Ntot/2.);
                          acceptance = .5/(1.-z);
                      } else {
                          z = 1. - std::exp(-(r*Ntot - N1)/4.)/2.;
                          acceptance = .25/z/z;
                      }
                  }
                  kt2 = 2*(1.-z)*z*p.e()/tauf;
                  // reject cases where qt2>kt2 for mode=0
                  double Num=0., Den = 0.;
                  for (int j=0; j<Ncolls; j++) {
   			if (ts[j]<0) continue;

                      double q2 = qt2s[j], t = ts[j];
                      if (kt2>q2) Num += q2*(1.-std::cos(t/tauf));
                      Den += q2*(1.-std::cos(t/tauf));
                  }
                  if (Num/Den < dist(gen)) continue;
              }
              else {
                  bool ok=false;
                  double minimum_q2 = 2*emin/tauf;
                  for (int j=0; j<Ncolls; j++){
                      if (qt2s[j]>minimum_q2 && ts[j]>0) ok=true;
                  }
                  if (!ok) continue;
                  while(acceptance<dist(gen)){
                      z = zmin * std::pow(zmax/zmin, dist(gen));
                      kt2 = 2*(1.-z)*z*p.e()/tauf;
                      double num = 0.;
                      for (int j=0; j<Ncolls; j++)
                          if (kt2<qt2s[j] && ts[j]>0)
                              num += 1.;
                      acceptance = num / Ncolls;
                  }
              }

              // step 3 correct for the splitting function
              // correct for splitting function
              if (p.id()==21 && (1+std::pow(1.-z,3))/2.<dist(gen) ) continue;
              if (p.id()!=21 && (1+std::pow(1.-z,2))/2.<dist(gen) ) continue;

              // finally, sample phikT2 and compute the deflection of the hard parton lt2
              if (mode==0) {
		  // no particular angular structure in the H-T expansion
                  phik = 2*M_PI*dist(gen);
                  lt2 = kt2;
              }
              else {
                  double Psum = 0.;
                  std::vector<double> dP;
                  dP.resize(Ncolls);
                  for (int j=0; j<Ncolls; j++){
                      if (kt2<qt2s[j]) Psum += (1.-std::cos(ts[j]/tauf));
                      dP[j] = Psum;
                  }
                  for (int j=0; j<Ncolls; j++) dP[j]/=Psum;
                  double rc = dist(gen);
                  int choice = -1;
                  for (int j=0; j<Ncolls; j++){
                      if ( rc<dP[j] ) {
                          choice = j;
                          break;
                      }
                  }
                  // sample phik ~ (1+delta cos) / (1+delta^2 + 2 delta cos)
                  double delta = std::sqrt(kt2/qt2s[choice]);
                  acceptance = 0.;
                  while(acceptance < dist(gen)){
                      r = dist(gen);
                      dphiqk = 2.*std::atan(std::tan(M_PI/2.*r) * (delta+1)/(delta-1));
                      acceptance = (1+delta*std::cos(dphiqk))/2.;
                  }
                  phik = phis[choice] + ( (dist(gen)>.5)? dphiqk : (-dphiqk) );
		  // lt2 = |(K-q) + q| = |k-q|^2
                  lt2 = kt2 + qt2s[choice] + 2.*std::sqrt(kt2*qt2s[choice])*std::cos(dphiqk);
              }
              // The virtuality was not allowed to be overlap with the high-Q shower
              if (lt2>kt2max_now ) continue;

              // Now, there is a radiation,
              // Hard parton splits into p -> p-k & k
              double kt = std::sqrt(kt2), k0 = z*p.e();
              if (kt>k0) continue;
              double kz = std::sqrt(k0*k0-kt2);
              Vec4 kmu{kt*std::cos(phik), kt*std::sin(phik), kz, k0};
              kmu.rot(p.theta(), 0.);
              kmu.rot(0., p.phi());
              p.p(p.p()-kmu);
              p.e(std::sqrt(p.pAbs2()+p.m2()));
              kmu.e(std::sqrt(kmu.pAbs2()));

	      // the gluon can continue to collide 
               
	      std::vector<double> g_qt2s, g_ts, g_phis;
              double R0sq = event.Rx()*event.Rx() + event.Ry()*event.Ry() + event.Rz()*event.Rz();
              double V2 = vx*vx+vy*vy+vz*vz;
              double TwoRdotV = 2*(event.Rx()*vx + event.Ry()*vy + event.Rz()*vz);
              Coll.sample_all_qt2(21, kmu.e(), L, TA, xB, Q20, {R0sq, V2, TwoRdotV, A},
                                  g_qt2s, g_ts, g_phis);
              Vec4 Qtot{0.,0.,0.,0.};
	      double e0 = kmu.e();
	      for (int ig=0; ig<g_ts.size(); ig++){
		 double qt = std::sqrt(g_qt2s[ig]);
                 double phi = g_phis[ig];
		 Vec4 qmu{qt*std::cos(phi), qt*std::sin(phi), -qt*qt/4./e0, 0.0};
                 Qtot = Qtot + qmu; 
	      }
	      Qtot.rot(kmu.theta(), 0.);
              Qtot.rot(0., kmu.phi());
	      kmu = kmu + Qtot;
	      kmu.e(std::sqrt(kmu.pAbs2()));
	      kmu = kmu*e0/kmu.e();

              // update the color if it is a hard gluon
              // first, the spliting process
              int k_col, k_acol;
              // if the gluon forms inside the nuclei,
	      // we consider it will lose color correlation with the original parton, 
	      // and form a new string with beam remnant 
              {
                  Particle gluon = Particle(21, 201, i, 0, 0, 0,
                                  event.nextColTag(), event.nextColTag(),
                                  kmu, 0.0, 0);
                                  int qid, diqid;
                                  double mq, mdiq, mn;
                                  if (dist(gen) < ZoverA) { // diquark from a proton
                                      diqid = 2101; mdiq = 0.57933;
                                      qid = 2; mq = 0.33;
                                      mn = 0.93847;
                                      if (dist(gen) < 2./3.) { // take away a u
                                          qid = 2;
                                          if (dist(gen) < .75) diqid = 2101;
                                          else diqid = 2103;
                                      } else { // take away the d
                                          qid = 1;
                                          diqid = 2203;
                                      }
                                 }
                                 else { // diquark from a neutron
                                      diqid = 2101; mdiq = 0.57933;
                                      qid = 1; mq = 0.33;
                                      mn = 0.93957;
                                      if (dist(gen) < 2./3.) { // take away a d
                                          qid = 1;
                                          if (dist(gen) < .75) diqid = 2101;
                                          else diqid = 2103;
                                      } else { // take away the u
                                          qid = 2;
                                          diqid = 1103;
                                      }
                                 }
                                 double pabs = std::sqrt((mn*mn-std::pow(mq+mdiq,2))*(mn*mn-std::pow(mq-mdiq,2))) / (2.*mn);
                                 double costheta = dist(gen)*2.-1.;
                                 double sintheta = std::sqrt(std::max(1.-costheta*costheta,1e-9));
                                 double rphi = 2*M_PI*dist(gen);
                                 double Nqz = pabs*costheta,
                                        Nqx = pabs*sintheta*std::cos(rphi),
                                        Nqy = pabs*sintheta*std::sin(rphi);
                                 Vec4 pq  { Nqx,  Nqy,  Nqz, 0},
                                      pdiq{-Nqx, -Nqy,  -Nqz, 0};
                                 pq.e(std::sqrt(pq.pAbs2()+mq*mq));
                                 pdiq.e(std::sqrt(pdiq.pAbs2()+mdiq*mdiq));
                                 Particle recolQ = Particle(qid, 201, i, 0, 0, 0,
                                                   gluon.acol(), 0,
                                                   pq, mq, 0);
                                 Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0,
                                                   0, gluon.col(),
                                                   pdiq, mdiq, 0);
                                 recoil_remnants.push_back(gluon);
                                 recoil_remnants.push_back(recolQ);
                                 recoil_remnants.push_back(recoldiQ);
              }
          }
          // Now handles recoil and remannts
          // if there are radiations, recoil goes to radiations
          // else: goes to the hard quark 
	  int Nrad = frag_gluons.size();
               for (int j=0; j<Ncolls; j++){
                   double qT = std::sqrt(qt2s[j]), phiq = phis[j];
                   double qx = qT*std::cos(phiq);
                   double qy = qT*std::sin(phiq);
                   double qz = -qT*qT/4./p.e();   
                   Vec4 qmu{qx, qy, qz, 0};
                   qmu.rot(p.theta(), 0.);
                   qmu.rot(0., p.phi());
                   p.p(p.p()+qmu);
                   p.e(std::sqrt(p.pAbs2()+p.m2()));
                   int q_col, q_acol;
                   // update color
                   if (p.id()==21){
                       if (std::rand()%2==0){
                           q_acol = p.acol();
                           p.acol(event.nextColTag());
                           q_col = p.acol();
                       } else {
                           q_col = p.col();
                           p.col(event.nextColTag());
                           q_acol = p.col();
                       }
                   }
                   else if(p.id()>0){
                       q_col = p.col();
                       p.col(event.nextColTag());
                       q_acol = p.col();
                   }
                   else {
                       q_acol = p.acol();
                       p.acol(event.nextColTag());
                       q_col = p.acol();
                   }
                   int qid, diqid;
                   double mq, mdiq, mn;
                   if (dist(gen) < ZoverA) { // diquark from a proton
                       diqid = 2101; mdiq = 0.57933;
                       qid = 2; mq = 0.33;
                       mn = 0.93847;
                       if (dist(gen) < 2./3.) { // take away a u
                           qid = 2;
                           if (dist(gen) < .75) diqid = 2101;
                           else diqid = 2103;
                       } else { // take away the d
                           qid = 1;
                           diqid = 2203;
                       }
                  }
                  else { // diquark from a neutron
                       diqid = 2101; mdiq = 0.57933;
                       qid = 1; mq = 0.33;
                       mn = 0.93957;
                       if (dist(gen) < 2./3.) { // take away a d
                           qid = 1;
                           if (dist(gen) < .75) diqid = 2101;
                           else diqid = 2103;
                       } else { // take away the u
                           qid = 2;
                           diqid = 1103;
                       }
                  }
                  double pabs = std::sqrt((mn*mn-std::pow(mq+mdiq,2))*(mn*mn-std::pow(mq-mdiq,2))) / (2.*mn);
                  double costheta = dist(gen)*2.-1.;
                  double sintheta = std::sqrt(std::max(1.-costheta*costheta,1e-9));
                  double rphi = 2*M_PI*dist(gen);
                  double Nqz = pabs*costheta,
                         Nqx = pabs*sintheta*std::cos(rphi),
                         Nqy = pabs*sintheta*std::sin(rphi);
                  Vec4 pq  { Nqx,  Nqy,  Nqz, 0},
                       pdiq{-Nqx, -Nqy,  -Nqz, 0};
                  // decide which object takes the recoil
                  if (std::rand()%2==0){
                      pq = pq - qmu;
                  } else {
                      pdiq = pdiq - qmu;
                  }
                  pq.e(std::sqrt(pq.pAbs2()+mq*mq));
                  pdiq.e(std::sqrt(pdiq.pAbs2()+mdiq*mdiq));
                  Particle recolQ = Particle(qid, 201, i, 0, 0, 0,
                                    q_col, 0,
                                    pq, mq, 0);
                  Particle recoldiQ = Particle(diqid, 201, i, 0, 0, 0,
                                    0, q_acol,
                                    pdiq, mdiq, 0);

                  recoil_remnants.push_back(recolQ);
                  recoil_remnants.push_back(recoldiQ);
           }
           for (auto & p : frag_gluons) new_particles.push_back(p);
           for (auto & p : recoil_remnants) new_particles.push_back(p);
        }

        
        for (auto & p : new_particles)
            event.append(p.id(), 201, p.col(), p.acol(),
                         p.px(), p.py(), p.pz(), p.e(), p.m());
}