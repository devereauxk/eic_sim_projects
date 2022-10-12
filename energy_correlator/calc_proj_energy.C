// calculate the projectile energy needed to use in eHIJING setting to produce certain energy combination in lab fram
// assumes fixed target in eHIJING

R__LOAD_LIBRARY(libeicsmear);

const Double_t Mp(0.9383); // in GeV/c^2
const Double_t Me(0.511E-3); // in GeV/c^2
const Double_t MAu(183.4343); // in GeV/c^2

void calc_proj_energy(double proj_lab_e = 10, double targ_lab_e = 100, double targ_atomic_mass = 197, targ_mass = MAu)
{
  // proj_lab is desired projectile energy in lab frame
  // trag_lab is desired target energy in lab frame
  // positive values
  // must be e+A where A is specified by target atomic mass
  // default is e+Au 10x100 GeV

  TLorentzVector proj_lab, targ_lab, proj_rest, targ_rest;
  proj_lab.SetXYZM(0,0,-proj_lab_e,Me);
  targ_lab.SetXYZM(0,0,targ_lab_e * targ_atomic_mass,targ_mass);

  TVector3 boost_vec = targ_lab.BoostVector();

  proj_rest = proj_lab; proj_rest.Boost(-boost_vec);
  targ_rest = targ_lab; targ_rest.Boost(-boost_vec);

  cout<<"projectile in target rest frame:"<<endl;
  proj_rest.Print();
  cout<<"target in target rest frame:"<<endl;
  targ_rest.Print();

  cout<<"set projectile energy to "<<proj_rest.E()<<" GeV"<<endl;

}
