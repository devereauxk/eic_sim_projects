// calculate the projectile energy needed to use in eHIJING setting to produce certain energy combination in lab fram
// assumes fixed target in eHIJING

R__LOAD_LIBRARY(libeicsmear);

const Double_t Mp(0.9383);
const Double_t Me(0.511E-3);

void calc_proj_energy(double proj_lab_e = 10, double targ_lab_e = 100, double proj_mass = Me, double targ_mass = Mp)
{
  // proj_lab is desired projectile energy in lab frame
  // trag_lab is desired target energy in lab frame
  // positive values
  // default is 10x100 GeV e+p

  TLorentzVector proj_lab, targ_lab, proj_rest, targ_rest;
  proj_lab.SetXYZM(0,0,-proj_lab_e,proj_mass);
  targ_lab.SetXYZM(0,0,targ_lab_e,targ_mass);

  TVector3 boost_vec = targ_lab.BoostVector();

  proj_rest = proj_lab; proj_rest.Boost(-boost_vec);
  targ_rest = targ_lab; targ_rest.Boost(-boost_vec);

  cout<<"projectile in target rest frame:"<<endl;
  proj_rest.Print();
  cout<<"target in target rest frame:"<<endl;
  targ_rest.Print();

  cout<<"set projectile energy to "<<proj_rest.E()<<" GeV"<<endl;

}
