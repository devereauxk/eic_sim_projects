#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "TString.h"

using namespace std;
using namespace Pythia8;

const int verbosity = 0;

int get_seed()
{
  int fSeed;

  ifstream devrandom;
  devrandom.open("/dev/urandom",ios::binary);
  devrandom.read((char*)&fSeed,sizeof(fSeed));
  devrandom.close();
  if ( fSeed != -1 )
  {
    cout << " Got seed from /dev/urandom" << endl;
    fSeed = abs(fSeed)%900000000; // pythia take seed from -1 to 900000000
  }
  else fSeed = 0;
  cout << "seed is " << fSeed << endl;
  return fSeed;
}

int main()
{
  
  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC3::Pythia8ToHepMC3 topHepMC;

  // Specify file where HepMC events will be stored.
  HepMC3::WriterAscii ascii_io("hepmcout.dat");

  // Beam energies, minimal Q2, number of events to generate.
  double eProton   = 275.;
  double eElectron = 18.;
  double Q2min     = 10.0;
  int    nEvent    = 400000;

  // Generator. Shorthand for event.
  Pythia pythia;
  //Event& event = pythia.event;

  // Pick new random number seed for each run, based on dev/urandom
  string seed = to_string(get_seed());
  cout<<"string seed "<<seed<<endl;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"+seed);

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  pythia.readString("Beams:idB = 11");
  pythia.settings.parm("Beams:eB", eElectron);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Uncomment to allow charged current.
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);
  pythia.readString("PhaseSpace:mHatMin = force 0.1");
  pythia.readString("PhaseSpace:pTHatMindiverge = force 0.1");

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // Input PDF set
  pythia.readString("PDF:pSet=8");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Hadronization and radiation settings
  pythia.readString("HadronLevel:Decay = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("PromptPhoton:all = off");

  // turn off charm meson decays
  pythia.readString("411:mayDecay = no");
  pythia.readString("421:mayDecay = no");
  pythia.readString("10411:mayDecay = no");
  pythia.readString("10421:mayDecay = no");
  pythia.readString("413:mayDecay = no");
  pythia.readString("423:mayDecay = no");
  pythia.readString("10413:mayDecay = no");
  pythia.readString("10423:mayDecay = no");
  pythia.readString("20413:mayDecay = no");
  pythia.readString("20423:mayDecay = no");
  pythia.readString("415:mayDecay = no");
  pythia.readString("425:mayDecay = no");
  pythia.readString("431:mayDecay = no");
  pythia.readString("10431:mayDecay = no");
  pythia.readString("433:mayDecay = no");
  pythia.readString("10433:mayDecay = no");
  pythia.readString("20433:mayDecay = no");
  pythia.readString("435:mayDecay = no");

  // turn off charm baryon decay
  pythia.readString("4122:mayDecay = no");
  pythia.readString("4222:mayDecay = no");
  pythia.readString("4212:mayDecay = no");
  pythia.readString("4112:mayDecay = no");
  pythia.readString("4224:mayDecay = no");
  pythia.readString("4214:mayDecay = no");
  pythia.readString("4114:mayDecay = no");
  pythia.readString("4232:mayDecay = no");
  pythia.readString("4132:mayDecay = no");
  pythia.readString("4322:mayDecay = no");
  pythia.readString("4312:mayDecay = no");
  pythia.readString("4324:mayDecay = no");
  pythia.readString("4314:mayDecay = no");
  pythia.readString("4332:mayDecay = no");
  pythia.readString("4334:mayDecay = no");
  pythia.readString("4412:mayDecay = no");
  pythia.readString("4422:mayDecay = no");
  pythia.readString("4414:mayDecay = no");
  pythia.readString("4424:mayDecay = no");
  pythia.readString("4432:mayDecay = no");
  pythia.readString("4434:mayDecay = no");
  pythia.readString("4444:mayDecay = no");

  // Initialize.
  pythia.init();

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next()) continue;

    //if (pythia.event.size()<7) continue; // somehow not even enough entries for a LO DIS process

    if (iEvent<10)
       pythia.event.list();

    // Construct new empty HepMC event and fill it.
    // Default units are ( HepMC3::Units::GEV, HepMC3::Units::MM)
    // but can be changed in the GenEvent constructor.
    HepMC3::GenEvent hepmcevt;
    topHepMC.fill_next_event( pythia, &hepmcevt );

    // Write the HepMC event to file.
    ascii_io.write_event(hepmcevt);

    if (verbosity>1)
      pythia.process.list();

  } //End of event loop.

  //Statistics
  pythia.stat();

  // Done.
  return 0;
}
