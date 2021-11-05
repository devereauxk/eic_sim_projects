const int chargebin = 3; // 0: -, 1: +, 2:+/-

// const int etabin = 3;
// static double eta_lo[etabin] = {-3,-1,1};
// static double eta_hi[etabin] = {-1,1,3};

const int etabin = 1;
static double eta_lo[etabin] = {-3.5};
static double eta_hi[etabin] = {3.5};

// parent pt bin
const int pptbin = 5;
static double ppt_lo[pptbin] = {0.2, 0.5, 1, 2, 4};
static double ppt_hi[pptbin] = {0.5, 1, 2, 4, 10};

//===============================================
//            for events with Q2>1
//===============================================
// const int Q2bin = 6; // last bin inclusive
// static double Q2_lo[Q2bin] = {1E0, 5E0, 1E1, 5E1, 1E2, 1E0};
// static double Q2_hi[Q2bin] = {5E0, 1E1, 5E1, 1E2, 5E2, 5E2};

// const int xbin = 5; // last bin inclusive
// static double x_lo[xbin] = {1E-4, 1E-3, 1E-2, 1E-1, 1E-4};
// static double x_hi[xbin] = {1E-3, 1E-2, 1E-1, 1E0, 1E0};

//===============================================
//            for events with Q2>10
//===============================================
const int Q2bin = 4; // last bin inclusive
static double Q2_lo[Q2bin] = {1E1, 5E1, 1E2, 1E0};
static double Q2_hi[Q2bin] = {5E1, 1E2, 5E2, 5E2};

const int xbin = 4; // last bin inclusive
static double x_lo[xbin] = {1E-3, 1E-2, 1E-1, 1E-4};
static double x_hi[xbin] = {1E-2, 1E-1, 1E0, 1E0};

const char* system_name[2] = {"Pythia, e+p @ 10+100 GeV", "BeAGLE, e+Au @ 18+110 GeV"};
const char* system_abbr[2] = {"ep", "eAu"};