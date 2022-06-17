const int chargebin = 3; // 0: -, 1: +, 2:+/-

const int etabin = 4;
static double eta_lo[etabin] = {-10,-1.5,1.5,-10};
static double eta_hi[etabin] = {-1.5,1.5,10,10};
const int eta_color[etabin] = {kGreen+1, kBlue, kOrange+1, kRed};

// parent pt bin
const int pptbin = 5;
static double ppt_lo[pptbin] = {0.2, 0.5, 1, 2, 4};
static double ppt_hi[pptbin] = {0.5, 1, 2, 4, 10};

//===============================================
//            for events with Q2>1
//===============================================
const int Q2bin = 4; // last bin inclusive
static double Q2_lo[Q2bin] = {1E0, 5E0, 1E1, 1E0};
static double Q2_hi[Q2bin] = {5E0, 1E1, 5E1, 5E1};
const int Q2_color[Q2bin] = {kGreen+1, kBlue, kOrange+1, kRed};

const int xbin = 4; // last bin inclusive
//static double x_lo[xbin] = {1E-3, 1E-2, 1E-1, 1E-3};
//static double x_hi[xbin] = {1E-2, 1E-1, 1E0, 1E0};
static double x_lo[xbin] = {1E-3, 1E-2, 1E-1, 0.023};
static double x_hi[xbin] = {1E-2, 1E-1, 1E0, 0.8};
const int x_color[xbin] = {kGreen+1, kBlue, kOrange+1, kRed};

const int nubin = 4; // last bin inclusive
//static double nu_lo[nubin] = {0, 12, 17, 0};
//static double nu_hi[nubin] = {12, 17, 23.5, 25};
static double nu_lo[nubin] = {0, 12, 17, 6};
static double nu_hi[nubin] = {12, 17, 23.5, 25};
const int nu_color[nubin] = {kGreen+1, kBlue, kOrange+1, kRed};
