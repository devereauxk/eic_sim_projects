#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 4;
const char* sys_name[sys_bins] = {"e+C", "e+Cu", "e+Au", "e+Pb"};
const char* sys_abbr[sys_bins] = {"eC", "eCu", "eAu", "ePb"};
const int sys_color[sys_bins] = {kBlack, kRed, kOrange+1, kGreen+1};

const int energy_bins = 6;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

TProfile* prof_thickness_in_b[sys_bins] = {0};
TProfile* prof_Ninc_in_thickness[sys_bins] = {0};
TProfile* prof_Ninc_in_b[sys_bins] = {0};

static int cno = 0;

void plot_comparison()
{

}

void compare_plot_inc_hadron_diff_sys(const int energy_option, const char* outDir = "figs/")
{


}
