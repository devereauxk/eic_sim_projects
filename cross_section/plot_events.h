// Constants
const double ebeam = 10;
const double hbeam = 100;
const double alpha(1./137.036);
const double s_cm = 4.*ebeam*hbeam;
const double microbgev(1./(0.3894E3));

// Cross section (microbarn)
const double cs_Pythia6_eAu = 3.64279E-2; // Q2>10
const double cs_BeAGLE_eAu = 0.037279059622903579; // 3.16941E-2; // genShd =3, LO nPDF
const double cs_BeAGLE_eAu_NLO = 0.035348639250616293; // NLO
const double cs_BeAGLE_eAu_genShd1 = 0.032861777825967387; // 4.086525E-2; // genShd = 1, LO nPDF
const double cs_Pythia6_ep = 4.08586E-2; // 10050 setup
const double cs_BeAGLE_ep = 4.06891E-2;
const double cs_Pythia6_eAu_10042 = 3.6427976319158305E-002;
const double cs_Pythia6_eAu_10042_EPSLO = 3.6255883519827781E-002;

// Cross Section Bins
const int Q2bin = 7;
const double Q2_mid[Q2bin] = {12, 20, 35, 60, 120, 200, 300};
double Q2_lo[Q2bin] = {0};
double Q2_hi[Q2bin] = {0};

const int xbin = 15;
double x_lo[xbin] = {0};
double x_hi[xbin] = {0};
double x_mid[xbin] = {0};

class PlotXsec
{
  private:
    // evt info
    double Xsec_gen;
    double L_gen;

  public:
    const char* sys_name;
    const char* sys_latex;
    const char* sys_abbr;

    TH1D* h1d_nevt;
    TH1D* h1d_x_c[Q2bin];
    TH1D* h1d_x_e[Q2bin];

    TGraphErrors* g_cs_vs_x_e[Q2bin];
    TGraphErrors* g_cs_vs_x_c[Q2bin];

    PlotXsec(int _sys_option)
    {
      if (_sys_option==1)
      {
        sys_name = "Pythia e+p";
        sys_latex = "Pythia e+p";
        sys_abbr = "Pythia_ep";

        Xsec_gen = cs_Pythia6_ep;
      }
      else if (_sys_option==2)
      {
        sys_name = "Pythia e+Au";
        sys_latex = "Pythia e+Au";
        sys_abbr = "Pythia_eAu";

        Xsec_gen = cs_Pythia6_eAu;
      }
      else if (_sys_option==3)
      {
        sys_name = "BeAGLE e+p";
        sys_latex = "BeAGLE e+p";
        sys_abbr = "BeAGLE_ep";

        Xsec_gen = cs_BeAGLE_ep;
      }
      else if (_sys_option==4)
      {
        sys_name = "BeAGLE e+Au";
        sys_latex = "BeAGLE e+Au";
        sys_abbr = "BeAGLE_eAu";

        Xsec_gen = cs_BeAGLE_eAu;
      }
      else if (_sys_option==5)
      {
        sys_name = "Pythia e+Au";
        sys_latex = "Pythia e+Au (10042)";
        sys_abbr = "Pythia_eAu_10042";

        Xsec_gen = cs_Pythia6_eAu_10042;
      }
      else if (_sys_option==6)
      {
        sys_name = "BeAGLE e+Au";
        sys_latex = "BeAGLE e+Au (NLO)";
        sys_abbr = "BeAGLE_eAu_NLO";

        Xsec_gen = cs_BeAGLE_eAu_NLO;
      }
      else if (_sys_option==7)
      {
        sys_name = "Pythia e+Au";
        sys_latex = "Pythia e+Au (LO)";
        sys_abbr = "Pythia_eAu_LO";

        Xsec_gen = cs_Pythia6_eAu_10042_EPSLO;
      }
      else if (_sys_option==8)
      {
        sys_name = "BeAGLE e+Au";
        sys_latex = "BeAGLE e+Au (genShd=1)";
        sys_abbr = "BeAGLE_eAu_genShd1";

        Xsec_gen = cs_BeAGLE_eAu_genShd1;
      }
      else
      {
        cout << "Not recognized, abort..." << endl;
        exit(0);
      }

      cout << "Constructing plotting module to calculate reduced cross section for system " << sys_name << endl;
    }

    virtual ~PlotXsec() { };

    void ReadEvtHists(TFile* fin)
    {
      h1d_nevt = (TH1D*)fin->Get("h1d_nevt");

      L_gen = h1d_nevt->GetEntries()/Xsec_gen; 
      
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        h1d_x_e[iQ2] = (TH1D*)fin->Get(Form("h1d_x_e_Q2_%d",iQ2));
        h1d_x_c[iQ2] = (TH1D*)fin->Get(Form("h1d_x_c_Q2_%d",iQ2));
      }
    }

    void CalculateReducedXsec()
    {
      //============================
      // inclusive reduced Xsection
      //============================
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        double rcs[xbin] = {0};
        double rcs_err[xbin] = {0};

        for (int ibin = 1; ibin < h1d_x_e[iQ2]->GetNbinsX()+1; ++ibin)
        { // ibin = 0 is underflow bin, skip it
          double raw_count = h1d_x_e[iQ2]->GetBinContent(ibin);
          double raw_count_err = h1d_x_e[iQ2]->GetBinError(ibin);

          double Q2_bw = Q2_hi[iQ2]-Q2_lo[iQ2];
          raw_count /= h1d_x_e[iQ2]->GetBinWidth(ibin)*Q2_bw; // dN/dxdQ2
          raw_count_err /= h1d_x_e[iQ2]->GetBinWidth(ibin)*Q2_bw; // dN/dxdQ2

          // calculate cross section (convert to GeV-2 unit)
          double cs = raw_count/L_gen*microbgev;
          double cs_err = raw_count_err/L_gen*microbgev;

          // phase space factor
          x_mid[ibin-1] = h1d_x_e[iQ2]->GetBinCenter(ibin);
          x_lo[ibin-1] = h1d_x_e[iQ2]->GetBinLowEdge(ibin);
          x_hi[ibin-1] = h1d_x_e[iQ2]->GetBinLowEdge(ibin+1);
          double y = Q2_mid[iQ2]/(x_mid[ibin-1]*s_cm);
          double ps_factor = x_mid[ibin-1]*pow(Q2_mid[iQ2],2)/(2*TMath::Pi()*pow(alpha,2)*(1+pow(1-y,2)));

          rcs[ibin-1] = cs*ps_factor;
          rcs_err[ibin-1] = cs_err*ps_factor;
        }

        g_cs_vs_x_e[iQ2] = new TGraphErrors(xbin,x_mid,rcs,0,rcs_err);
        g_cs_vs_x_e[iQ2]->SetMarkerStyle(20);
        g_cs_vs_x_e[iQ2]->SetMarkerSize(0.7);

        for (int ix = xbin-1; ix >= 0; --ix)
        {
          double temp_y = Q2_mid[iQ2]/(x_mid[ix]*s_cm);
          if (temp_y<0.001 || temp_y>0.90) g_cs_vs_x_e[iQ2]->RemovePoint(ix);
        }
      }
      
      //============================
      //    charm reduced Xsection
      //============================
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        double rcs[xbin] = {0};
        double rcs_err[xbin] = {0};

        for (int ibin = 1; ibin < h1d_x_c[iQ2]->GetNbinsX()+1; ++ibin)
        { // ibin = 0 is underflow bin, skip it
          double raw_count = h1d_x_c[iQ2]->GetBinContent(ibin);
          double raw_count_err = h1d_x_c[iQ2]->GetBinError(ibin);

          double Q2_bw = Q2_hi[iQ2]-Q2_lo[iQ2];
          raw_count /= h1d_x_c[iQ2]->GetBinWidth(ibin)*Q2_bw; // dN/dxdQ2
          raw_count_err /= h1d_x_c[iQ2]->GetBinWidth(ibin)*Q2_bw; // dN/dxdQ2

          // calculate cross section (convert to GeV-2 unit)
          double cs = raw_count/L_gen*microbgev;
          double cs_err = raw_count_err/L_gen*microbgev;

          // phase space factor
          x_mid[ibin-1] = h1d_x_c[iQ2]->GetBinCenter(ibin);
          x_lo[ibin-1] = h1d_x_c[iQ2]->GetBinLowEdge(ibin);
          x_hi[ibin-1] = h1d_x_c[iQ2]->GetBinLowEdge(ibin+1);
          double y = Q2_mid[iQ2]/(x_mid[ibin-1]*s_cm);
          double ps_factor = x_mid[ibin-1]*pow(Q2_mid[iQ2],2)/(2*TMath::Pi()*pow(alpha,2)*(1+pow(1-y,2)));

          rcs[ibin-1] = cs*ps_factor;
          rcs_err[ibin-1] = cs_err*ps_factor;
        }

        g_cs_vs_x_c[iQ2] = new TGraphErrors(xbin,x_mid,rcs,0,rcs_err);
        g_cs_vs_x_c[iQ2]->SetMarkerStyle(20);
        g_cs_vs_x_c[iQ2]->SetMarkerSize(0.7);

        for (int ix = xbin-1; ix >= 0; --ix)
        {
          double temp_y = Q2_mid[iQ2]/(x_mid[ix]*s_cm);
          if (temp_y<0.001 || temp_y>0.90) g_cs_vs_x_c[iQ2]->RemovePoint(ix);
        }
      }
    }

    void WriteHadronHists(TFile* fout)
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        g_cs_vs_x_e[iQ2]->Write(Form("g_cs_vs_x_e_%s_Q2_%d",sys_abbr,iQ2));
        g_cs_vs_x_c[iQ2]->Write(Form("g_cs_vs_x_c_%s_Q2_%d",sys_abbr,iQ2));
      }
    }
};
