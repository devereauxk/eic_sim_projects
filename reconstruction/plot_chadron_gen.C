#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;
unsigned int verbosity = 0;

const int sys_bins = 10;
const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)"};
const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia"};

const int energy_bins = 6;
const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27x0 GeV"};
const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27_0"};

static int cno = 0;

const int MaxVarbin = 10;

class PlotHadron
{
  private:
    // evt info
    int hadron_id;
    const char* hadron_latex;
    const char* hadron_abbr;

    int sys_ep_option;
    int energy_ep_option;
    int sys_eA_option;
    int energy_eA_option;

    int var_option; // as function of x or nu
    int varbin;
    const char* var_latex;
    const char* var_abbr;
    double var_lo[MaxVarbin];
    double var_hi[MaxVarbin];

    TH1D* h1d_nevt_ep[Q2bin][MaxVarbin];
    TH1D* h1d_nevt_w_charm_ep[Q2bin][MaxVarbin];

    TH1D* h1d_nevt_eA[Q2bin][MaxVarbin];
    TH1D* h1d_nevt_w_charm_eA[Q2bin][MaxVarbin];

    // event diagnostic
    TGraphErrors* g_incl_eA_over_ep[Q2bin];
    TGraphErrors* g_charm_eA_over_ep[Q2bin];

    TH2D* h2d_hadron_pt_vs_eta_gen_ep[Q2bin][MaxVarbin];
    TH2D* h2d_hadron_z_vs_eta_gen_ep[Q2bin][MaxVarbin];

    TH2D* h2d_hadron_pt_vs_eta_gen_eA[Q2bin][MaxVarbin];
    TH2D* h2d_hadron_z_vs_eta_gen_eA[Q2bin][MaxVarbin];

    TH1D* h1d_hadron_z_in_eta_gen_ep[Q2bin][MaxVarbin][etabin];
    TH1D* h1d_hadron_z_in_eta_gen_eA[Q2bin][MaxVarbin][etabin];
    TH1D* h1d_hadron_z_in_eta_gen_ratio[Q2bin][MaxVarbin][etabin];

  public:
    PlotHadron(int _hadron_id, int _sys_ep_option, int _energy_ep_option, int _sys_eA_option, int _energy_eA_option, int _var_option)
    {
      cout << "Constructing plotting module to study hadronization of particle with ID " << _hadron_id << endl;
      hadron_id = abs(_hadron_id);
      if (hadron_id==211)
      {
        hadron_latex = "#pi^{#pm}";
        hadron_abbr = "pion";
      }
      else if (hadron_id==111)
      {
        hadron_latex = "#pi^{0}";
        hadron_abbr = "pi0";
      }
      else if (hadron_id==321)
      {
        hadron_latex = "K^{#pm}";
        hadron_abbr = "kaon";
      }
      else if (hadron_id==2212)
      {
        hadron_latex = "p";
        hadron_abbr = "proton";
      }
      else if (hadron_id==421)
      {
        hadron_latex = "D^{0}";
        hadron_abbr = "D0";
      }
      else if (hadron_id==4122)
      {
        hadron_latex = "#Lambda_{c}";
        hadron_abbr = "Lc";
      }
      else
      {
        hadron_latex = "";
        hadron_abbr = "";
      }

      sys_ep_option = _sys_ep_option;
      energy_ep_option = _energy_ep_option;
      cout << "Equivalent \"e+p\" system (demoninator): " << sys_name[sys_ep_option] << " @ " << energy_name[energy_ep_option] << endl;

      sys_eA_option = _sys_eA_option;
      energy_eA_option = _energy_eA_option;
      cout << "Equivalent \"e+A\" system (numerator): " << sys_name[sys_eA_option] << " @ " << energy_name[energy_eA_option] << endl;

      var_option = _var_option;
      cout << "Making plot in Q2 and ";
      if (var_option==0)
      {
        varbin = xbin;
        var_latex = "x";
        var_abbr = "x";
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          var_lo[ivar] = x_lo[ivar];
          var_hi[ivar] = x_hi[ivar];
        }
      }
      else
      {
        varbin = nubin;
        var_latex = "#nu";
        var_abbr = "nu";
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          var_lo[ivar] = nu_lo[ivar];
          var_hi[ivar] = nu_hi[ivar];
        }
      }

      cout << var_abbr << " bins" << endl;
    }

    virtual ~PlotHadron() { };

    void ReadEvtHists(TFile* fin_ep, TFile* fin_eA)
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          // e+p
          h1d_nevt_ep[iQ2][ivar] = (TH1D*)fin_ep->Get(Form("h1d_nevt_in_%s_%d_%d",var_abbr,iQ2,ivar));
          h1d_nevt_ep[iQ2][ivar]->SetName(Form("h1d_nevt_ep_Q2%s_%s_%s%d%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],var_abbr,iQ2,ivar));

          h1d_nevt_w_charm_ep[iQ2][ivar] = (TH1D*)fin_ep->Get(Form("h1d_nevt_w_charm_in_%s_%d_%d",var_abbr,iQ2,ivar));
          h1d_nevt_w_charm_ep[iQ2][ivar]->SetName(Form("h1d_nevt_w_charm_ep_Q2%s_%s_%s%d%d",sys_abbr[sys_ep_option],energy_abbr[energy_ep_option],var_abbr,iQ2,ivar));

          cout << "ep # of inclusive events vs charm events is " << h1d_nevt_ep[iQ2][ivar]->Integral() << " vs " << h1d_nevt_w_charm_ep[iQ2][ivar]->Integral() << " in (Q2," << var_abbr << ") bin (" << iQ2 << ", " << ivar << ")" <<endl;

          // e+A
          h1d_nevt_eA[iQ2][ivar] = (TH1D*)fin_eA->Get(Form("h1d_nevt_in_%s_%d_%d",var_abbr,iQ2,ivar));
          h1d_nevt_eA[iQ2][ivar]->SetName(Form("h1d_nevt_eA_Q2%s_%s_%s%d%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],var_abbr,iQ2,ivar));

          h1d_nevt_w_charm_eA[iQ2][ivar] = (TH1D*)fin_eA->Get(Form("h1d_nevt_w_charm_in_%s_%d_%d",var_abbr,iQ2,ivar));
          h1d_nevt_w_charm_eA[iQ2][ivar]->SetName(Form("h1d_nevt_w_charm_eA_Q2%s_%s_%s%d%d",sys_abbr[sys_eA_option],energy_abbr[energy_eA_option],var_abbr,iQ2,ivar));

          cout << "eA # of inclusive events vs charm events is " << h1d_nevt_eA[iQ2][ivar]->Integral() << " vs " << h1d_nevt_w_charm_eA[iQ2][ivar]->Integral() << " in (Q2," << var_abbr << ") bin (" << iQ2 << ", " << ivar << ")" <<endl;
        }
      }
    }

    void ReadHadronHists(TFile* fin_ep, TFile* fin_eA)
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          // e+p
          h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar] = (TH2D*)fin_ep->Get(Form("h2d_hadron_%d_z_vs_eta_gen_in_%s_%d_%d",hadron_id,var_abbr,iQ2,ivar));
          h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->SetName(Form("h2d_hadron_%d_z_vs_eta_gen_ep_Q2%d_%s%d",hadron_id,iQ2,var_abbr,ivar));

          h2d_hadron_pt_vs_eta_gen_ep[iQ2][ivar] = (TH2D*)fin_ep->Get(Form("h2d_hadron_%d_pt_vs_eta_gen_in_%s_%d_%d",hadron_id,var_abbr,iQ2,ivar));
          h2d_hadron_pt_vs_eta_gen_ep[iQ2][ivar]->SetName(Form("h2d_hadron_%d_pt_vs_eta_gen_ep_Q2%d_%s%d",hadron_id,iQ2,var_abbr,ivar));

          // e+A
          h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar] = (TH2D*)fin_eA->Get(Form("h2d_hadron_%d_z_vs_eta_gen_in_%s_%d_%d",hadron_id,var_abbr,iQ2,ivar));
          h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar]->SetName(Form("h2d_hadron_%d_z_vs_eta_gen_eA_Q2%d_%s%d",hadron_id,iQ2,var_abbr,ivar));

          h2d_hadron_pt_vs_eta_gen_eA[iQ2][ivar] = (TH2D*)fin_eA->Get(Form("h2d_hadron_%d_pt_vs_eta_gen_in_%s_%d_%d",hadron_id,var_abbr,iQ2,ivar));
          h2d_hadron_pt_vs_eta_gen_eA[iQ2][ivar]->SetName(Form("h2d_hadron_%d_pt_vs_eta_gen_eA_Q2%d_%s%d",hadron_id,iQ2,var_abbr,ivar));
        }
      }
    }

    void SetNorm(const int option = 0)
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          float norm_ep = 1, norm_eA = 1;
          if (option==0)
          {
            if (iQ2==0 && ivar==0) cout << "Normlaize using inclusive DIS events" << endl;
            if (h1d_nevt_ep[iQ2][ivar]->Integral()>0 && h1d_nevt_eA[iQ2][ivar]->Integral()>0)
            {
              norm_ep = 1.0/h1d_nevt_ep[iQ2][ivar]->Integral();
              norm_eA = 1.0/h1d_nevt_eA[iQ2][ivar]->Integral();
            }
            else
            {
              cout << "Cannot obtain normalization factor because no event found in this (Q2, "<<var_abbr<<") bin (" << iQ2 << ", " << ivar << ")" << endl;
            }
          }
          else
          {
            if (iQ2==0 && ivar==0) cout << "Normlaize using events containing charm" << endl;
            if (h1d_nevt_w_charm_ep[iQ2][ivar]->Integral()>0 && h1d_nevt_w_charm_eA[iQ2][ivar]->Integral()>0)
            {
              norm_ep = 1.0/h1d_nevt_w_charm_ep[iQ2][ivar]->Integral();
              norm_eA = 1.0/h1d_nevt_w_charm_eA[iQ2][ivar]->Integral();
            }
            else
            {
              cout << "Cannot obtain normalization factor because no event found in this (Q2, "<<var_abbr<<") bin (" << iQ2 << ", " << ivar << ")" << endl;
            }
          }

          // e+p
          h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->Scale(norm_ep);
          h2d_hadron_pt_vs_eta_gen_ep[iQ2][ivar]->Scale(norm_ep);

          // e+A
          h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar]->Scale(norm_eA);
          h2d_hadron_pt_vs_eta_gen_eA[iQ2][ivar]->Scale(norm_eA);
        }
      }
    }

    void PlotEvents()
    {
      double var_mid[MaxVarbin] = {0}; // last varbin is inclusive bin, skip it here
      double var_mid_err[MaxVarbin] = {0};
      for (int ivar = 0; ivar < varbin-1; ++ivar)
      {
        if (var_option==0)
        {
          double logx_lo = log10(x_lo[ivar]);
          double logx_hi = log10(x_hi[ivar]);
          double logx_mid = 0.5*(logx_lo+logx_hi);
          double logx_mid_err = 0.5*(logx_hi-logx_lo);

          var_mid[ivar] = logx_mid;
          var_mid_err[ivar] = logx_mid_err;
        }
        else
        {
          var_mid[ivar] = 0.5*(nu_lo[ivar]+nu_hi[ivar]);
          var_mid_err[ivar] = 0.5*(nu_hi[ivar]-nu_lo[ivar]);
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        double incl_ratio[MaxVarbin] = {0};
        double incl_ratio_err[MaxVarbin] = {0};
        double charm_ratio[MaxVarbin] = {0};
        double charm_ratio_err[MaxVarbin] = {0};
        for (int ivar = 0; ivar < varbin-1; ++ivar)
        {
          incl_ratio[ivar] = (double)h1d_nevt_eA[iQ2][ivar]->GetEntries()/h1d_nevt_ep[iQ2][ivar]->GetEntries();
          incl_ratio_err[ivar] = incl_ratio[ivar]*sqrt(1.0/h1d_nevt_eA[iQ2][ivar]->GetEntries()+1.0/h1d_nevt_ep[iQ2][ivar]->GetEntries());
          if (h1d_nevt_ep[iQ2][ivar]->GetEntries()<=0)
          {
            incl_ratio[ivar] = -9999;
            incl_ratio_err[ivar] = 0;
          }

          charm_ratio[ivar] = (double)h1d_nevt_w_charm_eA[iQ2][ivar]->GetEntries()/h1d_nevt_w_charm_ep[iQ2][ivar]->GetEntries();
          charm_ratio_err[ivar] = charm_ratio[ivar]*sqrt(1.0/h1d_nevt_w_charm_eA[iQ2][ivar]->GetEntries()+1.0/h1d_nevt_w_charm_ep[iQ2][ivar]->GetEntries());
          if (h1d_nevt_w_charm_ep[iQ2][ivar]->GetEntries()<=0)
          {
            charm_ratio[ivar] = -9999;
            charm_ratio_err[ivar] = 0;
          }

          cout << var_mid[ivar] << " inclusive ratio " << incl_ratio[ivar] << " charm ratio " << charm_ratio[ivar] << endl;
        }

        g_incl_eA_over_ep[iQ2] = new TGraphErrors(varbin-1,var_mid,incl_ratio,var_mid_err,incl_ratio_err);
        g_charm_eA_over_ep[iQ2] = new TGraphErrors(varbin-1,var_mid,charm_ratio,var_mid_err,charm_ratio_err);

        mcs(cno++);
        {
          float plot_xrange_lo = var_mid[0]-var_mid_err[0];
          float plot_xrange_hi = var_mid[varbin-2]+var_mid_err[varbin-2];

          float plot_yrange_diff = fabs(TMath::MaxElement(varbin-1,incl_ratio)-TMath::MaxElement(varbin-1,charm_ratio));
          float plot_yrange_lo = TMath::MaxElement(varbin-1,incl_ratio)-3*plot_yrange_diff;
          float plot_yrange_hi = TMath::MaxElement(varbin-1,incl_ratio)+3*plot_yrange_diff;

          TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
          htemp->Draw();
          if (var_option==0) htemp->GetXaxis()->SetTitle("log(x)");
          else htemp->GetXaxis()->SetTitle("#nu");
          htemp->GetYaxis()->SetTitle(Form("N_{evt}^{%s @ %s}/N_{evt}^{%s @ %s}",sys_name[sys_eA_option],energy_name[energy_eA_option],sys_name[sys_ep_option],energy_name[energy_ep_option]));
          myhset(htemp,1.2,1.6,0.05,0.045);

          TLegend* leg = new TLegend(0.2,0.72,0.83,0.82);
          leg->SetBorderSize(0);
          leg->SetTextSize(0.03);
          leg->SetFillStyle(0);
          leg->SetMargin(0.1);

          g_incl_eA_over_ep[iQ2]->SetMarkerSize(0.7);
          g_incl_eA_over_ep[iQ2]->SetMarkerStyle(21);
          g_incl_eA_over_ep[iQ2]->Draw("psame");

          g_charm_eA_over_ep[iQ2]->SetMarkerSize(0.7);
          g_charm_eA_over_ep[iQ2]->SetMarkerStyle(20);
          g_charm_eA_over_ep[iQ2]->SetMarkerColor(kRed);
          g_charm_eA_over_ep[iQ2]->SetLineColor(kRed);
          g_charm_eA_over_ep[iQ2]->Draw("psame");

          leg->AddEntry(g_incl_eA_over_ep[iQ2],"inclusive events","p");
          leg->AddEntry(g_charm_eA_over_ep[iQ2],"charm events","p");

          leg->Draw("same");

          TLatex* tl = new TLatex();
          tl->SetTextAlign(11);
          tl->SetTextSize(0.03);
          tl->DrawLatexNDC(0.21,0.85,Form("%.0f < Q^{2} < %.0f GeV^{2}, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2]));

          gROOT->ProcessLine( Form("cc%d->Print(\"figs/event_ratio_%d.pdf\")", cno-1, iQ2) );

          delete htemp;
          delete leg;
          delete tl;
        }
      }
    }

    void Plot2D()
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          mclogz(cno++);
          {
            float plot_xrange_lo = 0;
            float plot_xrange_hi = 1;

            float plot_yrange_lo = eta_lo[0];
            float plot_yrange_hi = eta_hi[etabin-1];

            TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
            htemp->GetYaxis()->SetTitle(Form("#eta^{%s}",hadron_latex));
            myhset(htemp,1.2,1.6,0.05,0.045);

            h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->Draw("samecolz");

            TLatex* tl = new TLatex();
            tl->SetTextAlign(11);
            tl->SetTextSize(0.03);
            tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
            tl->DrawLatexNDC(0.21,0.80,Form("%.0f < Q^{2} < %.0f GeV^{2}, %.0e < %s < %.0e, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2],var_lo[ivar],var_latex,var_hi[ivar]));

            gROOT->ProcessLine( Form("cc%d->Print(\"figs/ep_%s_z_vs_eta_%d_%d.pdf\")", cno-1, hadron_abbr, iQ2, ivar) );

            delete htemp;
            delete tl;
          }

          mclogz(cno++);
          {
            float plot_xrange_lo = 0;
            float plot_xrange_hi = 10;

            float plot_yrange_lo = eta_lo[0];
            float plot_yrange_hi = eta_hi[etabin-1];

            TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle(Form("p_{T}^{%s}",hadron_latex));
            htemp->GetYaxis()->SetTitle(Form("#eta^{%s}",hadron_latex));
            myhset(htemp,1.2,1.6,0.05,0.045);

            h2d_hadron_pt_vs_eta_gen_ep[iQ2][ivar]->Draw("samecolz");

            TLatex* tl = new TLatex();
            tl->SetTextAlign(11);
            tl->SetTextSize(0.03);
            tl->DrawLatexNDC(0.21,0.85,Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]));
            tl->DrawLatexNDC(0.21,0.80,Form("%.0f < Q^{2} < %.0f GeV^{2}, %.0e < %s < %.0e, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2],var_lo[ivar],var_latex,var_hi[ivar]));

            gROOT->ProcessLine( Form("cc%d->Print(\"figs/ep_%s_pt_vs_eta_%d_%d.pdf\")", cno-1, hadron_abbr, iQ2, ivar) );

            delete htemp;
            delete tl;
          }
        }
      }
    }

    void SliceInEta()
    { // make sure slice_lo/hi are double before slicing
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          for (int ieta = 0; ieta < etabin; ++ieta)
          {
            h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
            h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta] = (TH1D*)h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->ProjectionX( Form("h1d_hadron_%d_z_in_eta_gen_ep_Q2%d_x%d_eta%d",hadron_id,iQ2,ivar,ieta) );

            h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar]->GetYaxis()->SetRangeUser(eta_lo[ieta], eta_hi[ieta]);
            h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta] = (TH1D*)h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar]->ProjectionX( Form("h1d_hadron_%d_z_in_eta_gen_eA_Q2%d_x%d_eta%d",hadron_id,iQ2,ivar,ieta) );

            h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->Rebin();
            h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->Rebin();

            h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->Rebin();
            h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->Rebin();
          }

          // set back to the original range
          h2d_hadron_z_vs_eta_gen_ep[iQ2][ivar]->GetYaxis()->SetRangeUser(eta_lo[0],eta_hi[etabin-1]);
          h2d_hadron_z_vs_eta_gen_eA[iQ2][ivar]->GetYaxis()->SetRangeUser(eta_lo[0],eta_hi[etabin-1]);

          // set ratio hists
          for (int ieta = 0; ieta < etabin; ++ieta)
          {
            h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta] = (TH1D*)h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->Clone();
            h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetName( Form("h1d_hadron_%d_z_in_eta_gen_ratio_Q2%d_x%d_eta%d",hadron_id,iQ2,ivar,ieta) );
            h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Divide(h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]);
          }
        }
      }
    }

    void Plot1D()
    {
      TGaxis::SetMaxDigits(3);

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        if (iQ2<Q2bin-1) continue;

        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          if (ivar<varbin-1) continue;

          for (int ieta = 0; ieta < etabin; ++ieta)
          {
            if (ieta<etabin-1) continue;

            mcs(cno++);
            {
              float plot_xrange_lo = 0;
              float plot_xrange_hi = 1;

              float plot_yrange_lo = 0.5*h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->GetMinimum();
              float plot_yrange_hi = 1.8*h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->GetMaximum();

              TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
              htemp->Draw();
              htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
              htemp->GetYaxis()->SetTitle(Form("Normalized N^{%s}",hadron_latex));
              myhset(htemp,1.2,1.6,0.05,0.045);

              TLegend* leg = new TLegend(0.2,0.69,0.83,0.78);
              leg->SetBorderSize(0);
              leg->SetTextSize(0.03);
              leg->SetFillStyle(0);
              leg->SetMargin(0.1);

              h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->SetMarkerSize(0.7);
              h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->SetLineColor(kBlack);
              h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->SetMarkerColor(kBlack);
              h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->Draw("same");
              leg->AddEntry(h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta],Form("%s @ %s",sys_name[sys_ep_option],energy_name[energy_ep_option]),"lp");

              h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->SetMarkerSize(0.7);
              h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->SetLineColor(kRed);
              h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->SetMarkerColor(kRed);
              h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->Draw("same");
              leg->AddEntry(h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta],Form("%s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option]),"lp");

              leg->Draw("same");

              TLatex* tl = new TLatex();
              tl->SetTextAlign(11);
              tl->SetTextSize(0.03);
              tl->DrawLatexNDC(0.21,0.85,Form("%.0f < Q^{2} < %.0f GeV^{2}, %.0e < %s < %.0e, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2],var_lo[ivar],var_latex,var_hi[ivar]));
              tl->DrawLatexNDC(0.21,0.80,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));

              gROOT->ProcessLine( Form("cc%d->Print(\"figs/%s_z_in_eta_%d_%d_%d.pdf\")", cno-1, hadron_abbr, iQ2, ivar, ieta) );

              delete htemp;
              delete leg;
              delete tl;
            }

            mcs(cno++);
            {
              float plot_xrange_lo = 0;
              float plot_xrange_hi = 1;

              float plot_yrange_lo = 0.;
              float plot_yrange_hi = 2;

              TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
              htemp->Draw();
              htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
              htemp->GetYaxis()->SetTitle(Form("R^{%s}_{eA}",hadron_latex));
              myhset(htemp,1.2,1.6,0.05,0.05);

              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerSize(0.7);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetLineColor(kRed);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerColor(kRed);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Draw("same");

              TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
              l1.SetLineStyle(7);
              l1.SetLineColor(kGray+2);
              l1.Draw("same");

              TLatex* tl = new TLatex();
              tl->SetTextAlign(11);
              tl->SetTextSize(0.03);
              tl->DrawLatexNDC(0.21,0.85,Form("%.0f < Q^{2} < %.0f GeV^{2}, %.0e < %s < %.0e, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2],var_lo[ivar],var_latex,var_hi[ivar]));
              tl->DrawLatexNDC(0.21,0.80,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
              tl->DrawLatexNDC(0.21,0.75,Form("%s @ %s / %s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option],sys_name[sys_ep_option],energy_name[energy_ep_option]));

              gROOT->ProcessLine( Form("cc%d->Print(\"figs/%s_ratio_z_in_eta_%d_%d_%d.pdf\")", cno-1, hadron_abbr, iQ2, ivar, ieta) );

              delete htemp;
              delete tl;
            }
          }
        }
      }
    }

    void PlotComparison()
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        if (iQ2<Q2bin-1) continue;

        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          if (ieta<etabin-1) continue;

          mcs(cno++);
          {
            float plot_xrange_lo = 0;
            float plot_xrange_hi = 1;

            float plot_yrange_lo = 0.;
            float plot_yrange_hi = 2;

            TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
            htemp->GetYaxis()->SetTitle(Form("R^{%s}_{eA}",hadron_latex));
            myhset(htemp,1.2,1.6,0.05,0.05);

            TLegend* leg = new TLegend(0.21,0.17,0.51,0.29);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->SetFillStyle(0);
            leg->SetMargin(0.1);

            for (int ivar = 0; ivar < varbin-1; ++ivar)
            { // excluding inclusive bin
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerColor(x_color[ivar]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetLineColor(x_color[ivar]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerSize(0.5);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Draw("same");
              leg->AddEntry(h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta],Form("%.1e < %s < %.1e",var_lo[ivar],var_latex,var_hi[ivar]),"p");
            }
            leg->Draw("same");

            TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
            l1.SetLineStyle(7);
            l1.SetLineColor(kGray+2);
            l1.Draw("same");

            TLatex* tl = new TLatex();
            tl->SetTextAlign(11);
            tl->SetTextSize(0.03);
            tl->DrawLatexNDC(0.21,0.85,Form("%.0f < Q^{2} < %.0f GeV^{2}, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2]));
            tl->DrawLatexNDC(0.21,0.80,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
            tl->DrawLatexNDC(0.21,0.75,Form("%s @ %s / %s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option],sys_name[sys_ep_option],energy_name[energy_ep_option]));

            gROOT->ProcessLine( Form("cc%d->Print(\"figs/%s_ratio_z_in_eta_diff_%s_%d_%d.pdf\")", cno-1, hadron_abbr, var_abbr, iQ2, ieta) );

            delete htemp;
            delete leg;
            delete tl;
          }
        }
      }

      for (int ivar = 0; ivar < varbin; ++ivar)
      {
        if (ivar<varbin-1) continue;

        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          if (ieta<etabin-1) continue;

          mcs(cno++);
          {
            float plot_xrange_lo = 0;
            float plot_xrange_hi = 1;

            float plot_yrange_lo = 0.;
            float plot_yrange_hi = 2;

            TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
            htemp->GetYaxis()->SetTitle(Form("R^{%s}_{eA}",hadron_latex));
            myhset(htemp,1.2,1.6,0.05,0.05);

            TLegend* leg = new TLegend(0.21,0.17,0.51,0.29);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->SetFillStyle(0);
            leg->SetMargin(0.1);

            for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
            { // excluding inclusive bin
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerColor(Q2_color[iQ2]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetLineColor(Q2_color[iQ2]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerSize(0.5);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Draw("same");
              leg->AddEntry(h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta],Form("%.1e < Q^{2} < %.1e",Q2_lo[iQ2],Q2_hi[iQ2]),"p");
            }
            leg->Draw("same");

            TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
            l1.SetLineStyle(7);
            l1.SetLineColor(kGray+2);
            l1.Draw("same");

            TLatex* tl = new TLatex();
            tl->SetTextAlign(11);
            tl->SetTextSize(0.03);
            tl->DrawLatexNDC(0.21,0.85,Form("%.0e < %s < %.0e, 0.05 < y < 0.8",var_lo[ivar],var_latex,var_hi[ivar]));
            tl->DrawLatexNDC(0.21,0.80,Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]));
            tl->DrawLatexNDC(0.21,0.75,Form("%s @ %s / %s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option],sys_name[sys_ep_option],energy_name[energy_ep_option]));

            gROOT->ProcessLine( Form("cc%d->Print(\"figs/%s_ratio_z_in_eta_diff_Q2_%d_%d.pdf\")", cno-1, hadron_abbr, ivar, ieta) );

            delete htemp;
            delete leg;
            delete tl;
          }
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        if (iQ2<Q2bin-1) continue;

        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          if (ivar<varbin-1) continue;

          mcs(cno++);
          {
            float plot_xrange_lo = 0;
            float plot_xrange_hi = 1;

            float plot_yrange_lo = 0.;
            float plot_yrange_hi = 2;

            TH2F* htemp = new TH2F("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
            htemp->Draw();
            htemp->GetXaxis()->SetTitle(Form("z^{%s}",hadron_latex));
            htemp->GetYaxis()->SetTitle(Form("R^{%s}_{eA}",hadron_latex));
            myhset(htemp,1.2,1.6,0.05,0.05);

            TLegend* leg = new TLegend(0.21,0.17,0.51,0.29);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);
            leg->SetFillStyle(0);
            leg->SetMargin(0.1);

            for (int ieta = 0; ieta < etabin-1; ++ieta)
            { // excluding inclusive bin
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerColor(eta_color[ieta]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetLineColor(eta_color[ieta]);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerSize(0.5);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->SetMarkerStyle(21);
              h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Draw("same");
              leg->AddEntry(h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta],Form("%.1f < #eta < %.1f",eta_lo[ieta],eta_hi[ieta]),"p");
            }
            leg->Draw("same");

            TLine l1(plot_xrange_lo,1,plot_xrange_hi,1);
            l1.SetLineStyle(7);
            l1.SetLineColor(kGray+2);
            l1.Draw("same");

            TLatex* tl = new TLatex();
            tl->SetTextAlign(11);
            tl->SetTextSize(0.03);
            tl->DrawLatexNDC(0.21,0.85,Form("%.0f < Q^{2} < %.0f GeV^{2}, %.0e < %s < %.0e, 0.05 < y < 0.8",Q2_lo[iQ2],Q2_hi[iQ2],var_lo[ivar],var_latex,var_hi[ivar]));
            tl->DrawLatexNDC(0.21,0.80,Form("%s @ %s / %s @ %s",sys_name[sys_eA_option],energy_name[energy_eA_option],sys_name[sys_ep_option],energy_name[energy_ep_option]));

            gROOT->ProcessLine( Form("cc%d->Print(\"figs/%s_ratio_z_in_eta_diff_eta_%d.pdf\")", cno-1, hadron_abbr, ivar) );

            delete htemp;
            delete leg;
            delete tl;
          }
        }
      }
    }

    void WriteHadronHists(TFile* fout)
    {
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ivar = 0; ivar < varbin; ++ivar)
        {
          for (int ieta = 0; ieta < etabin; ++ieta)
          {
            h1d_hadron_z_in_eta_gen_ep[iQ2][ivar][ieta]->Write();
            h1d_hadron_z_in_eta_gen_eA[iQ2][ivar][ieta]->Write();
            h1d_hadron_z_in_eta_gen_ratio[iQ2][ivar][ieta]->Write();
          }
        }
      }
    }
};

void plot_chadron_gen(const char* inFile_ep = "hists_gen_ep.root", const int sys_ep_option = 0, const int energy_ep_option = 0, const char* inFile_eA = "hists_gen_eA.root", const int sys_eA_option = 0, const int energy_eA_option = 0, const char* outFile = "hists_gen.root")
{
  mcs(-1);

  TFile* fin_ep = new TFile(inFile_ep,"READ");
  TFile* fin_eA = new TFile(inFile_eA,"READ");

  int var_option = 0; // 0: x, 1: nu
  PlotHadron plot_pion(211,sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option,var_option);
  plot_pion.ReadEvtHists(fin_ep,fin_eA);
  plot_pion.ReadHadronHists(fin_ep,fin_eA);
  plot_pion.SetNorm(0);
  // plot_pion.Plot2D();
  plot_pion.SliceInEta();
  plot_pion.Plot1D();
  plot_pion.PlotComparison();

  PlotHadron plot_kaon(321,sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option,var_option);
  plot_kaon.ReadEvtHists(fin_ep,fin_eA);
  plot_kaon.ReadHadronHists(fin_ep,fin_eA);
  plot_kaon.SetNorm(0);
  // plot_kaon.Plot2D();
  plot_kaon.SliceInEta();
  plot_kaon.Plot1D();
  plot_kaon.PlotComparison();

  // PlotHadron plot_proton(2212,sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option,var_option);
  // plot_proton.ReadEvtHists(fin_ep,fin_eA);
  // plot_proton.ReadHadronHists(fin_ep,fin_eA);
  // plot_proton.SetNorm(0);
  // // plot_proton.Plot2D();
  // plot_proton.SliceInEta();
  // plot_proton.Plot1D();
  // plot_proton.PlotComparison();

  PlotHadron plot_D0(421,sys_ep_option,energy_ep_option,sys_eA_option,energy_eA_option,var_option);
  plot_D0.ReadEvtHists(fin_ep,fin_eA);
  plot_D0.ReadHadronHists(fin_ep,fin_eA);
  plot_D0.PlotEvents();
  plot_D0.SetNorm(1);
  // plot_D0.Plot2D();
  plot_D0.SliceInEta();
  plot_D0.Plot1D();
  plot_D0.PlotComparison();

  TFile* fout = new TFile(outFile,"recreate");
  plot_pion.WriteHadronHists(fout);
  plot_kaon.WriteHadronHists(fout);
  // plot_D0.WriteHadronHists(fout);
  fout->Write();
}
