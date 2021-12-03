R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

static int cno = 0;

// TODO mcs and myhset alias

void plot_1D(const int ismear = 4, const int ibfield = 1, const int ipid = 1, const int isys = 0)
{
  TGaxis::SetMaxDigits(3);
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    for (int iz = 0; iz < zbin; ++iz)
    {
      mcs(cno++);
      {
        float plot_xrange_lo = 1.5;
        float plot_xrange_hi = 2.2;

        fg1d_Kpimass_vs_z[ismear][ibfield][1][ieta][iz]->GetXaxis()->SetRangeUser(plot_xrange_lo,plot_xrange_hi);
        float plot_yrange_lo = 0;
        float plot_yrange_hi = 2.0*fg1d_Kpimass_vs_z[ismear][ibfield][1][ieta][iz]->GetMaximum();

        TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
        htemp.Draw();
        htemp.GetXaxis()->SetTitle("m_{K#pi}");
        htemp.GetYaxis()->SetTitle("Counts");
        myhset(&htemp,1.2,1.6,0.05,0.045);

        TLegend leg(0.2,0.65,0.83,0.80);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.03);
        leg.SetFillStyle(0);
        leg.SetMargin(0.1);

        fg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->SetLineColor(kBlue);
        bg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->SetLineColor(kRed);
        sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->SetLineColor(kGreen+1);
        fg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Draw("hsame");
        bg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Draw("hsame");
        // sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Draw("hsame");

        float temp_mean = -9999;
        float temp_sigma = 0;
        sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Fit("gaus","0R","",1.8,1.95);
        TF1* gaus = sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->GetFunction("gaus");

        if (gaus!=NULL)
        {
          temp_mean = gaus->GetParameter(1);
          temp_sigma = gaus->GetParameter(2);

          TF1* func_peak = new TF1("func_peak","[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*pow(x-[4],2)/pow([5],2))",temp_mean-3*temp_sigma,temp_mean+3*temp_sigma);
          func_peak->SetLineColor(pid_color[ipid]);
          func_peak->FixParameter(3,gaus->GetParameter(0));
          func_peak->FixParameter(4,gaus->GetParameter(1));
          func_peak->FixParameter(5,gaus->GetParameter(2));
          fg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Fit(func_peak,"R","",temp_mean-8*temp_sigma,temp_mean+8*temp_sigma);
        }

        float int_range_lo = sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->FindBin(temp_mean-3*temp_sigma);
        float int_range_hi = sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->FindBin(temp_mean+3*temp_sigma);
        float N_SG = sg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Integral(int_range_lo, int_range_hi);
        float N_BG = bg1d_Kpimass_vs_z[ismear][ibfield][ipid][ieta][iz]->Integral(int_range_lo, int_range_hi);

        // TODO ANSWER

        TLatex* tl = new TLatex();
        tl->SetTextAlign(11);
        tl->SetTextSize(0.03);
        tl->DrawLatexNDC(0.21,0.85,Form("%s, L = 10fb^{-1}, Q^{2} > 10 GeV^{2}",system_name[isys]));

        tl->SetTextSize(0.03);
        // tl->SetTextColor(kRed);
        tl->DrawLatexNDC(0.21,0.80,Form("#eta #subseteq [%.1f, %.1f), z #subseteq [%.2f, %.2f)",eta_lo[ieta],eta_hi[ieta],z_lo[iz],z_hi[iz]));

        gROOT->ProcessLine( Form("cc%d->Print(\"figs/kpimass_vs_z_%d_%d.pdf\")", cno-1, ieta, iz) ); // TODO change 'figs' to other folder
      }
    }
  }
}
