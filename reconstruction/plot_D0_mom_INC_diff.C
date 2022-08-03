R__LOAD_LIBRARY(libeicsmear);
#include "bins.h"
#include "TStyle.h"
#include "TGraphErrors.h"
using namespace std;

const int sys_bins = 2;
const char* fin_dirs[sys_bins] = {"./BeAGLE_v102/eAu_10_100_qhat0_nlo/outForPythiaMode/", "./BeAGLE_v102/eAu_10_100_tauforOff_qhat0_nlo/outForPythiaMode/"};
const char* sys_name[sys_bins] = {"e+Au INC on", "e+Au INC off"};
const int sys_color[sys_bins] = {kBlack, kRed};

TH1D* D0_p[sys_bins] = {0};
TH1D* D0_pt[sys_bins] = {0};

static int cno = 0;

void standardLatex()
{
  TLatex* tl = new TLatex();
  tl->SetTextAlign(11);
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.50,0.84,"Q^{2} > 10GeV^{2}, 0.05 < y < 0.8");
}

void plot_D0_mom_INC_diff()
{
  mcs(-1);

  for (int isys = 0; isys < sys_bins; isys++)
  {
    TFile* f = new TFile(Form("%smerged.root", fin_dirs[isys]), "READ");

    TTree *tree = (TTree*)f->Get("EICTree");
    erhic::EventBeagle *event(NULL);
    tree->SetBranchAddress("event",&event);

    // initialize histograms
    D0_p[isys] = new TH1D(Form("D0_p_%d", isys), Form("D0_p_%d", isys), 100, 0, 50);
    D0_pt[isys] = new TH1D(Form("D0_pt_%d", isys), Form("D0_pt_%d", isys), 100, 0, 10);

    // fill histograms
    for(int ievt = 0; ievt < tree->GetEntries(); ievt++)
    {
      tree->GetEntry(ievt);
      for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = event->GetTrack(ipart);
        if (abs(part->Id()) == 421)
        {
          D0_p[isys]->Fill(part->GetP());
          D0_pt[isys]->Fill(part->GetPt());
        }
      }
    }

    // print histograms
    mcs(cno++);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 30;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1400000;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
      htemp.GetYaxis()->SetTitle("normalized counts");
      myhset(&htemp,1.2,1.6,0.05,0.05);

      TLegend leg(0.55,0.66,0.84,0.82);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      D0_p[isys]->SetMarkerColor(sys_color[isys]);
      D0_p[isys]->SetLineColor(sys_color[isys]);
      D0_p[isys]->Draw("hsame");
      leg.AddEntry(D0_p[isys],Form("%s", sys_name[isys]), "l");

      leg.Draw("hsame");

      standardLatex();

      gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_p.pdf\")", cno-1, fin_dirs[isys]) );
    }
    mcs(cno++);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 6;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1400000;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
      htemp.GetYaxis()->SetTitle("normalized counts");
      myhset(&htemp,1.2,1.6,0.05,0.05);

      TLegend leg(0.55,0.66,0.84,0.82);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.03);
      leg.SetFillStyle(0);
      leg.SetMargin(0.1);

      D0_pt[isys]->SetMarkerColor(sys_color[isys]);
      D0_pt[isys]->SetLineColor(sys_color[isys]);
      D0_pt[isys]->Draw("hsame");
      leg.AddEntry(D0_pt[isys],Form("%s", sys_name[isys]), "l");

      leg.Draw("hsame");

      standardLatex();

      gROOT->ProcessLine( Form("cc%d->Print(\"%sD0_pt.pdf\")", cno-1, fin_dirs[isys]) );
    }

  }

  // normalize histograms to have same integral as isys ones
  for(int isys = 1; isys < sys_bins; isys++)
  {
    D0_p[isys]->Scale(D0_p[0]->Integral() / D0_p[isys]->Integral());
    D0_pt[isys]->Scale(D0_pt[0]->Integral() / D0_pt[isys]->Integral());
  }

  // print compare histograms
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 30;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 14000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}} [GeV]");
    htemp.GetYaxis()->SetTitle("normalized counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    for (int isys = 0; isys < sys_bins; isys++)
    {
      D0_p[isys]->SetMarkerColor(sys_color[isys]);
      D0_p[isys]->SetLineColor(sys_color[isys]);
      D0_p[isys]->Draw("hsame");
      leg.AddEntry(D0_p[isys],Form("%s", sys_name[isys]), "l");
    }

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"D0_p_diff.pdf\")", cno-1) );
  }
  mcs(cno++);
  {
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 6;

    float plot_yrange_lo = 0;
    float plot_yrange_hi = 14000;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p^{D^{0}}_{T} [GeV]");
    htemp.GetYaxis()->SetTitle("normalized counts");
    myhset(&htemp,1.2,1.6,0.05,0.05);

    TLegend leg(0.55,0.66,0.84,0.82);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.03);
    leg.SetFillStyle(0);
    leg.SetMargin(0.1);

    for (int isys = 0; isys < sys_bins; isys++)
    {
      D0_pt[isys]->SetMarkerColor(sys_color[isys]);
      D0_pt[isys]->SetLineColor(sys_color[isys]);
      D0_pt[isys]->Draw("hsame");
      leg.AddEntry(D0_pt[isys],Form("%s", sys_name[isys]), "l");
    }

    leg.Draw("hsame");

    standardLatex();

    gROOT->ProcessLine( Form("cc%d->Print(\"D0_pt_diff.pdf\")", cno-1) );
  }

}
