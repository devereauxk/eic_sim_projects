void mc(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
  	TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
  	gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  	gStyle->SetMarkerSize(1.6);
  	c->cd();
    sprintf(name,"mpad%d",number);	
  	mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
  	mpad->SetTopMargin(ptop);
  	mpad->SetBottomMargin(pbot);
  	mpad->SetLeftMargin(pleft);
  	mpad->SetRightMargin(pright);
  	mpad->SetLogy(1);
  	mpad->Draw();
  	mpad->cd();
	c->Modified();
	c->Update();
  }
  else {
  	TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
	tmp->Modified();
	tmp->Update();	
  }
  gSystem->ProcessEvents();
}
void mcs(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogz(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.12, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.97,0.97,0,0,0);
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogz(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogx(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogx(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogy(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogy(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogxy(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogx(1);
    mpad->SetLogy(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogyz(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.12, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.97,0.97,0,0,0);
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogy(1);
    mpad->SetLogz(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogxz(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.12, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.97,0.97,0,0,0);
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogx(1);
    mpad->SetLogz(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogxyz(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.12, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.97,0.97,0,0,0);
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogx(1);
    mpad->SetLogy(1);
    mpad->SetLogz(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mpadlog(int number, int x, int y, int z){
	if(x == 0 || x == 1) gROOT->ProcessLine( Form("mpad%d->SetLogx(%d)", number, x) );			
	if(y == 0 || y == 1) gROOT->ProcessLine( Form("mpad%d->SetLogy(%d)", number, y) );			
	if(z == 0 || z == 1) gROOT->ProcessLine( Form("mpad%d->SetLogz(%d)", number, z) );			
}
// void hset(TH1F *href, TString x, TString y, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
//   href->GetXaxis()->CenterTitle(1);
//   href->GetYaxis()->CenterTitle(1);
//   href->GetXaxis()->SetTitleOffset(xtitoff);
//   href->GetYaxis()->SetTitleOffset(ytitoff);
//   href->GetXaxis()->SetTitleSize(xtitsize);
//   href->GetYaxis()->SetTitleSize(ytitsize);
//   href->GetXaxis()->SetLabelOffset(0.01);
//   href->GetYaxis()->SetLabelOffset(0.001);
//   href->GetXaxis()->SetLabelSize(0.045);
//   href->GetYaxis()->SetLabelSize(0.045);
//   href->GetXaxis()->SetNdivisions(505);
//   href->GetYaxis()->SetNdivisions(505);
//   href->GetXaxis()->SetTitle(x);
//   href->GetYaxis()->SetTitle(y);
// }
// void hset(TH2D *href, TString x, TString y, double xtitoff = 1.,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
//   href->GetXaxis()->CenterTitle(1);
//   href->GetYaxis()->CenterTitle(1);
//   href->GetXaxis()->SetTitleOffset(xtitoff);
//   href->GetYaxis()->SetTitleOffset(ytitoff);
//   href->GetXaxis()->SetTitleSize(xtitsize);
//   href->GetYaxis()->SetTitleSize(ytitsize);
//   href->GetXaxis()->SetLabelOffset(0.01);
//   href->GetYaxis()->SetLabelOffset(0.001);
//   href->GetXaxis()->SetLabelSize(0.045);
//   href->GetYaxis()->SetLabelSize(0.045);
//   href->GetXaxis()->SetNdivisions(505);
//   href->GetYaxis()->SetNdivisions(505);
//   href->GetXaxis()->SetTitle(x);
//   href->GetYaxis()->SetTitle(y);
// }
void myhset(TH1F *href, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
  href->GetXaxis()->CenterTitle(1);
  href->GetYaxis()->CenterTitle(1);
  href->GetXaxis()->SetTitleOffset(xtitoff);
  href->GetYaxis()->SetTitleOffset(ytitoff);
  href->GetXaxis()->SetTitleSize(xtitsize);
  href->GetYaxis()->SetTitleSize(ytitsize);
  href->GetXaxis()->SetLabelOffset(0.01);
  href->GetYaxis()->SetLabelOffset(0.001);
  href->GetXaxis()->SetLabelSize(0.5);
  href->GetYaxis()->SetLabelSize(0.5);
  href->GetXaxis()->SetNdivisions(505);
  href->GetYaxis()->SetNdivisions(505);
}
void myhset(TH1D *href, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
  href->GetXaxis()->CenterTitle(1);
  href->GetYaxis()->CenterTitle(1);
  href->GetXaxis()->SetTitleOffset(xtitoff);
  href->GetYaxis()->SetTitleOffset(ytitoff);
  href->GetXaxis()->SetTitleSize(xtitsize);
  href->GetYaxis()->SetTitleSize(ytitsize);
  href->GetXaxis()->SetLabelOffset(0.01);
  href->GetYaxis()->SetLabelOffset(0.001);
  href->GetXaxis()->SetLabelSize(0.045);
  href->GetYaxis()->SetLabelSize(0.045);
  href->GetXaxis()->SetNdivisions(505);
  href->GetYaxis()->SetNdivisions(505);
}
void myhset(TH2F *href, double xtitoff = 1.0,  double ytitoff = 1.4, double xtitsize = 0.06,  double ytitsize = 0.06){
  // href->GetXaxis()->CenterTitle(1);
  // href->GetYaxis()->CenterTitle(1);
  href->GetXaxis()->SetTitleOffset(xtitoff);
  href->GetYaxis()->SetTitleOffset(ytitoff);
  href->GetXaxis()->SetTitleSize(xtitsize);
  href->GetYaxis()->SetTitleSize(ytitsize);
  href->GetXaxis()->SetLabelOffset(0.01);
  href->GetYaxis()->SetLabelOffset(0.01);
  href->GetXaxis()->SetLabelSize(0.045);
  href->GetYaxis()->SetLabelSize(0.045);
  href->GetXaxis()->SetNdivisions(505);
  href->GetYaxis()->SetNdivisions(505);
}
void myhset(TH2D *href, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
  href->GetXaxis()->CenterTitle(1);
  href->GetYaxis()->CenterTitle(1);
  href->GetXaxis()->SetTitleOffset(xtitoff);
  href->GetYaxis()->SetTitleOffset(ytitoff);
  href->GetXaxis()->SetTitleSize(xtitsize);
  href->GetYaxis()->SetTitleSize(ytitsize);
  href->GetXaxis()->SetLabelOffset(0.01);
  href->GetYaxis()->SetLabelOffset(0.01);
  href->GetXaxis()->SetLabelSize(0.045);
  href->GetYaxis()->SetLabelSize(0.045);
  href->GetXaxis()->SetNdivisions(505);
  href->GetYaxis()->SetNdivisions(505);
}
void myhset(TH3F *href, double xtitoff = 1.6,  double ytitoff = 1.6, double ztitoff = 1.4, double xtitsize = 0.06, double ytitsize = 0.06, double ztitsize = 0.06){
  href->GetXaxis()->SetTitleOffset(xtitoff);
  href->GetYaxis()->SetTitleOffset(ytitoff);
  href->GetZaxis()->SetTitleOffset(ztitoff);
  href->GetXaxis()->SetTitleSize(xtitsize);
  href->GetYaxis()->SetTitleSize(ytitsize);
  href->GetZaxis()->SetTitleSize(ztitsize);
  href->GetXaxis()->SetLabelOffset(0.01);
  href->GetYaxis()->SetLabelOffset(0.01);
  href->GetZaxis()->SetLabelOffset(0.01);
  href->GetXaxis()->SetLabelSize(0.045);
  href->GetYaxis()->SetLabelSize(0.045);
  href->GetZaxis()->SetLabelSize(0.045);
  href->GetXaxis()->SetNdivisions(505);
  href->GetYaxis()->SetNdivisions(505);
  href->GetZaxis()->SetNdivisions(505);
}
void mygraph(TGraph* gr, int mcolor = 2, double msize = 0.5, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.045,  double ytitsize = 0.045)
{
  gr->SetMarkerColor(mcolor);
  gr->SetMarkerSize(msize);
  gr->GetXaxis()->CenterTitle(1);
  gr->GetYaxis()->CenterTitle(1);
  gr->GetXaxis()->SetTitleOffset(xtitoff);
  gr->GetYaxis()->SetTitleOffset(ytitoff);
  gr->GetXaxis()->SetTitleSize(xtitsize);
  gr->GetYaxis()->SetTitleSize(ytitsize);
}
void mygrapherrors(TGraphErrors* gr, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.045,  double ytitsize = 0.045)
{
  gr->GetXaxis()->CenterTitle(1);
  gr->GetYaxis()->CenterTitle(1);
  gr->GetXaxis()->SetTitleOffset(xtitoff);
  gr->GetYaxis()->SetTitleOffset(ytitoff);
  gr->GetXaxis()->SetTitleSize(xtitsize);
  gr->GetYaxis()->SetTitleSize(ytitsize);
  gr->GetXaxis()->SetLabelOffset(0.01);
  gr->GetYaxis()->SetLabelOffset(0.01);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
}
void mygrapherrors(TGraphAsymmErrors* gr, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.045,  double ytitsize = 0.045)
{
  gr->GetXaxis()->CenterTitle(1);
  gr->GetYaxis()->CenterTitle(1);
  gr->GetXaxis()->SetTitleOffset(xtitoff);
  gr->GetYaxis()->SetTitleOffset(ytitoff);
  gr->GetXaxis()->SetTitleSize(xtitsize);
  gr->GetYaxis()->SetTitleSize(ytitsize);
}
void mymultigr(TMultiGraph* mg, double xtitoff = 1.0,  double ytitoff = 1.0, double xtitsize = 0.05,  double ytitsize = 0.05)
{
  mg->GetXaxis()->CenterTitle(1);
  mg->GetYaxis()->CenterTitle(1);
  mg->GetXaxis()->SetLabelSize(0.045);
  mg->GetYaxis()->SetLabelSize(0.045);
  mg->GetXaxis()->SetTitleOffset(xtitoff);
  mg->GetYaxis()->SetTitleOffset(ytitoff);
  mg->GetXaxis()->SetTitleSize(xtitsize);
  mg->GetYaxis()->SetTitleSize(ytitsize);
}
void mytp2d(TProfile2D *tp2d, double xtitoff = 1.0,  double ytitoff = 1., double xtitsize = 0.06,  double ytitsize = 0.07){
  tp2d->GetXaxis()->CenterTitle(1);
  tp2d->GetYaxis()->CenterTitle(1);
  tp2d->GetXaxis()->SetTitleOffset(xtitoff);
  tp2d->GetYaxis()->SetTitleOffset(ytitoff);
  tp2d->GetXaxis()->SetTitleSize(xtitsize);
  tp2d->GetYaxis()->SetTitleSize(ytitsize);
  tp2d->GetXaxis()->SetLabelOffset(0.01);
  tp2d->GetYaxis()->SetLabelOffset(0.01);
  tp2d->GetXaxis()->SetLabelSize(0.045);
  tp2d->GetYaxis()->SetLabelSize(0.045);
  tp2d->GetXaxis()->SetNdivisions(505);
  tp2d->GetYaxis()->SetNdivisions(505);
}