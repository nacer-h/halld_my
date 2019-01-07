#include <TStyle.h>
#include <TH2.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TLegend.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include <algorithm>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <Riostream.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TPaletteAxis.h"
#include "TFile.h"
#include <TProfile.h>
#include <TProfile2D.h>

void p2pifinal(TString particle1, TString particle2)
{
  gStyle->SetOptStat(0);

  TFile *fpid1 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root", particle1.Data(), particle1.Data()));
  TFile *fpid2 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root", particle2.Data(), particle2.Data()));
  TFile *fpid3 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_sep/pid%s%s.root", particle1.Data(), particle2.Data()));
  TFile *fpid4 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_misid/pid%s%s.root", particle1.Data(), particle2.Data()));
  TFile *outputfig = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_final/pid%s%s.root",particle1.Data(),particle2.Data()),"UPDATE");  // file to save in all the plots.

  ofstream foutput0 ("pid0.txt");
  ofstream foutput20 ("pid20.txt");
  ofstream foutput455 ("pid455.txt");

  // TH2F *hptheta[190];
  TH2F *hptheta;

  TH2F *hpthetantracksa[20][20];
  TH2F *hpthetantracksb[20][20];

  TH2F *hpthetanhitsa[20][20];
  TH2F *hpthetanhitsb[20][20];
  TH2F *hpthetameana[20][20];
  TH2F *hpthetameanb[20][20];

  // TH2F *hpthetaresa[190];
  TH2F *hpthetaresa[20][20];
  TH2F *hpthetabestresa = new TH2F("hpthetabestresa", Form("Momentum vs. #theta vs. optimal resolution (%s); #theta [^{0}]; P [GeV/c]; resolution [keV/cm]",particle1.Data()), 20,60,80,20,0.3,1);
  TH2F *htextresa = new TH2F("htextresa","htextresa",20,60,80,20,0.3,1);
  TH2F *htextresa2 = new TH2F("htextresa2","htextresa2",20,60,80,20,0.3,1);
  TH2F *hpthetatruncresa = new TH2F("hpthetatruncresa", Form("Momentum vs. #theta vs. optimal truncation low (%s resolution); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data()), 20,60,80,20,0.3,1);
  TH2F *hpthetatruncresa2 = new TH2F("hpthetatruncresa2", Form("Momentum vs. #theta vs. optimal truncation high (%s resolution); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data()), 20,60,80,20,0.3,1);
  TH1F *htruncresa = new TH1F("htruncresa", Form("optimal truncation low (based on resolution %s); Truncation [%%]; Counts",particle1.Data()), 100,2,102);
  TH1F *htruncresa2 = new TH1F("htruncresa2", Form("optimal truncation high (based on resolution %s); Truncation [%%]; Counts",particle1.Data()), 100,2,102);

  // TH2F *hpthetaresb[190];
  TH2F *hpthetaresb[20][20];
  TH2F *hpthetabestresb = new TH2F("hpthetabestresb", Form("Momentum vs. #theta vs. optimal resolution (%s); #theta [^{0}]; P [GeV/c]; resolution [keV/cm]",particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *htextresb = new TH2F("htextresb","htextresb",20,60,80,20,0.3,1);
  TH2F *htextresb2 = new TH2F("htextresb2","htextresb2",20,60,80,20,0.3,1);
  TH2F *hpthetatruncresb = new TH2F("hpthetatruncresb", Form("Momentum vs. #theta vs. optimal truncation low (%s resolution); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *hpthetatruncresb2 = new TH2F("hpthetatruncresb2", Form("Momentum vs. #theta vs. optimal truncation high (%s resolution); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle2.Data()), 20,60,80,20,0.3,1);
  TH1F *htruncresb = new TH1F("htruncresb", Form("optimal truncation low (based on resolution %s); Truncation [%%]; Counts",particle2.Data()), 100,2,102);
  TH1F *htruncresb2 = new TH1F("htruncresb2", Form("optimal truncation high (based on resolution %s); Truncation [%%]; Counts",particle2.Data()), 100,2,102);

  // TH2F *hpthetasep[190];
  TH2F *hpthetasep[20][20];
  TH2F *hpthetabestsep = new TH2F("hpthetabestsep", Form("Momentum vs. #theta vs. optimal separation power (%s , %s); #theta [^{0}]; P [GeV/c]; separation power",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *htextsep = new TH2F("htextsep","htextsep",20,60,80,20,0.3,1);
  TH2F *htextsep2 = new TH2F("htextsep2","htextsep2",20,60,80,20,0.3,1);
  TH2F *hpthetatruncsep = new TH2F("hpthetatruncsep", Form("Momentum vs. #theta vs. optimal truncation low (separation power (%s , %s)); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *hpthetatruncsep2 = new TH2F("hpthetatruncsep2", Form("Momentum vs. #theta vs. optimal truncation high (separation power (%s , %s)); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH1F *htruncsep = new TH1F("htruncsep", Form("optimal truncation low (based on separation power (%s , %s)); Truncation [%%]; Counts",particle1.Data(),particle2.Data()), 100,2,102);
  TH1F *htruncsep2 = new TH1F("htruncsep2", Form("optimal truncation high (based on separation power (%s , %s)); Truncation [%%]; Counts",particle1.Data(),particle2.Data()), 100,2,102);

  // TH2F* hpthetamisid[190];
  TH2F *hpthetamisid[20][20];
  TH2F *hpthetabestmisid = new TH2F("hpthetabestmisid", Form("Momentum vs. #theta vs. optimal mis-id (%s , %s); #theta [^{0}]; P [GeV/c]; mis-id",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *htextmisid = new TH2F("htextmisid","htextmisid",20,60,80,20,0.3,1);
  TH2F *htextmisid2 = new TH2F("htextsep2","htextmisid2",20,60,80,20,0.3,1);
  TH2F *hpthetatruncmisid = new TH2F("hpthetatruncmisid", Form("Momentum vs. #theta vs. optimal truncation low (mis-id (%s , %s)); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH2F *hpthetatruncmisid2 = new TH2F("hpthetatruncmisid2", Form("Momentum vs. #theta vs. optimal truncation high (mis-id (%s , %s)); #theta [^{0}]; P [GeV/c]; truncation [%%]",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);
  TH1F *htruncmisid = new TH1F("htruncmisid", Form("optimal truncation low (based on mis-id (%s , %s)); Truncation [%%]; Counts",particle1.Data(),particle2.Data()), 100,2,102);
  TH1F *htruncmisid2 = new TH1F("htruncmisid2", Form("optimal truncation high (based on mis-id (%s , %s)); Truncation [%%]; Counts",particle1.Data(),particle2.Data()), 100,2,102);

  // TGraph *grsepp = new TGraph();

  TMultiGraph *mgnhitsap = new TMultiGraph();
  TGraphErrors *ergrnhitsap0 = new TGraphErrors();
  TGraphErrors *ergrnhitsap5 = new TGraphErrors();
  TGraphErrors *ergrnhitsap20 = new TGraphErrors();
  TGraphErrors *ergrnhitsap50 = new TGraphErrors();
  TGraphErrors *ergrnhitsap455 = new TGraphErrors();
  TGraphErrors *ergrnhitsap5030 = new TGraphErrors();
  TGraphErrors *ergrnhitsap6510 = new TGraphErrors();

  TMultiGraph *mgnhitsbp = new TMultiGraph();
  TGraphErrors *ergrnhitsbp0 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp5 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp20 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp50 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp455 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp5030 = new TGraphErrors();
  TGraphErrors *ergrnhitsbp6510 = new TGraphErrors();

  TMultiGraph *mgmeanap = new TMultiGraph();
  TGraphErrors *ergrmeanap0 = new TGraphErrors();
  TGraphErrors *ergrmeanap5 = new TGraphErrors();
  TGraphErrors *ergrmeanap20 = new TGraphErrors();
  TGraphErrors *ergrmeanap50 = new TGraphErrors();
  TGraphErrors *ergrmeanap455 = new TGraphErrors();

  TMultiGraph *mgmeanbp = new TMultiGraph();
  TGraphErrors *ergrmeanbp0 = new TGraphErrors();
  TGraphErrors *ergrmeanbp5 = new TGraphErrors();
  TGraphErrors *ergrmeanbp20 = new TGraphErrors();
  TGraphErrors *ergrmeanbp50 = new TGraphErrors();
  TGraphErrors *ergrmeanbp455 = new TGraphErrors();

  TMultiGraph *mgergrresap = new TMultiGraph();
  TGraphErrors *ergrresap0 = new TGraphErrors();
  TGraphErrors *ergrresap5 = new TGraphErrors();
  TGraphErrors *ergrresap20 = new TGraphErrors();
  TGraphErrors *ergrresap50 = new TGraphErrors();
  TGraphErrors *ergrresap455 = new TGraphErrors();
  TGraphErrors *ergrresap5030 = new TGraphErrors();
  TGraphErrors *ergrresap6510 = new TGraphErrors();

  TMultiGraph *mgergrresbp = new TMultiGraph();
  TGraphErrors *ergrresbp0 = new TGraphErrors();
  TGraphErrors *ergrresbp5 = new TGraphErrors();
  TGraphErrors *ergrresbp20 = new TGraphErrors();
  TGraphErrors *ergrresbp50 = new TGraphErrors();
  TGraphErrors *ergrresbp455 = new TGraphErrors();
  TGraphErrors *ergrresbp5030 = new TGraphErrors();
  TGraphErrors *ergrresbp6510 = new TGraphErrors();

  TMultiGraph *mgergrsepp = new TMultiGraph();
  TGraphErrors *ergrsepp0 = new TGraphErrors();
  TGraphErrors *ergrsepp5 = new TGraphErrors();
  TGraphErrors *ergrsepp20 = new TGraphErrors();
  TGraphErrors *ergrsepp50 = new TGraphErrors();
  TGraphErrors *ergrsepp455 = new TGraphErrors();
  TGraphErrors *ergrsepp5030 = new TGraphErrors();
  TGraphErrors *ergrsepp6510 = new TGraphErrors();

  TMultiGraph *mgergrmisidp = new TMultiGraph();
  TGraphErrors *ergrmisidp0 = new TGraphErrors();
  TGraphErrors *ergrmisidp5 = new TGraphErrors();
  TGraphErrors *ergrmisidp20 = new TGraphErrors();
  TGraphErrors *ergrmisidp50 = new TGraphErrors();
  TGraphErrors *ergrmisidp455 = new TGraphErrors();
  TGraphErrors *ergrmisidp5030 = new TGraphErrors();
  TGraphErrors *ergrmisidp6510 = new TGraphErrors();

  TMultiGraph *mgnotresap = new TMultiGraph();
  TGraph *grnotresap0 = new TGraph();
  TGraph *grnotresap5 = new TGraph();
  TGraph *grnotresap20 = new TGraph();
  TGraph *grnotresap50 = new TGraph();
  TGraph *grnotresap455 = new TGraph();
  TGraph *grnotresap5030 = new TGraph();
  TGraph *grnotresap6510 = new TGraph();

  TMultiGraph *mgnotresbp = new TMultiGraph();
  TGraph *grnotresbp0 = new TGraph();
  TGraph *grnotresbp5 = new TGraph();
  TGraph *grnotresbp20 = new TGraph();
  TGraph *grnotresbp50 = new TGraph();
  TGraph *grnotresbp455 = new TGraph();
  TGraph *grnotresbp5030 = new TGraph();
  TGraph *grnotresbp6510 = new TGraph();

  TMultiGraph *mgnotsepp = new TMultiGraph();
  TGraph *grnotsepp0 = new TGraph();
  TGraph *grnotsepp5 = new TGraph();
  TGraph *grnotsepp20 = new TGraph();
  TGraph *grnotsepp50 = new TGraph();
  TGraph *grnotsepp455 = new TGraph();
  TGraph *grnotsepp5030 = new TGraph();
  TGraph *grnotsepp6510 = new TGraph();

  TMultiGraph *mgnotmisidp = new TMultiGraph();
  TGraph *grnotmisidp0 = new TGraph();
  TGraph *grnotmisidp5 = new TGraph();
  TGraph *grnotmisidp20 = new TGraph();
  TGraph *grnotmisidp50 = new TGraph();
  TGraph *grnotmisidp455 = new TGraph();
  TGraph *grnotmisidp5030 = new TGraph();
  TGraph *grnotmisidp6510 = new TGraph();

  TMultiGraph *mgsepp = new TMultiGraph();
  TGraphErrors *grsepp5 = new TGraphErrors();
  TGraphErrors *grsepp20 = new TGraphErrors();
  TGraphErrors *grsepp50 = new TGraphErrors();
  TGraphErrors *grsepp455 = new TGraphErrors();

  TGraph *grtruncresacorr = new TGraph();
  TGraph *grtruncresbcorr = new TGraph();
  TGraph *grtruncsepcorr = new TGraph();
  TGraph *grtruncmisidcorr = new TGraph();
  TH2D *htruncresacorr = new TH2D("htruncresacorr","low truncation vs. high truncation (based on resolution of Protons);high dE/dx truncation [%];low dE/dx truncation [%]",100,0,80,100,0,80);
  TH2D *htruncresbcorr = new TH2D("htruncresbcorr","low truncation vs. high truncation (based on resolution of #pi^{+});high dE/dx truncation [%];low dE/dx truncation [%]",100,0,80,100,0,80);
  TH2D *htruncsepcorr = new TH2D("htruncsepcorr","low truncation vs. high truncation (based on separation power ofProtons & #pi^{+});high dE/dx truncation [%];low dE/dx truncation [%]",100,0,80,100,0,80);
  TH2D *htruncmisidcorr = new TH2D("htruncmisidcorr","low truncation vs. high truncation (based on Mis-ID ofProtons & #pi^{+});high dE/dx truncation [%];low dE/dx truncation [%]",100,0,80,100,0,80);

  hptheta = (TH2F*)fpid1->Get("hptheta");
  cout <<" ***** hptheta = "<<hptheta<<endl;

  // int itrunca = 0;
  for(int trunclevel = 0; trunclevel < 100; trunclevel+=5)
  {
    for(int trunclevel2 = 0; trunclevel2 < 100; trunclevel2+=5)
    {
      // int trunclevel = 0;
      if(trunclevel + trunclevel2 > 81) continue;
      printf("------------------ trunclevel= %i | trunclevel2= %i ---------------\n",trunclevel,trunclevel2);

      TString hpthetantracksaname = Form("hpthetantracks_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
      hpthetantracksa[trunclevel/5][trunclevel2/5] = (TH2F*)fpid1->Get(hpthetantracksaname);
      // hpthetantracksa[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetantracksa[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetantracksaname<<" ***** hpthetantracksa = "<<hpthetantracksa[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetantracksbname = Form("hpthetantracks_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
      hpthetantracksb[trunclevel/5][trunclevel2/5] = (TH2F*)fpid2->Get(hpthetantracksbname);
      // hpthetantracksb[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetantracksb[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetantracksbname<<" ***** hpthetantracksb = "<<hpthetantracksb[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetanhitsaname = Form("hpthetanhits_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
      hpthetanhitsa[trunclevel/5][trunclevel2/5] = (TH2F*)fpid1->Get(hpthetanhitsaname);
      // hpthetanhitsa[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetanhitsa[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetanhitsaname<<" ***** hpthetanhitsa = "<<hpthetanhitsa[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetanhitsbname = Form("hpthetanhits_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
      hpthetanhitsb[trunclevel/5][trunclevel2/5] = (TH2F*)fpid2->Get(hpthetanhitsbname);
      // hpthetanhitsb[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetanhitsb[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetanhitsbname<<" ***** hpthetanhitsb = "<<hpthetanhitsb[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetameananame = Form("hpthetamean_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
      hpthetameana[trunclevel/5][trunclevel2/5] = (TH2F*)fpid3->Get(hpthetameananame);
      // hpthetameana[itrunca] = (TH2F*)fpid1->Get(hpthetameananame);
      // hpthetameana[itrunca]->Draw("colz");
      // hpthetameana[itrunca]->Write();
      cout <<hpthetameananame<<" ***** hpthetameana = "<<hpthetameana[trunclevel/5][trunclevel2/5]<<endl;
      // cout <<hpthetameananame<<" ***** hpthetameana = "<<hpthetameana[itrunca]<<endl;

      TString hpthetameanbname = Form("hpthetamean_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
      hpthetameanb[trunclevel/5][trunclevel2/5] = (TH2F*)fpid3->Get(hpthetameanbname);
      // hpthetameanb[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetameanb[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetameanbname<<" ***** hpthetameanb = "<<hpthetameanb[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetaresaname = Form("hpthetares_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
      hpthetaresa[trunclevel/5][trunclevel2/5] = (TH2F*)fpid3->Get(hpthetaresaname);
      // hpthetaresa[itrunca] = (TH2F*)fpid1->Get(hpthetaresaname);
      // hpthetaresa[itrunca]->Draw("colz");
      // hpthetaresa[itrunca]->Write();
      cout <<hpthetaresaname<<" ***** hpthetaresa = "<<hpthetaresa[trunclevel/5][trunclevel2/5]<<endl;
      // cout <<hpthetaresaname<<" ***** hpthetaresa = "<<hpthetaresa[itrunca]<<endl;

      TString hpthetaresbname = Form("hpthetares_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
      hpthetaresb[trunclevel/5][trunclevel2/5] = (TH2F*)fpid3->Get(hpthetaresbname);
      // hpthetaresb[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetaresb[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetaresbname<<" ***** hpthetaresb = "<<hpthetaresb[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetasepname = Form("hpthetasep_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
      hpthetasep[trunclevel/5][trunclevel2/5] = (TH2F*)fpid3->Get(hpthetasepname);
      // hpthetasep[trunclevel/5][trunclevel2/5]->GetZaxis()->SetRangeUser(0,6);
      // hpthetasep[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetasep[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetasepname<<" ***** hpthetasep = "<<hpthetasep[trunclevel/5][trunclevel2/5]<<endl;

      TString hpthetamisidname = Form("hpthetamisid_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
      hpthetamisid[trunclevel/5][trunclevel2/5] = (TH2F*)fpid4->Get(hpthetamisidname);
      // hpthetamisid[trunclevel/5][trunclevel2/5]->Draw("colz");
      // hpthetamisid[trunclevel/5][trunclevel2/5]->Write();
      cout <<hpthetamisidname<<" ***** hpthetamisid = "<<hpthetamisid[trunclevel/5][trunclevel2/5]<<endl;
    }
  }

  double ntracksap0[20];
  double ntracksap5[20];
  double ntracksap20[20];
  double ntracksap50[20];
  double ntracksap455[20];

  double ntracksbp0[20];
  double ntracksbp5[20];
  double ntracksbp20[20];
  double ntracksbp50[20];
  double ntracksbp455[20];

  double nhitsap0[20];
  double nhitsap5[20];
  double nhitsap20[20];
  double nhitsap50[20];
  double nhitsap455[20];
  double nhitsap5030[20];
  double nhitsap6510[20];

  double nhitsbp0[20];
  double nhitsbp5[20];
  double nhitsbp20[20];
  double nhitsbp50[20];
  double nhitsbp455[20];
  double nhitsbp5030[20];
  double nhitsbp6510[20];

  double meanap0[20];
  double meanap5[20];
  double meanap20[20];
  double meanap50[20];
  double meanap455[20];

  double meanbp0[20];
  double meanbp5[20];
  double meanbp20[20];
  double meanbp50[20];
  double meanbp455[20];

  double resap0[20];
  double resap5[20];
  double resap20[20];
  double resap50[20];
  double resap455[20];
  double resap5030[20];
  double resap6510[20];

  double resbp0[20];
  double resbp5[20];
  double resbp20[20];
  double resbp50[20];
  double resbp455[20];
  double resbp5030[20];
  double resbp6510[20];

  double sepp0[20];
  double sepp5[20];
  double sepp20[20];
  double sepp50[20];
  double sepp455[20];
  double sepp5030[20];
  double sepp6510[20];

  double misidp0[20];
  double misidp5[20];
  double misidp20[20];
  double misidp50[20];
  double misidp455[20];
  double misidp5030[20];
  double misidp6510[20];

  int itruncresa = 0;
  int itruncresb = 0;
  int itruncsep = 0;
  int itruncmisid = 0;

  double sumbins = 0;

  int nseptheta = 0;
  foutput0<<"Left-trunc[%] Right-trunc[%] momentum [GeV/c] theta[deg] #hits(Protons) #hits(Pi+) #tracks(Protons) #tracks(Pi+) separation power"<<endl;
  foutput20<<"Left-trunc[%] Right-trunc[%] momentum [GeV/c] theta[deg] #hits(Protons) #hits(Pi+) #tracks(Protons) #tracks(Pi+) separation power"<<endl;
  foutput455<<"Left-trunc[%] Right-trunc[%] momentum [GeV/c] theta[deg] #hits(Protons) #hits(Pi+) #tracks(Protons) #tracks(Pi+) separation power"<<endl;

  for(int ip=0; ip<20; ++ip)
  {
    int notresa0 = 0;
    int notresa5 = 0;
    int notresa20 = 0;
    int notresa50 = 0;
    int notresa455 = 0;
    int notresa5030 = 0;
    int notresa6510 = 0;

    int notresb0 = 0;
    int notresb5 = 0;
    int notresb20 = 0;
    int notresb50 = 0;
    int notresb455 = 0;
    int notresb5030 = 0;
    int notresb6510 = 0;

    int notsep0 = 0;
    int notsep5 = 0;
    int notsep20 = 0;
    int notsep50 = 0;
    int notsep455 = 0;
    int notsep5030 = 0;
    int notsep6510 = 0;

    int notmisid0 = 0;
    int notmisid5 = 0;
    int notmisid20 = 0;
    int notmisid50 = 0;
    int notmisid455 = 0;
    int notmisid5030 = 0;
    int notmisid6510 = 0;

    for(int itheta=0; itheta<20; ++itheta)
    {
      double resa = 1000;
      double resb = 1000;
      double sep = -1;
      double misid   = 1000;
      int truncresa =1000;
      int truncresa2 = 1000;
      int truncresb =1000;
      int truncresb2 = 1000;
      int truncsep = 0;
      int truncsep2 = 0;
      int truncmisid = 1000;
      int truncmisid2 = 1000;
      // int itruncat = 0;

      for(int trunclevel=0; trunclevel<100; trunclevel+=5)
      {
        for(int trunclevel2=0; trunclevel2<100; trunclevel2+=5)
        {
          // int trunclevel = 0;
          if(trunclevel + trunclevel2 > 81) continue;

          printf("_______________________________ ip= %i | itheta= %i | trunclevel= %i | trunclevel2= %i _________________________\n",ip,itheta,trunclevel,trunclevel2);

          if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0 && hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1) < resa)
          {
            resa = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            truncresa = trunclevel;
            truncresa2 = trunclevel2;
            printf("****** resa= %f | truncresa = %i | truncresa2 = %i ******\n",resa,truncresa,truncresa2);
          }
          if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0 && hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1) < resb)
          {
            resb = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            truncresb = trunclevel;
            truncresb2 = trunclevel2;
            printf("****** resb= %f | truncresb = %i | truncresb2 = %i ******\n",resb,truncresb,truncresb2);
          }

          if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1) > 0 && hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1) > sep)
          {
            sep = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            truncsep = trunclevel;
            truncsep2 = trunclevel2;
            printf("******* sep = %f | truncsep = %i | truncsep2 = %i ********\n", sep, truncsep,truncsep2);
          }

          if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8 && hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1) < misid)
          {
            misid = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            truncmisid = trunclevel;
            truncmisid2 = trunclevel2;
            printf("****** misid= %f| truncmisid = %i | truncmisid2 = %i ******\n", misid, truncmisid, truncmisid2);
          }

          if(trunclevel == 0 && trunclevel2 == 0)
          {
            // if(hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksap0[itheta] = hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksbp0[itheta] = hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap0[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)

            nhitsbp0[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanap0[itheta] = hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanbp0[itheta] = hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap0[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp0[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp0[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp0[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa0;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb0;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep0;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid0;

            foutput0<<"    "<<trunclevel<<"              "<<trunclevel2<<"           "<<hptheta->GetYaxis()->GetBinCenter(ip+1)<<"           "<<
            hptheta->GetXaxis()->GetBinCenter(itheta+1)<<"         "<<nhitsap0[itheta]<<"      "<<nhitsbp0[itheta]<<"       "<<ntracksap0[itheta]<<
            "       "<<ntracksbp0[itheta]<<"         "<<sepp0[itheta]<<endl;
            printf("++++++++++++++++++++++++++++++++++++++++ misidp5[itheta] = %f\n", misidp5[itheta]);
          }

          if(trunclevel == 0 && trunclevel2 == 5)
          {
            // if(hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksap5[itheta] = hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksbp5[itheta] = hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap5[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp5[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanap5[itheta] = hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanbp5[itheta] = hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap5[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp5[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp5[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp5[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa5;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb5;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep5;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid5;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp5[itheta] = %f\n", misidp5[itheta]);
          }
          if(trunclevel == 0 && trunclevel2 == 20)
          {
            // if(hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksap20[itheta] = hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksbp20[itheta] = hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap20[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp20[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanap20[itheta] = hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanbp20[itheta] = hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap20[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp20[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp20[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp20[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa20;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb20;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep20;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid20;

            foutput20<<"    "<<trunclevel<<"              "<<trunclevel2<<"           "<<hptheta->GetYaxis()->GetBinCenter(ip+1)<<"           "<<
            hptheta->GetXaxis()->GetBinCenter(itheta+1)<<"         "<<nhitsap20[itheta]<<"      "<<nhitsbp20[itheta]<<"       "<<ntracksap20[itheta]<<
            "       "<<ntracksbp20[itheta]<<"         "<<sepp20[itheta]<<endl;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp20[itheta] = %f\n", misidp20[itheta]);
          }
          if(trunclevel == 0 && trunclevel2 == 50)
          {
            // if(hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksap50[itheta] = hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksbp50[itheta] = hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap50[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp50[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanap50[itheta] = hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanbp50[itheta] = hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap50[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp50[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp50[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp50[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa50;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb50;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep50;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid50;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp5[itheta] = %f\n", misidp50[itheta]);
          }
          if(trunclevel == 45 && trunclevel2 == 5)
          {
            // if(hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksap455[itheta] = hpthetantracksa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            ntracksbp455[itheta] = hpthetantracksb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap455[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp455[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanap455[itheta] = hpthetameana[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            meanbp455[itheta] = hpthetameanb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap455[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp455[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp455[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp455[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa455;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb455;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep455;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid455;

            foutput455<<"    "<<trunclevel<<"              "<<trunclevel2<<"           "<<hptheta->GetYaxis()->GetBinCenter(ip+1)<<"           "<<
            hptheta->GetXaxis()->GetBinCenter(itheta+1)<<"         "<<nhitsap455[itheta]<<"      "<<nhitsbp455[itheta]<<"       "<<ntracksap455[itheta]<<
            "       "<<ntracksbp455[itheta]<<"         "<<sepp455[itheta]<<endl;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp455[itheta] = %f\n", misidp455[itheta]);
          }
          if(trunclevel == 50 && trunclevel2 == 30)
          {
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap5030[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp5030[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap5030[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp5030[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp5030[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp5030[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa5030;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb5030;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep5030;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid5030;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp5030[itheta] = %f\n", misidp5030[itheta]);
          }
          if(trunclevel == 65 && trunclevel2 == 10)
          {
            // if(hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsap6510[itheta] = hpthetanhitsa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            // if(hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            nhitsbp6510[itheta] = hpthetanhitsb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resap6510[itheta] = hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            resbp6510[itheta] = hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>0)
            sepp6510[itheta] = hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)>1e-8)
            misidp6510[itheta] = hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1);

            if(hpthetaresa[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresa6510;
            if(hpthetaresb[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notresb6510;
            if(hpthetasep[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notsep6510;
            if(hpthetamisid[trunclevel/5][trunclevel2/5]->GetBinContent(itheta+1, ip+1)==0) ++notmisid6510;

            printf("++++++++++++++++++++++++++++++++++++++++ misidp6510[itheta] = %f\n", misidp6510[itheta]);
          }
          // itruncat++;
        }
      }

      double binentries = hptheta->GetBinContent(itheta+1, ip+1);
      sumbins += binentries;

      if(resa>0 && resa<1000)
      {
        hpthetabestresa->SetBinContent(itheta+1, ip+1, resa);
        htextresa->SetBinContent(itheta+1, ip+1, truncresa);
        htextresa2->SetBinContent(itheta+1, ip+1, truncresa2);
        hpthetatruncresa->SetBinContent(itheta+1, ip+1, truncresa);
        hpthetatruncresa2->SetBinContent(itheta+1, ip+1, truncresa2);
        grtruncresacorr->SetPoint(itruncresa, truncresa,truncresa2);
        htruncresacorr->Fill(truncresa,truncresa2);
        htruncresa->Fill(truncresa,binentries);
        htruncresa2->Fill(truncresa2,binentries);
        itruncresa++;
      }
      if(resb>0 && resb<1000)
      {
        hpthetabestresb->SetBinContent(itheta+1, ip+1, resb);
        htextresb->SetBinContent(itheta+1, ip+1, truncresb);
        htextresb2->SetBinContent(itheta+1, ip+1, truncresb2);
        hpthetatruncresb->SetBinContent(itheta+1, ip+1, truncresb);
        hpthetatruncresb2->SetBinContent(itheta+1, ip+1, truncresb2);
        grtruncresbcorr->SetPoint(itruncresb, truncresb,truncresb2);
        htruncresbcorr->Fill(truncresb,truncresb2);
        htruncresb->Fill(truncresb,binentries);
        htruncresb2->Fill(truncresb2,binentries);
        itruncresb++;
      }
      if(sep>0)
      {
        //sepp = sep;
        hpthetabestsep->SetBinContent(itheta+1, ip+1, sep);
        htextsep->SetBinContent(itheta+1, ip+1, truncsep);
        htextsep2->SetBinContent(itheta+1, ip+1, truncsep2);
        hpthetatruncsep->SetBinContent(itheta+1, ip+1, truncsep);
        hpthetatruncsep2->SetBinContent(itheta+1, ip+1, truncsep2);
        grtruncsepcorr->SetPoint(itruncsep, truncsep,truncsep2);
        htruncsepcorr->Fill(truncsep,truncsep2);
        htruncsep->Fill(truncsep,binentries);
        htruncsep2->Fill(truncsep2,binentries);
        itruncsep++;
      }
      if(misid>0 && misid<1000)
      {
        hpthetabestmisid->SetBinContent(itheta+1, ip+1, misid);
        htextmisid->SetBinContent(itheta+1, ip+1, truncmisid);
        htextmisid2->SetBinContent(itheta+1, ip+1, truncmisid2);
        hpthetatruncmisid->SetBinContent(itheta+1, ip+1, truncmisid);
        hpthetatruncmisid2->SetBinContent(itheta+1, ip+1, truncmisid2);
        grtruncmisidcorr->SetPoint(itruncmisid, truncmisid,truncmisid2);
        htruncmisidcorr->Fill(truncmisid,truncmisid2);
        htruncmisid->Fill(truncmisid,binentries);
        htruncmisid2->Fill(truncmisid2,binentries);
        itruncmisid++;
      }
    }

    double mom = hptheta->GetYaxis()->GetBinCenter(ip+1);

    ergrnhitsap0->SetPoint(ip,mom,TMath::Mean(20,nhitsap0));
    ergrnhitsap0->SetPointError(ip,0,TMath::RMS(20,nhitsap0));
    ergrnhitsap5->SetPoint(ip,mom-0.003,TMath::Mean(20,nhitsap5));
    ergrnhitsap5->SetPointError(ip,0,TMath::RMS(20,nhitsap5));
    ergrnhitsap20->SetPoint(ip,mom-0.006,TMath::Mean(20,nhitsap20));
    ergrnhitsap20->SetPointError(ip,0,TMath::RMS(20,nhitsap20));
    ergrnhitsap50->SetPoint(ip,mom-0.009,TMath::Mean(20,nhitsap50));
    ergrnhitsap50->SetPointError(ip,0,TMath::RMS(20,nhitsap50));
    ergrnhitsap455->SetPoint(ip,mom-0.012,TMath::Mean(20,nhitsap455));
    ergrnhitsap455->SetPointError(ip,0,TMath::RMS(20,nhitsap455));
    ergrnhitsap5030->SetPoint(ip,mom-0.015,TMath::Mean(20,nhitsap5030));
    ergrnhitsap5030->SetPointError(ip,0,TMath::RMS(20,nhitsap5030));
    ergrnhitsap6510->SetPoint(ip,mom-0.018,TMath::Mean(20,nhitsap6510));
    ergrnhitsap6510->SetPointError(ip,0,TMath::RMS(20,nhitsap6510));

    ergrnhitsbp0->SetPoint(ip,mom,TMath::Mean(20,nhitsbp0));
    ergrnhitsbp0->SetPointError(ip,0,TMath::RMS(20,nhitsbp0));
    ergrnhitsbp5->SetPoint(ip,mom-0.003,TMath::Mean(20,nhitsbp5));
    ergrnhitsbp5->SetPointError(ip,0,TMath::RMS(20,nhitsbp5));
    ergrnhitsbp20->SetPoint(ip,mom-0.006,TMath::Mean(20,nhitsbp20));
    ergrnhitsbp20->SetPointError(ip,0,TMath::RMS(20,nhitsbp20));
    ergrnhitsbp50->SetPoint(ip,mom-0.009,TMath::Mean(20,nhitsbp50));
    ergrnhitsbp50->SetPointError(ip,0,TMath::RMS(20,nhitsbp50));
    ergrnhitsbp455->SetPoint(ip,mom-0.012,TMath::Mean(20,nhitsbp455));
    ergrnhitsbp455->SetPointError(ip,0,TMath::RMS(20,nhitsbp455));
    ergrnhitsbp5030->SetPoint(ip,mom-0.015,TMath::Mean(20,nhitsbp5030));
    ergrnhitsbp5030->SetPointError(ip,0,TMath::RMS(20,nhitsbp5030));
    ergrnhitsbp6510->SetPoint(ip,mom-0.018,TMath::Mean(20,nhitsbp6510));
    ergrnhitsbp6510->SetPointError(ip,0,TMath::RMS(20,nhitsbp6510));

    ergrmeanap0->SetPoint(ip,mom,TMath::Mean(20,meanap0));
    ergrmeanap0->SetPointError(ip,0,TMath::RMS(20,meanap0));
    ergrmeanap5->SetPoint(ip,mom+0.004,TMath::Mean(20,meanap5));
    ergrmeanap5->SetPointError(ip,0,TMath::RMS(20,meanap5));
    ergrmeanap20->SetPoint(ip,mom+0.009,TMath::Mean(20,meanap20));
    ergrmeanap20->SetPointError(ip,0,TMath::RMS(20,meanap20));
    ergrmeanap50->SetPoint(ip,mom+0.014,TMath::Mean(20,meanap50));
    ergrmeanap50->SetPointError(ip,0,TMath::RMS(20,meanap50));
    ergrmeanap455->SetPoint(ip,mom+0.019,TMath::Mean(20,meanap455));
    ergrmeanap455->SetPointError(ip,0,TMath::RMS(20,meanap455));

    ergrmeanbp0->SetPoint(ip,mom,TMath::Mean(20,meanbp0));
    ergrmeanbp0->SetPointError(ip,0,TMath::RMS(20,meanbp0));
    ergrmeanbp5->SetPoint(ip,mom+0.004,TMath::Mean(20,meanbp5));
    ergrmeanbp5->SetPointError(ip,0,TMath::RMS(20,meanbp5));
    ergrmeanbp20->SetPoint(ip,mom+0.009,TMath::Mean(20,meanbp20));
    ergrmeanbp20->SetPointError(ip,0,TMath::RMS(20,meanbp20));
    ergrmeanbp50->SetPoint(ip,mom+0.014,TMath::Mean(20,meanbp50));
    ergrmeanbp50->SetPointError(ip,0,TMath::RMS(20,meanbp50));
    ergrmeanbp455->SetPoint(ip,mom+0.019,TMath::Mean(20,meanbp455));
    ergrmeanbp455->SetPointError(ip,0,TMath::RMS(20,meanbp455));

    ergrresap0->SetPoint(ip,mom,TMath::Mean(20,resap0));
    ergrresap0->SetPointError(ip,0,TMath::RMS(20,resap0));
    ergrresap5->SetPoint(ip,mom-0.003,TMath::Mean(20,resap5));
    ergrresap5->SetPointError(ip,0,TMath::RMS(20,resap5));
    ergrresap20->SetPoint(ip,mom-0.006,TMath::Mean(20,resap20));
    ergrresap20->SetPointError(ip,0,TMath::RMS(20,resap20));
    ergrresap50->SetPoint(ip,mom-0.009,TMath::Mean(20,resap50));
    ergrresap50->SetPointError(ip,0,TMath::RMS(20,resap50));
    ergrresap455->SetPoint(ip,mom-0.012,TMath::Mean(20,resap455));
    ergrresap455->SetPointError(ip,0,TMath::RMS(20,resap455));
    ergrresap5030->SetPoint(ip,mom-0.015,TMath::Mean(20,resap5030));
    ergrresap5030->SetPointError(ip,0,TMath::RMS(20,resap5030));
    ergrresap6510->SetPoint(ip,mom-0.018,TMath::Mean(20,resap6510));
    ergrresap6510->SetPointError(ip,0,TMath::RMS(20,resap6510));

    ergrresbp0->SetPoint(ip,mom,TMath::Mean(20,resbp0));
    ergrresbp0->SetPointError(ip,0,TMath::RMS(20,resbp0));
    ergrresbp5->SetPoint(ip,mom-0.003,TMath::Mean(20,resbp5));
    ergrresbp5->SetPointError(ip,0,TMath::RMS(20,resbp5));
    ergrresbp20->SetPoint(ip,mom-0.006,TMath::Mean(20,resbp20));
    ergrresbp20->SetPointError(ip,0,TMath::RMS(20,resbp20));
    ergrresbp50->SetPoint(ip,mom-0.009,TMath::Mean(20,resbp50));
    ergrresbp50->SetPointError(ip,0,TMath::RMS(20,resbp50));
    ergrresbp455->SetPoint(ip,mom-0.012,TMath::Mean(20,resbp455));
    ergrresbp455->SetPointError(ip,0,TMath::RMS(20,resbp455));
    ergrresbp5030->SetPoint(ip,mom-0.015,TMath::Mean(20,resbp5030));
    ergrresbp5030->SetPointError(ip,0,TMath::RMS(20,resbp5030));
    ergrresbp6510->SetPoint(ip,mom-0.018,TMath::Mean(20,resbp6510));
    ergrresbp6510->SetPointError(ip,0,TMath::RMS(20,resbp6510));

    ergrsepp0->SetPoint(ip,mom,TMath::Mean(20,sepp0));
    ergrsepp0->SetPointError(ip,0,TMath::RMS(20,sepp0));
    ergrsepp5->SetPoint(ip,mom-0.003,TMath::Mean(20,sepp5));
    ergrsepp5->SetPointError(ip,0,TMath::RMS(20,sepp5));
    ergrsepp20->SetPoint(ip,mom-0.006,TMath::Mean(20,sepp20));
    ergrsepp20->SetPointError(ip,0,TMath::RMS(20,sepp20));
    ergrsepp50->SetPoint(ip,mom-0.009,TMath::Mean(20,sepp50));
    ergrsepp50->SetPointError(ip,0,TMath::RMS(20,sepp50));
    ergrsepp455->SetPoint(ip,mom-0.012,TMath::Mean(20,sepp455));
    ergrsepp455->SetPointError(ip,0,TMath::RMS(20,sepp455));
    ergrsepp5030->SetPoint(ip,mom-0.015,TMath::Mean(20,sepp5030));
    ergrsepp5030->SetPointError(ip,0,TMath::RMS(20,sepp5030));
    ergrsepp6510->SetPoint(ip,mom-0.018,TMath::Mean(20,sepp6510));
    ergrsepp6510->SetPointError(ip,0,TMath::RMS(20,sepp6510));

    ergrmisidp0->SetPoint(ip,mom,TMath::Mean(20,misidp0));
    ergrmisidp0->SetPointError(ip,0,TMath::RMS(20,misidp0));
    ergrmisidp5->SetPoint(ip,mom-0.003,TMath::Mean(20,misidp5));
    ergrmisidp5->SetPointError(ip,0,TMath::RMS(20,misidp5));
    ergrmisidp20->SetPoint(ip,mom-0.006,TMath::Mean(20,misidp20));
    ergrmisidp20->SetPointError(ip,0,TMath::RMS(20,misidp20));
    ergrmisidp50->SetPoint(ip,mom-0.009,TMath::Mean(20,misidp50));
    ergrmisidp50->SetPointError(ip,0,TMath::RMS(20,misidp50));
    ergrmisidp455->SetPoint(ip,mom-0.012,TMath::Mean(20,misidp455));
    ergrmisidp455->SetPointError(ip,0,TMath::RMS(20,misidp455));
    ergrmisidp5030->SetPoint(ip,mom-0.015,TMath::Mean(20,misidp5030));
    ergrmisidp5030->SetPointError(ip,0,TMath::RMS(20,misidp5030));
    ergrmisidp6510->SetPoint(ip,mom-0.018,TMath::Mean(20,misidp6510));
    ergrmisidp6510->SetPointError(ip,0,TMath::RMS(20,misidp6510));

    grnotresap0->SetPoint(ip,mom-0.003,notresa0);
    grnotresap5->SetPoint(ip,mom+0.003,notresa5);
    grnotresap20->SetPoint(ip,mom+0.006,notresa20);
    grnotresap50->SetPoint(ip,mom+0.009,notresa50);
    grnotresap455->SetPoint(ip,mom+0.012,notresa455);
    grnotresap5030->SetPoint(ip,mom+0.015,notresa5030);
    grnotresap6510->SetPoint(ip,mom+0.018,notresa6510);

    grnotresbp0->SetPoint(ip,mom-0.003,notresb0);
    grnotresbp5->SetPoint(ip,mom+0.003,notresb5);
    grnotresbp20->SetPoint(ip,mom+0.006,notresb20);
    grnotresbp50->SetPoint(ip,mom+0.009,notresb50);
    grnotresbp455->SetPoint(ip,mom+0.012,notresb455);
    grnotresbp5030->SetPoint(ip,mom+0.015,notresb5030);
    grnotresbp6510->SetPoint(ip,mom+0.018,notresb6510);

    grnotsepp0->SetPoint(ip,mom-0.003,notsep0);
    grnotsepp5->SetPoint(ip,mom+0.003,notsep5);
    grnotsepp20->SetPoint(ip,mom+0.006,notsep20);
    grnotsepp50->SetPoint(ip,mom+0.009,notsep50);
    grnotsepp455->SetPoint(ip,mom+0.012,notsep455);
    grnotsepp5030->SetPoint(ip,mom+0.015,notsep5030);
    grnotsepp6510->SetPoint(ip,mom+0.018,notsep6510);

    grnotmisidp0->SetPoint(ip,mom-0.003,notmisid0);
    grnotmisidp5->SetPoint(ip,mom+0.003,notmisid5);
    grnotmisidp20->SetPoint(ip,mom+0.006,notmisid20);
    grnotmisidp50->SetPoint(ip,mom+0.009,notmisid50);
    grnotmisidp455->SetPoint(ip,mom+0.012,notmisid455);
    grnotmisidp5030->SetPoint(ip,mom+0.015,notmisid5030);
    grnotmisidp6510->SetPoint(ip,mom+0.018,notmisid6510);

    grsepp5->SetPoint(ip,mom+0.009,TMath::MaxElement(20,sepp5));
    grsepp20->SetPoint(ip,mom+0.014,TMath::MaxElement(20,sepp20));
    grsepp50->SetPoint(ip,mom+0.019,TMath::MaxElement(20,sepp50));
    grsepp455->SetPoint(ip,mom,TMath::MaxElement(20,sepp455));


    cout<<"****** notsep0 = "<<notsep0<<" | notsep5 = "<<notsep5<<" | notsep20 = "<<notsep20<<" | notsep50 = "<<notsep50<<" | notsep455 = "<<notsep455<<endl;
  }

  TCanvas *cpthetabestresa = new TCanvas("cpthetabestresa","cpthetabestresa",700,500);
  cpthetabestresa->cd();
  gStyle->SetPalette(1);
  hpthetabestresa->GetXaxis()->SetLabelSize(0.05);
  hpthetabestresa->GetXaxis()->SetTitleSize(0.05);
  hpthetabestresa->GetXaxis()->SetTitleOffset(0.95);
  hpthetabestresa->GetYaxis()->SetLabelSize(0.05);
  hpthetabestresa->GetYaxis()->SetTitleSize(0.05);
  hpthetabestresa->GetYaxis()->SetTitleOffset(0.95);
  hpthetabestresa->GetZaxis()->SetLabelSize(0.05);
  hpthetabestresa->GetZaxis()->SetTitleSize(0.05);
  hpthetabestresa->GetZaxis()->SetTitleOffset(0.95);
  // hpthetabestresa->GetZaxis()->SetRangeUser(0,2);
  hpthetabestresa->GetZaxis()->SetNdivisions(515);
  htextresa->SetMarkerSize(1.);
  htextresa2->SetMarkerSize(1.);
  htextresa2->SetMarkerColor(0);
  hpthetabestresa->Draw("colz");
  htextresa->SetBarOffset(0.21);
  htextresa->Draw("TEXT SAME");
  htextresa2->SetBarOffset(-1.*0.21);
  htextresa2->Draw("TEXT SAME");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpthetabestresa->GetListOfFunctions()->FindObject("palette");
  cpthetabestresa->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetabestresa->Write();
  cpthetabestresa->Modified();
  cpthetabestresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestres_%s.eps",particle1.Data()),"eps");
  cpthetabestresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestres_%s.root",particle1.Data()),"root");

  TCanvas *cpthetatruncresa = new TCanvas("cpthetatruncresa","cpthetatruncresa",700,500);
  cpthetatruncresa->cd();
  hpthetatruncresa->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncresa->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncresa->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncresa->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncresa->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncresa->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncresa->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncresa->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncresa->GetZaxis()->SetTitleOffset(0.95);
  // hpthetatruncresa->GetZaxis()->SetRangeUser(0,2);
  hpthetatruncresa->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncresa->GetZaxis()->SetNdivisions(515);
  hpthetatruncresa->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncresa->GetListOfFunctions()->FindObject("palette");
  cpthetatruncresa->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetatruncresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres_%s.eps",particle1.Data()),"eps");
  cpthetatruncresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres_%s.root",particle1.Data()),"root");

  TCanvas *cpthetatruncresa2 = new TCanvas("cpthetatruncresa2","cpthetatruncresa2",700,500);
  cpthetatruncresa2->cd();
  hpthetatruncresa2->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncresa2->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncresa2->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncresa2->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncresa2->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncresa2->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncresa2->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncresa2->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncresa2->GetZaxis()->SetTitleOffset(0.95);
  // hpthetatruncresa2->GetZaxis()->SetRangeUser(0,2);
  hpthetatruncresa2->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncresa2->GetZaxis()->SetNdivisions(515);
  hpthetatruncresa2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncresa2->GetListOfFunctions()->FindObject("palette");
  cpthetatruncresa2->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetatruncresa2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres2_%s.eps",particle2.Data()),"eps");
  cpthetatruncresa2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres2_%s.root",particle2.Data()),"root");

  TCanvas *cpthetabestresb = new TCanvas("cpthetabestresb","cpthetabestresb",700,500);
  cpthetabestresb->cd();
  gStyle->SetPalette(1);
  hpthetabestresb->GetXaxis()->SetLabelSize(0.05);
  hpthetabestresb->GetXaxis()->SetTitleSize(0.05);
  hpthetabestresb->GetXaxis()->SetTitleOffset(0.95);
  hpthetabestresb->GetYaxis()->SetLabelSize(0.05);
  hpthetabestresb->GetYaxis()->SetTitleSize(0.05);
  hpthetabestresb->GetYaxis()->SetTitleOffset(0.95);
  hpthetabestresb->GetZaxis()->SetLabelSize(0.05);
  hpthetabestresb->GetZaxis()->SetTitleSize(0.05);
  hpthetabestresb->GetZaxis()->SetTitleOffset(0.95);
  // hpthetabestresb->GetZaxis()->SetRangeUser(0,2);
  hpthetabestresb->GetZaxis()->SetNdivisions(515);
  htextresb->SetMarkerSize(1.);
  htextresb2->SetMarkerSize(1.);
  htextresb2->SetMarkerColor(4);
  hpthetabestresb->Draw("colz");
  htextresb->SetBarOffset(0.21);
  htextresb->Draw("TEXT SAME");
  htextresb2->SetBarOffset(-1.*0.21);
  htextresb2->Draw("TEXT SAME");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetabestresb->GetListOfFunctions()->FindObject("palette");
  cpthetabestresb->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetabestresb->Write();
  cpthetabestresb->Modified();
  cpthetabestresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestres_%s.eps",particle2.Data()),"eps");
  cpthetabestresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestres_%s.root",particle2.Data()),"root");

  TCanvas *cpthetatruncresb = new TCanvas("cpthetatruncresb","cpthetatruncresb",700,500);
  cpthetatruncresb->cd();
  hpthetatruncresb->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncresb->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncresb->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncresb->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncresb->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncresb->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncresb->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncresb->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncresb->GetZaxis()->SetTitleOffset(0.95);
  // hpthetatruncresb->GetZaxis()->SetRangeUser(0,2);
  hpthetatruncresb->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncresb->GetZaxis()->SetNdivisions(515);
  hpthetatruncresb->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncresb->GetListOfFunctions()->FindObject("palette");
  cpthetatruncresb->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetatruncresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres_%s.eps",particle2.Data()),"eps");
  cpthetatruncresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres_%s.root",particle2.Data()),"root");

  TCanvas *cpthetatruncresb2 = new TCanvas("cpthetatruncresb2","cpthetatruncresb2",700,500);
  cpthetatruncresb2->cd();
  hpthetatruncresb2->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncresb2->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncresb2->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncresb2->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncresb2->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncresb2->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncresb2->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncresb2->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncresb2->GetZaxis()->SetTitleOffset(0.95);
  // hpthetatruncresb2->GetZaxis()->SetRangeUser(0,2);
  hpthetatruncresb2->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncresb2->GetZaxis()->SetNdivisions(515);
  hpthetatruncresb2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncresb2->GetListOfFunctions()->FindObject("palette");
  cpthetatruncresb2->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetatruncresb2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres2_%s.eps",particle2.Data()),"eps");
  cpthetatruncresb2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncres2_%s.root",particle2.Data()),"root");

  TCanvas *cpthetabestsep = new TCanvas("cpthetabestsep","cpthetabestsep",700,500);
  cpthetabestsep->cd();
  gStyle->SetPalette(1);
  hpthetabestsep->SetTitleSize(0.05,"XYZ");
  hpthetabestsep->SetLabelSize(0.05,"XYZ");
  hpthetabestsep->SetTitleOffset(0.95,"XYZ");
  hpthetabestsep->GetZaxis()->SetRangeUser(0,7);
  hpthetabestsep->GetZaxis()->SetNdivisions(515);
  htextsep->SetMarkerSize(1.);
  htextsep2->SetMarkerSize(1.);
  htextsep2->SetMarkerColor(4);
  hpthetabestsep->Draw("colz");
  htextsep->SetBarOffset(0.21);
  htextsep->Draw("TEXT SAME");
  htextsep2->SetBarOffset(-1.*0.21);
  htextsep2->Draw("TEXT SAME");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetabestsep->GetListOfFunctions()->FindObject("palette");
  cpthetabestsep->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetabestsep->Write();
  cpthetabestsep->Modified();
  // hpthetabestsep->Write();
  cpthetabestsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestsep_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetabestsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestsep_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cpthetatruncsep = new TCanvas("cpthetatruncsep","cpthetatruncsep",700,500);
  cpthetatruncsep->cd();
  hpthetatruncsep->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncsep->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncsep->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncsep->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncsep->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncsep->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncsep->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncsep->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncsep->GetZaxis()->SetTitleOffset(0.95);
  hpthetatruncsep->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncsep->GetZaxis()->SetNdivisions(515);
  hpthetatruncsep->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncsep->GetListOfFunctions()->FindObject("palette");
  cpthetatruncsep->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetatruncsep->Write();
  cpthetatruncsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncsep_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetatruncsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncsep_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cpthetatruncsep2 = new TCanvas("cpthetatruncsep2","cpthetatruncsep2",700,500);
  cpthetatruncsep2->cd();
  hpthetatruncsep2->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncsep2->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncsep2->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncsep2->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncsep2->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncsep2->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncsep2->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncsep2->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncsep2->GetZaxis()->SetTitleOffset(0.95);
  hpthetatruncsep2->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncsep2->GetZaxis()->SetNdivisions(515);
  hpthetatruncsep2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncsep2->GetListOfFunctions()->FindObject("palette");
  cpthetatruncsep2->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetatruncsep2->Write();
  cpthetatruncsep2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncsep2_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetatruncsep2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncsep2_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cpthetabestmisid = new TCanvas("cpthetabestmisid","cpthetabestmisid",700,500);
  cpthetabestmisid->cd();
  hpthetabestmisid->GetXaxis()->SetLabelSize(0.05);
  hpthetabestmisid->GetXaxis()->SetTitleSize(0.05);
  hpthetabestmisid->GetXaxis()->SetTitleOffset(0.95);
  hpthetabestmisid->GetYaxis()->SetLabelSize(0.05);
  hpthetabestmisid->GetYaxis()->SetTitleSize(0.05);
  hpthetabestmisid->GetYaxis()->SetTitleOffset(0.95);
  hpthetabestmisid->GetZaxis()->SetLabelSize(0.05);
  hpthetabestmisid->GetZaxis()->SetTitleSize(0.05);
  hpthetabestmisid->GetZaxis()->SetTitleOffset(0.95);
  // hpthetabestmisid->GetZaxis()->SetRangeUser(0.000001,1);
  // gPad->SetLogz();
  hpthetabestmisid->GetZaxis()->SetNdivisions(515);
  htextmisid->SetMarkerSize(1.);
  htextmisid2->SetMarkerSize(1.);
  htextmisid2->SetMarkerColor(0);
  hpthetabestmisid->Draw("colz");
  htextmisid->SetBarOffset(0.21);
  htextmisid->Draw("TEXT SAME");
  htextmisid2->SetBarOffset(-1.*0.21);
  htextmisid2->Draw("TEXT SAME");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetabestmisid->GetListOfFunctions()->FindObject("palette");
  cpthetabestmisid->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetabestmisid->Write();
  cpthetabestmisid->Modified();
  // hpthetabestmisid->Write();
  cpthetabestmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestmisid_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetabestmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetabestmisid_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cpthetatruncmisid = new TCanvas("cpthetatruncmisid","cpthetatruncmisid",700,500);
  cpthetatruncmisid->cd();
  hpthetatruncmisid->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncmisid->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncmisid->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncmisid->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncmisid->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncmisid->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncmisid->GetZaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncmisid->GetZaxis()->SetNdivisions(515);
  hpthetatruncmisid->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncmisid->GetListOfFunctions()->FindObject("palette");
  cpthetatruncmisid->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetatruncmisid->Write();
  cpthetatruncmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncmisid_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetatruncmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncmisid_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cpthetatruncmisid2 = new TCanvas("cpthetatruncmisid2","cpthetatruncmisid2",700,500);
  cpthetatruncmisid2->cd();
  hpthetatruncmisid2->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncmisid2->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncmisid2->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid2->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncmisid2->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncmisid2->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid2->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncmisid2->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncmisid2->GetZaxis()->SetTitleOffset(0.95);
  hpthetatruncmisid2->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncmisid2->GetZaxis()->SetNdivisions(515);
  hpthetatruncmisid2->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncmisid2->GetListOfFunctions()->FindObject("palette");
  cpthetatruncmisid2->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  // hpthetatruncmisid2->Write();
  cpthetatruncmisid2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncmisid2_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cpthetatruncmisid2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cpthetatruncmisid2_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncresacorr = new TCanvas("ctruncresacorr","ctruncresacorr",700,500);
  ctruncresacorr->cd();
  htruncresacorr->Draw("colz");
  ctruncresacorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncrescorr_%s.eps",particle1.Data()),"eps");
  ctruncresacorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncrescorr_%s.root",particle1.Data()),"root");

  TCanvas *ctruncresbcorr = new TCanvas("ctruncresbcorr","ctruncresbcorr",700,500);
  ctruncresbcorr->cd();
  htruncresbcorr->Draw("colz");
  ctruncresbcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncrescorr_%s.eps",particle2.Data()),"eps");
  ctruncresbcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncrescorr_%s.root",particle2.Data()),"root");

  TCanvas *ctruncsepcorr = new TCanvas("ctruncsepcorr","ctruncsepcorr",700,500);
  ctruncsepcorr->cd();
  htruncsepcorr->Draw("colz");
  ctruncsepcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsepcorr_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncsepcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsepcorr_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncmisidcorr = new TCanvas("ctruncmisidcorr","ctruncmisidcorr",700,500);
  ctruncmisidcorr->cd();
  htruncmisidcorr->Draw("colz");
  ctruncmisidcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisidcorr_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncmisidcorr->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisidcorr_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncresa = new TCanvas("ctruncresa","ctruncresa",700,500);
  ctruncresa->cd();
  htruncresa->GetXaxis()->SetNdivisions(505);
  htruncresa->GetXaxis()->SetLabelSize(0.05);
  htruncresa->GetXaxis()->SetTitleSize(0.05);
  htruncresa->GetXaxis()->SetTitleOffset(0.95);
  htruncresa->GetYaxis()->SetLabelSize(0.05);
  htruncresa->GetYaxis()->SetTitleSize(0.05);
  htruncresa->GetYaxis()->SetTitleOffset(0.95);
  htruncresa->GetXaxis()->SetNdivisions(515);
  htruncresa->SetFillColor(kBlue);
  // htruncresa->Fit("gaus");
  htruncresa->Scale(1./sumbins);
  htruncresa->Draw();
  ctruncresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres_%s.eps",particle1.Data()),"eps");
  ctruncresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres_%s.root",particle1.Data()),"root");

  TCanvas *ctruncresa2 = new TCanvas("ctruncresa2","ctruncresa2",700,500);
  ctruncresa2->cd();
  htruncresa2->GetXaxis()->SetNdivisions(505);
  htruncresa2->GetXaxis()->SetLabelSize(0.05);
  htruncresa2->GetXaxis()->SetTitleSize(0.05);
  htruncresa2->GetXaxis()->SetTitleOffset(0.95);
  htruncresa2->GetYaxis()->SetLabelSize(0.05);
  htruncresa2->GetYaxis()->SetTitleSize(0.05);
  htruncresa2->GetYaxis()->SetTitleOffset(0.95);
  htruncresa2->GetXaxis()->SetNdivisions(515);
  htruncresa2->SetFillColor(kBlue);
  // htruncresa2->Fit("gaus");
  htruncresa2->Scale(1./sumbins);
  htruncresa2->Draw();
  ctruncresa2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres2_%s.eps",particle1.Data()),"eps");
  ctruncresa2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres2_%s.root",particle1.Data()),"root");

  TCanvas *ctruncresb = new TCanvas("ctruncresb","ctruncresb",700,500);
  ctruncresb->cd();
  htruncresb->GetXaxis()->SetNdivisions(505);
  htruncresb->GetXaxis()->SetLabelSize(0.05);
  htruncresb->GetXaxis()->SetTitleSize(0.05);
  htruncresb->GetXaxis()->SetTitleOffset(0.95);
  htruncresb->GetYaxis()->SetLabelSize(0.05);
  htruncresb->GetYaxis()->SetTitleSize(0.05);
  htruncresb->GetYaxis()->SetTitleOffset(0.95);
  htruncresb->GetXaxis()->SetNdivisions(515);
  htruncresb->SetFillColor(kBlue);
  // htruncresb->Fit("gaus");
  htruncresb->Scale(1./sumbins);
  htruncresb->Draw();
  ctruncresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres_%s.eps",particle2.Data()),"eps");
  ctruncresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres_%s.root",particle2.Data()),"root");

  TCanvas *ctruncresb2 = new TCanvas("ctruncresb2","ctruncresb2",700,500);
  ctruncresb2->cd();
  htruncresb2->GetXaxis()->SetNdivisions(505);
  htruncresb2->GetXaxis()->SetLabelSize(0.05);
  htruncresb2->GetXaxis()->SetTitleSize(0.05);
  htruncresb2->GetXaxis()->SetTitleOffset(0.95);
  htruncresb2->GetYaxis()->SetLabelSize(0.05);
  htruncresb2->GetYaxis()->SetTitleSize(0.05);
  htruncresb2->GetYaxis()->SetTitleOffset(0.95);
  htruncresb2->GetXaxis()->SetNdivisions(515);
  htruncresb2->SetFillColor(kBlue);
  // htruncresb2->Fit("gaus");
  htruncresb2->Scale(1./sumbins);
  htruncresb2->Draw();
  ctruncresb2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres2_%s.eps",particle2.Data()),"eps");
  ctruncresb2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncres2_%s.root",particle2.Data()),"root");

  TCanvas *ctruncsep = new TCanvas("ctruncsep","ctruncsep",700,500);
  ctruncsep->cd();
  htruncsep->GetXaxis()->SetLabelSize(0.05);
  htruncsep->GetXaxis()->SetTitleSize(0.05);
  htruncsep->GetXaxis()->SetTitleOffset(0.95);
  htruncsep->GetYaxis()->SetLabelSize(0.05);
  htruncsep->GetYaxis()->SetTitleSize(0.05);
  htruncsep->GetYaxis()->SetTitleOffset(0.95);
  htruncsep->GetXaxis()->SetNdivisions(515);
  // htruncsep->Fit("gaus");
  htruncsep->SetFillColor(kBlue);
  htruncsep->Scale(1./sumbins);
  htruncsep->Draw();
  ctruncsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsep_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncsep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsep_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncsep2 = new TCanvas("ctruncsep2","ctruncsep2",700,500);
  ctruncsep2->cd();
  htruncsep2->GetXaxis()->SetLabelSize(0.05);
  htruncsep2->GetXaxis()->SetTitleSize(0.05);
  htruncsep2->GetXaxis()->SetTitleOffset(0.95);
  htruncsep2->GetYaxis()->SetLabelSize(0.05);
  htruncsep2->GetYaxis()->SetTitleSize(0.05);
  htruncsep2->GetYaxis()->SetTitleOffset(0.95);
  htruncsep2->GetXaxis()->SetNdivisions(515);
  // htruncsep2->Fit("gaus");
  htruncsep2->SetFillColor(kBlue);
  htruncsep2->Scale(1./sumbins);
  htruncsep2->Draw();
  ctruncsep2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsep2_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncsep2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncsep2_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncmisid = new TCanvas("ctruncmisid","ctruncmisid",700,500);
  ctruncmisid->cd();
  htruncmisid->GetXaxis()->SetLabelSize(0.05);
  htruncmisid->GetXaxis()->SetTitleSize(0.05);
  htruncmisid->GetXaxis()->SetTitleOffset(0.95);
  htruncmisid->GetYaxis()->SetLabelSize(0.05);
  htruncmisid->GetYaxis()->SetTitleSize(0.05);
  htruncmisid->GetYaxis()->SetTitleOffset(0.95);
  htruncmisid->GetXaxis()->SetNdivisions(515);
  // htruncmisid->Fit("gaus");
  htruncmisid->SetFillColor(kBlue);
  htruncmisid->Scale(1./sumbins);
  htruncmisid->Draw();
  ctruncmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisid_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisid_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *ctruncmisid2 = new TCanvas("ctruncmisid2","ctruncmisid2",700,500);
  ctruncmisid2->cd();
  htruncmisid2->GetXaxis()->SetLabelSize(0.05);
  htruncmisid2->GetXaxis()->SetTitleSize(0.05);
  htruncmisid2->GetXaxis()->SetTitleOffset(0.95);
  htruncmisid2->GetYaxis()->SetLabelSize(0.05);
  htruncmisid2->GetYaxis()->SetTitleSize(0.05);
  htruncmisid2->GetYaxis()->SetTitleOffset(0.95);
  htruncmisid2->GetXaxis()->SetNdivisions(515);
  // htruncmisid2->Fit("gaus");
  htruncmisid2->SetFillColor(kBlue);
  htruncmisid2->Scale(1./sumbins);
  htruncmisid2->Draw();
  ctruncmisid2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisid2_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  ctruncmisid2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/ctruncmisid2_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgmeanap = new TCanvas("cmgmeanap","cmgmeanap",700,500);
  TLegend *lmgmeanap = new TLegend(0.60,0.80,0.95,0.99);
  lmgmeanap->SetTextSize(0.04);
  lmgmeanap->SetBorderSize(0);
  cmgmeanap->cd();
  cmgmeanap->SetGrid();
  ergrmeanap0->SetMarkerColor(1);
  ergrmeanap0->SetMarkerSize(1.5);
  ergrmeanap0->SetLineColor(1);
  ergrmeanap0->SetMarkerStyle(20);
  mgmeanap->Add(ergrmeanap0);
  ergrmeanap5->SetMarkerColor(kBlue);
  ergrmeanap5->SetMarkerSize(1.5);
  ergrmeanap5->SetLineColor(kBlue);
  ergrmeanap5->SetMarkerStyle(20);
  mgmeanap->Add(ergrmeanap5);
  ergrmeanap20->SetMarkerColor(kRed);
  ergrmeanap20->SetMarkerSize(1.5);
  ergrmeanap20->SetLineColor(kRed);
  ergrmeanap20->SetMarkerStyle(20);
  mgmeanap->Add(ergrmeanap20);
  ergrmeanap50->SetMarkerColor(kMagenta);
  ergrmeanap50->SetMarkerSize(1.5);
  ergrmeanap50->SetLineColor(kMagenta);
  ergrmeanap50->SetMarkerStyle(20);
  mgmeanap->Add(ergrmeanap50);
  ergrmeanap455->SetMarkerColor(28);
  ergrmeanap455->SetMarkerSize(1.5);
  ergrmeanap455->SetMarkerStyle(20);
  mgmeanap->Add(ergrmeanap455);
  mgmeanap->SetTitle("mean Energy loss in CDC vs. momentum (Protons); p [GeV/c]; #mu_{<dE/dx>} [keV/cm]");
  mgmeanap->Draw("ALP");
  mgmeanap->SetMinimum(0.);
  cmgmeanap->Modified();
  lmgmeanap->AddEntry(ergrmeanap0, "truncation (0%,0%)","p");
  lmgmeanap->AddEntry(ergrmeanap5, "truncation (0%,5%)","p");
  lmgmeanap->AddEntry(ergrmeanap20, "truncation (0%,20%)","p");
  lmgmeanap->AddEntry(ergrmeanap50, "truncation (0%,50%)","p");
  lmgmeanap->AddEntry(ergrmeanap455, "truncation (45%,5%)","p");
  lmgmeanap->Draw();
  cmgmeanap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgmeanp_%s.eps",particle1.Data()),"eps");
  cmgmeanap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgmeanp_%s.root",particle1.Data()),"root");

  TCanvas *cmgmeanbp = new TCanvas("cmgmeanbp","cmgmeanbp",700,500);
  TLegend *lmgmeanbp = new TLegend(0.60,0.80,0.95,0.99);
  lmgmeanbp->SetTextSize(0.04);
  lmgmeanbp->SetBorderSize(0);
  cmgmeanbp->cd();
  cmgmeanbp->SetGrid();
  ergrmeanbp0->SetMarkerColor(1);
  ergrmeanbp0->SetMarkerSize(1.5);
  ergrmeanbp0->SetLineColor(1);
  ergrmeanbp0->SetMarkerStyle(20);
  mgmeanbp->Add(ergrmeanbp0);
  ergrmeanbp5->SetMarkerColor(kBlue);
  ergrmeanbp5->SetMarkerSize(1.5);
  ergrmeanbp5->SetLineColor(kBlue);
  ergrmeanbp5->SetMarkerStyle(20);
  mgmeanbp->Add(ergrmeanbp5);
  ergrmeanbp20->SetMarkerColor(kRed);
  ergrmeanbp20->SetMarkerSize(1.5);
  ergrmeanbp20->SetLineColor(kRed);
  ergrmeanbp20->SetMarkerStyle(20);
  mgmeanbp->Add(ergrmeanbp20);
  ergrmeanbp50->SetMarkerColor(kMagenta);
  ergrmeanbp50->SetMarkerSize(1.5);
  ergrmeanbp50->SetLineColor(kMagenta);
  ergrmeanbp50->SetMarkerStyle(20);
  mgmeanbp->Add(ergrmeanbp50);
  ergrmeanbp455->SetMarkerColor(28);
  ergrmeanbp455->SetMarkerSize(1.5);
  ergrmeanbp455->SetMarkerStyle(20);
  mgmeanbp->Add(ergrmeanbp455);
  mgmeanbp->SetTitle("mean Energy loss in CDC vs. momentum (#pi^{+}); p [GeV/c]; #mu_{<dE/dx>} [keV/cm]");
  mgmeanbp->Draw("ALP");
  mgmeanbp->SetMinimum(0.);
  cmgmeanbp->Modified();
  lmgmeanbp->AddEntry(ergrmeanbp0, "truncation (0%,0%)","p");
  lmgmeanbp->AddEntry(ergrmeanbp5, "truncation (0%,5%)","p");
  lmgmeanbp->AddEntry(ergrmeanbp20, "truncation (0%,20%)","p");
  lmgmeanbp->AddEntry(ergrmeanbp50, "truncation (0%,50%)","p");
  lmgmeanbp->AddEntry(ergrmeanbp455, "truncation (45%,5%)","p");
  lmgmeanbp->Draw();
  cmgmeanbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgmeanp_%s.eps",particle2.Data()),"eps");
  cmgmeanbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgmeanp_%s.root",particle2.Data()),"root");

  TCanvas *cmgergrresap = new TCanvas("cmgergrresap","cmgergrresap",700,500);
  TLegend *lmgergrresap = new TLegend(0.60,0.75,0.95,0.99);
  lmgergrresap->SetTextSize(0.04);
  lmgergrresap->SetBorderSize(0);
  cmgergrresap->cd();
  cmgergrresap->SetGrid();
  ergrresap0->SetMarkerColor(1);
  ergrresap0->SetMarkerSize(1.5);
  ergrresap0->SetLineColor(1);
  ergrresap0->SetMarkerStyle(20);
  ergrresap0->Draw("ALP");
  mgergrresap->Add(ergrresap0);
  ergrresap5->SetMarkerColor(kBlue);
  ergrresap5->SetMarkerSize(1.5);
  ergrresap5->SetLineColor(kBlue);
  ergrresap5->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap5);
  ergrresap20->SetMarkerColor(kRed);
  ergrresap20->SetMarkerSize(1.5);
  ergrresap20->SetLineColor(kRed);
  ergrresap20->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap20);
  ergrresap50->SetMarkerColor(kMagenta);
  ergrresap50->SetMarkerSize(1.5);
  ergrresap50->SetLineColor(kMagenta);
  ergrresap50->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap50);
  ergrresap455->SetMarkerColor(28);
  ergrresap455->SetMarkerSize(1.5);
  ergrresap455->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap455);
  ergrresap5030->SetMarkerColor(7);
  ergrresap5030->SetLineColor(7);
  ergrresap5030->SetFillColor(7);
  ergrresap5030->SetFillStyle(3004);
  ergrresap5030->SetMarkerSize(1.5);
  ergrresap5030->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap5030);
  ergrresap6510->SetMarkerColor(8);
  ergrresap6510->SetLineColor(8);
  ergrresap6510->SetFillColor(8);
  ergrresap6510->SetFillStyle(3004);
  ergrresap6510->SetMarkerSize(1.5);
  ergrresap6510->SetMarkerStyle(20);
  mgergrresap->Add(ergrresap6510);
  mgergrresap->SetTitle("resoluion of Energy loss in the CDC vs. momentum (protons); p [GeV/c]; #sigma_{<dE/dx>} [keV/cm]");
  mgergrresap->Draw("ALP");
  mgergrresap->SetMinimum(0.);
  cmgergrresap->Modified();
  lmgergrresap->AddEntry(ergrresap0, "truncation (0%,0%)","p");
  lmgergrresap->AddEntry(ergrresap5, "truncation (0%,5%)","p");
  lmgergrresap->AddEntry(ergrresap20, "truncation (0%,20%)","p");
  lmgergrresap->AddEntry(ergrresap50, "truncation (0%,50%)","p");
  lmgergrresap->AddEntry(ergrresap455, "truncation (45%,5%)","p");
  lmgergrresap->AddEntry(ergrresap5030, "truncation (50%,30%)","p");
  lmgergrresap->AddEntry(ergrresap6510, "truncation (65%,10%)","p");
  lmgergrresap->Draw();
  cmgergrresap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrresp_%s.eps",particle1.Data()),"eps");
  cmgergrresap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrresp_%s.root",particle1.Data()),"root");

  TCanvas *cmgergrresbp = new TCanvas("cmgergrresbp","cmgergrresbp",700,500);
  TLegend *lmgergrresbp = new TLegend(0.60,0.75,0.95,0.99);
  lmgergrresbp->SetTextSize(0.04);
  lmgergrresbp->SetBorderSize(0);
  cmgergrresbp->cd();
  cmgergrresbp->SetGrid();
  ergrresbp0->SetMarkerColor(1);
  ergrresbp0->SetMarkerSize(1.5);
  ergrresbp0->SetLineColor(1);
  ergrresbp0->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp0);
  ergrresbp5->SetMarkerColor(kBlue);
  ergrresbp5->SetMarkerSize(1.5);
  ergrresbp5->SetLineColor(kBlue);
  ergrresbp5->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp5);
  ergrresbp20->SetMarkerColor(kRed);
  ergrresbp20->SetMarkerSize(1.5);
  ergrresbp20->SetLineColor(kRed);
  ergrresbp20->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp20);
  ergrresbp50->SetMarkerColor(kMagenta);
  ergrresbp50->SetMarkerSize(1.5);
  ergrresbp50->SetLineColor(kMagenta);
  ergrresbp50->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp50);
  ergrresbp455->SetMarkerColor(28);
  ergrresbp455->SetMarkerSize(1.5);
  ergrresbp455->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp455);
  ergrresbp5030->SetMarkerColor(7);
  ergrresbp5030->SetLineColor(7);
  ergrresbp5030->SetFillColor(7);
  ergrresbp5030->SetFillStyle(3004);
  ergrresbp5030->SetMarkerSize(1.5);
  ergrresbp5030->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp5030);
  ergrresbp6510->SetMarkerColor(8);
  ergrresbp6510->SetLineColor(8);
  ergrresbp6510->SetFillColor(8);
  ergrresbp6510->SetFillStyle(3004);
  ergrresbp6510->SetMarkerSize(1.5);
  ergrresbp6510->SetMarkerStyle(20);
  mgergrresbp->Add(ergrresbp6510);
  mgergrresbp->SetTitle("resoluion of Energy loss in the CDC vs. momentum (#pi^{+}); p [GeV/c]; #sigma_{<dE/dx>} [keV/cm]");
  mgergrresbp->Draw("ALP");
  mgergrresbp->SetMinimum(0.);
  cmgergrresbp->Modified();
  lmgergrresbp->AddEntry(ergrresbp0, "truncation (0%,0%)","p");
  lmgergrresbp->AddEntry(ergrresbp5, "truncation (0%,5%)","p");
  lmgergrresbp->AddEntry(ergrresbp20, "truncation (0%,20%)","p");
  lmgergrresbp->AddEntry(ergrresbp50, "truncation (0%,50%)","p");
  lmgergrresbp->AddEntry(ergrresbp455, "truncation (45%,5%)","p");
  lmgergrresbp->AddEntry(ergrresbp5030, "truncation (50%,30%)","p");
  lmgergrresbp->AddEntry(ergrresbp6510, "truncation (65%,10%)","p");
  lmgergrresbp->Draw();
  cmgergrresbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrresp_%s.eps",particle2.Data()),"eps");
  cmgergrresbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrresp_%s.root",particle2.Data()),"root");

  TCanvas *cmgergrsepp = new TCanvas("cmgergrsepp","cmgergrsepp",700,500);
  TLegend *lmgergrsepp = new TLegend(0.60,0.75,0.95,0.99);
  lmgergrsepp->SetTextSize(0.04);
  lmgergrsepp->SetBorderSize(0);
  cmgergrsepp->cd();
  cmgergrsepp->SetGrid();
  ergrsepp0->SetMarkerColor(1);
  ergrsepp0->SetMarkerSize(1.5);
  ergrsepp0->SetLineColor(1);
  ergrsepp0->SetMarkerStyle(20);
  // mgergrsepp->Add(ergrsepp0);
  ergrsepp5->SetMarkerColor(kBlue);
  ergrsepp5->SetLineColor(kBlue);
  ergrsepp5->SetFillColor(kBlue);
  ergrsepp5->SetFillStyle(3004);
  ergrsepp5->SetMarkerSize(1.5);
  ergrsepp5->SetMarkerStyle(20);
  mgergrsepp->Add(ergrsepp5);
  ergrsepp20->SetMarkerColor(kRed);
  ergrsepp20->SetLineColor(kRed);
  ergrsepp20->SetFillColor(kRed);
  ergrsepp20->SetFillStyle(3004);
  ergrsepp20->SetMarkerSize(1.5);
  ergrsepp20->SetMarkerStyle(20);
  mgergrsepp->Add(ergrsepp20);
  ergrsepp50->SetMarkerColor(kMagenta);
  ergrsepp50->SetLineColor(kMagenta);
  ergrsepp50->SetFillColor(kMagenta);
  ergrsepp50->SetFillStyle(3004);
  ergrsepp50->SetMarkerSize(1.5);
  ergrsepp50->SetMarkerStyle(20);
  mgergrsepp->Add(ergrsepp50);
  ergrsepp455->SetMarkerColor(28);
  ergrsepp455->SetLineColor(28);
  ergrsepp455->SetFillColor(28);
  ergrsepp455->SetFillStyle(3004);
  ergrsepp455->SetMarkerSize(1.5);
  ergrsepp455->SetMarkerStyle(20);
  // mgergrsepp->Add(ergrsepp455);
  // ergrsepp5030->SetMarkerColor(7);
  // ergrsepp5030->SetLineColor(7);
  // ergrsepp5030->SetFillColor(7);
  // ergrsepp5030->SetFillStyle(3004);
  // ergrsepp5030->SetMarkerSize(1.5);
  // ergrsepp5030->SetMarkerStyle(20);
  // mgergrsepp->Add(ergrsepp5030);
  // ergrsepp6510->SetMarkerColor(8);
  // ergrsepp6510->SetLineColor(8);
  // ergrsepp6510->SetFillColor(8);
  // ergrsepp6510->SetFillStyle(3004);
  // ergrsepp6510->SetMarkerSize(1.5);
  // ergrsepp6510->SetMarkerStyle(20);
  // mgergrsepp->Add(ergrsepp6510);
  mgergrsepp->SetTitle("separation power vs. momentum; p [GeV/c]; separation power");
  mgergrsepp->Draw("ALP");
  mgergrsepp->SetMinimum(0.);
  cmgergrsepp->Modified();
  // lmgergrsepp->AddEntry(ergrsepp0, "truncation (0%,0%)","p");
  lmgergrsepp->AddEntry(ergrsepp5, "truncation (0%,5%)","p");
  lmgergrsepp->AddEntry(ergrsepp20, "truncation (0%,20%)","p");
  lmgergrsepp->AddEntry(ergrsepp50, "truncation (0%,50%)","p");
  // lmgergrsepp->AddEntry(ergrsepp455, "truncation (45%,5%)","p");
  // lmgergrsepp->AddEntry(ergrsepp5030, "truncation (50%,30%)","p");
  // lmgergrsepp->AddEntry(ergrsepp6510, "truncation (65%,10%)","p");
  lmgergrsepp->Draw();
  cmgergrsepp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrsepp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgergrsepp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrsepp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgergrmisidp = new TCanvas("cmgergrmisidp","cmgergrmisidp",700,500);
  TLegend *lmgergrmisidp = new TLegend(0.60,0.75,0.95,0.99);
  lmgergrmisidp->SetTextSize(0.04);
  lmgergrmisidp->SetBorderSize(0);
  cmgergrmisidp->cd();
  cmgergrmisidp->SetGrid();
  ergrmisidp0->SetMarkerColor(1);
  ergrmisidp0->SetMarkerSize(1.5);
  ergrmisidp0->SetLineColor(1);
  ergrmisidp0->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp0);
  ergrmisidp5->SetMarkerColor(kBlue);
  ergrmisidp5->SetMarkerSize(1.5);
  ergrmisidp5->SetLineColor(kBlue);
  ergrmisidp5->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp5);
  ergrmisidp20->SetMarkerColor(kRed);
  ergrmisidp20->SetMarkerSize(1.5);
  ergrmisidp20->SetLineColor(kRed);
  ergrmisidp20->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp20);
  ergrmisidp50->SetMarkerColor(kMagenta);
  ergrmisidp50->SetMarkerSize(1.5);
  ergrmisidp50->SetLineColor(kMagenta);
  ergrmisidp50->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp50);
  ergrmisidp455->SetMarkerColor(28);
  ergrmisidp455->SetMarkerSize(1.5);
  ergrmisidp455->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp455);
  ergrmisidp5030->SetMarkerColor(7);
  ergrmisidp5030->SetLineColor(7);
  ergrmisidp5030->SetFillColor(7);
  ergrmisidp5030->SetFillStyle(3004);
  ergrmisidp5030->SetMarkerSize(1.5);
  ergrmisidp5030->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp5030);
  ergrmisidp6510->SetMarkerColor(8);
  ergrmisidp6510->SetLineColor(8);
  ergrmisidp6510->SetFillColor(8);
  ergrmisidp6510->SetFillStyle(3004);
  ergrmisidp6510->SetMarkerSize(1.5);
  ergrmisidp6510->SetMarkerStyle(20);
  mgergrmisidp->Add(ergrmisidp6510);
  mgergrmisidp->SetTitle("mis-ID vs. momentum; p [GeV/c]; mis-ID");
  mgergrmisidp->Draw("ALP");
  mgergrmisidp->SetMinimum(0.);
  cmgergrmisidp->Modified();
  lmgergrmisidp->AddEntry(ergrmisidp0, "truncation (0%,0%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp5, "truncation (0%,5%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp20, "truncation (0%,20%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp50, "truncation (0%,50%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp455, "truncation (45%,5%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp5030, "truncation (50%,30%)","p");
  lmgergrmisidp->AddEntry(ergrmisidp6510, "truncation (65%,10%)","p");
  lmgergrmisidp->Draw();
  cmgergrmisidp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrmisidp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgergrmisidp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgergrmisidp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgnhitsap = new TCanvas("cmgnhitsap","cmgnhitsap",700,500);
  TLegend *lmgnhitsap = new TLegend(0.60,0.75,0.95,0.99);
  lmgnhitsap->SetTextSize(0.04);
  lmgnhitsap->SetBorderSize(0);
  cmgnhitsap->cd();
  cmgnhitsap->SetGrid();
  ergrnhitsap0->SetMarkerColor(1);
  ergrnhitsap0->SetMarkerSize(1.5);
  ergrnhitsap0->SetLineColor(1);
  ergrnhitsap0->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap0);
  ergrnhitsap5->SetMarkerColor(kBlue);
  ergrnhitsap5->SetMarkerSize(1.5);
  ergrnhitsap5->SetLineColor(kBlue);
  ergrnhitsap5->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap5);
  ergrnhitsap20->SetMarkerColor(kRed);
  ergrnhitsap20->SetMarkerSize(1.5);
  ergrnhitsap20->SetLineColor(kRed);
  ergrnhitsap20->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap20);
  ergrnhitsap50->SetMarkerColor(kMagenta);
  ergrnhitsap50->SetMarkerSize(1.5);
  ergrnhitsap50->SetLineColor(kMagenta);
  ergrnhitsap50->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap50);
  ergrnhitsap455->SetMarkerColor(28);
  ergrnhitsap455->SetMarkerSize(1.5);
  ergrnhitsap455->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap455);
  ergrnhitsap5030->SetMarkerColor(7);
  ergrnhitsap5030->SetLineColor(7);
  ergrnhitsap5030->SetFillColor(7);
  ergrnhitsap5030->SetFillStyle(3004);
  ergrnhitsap5030->SetMarkerSize(1.5);
  ergrnhitsap5030->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap5030);
  ergrnhitsap6510->SetMarkerColor(8);
  ergrnhitsap6510->SetLineColor(8);
  ergrnhitsap6510->SetFillColor(8);
  ergrnhitsap6510->SetFillStyle(3004);
  ergrnhitsap6510->SetMarkerSize(1.5);
  ergrnhitsap6510->SetMarkerStyle(20);
  mgnhitsap->Add(ergrnhitsap6510);
  mgnhitsap->SetTitle("# hits vs. momentum (Protons); p [GeV/c]; # nhits");
  mgnhitsap->Draw("ALP");
  mgnhitsap->SetMinimum(0.);
  cmgnhitsap->Modified();
  lmgnhitsap->AddEntry(ergrnhitsap0, "truncation (0%,0%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap5, "truncation (0%,5%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap20, "truncation (0%,20%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap50, "truncation (0%,50%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap455, "truncation (45%,5%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap5030, "truncation (50%,30%)","p");
  lmgnhitsap->AddEntry(ergrnhitsap6510, "truncation (65%,10%)","p");
  lmgnhitsap->Draw();
  cmgnhitsap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnhitsp_%s.eps",particle1.Data()),"eps");
  cmgnhitsap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnhitsp_%s.root",particle1.Data()),"root");

  TCanvas *cmgnhitsbp = new TCanvas("cmgnhitsbp","cmgnhitsbp",700,500);
  TLegend *lmgnhitsbp = new TLegend(0.60,0.75,0.95,0.99);
  lmgnhitsbp->SetTextSize(0.04);
  lmgnhitsbp->SetBorderSize(0);
  cmgnhitsbp->cd();
  cmgnhitsbp->SetGrid();
  ergrnhitsbp0->SetMarkerColor(1);
  ergrnhitsbp0->SetMarkerSize(1.5);
  ergrnhitsbp0->SetLineColor(1);
  ergrnhitsbp0->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp0);
  ergrnhitsbp5->SetMarkerColor(kBlue);
  ergrnhitsbp5->SetMarkerSize(1.5);
  ergrnhitsbp5->SetLineColor(kBlue);
  ergrnhitsbp5->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp5);
  ergrnhitsbp20->SetMarkerColor(kRed);
  ergrnhitsbp20->SetMarkerSize(1.5);
  ergrnhitsbp20->SetLineColor(kRed);
  ergrnhitsbp20->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp20);
  ergrnhitsbp50->SetMarkerColor(kMagenta);
  ergrnhitsbp50->SetMarkerSize(1.5);
  ergrnhitsbp50->SetLineColor(kMagenta);
  ergrnhitsbp50->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp50);
  ergrnhitsbp455->SetMarkerColor(28);
  ergrnhitsbp455->SetMarkerSize(1.5);
  ergrnhitsbp455->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp455);
  ergrnhitsbp5030->SetMarkerColor(7);
  ergrnhitsbp5030->SetLineColor(7);
  ergrnhitsbp5030->SetFillColor(7);
  ergrnhitsbp5030->SetFillStyle(3004);
  ergrnhitsbp5030->SetMarkerSize(1.5);
  ergrnhitsbp5030->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp5030);
  ergrnhitsbp6510->SetMarkerColor(8);
  ergrnhitsbp6510->SetLineColor(8);
  ergrnhitsbp6510->SetFillColor(8);
  ergrnhitsbp6510->SetFillStyle(3004);
  ergrnhitsbp6510->SetMarkerSize(1.5);
  ergrnhitsbp6510->SetMarkerStyle(20);
  mgnhitsbp->Add(ergrnhitsbp6510);
  mgnhitsbp->SetTitle("# hits vs. momentum (#pi^{-}); p [GeV/c]; # nhits");
  mgnhitsbp->Draw("ALP");
  mgnhitsbp->SetMinimum(0.);
  cmgnhitsbp->Modified();
  lmgnhitsbp->AddEntry(ergrnhitsbp0, "truncation (0%,0%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp5, "truncation (0%,5%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp20, "truncation (0%,20%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp50, "truncation (0%,50%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp455, "truncation (45%,5%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp5030, "truncation (50%,30%)","p");
  lmgnhitsbp->AddEntry(ergrnhitsbp6510, "truncation (65%,10%)","p");
  lmgnhitsbp->Draw();
  cmgnhitsbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnhitsp_%s.eps",particle2.Data()),"eps");
  cmgnhitsbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnhitsp_%s.root",particle2.Data()),"root");

  TCanvas *cmgnotresap = new TCanvas("cmgnotresap","cmgnotresap",700,500);
  TLegend *lmgnotresap = new TLegend(0.60,0.75,0.95,0.99);
  lmgnotresap->SetTextSize(0.04);
  lmgnotresap->SetBorderSize(0);
  cmgnotresap->cd();
  cmgnotresap->SetGrid();
  grnotresap0->SetMarkerColor(1);
  grnotresap0->SetMarkerSize(1.5);
  grnotresap0->SetLineColor(1);
  grnotresap0->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap0);
  grnotresap5->SetMarkerColor(kBlue);
  grnotresap5->SetLineColor(kBlue);
  grnotresap5->SetFillColor(kBlue);
  grnotresap5->SetFillStyle(3004);
  grnotresap5->SetMarkerSize(1.5);
  grnotresap5->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap5);
  grnotresap20->SetMarkerColor(kRed);
  grnotresap20->SetLineColor(kRed);
  grnotresap20->SetFillColor(kRed);
  grnotresap20->SetFillStyle(3004);
  grnotresap20->SetMarkerSize(1.5);
  grnotresap20->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap20);
  grnotresap50->SetMarkerColor(kMagenta);
  grnotresap50->SetLineColor(kMagenta);
  grnotresap50->SetFillColor(kMagenta);
  grnotresap50->SetFillStyle(3004);
  grnotresap50->SetMarkerSize(1.5);
  grnotresap50->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap50);
  grnotresap455->SetMarkerColor(28);
  grnotresap455->SetLineColor(28);
  grnotresap455->SetFillColor(28);
  grnotresap455->SetFillStyle(3004);
  grnotresap455->SetMarkerSize(1.5);
  grnotresap455->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap455);
  grnotresap5030->SetMarkerColor(7);
  grnotresap5030->SetLineColor(7);
  grnotresap5030->SetFillColor(7);
  grnotresap5030->SetFillStyle(3004);
  grnotresap5030->SetMarkerSize(1.5);
  grnotresap5030->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap5030);
  grnotresap6510->SetMarkerColor(8);
  grnotresap6510->SetLineColor(8);
  grnotresap6510->SetFillColor(8);
  grnotresap6510->SetFillStyle(3004);
  grnotresap6510->SetMarkerSize(1.5);
  grnotresap6510->SetMarkerStyle(20);
  mgnotresap->Add(grnotresap6510);
  mgnotresap->SetTitle("Number of bins in #theta (#sigma_{dE/dx(Proton)} = 0) vs. momentum; p [GeV/c]; number of bins in #theta (#sigma_{dE/dx(Proton)} = 0)");
  mgnotresap->Draw("ALP");
  mgnotresap->SetMinimum(0.);
  cmgnotresap->Modified();
  lmgnotresap->AddEntry(grnotresap0, "truncation (0%,0%)","p");
  lmgnotresap->AddEntry(grnotresap5, "truncation (0%,5%)","p");
  lmgnotresap->AddEntry(grnotresap20, "truncation (0%,20%)","p");
  lmgnotresap->AddEntry(grnotresap50, "truncation (0%,50%)","p");
  lmgnotresap->AddEntry(grnotresap455, "truncation (45%,5%)","p");
  lmgnotresap->AddEntry(grnotresap5030, "truncation (50%,30%)","p");
  lmgnotresap->AddEntry(grnotresap6510, "truncation (65%,10%)","p");
  lmgnotresap->Draw();
  cmgnotresap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotresap_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgnotresap->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotresap_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgnotresbp = new TCanvas("cmgnotresbp","cmgnotresbp",700,500);
  TLegend *lmgnotresbp = new TLegend(0.60,0.75,0.95,0.99);
  lmgnotresbp->SetTextSize(0.04);
  lmgnotresbp->SetBorderSize(0);
  cmgnotresbp->cd();
  cmgnotresbp->SetGrid();
  grnotresbp0->SetMarkerColor(1);
  grnotresbp0->SetMarkerSize(1.5);
  grnotresbp0->SetLineColor(1);
  grnotresbp0->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp0);
  grnotresbp5->SetMarkerColor(kBlue);
  grnotresbp5->SetLineColor(kBlue);
  grnotresbp5->SetFillColor(kBlue);
  grnotresbp5->SetFillStyle(3004);
  grnotresbp5->SetMarkerSize(1.5);
  grnotresbp5->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp5);
  grnotresbp20->SetMarkerColor(kRed);
  grnotresbp20->SetLineColor(kRed);
  grnotresbp20->SetFillColor(kRed);
  grnotresbp20->SetFillStyle(3004);
  grnotresbp20->SetMarkerSize(1.5);
  grnotresbp20->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp20);
  grnotresbp50->SetMarkerColor(kMagenta);
  grnotresbp50->SetLineColor(kMagenta);
  grnotresbp50->SetFillColor(kMagenta);
  grnotresbp50->SetFillStyle(3004);
  grnotresbp50->SetMarkerSize(1.5);
  grnotresbp50->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp50);
  grnotresbp455->SetMarkerColor(28);
  grnotresbp455->SetLineColor(28);
  grnotresbp455->SetFillColor(28);
  grnotresbp455->SetFillStyle(3004);
  grnotresbp455->SetMarkerSize(1.5);
  grnotresbp455->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp455);
  grnotresbp5030->SetMarkerColor(7);
  grnotresbp5030->SetLineColor(7);
  grnotresbp5030->SetFillColor(7);
  grnotresbp5030->SetFillStyle(3004);
  grnotresbp5030->SetMarkerSize(1.5);
  grnotresbp5030->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp5030);
  grnotresbp6510->SetMarkerColor(8);
  grnotresbp6510->SetLineColor(8);
  grnotresbp6510->SetFillColor(8);
  grnotresbp6510->SetFillStyle(3004);
  grnotresbp6510->SetMarkerSize(1.5);
  grnotresbp6510->SetMarkerStyle(20);
  mgnotresbp->Add(grnotresbp6510);
  mgnotresbp->SetTitle("Number of bins in #theta (#sigma_{dE/dx(#pi^{+})} = 0) vs. momentum; p [GeV/c]; number of bins in #theta (#sigma_{dE/dx(#pi^{+})} = 0)");
  mgnotresbp->Draw("ALP");
  mgnotresbp->SetMinimum(0.);
  cmgnotresbp->Modified();
  lmgnotresbp->AddEntry(grnotresbp0, "truncation (0%,0%)","p");
  lmgnotresbp->AddEntry(grnotresbp5, "truncation (0%,5%)","p");
  lmgnotresbp->AddEntry(grnotresbp20, "truncation (0%,20%)","p");
  lmgnotresbp->AddEntry(grnotresbp50, "truncation (0%,50%)","p");
  lmgnotresbp->AddEntry(grnotresbp455, "truncation (45%,5%)","p");
  lmgnotresbp->AddEntry(grnotresbp5030, "truncation (50%,30%)","p");
  lmgnotresbp->AddEntry(grnotresbp6510, "truncation (65%,10%)","p");
  lmgnotresbp->Draw();
  cmgnotresbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotresbp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgnotresbp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotresbp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgnotsepp = new TCanvas("cmgnotsepp","cmgnotsepp",700,500);
  TLegend *lmgnotsepp = new TLegend(0.60,0.75,0.95,0.99);
  lmgnotsepp->SetTextSize(0.04);
  lmgnotsepp->SetBorderSize(0);
  cmgnotsepp->cd();
  cmgnotsepp->SetGrid();
  grnotsepp0->SetMarkerColor(1);
  grnotsepp0->SetMarkerSize(1.5);
  grnotsepp0->SetLineColor(1);
  grnotsepp0->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp0);
  grnotsepp5->SetMarkerColor(kBlue);
  grnotsepp5->SetLineColor(kBlue);
  grnotsepp5->SetFillColor(kBlue);
  grnotsepp5->SetFillStyle(3004);
  grnotsepp5->SetMarkerSize(1.5);
  grnotsepp5->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp5);
  grnotsepp20->SetMarkerColor(kRed);
  grnotsepp20->SetLineColor(kRed);
  grnotsepp20->SetFillColor(kRed);
  grnotsepp20->SetFillStyle(3004);
  grnotsepp20->SetMarkerSize(1.5);
  grnotsepp20->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp20);
  grnotsepp50->SetMarkerColor(kMagenta);
  grnotsepp50->SetLineColor(kMagenta);
  grnotsepp50->SetFillColor(kMagenta);
  grnotsepp50->SetFillStyle(3004);
  grnotsepp50->SetMarkerSize(1.5);
  grnotsepp50->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp50);
  grnotsepp455->SetMarkerColor(28);
  grnotsepp455->SetLineColor(28);
  grnotsepp455->SetFillColor(28);
  grnotsepp455->SetFillStyle(3004);
  grnotsepp455->SetMarkerSize(1.5);
  grnotsepp455->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp455);
  grnotsepp5030->SetMarkerColor(7);
  grnotsepp5030->SetLineColor(7);
  grnotsepp5030->SetFillColor(7);
  grnotsepp5030->SetFillStyle(3004);
  grnotsepp5030->SetMarkerSize(1.5);
  grnotsepp5030->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp5030);
  grnotsepp6510->SetMarkerColor(8);
  grnotsepp6510->SetLineColor(8);
  grnotsepp6510->SetFillColor(8);
  grnotsepp6510->SetFillStyle(3004);
  grnotsepp6510->SetMarkerSize(1.5);
  grnotsepp6510->SetMarkerStyle(20);
  mgnotsepp->Add(grnotsepp6510);
  mgnotsepp->SetTitle("Number of bins in #theta (separation power = 0) vs. momentum; p [GeV/c]; number of bins in #theta (separation power = 0)");
  mgnotsepp->Draw("ALP");
  mgnotsepp->SetMinimum(0.);
  cmgnotsepp->Modified();
  lmgnotsepp->AddEntry(grnotsepp0, "truncation (0%,0%)","p");
  lmgnotsepp->AddEntry(grnotsepp5, "truncation (0%,5%)","p");
  lmgnotsepp->AddEntry(grnotsepp20, "truncation (0%,20%)","p");
  lmgnotsepp->AddEntry(grnotsepp50, "truncation (0%,50%)","p");
  lmgnotsepp->AddEntry(grnotsepp455, "truncation (45%,5%)","p");
  lmgnotsepp->AddEntry(grnotsepp5030, "truncation (50%,30%)","p");
  lmgnotsepp->AddEntry(grnotsepp6510, "truncation (65%,10%)","p");
  lmgnotsepp->Draw();
  cmgnotsepp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotsepp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgnotsepp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotsepp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cmgnotmisidp = new TCanvas("cmgnotmisidp","cmgnotmisidp",700,500);
  TLegend *lmgnotmisidp = new TLegend(0.60,0.75,0.95,0.99);
  lmgnotmisidp->SetTextSize(0.04);
  lmgnotmisidp->SetBorderSize(0);
  cmgnotmisidp->cd();
  cmgnotmisidp->SetGrid();
  grnotmisidp0->SetMarkerColor(1);
  grnotmisidp0->SetMarkerSize(1.5);
  grnotmisidp0->SetLineColor(1);
  grnotmisidp0->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp0);
  grnotmisidp5->SetMarkerColor(kBlue);
  grnotmisidp5->SetLineColor(kBlue);
  grnotmisidp5->SetFillColor(kBlue);
  grnotmisidp5->SetFillStyle(3004);
  grnotmisidp5->SetMarkerSize(1.5);
  grnotmisidp5->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp5);
  grnotmisidp20->SetMarkerColor(kRed);
  grnotmisidp20->SetLineColor(kRed);
  grnotmisidp20->SetFillColor(kRed);
  grnotmisidp20->SetFillStyle(3004);
  grnotmisidp20->SetMarkerSize(1.5);
  grnotmisidp20->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp20);
  grnotmisidp50->SetMarkerColor(kMagenta);
  grnotmisidp50->SetLineColor(kMagenta);
  grnotmisidp50->SetFillColor(kMagenta);
  grnotmisidp50->SetFillStyle(3004);
  grnotmisidp50->SetMarkerSize(1.5);
  grnotmisidp50->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp50);
  grnotmisidp455->SetMarkerColor(28);
  grnotmisidp455->SetLineColor(28);
  grnotmisidp455->SetFillColor(28);
  grnotmisidp455->SetFillStyle(3004);
  grnotmisidp455->SetMarkerSize(1.5);
  grnotmisidp455->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp455);
  grnotmisidp5030->SetMarkerColor(7);
  grnotmisidp5030->SetLineColor(7);
  grnotmisidp5030->SetFillColor(7);
  grnotmisidp5030->SetFillStyle(3004);
  grnotmisidp5030->SetMarkerSize(1.5);
  grnotmisidp5030->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp5030);
  grnotmisidp6510->SetMarkerColor(8);
  grnotmisidp6510->SetLineColor(8);
  grnotmisidp6510->SetFillColor(8);
  grnotmisidp6510->SetFillStyle(3004);
  grnotmisidp6510->SetMarkerSize(1.5);
  grnotmisidp6510->SetMarkerStyle(20);
  mgnotmisidp->Add(grnotmisidp6510);
  mgnotmisidp->SetTitle("Number of bins in #theta (Mis-ID = 0) vs. momentum; p [GeV/c]; number of bins in #theta (Mis-ID = 0)");
  mgnotmisidp->Draw("ALP");
  mgnotmisidp->SetMinimum(0.);
  cmgnotmisidp->Modified();
  lmgnotmisidp->AddEntry(grnotmisidp0, "truncation (0%,0%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp5, "truncation (0%,5%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp20, "truncation (0%,20%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp50, "truncation (0%,50%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp455, "truncation (45%,5%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp5030, "truncation (50%,30%)","p");
  lmgnotmisidp->AddEntry(grnotmisidp6510, "truncation (65%,10%)","p");
  lmgnotmisidp->Draw();
  cmgnotmisidp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotmisidp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cmgnotmisidp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cmgnotmisidp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  TCanvas *cdedxsepcomp = new TCanvas("cdedxsepcomp","cdedxsepcomp",700,500);
  TLegend *ldedxsepcomp = new TLegend(0.50,0.65,0.95,0.95);
  ldedxsepcomp->SetTextSize(0.03);
  cdedxsepcomp->cd();
  TH1F *hdedxsep_protonflat_5_5_0_0 = (TH1F*) fpid3->Get("hdedxsep_protonflat_5_5_0_0");
  hdedxsep_protonflat_5_5_0_0->GetXaxis()->SetRangeUser(0,20);
  cout <<" ***** hdedxsep_protonflat_5_5_0_0 = "<<hdedxsep_protonflat_5_5_0_0<<endl;
  hdedxsep_protonflat_5_5_0_0->SetFillColor(8);
  hdedxsep_protonflat_5_5_0_0->SetLineWidth(0);
  hdedxsep_protonflat_5_5_0_0->Draw();
  ldedxsepcomp->AddEntry(hdedxsep_protonflat_5_5_0_0, "truncation (0%,0%), S = 5.237","f");
  TH1F *hdedxsep_protonflat_5_5_0_5 = (TH1F*) fpid3->Get("hdedxsep_protonflat_5_5_0_5");
  cout <<" ***** hdedxsep_protonflat_5_5_0_5 = "<<hdedxsep_protonflat_5_5_0_5<<endl;
  hdedxsep_protonflat_5_5_0_5->SetLineColor(6);
  hdedxsep_protonflat_5_5_0_5->SetLineWidth(2);
  hdedxsep_protonflat_5_5_0_5->Draw("same");
  ldedxsepcomp->AddEntry(hdedxsep_protonflat_5_5_0_5, "truncation (0%,5%), S = 5.344","l");
  TH1F *hdedxsep_protonflat_5_5_45_5 = (TH1F*) fpid3->Get("hdedxsep_protonflat_5_5_45_5");
  cout <<" ***** hdedxsep_protonflat_5_5_45_5 = "<<hdedxsep_protonflat_5_5_45_5<<endl;
  hdedxsep_protonflat_5_5_45_5->SetLineColor(28);
  hdedxsep_protonflat_5_5_45_5->SetLineWidth(2);
  hdedxsep_protonflat_5_5_45_5->Draw("same");
  ldedxsepcomp->AddEntry(hdedxsep_protonflat_5_5_45_5, "truncation (45%,5%), S = 6.450","l");
  TH1F *hdedxsep_pip_5_5_0_0 = (TH1F*) fpid3->Get("hdedxsep_pip_5_5_0_0");
  cout <<" ***** hdedxsep_pip_5_5_0_0 = "<<hdedxsep_pip_5_5_0_0<<endl;
  hdedxsep_pip_5_5_0_0->SetFillColor(8);
  hdedxsep_pip_5_5_0_0->SetLineWidth(0);
  hdedxsep_pip_5_5_0_0->Draw("same");
  TH1F *hdedxsep_pip_5_5_0_5 = (TH1F*) fpid3->Get("hdedxsep_pip_5_5_0_5");
  cout <<" ***** hdedxsep_pip_5_5_0_5 = "<<hdedxsep_pip_5_5_0_5<<endl;
  hdedxsep_pip_5_5_0_5->SetLineColor(6);
  hdedxsep_pip_5_5_0_5->SetLineWidth(2);
  hdedxsep_pip_5_5_0_5->Draw("same");
  TH1F *hdedxsep_pip_5_5_45_5 = (TH1F*) fpid3->Get("hdedxsep_pip_5_5_45_5");
  cout <<" ***** hdedxsep_pip_5_5_45_5 = "<<hdedxsep_pip_5_5_45_5<<endl;
  hdedxsep_pip_5_5_45_5->SetLineColor(28);
  hdedxsep_pip_5_5_45_5->SetLineWidth(2);
  hdedxsep_pip_5_5_45_5->Draw("same");
  ldedxsepcomp->Draw();
  cdedxsepcomp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cdedxsepcomp_%s_%s.eps",particle1.Data(),particle2.Data()),"eps");
  cdedxsepcomp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_final/cdedxsepcomp_%s_%s.root",particle1.Data(),particle2.Data()),"root");

  outputfig->Print();
  // outputfig->Close();
  foutput0.close();
  foutput20.close();
  foutput455.close();

}
