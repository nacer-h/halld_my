#define piddedx_cxx
// The class definition in piddedx.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("piddedx.C")
// Root > T->Process("piddedx.C","some options")
// Root > T->Process("piddedx.C+")
//

#include "piddedx.h"
#include <TH2.h>
#include <TStyle.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TLegend.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include <algorithm>
#include <TLine.h>
#include <TGraph.h>
#include <Riostream.h>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TPaletteAxis.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TProfile.h>
#include <TProfile2D.h>
#include <TPave.h>
#include <TBox.h>

bool dofit;
int trunclevel = 0;

TFile *outputfig;
// TH1F *hdedx;
TH2F *hdedxp;
TH2F *hddedxp;
TH2F *hdedxtheta;
TF1 *func_dedx_p1;
TF1 *func_dedx_p2;
TF1 *func_dEdxCut_SelectHeavy;
TH2F *hprot;
TH3F *h3;
TF1 *func_dedx_ex1;
TF1 *func_dedx_ex2;
TH2F *hptheta;
TH2F *hpthetaprot;
TH1F *hdedxres;
TF1 *fdedxsig;
TF1 *fdedxbkg;
TF1 *fdedx;
TH1 *hchi2ndf;

int ntheta, np;
double pmin, pmax, pstep, pcut, thetamin, thetamax, thetastep, thetacut;
TH1F *hdedxav[20][20];
double res, sep;
TH2F *hpthetares;
TH2F *hpthetasep;

double bethblock(double *x, double *par)
{
  double p = x[0];
  double m = par[0];
  double A = par[1];
  double I = par[2]; // Mean excitation of the medium.
  double F = 1; // = z²*(Z/A) ,z: charge of incident particle, Z and A: charge number and atomic mass of the medium.
  double D = 0; // Density correction (depend on B and G).
  double k = 0.3071;  // MeV.cm².g¯1
  double me = 0.00051099892; // GeV , electron mass.
  double c = 29979245800; // cm.s¯1 , celerity of light in vacuum.
  double B = p/sqrt(m*m+p*p); // Velocity.
  double G = 1./sqrt(1-B*B); // Lorentz factor.
  // double T = (2*me*B*B*G*G)/(1+2*G*me/m+(me*me/m*m)); // Maximum energy transfer in single collision.

  return A*k/(B*B)*(0.5*log(B*B*G*G*(2.*me/I))-(B*B));
  //  [0]/(x*x/(x*x+[1]*[1]))*(0.5*log([2]*x*x/[1]*[1])-(x*x/(x*x+[1]*[1])))
  //  (A*(-k)*F*(1./(B*B)))*(0.5*log((2.*me*B*B*G*G*T)/I)-(B*B)-(D/2.));
}

double mylan(double *x, double *par)
{
  double lambda = (x[0]-par[2])/par[1];
  return par[0]*exp(-0.5*(lambda+exp(-1*lambda)));
}

void piddedx::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  dofit =false;
  if (option!="fit")
  trunclevel=option.Atoi();
  else
  {
    trunclevel = 0;
    dofit = true;
  }

  // hdedx = new TH1F("hdedx","",100,0,2);

  outputfig = new TFile("fig/pid.root","UPDATE");  // file to save in all the plots.

  TString hnamededxp="hdedxp_";
  hnamededxp+=(int)(trunclevel);
  hdedxp = new TH2F(hnamededxp, "dE/dx Vs p;p [Gev/c];dE/dx [keV/cm]", 300, 0, 3, 300, 0, 12);
  hchi2ndf = new TH1F("hchi2ndf", "counts;#chi^{2}/NDF", 100, 0, 30);
  hdedxtheta = new TH2F("hdedxtheta","dE/dx Vs #theta;#theta [deg];dE/dx [keV/cm]",300, 0, 100, 300, 0, 5);
  func_dedx_p1 = new TF1("func_dedx_p1",bethblock, 0, 10, 3);
  func_dedx_p1->SetParameters(0.938, 0.10, 9e-07);
  func_dedx_p2 = new TF1("func_dedx_p2",bethblock, 0, 10, 3);
  func_dedx_p2->SetParameters(0.938, 0.25, 1.60585e-06);
  hddedxp = new TH2F("hddedxp", "(dE/dx_{proton} - dE/dx_{Measured}) Vs p;p [Gev/c];#Delta dE/dx [keV/cm]", 300, 0.4, 0.8, 300, -9, 4);
  hdedxres = new TH1F("hdedxres", "hdedxres;#sigma", 300, 0, 3);
  h3 = new TH3F("h3","dE/dx vs. p vs. #theta;p [Gev/c];#theta [deg];dE/dx [keV/cm]",300,0,3, 300,0,180, 300,0,14);
  func_dedx_ex1 = new TF1("func_dedx_ex1","pol3");
  func_dedx_ex2 = new TF1("func_dedx_ex2","pol3");
  func_dEdxCut_SelectHeavy = new TF1("func_dEdxCut_SelectHeavy", "10*exp(-1*[0]*x + [1]) + 10*[2]", 0, 12);
  if(trunclevel==0) {func_dEdxCut_SelectHeavy->SetParameters(4.5, 0.3, 0.31); func_dedx_ex2->SetParameters(25.357149, -60.331129, 60.140292, -21.611668);}
  if(trunclevel==5) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.28); func_dedx_ex2->SetParameters(23.809325, -57.130155, 57.620590, -21.038446);}
  if(trunclevel==10) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.28); func_dedx_ex2->SetParameters(23.304850, -56.430713, 57.444653, -21.139907);}
  if(trunclevel==15) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.26); func_dedx_ex2->SetParameters(21.738962, -50.684162, 49.432553, -17.351713);}
  if(trunclevel==20) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.26); func_dedx_ex2->SetParameters(21.242925, -49.344271, 47.607085, -16.362078);}
  if(trunclevel==25) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.24); func_dedx_ex2->SetParameters(20.604513, -47.936337, 46.439757, -16.168312);}
  if(trunclevel==30) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.24); func_dedx_ex2->SetParameters(20.558956, -48.971782, 48.634803, -17.386550);}
  if(trunclevel==35) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.23); func_dedx_ex2->SetParameters(19.812376, -46.578144, 45.301752, -15.705259);}
  if(trunclevel==40) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.23); func_dedx_ex2->SetParameters(18.856370, -43.087046, 40.292455, -13.175311);}
  if(trunclevel==45) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.22); func_dedx_ex2->SetParameters(18.734271, -43.977406, 42.472483, -14.523108);}
  if(trunclevel==50) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.22); func_dedx_ex2->SetParameters(18.636307, -44.703826, 44.450140, -15.818927);}
  if(trunclevel==55) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.20); func_dedx_ex2->SetParameters(17.893084, -42.682331, 41.846847, -14.580762);}
  if(trunclevel==60) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.20); func_dedx_ex2->SetParameters(16.340613, -36.335950, 32.198897, -9.547386);}
  if(trunclevel==65) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.18); func_dedx_ex2->SetParameters(15.989939, -37.064222, 35.204961, -11.875087);}
  if(trunclevel==70) {func_dEdxCut_SelectHeavy->SetParameters(5, 0.25, 0.18); func_dedx_ex2->SetParameters(15.012380, -33.353548, 29.303274, -8.500944);}
  if(trunclevel==75) {func_dEdxCut_SelectHeavy->SetParameters(6.3, 0.18, 0.18); func_dedx_ex2->SetParameters(13.401793, -30.079924, 27.710779, -8.891197);}
  if(trunclevel==80) {func_dEdxCut_SelectHeavy->SetParameters(6.3, 0.18, 0.18); func_dedx_ex2->SetParameters(13.172396, -30.957738, 30.095553, -10.181711);}
  if(trunclevel==85) {func_dEdxCut_SelectHeavy->SetParameters(13, 0.3, 0.19); func_dedx_ex2->SetParameters(10.067001, -21.161270, 18.804877, -6.158741);}
  if(trunclevel==90) {func_dEdxCut_SelectHeavy->SetParameters(13, 0.3, 0.17); func_dedx_ex2->SetParameters(9.876912, -22.159733, 20.484941, -6.769080);}
  if(trunclevel==95) {func_dEdxCut_SelectHeavy->SetParameters(13, 0.3, 0.14); func_dedx_ex2->SetParameters(10.224370, -26.325900, 28.968147, -11.451508);}

  hprot = new TH2F("hprot", "dE/dx Vs p;p [Gev/c];dE/dx [keV/cm]", 300, 0.4, 0.8, 300, 0.0, 15);
  hptheta = new TH2F("hptheta","Momentum Vs #theta; #theta [deg]; P [GeV/c]",100, 0, 180, 100, 0, 3);
  hpthetaprot = new TH2F("hpthetaprot","Momentum Vs #theta (Protons); #theta [deg]; P [GeV/c]",100, 0, 180, 100, 0.4, 0.8);
  TString hnamepthetares="hpthetares_";
  hnamepthetares+=(int)(trunclevel);
  hpthetares = new TH2F(hnamepthetares, "Momentum vs. #theta vs.resolution;#theta [deg]; P [GeV/c]; resolution", 20,20,60,20,0.4,0.8);
  TString hnamepthetasep="hpthetasep_";
  hnamepthetasep+=(int)(trunclevel);
  hpthetasep = new TH2F(hnamepthetasep, "Momentum vs. #theta vs. separation power;#theta [deg]; P [GeV/c]; separation power", 20,20,60,20,0.4,0.8);

  fdedxsig = new TF1("fdedxsig","gaus",0.35,1);
  fdedxbkg = new TF1("fdedxbkg", "landau", 0, 0.35);
  fdedx = new TF1("fdedx", "landau(0)+gaus(3)", 0, 1.);

  thetamin = hpthetares->GetXaxis()->GetBinLowEdge(1);
  thetamax = hpthetares->GetXaxis()->GetBinUpEdge(20);
  ntheta = hpthetares->GetXaxis()->GetNbins();
  thetastep = (thetamax-thetamin)/ntheta;
  pmin = hpthetares->GetYaxis()->GetBinLowEdge(1);
  pmax = hpthetares->GetYaxis()->GetBinUpEdge(20);
  np = hpthetares->GetYaxis()->GetNbins();
  pstep = (pmax-pmin)/np;
  for(int ip=0; ip<np; ++ip)
  {
    for(int itheta=0; itheta<ntheta; ++itheta)
    {
      TString hnamededxav="hdedxav_";
      hnamededxav+=ip;
      hnamededxav+="_";
      hnamededxav+=itheta;
      hnamededxav+="_";
      hnamededxav+=(int)(trunclevel);
      // hdedxav[ip][itheta]->SetMarkerStyle(6);
      hdedxav[ip][itheta] = new TH1F(hnamededxav, ";#DeltadE/dx [keV/cm];Counts", 100, -8, 8);
    }
  }
}

void piddedx::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

Bool_t piddedx::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either cdcdedxchi2::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  GetEntry(entry);
  if(entry%1000==0) cout <<entry<<endl;
  // if(entry==30) hdedx->Fill(dedx[i]);
  if (dofit)
  {
    // hdedx->Reset();
    // for (int i=0;i<nhits;++i)
    // {
    //   if(entry>30) hdedx->Fill(dedx[i]*1e5);
    // }
    // hdedx->Fit("landau","l0q");
    // double fitpeak = hdedx->GetFunction("landau")->GetParameter(1);
    // hdedxp->Fill(p, fitpeak);
  }
  else
  {
    sort(dedx, dedx+nhits);
    int ntrunc = nhits * (100-trunclevel)/100;
    double avgdedx=0;
    double sumde = 0;
    double sumdx = 0;
    for (int i=0;i<ntrunc;++i)
    {
      sumde+= de[i];
      sumdx+= dx[i];
      // hdedx->Fill(dedx[i]*1e5);
    }

    avgdedx = sumde/sumdx;

    hdedxp->Fill(p, avgdedx*1e6);
    hdedxtheta->Fill(theta, avgdedx*1e6);
    h3->Fill(p, theta, avgdedx*1e6);

    // focus on keeping protons
    if(avgdedx*1e6 > func_dEdxCut_SelectHeavy->Eval(p) && p > 0.40 && p < 0.80)
    {
      hprot->Fill(p, avgdedx*1e6);
      hpthetaprot->Fill(theta, p);
    }
    // hprot->Fill(p, avgdedx*1e5 - func_dedx_ex2->Eval(p));

    hddedxp->Fill(p, avgdedx*1e6 - func_dedx_ex2->Eval(p));
    hptheta->Fill(theta, p);

    for(int ip=0; ip<np; ++ip)
    {
      pcut = pmin+ip*pstep;

      for(int itheta=0; itheta<ntheta; ++itheta)
      {
        thetacut = thetamin+itheta*thetastep;
        double dp = p-pcut, dtht = theta-thetacut;
        //if(abs(p-pcut)<pstep && abs(theta-thetacut)<thetastep)
        if (dp>0 && dp<pstep && dtht>0 && dtht<thetastep)
        {
          // hdedxav[ip][itheta]->Fill(avgdedx*1e5 - func_dedx_ex2->Eval(p));
          hdedxav[ip][itheta]->Fill(avgdedx*1e6 - func_dedx_ex2->Eval(p));
        }
      }
    }
    // }

  }

  return kTRUE;
}

void piddedx::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void piddedx::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  //
  // TCanvas *cdedx = new TCanvas("cdedx","cdedx",700,500);
  // cdedx->cd();
  // hdedx->Draw();
  // TLine *ldedx = new TLine(hdedx->GetMean(), 0, hdedx->GetMean(), 2);
  // ldedx->SetLineStyle(9);
  // ldedx->SetLineColor(kRed);
  // ldedx->SetLineWidth(2);
  // ldedx->Draw();

  gStyle->SetOptStat(0);

  TString cname1 = "cdedxav_0_4_";
  cname1 +=(int)(trunclevel);
  TCanvas *cdedxav1 = new TCanvas(cname1,cname1,1500,1500);
  cdedxav1->Divide(10,10,0.00001,0.00001);
  TString cname2 = "cdedxav_5_9_";
  cname2 +=(int)(trunclevel);
  TCanvas *cdedxav2 = new TCanvas(cname2,cname2,1500,1500);
  cdedxav2->Divide(10,10,0.00001,0.00001);
  TString cname3 = "cdedxav_10_14_";
  cname3 +=(int)(trunclevel);
  TCanvas *cdedxav3 = new TCanvas(cname3,cname3,1500,1500);
  cdedxav3->Divide(10,10,0.00001,0.00001);
  TString cname4 = "cdedxav_15_19_";
  cname4 +=(int)(trunclevel);
  TCanvas *cdedxav4 = new TCanvas(cname4,cname4,1500,1500);
  cdedxav4->Divide(10,10,0.00001,0.00001);

  int k1 = 0;
  int k2 = 0;
  int k3 = 0;
  int k4 = 0;
  for(int ip=0; ip<np; ++ip)
  {
    for(int itheta=0; itheta<ntheta; ++itheta)
    {
      if(ip <= 4)
      {
        k1++;
        cdedxav1->cd(k1);
        hdedxav[ip][itheta]->Draw();
      }
      if(ip >= 5 && ip <= 9)
      {
        k2++;
        cdedxav2->cd(k2);
        hdedxav[ip][itheta]->Draw();
      }
      if(ip >= 10 && ip <= 14)
      {
        k3++;
        cdedxav3->cd(k3);
        hdedxav[ip][itheta]->Draw();
      }
      if(ip >= 15 && ip <= 19)
      {
        k4++;
        cdedxav4->cd(k4);
        hdedxav[ip][itheta]->Draw();
      }

      // double binmax = 0;
      // int bestbin = -1;
      // for (int i=35;i<90;++i)
      // {
      //   if (hdedxav[ip][itheta]->GetBinContent(i) > binmax)
      //   {binmax = hdedxav[ip][itheta]->GetBinContent(i); bestbin=i;}
      // }
      // double xmin = hdedxav[ip][itheta]->GetBinCenter(bestbin-7);
      // double xmax = hdedxav[ip][itheta]->GetBinCenter(bestbin+7);

      TSpectrum *s = new TSpectrum();
      int nfound = s->Search(hdedxav[ip][itheta],3,"",0.01);
      float *xpeaks = s->GetPositionX();
      int binpeak0 = -1;
      int binpeak1 = -1;
      double xpeak0 = -1;
      double xpeak1 = -1;
      double dxpeaks = -1;

      for (int p=0;p<nfound;p++)
      {
        xpeak0 = xpeaks[0];
        xpeak1 = xpeaks[1];
        if(xpeak0 > xpeak1) {double tmp=xpeak0;xpeak0=xpeak1;xpeak1=tmp;}
        binpeak0 = hdedxav[ip][itheta]->GetXaxis()->FindBin(xpeak0);
        binpeak1 = hdedxav[ip][itheta]->GetXaxis()->FindBin(xpeak1);
        dxpeaks = xpeak1 - xpeak0;
      }
      double xmin0 = hdedxav[ip][itheta]->GetBinCenter(binpeak0-2);
      double xmax0 = xpeak0 + dxpeaks*.5; //hdedxav[ip][itheta]->GetBinCenter(binpeak0+10);
      double xmin1 = xpeak1 - dxpeaks*.4; //hdedxav[ip][itheta]->GetBinCenter(binpeak1-10);
      double xmax1 = hdedxav[ip][itheta]->GetBinCenter(binpeak1+9);

      double par[6];
      // fdedxbkg->SetParLimits(1,xpeak0*0.9,xpeak0*1.1);
      fdedxbkg->SetLineColor(4);
      hdedxav[ip][itheta]->Fit(fdedxbkg,"","",xmin0,xmax0);
      fdedxbkg->GetParameters(&par[0]);
      // fdedxsig->SetParLimits(1,xpeak1*0.9,xpeak1*1.1);
      hdedxav[ip][itheta]->Fit(fdedxsig,"+","",xmin1,xmax1);
      fdedxsig->GetParameters(&par[3]);
      fdedx->SetParameters(par);
      // fdedx->SetParLimits(1,par[1]*.9,par[1]*1.1);
      // fdedx->SetParLimits(4,par[4]*.9,par[4]*1.1);
      hdedxav[ip][itheta]->Fit(fdedx,"R","",xmin0,xmax1);

      hdedxav[ip][itheta]->Write();
      res= par[5]; // hdedxav[ip][itheta]->GetFunction("fdedx")->GetParameter(4);
      sep= abs(par[1]-par[4])/(sqrt(par[2]*par[2]+par[5]*par[5]));
      // double chi2ndf = fdedx->GetChisquare()/fdedx->GetNDF();
      // hchi2ndf->Fill(chi2ndf);
      // if(chi2ndf < 10)
      hpthetares->SetBinContent(itheta+1, ip+1, res);
      hpthetasep->SetBinContent(itheta+1, ip+1, sep);

      // hpthetares->Fill(res[ip][itheta]);
      // printf("************* ip= %d | itheta= %d | res = %f | separation power = %f\n",ip,itheta,res,sep);
      // printf("#peaks=%d|\n binpeak0=%d| xpeak0=%f| xmin0=%f| xmax0=%f|\n binpeak1=%d| xpeak1=%f| xmin1=%f| xmax1=%f\n",nfound,binpeak0,xpeak0,xmin0,xmax0,binpeak1,xpeak1,xmin1,xmax1);
    }
  }
  cdedxav1->Print(Form("fig/%s.eps",cname1.Data()),"eps");
  cdedxav2->Print(Form("fig/%s.eps",cname2.Data()),"eps");
  cdedxav3->Print(Form("fig/%s.eps",cname3.Data()),"eps");
  cdedxav4->Print(Form("fig/%s.eps",cname4.Data()),"eps");
  cdedxav1->Print(Form("fig/%s.root",cname1.Data()),"root");
  cdedxav2->Print(Form("fig/%s.root",cname2.Data()),"root");
  cdedxav3->Print(Form("fig/%s.root",cname3.Data()),"root");
  cdedxav4->Print(Form("fig/%s.root",cname4.Data()),"root");

  //
  // TCanvas *cchi2ndf = new TCanvas("cchi2ndf","cchi2ndf",700,500);
  // cchi2ndf->cd();
  // hchi2ndf->Draw();

  TString cnamepthetasep = "cpthetasep_";
  cnamepthetasep +=(int)(trunclevel);
  TCanvas *cpthetasep= new TCanvas(cnamepthetasep,cnamepthetasep,700,500);
  cpthetasep->cd();
  gStyle->SetPalette(1);
  hpthetasep->GetZaxis()->SetRangeUser(0,7);
  hpthetasep->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpthetasep->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetasep->Modified();
  hpthetasep->Write();
  cpthetasep->Print(Form("fig/%s.eps",cnamepthetasep.Data()),"eps");
  cpthetasep->Print(Form("fig/%s.root",cnamepthetasep.Data()),"root");

  TString cnamepthetares = "cpthetares_";
  cnamepthetares +=(int)(trunclevel);
  TCanvas *cpthetares= new TCanvas(cnamepthetares,cnamepthetares,700,500);
  cpthetares->cd();
  // gStyle->SetPalette(1);
  hpthetares->GetZaxis()->SetRangeUser(0,2);
  hpthetares->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetares->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetares->Modified();
  hpthetares->Write();
  cpthetares->Print(Form("fig/%s.eps",cnamepthetares.Data()),"eps");
  cpthetares->Print(Form("fig/%s.root",cnamepthetares.Data()),"root");

  TCanvas *cptheta=new TCanvas("cptheta","cptheta",700,500);
  cptheta->cd();
  hptheta->Draw("colz");
  TLine *lptheta_x1 = new TLine(20, 0.40, 20, 0.80);
  TLine *lptheta_x2 = new TLine(60, 0.40, 60, 0.80);
  TLine *lptheta_y1 = new TLine(20, 0.40, 60, 0.40);
  TLine *lptheta_y2 = new TLine(20, 0.80, 60, 0.80);
  // lptheta_x1->SetLineColor(kRed);
  // lptheta_x2->SetLineColor(kRed);
  // lptheta_y1->SetLineColor(kRed);
  // lptheta_y2->SetLineColor(kRed);
  // lptheta_x1->SetLineStyle(9);
  // lptheta_x2->SetLineStyle(9);
  // lptheta_y1->SetLineStyle(9);
  // lptheta_y2->SetLineStyle(9);
  lptheta_x1->SetLineWidth(2);
  lptheta_x2->SetLineWidth(2);
  lptheta_y1->SetLineWidth(2);
  lptheta_y2->SetLineWidth(2);
  lptheta_x1->Draw();
  lptheta_x2->Draw();
  lptheta_y1->Draw();
  lptheta_y2->Draw();
  hptheta->Write();
  cptheta->Print("fig/cptheta.eps","eps");
  cptheta->Print("fig/cptheta.root","root");

  TCanvas *cpthetaprot=new TCanvas("cpthetaprot","cpthetaprot",700,500);
  cpthetaprot->cd();
  hpthetaprot->Draw("colz");
  hpthetaprot->Write();
  cpthetaprot->Print("fig/cpthetaprot.eps","eps");
  cpthetaprot->Print("fig/cpthetaprot.root","root");

  TString cnameprot = "cprot_";
  cnameprot +=(int)(trunclevel);
  TCanvas *cprot=new TCanvas(cnameprot,cnameprot,30,30,1500,500);
  cprot->Divide(3,1);
  cprot->cd(1);
  hprot->Draw("col");
  hprot->FitSlicesY();
  cprot->cd(2);
  TH1F *hprot_1 = (TH1F*)gDirectory->Get("hprot_1");
  hprot_1->SetTitle("dE/dx_{proton}");
  hprot_1->SetMinimum(1);
  hprot_1->SetMaximum(10);
  hprot_1->Draw();
  hprot_1->Fit("pol3");
  func_dedx_ex1 = hprot_1->GetFunction("pol3");
  double para0 = func_dedx_ex1->GetParameter(0);
  double para1 = func_dedx_ex1->GetParameter(1);
  double para2 = func_dedx_ex1->GetParameter(2);
  double para3 = func_dedx_ex1->GetParameter(3);
  printf("func_dedx_ex2->SetParameters(%f, %f, %f, %f);\n",para0, para1, para2, para3);
  cprot->cd(3);
  hddedxp->Draw("col");
  TH1F *hprot_2 = (TH1F*)gDirectory->Get("hprot_2");
  hprot_2->SetTitle("resolution");
  hprot_2->SetMinimum(0.05);
  hprot_2->SetMaximum(0.25);
  // hprot_2->Draw();
  cprot->Print(Form("fig/%s.eps",cnameprot.Data()),"eps");
  cprot->Print(Form("fig/%s.root",cnameprot.Data()),"root");

  TString cnamededxp = "cdedxp_";
  cnamededxp +=(int)(trunclevel);
  TCanvas *cdedxp=new TCanvas(cnamededxp,cnamededxp,700,500);
  TLine *l4_1 = new TLine(0.40, 0, 0.40, 10);
  TLine *l4_2 = new TLine(0.80, 0, 0.80, 10);
  l4_1->SetLineColor(kRed);
  l4_1->SetLineStyle(9);
  l4_1->SetLineWidth(2);
  l4_2->SetLineColor(kRed);
  l4_2->SetLineStyle(9);
  l4_2->SetLineWidth(2);
  cdedxp->cd();
  hdedxp->Draw("colz");
  func_dEdxCut_SelectHeavy->Draw("same");
  l4_1->Draw();
  l4_2->Draw();
  hdedxp->Write();
  // func_dedx_p1->Draw("same");
  cdedxp->Print(Form("fig/%s.eps",cnamededxp.Data()),"eps");
  cdedxp->Print(Form("fig/%s.root",cnamededxp.Data()),"root");

  TString cnamededxtheta = "cdedxtheta_";
  cnamededxtheta +=(int)(trunclevel);
  TCanvas *cdedxtheta = new TCanvas(cnamededxtheta,cnamededxtheta,700,500);
  // c2->Divide(3,1);
  cdedxtheta->cd();
  hdedxtheta->Draw("colz");
  hdedxtheta->Write();
  cdedxtheta->Print(Form("fig/%s.eps",cnamededxtheta.Data()),"eps");
  cdedxtheta->Print(Form("fig/%s.root",cnamededxtheta.Data()),"root");

  outputfig->Print();
  outputfig->Close();
}
