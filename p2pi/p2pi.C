#define p2pi_cxx
#include "p2pi.h"
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
#include <vector>

// Splits a TString into vector of strings; separation character contained in delim

typedef std::vector<TString> StrVec;

int SplitString(TString s, TString delim, StrVec &toks)
{
  toks.clear();
  TObjArray *tok = s.Tokenize(delim);
  int N = tok->GetEntries();
  for (int i=0;i<N;++i)
  {
    TString token = (((TObjString*)tok->At(i))->String()).Strip(TString::kBoth);
    toks.push_back(token);
  }
  return toks.size();
}

int trunclevel = 0;
int trunclevel2 = 0;
TString particle;

TFile *outputfig;
TH2F *hdedxp;
TH2F *hdedxtheta;
TH2F *hpthetall;
TH2F *hptheta;
TH1 *hdedxproj[200];
TGraph *grdedxexp = new TGraph();
TF1 *fdedxexp;
TH2F *hddedxp;
TH1F *hddedxav[20][20];
TH1F *hntracks[20][20];
TH1F *hnhits[20][20];
TH1F *hntracksall;
TH1F *hnhitsall;
int ntheta, np;
double pmin, pmax, pstep, pcut, thetamin, thetamax, thetastep, thetacut;
TH2F *hpthetantracks;
TH2F *hpthetanhits;

double res;
TProfile *profdedxtheta[20];

void p2pi::Begin(TTree * /*tree*/)
{
  TString option = GetOption();
  // trunclevel=option.Atoi();
  // TString truncstr = option(0,option.Index(":"));
  // trunclevel = truncstr.Atoi();
  // particle = option(option.Index(":")+1,1000);

  StrVec v;
  SplitString(option,":",v);
  trunclevel  = v[0].Atoi();
  trunclevel2 = v[1].Atoi();
  particle    = v[2];

  for (unsigned int i=0;i<v.size(); ++i)
  {
    cout <<v[i]<<endl;
  }

  outputfig = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root",particle.Data(),particle.Data()),"UPDATE");  // file to save in all the plots.
  hdedxp = new TH2F("hdedxp", "dE/dx Vs p;p [Gev/c];dE/dx [keV/cm]", 300, 0, 3, 300, 0, 30);
  hdedxtheta = new TH2F("hdedxtheta","dE/dx Vs #theta;#theta [deg];dE/dx [keV/cm]",300, 0, 100, 300, 0, 30);
  hpthetall = new TH2F("hpthetall","Momentum Vs #theta; #theta [deg]; P [GeV/c]",300, 0, 100, 300, 0.3, 1);
  fdedxexp = new TF1("fdedxexp","pol6",0.3,1.);
  if(trunclevel2==0) fdedxexp->SetParameters(97.912970, -635.688032, 1970.100727, -3387.329892,3322.266135, -1742.237229, 379.620155);
  if(trunclevel2==5) fdedxexp->SetParameters(96.757570, -630.127825, 1938.303650, -3301.260749,3207.641088, -1668.608514, 361.403850);
  if(trunclevel2==10) fdedxexp->SetParameters(97.519521, -655.021335, 2080.263153, -3663.469347,3682.344380, -1981.192432, 443.460578);
  if(trunclevel2==15) fdedxexp->SetParameters(100.175880, -690.652518, 2241.428559, -4026.414018,4121.025784, -2253.541192, 511.689073);
  if(trunclevel2==20) fdedxexp->SetParameters(98.081956, -674.890423, 2182.782070, -3907.282612,3984.826159, -2171.071283, 491.093684);
  if(trunclevel2==25) fdedxexp->SetParameters(97.630921, -675.990562, 2195.330749, -3945.097093,4039.269122, -2209.665159, 501.916230);
  if(trunclevel2==30) fdedxexp->SetParameters(96.462791, -670.891592, 2185.608209, -3940.429121,4048.271630, -2222.350229, 506.575117);
  if(trunclevel2==35) fdedxexp->SetParameters(96.033166, -670.342411, 2185.280328, -3939.116174,4044.517730, -2218.483368, 505.224089);
  if(trunclevel2==40) fdedxexp->SetParameters(95.407153, -669.338055, 2188.545305, -3953.096578,4063.887870, -2230.220756, 507.812718);
  if(trunclevel2==45) fdedxexp->SetParameters(95.965496, -681.833137, 2253.058404, -4111.109156,4267.994616, -2364.473612, 543.274364);
  if(trunclevel2==50) fdedxexp->SetParameters(95.753097, -681.702537, 2251.155222, -4101.589134,4250.713278, -2350.859821, 539.318370);
  if(trunclevel2==55) fdedxexp->SetParameters(90.776002, -639.196106, 2091.707970, -3784.194916,3900.359845, -2147.989616, 491.187217);
  if(trunclevel2==60) fdedxexp->SetParameters(91.302883, -649.038669, 2137.416696, -3886.263293,4021.733943, -2222.066508, 509.463303);
  if(trunclevel2==65) fdedxexp->SetParameters(91.517852, -658.946319, 2193.334256, -4027.758019,4206.880103, -2344.370503, 541.782645);
  if(trunclevel2==70) fdedxexp->SetParameters(91.186367, -659.314126, 2196.996761, -4036.345113,4217.470350, -2351.522476, 543.856313);
  if(trunclevel2==75) fdedxexp->SetParameters(91.009946, -661.567038, 2208.353965, -4057.854555,4236.035187, -2357.826943, 544.070573);
  if(trunclevel2==80) fdedxexp->SetParameters(90.928768, -666.696617, 2236.549201, -4124.589796,4317.710773, -2408.544160, 556.742123);
  if(trunclevel2==85) fdedxexp->SetParameters(85.253717, -615.437423, 2033.171490, -3698.544140,3824.973552, -2110.637659, 483.169492);
  if(trunclevel2==90) fdedxexp->SetParameters(67.131048, -443.881399, 1356.027197, -2298.646349,2229.332182, -1159.441512, 251.258291);
  if(trunclevel2==95) fdedxexp->SetParameters(58.306608, -366.997876, 1069.196578, -1733.089614,1610.681333, -803.685116, 167.226217);
  hddedxp = new TH2F("hddedxp", "#DeltadE/dx= dEdx_{measured} - dEdx_{expected} Vs p;p [Gev/c];#DeltadE/dx [keV/cm]", 300, 0.3, 1., 300, -10, 30);

  hptheta = new TH2F("hptheta", "Momentum Vs #theta; #theta [deg]; P [GeV/c]",20, 60, 80, 20, 0.3, 1);

  TString hnamepthetantracks=Form("hpthetantracks_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  hpthetantracks = new TH2F(hnamepthetantracks, "Momentum vs. #theta vs. number of tracks;#theta [deg]; P [GeV/c]; # tracks", 20,60,80,20,0.30,1);
  TString hnamentracksall=Form("hntracksall_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  hntracksall = new TH1F(hnamentracksall, "number of tracks after truncation ;# tracks;Counts", 20, 0, 40);

  TString hnamepthetanhits=Form("hpthetanhits_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  hpthetanhits = new TH2F(hnamepthetanhits, Form("Momentum vs. #theta vs. number of hits (%s);#theta [deg]; P [GeV/c]; # hits",particle.Data()), 20,60,80,20,0.30,1);
  TString hnamenhitsall=Form("hnhitsall_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  hnhitsall = new TH1F(hnamenhitsall, Form("number of hits after truncation (%s);# hits;Counts",particle.Data()), 20, 0, 40);

  thetamin = hptheta->GetXaxis()->GetBinLowEdge(1);
  thetamax = hptheta->GetXaxis()->GetBinUpEdge(20);
  ntheta = hptheta->GetXaxis()->GetNbins();
  thetastep = (thetamax-thetamin)/ntheta;
  pmin = hptheta->GetYaxis()->GetBinLowEdge(1);
  pmax = hptheta->GetYaxis()->GetBinUpEdge(20);
  np = hptheta->GetYaxis()->GetNbins();
  pstep = (pmax-pmin)/np;
  for(int ip=0; ip<np; ++ip)
  {
    for(int itheta=0; itheta<ntheta; ++itheta)
    {
      TString hnameddedxav=Form("hddedxav_%s_%d_%d_%d_%d",particle.Data(),ip,itheta,trunclevel,trunclevel2);
      if(particle == "protonflat")
      hddedxav[ip][itheta] = new TH1F(hnameddedxav, Form("Energy loss in CDC (%s);#DeltadE/dx [keV/cm];Counts",particle.Data()), 100, -10, 30);
      if(particle == "proton" || particle == "pip" || particle == "pim")
      hddedxav[ip][itheta] = new TH1F(hnameddedxav, Form("Energy loss in CDC (%s);dE/dx [keV/cm];Counts",particle.Data()), 100, 0, 30);

      TString hnamentracks=Form("hntracks_%s_%d_%d_%d_%d",particle.Data(),ip,itheta,trunclevel,trunclevel2);
      hntracks[ip][itheta] = new TH1F(hnamentracks, Form("Number of tracks (%s);# tracks;Counts",particle.Data()), 20, 0, 20);

      TString hnamenhits=Form("hnhits_%s_%d_%d_%d_%d",particle.Data(),ip,itheta,trunclevel,trunclevel2);
      hnhits[ip][itheta] = new TH1F(hnamenhits, Form("Number of hits (%s);# hits;Counts",particle.Data()), 20, 0, 40);
    }
    TString hnameprofdedxtheta=Form("hprofdedxtheta_%s_%d_%d_%d",particle.Data(),ip,trunclevel,trunclevel2);
    profdedxtheta[ip] = new TProfile(hnameprofdedxtheta,Form("dE/dx vs #theta (%s);#theta [deg];dE/dx [keV/cm]",particle.Data()),100,60,80,0,20);
  }

}

void p2pi::SlaveBegin(TTree * /*tree*/)
{
  TString option = GetOption();
}

Bool_t p2pi::Process(Long64_t entry)
{
  GetEntry(entry);

  // if(nhits<10) return kFALSE;
  // if(entry%1000000==0) cout <<entry<<endl;
  // sort(dedx, dedx+nhits);
  // if (nhits>10) cout <<"dedx= "<<dedx[0]<<" "<<dedx[1]<<" "<<dedx[2]<<" "<<dedx[3]<<" "<<dedx[4]<<" "<<dedx[5]<<" "<<dedx[6]<<" "<<dedx[7]<<" "<<dedx[8]<<" "<<dedx[9]<<" "<<dedx[10]<<" "<<endl;
  // if (nhits>10) cout <<"de/dx= "<<de[0]/dx[0]<<" "<<de[1]/dx[1]<<" "<<de[2]/dx[2]<<" "<<de[3]/dx[3]<<" "<<de[4]/dx[4]<<" "<<de[5]/dx[5]<<" "<<de[6]/dx[6]<<" "<<de[7]/dx[7]<<" "<<de[8]/dx[8]<<" "<<de[9]/dx[9]<<" "<<de[10]/dx[10]<<" "<<endl;
  // if (nhits>10) cout <<"dx= "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<" "<<dx[3]<<" "<<dx[4]<<" "<<dx[5]<<" "<<dx[6]<<" "<<dx[7]<<" "<<dx[8]<<" "<<dx[9]<<" "<<dx[10]<<" "<<endl;
  // if(entry%1000000==0) printf("ntracks = %d | nhits = %d\n", ntracks, nhits);
  int ntrunc2 = nhits * (100.-(double)trunclevel2)/100.;
  // sort(dedx, dedx+ntrunc2, std::greater<int>());
  int ntrunc = nhits * (double)trunclevel/100.;

  // if (nhits>10) cout<<"dedx= "<<dedx[0]<<" | ntrunc2 = "<<ntrunc2<<endl;
  // if (nhits>10) cout<<"dedx= "<<dedx[0]<<" | ntrunc = "<<ntrunc<<endl;
  // if(ntrunc2-ntrunc<10) return kFALSE;

  double avgdedx=0;
  double sumde = 0;
  double sumdx = 0;
  for (int i=ntrunc;i<ntrunc2;++i)
  {
    sumde+= de[i];
    sumdx+= dx[i];
  }
  if(entry%1000000==0) printf("entry = %d | ntracks=%d | nhits=%d | min=%d | max=%d | dh = %d\n",entry, ntracks, nhits,ntrunc, ntrunc2, (ntrunc2-ntrunc));
  avgdedx = sumde/sumdx;

  hntracksall->Fill(ntracks);
  hnhitsall->Fill(ntrunc2-ntrunc);
  hdedxp->Fill(p, avgdedx*1e6);
  hdedxtheta->Fill(theta, avgdedx*1e6);
  hpthetall->Fill(theta, p);
  hptheta->Fill(theta, p);
  hddedxp->Fill(p,avgdedx*1e6 - fdedxexp->Eval(p));

  for(int ip=0; ip<np; ++ip)
  {
    pcut = pmin+ip*pstep;
    double dp = p-pcut;
    for(int itheta=0; itheta<ntheta; ++itheta)
    {
      thetacut = thetamin+itheta*thetastep;
      double dtht = theta-thetacut;
      //if(abs(p-pcut)<pstep && abs(theta-thetacut)<thetastep)
      if (dp>0 && dp<pstep && dtht>0 && dtht<thetastep)
      {
        // hddedxav[ip][itheta]->Fill(avgdedx*1e5 - func_dedx_ex2->Eval(p));
        if(particle == "protonflat")
        hddedxav[ip][itheta]->Fill(avgdedx*1e6 - fdedxexp->Eval(p));
        if(particle == "proton" || particle == "pip" || particle == "pim")
        hddedxav[ip][itheta]->Fill(avgdedx*1e6);

        hntracks[ip][itheta]->Fill(ntracks);
        hnhits[ip][itheta]->Fill(ntrunc2-ntrunc);
      }
    }
    if(dp>0 && dp<pstep) profdedxtheta[ip]->Fill(theta, avgdedx*1e6);
  }

  return kTRUE;
}

void p2pi::SlaveTerminate()
{

}

void p2pi::Terminate()
{
  gStyle->SetOptStat(0);

  TString cname1 = Form("cddedxav_0_4_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cddedxav1 = new TCanvas(cname1,cname1,1500,1500);
  cddedxav1->Divide(10,10,0.00001,0.00001);
  TString cname2 = Form("cddedxav_5_9_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cddedxav2 = new TCanvas(cname2,cname2,1500,1500);
  cddedxav2->Divide(10,10,0.00001,0.00001);
  TString cname3 = Form("cddedxav_10_14_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cddedxav3 = new TCanvas(cname3,cname3,1500,1500);
  cddedxav3->Divide(10,10,0.00001,0.00001);
  TString cname4 = Form("cddedxav_15_19_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cddedxav4 = new TCanvas(cname4,cname4,1500,1500);
  cddedxav4->Divide(10,10,0.00001,0.00001);

  TString cnameprofdedxtheta = Form("cprofdedxtheta_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cprofdedxtheta = new TCanvas(cnameprofdedxtheta,cnameprofdedxtheta,1500,1500);
  cprofdedxtheta->Divide(4,5,0.00001,0.00001);

  int k1 = 0;
  int k2 = 0;
  int k3 = 0;
  int k4 = 0;
  int k0 = 0;

  for(int ip=0; ip<np; ++ip)
  {
    for(int itheta=0; itheta<ntheta; ++itheta)
    {
      if(ip <= 4)
      {
        k1++;
        cddedxav1->cd(k1);
        hddedxav[ip][itheta]->Draw();
      }
      if(ip >= 5 && ip <= 9)
      {
        k2++;
        cddedxav2->cd(k2);
        hddedxav[ip][itheta]->Draw();
      }
      if(ip >= 10 && ip <= 14)
      {
        k3++;
        cddedxav3->cd(k3);
        hddedxav[ip][itheta]->Draw();
      }
      if(ip >= 15 && ip <= 19)
      {
        k4++;
        cddedxav4->cd(k4);
        hddedxav[ip][itheta]->Draw();
      }

      hddedxav[ip][itheta]->Write();
      hpthetantracks->SetBinContent(itheta+1, ip+1, hntracks[ip][itheta]->GetMean());
      hpthetanhits->SetBinContent(itheta+1, ip+1, hnhits[ip][itheta]->GetMean());
    }

    k0++;
    cprofdedxtheta->cd(k0);
    cprofdedxtheta->SetGrid();
    profdedxtheta[ip]->GetYaxis()->SetRangeUser(0.,20.);
    profdedxtheta[ip]->Draw();
  }

  cddedxav1->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cname1.Data()),"eps");
  cddedxav1->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cname1.Data()),"root");
  cddedxav2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cname2.Data()),"eps");
  cddedxav2->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cname2.Data()),"root");
  cddedxav3->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cname3.Data()),"eps");
  cddedxav3->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cname3.Data()),"root");
  cddedxav4->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cname4.Data()),"eps");
  cddedxav4->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cname4.Data()),"root");

  TString cnameddedxexp = "cddedxexp_";
  cnameddedxexp +=(int)(trunclevel);
  cnameddedxexp +="_";
  cnameddedxexp +=(int)(trunclevel2);
  TCanvas *cddedxexp = new TCanvas(cnameddedxexp,cnameddedxexp,700,500);
  cddedxexp->cd();
  hddedxp->Draw("colz");
  cddedxexp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnameddedxexp.Data()),"eps");
  cddedxexp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnameddedxexp.Data()),"root");

  TString cnamededxexp = "cdedxexp_";
  cnamededxexp +=(int)(trunclevel);
  cnamededxexp +="_";
  cnamededxexp +=(int)(trunclevel2);
  TCanvas *cdedxexp = new TCanvas(cnamededxexp,cnamededxexp,700,500);

  TString cnamededxproj = "cdedxproj_";
  cnamededxproj +=(int)(trunclevel);
  cnamededxproj +="_";
  cnamededxproj +=(int)(trunclevel2);
  TCanvas *cdedxproj = new TCanvas(cnamededxproj,cnamededxproj,1500,1500);
  cdedxproj->Divide(8,9,0.00001,0.00001);
  int k5=0;

  TF1 *fgaus = new TF1("fgaus","gaus",0.3,1);
  TSpectrum *s1 = new TSpectrum();
  for (int i=29;i<100;i++)
  {
    ++k5;
    cdedxproj->cd(k5);
    hdedxproj[i] = hdedxp->ProjectionY(Form("bin%d",i+1),i+1,i+2);
    hdedxproj[i]->Draw();

    int nfound1 = s1->Search(hdedxproj[i],13,"",1);
    float *xpeaks1 = s1->GetPositionX();
    int binpeak10 = -1;
    double xpeak10 = -1;
    double xpeak11 = -1;
    for (int p1=0;p1<nfound1;p1++)
    {
      xpeak10 = xpeaks1[0];
      xpeak11 = xpeaks1[1];
      binpeak10 = hdedxproj[i]->GetXaxis()->FindBin(xpeak10);
    }
    double xmin1 = hdedxproj[i]->GetBinCenter(binpeak10-13);
    double xmax1 = hdedxproj[i]->GetBinCenter(binpeak10+13);
    fgaus->SetParameter(1,xpeak10);
    hdedxproj[i]->Fit(fgaus,"R","",xmin1,xmax1);
    double param1 = fgaus->GetParameter(1);
    grdedxexp->SetPoint(i,hdedxp->GetXaxis()->GetBinCenter(i),param1);
  }

  grdedxexp->Fit(fdedxexp,"R","",0.3,1);
  if(trunclevel2==90) grdedxexp->Fit(fdedxexp,"R","",0.35,1);
  if(trunclevel2==95) grdedxexp->Fit(fdedxexp,"R","",0.45,0.90);
  double para0 = fdedxexp->GetParameter(0);
  double para1 = fdedxexp->GetParameter(1);
  double para2 = fdedxexp->GetParameter(2);
  double para3 = fdedxexp->GetParameter(3);
  double para4 = fdedxexp->GetParameter(4);
  double para5 = fdedxexp->GetParameter(5);
  double para6 = fdedxexp->GetParameter(6);
  printf("fdedxexp->SetParameters(%f, %f, %f, %f,%f, %f, %f);\n",para0, para1, para2, para3, para4, para5, para6);

  cdedxexp->cd();
  grdedxexp->SetMarkerStyle(20);
  grdedxexp->Draw("AP");
  fdedxexp->Write(Form("%s",cnamededxexp.Data()));
  cdedxexp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamededxexp.Data()),"eps");
  cdedxexp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamededxexp.Data()),"root");

  cdedxproj->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamededxproj.Data()),"eps");
  cdedxproj->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamededxproj.Data()),"root");

  TString cnamededxp = Form("cdedxp_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cdedxp=new TCanvas(cnamededxp,cnamededxp,700,500);
  cdedxp->cd();
  hdedxp->Draw("colz");
  cdedxp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamededxp.Data()),"eps");
  cdedxp->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamededxp.Data()),"root");

  TString cnamededxtheta = Form("cdedxtheta_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cdedxtheta = new TCanvas(cnamededxtheta,cnamededxtheta,700,500);
  cdedxtheta->cd();
  hdedxtheta->Draw("colz");
  cdedxtheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamededxtheta.Data()),"eps");
  cdedxtheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamededxtheta.Data()),"root");

  TCanvas *cpthetall=new TCanvas("cpthetall","cpthetall",700,500);
  cpthetall->cd();
  hpthetall->Draw("colz");
  cpthetall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/cpthetall.eps",particle.Data()),"eps");
  cpthetall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/cpthetall.root",particle.Data()),"root");

  TCanvas *cptheta= new TCanvas("cptheta","cptheta",700,500);
  cptheta->cd();
  hptheta->Draw("colz");
  hptheta->Write();
  cptheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/cptheta.eps",particle.Data()),"eps");
  cptheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/cptheta.root",particle.Data()),"root");

  cprofdedxtheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnameprofdedxtheta.Data()),"eps");
  cprofdedxtheta->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnameprofdedxtheta.Data()),"root");

  TString cnamepthetantracks = Form("cpthetantracks_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetantracks= new TCanvas(cnamepthetantracks,cnamepthetantracks,700,500);
  cpthetantracks->cd();
  hpthetantracks->Draw("colz");
  cpthetantracks->Modified();
  hpthetantracks->Write();
  cpthetantracks->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamepthetantracks.Data()),"eps");
  cpthetantracks->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamepthetantracks.Data()),"root");

  TString cnamentracksall = Form("cntracksall_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cntracksall = new TCanvas(cnamentracksall,cnamentracksall,700,500);
  cntracksall->cd();
  hntracksall->Draw();
  hntracksall->Write();
  cntracksall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamentracksall.Data()),"eps");
  cntracksall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamentracksall.Data()),"root");

  TString cnamepthetanhits = Form("cpthetanhits_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetanhits= new TCanvas(cnamepthetanhits,cnamepthetanhits,700,500);
  cpthetanhits->cd();
  hpthetanhits->Draw("colz");
  cpthetanhits->Modified();
  hpthetanhits->Write();
  cpthetanhits->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamepthetanhits.Data()),"eps");
  cpthetanhits->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamepthetanhits.Data()),"root");

  TString cnamenhitsall = Form("cnhitsall_%s_%d_%d",particle.Data(),trunclevel,trunclevel2);
  TCanvas *cnhitsall = new TCanvas(cnamenhitsall,cnamenhitsall,700,500);
  cnhitsall->cd();
  hnhitsall->Draw();
  hnhitsall->Write();
  cnhitsall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.eps",particle.Data(),cnamenhitsall.Data()),"eps");
  cnhitsall->Print(Form("/data.local/nacer/halld_my/p2pi/fig_%s/%s.root",particle.Data(),cnamenhitsall.Data()),"root");

  outputfig->Print();
  // outputfig->Close();
}
