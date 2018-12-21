#include <TStyle.h>
#include <TH2.h>
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
#include <TMultiGraph.h>
#include <Riostream.h>
#include <vector>
#include <string>
#include <fstream>
#include "TString.h"
#include "TPaletteAxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TProfile.h>
#include <TProfile2D.h>
#include <TPave.h>
#include <TBox.h>
#include <vector>
#include "TVectorD.h"

int prt_shiftHist(TH1F *hist, double double_shift)
{
  int bins=hist->GetXaxis()->GetNbins();
  double xmin=hist->GetXaxis()->GetBinLowEdge(1);
  double xmax=hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift=double_shift*(bins/(xmax-xmin));
  int shift=0;
  if(double_shift<0) shift=TMath::FloorNint(double_shift);
  if(double_shift>0) shift=TMath::CeilNint(double_shift);
  if(shift==0) return 0;
  if(shift>0){
    for(int i=1; i<=bins; i++){
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
    return 0;
  }
  if(shift<0){
    for(int i=bins; i>0; i--){
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }
    return 0;
  }
  return 1;
}

double getmisid(TF1* fsig, TF1 *fbkg, double &xcut, double eff=0.95)
{
  double xmin, xmax;
  fsig->GetRange(xmin,xmax);
  double integ = 0;
  double step = (xmax - xmin)/3;
  xcut = (xmax + xmin)/2;
  double normsig = fsig->Integral(xmin, xmax);
  if (normsig<1e-8) return 1.0;

  double diffinteg = 1;
  int maxloop = 100;
  while(diffinteg > 1e-5 && --maxloop>0)
  {
    //if(normsig < 1e-5) continue;
    //if(fsig->Integral(xcut, xmax) < 1e-5) continue;
    integ = fsig->Integral(xcut, xmax)/normsig;
    diffinteg = fabs(integ - eff);
    if(integ > eff) xcut += step;
    else xcut -= step;
    step /= 2;
  }

  cout <<"DIFFINTEG "<<diffinteg<<endl;

  // printf("****************** no problem up to here *****************\n");
  // if(fbkg->Integral(xmin,xmax) == 0) continue;
  double misid = fbkg->Integral(xcut,xmax)/fbkg->Integral(xmin,xmax);

  printf("integ = %f | misid = %f | xcut = %f xmin = %f | xmax = %f | normsig = %f\n", integ, misid,xcut,xmin,xmax,normsig);

  return misid;
}

void p2pimisid(int trunclevel, int trunclevel2, TString particle1, TString particle2)
{
  TFile *fpid = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root",particle1.Data(),particle1.Data()));
  TFile *fpid1 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root",particle2.Data(),particle2.Data()));
  TFile *outputfig = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_misid/pid%s%s.root",particle1.Data(),particle2.Data()),"UPDATE");  // file to save in all the plots.

  // gStyle->SetOptStat(0);

  TString hnamepthetamisid = Form("hpthetamisid_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TH2F *hpthetamisid = new TH2F(hnamepthetamisid, Form("Momentum vs. #theta vs. mis-id (%s , %s);#theta [deg]; P [GeV/c]; mis-id",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);

  TH2F *hptheta = (TH2F*)fpid->Get("hptheta");
  cout <<"  hptheta =  " <<hptheta<<endl;

  TString hnamemisid = Form("hmisid_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TH1F *hmisid = new TH1F(hnamemisid,Form(";Mis-ID (%s , %s)",particle1.Data(),particle2.Data()), 50,0,1);

  TString cnamededxexp = Form("cdedxexp_%d_%d",trunclevel,trunclevel2);
  TF1 *fdedxexp = (TF1*) fpid->Get(cnamededxexp);
  cout <<cnamededxexp<<" ***** fdedxexp = "<<fdedxexp<<endl;

  TString cnamededxmisid = Form("cdedxmisid_%d_%d",trunclevel,trunclevel2);
  TCanvas *cdedxmisid=new TCanvas(cnamededxmisid,cnamededxmisid,700,500);

  double misid = 1000;
  double sumbins = 0;

  for (int ip=0;ip<20;++ip)
  {
    for (int itheta=0;itheta<20;++itheta)
    {
      TString hname = Form("hddedxav_%s_%d_%d_%d_%d",particle1.Data(),ip,itheta,trunclevel,trunclevel2);
      TH1F *hddedxav = (TH1F*) fpid->Get(hname);
      double integ = hddedxav->Integral(1,100);
      if(integ < 300) continue;
      cout <<hname<<" ***** hddedxav = "<<hddedxav<<endl;

      TString hname1 = Form("hddedxav_%s_%d_%d_%d_%d",particle2.Data(),ip,itheta,trunclevel,trunclevel2);
      TH1F *hddedxav1 = (TH1F*) fpid1->Get(hname1);
      double integ1 = hddedxav1->Integral(1,100);
      if(integ1 < 300) continue;
      cout <<hname1<<" ***** hddedxav1 = "<<hddedxav1<<endl;

      double mom = hptheta->GetYaxis()->GetBinCenter(ip+1);
      double shiftproton = fdedxexp->Eval(mom);
      if(particle1 == "protonflat") prt_shiftHist(hddedxav,-1.*shiftproton);
      double xmin = hddedxav->GetXaxis()->GetXmin();
      double xmax = hddedxav->GetXaxis()->GetXmax();
      double xpeak  = hddedxav->GetBinCenter(hddedxav->GetMaximumBin());
      double xmin1 = hddedxav1->GetXaxis()->GetXmin();
      double xmax1 = hddedxav1->GetXaxis()->GetXmax();
      double xpeak1 = hddedxav1->GetBinCenter(hddedxav1->GetMaximumBin());

      TString gexp1 = "[0]*((((x-[1])/[2])<=-[3])*(exp(0.5*[3]^2+[3]*((x-[1])/[2]))) + (((x-[1])/[2])>-[3] && ((x-[1])/[2])<=[4])*(exp(-0.5*(x-[1])^2/[2]^2)) + (((x-[1])/[2])>[4])*(exp(0.5*[4]^2-[4]*((x-[1])/[2]))))";

      TF1 fsig("fsig",gexp1);
      fsig.SetLineColor(2); fsig.SetLineStyle(9);
      fsig.SetRange(xmin,xmax);
      fsig.SetParameters(hddedxav->GetMaximum(),xpeak,0.5,5,1);
      hddedxav->Fit("fsig","R","",xmin,xmax);
      fsig.SetParLimits(1,xpeak*0.9,xpeak*1.1);
      TF1 fbkg("fbkg",gexp1);
      fbkg.SetLineColor(4); fbkg.SetLineStyle(7);
      fbkg.SetRange(xmin1,xmax1);
      fbkg.SetParameters(hddedxav1->GetMaximum(),xpeak1,0.5,5,1);
      hddedxav1->Fit("fbkg","R","",xmin1,xmax1);
      fbkg.SetParLimits(1,xpeak1*0.9,xpeak1*1.1);

      printf("########################### shiftproton = %f | mean = %f | xpeak = %f \n",shiftproton,hddedxav->GetMean(),xpeak);

      double xcut = 0.0;
      misid = getmisid(&fsig, &fbkg, xcut);
      hpthetamisid->SetBinContent(itheta+1, ip+1, misid);

      TLine *lxcut = new TLine(xcut, hddedxav->GetMinimum(), xcut, hddedxav->GetMaximum());
      lxcut->SetLineWidth(2);
      lxcut->SetLineStyle(9);
      lxcut->SetLineColor(8);

      // plot average dE/dx for proton & pi+ (expGausexp fit)
      cdedxmisid->cd();
      hddedxav->SetTitle(Form("[ip=%i, itheta=%i]: p= %0.3f GeV/c, #theta= %0.2f^{0}, truncation low= %i%%, truncation high= %i%%, misid= %f; dE/dx [keV/cm]",ip,itheta,hptheta->GetYaxis()->GetBinCenter(ip+1),hptheta->GetXaxis()->GetBinCenter(itheta+1),trunclevel,trunclevel2,misid));
      // hddedxav->GetXaxis()->SetRangeUser(0,30);
      hddedxav->Draw();
      // hddedxav->Write();
      hddedxav->SetName(Form("hdedxmisid_%s_%d_%d_%d_%d",particle1.Data(),ip,itheta,trunclevel,trunclevel2));
      hddedxav->Write(Form("hdedxmisid_%s_%d_%d_%d_%d",particle1.Data(),ip,itheta,trunclevel,trunclevel2),TObject::kWriteDelete);
      hddedxav1->SetTitle(Form("[ip=%i, itheta=%i]: p= %0.3f GeV/c, #theta= %0.2f^{0}, truncation low= %i%%, truncation high= %i%%, misid= %f; dE/dx [keV/cm]",ip,itheta,hptheta->GetYaxis()->GetBinCenter(ip+1),hptheta->GetXaxis()->GetBinCenter(itheta+1),trunclevel,trunclevel2,misid));
      hddedxav1->Draw("same");
      // hddedxav1->Write();
      hddedxav1->SetName(Form("hdedxmisid_%s_%d_%d_%d_%d",particle2.Data(),ip,itheta,trunclevel,trunclevel2));
      hddedxav1->Write(Form("hdedxmisid_%s_%d_%d_%d_%d",particle2.Data(),ip,itheta,trunclevel,trunclevel2),TObject::kWriteDelete);
      lxcut->Draw();
      // lxcut->Write();
      // fsig.SetTitle(Form("[ip=%i, itheta=%i]: p= %0.3f GeV/c, #theta= %0.2f^{0}, truncation low= %i%%, truncation high= %i%%, misid= %f; dE/dx [keV/cm]",ip,itheta,hptheta->GetYaxis()->GetBinCenter(ip+1),hptheta->GetYaxis()->GetBinCenter(itheta+1),trunclevel,trunclevel2,misid));
      // fsig.Draw();
      // fbkg.Draw("same");
      cdedxmisid->Write();
      cdedxmisid->Update();

      double binentries = hptheta->GetBinContent(itheta+1, ip+1);
      sumbins += binentries;
      hmisid->Fill(misid,binentries);

      printf("particle1 = %s | particle2 = %s | trunclevel = %i | trunclevel2 = %i | ip = %i | itheta = %i | mom = %f | misid = %f\n",particle1.Data(), particle2.Data(), trunclevel, trunclevel2, ip, itheta,mom, misid);
    }
  }

  hmisid->Scale(1./sumbins);

  TString cnamemisid = Form("cmisid_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cmisid= new TCanvas(cnamemisid,cnamemisid,800,500);
  cmisid->cd();
  hmisid->Draw();
  hmisid->Write();
  cmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_misid/%s.eps",cnamemisid.Data()),"eps");
  cmisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_misid/%s.root",cnamemisid.Data()),"root");

  TString cnamepthetamisid = Form("cpthetamisid_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetamisid= new TCanvas(cnamepthetamisid,cnamepthetamisid,800,500);
  cpthetamisid->cd();
  gStyle->SetPalette(1);
  // hpthetamisid->GetZaxis()->SetRangeUser(0.00001,1);
  // if(trunclevel != 95) gPad->SetLogz();
  hpthetamisid->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpthetamisid->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetamisid->Modified();
  hpthetamisid->Write();
  cpthetamisid->Update();
  cpthetamisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_misid/%s.eps",cnamepthetamisid.Data()),"eps");
  cpthetamisid->Print(Form("/data.local/nacer/halld_my/p2pi/fig_misid/%s.root",cnamepthetamisid.Data()),"root");

  outputfig->Print();
  // outputfig->Close();
}
