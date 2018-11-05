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
  if(shift>0)
  {
    for(int i=1; i<=bins; i++)
    {
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
    return 0;
  }
  if(shift<0)
  {
    for(int i=bins; i>0; i--)
    {
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }
    return 0;
  }
  return 1;
}

void p2pisep(int trunclevel, int trunclevel2, TString particle1, TString particle2)
{
  TFile *fpid = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root",particle1.Data(),particle1.Data()));
  TFile *fpid1 = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_%s/pid%s.root",particle2.Data(),particle2.Data()));
  TFile *outputfig = new TFile(Form("/data.local/nacer/halld_my/p2pi/fig_sep/pid%s%s.root",particle1.Data(),particle2.Data()),"UPDATE");  // file to save in all the plots.

  // gStyle->SetOptStat(0);

  TString hnameresa=Form("hres_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TH1F *hresa = new TH1F(hnameresa,Form(";resolution [keV/cm] (%s)",particle1.Data()), 50,0,1);
  TString hnameresb = Form("hres_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TH1F *hresb = new TH1F(hnameresb,Form(";resolution [keV/cm] (%s)",particle2.Data()), 50,0,1);

  TString hnamesep=Form("hsep_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TH1F *hsep = new TH1F(hnamesep,Form(";separation power (%s , %s)",particle1.Data(),particle2.Data()), 50,0,20);

  TH2F *hptheta = (TH2F*)fpid->Get("hptheta");
  cout <<" ****** hptheta =  " <<hptheta<<endl;

  TString hnamepthetameana = Form("hpthetamean_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TH2F *hpthetameana = new TH2F(hnamepthetameana,Form("Momentum vs. #theta vs. mean (%s);#theta [deg]; P [GeV/c]; #mu [keV/cm]",particle1.Data()), 20,60,80,20,0.3,1);

  TString hnamepthetameanb = Form("hpthetamean_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TH2F *hpthetameanb = new TH2F(hnamepthetameanb, Form("Momentum vs. #theta vs. mean (%s);#theta [deg]; P [GeV/c]; #mu [keV/cm]",particle2.Data()), 20,60,80,20,0.3,1);

  TString hnamepthetaresa = Form("hpthetares_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TH2F *hpthetaresa = new TH2F(hnamepthetaresa,Form("Momentum vs. #theta vs. resolution (%s);#theta [deg]; P [GeV/c]; #sigma [keV/cm]",particle1.Data()), 20,60,80,20,0.3,1);

  TString hnamepthetaresb = Form("hpthetares_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TH2F *hpthetaresb = new TH2F(hnamepthetaresb, Form("Momentum vs. #theta vs. resolution (%s);#theta [deg]; P [GeV/c]; #sigma [keV/cm]",particle2.Data()), 20,60,80,20,0.3,1);

  TString hnamepthetasep = Form("hpthetasep_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TH2F *hpthetasep = new TH2F(hnamepthetasep, Form("Momentum vs. #theta vs. separation power (%s , %s);#theta [deg]; P [GeV/c]; separation power",particle1.Data(),particle2.Data()), 20,60,80,20,0.3,1);

  TString cnamededxexp = Form("cdedxexp_%d_%d",trunclevel,trunclevel2);
  TF1 *fdedxexp = (TF1*) fpid->Get(cnamededxexp);
  cout <<cnamededxexp<<" ***** fdedxexp = "<<fdedxexp<<endl;

  TString cnamededxsep = Form("cdedxsep_%d_%d",trunclevel,trunclevel2);
  TCanvas *cdedxsep=new TCanvas(cnamededxsep,cnamededxsep,700,500);

  double mu1, mu2, sigma1, sigma2, sep;
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

      TF1 *fdedxsig  = new TF1("fdedxsig","gaus",xpeak-hddedxav->GetRMS()*1.5,xpeak+hddedxav->GetRMS()*1.5);
      fdedxsig->SetLineColor(2);
      hddedxav->Fit(fdedxsig,"R","",xpeak-hddedxav->GetRMS()*1.5,xpeak+hddedxav->GetRMS()*1.5);
      TF1 *fdedxsig1  = new TF1("fdedxsig1","gaus",xpeak1-hddedxav1->GetRMS()*1.5,xpeak1+hddedxav1->GetRMS()*1.5);
      fdedxsig1->SetLineColor(4);
      hddedxav1->Fit(fdedxsig1,"R","",xpeak1-hddedxav1->GetRMS()*1.5,xpeak1+hddedxav1->GetRMS()*1.5);

      mu1 = fdedxsig->GetParameter(1);
      mu2 = fdedxsig1->GetParameter(1);
      sigma1 = fdedxsig->GetParameter(2);
      sigma2 = fdedxsig1->GetParameter(2);

      hpthetaresa->SetBinContent(itheta+1, ip+1, sigma1);
      hpthetaresb->SetBinContent(itheta+1, ip+1, sigma2);
      hpthetameana->SetBinContent(itheta+1, ip+1, mu1);
      hpthetameanb->SetBinContent(itheta+1, ip+1, mu2);
      sep= abs(mu1-mu2)/(sqrt(sigma1*sigma1+sigma2*sigma2));
      hpthetasep->SetBinContent(itheta+1, ip+1, sep);

      printf("########################### shiftproton = %f | hddedxav->GetMean() = %f | mu1 = %f\n",shiftproton,hddedxav->GetMean(),mu1);

      // plot average dE/dx for proton & pi+ (Gaus fit)
      cdedxsep->cd();
      hddedxav->SetTitle(Form("(ip=%i, itheta=%i): p= %0.3f GeV/c, #theta= %0.2f^{0}, truncation low= %i%%, truncation high= %i%%, sep= %f; dE/dx [keV/cm]",ip,itheta,hptheta->GetYaxis()->GetBinCenter(ip+1),hptheta->GetXaxis()->GetBinCenter(itheta+1),trunclevel,trunclevel2,sep));
      hddedxav->GetXaxis()->SetRangeUser(0,30);
      hddedxav->Draw();
      // hddedxav->Write();
      hddedxav->SetName(Form("hdedxsep_%s_%d_%d_%d_%d",particle1.Data(),ip,itheta,trunclevel,trunclevel2));
      hddedxav->Write(Form("hdedxsep_%s_%d_%d_%d_%d",particle1.Data(),ip,itheta,trunclevel,trunclevel2),TObject::kWriteDelete);
      hddedxav1->SetTitle(Form("(ip=%i, itheta=%i): p= %0.3f GeV/c, #theta= %0.2f^{0}, truncation low= %i%%, truncation high= %i%%, sep= %f; dE/dx [keV/cm]",ip,itheta,hptheta->GetYaxis()->GetBinCenter(ip+1),hptheta->GetXaxis()->GetBinCenter(itheta+1),trunclevel,trunclevel2,sep));
      hddedxav1->Draw("same");
      // hddedxav1->Write();
      hddedxav1->SetName(Form("hdedxsep_%s_%d_%d_%d_%d",particle2.Data(),ip,itheta,trunclevel,trunclevel2));
      hddedxav1->Write(Form("hdedxsep_%s_%d_%d_%d_%d",particle2.Data(),ip,itheta,trunclevel,trunclevel2),TObject::kWriteDelete);
      cdedxsep->Write();
      cdedxsep->Update();

      double binentries = hptheta->GetBinContent(itheta+1, ip+1);
      sumbins += binentries;
      hresa->Fill(sigma1,binentries);
      hresb->Fill(sigma2,binentries);
      hsep->Fill(sep,binentries);

      printf("particle1 = %s | particle2 = %s | trunclevel = %i | trunclevel2 = %i | ip = %i | itheta = %i | mu1 = %f | mu2 = %f | sigma1 = %f | sigma2 = %f | mom = %f | sep = %f\n",
      particle1.Data(), particle2.Data(), trunclevel, trunclevel2, ip, itheta, mu1, mu2, sigma1, sigma2, mom, sep);
    }
  }

  hresa->Scale(1./sumbins);
  hresb->Scale(1./sumbins);
  hsep->Scale(1./sumbins);

  TString cnameresa = Form("cres_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TCanvas *cresa= new TCanvas(cnameresa,cnameresa,800,500);
  cresa->cd();
  hresa->Draw();
  hresa->Write();
  cresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnameresa.Data()),"eps");
  cresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnameresa.Data()),"root");

  TString cnameresb = Form("cres_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cresb= new TCanvas(cnameresb,cnameresb,800,500);
  cresb->cd();
  hresb->Draw();
  hresb->Write();
  cresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnameresb.Data()),"eps");
  cresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnameresb.Data()),"root");

  TString cnamesep = Form("csep_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TCanvas *csep= new TCanvas(cnamesep,cnamesep,800,500);
  csep->cd();
  hsep->Write();
  hsep->Draw();
  csep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_sep/%s.eps",cnamesep.Data()),"eps");
  csep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_sep/%s.root",cnamesep.Data()),"root");

  TString cnamepthetameana = Form("cpthetamean_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetameana= new TCanvas(cnamepthetameana,cnamepthetameana,800,500);
  cpthetameana->cd();
  gStyle->SetPalette(1);
  // hpthetameana->GetZaxis()->SetRangeUser(0,7);
  hpthetameana->Draw("colz");
  hpthetameana->Write();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpthetameana->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetameana->Modified();
  hpthetameana->Write();
  cpthetameana->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnamepthetameana.Data()),"eps");
  cpthetameana->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnamepthetameana.Data()),"root");

  TString cnamepthetameanb = Form("cpthetamean_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetameanb= new TCanvas(cnamepthetameanb,cnamepthetameanb,800,500);
  cpthetameanb->cd();
  gStyle->SetPalette(1);
  // hpthetameanb->GetZaxis()->SetRangeUser(0,7);
  hpthetameanb->Draw("colz");
  hpthetameanb->Write();
  gPad->Update();
  palette = (TPaletteAxis*)hpthetameanb->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetameanb->Modified();
  hpthetameanb->Write();
  cpthetameanb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnamepthetameanb.Data()),"eps");
  cpthetameanb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnamepthetameanb.Data()),"root");

  TString cnamepthetaresa = Form("cpthetares_%s_%d_%d",particle1.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetaresa= new TCanvas(cnamepthetaresa,cnamepthetaresa,800,500);
  cpthetaresa->cd();
  gStyle->SetPalette(1);
  // hpthetaresa->GetZaxis()->SetRangeUser(0,7);
  hpthetaresa->Draw("colz");
  hpthetaresa->Write();
  gPad->Update();
  palette = (TPaletteAxis*)hpthetaresa->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetaresa->Modified();
  hpthetaresa->Write();
  cpthetaresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnamepthetaresa.Data()),"eps");
  cpthetaresa->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnamepthetaresa.Data()),"root");

  TString cnamepthetaresb = Form("cpthetares_%s_%d_%d",particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetaresb= new TCanvas(cnamepthetaresb,cnamepthetaresb,800,500);
  cpthetaresb->cd();
  gStyle->SetPalette(1);
  // hpthetaresb->GetZaxis()->SetRangeUser(0,7);
  hpthetaresb->Draw("colz");
  hpthetaresb->Write();
  gPad->Update();
  palette = (TPaletteAxis*)hpthetaresb->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetaresb->Modified();
  hpthetaresb->Write();
  cpthetaresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.eps",cnamepthetaresb.Data()),"eps");
  cpthetaresb->Print(Form("/data.local/nacer/halld_my/p2pi/fig_res/%s.root",cnamepthetaresb.Data()),"root");

  TString cnamepthetasep = Form("cpthetasep_%s_%s_%d_%d",particle1.Data(),particle2.Data(),trunclevel,trunclevel2);
  TCanvas *cpthetasep= new TCanvas(cnamepthetasep,cnamepthetasep,800,500);
  cpthetasep->cd();
  gStyle->SetPalette(1);
  // hpthetasep->GetZaxis()->SetRangeUser(0,7);
  hpthetasep->Draw("colz");
  hpthetasep->Write();
  gPad->Update();
  palette = (TPaletteAxis*)hpthetasep->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  cpthetasep->Modified();
  hpthetasep->Write();
  cpthetasep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_sep/%s.eps",cnamepthetasep.Data()),"eps");
  cpthetasep->Print(Form("/data.local/nacer/halld_my/p2pi/fig_sep/%s.root",cnamepthetasep.Data()),"root");

  outputfig->Print();
  // outputfig->Close();
}
