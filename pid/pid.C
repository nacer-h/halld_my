
#include <TROOT.h>
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

void pid()
{
  TFile *fpid = new TFile("fig/pid.root");
  TFile *fpid1 = new TFile("fig/pid1.root","UPDATE");

  TH2F* hpthetares[20];
  TH2F* hpthetasep[20];
  TH2F* hpthetamisid[20];
  TCanvas *cpthetaresall = new TCanvas("cpthetaresall","cpthetaresall",1500,1500);
  TCanvas *cpthetasepall = new TCanvas("cpthetasepall","cpthetasepall",1500,1500);
  TCanvas *cpthetamisidall = new TCanvas("cpthetamisidall","cpthetamisidall",1500,1500);

  cpthetaresall->Divide(5,4,0.00001,0.00001);
  cpthetasepall->Divide(5,4,0.00001,0.00001);
  cpthetamisidall->Divide(5,4,0.00001,0.00001);

  TH2F *hpthetabestres = new TH2F("hpthetabestres", "Momentum vs. #theta vs. optimal resolution; #theta [deg]; P [GeV/c]; resolution", 20,20,60,20,0.4,0.8);
  TH2F *hpthetabestsep = new TH2F("hpthetabestsep", "Momentum vs. #theta vs. optimal separation power; #theta [deg]; P [GeV/c]; separation power", 20,20,60,20,0.4,0.8);
  TH2F *hpthetabestmisid = new TH2F("hpthetabestmisid", "Momentum vs. #theta vs. optimal mis-id; #theta [deg]; P [GeV/c]; mis-id", 20,20,60,20,0.4,0.8);

  TH2F *hpthetatruncres = new TH2F("hpthetatruncres", "Momentum vs. #theta vs. optimal truncation (resolution); #theta [deg]; P [GeV/c]; truncation [%]", 20,20,60,20,0.4,0.8);
  TH2F *hpthetatruncsep = new TH2F("hpthetatruncsep", "Momentum vs. #theta vs. optimal truncation (separation power); #theta [deg]; P [GeV/c]; truncation [%]", 20,20,60,20,0.4,0.8);
  TH2F *hpthetatruncmisid = new TH2F("hpthetatruncmisid", "Momentum vs. #theta vs. optimal truncation (mis-id); #theta [deg]; P [GeV/c]; truncation [%]", 20,20,60,20,0.4,0.8);

  // TGraph *grestrunc = new TGraph("grestrunc","classifiers evolution vs. truncation for (GeV/c, °); truncation [%%]; #sigma/S/mis-id",20,0,100,20,0,0.2);
  TGraph *grestrunc = new TGraph(20);
  TGraph *gseptrunc = new TGraph(20);
  TGraph *gmisidtrunc = new TGraph(20);
  TMultiGraph *mg = new TMultiGraph();

  TH1F *htruncres = new TH1F("htruncres", "optimal truncation (based on resolution); Truncation [%]; Counts", 20,2,102);
  TH1F *htruncsep = new TH1F("htruncsep", "optimal truncation (based on separation power); Truncation [%]; Counts", 20,2,102);
  TH1F *htruncmisid = new TH1F("htruncmisid", "optimal truncation (based on mis-id); Truncation [%]; Counts", 20,2,102);

  // gPad->SetRightMargin(0.1461318);
  gStyle->SetOptStat(0);
  // gStyle->SetLabelSize(0.05);
  // gStyle->SetTitleSize(0.05);
  // gStyle->SetTitleOffset(0.95);

  int k = 0;
  for (int i=0;i<20;++i)
  {
    hpthetares[i] = (TH2F*)fpid->Get(Form("hpthetares_%d",k));
    cpthetaresall->cd(i+1);
    hpthetares[i]->SetTitle(Form("truncation = %d%%",(int)k));
    // hpthetares[i]->GetZaxis()->SetRangeUser(0,2);
    hpthetares[i]->Draw("colz");

    hpthetasep[i] = (TH2F*)fpid->Get(Form("hpthetasep_%d",k));
    cpthetasepall->cd(i+1);
    hpthetasep[i]->SetTitle(Form("truncation = %d%%",(int)k));
    // hpthetasep[i]->GetZaxis()->SetRangeUser(0,6);
    hpthetasep[i]->Draw("colz");

    hpthetamisid[i] = (TH2F*)fpid1->Get(Form("hpthetamisid_%d",k));
    cpthetamisidall->cd(i+1);
    hpthetamisid[i]->SetTitle(Form("truncation = %d%%",(int)k));
    // hpthetamisid[i]->GetZaxis()->SetRangeUser(0.000001,1);
    // gPad->SetLogz();
    hpthetamisid[i]->Draw("colz");

    k+=5;
  }
  cpthetaresall->Print("fig/cpthetaresall.eps","eps");
  cpthetasepall->Print("fig/cpthetasepall.eps","eps");
  cpthetamisidall->Print("fig/cpthetamisidall.eps","eps");
  cpthetaresall->Print("fig/cpthetaresall.root","root");
  cpthetasepall->Print("fig/cpthetasepall.root","root");
  cpthetamisidall->Print("fig/cpthetamisidall.root","root");

  for(int ip=0; ip<20; ++ip)
  {
    for(int itheta=0; itheta<20; ++itheta)
    {
      double res   = 1000;
      double sep = -1;
      double misid   = 1000;

      double truncres = 1000;
      double truncsep = 1000;
      double truncmisid = 1000;

      int k1 = 0;
      double k2 = 0;
      for (int i=0;i<20;i++)
      {
        if (hpthetares[i]->GetBinContent(itheta+1, ip+1) < res)
        {
          res = hpthetares[i]->GetBinContent(itheta+1, ip+1);
          truncres = k1;
          printf("****** ip= %d | itheta= %d | res= %f | truncres = %.0f ******\n", ip, itheta, res, truncres);
        }

        if (hpthetasep[i]->GetBinContent(itheta+1, ip+1) > sep)
        {
          sep = hpthetasep[i]->GetBinContent(itheta+1, ip+1);
          truncsep = k1;
          printf("------ ip= %d | itheta= %d | sep = %f | truncsep = %.0f ------\n", ip, itheta, sep, truncsep);
        }

        if (hpthetamisid[i]->GetBinContent(itheta+1, ip+1) < misid && hpthetamisid[i]->GetBinContent(itheta+1, ip+1)>1e-6)
        {
          misid = hpthetamisid[i]->GetBinContent(itheta+1, ip+1);
          truncmisid = k1;
          printf("****** ip= %d | itheta= %d | misid= %f | truncmisid = %.0f ******\n", ip, itheta, misid, truncmisid);
        }

        k1+=5;
        if(ip == 15 && itheta == 10)
        {
          grestrunc->SetPoint(i, k2, hpthetares[i]->GetBinContent(itheta+1, ip+1));
          gseptrunc->SetPoint(i, k2, hpthetasep[i]->GetBinContent(itheta+1, ip+1));
          gmisidtrunc->SetPoint(i, k2, hpthetamisid[i]->GetBinContent(itheta+1, ip+1)*1000);
          k2+=5;
        }
      }

      htruncres->Fill(truncres);
      htruncsep->Fill(truncsep);
      htruncmisid->Fill(truncmisid);

      hpthetabestres->SetBinContent(itheta+1, ip+1, res);
      hpthetatruncres->SetBinContent(itheta+1, ip+1, truncres);

      hpthetabestsep->SetBinContent(itheta+1, ip+1, sep);
      hpthetatruncsep->SetBinContent(itheta+1, ip+1, truncsep);

      hpthetabestmisid->SetBinContent(itheta+1, ip+1, misid);
      hpthetatruncmisid->SetBinContent(itheta+1, ip+1, truncmisid);

    }
  }

  TCanvas *ctruncres = new TCanvas("ctruncres","ctruncres",700,500);
  ctruncres->cd();
  htruncres->GetXaxis()->SetNdivisions(505);
  htruncres->GetXaxis()->SetLabelSize(0.05);
  htruncres->GetXaxis()->SetTitleSize(0.05);
  htruncres->GetXaxis()->SetTitleOffset(0.95);
  htruncres->GetYaxis()->SetLabelSize(0.05);
  htruncres->GetYaxis()->SetTitleSize(0.05);
  htruncres->GetYaxis()->SetTitleOffset(0.95);
  htruncres->Fit("gaus");
  htruncres->Draw();
  ctruncres->Print("fig/ctruncres.eps","eps");
  ctruncres->Print("fig/ctruncres.root","root");

  TCanvas *ctruncsep = new TCanvas("ctruncsep","ctruncsep",700,500);
  ctruncsep->cd();
  htruncsep->GetXaxis()->SetNdivisions(505);
  htruncsep->GetXaxis()->SetLabelSize(0.05);
  htruncsep->GetXaxis()->SetTitleSize(0.05);
  htruncsep->GetXaxis()->SetTitleOffset(0.95);
  htruncsep->GetYaxis()->SetLabelSize(0.05);
  htruncsep->GetYaxis()->SetTitleSize(0.05);
  htruncsep->GetYaxis()->SetTitleOffset(0.95);
  htruncsep->Fit("gaus");
  htruncsep->Draw();
  ctruncsep->Print("fig/ctruncsep.eps","eps");
  ctruncsep->Print("fig/ctruncsep.root","root");

  TCanvas *ctruncmisid = new TCanvas("ctruncmisid","ctruncmisid",700,500);
  ctruncmisid->cd();
  htruncmisid->GetXaxis()->SetNdivisions(505);
  htruncmisid->GetXaxis()->SetLabelSize(0.05);
  htruncmisid->GetXaxis()->SetTitleSize(0.05);
  htruncmisid->GetXaxis()->SetTitleOffset(0.95);
  htruncmisid->GetYaxis()->SetLabelSize(0.05);
  htruncmisid->GetYaxis()->SetTitleSize(0.05);
  htruncmisid->GetYaxis()->SetTitleOffset(0.95);
  htruncmisid->Fit("gaus");
  htruncmisid->Draw();
  ctruncmisid->Print("fig/ctruncmisid.eps","eps");
  ctruncmisid->Print("fig/ctruncmisid.root","root");

  TCanvas *ctrunc = new TCanvas("ctrunc","ctrunc",700,500);
  ctrunc->cd();
  gStyle->SetHistMinimumZero();
  TLegend *ltrunc = new TLegend(0.80,0.80,0.90,0.90);
  ltrunc->SetTextSize(0.04);
  ltrunc->SetBorderSize(0);
  // grestrunc->SetMarkerSize(2);
  grestrunc->SetMarkerStyle(20);
  grestrunc->SetMarkerColor(4);
  grestrunc->SetLineColor(4);
  // gseptrunc->SetMarkerSize(2);
  gseptrunc->SetMarkerStyle(21);
  gseptrunc->SetMarkerColor(2);
  gseptrunc->SetLineColor(2);
  // gmisidtrunc->SetMarkerSize(2);
  gmisidtrunc->SetMarkerStyle(22);
  gmisidtrunc->SetMarkerColor(1);
  gmisidtrunc->SetLineColor(1);
  mg->Add(grestrunc);
  mg->Add(gseptrunc);
  mg->Add(gmisidtrunc);
  mg->SetTitle("classifiers evolution vs. truncation; truncation [%]; Classifiers");
  mg->Draw("AP");
  ltrunc->AddEntry(grestrunc, "#sigma","p");
  ltrunc->AddEntry(gseptrunc, "S","p");
  ltrunc->AddEntry(gmisidtrunc, "mis-id","p");
  ltrunc->Draw();
  ctrunc->Modified();
  ctrunc->Print("fig/ctrunc.eps","eps");

  // TFile *outputfig = new TFile("fig/pid.root","RECREATE");
  TCanvas *cpthetabestres = new TCanvas("cpthetabestres","cpthetabestres",700,500);
  cpthetabestres->cd();
  gStyle->SetPalette(1);
  hpthetabestres->GetXaxis()->SetLabelSize(0.05);
  hpthetabestres->GetXaxis()->SetTitleSize(0.05);
  hpthetabestres->GetXaxis()->SetTitleOffset(0.95);
  hpthetabestres->GetYaxis()->SetLabelSize(0.05);
  hpthetabestres->GetYaxis()->SetTitleSize(0.05);
  hpthetabestres->GetYaxis()->SetTitleOffset(0.95);
  hpthetabestres->GetZaxis()->SetLabelSize(0.05);
  hpthetabestres->GetZaxis()->SetTitleSize(0.05);
  hpthetabestres->GetZaxis()->SetTitleOffset(0.95);
  // hpthetabestres->GetZaxis()->SetRangeUser(0,2);
  hpthetabestres->GetZaxis()->SetNdivisions(515);
  hpthetabestres->Draw("colz");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hpthetabestres->GetListOfFunctions()->FindObject("palette");
  cpthetabestres->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  hpthetabestres->Write();
  cpthetabestres->Modified();
  hpthetabestres->Write();
  cpthetabestres->Print("fig/cpthetabestres.eps","eps");
  cpthetabestres->Print("fig/cpthetabestres.root","root");

  TCanvas *cpthetatruncres = new TCanvas("cpthetatruncres","cpthetatruncres",700,500);
  cpthetatruncres->cd();
  hpthetatruncres->GetXaxis()->SetLabelSize(0.05);
  hpthetatruncres->GetXaxis()->SetTitleSize(0.05);
  hpthetatruncres->GetXaxis()->SetTitleOffset(0.95);
  hpthetatruncres->GetYaxis()->SetLabelSize(0.05);
  hpthetatruncres->GetYaxis()->SetTitleSize(0.05);
  hpthetatruncres->GetYaxis()->SetTitleOffset(0.95);
  hpthetatruncres->GetZaxis()->SetLabelSize(0.05);
  hpthetatruncres->GetZaxis()->SetTitleSize(0.05);
  hpthetatruncres->GetZaxis()->SetTitleOffset(0.95);
  // hpthetatruncres->GetZaxis()->SetRangeUser(0,2);
  hpthetatruncres->GetZaxis()->SetRangeUser(0,95);
  hpthetatruncres->GetZaxis()->SetNdivisions(515);
  hpthetatruncres->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetatruncres->GetListOfFunctions()->FindObject("palette");
  cpthetatruncres->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  hpthetatruncres->Write();
  cpthetatruncres->Print("fig/cpthetatruncres.eps","eps");
  cpthetatruncres->Print("fig/cpthetatruncres.root","root");

  TCanvas *cpthetabestsep = new TCanvas("cpthetabestsep","cpthetabestsep",700,500);
  cpthetabestsep->cd();
  gStyle->SetPalette(1);
  hpthetabestsep->SetTitleSize(0.05,"XYZ");
  hpthetabestsep->SetLabelSize(0.05,"XYZ");
  hpthetabestsep->SetTitleOffset(0.95,"XYZ");
  // hpthetabestsep->GetZaxis()->SetRangeUser(0,7);
  hpthetabestsep->GetZaxis()->SetNdivisions(515);
  hpthetabestsep->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetabestsep->GetListOfFunctions()->FindObject("palette");
  cpthetabestsep->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  hpthetabestsep->Write();
  cpthetabestsep->Modified();
  hpthetabestsep->Write();
  cpthetabestsep->Print("fig/cpthetabestsep.eps","eps");
  cpthetabestsep->Print("fig/cpthetabestsep.root","root");

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
  hpthetatruncsep->Write();
  cpthetatruncsep->Print("fig/cpthetatruncsep.eps","eps");
  cpthetatruncsep->Print("fig/cpthetatruncsep.root","root");

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
  gPad->SetLogz();
  hpthetabestmisid->GetZaxis()->SetNdivisions(515);
  hpthetabestmisid->Draw("colz");
  gPad->Update();
  palette = (TPaletteAxis*)hpthetabestmisid->GetListOfFunctions()->FindObject("palette");
  cpthetabestmisid->SetRightMargin(0.1461318);
  palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
  hpthetabestmisid->Write();
  cpthetabestmisid->Modified();
  hpthetabestmisid->Write();
  cpthetabestmisid->Print("fig/cpthetabestmisid.eps","eps");
  cpthetabestmisid->Print("fig/cpthetabestmisid.root","root");

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
  hpthetatruncmisid->Write();
  cpthetatruncmisid->Print("fig/cpthetatruncmisid.eps","eps");
  cpthetatruncmisid->Print("fig/cpthetatruncmisid.root","root");

}
