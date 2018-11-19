//  File: y2175.C # Author: Nacer # Date: 2.9.2018 # Email: a.hamdi@gsi.de # Description: Macro to study Y(2175).                                                      
                                                                                                                                                                                   
#include "TCanvas.h"                                                                                                                                                               
#include "TFile.h"                                                                                                                                                                 
#include "TGraph.h"                                                                                                                                                                
#include "TString.h"                                                                                                                                                               
#include "TTree.h"                                                                                                                                                                 
#include "TH1.h"                                                                                                                                                                   
#include "TH2.h"                                                                                                                                                                   
#include "THStack.h"                                                                                                                                                               
#include "TMath.h"                                                                                                                                                                 
#include "TStyle.h"                                                                                                                                                                
#include "vector"                                                                                                                                                                  
#include "TF1.h"                                                                                                                                                                   
#include "TFormula.h"                                                                                                                                                              
#include "TLatex.h"                                                                                                                                                                
#include "TLegend.h"                                                                                                                                                               
#include "TLegendEntry.h"                                                                                                                                                          
#include "iostream"                                                                                                                                                                
#include "iomanip"                                                                                                                                                                 
#include "TLorentzVector.h"                                                                                                                                                        
#include "TGraphErrors.h"                                                                                                                                                          
#include "TAxis.h"                                                                                                                                                                 
#include "TROOT.h"                                                                                                                                                                 
#include "TSystem.h"                                                                                                                                                               
#include "TApplication.h"                                                                                                                                                          
#include "TPaveStats.h"                                                                                                                                                            
using namespace std;

void y2175(TString name)
{
  TFile *f = NULL;
  if(name == "sim") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_genr8_17v03_mm2cut.root");
  if(name == "data") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_17v21_mm2cut.root");
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/y2175.root","UPDATE");

  //gStyle->SetLabelSize(0.1,"xyz"); // size of axis values
  //gStyle->SetTitleSize(0.1,"XYZ"); // size of axis titles
  //gStyle->SetTitleOffset(0.95,"X");

  // ******* Tagger accidentals
  TCanvas *c_TaggerAccidentals = new TCanvas("c_TaggerAccidentals","c_TaggerAccidentals",600,400);
  c_TaggerAccidentals->cd();
  TH1F *h_TaggerAccidentals = (TH1F*) f->Get("h_TaggerAccidentals");
  h_TaggerAccidentals->Draw();
  c_TaggerAccidentals->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.root",name.Data()), "root");
  c_TaggerAccidentals->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.eps",name.Data()), "eps");

  TCanvas *c_TaggerAccidentals_postcut = new TCanvas("c_TaggerAccidentals_postcut","c_TaggerAccidentals_postcut",600,400);
  c_TaggerAccidentals_postcut->cd();
  TH1F *h_TaggerAccidentals_postcut = (TH1F*) f->Get("h_TaggerAccidentals_postcut");
  h_TaggerAccidentals_postcut->Draw();
  c_TaggerAccidentals_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.root",name.Data()), "root");
  c_TaggerAccidentals_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.eps",name.Data()), "eps");

    // ********  Beam Energy
  TCanvas *c_BeamEnergy = new TCanvas("c_BeamEnergy","c_BeamEnergy",600,400);
  c_BeamEnergy->cd();
  TH1F *h_BeamEnergy = (TH1F*) f->Get("BeamEnergy");
  h_BeamEnergy->Draw();
  c_BeamEnergy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_BeamEnergy.root",name.Data()), "root");
  c_BeamEnergy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_BeamEnergy.eps",name.Data()), "eps");

  // ********  Missing Mass against (Pi+, Pi-, Proton)
  TCanvas *c_kkmm2 = new TCanvas("c_kkmm2","c_kkmm2",1500,800);
  c_kkmm2->Divide(4,4);
  // c_kkmm2->Divide(1,3);
    for(int j=0; j<16; j++)
    {
      c_kkmm2->cd(j+1);
      TH1F *h_kkmm2_KinFit = (TH1F*) f->Get(Form("h_kkmm2_KinFit_%d",j));
      cout<<"h_kkmm2_KinFit = "<<h_kkmm2_KinFit<<endl;
      h_kkmm2_KinFit->SetLineColor(kBlue);
      TH1F *h_kkmm2_Measured = (TH1F*) f->Get(Form("h_kkmm2_Measured_%d",j));
      cout<<"h_kkmm2_Measured = "<<h_kkmm2_Measured<<endl;
      h_kkmm2_Measured->SetLineColor(kRed);
      h_kkmm2_KinFit->Draw();
      h_kkmm2_Measured->Draw("same");
    }

  c_kkmm2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm2.root",name.Data()), "root");
  c_kkmm2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm2.eps",name.Data()), "eps");

  // ******* Phi(1020) Mass: Measured & KinFit
  TCanvas *c_PhiMass = new TCanvas("c_PhiMass","c_PhiMass",600,400);
  c_PhiMass->cd();
  TLegend *l_PhiMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_PhiMass->SetTextSize(0.02);                                                                                                                                                 
  l_PhiMass->SetBorderSize(0); 
  c_PhiMass->SetGrid();
  TH1F *h_PhiMass_Measured = (TH1F*) f->Get("PhiMass_Measured");
  cout<<"h_PhiMass_Measured = "<<h_PhiMass_Measured<<endl;
  h_PhiMass_Measured->SetMinimum(0.);
  h_PhiMass_Measured->SetLineWidth(3);
  h_PhiMass_Measured->SetLineColor(kBlack);
  h_PhiMass_Measured->Draw();
  TH1F *h_PhiMass_KinFit = (TH1F*) f->Get("PhiMass_KinFit");
  cout<<"h_PhiMass_KinFit = "<<h_PhiMass_KinFit<<endl;
  h_PhiMass_KinFit->SetLineWidth(3);
  h_PhiMass_KinFit->SetLineColor(kBlue);
  h_PhiMass_KinFit->Draw("same");
  TH1F *h_PhiMass_beambunchcut = (TH1F*) f->Get("h_PhiMass_beambunchcut");
  cout<<"h_PhiMass_beambunchcut = "<<h_PhiMass_beambunchcut<<endl;
  h_PhiMass_beambunchcut->SetLineWidth(3);
  h_PhiMass_beambunchcut->SetLineColor(kRed);
  h_PhiMass_beambunchcut->Draw("same");
  l_PhiMass->AddEntry(h_PhiMass_Measured, "Measured","l"); 
  l_PhiMass->AddEntry(h_PhiMass_KinFit, "KinFit","l");
  l_PhiMass->AddEntry(h_PhiMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_PhiMass->Draw();
  c_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.root",name.Data()), "root");
  c_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  TCanvas *c_PhiMass_beambunchcut = new TCanvas("c_PhiMass_beambunchcut","c_PhiMass_beambunchcut",600,400);
  c_PhiMass_beambunchcut->cd();
  h_PhiMass_beambunchcut->Rebin(4);
  h_PhiMass_beambunchcut->SetMarkerSize(0.7);
  h_PhiMass_beambunchcut->SetLineColor(kRed);
  h_PhiMass_beambunchcut->Draw();
  c_PhiMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.root",name.Data()), "root");
  c_PhiMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.eps",name.Data()), "eps");

  TCanvas *c_PhiMass_nobeambunchcut = new TCanvas("c_PhiMass_nobeambunchcut","c_PhiMass_nobeambunchcut",600,400);
  c_PhiMass_nobeambunchcut->cd();
  TH1F *h_PhiMass_nobeambunchcut = (TH1F*) f->Get("h_PhiMass_nobeambunchcut");
  cout<<"h_PhiMass_nobeambunchcut = "<<h_PhiMass_nobeambunchcut<<endl;
  h_PhiMass_nobeambunchcut->Rebin(4);
  h_PhiMass_nobeambunchcut->SetMarkerSize(0.7);
  h_PhiMass_nobeambunchcut->Draw();
  c_PhiMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_nobeambunchcut.root",name.Data()), "root");
  c_PhiMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_nobeambunchcut.eps",name.Data()), "eps");
  
   // ******* fo(980) Mass: Measured & KinFit
  TCanvas *c_foMass = new TCanvas("c_foMass","c_foMass",600,400);
  c_foMass->cd();
  TLegend *l_foMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_foMass->SetTextSize(0.02);                                                                                                                                                 
  l_foMass->SetBorderSize(0); 
  c_foMass->SetGrid();
  TH1F *h_foMass_Measured = (TH1F*) f->Get("foMass_Measured");
  cout<<"h_foMass_Measured = "<<h_foMass_Measured<<endl;
  h_foMass_Measured->SetMinimum(0.);
  h_foMass_Measured->SetLineWidth(3);
  h_foMass_Measured->SetLineColor(kBlack);
  h_foMass_Measured->Draw();
  TH1F *h_foMass_KinFit = (TH1F*) f->Get("foMass_KinFit");
  cout<<"h_foMass_KinFit = "<<h_foMass_KinFit<<endl;
  h_foMass_KinFit->SetLineWidth(3);
  h_foMass_KinFit->SetLineColor(kBlue);
  h_foMass_KinFit->Draw("same");
  TH1F *h_foMass_beambunchcut = (TH1F*) f->Get("h_foMass_beambunchcut");
  cout<<"h_foMass_beambunchcut = "<<h_foMass_beambunchcut<<endl;
  h_foMass_beambunchcut->SetLineWidth(3);
  h_foMass_beambunchcut->SetLineColor(kRed);
  h_foMass_beambunchcut->Draw("same");
  l_foMass->AddEntry(h_foMass_Measured, "Measured","l"); 
  l_foMass->AddEntry(h_foMass_KinFit, "KinFit","l");
  l_foMass->AddEntry(h_foMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_foMass->Draw();
  c_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.root",name.Data()), "root");
  c_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  TCanvas *c_foMass_beambunchcut = new TCanvas("c_foMass_beambunchcut","c_foMass_beambunchcut",600,400);
  c_foMass_beambunchcut->cd();
  h_foMass_beambunchcut->Rebin(4);
  h_foMass_beambunchcut->SetMarkerSize(0.7);
  h_foMass_beambunchcut->SetLineColor(kRed);
  h_foMass_beambunchcut->Draw();
  c_foMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_beambunchcut.root",name.Data()), "root");
  c_foMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_beambunchcut.eps",name.Data()), "eps");

  TCanvas *c_foMass_nobeambunchcut = new TCanvas("c_foMass_nobeambunchcut","c_foMass_nobeambunchcut",600,400);
  c_foMass_nobeambunchcut->cd();
  TH1F *h_foMass_nobeambunchcut = (TH1F*) f->Get("h_foMass_nobeambunchcut");
  cout<<"h_foMass_nobeambunchcut = "<<h_foMass_nobeambunchcut<<endl;
  h_foMass_nobeambunchcut->Rebin(4);
  h_foMass_nobeambunchcut->SetMarkerSize(0.7);
  h_foMass_nobeambunchcut->Draw();
  c_foMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_nobeambunchcut.root",name.Data()), "root");
  c_foMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_nobeambunchcut.eps",name.Data()), "eps");

  // ******* Y(2175) Mass: Measured & KinFit
  TCanvas *c_YMass = new TCanvas("c_YMass","c_YMass",600,400);
  c_YMass->cd();
  TLegend *l_YMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_YMass->SetTextSize(0.02);                                                                                                                                                 
  l_YMass->SetBorderSize(0); 
  c_YMass->SetGrid();
  TH1F *h_YMass_Measured = (TH1F*) f->Get("YMass_Measured");
  cout<<"h_YMass_Measured = "<<h_YMass_Measured<<endl;
  h_YMass_Measured->SetMinimum(0.);
  h_YMass_Measured->SetLineWidth(3);
  h_YMass_Measured->SetLineColor(kBlack);
  TH1F *h_YMass_KinFit = (TH1F*) f->Get("YMass_KinFit");
  cout<<"h_YMass_KinFit = "<<h_YMass_KinFit<<endl;
  h_YMass_KinFit->SetLineWidth(3);
  h_YMass_KinFit->SetLineColor(kBlue);
  TH1F *h_YMass_beambunchcut = (TH1F*) f->Get("h_YMass_beambunchcut");
  cout<<"h_YMass_beambunchcut = "<<h_YMass_beambunchcut<<endl;
  h_YMass_beambunchcut->SetLineWidth(3);
  h_YMass_beambunchcut->SetLineColor(kRed);
  h_YMass_Measured->Draw();
  h_YMass_KinFit->Draw("same");
  h_YMass_beambunchcut->Draw("same");
  l_YMass->AddEntry(h_YMass_Measured, "Measured","l"); 
  l_YMass->AddEntry(h_YMass_KinFit, "KinFit","l");
  l_YMass->AddEntry(h_YMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_YMass->Draw();
  c_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.root",name.Data()), "root");
  c_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  TCanvas *c_YMass_beambunchcut = new TCanvas("c_YMass_beambunchcut","c_YMass_beambunchcut",600,400);
  c_YMass_beambunchcut->cd();
  h_YMass_beambunchcut->Rebin(4);
  h_YMass_beambunchcut->SetMarkerSize(0.7);
  h_YMass_beambunchcut->SetLineColor(kRed);
  h_YMass_beambunchcut->Draw();
  c_YMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_beambunchcut.root",name.Data()), "root");
  c_YMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_beambunchcut.eps",name.Data()), "eps");

  TCanvas *c_YMass_nobeambunchcut = new TCanvas("c_YMass_nobeambunchcut","c_YMass_nobeambunchcut",600,400);
  c_YMass_nobeambunchcut->cd();
  TH1F *h_YMass_nobeambunchcut = (TH1F*) f->Get("h_YMass_nobeambunchcut");
  cout<<"h_YMass_nobeambunchcut = "<<h_YMass_nobeambunchcut<<endl;
  h_YMass_nobeambunchcut->Rebin(4);
  h_YMass_nobeambunchcut->SetMarkerSize(0.7);
  h_YMass_nobeambunchcut->Draw();
  c_YMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_nobeambunchcut.root",name.Data()), "root");
  c_YMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_nobeambunchcut.eps",name.Data()), "eps");

}