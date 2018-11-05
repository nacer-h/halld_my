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
  if(name == "sim") f = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_sim_recon17_1v03_mm2cut.root");
  if(name == "data") f = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_17v20_mm2cut.root");
  TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_y2175/y2175.root","UPDATE");

  //gStyle->SetLabelSize(0.1,"xyz"); // size of axis values
  //gStyle->SetTitleSize(0.1,"XYZ"); // size of axis titles
  //gStyle->SetTitleOffset(0.95,"X");

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
  c_PhiMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.root",name.Data()), "root");
  c_PhiMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.eps",name.Data()), "eps");

  TCanvas *c_PhiMass_beambunchcut = new TCanvas("c_PhiMass_beambunchcut","c_PhiMass_beambunchcut",600,400);
  c_PhiMass_beambunchcut->cd();
  h_PhiMass_beambunchcut->Rebin(6);
  h_PhiMass_beambunchcut->SetMarkerSize(0.7);
  h_PhiMass_beambunchcut->SetLineColor(kRed);
  h_PhiMass_beambunchcut->Draw();
  c_PhiMass_beambunchcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.root",name.Data()), "root");
  c_PhiMass_beambunchcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.eps",name.Data()), "eps");
  
}