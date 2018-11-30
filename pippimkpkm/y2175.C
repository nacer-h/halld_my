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

void y2175(TString name, TString cut)
{
  TFile *f = NULL;
  if(name == "sim") f = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_sig_17v03_%scut.root",cut.Data()));
  if(name == "data") f = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_17v21_%scut.root",cut.Data()));
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/y2175.root","UPDATE");

  //gStyle->SetLabelSize(0.1,"xyz"); // size of axis values
  //gStyle->SetTitleSize(0.1,"XYZ"); // size of axis titles
  //gStyle->SetTitleOffset(0.95,"X");

  // ******* Tagger accidentals
  TCanvas *c_TaggerAccidentals = new TCanvas("c_TaggerAccidentals","c_TaggerAccidentals",600,400);
  c_TaggerAccidentals->cd();
  TH1F *h_TaggerAccidentals = (TH1F*) f->Get("h_TaggerAccidentals");
  h_TaggerAccidentals->Draw();
  TLine *l1_TaggerAccidentals = new TLine(0.5*4.008, 0, 0.5*4.008, 3500000);
  l1_TaggerAccidentals->SetLineColor(kRed);
  l1_TaggerAccidentals->SetLineStyle(2);
  l1_TaggerAccidentals->SetLineWidth(2);
  l1_TaggerAccidentals->Draw();
  TLine *l2_TaggerAccidentals = new TLine(-0.5*4.008, 0, -0.5*4.008, 3500000);
  l2_TaggerAccidentals->SetLineColor(kRed);
  l2_TaggerAccidentals->SetLineStyle(2);
  l2_TaggerAccidentals->SetLineWidth(2);
  l2_TaggerAccidentals->Draw();
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
  TCanvas *c_kkmm = new TCanvas("c_kkmm","c_kkmm",1500,800);
  c_kkmm->Divide(4,4);
  // c_kkmm2->Divide(1,3);
    for(int j=0; j<16; j++)
    {
      c_kkmm->cd(j+1);
      TH1F *h_kkmm_KinFit = (TH1F*) f->Get(Form("h_kkmm_KinFit_%d",j));
      cout<<"h_kkmm_KinFit = "<<h_kkmm_KinFit<<endl;
      h_kkmm_KinFit->SetLineColor(kBlue);
      TH1F *h_kkmm_Measured = (TH1F*) f->Get(Form("h_kkmm_Measured_%d",j));
      cout<<"h_kkmm_Measured = "<<h_kkmm_Measured<<endl;
      h_kkmm_Measured->SetLineColor(kRed);
      h_kkmm_KinFit->Draw();
      h_kkmm_Measured->Draw("same");
    }

  c_kkmm->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm.root",name.Data()), "root");
  c_kkmm->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm.eps",name.Data()), "eps");

  if(name=="sim")
  {
    TCanvas *c_kkmm_Thrown = new TCanvas("c_kkmm_Thrown","c_kkmm_Thrown",600,400);
    c_kkmm_Thrown->cd();
    TH1F *h_kkmm_Thrown = (TH1F*) f->Get("h_kkmm_Thrown");
    cout<<"h_kkmm_Thrown = "<<h_kkmm_Thrown<<endl;
    h_kkmm_Thrown->Draw();
    c_kkmm_Thrown->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm_Thrown.root",name.Data()), "root");
    c_kkmm_Thrown->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kkmm_Thrown.eps",name.Data()), "eps");
  }

 
  // ********************************* Phi(1020) Mass: Measured & KinFit
  // Precut
  TCanvas *c_PhiMass_precut = new TCanvas("c_PhiMass_precut","c_PhiMass_precut",600,400);
  c_PhiMass_precut->cd();
  TLegend *l_PhiMass_precut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_PhiMass_precut->SetTextSize(0.02);                                                                                                                                                 
  l_PhiMass_precut->SetBorderSize(0); 
  c_PhiMass_precut->SetGrid();
  TH1F *h_PhiMass_Measured_precut = (TH1F*) f->Get("PhiMass_Measured");
  cout<<"h_PhiMass_Measured_precut = "<<h_PhiMass_Measured_precut<<endl;
  h_PhiMass_Measured_precut->SetMinimum(0.);
  h_PhiMass_Measured_precut->SetLineWidth(3);
  h_PhiMass_Measured_precut->SetLineColor(kBlack);
  TH1F *h_PhiMass_KinFit_precut = (TH1F*) f->Get("PhiMass_KinFit");
  cout<<"h_PhiMass_KinFit_precut = "<<h_PhiMass_KinFit_precut<<endl;
  h_PhiMass_KinFit_precut->SetLineWidth(3);
  h_PhiMass_KinFit_precut->SetLineColor(kBlue);
  h_PhiMass_Measured_precut->Draw();
  h_PhiMass_KinFit_precut->Draw("same");
  // TH1F *h_PhiMass_beambunchcut = (TH1F*) f->Get("h_PhiMass_beambunchcut");
  // cout<<"h_PhiMass_beambunchcut = "<<h_PhiMass_beambunchcut<<endl;
  // h_PhiMass_beambunchcut->SetLineWidth(3);
  // h_PhiMass_beambunchcut->SetLineColor(kRed);
  // h_PhiMass_beambunchcut->Draw("same");
  l_PhiMass_precut->AddEntry(h_PhiMass_Measured_precut, "Measured","l"); 
  l_PhiMass_precut->AddEntry(h_PhiMass_KinFit_precut, "KinFit","l");
  // l_PhiMass->AddEntry(h_PhiMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_PhiMass_precut->Draw();
  c_PhiMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_precut.root",name.Data()), "root");
  c_PhiMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_precut.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  // TCanvas *c_PhiMass_beambunchcut = new TCanvas("c_PhiMass_beambunchcut","c_PhiMass_beambunchcut",600,400);
  // c_PhiMass_beambunchcut->cd();
  // h_PhiMass_beambunchcut->Rebin(4);
  // h_PhiMass_beambunchcut->SetMarkerSize(0.7);
  // h_PhiMass_beambunchcut->SetLineColor(kRed);
  // h_PhiMass_beambunchcut->Draw();
  // c_PhiMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.root",name.Data()), "root");
  // c_PhiMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_beambunchcut.eps",name.Data()), "eps");

  // Postcut
  TCanvas *c_PhiMass_postcut = new TCanvas("c_PhiMass_postcut","c_PhiMass_postcut",600,400);
  c_PhiMass_postcut->cd();
  TLegend *l_PhiMass_postcut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_PhiMass_postcut->SetTextSize(0.02);                                                                                                                                                 
  l_PhiMass_postcut->SetBorderSize(0); 
  c_PhiMass_postcut->SetGrid();
  TH1F *h_PhiMass_Measured_postcut = (TH1F*) f->Get("h_PhiMass_chi2cut_Measured");
  cout<<"h_PhiMass_Measured_postcut = "<<h_PhiMass_Measured_postcut<<endl;
  h_PhiMass_Measured_postcut->SetMinimum(0.);
  h_PhiMass_Measured_postcut->SetLineWidth(3);
  h_PhiMass_Measured_postcut->SetLineColor(kBlack);
  TH1F *h_PhiMass_KinFit_postcut = (TH1F*) f->Get("h_PhiMass_chi2cut_0");
  cout<<"h_PhiMass_KinFit_postcut = "<<h_PhiMass_KinFit_postcut<<endl;
  h_PhiMass_KinFit_postcut->SetLineWidth(3);
  h_PhiMass_KinFit_postcut->SetLineColor(kBlue);
  h_PhiMass_KinFit_postcut->Draw();
  h_PhiMass_Measured_postcut->Draw("same");
  l_PhiMass_postcut->AddEntry(h_PhiMass_Measured_postcut, "Measured","l"); 
  l_PhiMass_postcut->AddEntry(h_PhiMass_KinFit_postcut, "KinFit","l");
  l_PhiMass_postcut->Draw();
  c_PhiMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcut.root",name.Data()), "root");
  c_PhiMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcut.eps",name.Data()), "eps");

  // TCanvas *c_PhiMass_nobeambunchcut = new TCanvas("c_PhiMass_nobeambunchcut","c_PhiMass_nobeambunchcut",600,400);
  // c_PhiMass_nobeambunchcut->cd();
  // TH1F *h_PhiMass_nobeambunchcut = (TH1F*) f->Get(Form("h_PhiMass_%scut_nobeambunchcut",cut.Data()));
  // cout<<"h_PhiMass_nobeambunchcut = "<<h_PhiMass_nobeambunchcut<<endl;
  // // h_PhiMass_nobeambunchcut->Rebin(4);
  // h_PhiMass_nobeambunchcut->SetMarkerSize(0.7);
  // h_PhiMass_nobeambunchcut->Draw();
  // c_PhiMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_%scut_nobeambunchcut.root",name.Data(),cut.Data()), "root");
  // c_PhiMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_%scut_nobeambunchcut.eps",name.Data(),cut.Data()), "eps");
  
   // ******* fo(980) Mass: Measured & KinFit
  TCanvas *c_foMass_precut = new TCanvas("c_foMass_precut","c_foMass_precut",600,400);
  c_foMass_precut->cd();
  TLegend *l_foMass_precut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_foMass_precut->SetTextSize(0.02);                                                                                                                                                 
  l_foMass_precut->SetBorderSize(0); 
  c_foMass_precut->SetGrid();
  TH1F *h_foMass_Measured_precut = (TH1F*) f->Get("foMass_Measured");
  cout<<"h_foMass_Measured_precut = "<<h_foMass_Measured_precut<<endl;
  h_foMass_Measured_precut->SetMinimum(0.);
  h_foMass_Measured_precut->SetLineWidth(3);
  h_foMass_Measured_precut->SetLineColor(kBlack);
  TH1F *h_foMass_KinFit_precut = (TH1F*) f->Get("foMass_KinFit");
  cout<<"h_foMass_KinFit_precut = "<<h_foMass_KinFit_precut<<endl;
  h_foMass_KinFit_precut->SetLineWidth(3);
  h_foMass_KinFit_precut->SetLineColor(kBlue);
  h_foMass_KinFit_precut->Draw();
  h_foMass_Measured_precut->Draw("same");
  // TH1F *h_foMass_beambunchcut = (TH1F*) f->Get("h_foMass_beambunchcut");
  // cout<<"h_foMass_beambunchcut = "<<h_foMass_beambunchcut<<endl;
  // h_foMass_beambunchcut->SetLineWidth(3);
  // h_foMass_beambunchcut->SetLineColor(kRed);
  // h_foMass_beambunchcut->Draw("same");
  l_foMass_precut->AddEntry(h_foMass_Measured_precut, "Measured","l"); 
  l_foMass_precut->AddEntry(h_foMass_KinFit_precut, "KinFit","l");
  // l_foMass->AddEntry(h_foMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_foMass_precut->Draw();
  c_foMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_precut.root",name.Data()), "root");
  c_foMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_precut.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  // TCanvas *c_foMass_beambunchcut = new TCanvas("c_foMass_beambunchcut","c_foMass_beambunchcut",600,400);
  // c_foMass_beambunchcut->cd();
  // h_foMass_beambunchcut->Rebin(4);
  // h_foMass_beambunchcut->SetMarkerSize(0.7);
  // h_foMass_beambunchcut->SetLineColor(kRed);
  // h_foMass_beambunchcut->Draw();
  // c_foMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_beambunchcut.root",name.Data()), "root");
  // c_foMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_beambunchcut.eps",name.Data()), "eps");

  // Postcut
  TCanvas *c_foMass_postcut = new TCanvas("c_foMass_postcut","c_foMass_postcut",600,400);
  c_foMass_postcut->cd();
  TLegend *l_foMass_postcut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_foMass_postcut->SetTextSize(0.02);                                                                                                                                                 
  l_foMass_postcut->SetBorderSize(0); 
  c_foMass_postcut->SetGrid();
  TH1F *h_foMass_Measured_postcut = (TH1F*) f->Get("h_foMass_chi2cut_Measured");
  cout<<"h_foMass_Measured_postcut = "<<h_foMass_Measured_postcut<<endl;
  h_foMass_Measured_postcut->SetMinimum(0.);
  h_foMass_Measured_postcut->SetLineWidth(3);
  h_foMass_Measured_postcut->SetLineColor(kBlack);
  TH1F *h_foMass_KinFit_postcut = (TH1F*) f->Get("h_foMass_chi2cut_0");
  cout<<"h_foMass_KinFit_postcut = "<<h_foMass_KinFit_postcut<<endl;
  h_foMass_KinFit_postcut->SetLineWidth(3);
  h_foMass_KinFit_postcut->SetLineColor(kBlue);
  h_foMass_KinFit_postcut->Draw();
  h_foMass_Measured_postcut->Draw("same");
  l_foMass_postcut->AddEntry(h_foMass_Measured_postcut, "Measured","l"); 
  l_foMass_postcut->AddEntry(h_foMass_KinFit_postcut, "KinFit","l");
  l_foMass_postcut->Draw();
  c_foMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcut.root",name.Data()), "root");
  c_foMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcut.eps",name.Data()), "eps");

  // TCanvas *c_foMass_nobeambunchcut = new TCanvas("c_foMass_nobeambunchcut","c_foMass_nobeambunchcut",600,400);
  // c_foMass_nobeambunchcut->cd();
  // TH1F *h_foMass_nobeambunchcut = (TH1F*) f->Get(Form("h_foMass_%scut_nobeambunchcut",cut.Data()));
  // cout<<"h_foMass_nobeambunchcut = "<<h_foMass_nobeambunchcut<<endl;
  // // h_foMass_nobeambunchcut->Rebin(4);
  // h_foMass_nobeambunchcut->SetMarkerSize(0.7);
  // h_foMass_nobeambunchcut->Draw();
  // c_foMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_%scut_nobeambunchcut.root",name.Data(),cut.Data()), "root");
  // c_foMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_%scut_nobeambunchcut.eps",name.Data(),cut.Data()), "eps");

  // ******* Y(2175) Mass: Measured & KinFit
  TCanvas *c_YMass_precut = new TCanvas("c_YMass_precut","c_YMass_precut",600,400);
  c_YMass_precut->cd();
  TLegend *l_YMass_precut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_YMass_precut->SetTextSize(0.02);                                                                                                                                                 
  l_YMass_precut->SetBorderSize(0); 
  c_YMass_precut->SetGrid();
  TH1F *h_YMass_Measured_precut = (TH1F*) f->Get("YMass_Measured");
  cout<<"h_YMass_Measured_precut = "<<h_YMass_Measured_precut<<endl;
  h_YMass_Measured_precut->SetMinimum(0.);
  h_YMass_Measured_precut->SetLineWidth(3);
  h_YMass_Measured_precut->SetLineColor(kBlack);
  TH1F *h_YMass_KinFit_precut = (TH1F*) f->Get("YMass_KinFit");
  cout<<"h_YMass_KinFit_precut = "<<h_YMass_KinFit_precut<<endl;
  h_YMass_KinFit_precut->SetLineWidth(3);
  h_YMass_KinFit_precut->SetLineColor(kBlue);
  // TH1F *h_YMass_beambunchcut = (TH1F*) f->Get("h_YMass_beambunchcut");
  // cout<<"h_YMass_beambunchcut = "<<h_YMass_beambunchcut<<endl;
  // h_YMass_beambunchcut->SetLineWidth(3);
  // h_YMass_beambunchcut->SetLineColor(kRed);
  h_YMass_KinFit_precut->Draw();
  h_YMass_Measured_precut->Draw("same");
  // h_YMass_beambunchcut->Draw("same");
  l_YMass_precut->AddEntry(h_YMass_Measured_precut, "Measured","l"); 
  l_YMass_precut->AddEntry(h_YMass_KinFit_precut, "KinFit","l");
  // l_YMass->AddEntry(h_YMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_YMass_precut->Draw();
  c_YMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_precut.root",name.Data()), "root");
  c_YMass_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_precut.eps",name.Data()), "eps");
  // Zoom in after beam bunch cut
  // TCanvas *c_YMass_beambunchcut = new TCanvas("c_YMass_beambunchcut","c_YMass_beambunchcut",600,400);
  // c_YMass_beambunchcut->cd();
  // h_YMass_beambunchcut->Rebin(4);
  // h_YMass_beambunchcut->SetMarkerSize(0.7);
  // h_YMass_beambunchcut->SetLineColor(kRed);
  // h_YMass_beambunchcut->Draw();
  // c_YMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_beambunchcut.root",name.Data()), "root");
  // c_YMass_beambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_beambunchcut.eps",name.Data()), "eps");

  // Postcut
  TCanvas *c_YMass_postcut = new TCanvas("c_YMass_postcut","c_YMass_postcut",600,400);
  c_YMass_postcut->cd();
  TLegend *l_YMass_postcut = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_YMass_postcut->SetTextSize(0.02);                                                                                                                                                 
  l_YMass_postcut->SetBorderSize(0); 
  c_YMass_postcut->SetGrid();
  TH1F *h_YMass_Measured_postcut = (TH1F*) f->Get("h_YMass_chi2cut_Measured");
  cout<<"h_YMass_Measured_postcut = "<<h_YMass_Measured_postcut<<endl;
  h_YMass_Measured_postcut->SetMinimum(0.);
  h_YMass_Measured_postcut->SetLineWidth(3);
  h_YMass_Measured_postcut->SetLineColor(kBlack);
  TH1F *h_YMass_KinFit_postcut = (TH1F*) f->Get("h_YMass_chi2cut_0");
  cout<<"h_YMass_KinFit_postcut = "<<h_YMass_KinFit_postcut<<endl;
  h_YMass_KinFit_postcut->SetLineWidth(3);
  h_YMass_KinFit_postcut->SetLineColor(kBlue);
  h_YMass_KinFit_postcut->Draw();
  h_YMass_Measured_postcut->Draw("same");
  l_YMass_postcut->AddEntry(h_YMass_Measured_postcut, "Measured","l"); 
  l_YMass_postcut->AddEntry(h_YMass_KinFit_postcut, "KinFit","l");
  l_YMass_postcut->Draw();
  c_YMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcut.root",name.Data()), "root");
  c_YMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcut.eps",name.Data()), "eps");

  // TCanvas *c_YMass_nobeambunchcut = new TCanvas("c_YMass_nobeambunchcut","c_YMass_nobeambunchcut",600,400);
  // c_YMass_nobeambunchcut->cd();
  // TH1F *h_YMass_nobeambunchcut = (TH1F*) f->Get(Form("h_YMass_%scut_nobeambunchcut",cut.Data()));
  // cout<<"h_YMass_nobeambunchcut = "<<h_YMass_nobeambunchcut<<endl;
  // // h_YMass_nobeambunchcut->Rebin(4);
  // h_YMass_nobeambunchcut->SetMarkerSize(0.7);
  // h_YMass_nobeambunchcut->Draw();
  // c_YMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_%scut_nobeambunchcut.root",name.Data(),cut.Data()), "root");
  // c_YMass_nobeambunchcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_%scut_nobeambunchcut.eps",name.Data(),cut.Data()), "eps");

  // ******** Phi vs. fo
  TCanvas *c_h2_PhiVsfoMass_Measured = new TCanvas("c_h2_PhiVsfoMass_Measured","c_h2_PhiVsfoMass_Measured",600,400);
  c_h2_PhiVsfoMass_Measured->cd();
  TH2D *h2_PhiVsfoMass_Measured = (TH2D*) f->Get("h2_PhiVsfoMass_Measured");
  h2_PhiVsfoMass_Measured->Draw("colz");
  c_h2_PhiVsfoMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_Measured.root",name.Data()), "root");
  c_h2_PhiVsfoMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_Measured.eps",name.Data()), "eps");

  TCanvas *c_h2_PhiVsfoMass_KinFit = new TCanvas("c_h2_PhiVsfoMass_KinFit","c_h2_PhiVsfoMass_KinFit",600,400);
  c_h2_PhiVsfoMass_KinFit->cd();
  TH2D *h2_PhiVsfoMass_KinFit = (TH2D*) f->Get("h2_PhiVsfoMass_KinFit");
  h2_PhiVsfoMass_KinFit->Draw("colz");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.eps",name.Data()), "eps");

  TCanvas *c_h2_PhiVsfoMass_postcut = new TCanvas("c_h2_PhiVsfoMass_postcut","c_h2_PhiVsfoMass_postcut",600,400);
  c_h2_PhiVsfoMass_postcut->cd();
  TH2D *h2_PhiVsfoMass_postcut = (TH2D*) f->Get("h2_PhiVsfoMass_postcut");
  h2_PhiVsfoMass_postcut->Draw("colz");
  c_h2_PhiVsfoMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.root",name.Data()), "root");
  c_h2_PhiVsfoMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.eps",name.Data()), "eps");

  // ******** Phi vs. Y
  TCanvas *c_h2_PhiVsYMass_Measured = new TCanvas("c_h2_PhiVsYMass_Measured","c_h2_PhiVsYMass_Measured",600,400);
  c_h2_PhiVsYMass_Measured->cd();
  TH2D *h2_PhiVsYMass_Measured = (TH2D*) f->Get("h2_PhiVsYMass_Measured");
  h2_PhiVsYMass_Measured->Draw("colz");
  c_h2_PhiVsYMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_Measured.root",name.Data()), "root");
  c_h2_PhiVsYMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_Measured.eps",name.Data()), "eps");

  TCanvas *c_h2_PhiVsYMass_KinFit = new TCanvas("c_h2_PhiVsYMass_KinFit","c_h2_PhiVsYMass_KinFit",600,400);
  c_h2_PhiVsYMass_KinFit->cd();
  TH2D *h2_PhiVsYMass_KinFit = (TH2D*) f->Get("h2_PhiVsYMass_KinFit");
  h2_PhiVsYMass_KinFit->Draw("colz");
  c_h2_PhiVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.eps",name.Data()), "eps");

  TCanvas *c_h2_PhiVsYMass_postcut = new TCanvas("c_h2_PhiVsYMass_postcut","c_h2_PhiVsYMass_postcut",600,400);
  c_h2_PhiVsYMass_postcut->cd();
  TH2D *h2_PhiVsYMass_postcut = (TH2D*) f->Get("h2_PhiVsYMass_postcut");
  h2_PhiVsYMass_postcut->Draw("colz");
  c_h2_PhiVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.root",name.Data()), "root");
  c_h2_PhiVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.eps",name.Data()), "eps");


  // ******** fo vs. Y
  TCanvas *c_h2_foVsYMass_Measured = new TCanvas("c_h2_foVsYMass_Measured","c_h2_foVsYMass_Measured",600,400);
  c_h2_foVsYMass_Measured->cd();
  TH2D *h2_foVsYMass_Measured = (TH2D*) f->Get("h2_foVsYMass_Measured");
  h2_foVsYMass_Measured->Draw("colz");
  c_h2_foVsYMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_Measured.root",name.Data()), "root");
  c_h2_foVsYMass_Measured->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_Measured.eps",name.Data()), "eps");

  TCanvas *c_h2_foVsYMass_KinFit = new TCanvas("c_h2_foVsYMass_KinFit","c_h2_foVsYMass_KinFit",600,400);
  c_h2_foVsYMass_KinFit->cd();
  TH2D *h2_foVsYMass_KinFit = (TH2D*) f->Get("h2_foVsYMass_KinFit");
  h2_foVsYMass_KinFit->Draw("colz");
  c_h2_foVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.root",name.Data()), "root");
  c_h2_foVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.eps",name.Data()), "eps");

  TCanvas *c_h2_foVsYMass_postcut = new TCanvas("c_h2_foVsYMass_postcut","c_h2_foVsYMass_postcut",600,400);
  c_h2_foVsYMass_postcut->cd();
  TH2D *h2_foVsYMass_postcut = (TH2D*) f->Get("h2_foVsYMass_postcut");
  h2_foVsYMass_postcut->Draw("colz");
  c_h2_foVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.root",name.Data()), "root");
  c_h2_foVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.eps",name.Data()), "eps");
}