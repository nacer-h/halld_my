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

void y2175(TString name)//, TString cut)
{
  TFile *f = NULL;
  if(name == "sim") f = new TFile("/Users/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v03.root");
  if(name == "data") f = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm_17v21.root");
  TFile *fps = new TFile("/Users/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_y2175/y2175.root","UPDATE");
  TTree *t=(TTree*)f->Get("ntp");

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  //gStyle->SetLineWidth(2);
  //gStyle->SetHistLineWidth(3);
/*    
  // ******* Tagger accidentals
  TCanvas *c_TaggerAccidentals = new TCanvas("c_TaggerAccidentals","c_TaggerAccidentals",600,400);
  c_TaggerAccidentals->cd();
  TH1F *h_TaggerAccidentals = (TH1F*) f->Get("h_TaggerAccidentals");
  h_TaggerAccidentals->Draw();
  c_TaggerAccidentals->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.root",name.Data()), "root");
  c_TaggerAccidentals->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.eps",name.Data()), "eps");

  TCanvas *c_TaggerAccidentals_postcut = new TCanvas("c_TaggerAccidentals_postcut","c_TaggerAccidentals_postcut",600,400);
  c_TaggerAccidentals_postcut->cd();
  TH1F *h_TaggerAccidentals_postcut = (TH1F*) f->Get("h_TaggerAccidentals_postcut");
  h_TaggerAccidentals_postcut->Draw();
  c_TaggerAccidentals_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.root",name.Data()), "root");
  c_TaggerAccidentals_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.eps",name.Data()), "eps");

  // ********  Chi2 of KinFit
  TCanvas *c_kin_chisq = new TCanvas("c_kin_chisq","c_kin_chisq",600,400);
  c_kin_chisq->cd();
  TH1F *h_kin_chisq = new TH1F("h_kin_chisq", ";#chi^{2} of Kinematic Fit;Counts", 600, 0, 100);
  t->Project("h_kin_chisq","kin_chisq","w8*(abs(mm2)<0.06)");
  h_kin_chisq->Draw("hist");
  c_kin_chisq->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq.root",name.Data()), "root");
  c_kin_chisq->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq.eps",name.Data()), "eps");  

  // ********  Missing Mass Squared
  TCanvas *c_mm2 = new TCanvas("c_mm2","c_mm2",600,400);
  c_mm2->cd();
  TH1F *h_mm2 = new TH1F("h_mm2", ";MM^{2} (GeV/c^{2})^{2};Counts", 600, -0.06, 0.06);
  t->Project("h_mm2","mm2","w8*(kin_chisq<30)");
  h_mm2->Draw("hist");
  c_mm2->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2.root",name.Data()), "root");
  c_mm2->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2.eps",name.Data()), "eps");

  // ********  VertexZ
  // proton
  TCanvas *c_p_vertexz = new TCanvas("c_p_vertexz","c_p_vertexz",600,400);
  c_p_vertexz->cd();
  TH1F *h_p_vertexz = new TH1F("h_p_vertexz", ";Vertex-Z (cm) (Protons); Counts", 300, 45, 85);
  t->Project("h_p_vertexz","p_x4_kin.Z()","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_p_vertexz->Draw("hist");
  c_p_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexz.root",name.Data()), "root");
  c_p_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexz.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_vertexz = new TCanvas("c_pip_vertexz","c_pip_vertexz",600,400);
  c_pip_vertexz->cd();
  TH1F *h_pip_vertexz = new TH1F("h_pip_vertexz", ";Vertex-Z (cm) (#pi^{+}); Counts", 300, 45, 85);
  t->Project("h_pip_vertexz","pip_x4_kin.Z()","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pip_vertexz->Draw("hist");
  c_pip_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexz.root",name.Data()), "root");
  c_pip_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexz.eps",name.Data()), "eps");

 // pim
  TCanvas *c_pim_vertexz = new TCanvas("c_pim_vertexz","c_pim_vertexz",600,400);
  c_pim_vertexz->cd();
  TH1F *h_pim_vertexz = new TH1F("h_pim_vertexz", ";Vertex-Z (cm) (#pi^{-}); Counts", 300, 45, 85);
  t->Project("h_pim_vertexz","pim_x4_kin.Z()","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pim_vertexz->Draw("hist");
  c_pim_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexz.root",name.Data()), "root");
  c_pim_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexz.eps",name.Data()), "eps"); 

 // kp
  TCanvas *c_kp_vertexz = new TCanvas("c_kp_vertexz","c_kp_vertexz",600,400);
  c_kp_vertexz->cd();
  TH1F *h_kp_vertexz = new TH1F("h_kp_vertexz", ";Vertex-Z (cm) (K^{+}); Counts", 300, 45, 85);
  t->Project("h_kp_vertexz","kp_x4_kin.Z()","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_kp_vertexz->Draw("hist");
  c_kp_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexz.root",name.Data()), "root");
  c_kp_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexz.eps",name.Data()), "eps"); 

 // km
  TCanvas *c_km_vertexz = new TCanvas("c_km_vertexz","c_km_vertexz",600,400);
  c_km_vertexz->cd();
  TH1F *h_km_vertexz = new TH1F("h_km_vertexz", ";Vertex-Z (cm) (K^{-}); Counts", 300, 45, 85);
  t->Project("h_km_vertexz","km_x4_kin.Z()","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_km_vertexz->Draw("hist");
  c_km_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexz.root",name.Data()), "root");
  c_km_vertexz->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexz.eps",name.Data()), "eps");

   // ********  Vertex:sqrt(X^2+Y^2)
  // proton
  TCanvas *c_p_vertexxy = new TCanvas("c_p_vertexxy","c_p_vertexxy",600,400);
  c_p_vertexxy->cd();
  c_p_vertexxy->SetLogy();
  TH1F *h_p_vertexxy = new TH1F("h_p_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (Proton); Counts", 300, 0, 2);
  t->Project("h_p_vertexxy","sqrt((p_x4_kin.X()*p_x4_kin.X())+(p_x4_kin.Y()*p_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_p_vertexxy->Draw("hist");
  c_p_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexxy.root",name.Data()), "root");
  c_p_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexxy.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_vertexxy = new TCanvas("c_pip_vertexxy","c_pip_vertexxy",600,400);
  c_pip_vertexxy->cd();
  c_pip_vertexxy->SetLogy();
  TH1F *h_pip_vertexxy = new TH1F("h_pip_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (#pi^{+}); Counts", 300,  0, 2);
  t->Project("h_pip_vertexxy","sqrt((pip_x4_kin.X()*pip_x4_kin.X())+(pip_x4_kin.Y()*pip_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pip_vertexxy->Draw("hist");
  c_pip_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexxy.root",name.Data()), "root");
  c_pip_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexxy.eps",name.Data()), "eps");

 // pim
  TCanvas *c_pim_vertexxy = new TCanvas("c_pim_vertexxy","c_pim_vertexxy",600,400);
  c_pim_vertexxy->cd();
  c_pim_vertexxy->SetLogy();
  TH1F *h_pim_vertexxy = new TH1F("h_pim_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (#pi^{-}); Counts", 300,  0, 2);
  t->Project("h_pim_vertexxy","sqrt((pim_x4_kin.X()*pim_x4_kin.X())+(pim_x4_kin.Y()*pim_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pim_vertexxy->Draw("hist");
  c_pim_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexxy.root",name.Data()), "root");
  c_pim_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexxy.eps",name.Data()), "eps"); 

 // kp
  TCanvas *c_kp_vertexxy = new TCanvas("c_kp_vertexxy","c_kp_vertexxy",600,400);
  c_kp_vertexxy->cd();
  c_kp_vertexxy->SetLogy();
  TH1F *h_kp_vertexxy = new TH1F("h_kp_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (K^{+}); Counts", 300,  0, 2);
  t->Project("h_kp_vertexxy","sqrt((kp_x4_kin.X()*kp_x4_kin.X())+(kp_x4_kin.Y()*kp_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_kp_vertexxy->Draw("hist");
  c_kp_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexxy.root",name.Data()), "root");
  c_kp_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexxy.eps",name.Data()), "eps"); 

 // km
  TCanvas *c_km_vertexxy = new TCanvas("c_km_vertexxy","c_km_vertexxy",600,400);
  c_km_vertexxy->cd();
  c_km_vertexxy->SetLogy();
  TH1F *h_km_vertexxy = new TH1F("h_km_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (K^{-}); Counts", 300,  0, 2);
  t->Project("h_km_vertexxy","sqrt((km_x4_kin.X()*km_x4_kin.X())+(km_x4_kin.Y()*km_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_km_vertexxy->Draw("hist");
  c_km_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexxy.root",name.Data()), "root");
  c_km_vertexxy->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexxy.eps",name.Data()), "eps");  

  // ********  Missing Energy
  TCanvas *c_me = new TCanvas("c_me","c_me",600,400);
  c_me->cd();
  TH1F *h_me = new TH1F("h_me", ";Missing Energy (GeV);Counts", 600, -1.1, 1.1);
  t->Project("h_me","me","w8*(kin_chisq<30 && abs(mm2)<0.06)");
  h_me->Draw("hist");
  c_me->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_me.root",name.Data()), "root");
  c_me->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_me.eps",name.Data()), "eps"); 

  // ********  Total Transverse Momentum
  TCanvas *c_pt = new TCanvas("c_pt","c_pt",600,400);
  c_pt->cd();
  c_pt->SetLogy();
  TH1F *h_pt = new TH1F("h_pt", ";P_{t} total (GeV);Counts", 600, 0, 3);
  t->Project("h_pt","pt","w8*(kin_chisq<30 && abs(mm2)<0.06)");
  h_pt->Draw("hist");
  c_pt->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pt.root",name.Data()), "root");
  c_pt->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pt.eps",name.Data()), "eps"); 

  // ********  Beam Energy
  TCanvas *c_beam_e = new TCanvas("c_beam_e","c_beam_e",600,400);
  c_beam_e->cd();
  TH1F *h_beam_e = new TH1F("h_beam_e", ";Beam Energy (GeV);Counts", 100, 3.0, 11.6);
  t->Project("h_beam_e","beam_p4_kin.E()","w8*(kin_chisq<30 && abs(mm2)<0.015)");
  h_beam_e->Draw("hist");
  c_beam_e->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beam_e.root",name.Data()), "root");
  c_beam_e->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beam_e.eps",name.Data()), "eps");

  // ********  proton momentum transfer (-t)
  TCanvas *c_t_kin_precut = new TCanvas("c_t_kin_precut","c_t_kin_precut",600,400);
  c_t_kin_precut->cd();
  TH1F *h_t_kin_precut = new TH1F("h_t_kin_precut", ";-t (GeV/c^{2});Counts", 300, 0, 10);
  t->Project("h_t_kin_precut","-t_kin","w8*(kin_chisq<30 && abs(mm2)<0.015)");
  h_t_kin_precut->Draw("hist");
  c_t_kin_precut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin_precut.root",name.Data()), "root");
  c_t_kin_precut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin_precut.eps",name.Data()), "eps"); 

  // post PhiMass cut
  TCanvas *c_t_kin = new TCanvas("c_t_kin","c_t_kin",600,400);
  c_t_kin->cd();
  TH1F *h_t_kin = new TH1F("h_t_kin", ";-t (GeV/c^{2});Counts", 300, 0, 10);
  t->Project("h_t_kin","-t_kin","w8*(kin_chisq<30 && abs(mm2)<0.015 && kpkm_mf>1.005 && kpkm_mf<1.035)");
  h_t_kin->Draw("hist");
  c_t_kin->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin.root",name.Data()), "root");
  c_t_kin->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin.eps",name.Data()), "eps");

  // ********  Missing mass squared with Pions as Kions hypothesis
  TCanvas *c_mm2_piask = new TCanvas("c_mm2_piask","c_mm2_piask",600,400);
  c_mm2_piask->cd();
  TH1F *h_mm2_piask = new TH1F("h_mm2_piask", ";MM^{2} (Pions As Kaons) (GeV/c^{2})^{2};Counts", 300, -10, 4);
  t->Project("h_mm2_piask","mm2_piask","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_mm2_piask->Draw("hist");
  c_mm2_piask->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2_piask.root",name.Data()), "root");
  c_mm2_piask->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2_piask.eps",name.Data()), "eps");   

  TCanvas *c_mm2vsmm2_piask = new TCanvas("c_mm2vsmm2_piask","c_mm2vsmm2_piask",600,400);
  c_mm2vsmm2_piask->cd();
  TH2D *h_mm2vsmm2_piask = new TH2D("h_mm2vsmm2_piask", ";MM^{2} (GeV/c^{2})^{2};MM^{2} (Pions As Kaons) (GeV/c^{2})^{2}", 300, -0.06, 0.06, 300, -10, 5);
  t->Project("h_mm2vsmm2_piask","mm2_piask:mm2","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_mm2vsmm2_piask->Draw("colz");
  c_mm2vsmm2_piask->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2vsmm2_piask.root",name.Data()), "root");
  c_mm2vsmm2_piask->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2vsmm2_piask.eps",name.Data()), "eps");     
*/
/*
  // **************** P vs. theta
  // proton
  TCanvas *c_p_ptheta = new TCanvas("c_p_ptheta", "c_p_ptheta", 600, 400);
  c_p_ptheta->cd();
  TH2D *h_p_ptheta = new TH2D("h_p_ptheta", "proton ; #theta#circ ;p (GeV/c)", 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_p_ptheta", "p_p4_kin.P():p_p4_kin.Theta()*TMath::RadToDeg()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_p_ptheta = " << h_p_ptheta << endl;
  h_p_ptheta->Draw("colz");
  c_p_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_ptheta.root",name.Data()), "root");
  c_p_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_ptheta.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_ptheta = new TCanvas("c_pip_ptheta", "c_pip_ptheta", 600, 400);
  c_pip_ptheta->cd();
  TH2D *h_pip_ptheta = new TH2D("h_pip_ptheta", "#pi^{+}; #theta#circ ;p (GeV/c)", 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_pip_ptheta", "pip_p4_kin.P():pip_p4_kin.Theta()*TMath::RadToDeg()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pip_ptheta = " << h_pip_ptheta << endl;
  h_pip_ptheta->Draw("colz");
  c_pip_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_ptheta.root",name.Data()), "root");
  c_pip_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_ptheta.eps",name.Data()), "eps");

  // pim
  TCanvas *c_pim_ptheta = new TCanvas("c_pim_ptheta", "c_pim_ptheta", 600, 400);
  c_pim_ptheta->cd();
  TH2D *h_pim_ptheta = new TH2D("h_pim_ptheta", "#pi^{-}; #theta#circ ;p (GeV/c)", 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_pim_ptheta", "pim_p4_kin.P():pim_p4_kin.Theta()*TMath::RadToDeg()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pim_ptheta = " << h_pim_ptheta << endl;
  h_pim_ptheta->Draw("colz");
  c_pim_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_ptheta.root",name.Data()), "root");
  c_pim_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_ptheta.eps",name.Data()), "eps");

  // kp
  TCanvas *c_kp_ptheta = new TCanvas("c_kp_ptheta", "c_kp_ptheta", 600, 400);
  c_kp_ptheta->cd();
  TH2D *h_kp_ptheta = new TH2D("h_kp_ptheta", "K^{+}; #theta#circ ;p (GeV/c)", 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_kp_ptheta", "kp_p4_kin.P():kp_p4_kin.Theta()*TMath::RadToDeg()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_kp_ptheta = " << h_kp_ptheta << endl;
  h_kp_ptheta->Draw("colz");
  c_kp_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_ptheta.root",name.Data()), "root");
  c_kp_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_ptheta.eps",name.Data()), "eps");
  
  // km
  TCanvas *c_km_ptheta = new TCanvas("c_km_ptheta", "c_km_ptheta", 600, 400);
  c_km_ptheta->cd();
  TH2D *h_km_ptheta = new TH2D("h_km_ptheta", "K^{-}; #theta#circ ;p (GeV/c)", 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_km_ptheta", "km_p4_kin.P():km_p4_kin.Theta()*TMath::RadToDeg()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_km_ptheta = " << h_km_ptheta << endl;
  h_km_ptheta->Draw("colz");
  c_km_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_ptheta.root",name.Data()), "root");
  c_km_ptheta->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_ptheta.eps",name.Data()), "eps");
*/
/*
  // **************** dt vs. p
  // proton
  // TOF
  TCanvas *c_p_dttof = new TCanvas("c_p_dttof", "c_p_dttof", 600, 400);
  c_p_dttof->cd();
  TH2D *h_p_dttof = new TH2D("h_p_dttof", "proton; p (GeV/c); TOF #Delta T (ns)", 250, 0.0, 10.0, 500,-0.6, 0.6);
  t->Project("h_p_dttof", "p_dttof:p_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_p_dttof = " << h_p_dttof << endl;
  h_p_dttof->Draw("colz");
  c_p_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dttof.root",name.Data()), "root");
  c_p_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_p_dtbcal = new TCanvas("c_p_dtbcal", "c_p_dtbcal", 600, 400);
  c_p_dtbcal->cd();
  TH2D *h_p_dtbcal = new TH2D("h_p_dtbcal", "proton; p (GeV/c); BCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_p_dtbcal", "p_dtbcal:p_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_p_dtbcal = " << h_p_dtbcal << endl;
  h_p_dtbcal->Draw("colz");
  c_p_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtbcal.root",name.Data()), "root");
  c_p_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_p_dtfcal = new TCanvas("c_p_dtfcal", "c_p_dtfcal", 600, 400);
  c_p_dtfcal->cd();
  TH2D *h_p_dtfcal = new TH2D("h_p_dtfcal", "proton; p (GeV/c); FCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_p_dtfcal", "p_dtfcal:p_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_p_dtfcal = " << h_p_dtfcal << endl;
  h_p_dtfcal->Draw("colz");
  c_p_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtfcal.root",name.Data()), "root");
  c_p_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtfcal.eps",name.Data()), "eps");

  // pip
  // TOF
  TCanvas *c_pip_dttof = new TCanvas("c_pip_dttof", "c_pip_dttof", 600, 400);
  c_pip_dttof->cd();
  TH2D *h_pip_dttof = new TH2D("h_pip_dttof", "#pi^{+}; p (GeV/c); TOF #Delta T (ns)", 250, 0.0, 10.0, 500,-0.5, 0.5);
  t->Project("h_pip_dttof", "pip_dttof:pip_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pip_dttof = " << h_pip_dttof << endl;
  h_pip_dttof->Draw("colz");
  c_pip_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dttof.root",name.Data()), "root");
  c_pip_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_pip_dtbcal = new TCanvas("c_pip_dtbcal", "c_pip_dtbcal", 600, 400);
  c_pip_dtbcal->cd();
  TH2D *h_pip_dtbcal = new TH2D("h_pip_dtbcal", "#pi^{+}; p (GeV/c); BCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_pip_dtbcal", "pip_dtbcal:pip_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pip_dtbcal = " << h_pip_dtbcal << endl;
  h_pip_dtbcal->Draw("colz");
  c_pip_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtbcal.root",name.Data()), "root");
  c_pip_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_pip_dtfcal = new TCanvas("c_pip_dtfcal", "c_pip_dtfcal", 600, 400);
  c_pip_dtfcal->cd();
  TH2D *h_pip_dtfcal = new TH2D("h_pip_dtfcal", "#pi^{+}; p (GeV/c); FCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_pip_dtfcal", "pip_dtfcal:pip_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pip_dtfcal = " << h_pip_dtfcal << endl;
  h_pip_dtfcal->Draw("colz");
  c_pip_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtfcal.root",name.Data()), "root");
  c_pip_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtfcal.eps",name.Data()), "eps");

  // pim
  // TOF
  TCanvas *c_pim_dttof = new TCanvas("c_pim_dttof", "c_pim_dttof", 600, 400);
  c_pim_dttof->cd();
  TH2D *h_pim_dttof = new TH2D("h_pim_dttof", "#pi^{-}; p (GeV/c); TOF #Delta T (ns)", 250, 0.0, 10.0, 500,-0.5, 0.5);
  t->Project("h_pim_dttof", "pim_dttof:pim_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pim_dttof = " << h_pim_dttof << endl;
  h_pim_dttof->Draw("colz");
  c_pim_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dttof.root",name.Data()), "root");
  c_pim_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_pim_dtbcal = new TCanvas("c_pim_dtbcal", "c_pim_dtbcal", 600, 400);
  c_pim_dtbcal->cd();
  TH2D *h_pim_dtbcal = new TH2D("h_pim_dtbcal", "#pi^{-}; p (GeV/c); BCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_pim_dtbcal", "pim_dtbcal:pim_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pim_dtbcal = " << h_pim_dtbcal << endl;
  h_pim_dtbcal->Draw("colz");
  c_pim_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtbcal.root",name.Data()), "root");
  c_pim_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_pim_dtfcal = new TCanvas("c_pim_dtfcal", "c_pim_dtfcal", 600, 400);
  c_pim_dtfcal->cd();
  TH2D *h_pim_dtfcal = new TH2D("h_pim_dtfcal", "#pi^{-}; p (GeV/c); FCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_pim_dtfcal", "pim_dtfcal:pim_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_pim_dtfcal = " << h_pim_dtfcal << endl;
  h_pim_dtfcal->Draw("colz");
  c_pim_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtfcal.root",name.Data()), "root");
  c_pim_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtfcal.eps",name.Data()), "eps");

  // kp
  // TOF
  TCanvas *c_kp_dttof = new TCanvas("c_kp_dttof", "c_kp_dttof", 600, 400);
  c_kp_dttof->cd();
  TH2D *h_kp_dttof = new TH2D("h_kp_dttof", "K^{+}; p (GeV/c); TOF #Delta T (ns)", 250, 0.0, 10.0, 500,-0.3, 0.3);
  t->Project("h_kp_dttof", "kp_dttof:kp_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_kp_dttof = " << h_kp_dttof << endl;
  h_kp_dttof->Draw("colz");
  c_kp_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dttof.root",name.Data()), "root");
  c_kp_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_kp_dtbcal = new TCanvas("c_kp_dtbcal", "c_kp_dtbcal", 600, 400);
  c_kp_dtbcal->cd();
  TH2D *h_kp_dtbcal = new TH2D("h_kp_dtbcal", "K^{+}; p (GeV/c); BCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-0.75, 0.75);
  t->Project("h_kp_dtbcal", "kp_dtbcal:kp_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_kp_dtbcal = " << h_kp_dtbcal << endl;
  h_kp_dtbcal->Draw("colz");
  c_kp_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtbcal.root",name.Data()), "root");
  c_kp_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_kp_dtfcal = new TCanvas("c_kp_dtfcal", "c_kp_dtfcal", 600, 400);
  c_kp_dtfcal->cd();
  TH2D *h_kp_dtfcal = new TH2D("h_kp_dtfcal", "K^{+}; p (GeV/c); FCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-2.5, 2.5);
  t->Project("h_kp_dtfcal", "kp_dtfcal:kp_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_kp_dtfcal = " << h_kp_dtfcal << endl;
  h_kp_dtfcal->Draw("colz");
  c_kp_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtfcal.root",name.Data()), "root");
  c_kp_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtfcal.eps",name.Data()), "eps");

  // km
  // TOF
  TCanvas *c_km_dttof = new TCanvas("c_km_dttof", "c_km_dttof", 600, 400);
  c_km_dttof->cd();
  TH2D *h_km_dttof = new TH2D("h_km_dttof", "K^{-}; p (GeV/c); TOF #Delta T (ns)", 250, 0.0, 10.0, 500,-0.3, 0.3);
  t->Project("h_km_dttof", "km_dttof:km_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_km_dttof = " << h_km_dttof << endl;
  h_km_dttof->Draw("colz");
  c_km_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dttof.root",name.Data()), "root");
  c_km_dttof->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_km_dtbcal = new TCanvas("c_km_dtbcal", "c_km_dtbcal", 600, 400);
  c_km_dtbcal->cd();
  TH2D *h_km_dtbcal = new TH2D("h_km_dtbcal", "K^{-}; p (GeV/c); BCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-0.75, 0.75);
  t->Project("h_km_dtbcal", "km_dtbcal:km_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_km_dtbcal = " << h_km_dtbcal << endl;
  h_km_dtbcal->Draw("colz");
  c_km_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtbcal.root",name.Data()), "root");
  c_km_dtbcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_km_dtfcal = new TCanvas("c_km_dtfcal", "c_km_dtfcal", 600, 400);
  c_km_dtfcal->cd();
  TH2D *h_km_dtfcal = new TH2D("h_km_dtfcal", "K^{-}; p (GeV/c); FCAL #Delta T (ns)", 250, 0.0, 10.0, 500,-2.5, 2.5);
  t->Project("h_km_dtfcal", "km_dtfcal:km_p4_kin.P()", "kin_chisq<30 && abs(mm2)<0.06");
  cout << "h_km_dtfcal = " << h_km_dtfcal << endl;
  h_km_dtfcal->Draw("colz");
  c_km_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtfcal.root",name.Data()), "root");
  c_km_dtfcal->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtfcal.eps",name.Data()), "eps");  
*/
/*
  // ********************************* Phi(1020) Mass: Measured & KinFit *************************************
  TCanvas *c_PhiMass = new TCanvas("c_PhiMass","c_PhiMass",600,400);
  c_PhiMass->cd();
  TLegend *l_PhiMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_PhiMass->SetTextSize(0.02);                                                                                                                                                 
  l_PhiMass->SetBorderSize(0); 
  c_PhiMass->SetGrid();
  TH1F *h_PhiMass_Measured = new TH1F("h_PhiMass_Measured", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 0.98, 1.2);
  t->Project("h_PhiMass_Measured", "kpkm_m", "kpkm_uni && abs(rf_dt)<2.004");
  // (TH1F*) f->Get("PhiMass_Measured");
  cout<<"h_PhiMass_Measured = "<<h_PhiMass_Measured<<endl;
  h_PhiMass_Measured->SetMinimum(0.);
  h_PhiMass_Measured->SetLineColor(kBlack);
  TH1F *h_PhiMass_KinFit = new TH1F("h_PhiMass_KinFit", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 0.98, 1.2);
  t->Project("h_PhiMass_KinFit", "kpkm_mf", "kpkm_uni && abs(rf_dt)<2.004 && abs(mm2)<0.015 && kin_chisq<30");
  h_PhiMass_KinFit->SetLineColor(kBlue);
  TH1F *h_PhiMass_beambunchcut = new TH1F("h_PhiMass_beambunchcut", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 0.98, 1.2);
  t->Project("h_PhiMass_beambunchcut", "kpkm_mf", "w8*(kpkm_uni && abs(mm2)<0.06 && kin_chisq<100)");
  h_PhiMass_beambunchcut->SetLineColor(kRed);
  
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 0.98, 1.2);
  t->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni && abs(mm2)<0.015 && kin_chisq<30)");
  h_PhiMass_postcuts->SetLineColor(kMagenta);
  h_PhiMass_KinFit->Draw();
  h_PhiMass_Measured->Draw("same");
  h_PhiMass_beambunchcut->Draw("same hist");
  h_PhiMass_postcuts->Draw("same hist");
  l_PhiMass->AddEntry(h_PhiMass_Measured, "Measured + |#Delta T|<2ns","l"); 
  l_PhiMass->AddEntry(h_PhiMass_KinFit, "KinFit + |#Delta T|<2ns","l");
  l_PhiMass->AddEntry(h_PhiMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_PhiMass->AddEntry(h_PhiMass_postcuts, "1.005<m_{K^{+}K^{-}}<1.035","l");
  // l_PhiMass->Draw();
  c_PhiMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.root",name.Data()), "root");
  c_PhiMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.eps",name.Data()), "eps");
*/
  // +++ Post cut
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.98, 1.2);
  t->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");  
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts","c_PhiMass_postcuts",600,400);
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");
  c_PhiMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcuts.root",name.Data()), "root");
  c_PhiMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcuts.eps",name.Data()), "eps");
/*
  // ******* fo(980) Mass: Measured & KinFit ********************
  TCanvas *c_foMass = new TCanvas("c_foMass","c_foMass",600,400);
  c_foMass->cd();
  TLegend *l_foMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_foMass->SetTextSize(0.02);                                                                                                                                                 
  l_foMass->SetBorderSize(0); 
  c_foMass->SetGrid();
  TH1F *h_foMass_Measured = new TH1F("h_foMass_Measured", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.3, 1.3);
  t->Project("h_foMass_Measured", "pippim_m", "pippim_uni && abs(rf_dt)<2.004 && abs(mm2)<0.06 && kin_chisq<100");
  h_foMass_Measured->SetMinimum(0.);
  // h_foMass_Measured->SetLineWidth(3);
  h_foMass_Measured->SetLineColor(kBlack);
  TH1F *h_foMass_KinFit = new TH1F("h_foMass_KinFit", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.3, 1.3);
  t->Project("h_foMass_KinFit", "pippim_mf", "pippim_uni && abs(rf_dt)<2.004 && abs(mm2)<0.06 && kin_chisq<100");
  h_foMass_KinFit->SetLineColor(kBlue);
  TH1F *h_foMass_beambunchcut = new TH1F("h_foMass_beambunchcut", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.3, 1.3);
  t->Project("h_foMass_beambunchcut", "pippim_mf", "w8*(pippim_uni && abs(mm2)<0.06 && kin_chisq<100)");
  h_foMass_beambunchcut->SetLineColor(kRed);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.3, 1.3);
  t->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni && abs(mm2)<0.06 && kin_chisq<100 && kpkm_mf>1.005 && kpkm_mf<1.035)");
  h_foMass_postcuts->SetLineColor(kMagenta);
  h_foMass_KinFit->Draw();
  h_foMass_Measured->Draw("same");
  h_foMass_beambunchcut->Draw("hist same");
  h_foMass_postcuts->Draw("hist same");
  l_foMass->AddEntry(h_foMass_Measured, "Measured + |#Delta T|<2ns","l"); 
  l_foMass->AddEntry(h_foMass_KinFit, "KinFit + |#Delta T|<2ns","l");
  l_foMass->AddEntry(h_foMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_foMass->AddEntry(h_foMass_postcuts, "1.005<m_{K^{+}K^{-}}<1.035","l");
  l_foMass->Draw();
  c_foMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.root",name.Data()), "root");
  c_foMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.eps",name.Data()), "eps");
*/
  // +++ Post cut
  TCanvas *c_foMass_postcuts = new TCanvas("c_foMass_postcuts","c_foMass_postcuts",600,400);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 0.3, 1.3);
  t->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni)");
  c_foMass_postcuts->cd();
  h_foMass_postcuts->Draw("e");
  c_foMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcuts.root",name.Data()), "root");
  c_foMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcuts.eps",name.Data()), "eps");  
/*
  // ******* Y(2175) Mass: Measured & KinFit
  TCanvas *c_YMass = new TCanvas("c_YMass","c_YMass",600,400);
  c_YMass->cd();
  TLegend *l_YMass = new TLegend(0.7,0.7,0.89,0.89);                                                                                                                          
  l_YMass->SetTextSize(0.02);                                                                                                                                                 
  l_YMass->SetBorderSize(0); 
  c_YMass->SetGrid();
  TH1F *h_YMass_Measured = new TH1F("h_YMass_Measured", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.6, 3.2);
  t->Project("h_YMass_Measured", "kpkmpippim_m", "kpkmpippim_uni && abs(rf_dt)<2.004 && abs(mm2)<0.06 && kin_chisq<100");
  h_YMass_Measured->SetMinimum(0.);
  h_YMass_Measured->SetLineColor(kBlack);
  TH1F *h_YMass_KinFit = new TH1F("h_YMass_KinFit", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.6, 3.2);
  t->Project("h_YMass_KinFit", "kpkmpippim_mf", "kpkmpippim_uni && abs(rf_dt)<2.004 && abs(mm2)<0.06 && kin_chisq<100");
  h_YMass_KinFit->SetLineColor(kBlue);
  TH1F *h_YMass_beambunchcut = new TH1F("h_YMass_beambunchcut", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.6, 3.2);
  t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && abs(mm2)<0.06 && kin_chisq<100)");
  h_YMass_beambunchcut->SetLineColor(kRed);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_beambunchcut", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.6, 3.2);
  t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && abs(mm2)<0.06 && kin_chisq<100 && kpkm_mf>1.005 && kpkm_mf<1.035)");
  h_YMass_postcuts->SetLineColor(kMagenta);
  h_YMass_KinFit->Draw();
  h_YMass_Measured->Draw("same");
  h_YMass_beambunchcut->Draw("same hist");
  h_YMass_postcuts->Draw("same hist");
  l_YMass->AddEntry(h_YMass_Measured, "Measured + |#Delta T|<2ns","l"); 
  l_YMass->AddEntry(h_YMass_KinFit, "KinFit + |#Delta T|<2ns","l");
  l_YMass->AddEntry(h_YMass_beambunchcut, "KinFit + accidental subtracted","l");
  l_YMass->AddEntry(h_YMass_postcuts, "1.005<m_{K^{+}K^{-}}<1.035","l");
  l_YMass->Draw();
  c_YMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.root",name.Data()), "root");
  c_YMass->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.eps",name.Data()), "eps");
*/
  // +++ Post cut
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts","c_YMass_postcuts",600,400);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_beambunchcut", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni)");
  // TH1F *h_YMass_postcuts = new TH1F("h_YMass_beambunchcut", ";m_{#phif_{0}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pippim_mf-0.99)<0.1)"); 
  c_YMass_postcuts->cd();
  h_YMass_postcuts->Draw("e");//"hist"
  c_YMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcuts.root",name.Data()), "root");
  c_YMass_postcuts->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcuts.eps",name.Data()), "eps");

/*
  // ******** Phi vs. fo

  TCanvas *c_h2_PhiVsfoMass_KinFit = new TCanvas("c_h2_PhiVsfoMass_KinFit","c_h2_PhiVsfoMass_KinFit",600,400);
  c_h2_PhiVsfoMass_KinFit->cd();
  TH2D *h2_PhiVsfoMass_KinFit = (TH2D*) f->Get("h2_PhiVsfoMass_KinFit");
  h2_PhiVsfoMass_KinFit->Draw("colz");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_PhiVsfoMass_postcut = new TCanvas("c_h2_PhiVsfoMass_postcut","c_h2_PhiVsfoMass_postcut",600,400);
  // c_h2_PhiVsfoMass_postcut->cd();
  // TH2D *h2_PhiVsfoMass_postcut = (TH2D*) f->Get("h2_PhiVsfoMass_postcut");
  // h2_PhiVsfoMass_postcut->Draw("colz");
  // c_h2_PhiVsfoMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.root",name.Data()), "root");
  // c_h2_PhiVsfoMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.eps",name.Data()), "eps");

  // ******** Phi vs. Y

  TCanvas *c_h2_PhiVsYMass_KinFit = new TCanvas("c_h2_PhiVsYMass_KinFit","c_h2_PhiVsYMass_KinFit",600,400);
  c_h2_PhiVsYMass_KinFit->cd();
  TH2D *h2_PhiVsYMass_KinFit = (TH2D*) f->Get("h2_PhiVsYMass_KinFit");
  h2_PhiVsYMass_KinFit->Draw("colz");
  c_h2_PhiVsYMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsYMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_PhiVsYMass_postcut = new TCanvas("c_h2_PhiVsYMass_postcut","c_h2_PhiVsYMass_postcut",600,400);
  // c_h2_PhiVsYMass_postcut->cd();
  // TH2D *h2_PhiVsYMass_postcut = (TH2D*) f->Get("h2_PhiVsYMass_postcut");
  // h2_PhiVsYMass_postcut->Draw("colz");
  // c_h2_PhiVsYMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.root",name.Data()), "root");
  // c_h2_PhiVsYMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.eps",name.Data()), "eps");


  // ******** fo vs. Y

  TCanvas *c_h2_foVsYMass_KinFit = new TCanvas("c_h2_foVsYMass_KinFit","c_h2_foVsYMass_KinFit",600,400);
  c_h2_foVsYMass_KinFit->cd();
  TH2D *h2_foVsYMass_KinFit = (TH2D*) f->Get("h2_foVsYMass_KinFit");
  h2_foVsYMass_KinFit->Draw("colz");
  c_h2_foVsYMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.root",name.Data()), "root");
  c_h2_foVsYMass_KinFit->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_foVsYMass_postcut = new TCanvas("c_h2_foVsYMass_postcut","c_h2_foVsYMass_postcut",600,400);
  // c_h2_foVsYMass_postcut->cd();
  // TH2D *h2_foVsYMass_postcut = (TH2D*) f->Get("h2_foVsYMass_postcut");
  // h2_foVsYMass_postcut->Draw("colz");
  // c_h2_foVsYMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.root",name.Data()), "root");
  // c_h2_foVsYMass_postcut->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.eps",name.Data()), "eps");

*/ 
/*
  // +++ kppim
  TH1F *h_kppim = new TH1F("h_kppim", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kppim", "kppim_mf", "w8*(kppim_uni)");  
  TCanvas *c_kppim = new TCanvas("c_kppim","c_kppim",600,400);
  c_kppim->cd();
  h_kppim->Draw("hist");
  c_kppim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppim.root",name.Data()), "root");
  c_kppim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppim.eps",name.Data()), "eps");
 
  // +++ kmpip
  TH1F *h_kmpip = new TH1F("h_kmpip", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kmpip", "kmpip_mf", "w8*(kmpip_uni)");  
  TCanvas *c_kmpip = new TCanvas("c_kmpip","c_kmpip",600,400);
  c_kmpip->cd();
  h_kmpip->Draw("hist");
  c_kmpip->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpip.root",name.Data()), "root");
  c_kmpip->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpip.eps",name.Data()), "eps");
 
  // +++ pimp
  TH1F *h_pimp = new TH1F("h_pimp", ";m_{#pi^{-}p} [GeV/c^{2}];Counts", 200, 0.9, 2.5);
  t->Project("h_pimp", "pimp_mf", "w8*(pimp_uni)");  
  TCanvas *c_pimp = new TCanvas("c_pimp","c_pimp",600,400);
  c_pimp->cd();
  h_pimp->Draw("hist");
  c_pimp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pimp.root",name.Data()), "root");
  c_pimp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pimp.eps",name.Data()), "eps");
 
  // +++ kmp
  TH1F *h_kmp = new TH1F("h_kmp", ";m_{K^{-}p} [GeV/c^{2}];Counts", 200, 1.4, 3);
  t->Project("h_kmp", "kmp_mf", "w8*(kmp_uni)");  
  TCanvas *c_kmp = new TCanvas("c_kmp","c_kmp",600,400);
  c_kmp->cd();
  h_kmp->Draw("hist");
  c_kmp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmp.root",name.Data()), "root");
  c_kmp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmp.eps",name.Data()), "eps");

  // +++ kpkmp
  TH1F *h_kpkmp = new TH1F("h_kpkmp", ";m_{K^{+}K^{-}p} [GeV/c^{2}];Counts", 200, 2, 4);
  t->Project("h_kpkmp", "kpkmp_mf", "w8*(kpkmp_uni)");  
  TCanvas *c_kpkmp = new TCanvas("c_kpkmp","c_kpkmp",600,400);
  c_kpkmp->cd();
  h_kpkmp->Draw("hist");
  c_kpkmp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmp.root",name.Data()), "root");
  c_kpkmp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmp.eps",name.Data()), "eps");

  // +++ kpkmpip
  TH1F *h_kpkmpip = new TH1F("h_kpkmpip", ";m_{K^{+}K^{-}#pi^{+}} [GeV/c^{2}];Counts", 200, 1.2, 3.2);
  t->Project("h_kpkmpip", "kpkmpip_mf", "w8*(kpkmpip_uni)");  
  TCanvas *c_kpkmpip = new TCanvas("c_kpkmpip","c_kpkmpip",600,400);
  c_kpkmpip->cd();
  h_kpkmpip->Draw("hist");
  c_kpkmpip->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpip.root",name.Data()), "root");
  c_kpkmpip->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpip.eps",name.Data()), "eps");

  // +++ kpkmpim
  TH1F *h_kpkmpim = new TH1F("h_kpkmpim", ";m_{K^{+}K^{-}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.2, 3.2);
  t->Project("h_kpkmpim", "kpkmpim_mf", "w8*(kpkmpim_uni)");  
  TCanvas *c_kpkmpim = new TCanvas("c_kpkmpim","c_kpkmpim",600,400);
  c_kpkmpim->cd();
  h_kpkmpim->Draw("hist");
  c_kpkmpim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpim.root",name.Data()), "root");
  c_kpkmpim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpim.eps",name.Data()), "eps");

  // +++ ppippim
  TH1F *h_ppippim = new TH1F("h_ppippim", ";m_{p#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.3, 3.2);
  t->Project("h_ppippim", "ppippim_mf", "w8*(ppippim_uni)");  
  TCanvas *c_ppippim = new TCanvas("c_ppippim","c_ppippim",600,400);
  c_ppippim->cd();
  h_ppippim->Draw("hist");
  c_ppippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_ppippim.root",name.Data()), "root");
  c_ppippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_ppippim.eps",name.Data()), "eps");

  // +++ kppippim
  TH1F *h_kppippim = new TH1F("h_kppippim", ";m_{K^{+}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.9, 2.7);
  t->Project("h_kppippim", "kppippim_mf", "w8*(kppippim_uni)");  
  TCanvas *c_kppippim = new TCanvas("c_kppippim","c_kppippim",600,400);
  c_kppippim->cd();
  h_kppippim->Draw("hist");
  c_kppippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppippim.root",name.Data()), "root");
  c_kppippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppippim.eps",name.Data()), "eps");

  // +++ kmpippim
  TH1F *h_kmpippim = new TH1F("h_kmpippim", ";m_{K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.9, 2.7);
  t->Project("h_kmpippim", "kmpippim_mf", "w8*(kmpippim_uni)");  
  TCanvas *c_kmpippim = new TCanvas("c_kmpippim","c_kmpippim",600,400);
  c_kmpippim->cd();
  h_kmpippim->Draw("hist");
  c_kmpippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpippim.root",name.Data()), "root");
  c_kmpippim->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpippim.eps",name.Data()), "eps");

  // +++ kppimp
  TH1F *h_kppimp = new TH1F("h_kppimp", ";m_{K^{+}#pi^{-}p} [GeV/c^{2}];Counts", 200, 1.6, 3.5);
  t->Project("h_kppimp", "kppimp_mf", "w8*(kppimp_uni)");  
  TCanvas *c_kppimp = new TCanvas("c_kppimp","c_kppimp",600,400);
  c_kppimp->cd();
  h_kppimp->Draw("hist");
  c_kppimp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimp.root",name.Data()), "root");
  c_kppimp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimp.eps",name.Data()), "eps");

  // +++ kmpipp
  TH1F *h_kmpipp = new TH1F("h_kmpipp", ";m_{K^{-}#pi^{+}p} [GeV/c^{2}];Counts", 200, 1.6, 3.5);
  t->Project("h_kmpipp", "kmpipp_mf", "w8*(kmpipp_uni)");  
  TCanvas *c_kmpipp = new TCanvas("c_kmpipp","c_kmpipp",600,400);
  c_kmpipp->cd();
  h_kmpipp->Draw("hist");
  c_kmpipp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpipp.root",name.Data()), "root");
  c_kmpipp->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpipp.eps",name.Data()), "eps");
*/
                                          
}