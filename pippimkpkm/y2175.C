// File: y2175.C # Author: Nacer # Date: 2.9.2018 # Email: a.hamdi@gsi.de # Description: Macro to study Y(2175).                                                                                                                                                                                  
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
#include "TBox.h"
#include "TList.h"                                                                                                                            
using namespace std;

void y2175(TString name)//, TString cut)
{
  TFile *f = NULL;
  if(name == "mc_phifo") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_phifo_genr8_17v3_flat.root");
  if(name == "mc_phi2pi") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_phi2pi_genr8_17v3_flat.root");
  if(name == "mc_kkpipi") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_kkpipi_genr8_17v3_flat.root");
  if(name == "mc_bggen") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_bggen_17v3_flat.root");
  if(name == "mc_all") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_sim_17v3_flat.root");
  if(name == "data_16") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_16_flat.root");
  if(name == "data_17") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_17_flat.root");
  if(name == "data_18") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_18_flat.root");
  if(name == "data_all") f = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_all_flat.root");
  // TFile *fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/y2175.root","UPDATE");
  TTree *t=(TTree*)f->Get("ntp");

  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerSize(0.5);
  gStyle->SetMarkerStyle(8);
  // gStyle->SetMarkerSize(1);
  // gStyle->SetLineWidth(1);
  // gStyle->SetHistLineWidth(3);
  // TString cutlist = "(abs(kp_dttof)<0.135 || kp_dttof == -999) && (abs(kp_dtbcal)<0.4125 || kp_dtbcal == -999) && (abs(kp_dtfcal)<1 || kp_dtfcal == -999) && (abs(pip_dttof)<0.15 || pip_dttof == -999) && (abs(pip_dtbcal)<0.85 || pip_dtbcal == -999) && (abs(p_dttof)<0.27 || p_dttof == -999) && abs(mm2)<0.01 && kp_x4_kin.Z()<80 && kp_x4_kin.Z()>50";// abs(mm2)<0.01 && kp_x4_kin.Z()<80 && kp_x4_kin.Z()>50 &&

  // // Van Hove
  // TCanvas *c_vhthetaphi = new TCanvas("c_vhthetaphi", "c_vhthetaphi", 900,600);
  // c_vhthetaphi->cd();
  // TH2D *h2_vhthetaphi = new TH2D("h2_vhthetaphi", Form("Great Circle Plot (%s); #phi (rad) ;#theta (rad)",name.Data()), 260, 0.0, 6.3, 130, 0.0, 3.25);
  // t->Project("h2_vhthetaphi", "vhtheta:vhphi", "vhtheta>0.4 && vhtheta<1.5 && vhphi>1.55 && vhphi<2.2");//vhtheta>0.4 && vhtheta<1.5 && vhphi>1 && vhphi<2.2   //vhphi>1.55
  // cout << "h2_vhthetaphi = " << h2_vhthetaphi << endl;
  // h2_vhthetaphi->Draw("colz");
  // TF1 *gc_one = new TF1("gc_one", "TMath::ATan2(1,((2*TMath::Sqrt(2))*TMath::Sin(x)))", 0, 6.3);
  // gc_one->Draw("same");
  // TF1 *gc_two = new TF1("gc_two", "TMath::ATan2(1,((2*TMath::Sqrt(2))*TMath::Sin(x-((2*TMath::Pi())/3))))", 0, 6.3);
  // gc_two->Draw("same");
  // TF1 *gc_three = new TF1("gc_three", "TMath::ATan2(1,((2*TMath::Sqrt(2))*TMath::Sin(x-((4*TMath::Pi())/3))))", 0, 6.3);
  // gc_three->Draw("same");
  // TLine *gc_four = new TLine(0, 1.5708, 6.3, 1.5708);
  // gc_four->SetLineColor(2);
  // gc_four->SetLineWidth(2);
  // gc_four->Draw("same");
  // h2_vhthetaphi->Write(Form("h%s_vhthetaphi",name.Data()),TObject::kWriteDelete);
  // c_vhthetaphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_vhthetaphi.root",name.Data()), "root");
  // c_vhthetaphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_vhthetaphi.eps",name.Data()), "eps");

/*  
if(name == "data")
{
   // ********  Chi2pi of KinFit
  TCanvas *c_kin_chisq_chi2pi = new TCanvas("c_kin_chisq_chi2pi","c_kin_chisq_chi2pi",900,600);
  c_kin_chisq_chi2pi->cd();
  TH2D *h2_kin_chisq_chi2pi = new TH2D("h2_kin_chisq_chi2pi", "#chi^{2} of KinFit with different hypothesis;#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}};#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}", 100, 0, 2500, 100, 0, 25);
  t->Project("h2_kin_chisq_chi2pi","kin_chisq:chi2pi","w8");//("+cutlist+"&& abs(mm2)<0.01 && kp_x4_kin.Z()<80 && kp_x4_kin.Z()>50)
  h2_kin_chisq_chi2pi->Draw("colz");
  h2_kin_chisq_chi2pi->Write();
  TLine *l_kin_chisq_05chi2pi = new TLine(0, 0, 50, 25);
  l_kin_chisq_05chi2pi->SetLineColor(kRed);
  l_kin_chisq_05chi2pi->SetLineWidth(2);  
  l_kin_chisq_05chi2pi->Draw("same"); 
  TLine *l_kin_chisq_01chi2pi = new TLine(0, 0, 250, 25);
  l_kin_chisq_01chi2pi->SetLineColor(kMagenta);
  l_kin_chisq_01chi2pi->SetLineWidth(2);  
  l_kin_chisq_01chi2pi->Draw("same");  
  TLine *l_kin_chisq_001chi2pi = new TLine(0, 0, 2500, 25);
  l_kin_chisq_001chi2pi->SetLineColor(kGreen);
  l_kin_chisq_001chi2pi->SetLineWidth(2);  
  l_kin_chisq_001chi2pi->Draw("same"); 
  h2_kin_chisq_chi2pi->Write();
  c_kin_chisq_chi2pi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq_chi2pi.root",name.Data()), "root");
  c_kin_chisq_chi2pi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq_chi2pi.eps",name.Data()), "eps");  
}
  // ******* Tagger accidentals
  TCanvas *c_TaggerAccidentals = new TCanvas("c_TaggerAccidentals","c_TaggerAccidentals",900,600);
  c_TaggerAccidentals->cd();
  TH1F *h_TaggerAccidentals = new TH1F("h_TaggerAccidentals", Form("Vertex time - RF (%s);#Deltat_{Beam-RF}(ns)",name.Data()), 300, -18, 18);
  t->Project("h_TaggerAccidentals","rf_dt","");
  TLine *l_TaggerAccidentals1 = new TLine(2.004, 0, 2.004, h_TaggerAccidentals->GetMaximum());
  TLine *l_TaggerAccidentals2 = new TLine(-2.004, 0, -2.004, h_TaggerAccidentals->GetMaximum());
  TLine *l_TaggerAccidentals3 = new TLine(6.012, 0, 6.012, h_TaggerAccidentals->GetMaximum());
  TLine *l_TaggerAccidentals4 = new TLine(-6.012, 0, -6.012, h_TaggerAccidentals->GetMaximum());
  l_TaggerAccidentals1->SetLineColor(kRed);
  l_TaggerAccidentals2->SetLineColor(kRed);
  l_TaggerAccidentals3->SetLineColor(kBlue);
  l_TaggerAccidentals4->SetLineColor(kBlue);
  l_TaggerAccidentals1->SetLineStyle(10);
  l_TaggerAccidentals2->SetLineStyle(10);
  l_TaggerAccidentals3->SetLineStyle(10);
  l_TaggerAccidentals4->SetLineStyle(10);
  h_TaggerAccidentals->Draw();
  l_TaggerAccidentals1->Draw("same");
  l_TaggerAccidentals2->Draw("same");
  l_TaggerAccidentals3->Draw("same");
  l_TaggerAccidentals4->Draw("same");
  h_TaggerAccidentals->Write();
  c_TaggerAccidentals->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.root",name.Data()), "root");
  c_TaggerAccidentals->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals.eps",name.Data()), "eps");

  TCanvas *c_TaggerAccidentals_postcut = new TCanvas("c_TaggerAccidentals_postcut","c_TaggerAccidentals_postcut",900,600);
  c_TaggerAccidentals_postcut->cd();
  TH1F *h_TaggerAccidentals_postcut = new TH1F("h_TaggerAccidentals_postcut",  Form("Vertex time - RF (%s);#Deltat_{Beam-RF}(ns)",name.Data()), 300, -18, 18);
  t->Project("h_TaggerAccidentals_postcut","rf_dt","w8");
  h_TaggerAccidentals_postcut->Draw("hist");
  h_TaggerAccidentals_postcut->Write();
  c_TaggerAccidentals_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.root",name.Data()), "root");
  c_TaggerAccidentals_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_TaggerAccidentals_postcut.eps",name.Data()), "eps");
*/
/*
  // ********  Chi2 of KinFit
  TCanvas *c_kin_chisq = new TCanvas("c_kin_chisq","c_kin_chisq",900,600);
  c_kin_chisq->cd();
  TH1F *h_kin_chisq = new TH1F("h_kin_chisq", Form("%s;#chi^{2} of Kinematic Fit;Counts",name.Data()), 300, 0, 100);
  t->Project("h_kin_chisq","kin_chisq","w8*("+cutlist+"&& abs(mm2)<0.01 && kp_x4_kin.Z()<80 && kp_x4_kin.Z()>50)");
  TLine *l_kin_chisq = new TLine(25, 0, 25, h_kin_chisq->GetMaximum());
  l_kin_chisq->SetLineColor(kRed);
  l_kin_chisq->SetLineStyle(10);
  h_kin_chisq->Draw("hist");
  h_kin_chisq->Write(Form("h%s_kin_chisq",name.Data()),TObject::kWriteDelete);
  l_kin_chisq->Draw("same");

  //  //Copy h1 in a clone h1c. Set range and color for h1c
  //  TH1F *h_kin_chisq_side = (TH1F*)h_kin_chisq->Clone();
  //  h_kin_chisq_side->SetFillColor(2);
  //  h_kin_chisq_side->GetXaxis()->SetRange(h_kin_chisq->FindBin(50),h_kin_chisq->FindBin(70));
  //  h_kin_chisq_side->Draw("hist same");
  //  TH1F *h_kin_chisq_sub = (TH1F*)h_kin_chisq->Clone();
  //  h_kin_chisq_sub->SetFillColor(1);
  //  h_kin_chisq_sub->GetXaxis()->SetRange(h_kin_chisq->FindBin(0),h_kin_chisq->FindBin(20));
  //  h_kin_chisq_sub->Draw("hist same");

  h_kin_chisq->Write();
  c_kin_chisq->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq.root",name.Data()), "root");
  c_kin_chisq->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kin_chisq.eps",name.Data()), "eps");  


  // ********  Missing Mass Squared
  TCanvas *c_mm2 = new TCanvas("c_mm2","c_mm2",900,600);
  c_mm2->cd();
  // c_mm2->SetLogy();
  TH1F *h_mm2 = new TH1F("h_mm2", Form("%s; MM^{2} (GeV/c^{2})^{2}; Counts",name.Data()), 600, -0.1, 0.1);
  t->Project("h_mm2","mm2","w8*("+cutlist+")");
  TLine *l_mm21 = new TLine(0.01, 0, 0.01, h_mm2->GetMaximum());
  TLine *l_mm22 = new TLine(-0.01, 0, -0.01, h_mm2->GetMaximum());
  l_mm21->SetLineColor(kRed);
  l_mm22->SetLineColor(kRed);
  l_mm21->SetLineStyle(10);
  l_mm22->SetLineStyle(10);  
  h_mm2->Draw("hist");
  l_mm21->Draw("same");
  l_mm22->Draw("same");
  h_mm2->Write(Form("h%s_mm2",name.Data()),TObject::kWriteDelete);
  c_mm2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2.root",name.Data()), "root");
  c_mm2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2.eps",name.Data()), "eps");
*/
/*
  // ********  VertexZ
  // proton
  TCanvas *c_p_vertexz = new TCanvas("c_p_vertexz","c_p_vertexz",600,400);
  c_p_vertexz->cd();
  TH1F *h_p_vertexz = new TH1F("h_p_vertexz", Form("%s;Vertex-Z (cm) (Protons); Counts",name.Data()), 300, 45, 85);
  t->Project("h_p_vertexz","p_x4_kin.Z()","w8*("+cutlist+"&& abs(mm2)<0.01)");
  TLine *l_p_vertexz1 = new TLine(50, 0, 50, h_p_vertexz->GetMaximum());
  TLine *l_p_vertexz2 = new TLine(80, 0, 80, h_p_vertexz->GetMaximum());
  l_p_vertexz1->SetLineColor(kRed);
  l_p_vertexz2->SetLineColor(kRed);
  l_p_vertexz1->SetLineStyle(10);
  l_p_vertexz2->SetLineStyle(10);  
  h_p_vertexz->Draw("hist");
  l_p_vertexz1->Draw("same");
  l_p_vertexz2->Draw("same");
  h_p_vertexz->Write();
  c_p_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexz.root",name.Data()), "root");
  c_p_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexz.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_vertexz = new TCanvas("c_pip_vertexz","c_pip_vertexz",600,400);
  c_pip_vertexz->cd();
  TH1F *h_pip_vertexz = new TH1F("h_pip_vertexz", ";Vertex-Z (cm) (#pi^{+}); Counts", 300, 45, 85);
  t->Project("h_pip_vertexz","pip_x4_kin.Z()","w8");
  h_pip_vertexz->Draw("hist");
  c_pip_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexz.root",name.Data()), "root");
  c_pip_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexz.eps",name.Data()), "eps");

 // pim
  TCanvas *c_pim_vertexz = new TCanvas("c_pim_vertexz","c_pim_vertexz",600,400);
  c_pim_vertexz->cd();
  TH1F *h_pim_vertexz = new TH1F("h_pim_vertexz", ";Vertex-Z (cm) (#pi^{-}); Counts", 300, 45, 85);
  t->Project("h_pim_vertexz","pim_x4_kin.Z()","w8");
  h_pim_vertexz->Draw("hist");
  c_pim_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexz.root",name.Data()), "root");
  c_pim_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexz.eps",name.Data()), "eps"); 

 // kp
  TCanvas *c_kp_vertexz = new TCanvas("c_kp_vertexz","c_kp_vertexz",600,400);
  c_kp_vertexz->cd();
  TH1F *h_kp_vertexz = new TH1F("h_kp_vertexz", ";Vertex-Z (cm) (K^{+}); Counts", 300, 45, 85);
  t->Project("h_kp_vertexz","kp_x4_kin.Z()","w8");
  h_kp_vertexz->Draw("hist");
  c_kp_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexz.root",name.Data()), "root");
  c_kp_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexz.eps",name.Data()), "eps"); 

 // km
  TCanvas *c_km_vertexz = new TCanvas("c_km_vertexz","c_km_vertexz",600,400);
  c_km_vertexz->cd();
  TH1F *h_km_vertexz = new TH1F("h_km_vertexz", ";Vertex-Z (cm) (K^{-}); Counts", 300, 45, 85);
  t->Project("h_km_vertexz","km_x4_kin.Z()","w8");
  h_km_vertexz->Draw("hist");
  c_km_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexz.root",name.Data()), "root");
  c_km_vertexz->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexz.eps",name.Data()), "eps");
*/
/*
  // ********  Vertex:sqrt(X^2+Y^2)
  // proton
  TCanvas *c_p_vertexxy = new TCanvas("c_p_vertexxy","c_p_vertexxy",600,400);
  c_p_vertexxy->cd();
  c_p_vertexxy->SetLogy();
  TH1F *h_p_vertexxy = new TH1F("h_p_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (Proton); Counts", 300, 0, 2);
  t->Project("h_p_vertexxy","sqrt((p_x4_kin.X()*p_x4_kin.X())+(p_x4_kin.Y()*p_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_p_vertexxy->Draw("hist");
  c_p_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexxy.root",name.Data()), "root");
  c_p_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_vertexxy.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_vertexxy = new TCanvas("c_pip_vertexxy","c_pip_vertexxy",600,400);
  c_pip_vertexxy->cd();
  c_pip_vertexxy->SetLogy();
  TH1F *h_pip_vertexxy = new TH1F("h_pip_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (#pi^{+}); Counts", 300,  0, 2);
  t->Project("h_pip_vertexxy","sqrt((pip_x4_kin.X()*pip_x4_kin.X())+(pip_x4_kin.Y()*pip_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pip_vertexxy->Draw("hist");
  c_pip_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexxy.root",name.Data()), "root");
  c_pip_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_vertexxy.eps",name.Data()), "eps");

 // pim
  TCanvas *c_pim_vertexxy = new TCanvas("c_pim_vertexxy","c_pim_vertexxy",600,400);
  c_pim_vertexxy->cd();
  c_pim_vertexxy->SetLogy();
  TH1F *h_pim_vertexxy = new TH1F("h_pim_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (#pi^{-}); Counts", 300,  0, 2);
  t->Project("h_pim_vertexxy","sqrt((pim_x4_kin.X()*pim_x4_kin.X())+(pim_x4_kin.Y()*pim_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_pim_vertexxy->Draw("hist");
  c_pim_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexxy.root",name.Data()), "root");
  c_pim_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_vertexxy.eps",name.Data()), "eps"); 

 // kp
  TCanvas *c_kp_vertexxy = new TCanvas("c_kp_vertexxy","c_kp_vertexxy",600,400);
  c_kp_vertexxy->cd();
  c_kp_vertexxy->SetLogy();
  TH1F *h_kp_vertexxy = new TH1F("h_kp_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (K^{+}); Counts", 300,  0, 2);
  t->Project("h_kp_vertexxy","sqrt((kp_x4_kin.X()*kp_x4_kin.X())+(kp_x4_kin.Y()*kp_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_kp_vertexxy->Draw("hist");
  c_kp_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexxy.root",name.Data()), "root");
  c_kp_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_vertexxy.eps",name.Data()), "eps"); 

 // km
  TCanvas *c_km_vertexxy = new TCanvas("c_km_vertexxy","c_km_vertexxy",600,400);
  c_km_vertexxy->cd();
  c_km_vertexxy->SetLogy();
  TH1F *h_km_vertexxy = new TH1F("h_km_vertexxy", ";#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} [cm] (K^{-}); Counts", 300,  0, 2);
  t->Project("h_km_vertexxy","sqrt((km_x4_kin.X()*km_x4_kin.X())+(km_x4_kin.Y()*km_x4_kin.Y()))","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_km_vertexxy->Draw("hist");
  c_km_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexxy.root",name.Data()), "root");
  c_km_vertexxy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_vertexxy.eps",name.Data()), "eps");  
*/
/*
  // ********  Missing Energy
  TCanvas *c_me = new TCanvas("c_me","c_me",600,400);
  c_me->cd();
  TH1F *h_me = new TH1F("h_me", ";Missing Energy (GeV);Counts", 600, -3, 3);
  t->Project("h_me","me","w8");
  h_me->Draw("hist");
  h_me->Write(Form("h%s_me",name.Data()),TObject::kWriteDelete);
  c_me->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_me.root", name.Data()), "root");
  c_me->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_me.eps",name.Data()), "eps"); 

  // ********  Total Transverse Momentum
  TCanvas *c_pt = new TCanvas("c_pt","c_pt",600,400);
  c_pt->cd();
  c_pt->SetLogy();
  TH1F *h_pt = new TH1F("h_pt", ";P_{t} total (GeV);Counts", 600, 0, 3);
  t->Project("h_pt","pt","w8");
  h_pt->Draw("hist");
  h_pt->Write(Form("h%s_pt",name.Data()),TObject::kWriteDelete);
  c_pt->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pt.root",name.Data()), "root");
  c_pt->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pt.eps",name.Data()), "eps"); 
*/

  // ********  Beam Energy
  TCanvas *c_beam_e = new TCanvas("c_beam_e", "c_beam_e", 900, 600);
  c_beam_e->cd();
  TH1F *h_beam_e = new TH1F("h_beam_e", Form("%s;Beam Energy (GeV);Counts",name.Data()), 100, 3.0, 11.6);
  t->Project("h_beam_e","beam_p4_kin.E()","w8");
  h_beam_e->Draw("hist");
  h_beam_e->Write(Form("h%s_beam_e",name.Data()),TObject::kWriteDelete);
  c_beam_e->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beam_e.root", name.Data()), "root");
  c_beam_e->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beam_e.eps",name.Data()), "eps");

  // // ********  proton momentum transfer (-t)
  // TCanvas *c_t_kin_precut = new TCanvas("c_t_kin_precut","c_t_kin_precut",900,600);
  // c_t_kin_precut->cd();
  // TH1F *h_t_kin_precut = new TH1F("h_t_kin_precut", ";-t (GeV/c^{2});Counts", 300, 0, 10);
  // t->Project("h_t_kin_precut","-t_kin","w8*(kin_chisq<30 && abs(mm2)<0.015)");
  // h_t_kin_precut->Draw("hist");
  // c_t_kin_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin_precut.root",name.Data()), "root");
  // c_t_kin_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin_precut.eps",name.Data()), "eps"); 

  // post PhiMass cut
  TCanvas *c_t_kin = new TCanvas("c_t_kin","c_t_kin",900,600);
  c_t_kin->cd();
  c_t_kin->SetLogy();
  TH1F *h_t_kin = new TH1F("h_t_kin", Form("%s;-t (GeV/c)^{2});Counts",name.Data()), 300, 0, 10);
  t->Project("h_t_kin","-t_kin","w8*(kpkm_mf>1.005 && kpkm_mf<1.035)");
  h_t_kin->Fit("expo", "R", "", 0.5, 3);
  h_t_kin->Draw("e");
  h_t_kin->Write(Form("h%s_t_kin",name.Data()),TObject::kWriteDelete);
  c_t_kin->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin.root", name.Data()), "root");
  c_t_kin->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_t_kin.eps",name.Data()), "eps");
 
  // ********  Eg vs. -t
  TCanvas *c_beamevst = new TCanvas("c_beamevst", "c_beamevst", 900, 600);
  c_beamevst->cd();
  TH2D *h_beamevst = new TH2D("h_beamevst", Form("%s;E_{#gamma} [GeV];-t [GeV/c)^{2}]",name.Data()), 100, 3.0, 12, 100, 0, 10);
  t->Project("h_beamevst","-t_kin:beam_p4_kin.E()","");//w8*(kpkm_mf>1.005 && kpkm_mf<1.035)
  h_beamevst->Draw("colz");
  h_beamevst->Write(Form("h%s_beamevst",name.Data()),TObject::kWriteDelete);
  c_beamevst->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beamevst.root", name.Data()), "root");
  c_beamevst->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_beamevst.eps",name.Data()), "eps");

/*
  // ********  Missing mass squared with Pions as Kions hypothesis
  TCanvas *c_mm2_piask = new TCanvas("c_mm2_piask","c_mm2_piask",600,400);
  c_mm2_piask->cd();
  TH1F *h_mm2_piask = new TH1F("h_mm2_piask", ";MM^{2} (Pions As Kaons) (GeV/c^{2})^{2};Counts", 300, -10, 4);
  t->Project("h_mm2_piask","mm2_piask","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_mm2_piask->Draw("hist");
  c_mm2_piask->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2_piask.root",name.Data()), "root");
  c_mm2_piask->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2_piask.eps",name.Data()), "eps");   

  TCanvas *c_mm2vsmm2_piask = new TCanvas("c_mm2vsmm2_piask","c_mm2vsmm2_piask",600,400);
  c_mm2vsmm2_piask->cd();
  TH2D *h_mm2vsmm2_piask = new TH2D("h_mm2vsmm2_piask", ";MM^{2} (GeV/c^{2})^{2};MM^{2} (Pions As Kaons) (GeV/c^{2})^{2}", 300, -0.06, 0.06, 300, -10, 5);
  t->Project("h_mm2vsmm2_piask","mm2_piask:mm2","w8*(abs(mm2)<0.06 && kin_chisq<30)");
  h_mm2vsmm2_piask->Draw("colz");
  c_mm2vsmm2_piask->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2vsmm2_piask.root",name.Data()), "root");
  c_mm2vsmm2_piask->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_mm2vsmm2_piask.eps",name.Data()), "eps");     
*/
/*
  // **************** P vs. theta
  // proton
  TCanvas *c_p_ptheta = new TCanvas("c_p_ptheta", "c_p_ptheta", 600, 400);
  c_p_ptheta->cd();
  TH2D *h_p_ptheta = new TH2D("h_p_ptheta", Form("proton (%s); #theta#circ ;p (GeV/c)",name.Data()), 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_p_ptheta", "p_p4_kin.P():p_p4_kin.Theta()*TMath::RadToDeg()", "");
  cout << "h_p_ptheta = " << h_p_ptheta << endl;
  h_p_ptheta->Draw("colz");
  h_p_ptheta->Write(Form("h%s_p_ptheta",name.Data()),TObject::kWriteDelete);
  c_p_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_ptheta.root",name.Data()), "root");
  c_p_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_ptheta.eps",name.Data()), "eps");

  // pip
  TCanvas *c_pip_ptheta = new TCanvas("c_pip_ptheta", "c_pip_ptheta", 600, 400);
  c_pip_ptheta->cd();
  TH2D *h_pip_ptheta = new TH2D("h_pip_ptheta", Form("#pi^{+} (%s); #theta#circ ;p (GeV/c)",name.Data()), 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_pip_ptheta", "pip_p4_kin.P():pip_p4_kin.Theta()*TMath::RadToDeg()", "");
  cout << "h_pip_ptheta = " << h_pip_ptheta << endl;
  h_pip_ptheta->Draw("colz");
  h_pip_ptheta->Write(Form("h%s_pip_ptheta",name.Data()),TObject::kWriteDelete);
  c_pip_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_ptheta.root",name.Data()), "root");
  c_pip_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_ptheta.eps",name.Data()), "eps");

  // pim
  TCanvas *c_pim_ptheta = new TCanvas("c_pim_ptheta", "c_pim_ptheta", 600, 400);
  c_pim_ptheta->cd();
  TH2D *h_pim_ptheta = new TH2D("h_pim_ptheta", Form("#pi^{-} (%s); #theta#circ ;p (GeV/c)",name.Data()), 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_pim_ptheta", "pim_p4_kin.P():pim_p4_kin.Theta()*TMath::RadToDeg()", "");
  cout << "h_pim_ptheta = " << h_pim_ptheta << endl;
  h_pim_ptheta->Draw("colz");
  h_pim_ptheta->Write(Form("h%s_pim_ptheta",name.Data()),TObject::kWriteDelete);
  c_pim_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_ptheta.root",name.Data()), "root");
  c_pim_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_ptheta.eps",name.Data()), "eps");

  // kp
  TCanvas *c_kp_ptheta = new TCanvas("c_kp_ptheta", "c_kp_ptheta", 600, 400);
  c_kp_ptheta->cd();
  TH2D *h_kp_ptheta = new TH2D("h_kp_ptheta", Form("K^{+} (%s); #theta#circ ;p (GeV/c)",name.Data()), 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_kp_ptheta", "kp_p4_kin.P():kp_p4_kin.Theta()*TMath::RadToDeg()", "");
  cout << "h_kp_ptheta = " << h_kp_ptheta << endl;
  h_kp_ptheta->Draw("colz");
  h_kp_ptheta->Write(Form("h%s_kp_ptheta",name.Data()),TObject::kWriteDelete);
  c_kp_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_ptheta.root",name.Data()), "root");
  c_kp_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_ptheta.eps",name.Data()), "eps");
  
  // km
  TCanvas *c_km_ptheta = new TCanvas("c_km_ptheta", "c_km_ptheta", 600, 400);
  c_km_ptheta->cd();
  TH2D *h_km_ptheta = new TH2D("h_km_ptheta", Form("K^{-} (%s); #theta#circ ;p (GeV/c)",name.Data()), 140, 0.0, 140.0, 250, 0.0, 10.0);
  t->Project("h_km_ptheta", "km_p4_kin.P():km_p4_kin.Theta()*TMath::RadToDeg()", "");
  cout << "h_km_ptheta = " << h_km_ptheta << endl;
  h_km_ptheta->Draw("colz");
  h_km_ptheta->Write(Form("h%s_km_ptheta",name.Data()),TObject::kWriteDelete);
  c_km_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_ptheta.root",name.Data()), "root");
  c_km_ptheta->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_ptheta.eps",name.Data()), "eps");
*/
/*
  // **************** dt vs. p
  // proton
  // TOF
  TCanvas *c_p_dttof = new TCanvas("c_p_dttof", "c_p_dttof", 600, 400);
  c_p_dttof->cd();
  TH2D *h_p_dttof = new TH2D("h_p_dttof", Form("proton (%s); p (GeV/c); TOF #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.6, 0.6);
  t->Project("h_p_dttof", "p_dttof:p_p4_kin.P()","");
  cout << "h_p_dttof = " << h_p_dttof << endl;
  TLine *l_p_dttof1 = new TLine(0, 0.27, 8, 0.27);
  TLine *l_p_dttof2 = new TLine(0, -0.27, 8, -0.27); 
  l_p_dttof1->SetLineColor(kRed);
  l_p_dttof2->SetLineColor(kRed);
  l_p_dttof1->SetLineStyle(10);
  l_p_dttof2->SetLineStyle(10); 
  h_p_dttof->Draw("colz");
  l_p_dttof1->Draw("same");
  l_p_dttof2->Draw("same");
  h_p_dttof->Write(Form("h%s_p_dttof",name.Data()),TObject::kWriteDelete);
  c_p_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dttof.root",name.Data()), "root");
  c_p_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dttof.eps",name.Data()), "eps");

  // BCAL
  TCanvas *c_p_dtbcal = new TCanvas("c_p_dtbcal", "c_p_dtbcal", 600, 400);
  c_p_dtbcal->cd();
  TH2D *h_p_dtbcal = new TH2D("h_p_dtbcal", Form("proton (%s); p (GeV/c); BCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_p_dtbcal", "p_dtbcal:p_p4_kin.P()","");
  cout << "h_p_dtbcal = " << h_p_dtbcal << endl;
  h_p_dtbcal->Draw("colz");
  h_p_dtbcal->Write(Form("h%s_p_dtbcal",name.Data()),TObject::kWriteDelete);
  c_p_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtbcal.root",name.Data()), "root");
  c_p_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_p_dtfcal = new TCanvas("c_p_dtfcal", "c_p_dtfcal", 600, 400);
  c_p_dtfcal->cd();
  TH2D *h_p_dtfcal = new TH2D("h_p_dtfcal", Form("proton (%s); p (GeV/c); FCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_p_dtfcal", "p_dtfcal:p_p4_kin.P()","");
  cout << "h_p_dtfcal = " << h_p_dtfcal << endl;
  h_p_dtfcal->Draw("colz");
  h_p_dtfcal->Write(Form("h%s_p_dtfcal",name.Data()),TObject::kWriteDelete);
  c_p_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtfcal.root",name.Data()), "root");
  c_p_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_p_dtfcal.eps",name.Data()), "eps");

  // pip
  // TOF
  TCanvas *c_pip_dttof = new TCanvas("c_pip_dttof", "c_pip_dttof", 600, 400);
  c_pip_dttof->cd();
  TH2D *h_pip_dttof = new TH2D("h_pip_dttof", Form("#pi^{+} (%s); p (GeV/c); TOF #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.5, 0.5);
  t->Project("h_pip_dttof", "pip_dttof:pip_p4_kin.P()","");
  cout << "h_pip_dttof = " << h_pip_dttof << endl;
  TLine *l_pip_dttof1 = new TLine(0, 0.15, 8, 0.15);
  TLine *l_pip_dttof2 = new TLine(0, -0.15, 8, -0.15); 
  l_pip_dttof1->SetLineColor(kRed);
  l_pip_dttof2->SetLineColor(kRed);
  l_pip_dttof1->SetLineStyle(10);
  l_pip_dttof2->SetLineStyle(10); 
  h_pip_dttof->Draw("colz");
  l_pip_dttof1->Draw("same");
  l_pip_dttof2->Draw("same");
  h_pip_dttof->Write(Form("h%s_pip_dttof",name.Data()),TObject::kWriteDelete);
  c_pip_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dttof.root",name.Data()), "root");
  c_pip_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_pip_dtbcal = new TCanvas("c_pip_dtbcal", "c_pip_dtbcal", 600, 400);
  c_pip_dtbcal->cd();
  TH2D *h_pip_dtbcal = new TH2D("h_pip_dtbcal", Form("#pi^{+} (%s); p (GeV/c); BCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_pip_dtbcal", "pip_dtbcal:pip_p4_kin.P()","");
  cout << "h_pip_dtbcal = " << h_pip_dtbcal << endl;
  TLine *l_pip_dtbcal1 = new TLine(0, 0.85, 8, 0.85);
  TLine *l_pip_dtbcal2 = new TLine(0, -0.85, 8, -0.85); 
  l_pip_dtbcal1->SetLineColor(kRed);
  l_pip_dtbcal2->SetLineColor(kRed);
  l_pip_dtbcal1->SetLineStyle(10);
  l_pip_dtbcal2->SetLineStyle(10); 
  h_pip_dtbcal->Draw("colz");
  l_pip_dtbcal1->Draw("same");
  l_pip_dtbcal2->Draw("same");
  h_pip_dtbcal->Write(Form("h%s_pip_dtbcal",name.Data()),TObject::kWriteDelete);
  c_pip_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtbcal.root",name.Data()), "root");
  c_pip_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_pip_dtfcal = new TCanvas("c_pip_dtfcal", "c_pip_dtfcal", 600, 400);
  c_pip_dtfcal->cd();
  TH2D *h_pip_dtfcal = new TH2D("h_pip_dtfcal", Form("#pi^{+} (%s); p (GeV/c); FCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_pip_dtfcal", "pip_dtfcal:pip_p4_kin.P()","");
  cout << "h_pip_dtfcal = " << h_pip_dtfcal << endl;
  h_pip_dtfcal->Draw("colz");
  h_pip_dtfcal->Write(Form("h%s_pip_dtfcal",name.Data()),TObject::kWriteDelete);
  c_pip_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtfcal.root",name.Data()), "root");
  c_pip_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pip_dtfcal.eps",name.Data()), "eps");

  // pim
  // TOF
  TCanvas *c_pim_dttof = new TCanvas("c_pim_dttof", "c_pim_dttof", 600, 400);
  c_pim_dttof->cd();
  TH2D *h_pim_dttof = new TH2D("h_pim_dttof", Form("#pi^{-} (%s); p (GeV/c); TOF #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.5, 0.5);
  t->Project("h_pim_dttof", "pim_dttof:pim_p4_kin.P()","");
  cout << "h_pim_dttof = " << h_pim_dttof << endl;
  h_pim_dttof->Draw("colz");
  h_pim_dttof->Write(Form("h%s_pim_dttof",name.Data()),TObject::kWriteDelete);
  c_pim_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dttof.root",name.Data()), "root");
  c_pim_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_pim_dtbcal = new TCanvas("c_pim_dtbcal", "c_pim_dtbcal", 600, 400);
  c_pim_dtbcal->cd();
  TH2D *h_pim_dtbcal = new TH2D("h_pim_dtbcal", Form("#pi^{-} (%s); p (GeV/c); BCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-1.0, 1.0);
  t->Project("h_pim_dtbcal", "pim_dtbcal:pim_p4_kin.P()","");
  cout << "h_pim_dtbcal = " << h_pim_dtbcal << endl;
  h_pim_dtbcal->Draw("colz");
  h_pim_dtbcal->Write(Form("h%s_pim_dtbcal",name.Data()),TObject::kWriteDelete);
  c_pim_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtbcal.root",name.Data()), "root");
  c_pim_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_pim_dtfcal = new TCanvas("c_pim_dtfcal", "c_pim_dtfcal", 600, 400);
  c_pim_dtfcal->cd();
  TH2D *h_pim_dtfcal = new TH2D("h_pim_dtfcal", Form("#pi^{-} (%s); p (GeV/c); FCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-2.0, 2.0);
  t->Project("h_pim_dtfcal", "pim_dtfcal:pim_p4_kin.P()","");
  cout << "h_pim_dtfcal = " << h_pim_dtfcal << endl;
  h_pim_dtfcal->Draw("colz");
  h_pim_dtfcal->Write(Form("h%s_pim_dtfcal",name.Data()),TObject::kWriteDelete);
  c_pim_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtfcal.root",name.Data()), "root");
  c_pim_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pim_dtfcal.eps",name.Data()), "eps");

  // kp
  // TOF
  TCanvas *c_kp_dttof = new TCanvas("c_kp_dttof", "c_kp_dttof", 600, 400);
  c_kp_dttof->cd();
  TH2D *h_kp_dttof = new TH2D("h_kp_dttof", Form("K^{+} (%s); p (GeV/c); TOF #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.3, 0.3);
  t->Project("h_kp_dttof", "kp_dttof:kp_p4_kin.P()","");
  cout << "h_kp_dttof = " << h_kp_dttof << endl;
  TLine *l_kp_dttof1 = new TLine(0, 0.135, 8, 0.135);
  TLine *l_kp_dttof2 = new TLine(0, -0.135, 8, -0.135); 
  l_kp_dttof1->SetLineColor(kRed);
  l_kp_dttof2->SetLineColor(kRed);
  l_kp_dttof1->SetLineStyle(10);
  l_kp_dttof2->SetLineStyle(10);
  h_kp_dttof->Draw("colz");
  l_kp_dttof1->Draw("same");
  l_kp_dttof2->Draw("same");
  h_kp_dttof->Write(Form("h%s_kp_dttof",name.Data()),TObject::kWriteDelete);
  c_kp_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dttof.root",name.Data()), "root");
  c_kp_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_kp_dtbcal = new TCanvas("c_kp_dtbcal", "c_kp_dtbcal", 600, 400);
  c_kp_dtbcal->cd();
  TH2D *h_kp_dtbcal = new TH2D("h_kp_dtbcal", Form("K^{+} (%s); p (GeV/c); BCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.75, 0.75);
  t->Project("h_kp_dtbcal", "kp_dtbcal:kp_p4_kin.P()","");
  cout << "h_kp_dtbcal = " << h_kp_dtbcal << endl;
  TLine *l_kp_dtbcal1 = new TLine(0, 0.4125, 8, 0.4125);
  TLine *l_kp_dtbcal2 = new TLine(0, -0.4125, 8, -0.4125); 
  l_kp_dtbcal1->SetLineColor(kRed);
  l_kp_dtbcal2->SetLineColor(kRed);
  l_kp_dtbcal1->SetLineStyle(10);
  l_kp_dtbcal2->SetLineStyle(10);  
  h_kp_dtbcal->Draw("colz");
  l_kp_dtbcal1->Draw("same");
  l_kp_dtbcal2->Draw("same");
  h_kp_dtbcal->Write(Form("h%s_kp_dtbcal",name.Data()),TObject::kWriteDelete);
  c_kp_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtbcal.root",name.Data()), "root");
  c_kp_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_kp_dtfcal = new TCanvas("c_kp_dtfcal", "c_kp_dtfcal", 600, 400);
  c_kp_dtfcal->cd();
  TH2D *h_kp_dtfcal = new TH2D("h_kp_dtfcal", Form("K^{+} (%s); p (GeV/c); FCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-2.5, 2.5);
  t->Project("h_kp_dtfcal", "kp_dtfcal:kp_p4_kin.P()","");
  cout << "h_kp_dtfcal = " << h_kp_dtfcal << endl;
  TLine *l_kp_dtfcal1 = new TLine(0, 1.0, 8, 1.0);
  TLine *l_kp_dtfcal2 = new TLine(0, -1.0, 8, -1.0); 
  l_kp_dtfcal1->SetLineColor(kRed);
  l_kp_dtfcal2->SetLineColor(kRed);
  l_kp_dtfcal1->SetLineStyle(10);
  l_kp_dtfcal2->SetLineStyle(10);  
  h_kp_dtfcal->Draw("colz");
  l_kp_dtfcal1->Draw("same");
  l_kp_dtfcal2->Draw("same");
  h_kp_dtfcal->Write(Form("h%s_kp_dtfcal",name.Data()),TObject::kWriteDelete);
  c_kp_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtfcal.root",name.Data()), "root");
  c_kp_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kp_dtfcal.eps",name.Data()), "eps");

  // km
  // TOF
  TCanvas *c_km_dttof = new TCanvas("c_km_dttof", "c_km_dttof", 600, 400);
  c_km_dttof->cd();
  TH2D *h_km_dttof = new TH2D("h_km_dttof", Form("K^{-} (%s); p (GeV/c); TOF #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.3, 0.3);
  t->Project("h_km_dttof", "km_dttof:km_p4_kin.P()","");
  cout << "h_km_dttof = " << h_km_dttof << endl;
  h_km_dttof->Draw("colz");
  h_km_dttof->Write(Form("h%s_km_dttof",name.Data()),TObject::kWriteDelete);
  c_km_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dttof.root",name.Data()), "root");
  c_km_dttof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dttof.eps",name.Data()), "eps");
  // BCAL
  TCanvas *c_km_dtbcal = new TCanvas("c_km_dtbcal", "c_km_dtbcal", 600, 400);
  c_km_dtbcal->cd();
  TH2D *h_km_dtbcal = new TH2D("h_km_dtbcal", Form("K^{-} (%s); p (GeV/c); BCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-0.75, 0.75);
  t->Project("h_km_dtbcal", "km_dtbcal:km_p4_kin.P()","");
  cout << "h_km_dtbcal = " << h_km_dtbcal << endl;
  h_km_dtbcal->Draw("colz");
  h_km_dtbcal->Write(Form("h%s_km_dtbcal",name.Data()),TObject::kWriteDelete);
  c_km_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtbcal.root",name.Data()), "root");
  c_km_dtbcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtbcal.eps",name.Data()), "eps");
  // FCAL
  TCanvas *c_km_dtfcal = new TCanvas("c_km_dtfcal", "c_km_dtfcal", 600, 400);
  c_km_dtfcal->cd();
  TH2D *h_km_dtfcal = new TH2D("h_km_dtfcal", Form("K^{-} (%s); p (GeV/c); FCAL #Delta T (ns)",name.Data()), 250, 0.0, 10.0, 500,-2.5, 2.5);
  t->Project("h_km_dtfcal", "km_dtfcal:km_p4_kin.P()","");
  cout << "h_km_dtfcal = " << h_km_dtfcal << endl;
  h_km_dtfcal->Draw("colz");
  h_km_dtfcal->Write(Form("h%s_km_dtfcal",name.Data()),TObject::kWriteDelete);
  c_km_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtfcal.root",name.Data()), "root");
  c_km_dtfcal->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_km_dtfcal.eps",name.Data()), "eps");  
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
  c_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.root",name.Data()), "root");
  c_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass.eps",name.Data()), "eps");

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
  c_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.root",name.Data()), "root");
  c_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass.eps",name.Data()), "eps");

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
  c_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.root",name.Data()), "root");
  c_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass.eps",name.Data()), "eps");
*/
/*
  // ******** Phi vs. fo

  TCanvas *c_h2_PhiVsfoMass_KinFit = new TCanvas("c_h2_PhiVsfoMass_KinFit","c_h2_PhiVsfoMass_KinFit",600,400);
  c_h2_PhiVsfoMass_KinFit->cd();
  TH2D *h2_PhiVsfoMass_KinFit = (TH2D*) f->Get("h2_PhiVsfoMass_KinFit");
  h2_PhiVsfoMass_KinFit->Draw("colz");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsfoMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_PhiVsfoMass_postcut = new TCanvas("c_h2_PhiVsfoMass_postcut","c_h2_PhiVsfoMass_postcut",600,400);
  // c_h2_PhiVsfoMass_postcut->cd();
  // TH2D *h2_PhiVsfoMass_postcut = (TH2D*) f->Get("h2_PhiVsfoMass_postcut");
  // h2_PhiVsfoMass_postcut->Draw("colz");
  // c_h2_PhiVsfoMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.root",name.Data()), "root");
  // c_h2_PhiVsfoMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsfoMass_postcut.eps",name.Data()), "eps");

  // ******** Phi vs. Y

  TCanvas *c_h2_PhiVsYMass_KinFit = new TCanvas("c_h2_PhiVsYMass_KinFit","c_h2_PhiVsYMass_KinFit",600,400);
  c_h2_PhiVsYMass_KinFit->cd();
  TH2D *h2_PhiVsYMass_KinFit = (TH2D*) f->Get("h2_PhiVsYMass_KinFit");
  h2_PhiVsYMass_KinFit->Draw("colz");
  c_h2_PhiVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.root",name.Data()), "root");
  c_h2_PhiVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_PhiVsYMass_postcut = new TCanvas("c_h2_PhiVsYMass_postcut","c_h2_PhiVsYMass_postcut",600,400);
  // c_h2_PhiVsYMass_postcut->cd();
  // TH2D *h2_PhiVsYMass_postcut = (TH2D*) f->Get("h2_PhiVsYMass_postcut");
  // h2_PhiVsYMass_postcut->Draw("colz");
  // c_h2_PhiVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.root",name.Data()), "root");
  // c_h2_PhiVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_PhiVsYMass_postcut.eps",name.Data()), "eps");


  // ******** fo vs. Y

  TCanvas *c_h2_foVsYMass_KinFit = new TCanvas("c_h2_foVsYMass_KinFit","c_h2_foVsYMass_KinFit",600,400);
  c_h2_foVsYMass_KinFit->cd();
  TH2D *h2_foVsYMass_KinFit = (TH2D*) f->Get("h2_foVsYMass_KinFit");
  h2_foVsYMass_KinFit->Draw("colz");
  c_h2_foVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.root",name.Data()), "root");
  c_h2_foVsYMass_KinFit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_KinFit.eps",name.Data()), "eps");

  // TCanvas *c_h2_foVsYMass_postcut = new TCanvas("c_h2_foVsYMass_postcut","c_h2_foVsYMass_postcut",600,400);
  // c_h2_foVsYMass_postcut->cd();
  // TH2D *h2_foVsYMass_postcut = (TH2D*) f->Get("h2_foVsYMass_postcut");
  // h2_foVsYMass_postcut->Draw("colz");
  // c_h2_foVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.root",name.Data()), "root");
  // c_h2_foVsYMass_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_h2_foVsYMass_postcut.eps",name.Data()), "eps");
*/

  // +++ Phi Mass
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", Form("%s;m_{K^{+}K^{-}} [GeV/c^{2}];Counts",name.Data()), 100, 0.98, 1.2);
  t->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");//+cutlist+" && kin_chisq<25)"
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 900, 600);
  // TLine *lphicutmin = new TLine(1.005, 0, 1.005, h_PhiMass_postcuts->GetMaximum());
  // lphicutmin->SetLineColor(kRed);
  // lphicutmin->SetLineStyle(10);
  // TLine *lphicutmax = new TLine(1.035, 0, 1.035, h_PhiMass_postcuts->GetMaximum());
  // lphicutmax->SetLineColor(kRed);
  // lphicutmax->SetLineStyle(10);
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");
  // lphicutmin->Draw("same");
  // lphicutmax->Draw("same");
  
  // // Van Hove cut
  // TH1F *h_PhiMass_vhcut = new TH1F("h_PhiMass_vhcut", "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2];m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.98, 1.2);
  // h_PhiMass_vhcut->SetMarkerColor(kRed);
  // t->Project("h_PhiMass_vhcut", "kpkm_mf", "w8*(kpkm_uni && vhtheta>0.4 && vhtheta<1.5 && vhphi>1.55 && vhphi<2.2)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_PhiMass_vhcut->Draw("same");  
  // TLegend *l_PhiMass = new TLegend(0.55, 0.75, 0.85, 0.88);
  // l_PhiMass->SetFillColor(kWhite);
  // l_PhiMass->SetLineColor(kWhite);
  // l_PhiMass->AddEntry(h_PhiMass_postcuts, "#chi^{2}<25", "p");
  // l_PhiMass->AddEntry(h_PhiMass_vhcut, "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2]", "p");
  // l_PhiMass->Draw("same");

  // // Chi2pi cut
  // TH1F *h_PhiMass_05chi2pi = new TH1F("h_PhiMass_05chi2pi", "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.5*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}};m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.98, 1.2);
  // h_PhiMass_05chi2pi->SetMarkerColor(kRed);
  // t->Project("h_PhiMass_05chi2pi", "kpkm_mf", "w8*(kpkm_uni && kin_chisq<0.5*chi2pi)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_PhiMass_05chi2pi->Draw("same");
  // TH1F *h_PhiMass_01chi2pi = new TH1F("h_PhiMass_01chi2pi", "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.1*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}};m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.98, 1.2);
  // h_PhiMass_01chi2pi->SetMarkerColor(kMagenta);
  // t->Project("h_PhiMass_01chi2pi", "kpkm_mf", "w8*(kpkm_uni && kin_chisq<0.1*chi2pi)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_PhiMass_01chi2pi->Draw("same");
  // TH1F *h_PhiMass_001chi2pi = new TH1F("h_PhiMass_001chi2pi", "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.01*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}};m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.98, 1.2);
  // h_PhiMass_001chi2pi->SetMarkerColor(kGreen);
  // t->Project("h_PhiMass_001chi2pi", "kpkm_mf", "w8*(kpkm_uni && kin_chisq<0.01*chi2pi)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_PhiMass_001chi2pi->Draw("same");  
  // TLegend *l_PhiMass = new TLegend(0.7, 0.7, 0.86, 0.87);
  // l_PhiMass->SetFillColor(kWhite);
  // l_PhiMass->SetLineColor(kWhite);
  // l_PhiMass->AddEntry(h_PhiMass_postcuts, "#chi^{2}<25", "p");
  // l_PhiMass->AddEntry(h_PhiMass_05chi2pi, "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.5*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}}", "p");
  // l_PhiMass->AddEntry(h_PhiMass_01chi2pi, "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.1*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}}", "p");
  // l_PhiMass->AddEntry(h_PhiMass_001chi2pi, "#chi^{2}_{K^{+}K^{-}#pi^{+}#pi^{-}}<0.01*#chi^{2}_{#pi^{+}#pi^{-}#pi^{+}#pi^{-}}", "p");
  // l_PhiMass->Draw("same");
  // // h_PhiMass_sub->Add(h_PhiMass_postcuts, h_PhiMass_sideband, 1., -1.);


  c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcuts.root", name.Data()), "root");
  c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_PhiMass_postcuts.eps", name.Data()), "eps");

  // +++ fo Mass
  TCanvas *c_foMass_postcuts = new TCanvas("c_foMass_postcuts", "c_foMass_postcuts", 900, 600);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", Form("%s;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",name.Data()), 100, 0.3, 1.2);
  t->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pippim_mf-0.97)<0.1)");
  c_foMass_postcuts->cd();
  h_foMass_postcuts->Draw("e");

  // // Van Hove cut
  // TH1F *h_foMass_vhcut = new TH1F("h_foMass_vhcut", "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2];m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 0.3, 1.2);
  // h_foMass_vhcut->SetMarkerColor(kRed);
  // t->Project("h_foMass_vhcut", "pippim_mf", "w8*(pippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && vhtheta>0.4 && vhtheta<1.5 && vhphi>1.55 && vhphi<2.2)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_foMass_vhcut->Draw("same");  
  // TLegend *l_foMass = new TLegend(0.55, 0.75, 0.85, 0.88);
  // l_foMass->SetFillColor(kWhite);
  // l_foMass->SetLineColor(kWhite);
  // l_foMass->AddEntry(h_foMass_postcuts, "#chi^{2}<25", "p");
  // l_foMass->AddEntry(h_foMass_vhcut, "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2]", "p");
  // l_foMass->Draw("same");

  // // Chi2 side band subtraction
  // TH1F *h_foMass_sideband = new TH1F("h_foMass_sideband", "50<#chi^{2}<70;m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 0.3, 1.2);
  // h_foMass_sideband->SetLineColor(2);
  // t->Project("h_foMass_sideband", "pippim_mf", "w8*(pippim_uni && "+cutlist+" && kin_chisq>50 && kin_chisq<70 && kpkm_mf>1.005 && kpkm_mf<1.035)");
  // h_foMass_sideband->Draw("hist same");
  // TH1F *h_foMass_sub = new TH1F("h_foMass_sub", "side band subtracted;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 0.3, 1.2);
  // h_foMass_sub->Add(h_foMass_postcuts, h_foMass_sideband, 1., -1.);
  // h_foMass_sub->SetLineColor(4);
  // h_foMass_sub->Draw("hist same");
  // TLegend *l_foMass = new TLegend(0.7, 0.7, 0.86, 0.87);
  // l_foMass->SetFillColor(kWhite);
  // l_foMass->SetLineColor(kWhite);
  // l_foMass->AddEntry(h_foMass_postcuts, "#chi^{2}<25", "l");
  // l_foMass->AddEntry(h_foMass_sideband, "50<#chi^{2}<70", "l");
  // l_foMass->AddEntry(h_foMass_sub, "subtracted", "l");
  // l_foMass->Draw();

  c_foMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcuts.root", name.Data()), "root");
  c_foMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_foMass_postcuts.eps", name.Data()), "eps");

  // +++ Y Mass
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts", "c_YMass_postcuts", 900, 600);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_postcuts", Form("%s 0.87<f_{0}<1.07;m_{#phi f_{0}} [GeV/c^{2}];Counts",name.Data()), 100, 1.6, 3.2);
  TH1F *h_YMass_focut = new TH1F("h_YMass_focut", Form("%s;m_{#phi f_{0}} [GeV/c^{2}];Counts",name.Data()), 100, 1.6, 3.2);
  // TH1F *h_YMass_deltacut = new TH1F("h_YMass_deltacut", Form("%s 1<#Delta^{++}<1.4;m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts",name.Data()), 100, 1.6, 3.2);
  t->Project("h_YMass_postcuts", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035)");
  t->Project("h_YMass_focut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pippim_mf-0.97)<0.1)");
  // t->Project("h_YMass_deltacut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pipp_mf-1.2)>0.2)");
  // TH1F *h_YMass_postcuts = new TH1F("h_YMass_postcuts", ";m_{#phif_{0}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pippim_mf-0.99)<0.1)");
  c_YMass_postcuts->cd();
  h_YMass_postcuts->Draw("e"); //"hist"
  h_YMass_focut->SetLineColor(kBlue);
  h_YMass_focut->Draw("same"); //"hist"

  // // Van Hove cut
  // TH1F *h_YMass_vhcut = new TH1F("h_YMass_vhcut", "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // h_YMass_vhcut->SetMarkerColor(kRed);
  // t->Project("h_YMass_vhcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && vhtheta>0.4 && vhtheta<1.5 && vhphi>1.55 && vhphi<2.2)"); // && "+cutlist+" && kin_chisq>50 && kin_chisq<70
  // h_YMass_vhcut->Draw("same");  
  // TLegend *l_YMass = new TLegend(0.55, 0.75, 0.85, 0.88);
  // l_YMass->SetFillColor(kWhite);
  // l_YMass->SetLineColor(kWhite);
  // l_YMass->AddEntry(h_YMass_postcuts, "#chi^{2}<25", "p");
  // l_YMass->AddEntry(h_YMass_vhcut, "#theta_{VH} #in [0.4,1.5] & #phi_{VH} #in [1.55,2.2]", "p");
  // l_YMass->Draw("same");


  // // Chi2 side band subtraction
  // TH1F *h_YMass_sideband = new TH1F("h_YMass_sideband", "50<#chi^{2}<70;m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // h_YMass_sideband->SetLineColor(2);
  // t->Project("h_YMass_sideband", "kpkmpippim_mf", "w8*(kpkmpippim_uni && "+cutlist+" && kin_chisq>50 && kin_chisq<70 && kpkm_mf>1.005 && kpkm_mf<1.035)");
  // h_YMass_sideband->Draw("hist same");
  // TH1F *h_YMass_sub = new TH1F("h_YMass_sub", "side band subtracted;m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // h_YMass_sub->Add(h_YMass_postcuts, h_YMass_sideband, 1., -1.);
  // h_YMass_sub->SetLineColor(4);
  // h_YMass_sub->Draw("hist same");
  // TLegend *l_YMass = new TLegend(0.7, 0.7, 0.86, 0.87);
  // l_YMass->SetFillColor(kWhite);
  // l_YMass->SetLineColor(kWhite);
  // l_YMass->AddEntry(h_YMass_postcuts, "#chi^{2}<25", "l");
  // l_YMass->AddEntry(h_YMass_sideband, "50<#chi^{2}<70", "l");
  // l_YMass->AddEntry(h_YMass_sub, "subtracted", "l");
  // l_YMass->Draw();


  c_YMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcuts.root", name.Data()), "root");
  c_YMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_YMass_postcuts.eps", name.Data()), "eps");

/*
  //3 +++ kppim
  TH1F *h_kppim = new TH1F("h_kppim", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kppim", "kppim_mf", "w8*(kppim_uni)");  
  TCanvas *c_kppim = new TCanvas("c_kppim","c_kppim",900,600);
  c_kppim->cd();
  h_kppim->Draw();
  c_kppim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppim.root",name.Data()), "root");
  c_kppim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppim.eps",name.Data()), "eps");
 
  //4 +++ kppip
  TH1F *h_kppip = new TH1F("h_kppip", ";m_{K^{+}#pi^{+}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kppip", "kppip_mf", "w8*(kppip_uni)");  
  TCanvas *c_kppip = new TCanvas("c_kppip","c_kppip",900,600);
  c_kppip->cd();
  h_kppip->Draw();
  c_kppip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppip.root",name.Data()), "root");
  c_kppip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppip.eps",name.Data()), "eps");

  //5 +++ kpp
  TH1F *h_kpp = new TH1F("h_kpp", ";m_{K^{+}p} [GeV/c^{2}];Counts", 200, 1.4, 3.);//1.4, 3.2
  t->Project("h_kpp", "kpp_mf", "w8*(kpp_uni)");  
  TCanvas *c_kpp = new TCanvas("c_kpp","c_kpp",900,600);
  c_kpp->cd();
  h_kpp->Draw();//"hist"
  c_kpp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpp.root",name.Data()), "root");
  c_kpp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpp.eps",name.Data()), "eps");

  //6 +++ kmpim
  TH1F *h_kmpim = new TH1F("h_kmpim", ";m_{K^{-}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kmpim", "kmpim_mf", "w8*(kmpim_uni)");  
  TCanvas *c_kmpim = new TCanvas("c_kmpim","c_kmpim",900,600);
  c_kmpim->cd();
  h_kmpim->Draw();
  c_kmpim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpim.root",name.Data()), "root");
  c_kmpim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpim.eps",name.Data()), "eps");

  //7 +++ kmpip
  TH1F *h_kmpip = new TH1F("h_kmpip", ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Counts", 200, 0.6, 2);
  t->Project("h_kmpip", "kmpip_mf", "w8*(kmpip_uni)");  
  TCanvas *c_kmpip = new TCanvas("c_kmpip","c_kmpip",900,600);
  c_kmpip->cd();
  h_kmpip->Draw();
  c_kmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpip.root",name.Data()), "root");
  c_kmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpip.eps",name.Data()), "eps");

  //8 +++ kmp
  TH1F *h_kmp = new TH1F("h_kmp", ";m_{K^{-}p} [GeV/c^{2}];Counts", 200, 1.4, 3);
  t->Project("h_kmp", "kmp_mf", "w8*(kmp_uni)");  
  TCanvas *c_kmp = new TCanvas("c_kmp","c_kmp",900,600);
  c_kmp->cd();
  h_kmp->Draw();
  c_kmp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmp.root",name.Data()), "root");
  c_kmp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmp.eps",name.Data()), "eps");

  // 9 +++ pipp
  TH1F *h_pipp = new TH1F("h_pipp", ";m_{#pi^{+}p} [GeV/c^{2}];Counts", 200, 0.9, 2.5);
  t->Project("h_pipp", "pipp_mf", "w8*(pipp_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pipp_mf-1.2)>0.2)");  
  TCanvas *c_pipp = new TCanvas("c_pipp","c_pipp",900,600);
  c_pipp->cd();
  h_pipp->Draw();
  c_pipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pipp.root",name.Data()), "root");
  c_pipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pipp.eps",name.Data()), "eps");

  //10 +++ pimp
  TH1F *h_pimp = new TH1F("h_pimp", ";m_{#pi^{-}p} [GeV/c^{2}];Counts", 200, 0.9, 2.5);
  t->Project("h_pimp", "pimp_mf", "w8*(pimp_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pimp_mf-1.5)>0.2)");  
  TCanvas *c_pimp = new TCanvas("c_pimp","c_pimp",900,600);
  c_pimp->cd();
  h_pimp->Draw();
  c_pimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pimp.root",name.Data()), "root");
  c_pimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pimp.eps",name.Data()), "eps");

  //11 +++ kpkmp
  TH1F *h_kpkmp = new TH1F("h_kpkmp", ";m_{K^{+}K^{-}p} [GeV/c^{2}];Counts", 200, 2, 4);
  t->Project("h_kpkmp", "kpkmp_mf", "w8*(kpkmp_uni)");  
  TCanvas *c_kpkmp = new TCanvas("c_kpkmp","c_kpkmp",900,600);
  c_kpkmp->cd();
  h_kpkmp->Draw();
  c_kpkmp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmp.root",name.Data()), "root");
  c_kpkmp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmp.eps",name.Data()), "eps");

  //12 +++ kpkmpip
  TH1F *h_kpkmpip = new TH1F("h_kpkmpip", ";m_{K^{+}K^{-}#pi^{+}} [GeV/c^{2}];Counts", 200, 1.2, 3.2);
  t->Project("h_kpkmpip", "kpkmpip_mf", "w8*(kpkmpip_uni)");  
  TCanvas *c_kpkmpip = new TCanvas("c_kpkmpip","c_kpkmpip",900,600);
  c_kpkmpip->cd();
  h_kpkmpip->Draw();
  c_kpkmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpip.root",name.Data()), "root");
  c_kpkmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpip.eps",name.Data()), "eps");

  //13 +++ kpkmpim
  TH1F *h_kpkmpim = new TH1F("h_kpkmpim", ";m_{K^{+}K^{-}#pi^{-}} [GeV/c^{2}];Counts", 200, 1.2, 3.2);
  t->Project("h_kpkmpim", "kpkmpim_mf", "w8*(kpkmpim_uni)");  
  TCanvas *c_kpkmpim = new TCanvas("c_kpkmpim","c_kpkmpim",900,600);
  c_kpkmpim->cd();
  h_kpkmpim->Draw();
  c_kpkmpim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpim.root",name.Data()), "root");
  c_kpkmpim->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpim.eps",name.Data()), "eps");

  //14 +++ pippimp
  TH1F *h_pippimp = new TH1F("h_pippimp", ";m_{#pi^{+}#pi^{-}p} [GeV/c^{2}];Counts", 200, 1.3, 3.2);
  t->Project("h_pippimp", "pippimp_mf", "w8*(pippimp_uni)");  
  TCanvas *c_pippimp = new TCanvas("c_pippimp","c_pippimp",900,600);
  c_pippimp->cd();
  h_pippimp->Draw();
  c_pippimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimp.root",name.Data()), "root");
  c_pippimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimp.eps",name.Data()), "eps");

  //15 +++ pippimkp
  TH1F *h_pippimkp = new TH1F("h_pippimkp", ";m_{#pi^{+}#pi^{-}K^{+}} [GeV/c^{2}];Counts", 200, 0.9, 2.7);
  t->Project("h_pippimkp", "pippimkp_mf", "w8*(pippimkp_uni)");  
  TCanvas *c_pippimkp = new TCanvas("c_pippimkp","c_pippimkp",900,600);
  c_pippimkp->cd();
  h_pippimkp->Draw();
  c_pippimkp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimkp.root",name.Data()), "root");
  c_pippimkp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimkp.eps",name.Data()), "eps");

  //16 +++ pippimkm
  TH1F *h_pippimkm = new TH1F("h_pippimkm", ";m_{#pi^{+}#pi^{-}K^{-}} [GeV/c^{2}];Counts", 200, 0.9, 2.7);
  t->Project("h_pippimkm", "pippimkm_mf", "w8*(pippimkm_uni)");  
  TCanvas *c_pippimkm = new TCanvas("c_pippimkm","c_pippimkm",900,600);
  c_pippimkm->cd();
  h_pippimkm->Draw();
  c_pippimkm->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimkm.root",name.Data()), "root");
  c_pippimkm->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_pippimkm.eps",name.Data()), "eps");

  //17 +++ kppimp
  TH1F *h_kppimp = new TH1F("h_kppimp", ";m_{K^{+}#pi^{-}p} [GeV/c^{2}];Counts", 200, 1.6, 3.5);
  t->Project("h_kppimp", "kppimp_mf", "w8*(kppimp_uni)");  
  TCanvas *c_kppimp = new TCanvas("c_kppimp","c_kppimp",900,600);
  c_kppimp->cd();
  h_kppimp->Draw();
  c_kppimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimp.root",name.Data()), "root");
  c_kppimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimp.eps",name.Data()), "eps");

  //18 +++ kmpipp
  TH1F *h_kmpipp = new TH1F("h_kmpipp", ";m_{K^{-}#pi^{+}p} [GeV/c^{2}];Counts", 200, 1.6, 3.5);
  t->Project("h_kmpipp", "kmpipp_mf", "w8*(kmpipp_uni)");  
  TCanvas *c_kmpipp = new TCanvas("c_kmpipp","c_kmpipp",900,600);
  c_kmpipp->cd();
  h_kmpipp->Draw();
  c_kmpipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpipp.root",name.Data()), "root");
  c_kmpipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kmpipp.eps",name.Data()), "eps");

*/
/*  
  // Daltz analysis phipip vs. phipim

  TH2D *h2_phipipvsphipim_tot = new TH2D("h2_phipipvsphipim_tot", ";m^{2}_{#phi#pi^{+}} [GeV/c^{2}];m^{2}_{#phi#pi^{-}} [GeV/c^{2}]", 200, 0.9, 8, 200, 0.9, 8);
  t->Project("h2_phipipvsphipim_tot", "kpkmpip_mf^2:kpkmpim_mf^2", "w8*((kpkmpip_uni || kpkmpim_uni) && kpkm_mf>1.005 && kpkm_mf<1.035)"); 
  TCanvas *c_phipipvsphipim_tot = new TCanvas("c_phipipvsphipim_tot","c_phipipvsphipim_tot",900,600);
  c_phipipvsphipim_tot->cd();
  h2_phipipvsphipim_tot->Draw("colz");
  c_phipipvsphipim_tot->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_tot.root",name.Data()), "root");
  c_phipipvsphipim_tot->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_tot.eps",name.Data()), "eps");

  TH2D *h2_phipipvsphipim_1 = new TH2D("h2_phipipvsphipim_1", "1.6<m_{#phi#pi^{+}#pi^{-}}<2 GeV/c^{2};m^{2}_{#phi#pi^{+}} [GeV/c^{2}];m^{2}_{#phi#pi^{-}} [GeV/c^{2}]", 200, 0.9, 8, 200, 0.9, 8);
  t->Project("h2_phipipvsphipim_1", "kpkmpip_mf^2:kpkmpim_mf^2", "w8*((kpkmpip_uni || kpkmpim_uni) && kpkm_mf>1.005 && kpkm_mf<1.035 && kpkmpippim_mf>1.6 && kpkmpippim_mf<2)"); 
  TCanvas *c_phipipvsphipim_1 = new TCanvas("c_phipipvsphipim_1","c_phipipvsphipim_1",900,600);
  c_phipipvsphipim_1->cd();
  h2_phipipvsphipim_1->Draw("colz");
  c_phipipvsphipim_1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_1.root",name.Data()), "root");
  c_phipipvsphipim_1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_1.eps",name.Data()), "eps"); 

  TH2D *h2_phipipvsphipim_2 = new TH2D("h2_phipipvsphipim_2", "2<m_{#phi#pi^{+}#pi^{-}}<2.5 GeV/c^{2};m^{2}_{#phi#pi^{+}} [GeV/c^{2}];m^{2}_{#phi#pi^{-}} [GeV/c^{2}]", 200, 0.9, 8, 200, 0.9, 8);
  t->Project("h2_phipipvsphipim_2", "kpkmpip_mf^2:kpkmpim_mf^2", "w8*((kpkmpip_uni || kpkmpim_uni) && kpkm_mf>1.005 && kpkm_mf<1.035 && kpkmpippim_mf>2 && kpkmpippim_mf<2.5)"); 
  TCanvas *c_phipipvsphipim_2 = new TCanvas("c_phipipvsphipim_2","c_phipipvsphipim_2",900,600);
  c_phipipvsphipim_2->cd();
  h2_phipipvsphipim_2->Draw("colz");
  c_phipipvsphipim_2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_2.root",name.Data()), "root");
  c_phipipvsphipim_2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_2.eps",name.Data()), "eps");

  TH2D *h2_phipipvsphipim_3 = new TH2D("h2_phipipvsphipim_3", "2.5<m_{#phi#pi^{+}#pi^{-}}<3.2 GeV/c^{2};m^{2}_{#phi#pi^{+}} [GeV/c^{2}];m^{2}_{#phi#pi^{-}} [GeV/c^{2}]", 200, 0.9, 8, 200, 0.9, 8);
  t->Project("h2_phipipvsphipim_3", "kpkmpip_mf^2:kpkmpim_mf^2", "w8*((kpkmpip_uni || kpkmpim_uni) && kpkm_mf>1.005 && kpkm_mf<1.035 && kpkmpippim_mf>2.5 && kpkmpippim_mf<3.2)"); 
  TCanvas *c_phipipvsphipim_3 = new TCanvas("c_phipipvsphipim_3","c_phipipvsphipim_3",900,600);
  c_phipipvsphipim_3->cd();
  h2_phipipvsphipim_3->Draw("colz");
  c_phipipvsphipim_3->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_3.root",name.Data()), "root");
  c_phipipvsphipim_3->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_phipipvsphipim_3.eps",name.Data()), "eps");
*/
/*
  // Nine tiles analysis
  // +++ kppim vs. kmpip
  TH2D *h2_kppimvskmpip = new TH2D("h2_kppimvskmpip", ";m_{K^{+}#pi^{-}} [GeV/c^{2}];m_{K^{-}#pi^{+}} [GeV/c^{2}]", 200, 0.6, 2, 200, 0.6, 2);
  t->Project("h2_kppimvskmpip", "kppim_mf:kmpip_mf", "w8*(kppim_uni || kmpip_uni)"); 
  TCanvas *c_kppimvskmpip = new TCanvas("c_kppimvskmpip","c_kppimvskmpip",900,600);
  c_kppimvskmpip->cd();
  TBox *box1 = new TBox(0.7, 1.0, 0.8, 1.1);
  TBox *box2 = new TBox(0.85, 1.0, 0.95, 1.1);  
  TBox *box3 = new TBox(1.0, 1.0, 1.1, 1.1);
  TBox *box4 = new TBox(0.7, 0.85, 0.8, 0.95);
  TBox *box5 = new TBox(0.85, 0.85, 0.95, 0.95);
  TBox *box6 = new TBox(1.0, 0.85, 1.1, 0.95);
  TBox *box7 = new TBox(0.7, 0.7, 0.8, 0.8);
  TBox *box8 = new TBox(0.85, 0.7, 0.95, 0.8);
  TBox *box9 = new TBox(1.0, 0.7, 1.1, 0.8);
  box1->SetLineColor(kMagenta);
  box2->SetLineColor(kGreen);
  box3->SetLineColor(kMagenta);
  box4->SetLineColor(kGreen);
  box5->SetLineColor(kBlack);
  box6->SetLineColor(kGreen);
  box7->SetLineColor(kMagenta);
  box8->SetLineColor(kGreen);
  box9->SetLineColor(kMagenta);
  box1->SetFillStyle(0);
  box2->SetFillStyle(0);
  box3->SetFillStyle(0);
  box4->SetFillStyle(0);
  box5->SetFillStyle(0);
  box6->SetFillStyle(0);
  box7->SetFillStyle(0);
  box8->SetFillStyle(0);
  box9->SetFillStyle(0);
  box1->SetLineWidth(2);
  box2->SetLineWidth(2);
  box3->SetLineWidth(2);
  box4->SetLineWidth(2);
  box5->SetLineWidth(2);
  box6->SetLineWidth(2);
  box7->SetLineWidth(2);
  box8->SetLineWidth(2);
  box9->SetLineWidth(2);
  h2_kppimvskmpip->Draw("colz");
  box1->Draw("same");
  box2->Draw("same");
  box3->Draw("same");
  box4->Draw("same");
  box5->Draw("same");
  box6->Draw("same");
  box7->Draw("same");
  box8->Draw("same");
  box9->Draw("same");
  c_kppimvskmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimvskmpip.root",name.Data()), "root");
  c_kppimvskmpip->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kppimvskmpip.eps",name.Data()), "eps");
 
  // +++ K*K*
  TCanvas *c_2kstar = new TCanvas("c_2kstar","c_2kstar",1200,900);
  c_2kstar->Divide(3,3);
  TH1F *h_2kstar_1 = new TH1F("h_2kstar_1", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_1", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.7 && kppim_mf<0.8 && kmpip_mf>1.0 && kmpip_mf<1.1)");
  c_2kstar->cd(1);  
  h_2kstar_1->Draw("colz");
  TH1F *h_2kstar_2 = new TH1F("h_2kstar_2", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_2", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.85 && kppim_mf<0.95 && kmpip_mf>1.0 && kmpip_mf<1.1)"); 
  c_2kstar->cd(2);
  h_2kstar_2->Draw("colz");
  TH1F *h_2kstar_3 = new TH1F("h_2kstar_3", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_3", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>1.0 && kppim_mf<1.1 && kmpip_mf>1.0 && kmpip_mf<1.1)"); 
  c_2kstar->cd(3);
  h_2kstar_3->Draw("colz");
  TH1F *h_2kstar_4 = new TH1F("h_2kstar_4", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_4", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.7 && kppim_mf<0.8 && kmpip_mf>0.85 && kmpip_mf<0.95)");
  c_2kstar->cd(4);
  h_2kstar_4->Draw("colz");
  TH1F *h_2kstar_5 = new TH1F("h_2kstar_5", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_5", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.85 && kppim_mf<0.95 && kmpip_mf>0.85 && kmpip_mf<0.95)"); 
  c_2kstar->cd(5);
  h_2kstar_5->Draw("colz");
  TH1F *h_2kstar_6 = new TH1F("h_2kstar_6", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_6", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>1.0 && kppim_mf<1.1 && kmpip_mf>0.85 && kmpip_mf<0.95)");
  c_2kstar->cd(6);
  h_2kstar_6->Draw("colz");
  TH1F *h_2kstar_7 = new TH1F("h_2kstar_7", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_7", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.7 && kppim_mf<0.8 && kmpip_mf>0.7 && kmpip_mf<0.8)"); 
  c_2kstar->cd(7);
  h_2kstar_7->Draw("colz");
  TH1F *h_2kstar_8 = new TH1F("h_2kstar_8", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_8", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>0.85 && kppim_mf<0.95 && kmpip_mf>0.7 && kmpip_mf<0.8)");
  c_2kstar->cd(8);
  h_2kstar_8->Draw("colz");
  TH1F *h_2kstar_9 = new TH1F("h_2kstar_9", ";m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  t->Project("h_2kstar_9", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kppim_mf>1.0 && kppim_mf<1.1 && kmpip_mf>0.7 && kmpip_mf<0.8)");
  c_2kstar->cd(9);
  h_2kstar_9->Draw("colz");
  c_2kstar->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_2kstar.root",name.Data()), "root");
  c_2kstar->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_2kstar.eps",name.Data()), "eps");

  TCanvas *c_2kstar_sub = new TCanvas("c_2kstar_sub","c_2kstar_sub",900,600);
  TH1F *h_2kstar_sub = new TH1F("h_2kstar_sub", "side band subtracted;m_{K^{*} #bar{K}^{*}} [GeV/c^{2}];Counts", 100, 1.7, 3.0);
  // TList *list = new TList;
  // list->Add(h_2kstar_5, 1.);
  // list->Add(h_2kstar_2, -0.5);
  // list->Add(h_2kstar_4, -0.5);
  // list->Add(h_2kstar_6, -0.5);
  // list->Add(h_2kstar_8, -0.5);
  // list->Add(h_2kstar_1, 0.25);
  // list->Add(h_2kstar_3, 0.25);
  // list->Add(h_2kstar_7, 0.25);
  // list->Add(h_2kstar_9, 0.25);
  // c_2kstar_sub->cd();
  // h_2kstar_sub->Reset();
  // h_2kstar_sub->Merge(list);
  // h_2kstar_sub->Draw();
  // h_2kstar_sub->Add(h_2kstar_5, h_2kstar_2,h_2kstar_4,h_2kstar_6,h_2kstar_8, h_2kstar_1,h_2kstar_3,h_2kstar_7,h_2kstar_9, 1., -0.5,-0.5,-0.5,-0.5, 0.25,0.25,0.25,0.25);
  TH1F *h_2kstar_sub1 = new TH1F("h_2kstar_sub1", "h_2kstar_sub1", 100, 1.7, 3.0);
  h_2kstar_sub1->Add(h_2kstar_5, h_2kstar_2, 1., -0.5);
  TH1F *h_2kstar_sub2 = new TH1F("h_2kstar_sub2", "h_2kstar_sub2", 100, 1.7, 3.0);
  h_2kstar_sub2->Add(h_2kstar_sub1, h_2kstar_4, 1., -0.5);
  TH1F *h_2kstar_sub3 = new TH1F("h_2kstar_sub3", "h_2kstar_sub3", 100, 1.7, 3.0);
  h_2kstar_sub3->Add(h_2kstar_sub2, h_2kstar_6, 1., -0.5);
  TH1F *h_2kstar_sub4 = new TH1F("h_2kstar_sub4", "h_2kstar_sub4", 100, 1.7, 3.0);
  h_2kstar_sub4->Add(h_2kstar_sub3, h_2kstar_8, 1., -0.5);
  TH1F *h_2kstar_sub5 = new TH1F("h_2kstar_sub5", "h_2kstar_sub5", 100, 1.7, 3.0);
  h_2kstar_sub5->Add(h_2kstar_sub4, h_2kstar_1, 1., 0.25);
  TH1F *h_2kstar_sub6 = new TH1F("h_2kstar_sub6", "h_2kstar_sub6", 100, 1.7, 3.0);
  h_2kstar_sub6->Add(h_2kstar_sub5, h_2kstar_3, 1., 0.25);
  TH1F *h_2kstar_sub7 = new TH1F("h_2kstar_sub7", "h_2kstar_sub7", 100, 1.7, 3.0);
  h_2kstar_sub7->Add(h_2kstar_sub6, h_2kstar_7, 1., 0.25);

  h_2kstar_sub->Add(h_2kstar_sub7, h_2kstar_9, 1., 0.25);
  c_2kstar_sub->cd();
  h_2kstar_sub->Draw();
  c_2kstar_sub->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_2kstar_sub.root",name.Data()), "root");
  c_2kstar_sub->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_2kstar_sub.eps",name.Data()), "eps");  


  // +++ kpkmpippim vs. pipp
  TH2D *h2_kpkmpippimvspipp = new TH2D("h2_kpkmpippimvspipp", Form("%s;m_{#pi^{+}p} [GeV/c^{2}];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]",name.Data()), 100, 1.0, 3, 100, 1.6, 3.2);
  t->Project("h2_kpkmpippimvspipp", "kpkmpippim_mf:pipp_mf", "w8*((kpkmpippim_uni || pipp_uni) && kpkm_mf>1.005 && kpkm_mf<1.035)"); 
  TCanvas *c_kpkmpippimvspipp = new TCanvas("c_kpkmpippimvspipp","c_kpkmpippimvspipp",900,600);
  c_kpkmpippimvspipp->cd();
  h2_kpkmpippimvspipp->Draw("colz");
  c_kpkmpippimvspipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspipp.root",name.Data()), "root");
  c_kpkmpippimvspipp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspipp.eps",name.Data()), "eps");

  // // After VanHove cut
  // TH2D *h2_kpkmpippimvspipp_vhcut = new TH2D("h2_kpkmpippimvspipp_vhcut", Form("%s;m_{#pi^{+}p} [GeV/c^{2}];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]",name.Data()), 100, 1.0, 3, 100, 1.6, 3.2);
  // t->Project("h2_kpkmpippimvspipp_vhcut", "kpkmpippim_mf:pipp_mf", "w8*((kpkmpippim_uni || pipp_uni) && kpkm_mf>1.005 && kpkm_mf<1.035 && vhtheta>0.4 && vhtheta<1.5 && vhphi>1.55 && vhphi<2.2)"); 
  // TCanvas *c_kpkmpippimvspipp_vhcut = new TCanvas("c_kpkmpippimvspipp_vhcut","c_kpkmpippimvspipp_vhcut",900,600);
  // c_kpkmpippimvspipp_vhcut->cd();
  // h2_kpkmpippimvspipp_vhcut->Draw("colz");
  // c_kpkmpippimvspipp_vhcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspipp_vhcut.root",name.Data()), "root");
  // c_kpkmpippimvspipp_vhcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspipp_vhcut.eps",name.Data()), "eps");

  // +++ kpkmpippim vs. pimp
  TH2D *h2_kpkmpippimvspimp = new TH2D("h2_kpkmpippimvspimp", Form("%s;m_{#pi^{-}p} [GeV/c^{2}];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]",name.Data()), 100, 1.0, 3, 100, 1.6, 3.2);
  t->Project("h2_kpkmpippimvspimp", "kpkmpippim_mf:pimp_mf", "w8*((kpkmpippim_uni || pimp_uni) && kpkm_mf>1.005 && kpkm_mf<1.035)"); 
  TCanvas *c_kpkmpippimvspimp = new TCanvas("c_kpkmpippimvspimp","c_kpkmpippimvspimp",900,600);
  c_kpkmpippimvspimp->cd();
  h2_kpkmpippimvspimp->Draw("colz");
  c_kpkmpippimvspimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspimp.root",name.Data()), "root");
  c_kpkmpippimvspimp->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_y2175/c%s_kpkmpippimvspimp.eps",name.Data()), "eps");
*/
}