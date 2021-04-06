#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFormula.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include <fstream>
#include "TStyle.h"
#include <TMultiGraph.h>
#include <iostream>
#include "TSystem.h"
#include "scanphi.C"
#include <vector>
#include "pdgavg.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooStats/NumberCountingUtils.h"
#endif
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAbsArg.h"
#include "RooCmdArg.h"
#include "RooCurve.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooBreitWigner.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooExponential.h"
#include "RooStats/RooStatsUtils.h"
using namespace RooFit;
using namespace RooStats;

void xsec_phifo(TString name, int n2k=100, int n2pi=50, int ne=1, int nt=1)//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", name.Data())); // _chi100
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_phifo_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/thrown_phifo_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_11366_11555.root");
  if(name == "17") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_30274_31057.root");
  if(name == "18") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_40856_42577.root");
  if(name == "18l") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_50677_51768.root");

  TFile *outputfig = new TFile("output/fig_phifo/xsec_phifo.root", "UPDATE");
  
  // ofstream ofs_phifo("phifo_tab.txt", ofstream::out);

  ofstream test;
  test.open ("test.txt");

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.99, mkk_max = 1.2; //mkk_min = 0.98, mkk_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;
  double m2pi_min2 = 0.83, m2pi_max2 = 1.16; // m2pi_min2 = 0.83, m2pi_max2 = 1.16, 1.15, 1.17

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // ##############################################  fo ALL ##############################################

  FILE *table_xsec_phifo = fopen(Form("table_xsec_phifo_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{0}$ & $N_{MC}$ & $N_{Data}$ & $\\epsilon$ [ \\% ] & $\\sigma \\times BR_{f_{0}\\rightarrow\\pi^{+}\\pi^{-}}$ [nb] & $<\\sigma>$ [nb] \\\\ \n \\hline\n");

  FILE *table_xsec_phifo_sys = fopen(Form("table_xsec_phifo_sys_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phifo_sys,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{makecell} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\small\n \\centering\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n \\thead{Bkg deg\\\\(3,4,5)} & \\thead{Fit range\\\\(2; 3, 2.9, 3.1)} & \\thead{bining\\\\(90, 100, 110)} & \\thead{Accidentals\\\\(2, 4, 8)} & \\thead{$\\Delta T_{RF-TOF}$ proton\\\\($\\pm$0.2, $\\pm$0.3, $\\pm$0.4)} & \\thead{$\\chi^{2}$ KinFit\\\\(45, 55, 65)} & \\thead{$MM^{2}$\\\\($\\pm$0.025, $\\pm$0.035, $\\pm$0.045)} \\\\ \n \\hline\n");  

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  // h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
   TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  // h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");
/*
  // *********************** Phi(1020) MC *********************************
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", "; m_{K^{+}K^{-}} (GeV/c^{2}); Counts", 200, 1.005, 1.035);
  tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");//+cutlist+" && kin_chisq<25)"
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");

  // w.factory("BreitWigner::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004])");
  // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.001,0.1])"); //m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.001,0.01]
  w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-10000,10000], c1_PhiMass_mc[-10000,10000], c2_PhiMass_mc[-10000,10000]})");
  w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,1000000]*sig_PhiMass_mc, nbkg_PhiMass_mc[0,1000000]*bkg_PhiMass_mc)");//, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); //nsig[0,100000000]*sig2,
  w.var("m_PhiMass_mc")->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
  RooDataHist dh_PhiMass_mc("dh_PhiMass_mc", "dh_PhiMass_mc", *w.var("m_PhiMass_mc"), Import(*h_PhiMass_postcuts));
  RooPlot *fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title(" "));
  // fr_PhiMass_mc->SetTitleOffset(0.90, "X");
  // fr_PhiMass_mc->SetTitleSize(0.06, "XYZ");
  // fr_PhiMass_mc->SetLabelSize(0.06, "xyz");
  w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc);

  // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
  dh_PhiMass_mc.plotOn(fr_PhiMass_mc, RooFit::Name("ldh_PhiMass_mc"));
  w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig_PhiMass_mc")), LineColor(kRed), RooFit::Name("lsig_PhiMass_mc"));
  w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("bkg_PhiMass_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass_mc"));
  w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, RooFit::Name("lmodel_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->paramOn(fr_PhiMass_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_PhiMass_mc"), *w.var("nbkg_PhiMass_mc")))); //,*w.var("mean_PhiMass_mc"),*w.var("width_PhiMass_mc"),*w.var("sigma_PhiMass_mc"))));
  fr_PhiMass_mc->Draw();

  TLegend *l_phi_mc = new TLegend(0.2, 0.65, 0.4, 0.85);
  l_phi_mc->SetFillColor(kWhite);
  l_phi_mc->SetLineColor(kWhite);
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lmodel_PhiMass_mc"), "total", "l");
  l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lsig_PhiMass_mc"), "Voigtian", "l");
  l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lbkg_PhiMass_mc"), "pol 2nd", "l");
  l_phi_mc->Draw();

  TLatex lat_PhiMass_mc;
  lat_PhiMass_mc.SetTextSize(0.05);
  lat_PhiMass_mc.SetTextAlign(13); //align at top
  lat_PhiMass_mc.SetNDC();
  lat_PhiMass_mc.SetTextColor(kBlue);
  lat_PhiMass_mc.DrawLatex(0.62, 0.86, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_PhiMass_mc")->getVal(), w.var("nsig_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.80, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_PhiMass_mc")->getVal(), w.var("nbkg_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.74, Form("#mu = %0.3f#pm%0.3f", w.var("mean_PhiMass_mc")->getVal(), w.var("mean_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.68, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_PhiMass_mc")->getVal(), w.var("width_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.62, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_PhiMass_mc")->getVal(), w.var("sigma_PhiMass_mc")->getError()));
*/
  // *********************** f0(980) MC *********************************
  TCanvas *c_foMass_postcuts = new TCanvas("c_foMass_postcuts", "c_foMass_postcuts", 1000, 600);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 200, m2pi_min2-0.025, m2pi_max2+0.05); // 0.85, 1.05
  tmc->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)"); 
  // && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(pippim_uni && abs(rf_dt)<1.5*4.008, w4*(pippim_uni && abs(rf_dt)<2.5*4.008, w8*(pippim_uni
  c_foMass_postcuts->cd();
  h_foMass_postcuts->Draw("e");
  
  TF1 *fsb_fo_mc = new TF1("fsb_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol1(3)", m2pi_min2, m2pi_max2); // pol1(3)
  fsb_fo_mc->SetLineColor(2);
  fsb_fo_mc->SetParameters(h_foMass_postcuts->GetMaximum(), 0.98, 0.05, 1,1); // 1,1
  // fsb_fo_mc->SetParLimits(0, 0, h_foMass_postcuts->GetMaximum());
  fsb_fo_mc->SetParLimits(1, 0.95, 0.99);  // [0.95, 0.1]
  fsb_fo_mc->SetParLimits(2, 0.03, 0.09); //[0.001, 0.1]
  // fsb_fo_mc->FixParameter(1, fsb_fo_mc->GetParameter(1)); //0.988
  // fsb_fo_mc->FixParameter(2, fsb_fo_mc->GetParameter(2)); // 0.059, 0.070, 0.063, 0.062

  // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
  TF1 *fs_fo_mc = new TF1("fs_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2])", m2pi_min2-0.025, m2pi_max2+0.05);
  fs_fo_mc->SetLineColor(4);

  TF1 *fb_fo_mc = new TF1("fb_fo_mc", "pol1(3)", m2pi_min2-0.025, m2pi_max2+0.05); //pol1(3)
  fb_fo_mc->SetLineColor(28);
  fb_fo_mc->SetLineStyle(2);

  h_foMass_postcuts->Fit("fsb_fo_mc", "", "", m2pi_min2-0.025, m2pi_max2+0.05);
  double par_fo_mc[fsb_fo_mc->GetNpar()]; //6
  fsb_fo_mc->GetParameters(&par_fo_mc[0]);
  fs_fo_mc->SetParameters(&par_fo_mc[0]);
  fb_fo_mc->SetParameters(&par_fo_mc[3]); //3

  double N_fo_mc = fs_fo_mc->Integral(m2pi_min2, m2pi_max2) * n2pi / (m2pi_max - m2pi_min);
  double dN_fo_mc = N_fo_mc * fsb_fo_mc->GetParError(0) / fsb_fo_mc->GetParameter(0);

  fs_fo_mc->Draw("same");
  fb_fo_mc->Draw("same");

  TLatex lat_fo_mc;
  lat_fo_mc.SetTextSize(0.05);
  lat_fo_mc.SetTextAlign(13); //align at top
  lat_fo_mc.SetNDC();
  lat_fo_mc.SetTextColor(kBlue);
  lat_fo_mc.DrawLatex(0.2, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_fo_mc->GetChisquare() / fsb_fo_mc->GetNDF()));
  lat_fo_mc.DrawLatex(0.2, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_fo_mc, dN_fo_mc));
  lat_fo_mc.DrawLatex(0.2, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_fo_mc->GetParameter(1), fsb_fo_mc->GetParError(1)));
  lat_fo_mc.DrawLatex(0.2, 0.68, Form("#Gamma = %0.3f#pm%0.3f", fsb_fo_mc->GetParameter(2), fsb_fo_mc->GetParError(2)));

  // ======================================== fo vs. Phi(1020) ===============================================
  // root -l 'xsec_phifo.C+(100,50,1,1)'

  TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 900, 600); //1500, 800
  TGraphErrors *gphifo = scanphi("fo", "phifo", Form("%s", name.Data()), n2k, n2pi);
  cgphifo->cd();
  gphifo->Draw("AP");

  TF1 *fsb_fo_data = new TF1("fsb_fo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", m2pi_min2, m2pi_max2); 
  // [0]*TMath::BreitWigner(x,[1],[2]) + expo(3)+pol0(5); pol2(3)
  fsb_fo_data->SetParameters(1000, 0.970, 0.06, 6.7, 2.8, 1,1,1); // 80, 0.970, 0.06, 6.7, 2.8, 1,1,1
  fsb_fo_data->SetLineColor(2);
  fsb_fo_data->SetParLimits(1, 0.95, 0.99); // 1, 0.95, 0.99
  fsb_fo_data->SetParLimits(2, 0.03, 0.09); // 2, 0.03, 0.09
  // fsb_fo_data->SetParLimits(2, 0.04, 0.1);
  // fsb_fo_data->SetParLimits(1, 0.94, 0.98); // 1, 0.96, 0.99
  // fsb_fo_data->SetParLimits(2, 0.03, 0.085); // 0.05, 0.075, 0.04, 0.08

  if (name == "16")
  {
    fsb_fo_data->SetParLimits(2, 0.06, 0.07); // 2, 0.060, 0.070

    // fsb_fo_data->FixParameter(1, 0.977-0.01);
    // fsb_fo_data->FixParameter(2, 0.070+0.002);
  }
  // if (name == "17")
  // {
  //   // fsb_fo_data->FixParameter(1, 0.968-0.01);
  //   // fsb_fo_data->FixParameter(2, 0.067+0.002);
  // }
  // if (name == "18")
  // {
  //   // fsb_fo_data->FixParameter(1, 0.973-0.01);
  //   // fsb_fo_data->FixParameter(2, 0.059+0.002);
  // }
  // if (name == "18l")
  // {
  //   // fsb_fo_data->FixParameter(1, 0.968-0.01);
  //   // fsb_fo_data->FixParameter(2, 0.058+0.002);
  // }

  TF1 *fs_fo_data = new TF1("fs_fo_data", "[0]*TMath::BreitWigner(x,[1],[2])", m2pi_min2, m2pi_max2);
  fs_fo_data->SetLineColor(4);

  TF1 *fb_fo_data = new TF1("fb_fo_data", "pol2(3)", m2pi_min2, m2pi_max2); //pol2(3)
  fb_fo_data->SetLineColor(28);
  fb_fo_data->SetLineStyle(2);

  gphifo->Fit("fsb_fo_data", "R", "", m2pi_min2, m2pi_max2);
  double par_fo_data[fsb_fo_data->GetNpar()]; //fsb_fo_data->GetNpar()
  fsb_fo_data->GetParameters(&par_fo_data[0]);
  fs_fo_data->SetParameters(&par_fo_data[0]);
  fb_fo_data->SetParameters(&par_fo_data[3]); //3

  double N_fo_data = fs_fo_data->Integral(m2pi_min2, m2pi_max2) * n2pi / (m2pi_max - m2pi_min);
  double dN_fo_data = N_fo_data * fsb_fo_data->GetParError(0) / fsb_fo_data->GetParameter(0);

  fs_fo_data->Draw("same");
  fb_fo_data->Draw("same");

  TLatex lat_phifo;
  lat_phifo.SetTextSize(0.05);
  lat_phifo.SetTextAlign(13); //align at top
  lat_phifo.SetNDC();
  lat_phifo.SetTextColor(kBlue);
  lat_phifo.DrawLatex(0.6, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_fo_data->GetChisquare() / fsb_fo_data->GetNDF()));
  lat_phifo.DrawLatex(0.6, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_fo_data, dN_fo_data));
  lat_phifo.DrawLatex(0.6, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(1), fsb_fo_data->GetParError(1)));
  lat_phifo.DrawLatex(0.6, 0.68, Form("#Gamma = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(2), fsb_fo_data->GetParError(2)));
  // lat_phifo.DrawLatex(0.6, 0.62, Form("#sigma = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(3), fsb_fo_data->GetParError(3)));

  // +++++++++++++++++++ Systematic Errors
  // if(name == "16")
  double arr_sig_mean_16[3] = {28.95, 27.38, 28.33};  //mean: mu,-0.01, +0.01 -> 1% of 0.990.  
  double arr_sig_width_16[3] = {28.95, 27.88, 30.05}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_16[3] = {28.95, 28.46, 27.85};  //pol:2,3,4
  double arr_fit_range_16[3] = {28.95, 26.71, 28.33}; //fit:[0.83, 1.16], [0.83, 1.15], [0.83, 1.17] 
  double arr_binning_16[3] = {28.95, 28.60, 29.17};   //binning: 100,90,110
  double arr_beambun_16[3] = {28.95, 26.66, 27.90};   //bunches: w8, w2, w4
  double arr_p_dttof_16[3] = {29.94, 31.03, 30.40};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_16[3] = {29.94, 28.58, 29.95}; // kin_chisq<55, 45, 65
  double arr_mm2_16[3] = {29.94, 28.95, 28.90};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "17")
  double arr_sig_mean_17[3] = {10.36, 12.58, 11.05};   //mean: mu,-0.01, +0.01 -> 1% of 0.990.
  double arr_sig_width_17[3] = {10.36, 10.10, 10.69}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_17[3] = {10.36, 8.47, 9.52};   //pol:2,3,4
  double arr_fit_range_17[3] = {10.36, 12.28, 10.14}; //fit:[0.83, 1.16], [0.83, 1.15], [0.83, 1.17]
  double arr_binning_17[3] = {10.36, 7.91, 11.75};    //binning: 100,90,110
  double arr_beambun_17[3] = {10.36, 8.61, 9.86};    //bunches: w8, w2, w4
  double arr_p_dttof_17[3] = {9.51, 11.22, 14.84};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_17[3] = {9.51, 9.74, 10.01}; // kin_chisq<55, 45, 65
  double arr_mm2_17[3] = {9.51, 9.84, 12.11};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18")
  double arr_sig_mean_18[3] = {13.49, 15.45, 12.97};  //mean: mu,-0.01, +0.01 -> 1% of 0.990. 12.29
  double arr_sig_width_18[3] = {13.49, 13.02, 13.88}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_18[3] = {13.49, 12.52, 11.57};  //pol:2,3,4
  double arr_fit_range_18[3] = {13.49, 14.06, 13.19}; //fit:[0.83, 1.16], [0.83, 1.15], [0.83, 1.17]
  double arr_binning_18[3] = {13.49, 15.09, 11.56};   //binning: 100,90,110
  double arr_beambun_18[3] = {13.49, 16.29, 14.23};   //bunches: w8, w2, w4
  double arr_p_dttof_18[3] = {12.99, 15.95, 17.06};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18[3] = {12.99, 13.97, 19.97}; // kin_chisq<55, 45, 65
  double arr_mm2_18[3] = {12.99, 11.48, 15.36};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18l")
  double arr_sig_mean_18l[3] = {10.28, 10.60, 9.49};  //mean: mu,-0.01, +0.01 -> 1% of 0.990. 9.51
  double arr_sig_width_18l[3] = {10.28, 10.02, 10.66}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_18l[3] = {10.28, 9.47, 8.65};  //pol:2,3,4
  double arr_fit_range_18l[3] = {10.28, 12.49, 10.05}; //fit:[0.83, 1.16], [0.83, 1.15], [0.83, 1.17]
  double arr_binning_18l[3] = {10.28, 8.62, 9.74};   //binning: 100,90,110
  double arr_beambun_18l[3] = {10.28, 9.85, 8.68};   //bunches: w8, w2, w4
  double arr_p_dttof_18l[3] = {9.99, 9.53, 10.26};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18l[3] = {9.99, 10.24, 8.71}; // kin_chisq<55, 45, 65
  double arr_mm2_18l[3] = {9.99, 9.99, 9.94};       // abs(mm2)<0.035, 0.025, 0.045

  double arr_sig_mean[3], arr_sig_width[3], arr_bkg_poln[3], arr_fit_range[3], arr_binning[3], arr_beambun[3], arr_p_dttof[3], arr_kin_chisq[3], arr_mm2[3];

  if (name == "16")
  {
    memcpy(arr_sig_mean, arr_sig_mean_16, sizeof(arr_sig_mean));
    memcpy(arr_sig_width, arr_sig_width_16, sizeof(arr_sig_width));
    memcpy(arr_bkg_poln, arr_bkg_poln_16, sizeof(arr_bkg_poln));
    memcpy(arr_fit_range, arr_fit_range_16, sizeof(arr_fit_range));
    memcpy(arr_binning, arr_binning_16, sizeof(arr_binning));
    memcpy(arr_beambun, arr_beambun_16, sizeof(arr_beambun));
    memcpy(arr_p_dttof, arr_p_dttof_16, sizeof(arr_p_dttof));
    memcpy(arr_kin_chisq, arr_kin_chisq_16, sizeof(arr_kin_chisq));
    memcpy(arr_mm2, arr_mm2_16, sizeof(arr_mm2));
  }
  if (name == "17")
  {
    memcpy(arr_sig_mean, arr_sig_mean_17, sizeof(arr_sig_mean));
    memcpy(arr_sig_width, arr_sig_width_17, sizeof(arr_sig_width));
    memcpy(arr_bkg_poln, arr_bkg_poln_17, sizeof(arr_bkg_poln));
    memcpy(arr_fit_range, arr_fit_range_17, sizeof(arr_fit_range));
    memcpy(arr_binning, arr_binning_17, sizeof(arr_binning));
    memcpy(arr_beambun, arr_beambun_17, sizeof(arr_beambun));
    memcpy(arr_p_dttof, arr_p_dttof_17, sizeof(arr_p_dttof));
    memcpy(arr_kin_chisq, arr_kin_chisq_17, sizeof(arr_kin_chisq));
    memcpy(arr_mm2, arr_mm2_17, sizeof(arr_mm2));
  }
  if (name == "18")
  {
    memcpy(arr_sig_mean, arr_sig_mean_18, sizeof(arr_sig_mean));
    memcpy(arr_sig_width, arr_sig_width_18, sizeof(arr_sig_width));
    memcpy(arr_bkg_poln, arr_bkg_poln_18, sizeof(arr_bkg_poln));
    memcpy(arr_fit_range, arr_fit_range_18, sizeof(arr_fit_range));
    memcpy(arr_binning, arr_binning_18, sizeof(arr_binning));
    memcpy(arr_beambun, arr_beambun_18, sizeof(arr_beambun));
    memcpy(arr_p_dttof, arr_p_dttof_18, sizeof(arr_p_dttof));
    memcpy(arr_kin_chisq, arr_kin_chisq_18, sizeof(arr_kin_chisq));
    memcpy(arr_mm2, arr_mm2_18, sizeof(arr_mm2));
  }
  if (name == "18l")
  {
    memcpy(arr_sig_mean, arr_sig_mean_18l, sizeof(arr_sig_mean));
    memcpy(arr_sig_width, arr_sig_width_18l, sizeof(arr_sig_width));
    memcpy(arr_bkg_poln, arr_bkg_poln_18l, sizeof(arr_bkg_poln));
    memcpy(arr_fit_range, arr_fit_range_18l, sizeof(arr_fit_range));
    memcpy(arr_binning, arr_binning_18l, sizeof(arr_binning));
    memcpy(arr_beambun, arr_beambun_18l, sizeof(arr_beambun));
    memcpy(arr_p_dttof, arr_p_dttof_18l, sizeof(arr_p_dttof));
    memcpy(arr_kin_chisq, arr_kin_chisq_18l, sizeof(arr_kin_chisq));
    memcpy(arr_mm2, arr_mm2_18l, sizeof(arr_mm2));
  }

  // double err_sig_mean = TMath::RMS(3, arr_sig_mean) / arr_sig_mean[0];
  // double err_sig_width = TMath::RMS(3, arr_sig_width) / arr_sig_width[0];
  double err_bkg_poln = TMath::RMS(3, arr_bkg_poln) / arr_bkg_poln[0];
  double err_fit_range = TMath::RMS(3, arr_fit_range) / arr_fit_range[0];
  double err_binning = TMath::RMS(3, arr_binning) / arr_binning[0];
  double err_beambun = TMath::RMS(3, arr_beambun) / arr_beambun[0];
  double err_p_dttof = TMath::RMS(3, arr_p_dttof) / arr_p_dttof[0];
  double err_kin_chisq = TMath::RMS(3, arr_kin_chisq) / arr_kin_chisq[0];
  double err_mm2 = TMath::RMS(3, arr_mm2) / arr_mm2[0];

  double err_sys = TMath::Sqrt(err_bkg_poln * err_bkg_poln + err_fit_range * err_fit_range + err_binning * err_binning + err_beambun * err_beambun + err_p_dttof * err_p_dttof + err_kin_chisq * err_kin_chisq + err_mm2 * err_mm2);

  // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
  double eff_phifo = N_fo_mc / h_beame_tru->Integral(0, 100); // Efficiency = N_observed/N_generated
  double deff_phifo = eff_phifo * (dN_fo_data / N_fo_mc);

  // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
  double lumi_phifo = h_tagged_flux->Integral(0, 100); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
  double T = 1.273;                                    //Target thickness [b^{-1}]
  double br_phi = 0.492;
  double dbr_phi = 0.005;
  // double br_fo = ?;
  // double dbr_fo = ?;
  double xsec_phifo = 1e9 * N_fo_data / (eff_phifo * lumi_phifo * T * br_phi);
  double dxsec_phifo_stat = xsec_phifo * (dN_fo_data / N_fo_data);
  double dxsec_phifo_sys = xsec_phifo * TMath::Sqrt((err_sys*err_sys) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
  // double dxsec_phifo_sys = xsec_phifo * TMath::Sqrt((deff_phifo / eff_phifo) * (deff_phifo / eff_phifo) + (dbr_phi / br_phi) * (dbr_phi / br_phi)); //(err_sys / N_fo_data) * (err_sys / N_fo_data) +
  double dxsec_phifo_tot = TMath::Sqrt(dxsec_phifo_stat*dxsec_phifo_stat + dxsec_phifo_sys*dxsec_phifo_sys);

  // +++++++++++++++++++ PDG cross section average
  vector<double> xsec(3), dxsec(3);
  pair<double, double> xsec_pair;

  xsec = {10.36, 13.49, 10.28};
  dxsec = {6.93, 6.32, 3.00}; // 6.00, 4.55, 2.47
  xsec_pair = pdgavg(xsec, dxsec);
  double xsec_avg = xsec_pair.first;
  double dxsec_avg = xsec_pair.second;

  test << std::setprecision(2) << std::fixed;
  test << xsec_phifo << " | " << dxsec_phifo_tot << endl;

  fprintf(table_xsec_phifo, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f \\\\ \n", h_tagged_flux->GetXaxis()->GetXmin(), h_tagged_flux->GetXaxis()->GetXmax(), h_beame_tru->Integral(0, 100), N_fo_mc, dN_fo_mc, N_fo_data, dN_fo_data, eff_phifo * 100, deff_phifo * 100, xsec_phifo, dxsec_phifo_stat, dxsec_phifo_sys, xsec_avg, dxsec_avg);

  fprintf(table_xsec_phifo_sys, "%0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", err_bkg_poln * 100, err_fit_range * 100, err_binning * 100, err_beambun * 100, err_p_dttof * 100, err_kin_chisq * 100, err_mm2 * 100);

  // table_phi << std::setprecision(2) << std::fixed;
  // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_fo_mc<<"$\\pm$"<<dN_fo_mc<<"&"<<N_fo_data<<"$\\pm$"<<dN_fo_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phifo<<"$\\pm$"<<dxsec_phifo_stat<<"$\\pm$"<<dxsec_phifo_sys<<" \\\\"<< endl;

  // c_PhiMass_postcuts->Print(Form("output/fig_phifo/cmc_PhiMass_postcuts_fitted_%s.root", name.Data()), "root");
  // c_PhiMass_postcuts->Print(Form("output/fig_phifo/cmc_PhiMass_postcuts_fitted_%s.eps", name.Data()), "eps");
  // c_PhiMass_postcuts->Print(Form("output/fig_phifo/cmc_PhiMass_postcuts_fitted_%s.png", name.Data()), "png");
  c_foMass_postcuts->Print(Form("output/fig_phifo/cmc_foMass_postcuts_fitted_%s.root", name.Data()), "root");
  c_foMass_postcuts->Print(Form("output/fig_phifo/cmc_foMass_postcuts_fitted_%s.eps", name.Data()), "eps");
  c_foMass_postcuts->Print(Form("output/fig_phifo/cmc_foMass_postcuts_fitted_%s.png", name.Data()), "png");
  // cphifo1->Print(Form("output/fig_phifo/c1_phifo_%s.root", name.Data()), "root");
  // cphifo1->Print(Form("output/fig_phifo/c1_phifo_%s.eps", name.Data()), "eps");
  // cphifo1->Print(Form("output/fig_phifo/c1_phifo_%s.png", name.Data()), "png");
  // cphifo2->Print(Form("output/fig_phifo/c2_phifo_%s.root", name.Data()), "root");
  // cphifo2->Print(Form("output/fig_phifo/c2_phifo_%s.eps", name.Data()), "eps");
  // cphifo2->Print(Form("output/fig_phifo/c2_phifo_%s.png", name.Data()), "png");
  // cphifo->Print(Form("output/fig_phifo/c_phifo_%s.root", name.Data()), "root");
  // cphifo->Print(Form("output/fig_phifo/c_phifo_%s.eps", name.Data()), "eps");
  // cphifo->Print(Form("output/fig_phifo/c_phifo_%s.png", name.Data()), "png");
  cgphifo->Print(Form("output/fig_phifo/c_gphifo_%s.root", name.Data()), "root");
  cgphifo->Print(Form("output/fig_phifo/c_gphifo_%s.eps", name.Data()), "eps");
  cgphifo->Print(Form("output/fig_phifo/c_gphifo_%s.png", name.Data()), "png");

  // c_tagged_flux->Print(Form("output/fig_phifo/c%s_tagged_flux.root", name.Data()), "root");
  // c_tagged_flux->Print(Form("output/fig_phifo/c%s_tagged_flux.eps", name.Data()), "eps");
  // c_tagged_flux->Print(Form("output/fig_phifo/c%s_tagged_flux.png", name.Data()), "png");
  // c_beame_tru->Print(Form("output/fig_phifo/c_beame_tru_%s.root", name.Data()), "root");
  // c_beame_tru->Print(Form("output/fig_phifo/c_beame_tru_%s.eps", name.Data()), "eps");
  // c_beame_tru->Print(Form("output/fig_phifo/c_beame_tru_%s.png", name.Data()), "png");

  outputfig->Print();

  fprintf(table_xsec_phifo, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
  fclose(table_xsec_phifo);
  gSystem->Exec(Form("pdflatex table_xsec_phifo_%s.tex", name.Data()));

   fprintf(table_xsec_phifo_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
   fclose(table_xsec_phifo_sys);
   gSystem->Exec(Form("pdflatex table_xsec_phifo_sys_%s.tex", name.Data()));

   // cgphifo_width->Print("output/fig_phifo/c_gphifo_width.root", "root");
   // cgphifo_width->Print("output/fig_phifo/c_gphifo_width.eps", "eps");
   // cgphifo_mean->Print("output/fig_phifo/c_gphifo_mean.root", "root");
   // cgphifo_mean->Print("output/fig_phifo/c_gphifo_mean.eps", "eps");

   // cgnophifo->Print("output/fig_phifo/c_gnophifo.root", "root");
   // cgnophifo->Print("output/fig_phifo/c_gnophifo.eps", "eps");

}