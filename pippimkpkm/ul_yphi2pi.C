#include "TCanvas.h"
#include "TFile.h"
#include "TH2F.h"
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
#include "TLine.h"
#include "scanphi.C"
#include "profileLikelihoodScan.C"

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

void ul_yphi2pi(TString name, int n2k=100, int n2pi2k=50, int ne=1, int nt=1) // , TString cut="&& kin_chisq<30 && abs(mm2)<0.015"
{
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", name.Data()));
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_yphi2pi_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/thrown_yphi2pi_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_11366_11555.root");
  if(name == "17") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_30274_31057.root");
  if(name == "18") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_40856_42577.root");
  if(name == "18l") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_50677_51768.root");

  TFile *outputfig = new TFile("output/fig_yphi2pi/ul_yphi2pi.root","UPDATE");//UPDATE

  FILE *table_ul_yphi2pi = fopen(Form("table_ul_yphi2pi_%s.tex", name.Data()),"w");
  fprintf(table_ul_yphi2pi,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\usepackage{pbox}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections and upper limit}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & \\pbox{10cm}{$N_{generated}$\\\\(MC)} & \\pbox{10cm}{$N_{measured}$\\\\(MC)} & \\pbox{10cm}{$N_{measured}$\\\\(Data)} & $\\epsilon$ [$\\%$] & $\\sigma$ [nb] & $90\\%$ CL limit [nb] \\\\ \n \\hline\n");
  //  \\pbox{10cm}{$90\\%$ CL limit [nb]\\\\(Bayesian)} & 

  FILE *table_ul_yphi2pi_sys = fopen(Form("table_ul_yphi2pi_sys_%s.tex", name.Data()),"w");
  fprintf(table_ul_yphi2pi_sys,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{makecell} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\small\n \\centering\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\n \\hline\n \\thead{Sig mean\\\\$\\pm 0.01\\%$ PDG error} & \\thead{Sig width\\\\$\\pm 0.012\\%$ PDG error} & \\thead{Bkg deg\\\\(3,4,5)}& \\thead{Fit range\\\\(2; 3, 2.9, 3.1)} & \\thead{bining\\\\(90, 100, 110)} & \\thead{Accidentals\\\\(2, 4, 8)} & \\thead{$\\Delta T_{RF-TOF}$ proton\\\\($\\pm$0.2, $\\pm$0.3, $\\pm$0.4)} & \\thead{$\\chi^{2}$ KinFit\\\\(45, 55, 65)} & \\thead{$MM^{2}$\\\\($\\pm$0.025, $\\pm$0.035, $\\pm$0.045)} \\\\ \n \\hline\n");

  ofstream test;
  test.open ("test.txt");

  RooWorkspace w("w", kTRUE);

  double mkk_min    = 0.99, mkk_max    = 1.2;
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2;
  double m2pi2k_min2 = 2.0,  m2pi2k_max2 = 3.0; // 3.0, 2.9, 3.1

  // gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetLabelSize(0.07, "XYZ");

  // ++++++++++++++++++++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // ++++++++++++++++++++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");  

  // *********************** Y(2175) MC *********************************
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts", "c_YMass_postcuts", 900, 600);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_postcuts", ";m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 100, m2pi2k_min2, m2pi2k_max2); //[2, 2.5]
  tmc->Project("h_YMass_postcuts", "kpkmpippim_mf", "w8*(kpkmpippim_uni && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)"); // is_truecombo && kpkm_mf>1.005 && kpkm_mf<1.035
  // && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(kpkmpippim_uni && abs(rf_dt)<1.5*4.008, w4*(kpkmpippim_uni && abs(rf_dt)<2.5*4.008, w8*(kpkmpippim_uni
  c_YMass_postcuts->cd();
  h_YMass_postcuts->Draw("e"); //"hist"

  TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", m2pi2k_min2, m2pi2k_max2);//pol4(4)
  // TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  fsb_mc->SetLineColor(2);
  fsb_mc->SetParameters(h_YMass_postcuts->GetMaximum(), 2.195, 0.025, 0.083, 1,1,1,1,1);//2.175, 0.004, 0.002, 1,1,1,1,1
  // fsb_mc->SetParLimits(0, 0, 10000);
  fsb_mc->SetParLimits(1, 2.190, 2.196); // 1.018, 1.021
  // fsb_mc->SetParLimits(2, 0.01, 0.1);
  // fsb_mc->SetParLimits(3, 0.079-0.14, 0.079+0.14);
  fsb_mc->FixParameter(3, 0.083); //0.079
  if(name == "16") fsb_mc->FixParameter(1, 2.192);

  TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", m2pi2k_min2, m2pi2k_max2);
  fs_mc->SetLineColor(4);

  TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", m2pi2k_min2, m2pi2k_max2); //pol4(4)
  fb_mc->SetLineColor(28);
  fb_mc->SetLineStyle(2);

  h_YMass_postcuts->Fit("fsb_mc", "", "", m2pi2k_min2, m2pi2k_max2);
  double par_mc[fsb_mc->GetNpar()]; //6
  fsb_mc->GetParameters(&par_mc[0]);
  fs_mc->SetParameters(&par_mc[0]);
  fb_mc->SetParameters(&par_mc[4]); //3

  fs_mc->Draw("same");
  fb_mc->Draw("same");
  fsb_mc->Draw("same");

  double N_y_mc = fs_mc->Integral(m2pi2k_min2, m2pi2k_max2) / h_YMass_postcuts->GetBinWidth(1);
  double dN_y_mc = N_y_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

  TLegend *l_y_mc = new TLegend(0.55, 0.27, 0.66, 0.45);
  l_y_mc->SetTextSize(0.05);
  l_y_mc->SetFillColor(kWhite);
  l_y_mc->SetLineColor(kWhite);
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  l_y_mc->AddEntry("fsb_mc", "Total", "l");
  l_y_mc->AddEntry("fs_mc", "Signal", "l");
  l_y_mc->AddEntry("fb_mc", "Bkg", "l");
  l_y_mc->Draw();

  TLatex lat_mc;
  lat_mc.SetTextSize(0.05);
  lat_mc.SetTextAlign(13); //align at top
  lat_mc.SetNDC();
  lat_mc.SetTextColor(kBlue);
  lat_mc.DrawLatex(0.55, 0.82, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
  lat_mc.DrawLatex(0.55, 0.75, Form("N_{sig} = %0.2f#pm%0.2f", N_y_mc, dN_y_mc));
  lat_mc.DrawLatex(0.55, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
  lat_mc.DrawLatex(0.55, 0.61, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
  lat_mc.DrawLatex(0.55, 0.54, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));
  // lat_mc.Draw("same");

  c_YMass_postcuts->Print(Form("output/fig_yphi2pi/cmc_YMass_postcuts_fitted_%s.root",name.Data()), "root");
  c_YMass_postcuts->Print(Form("output/fig_yphi2pi/cmc_YMass_postcuts_fitted_%s.eps",name.Data()), "eps");
  c_YMass_postcuts->Print(Form("output/fig_yphi2pi/cmc_YMass_postcuts_fitted_%s.png",name.Data()), "png");

  // cout<<"=============================  no problem up to here ! ========================"<<endl;

  // ======================================== Y vs. Phi(1020) ===============================================

  TCanvas *cgphiy = new TCanvas("cgphiy", "cgphiy", 900, 600); //1500, 800
  TGraphErrors *gphiy = scanphi("y", "yphi2pi", Form("%s", name.Data()), n2k, n2pi2k);

  cgphiy->cd();
  gphiy->Draw("ALP");

  double wid = gphiy->GetX()[1] - gphiy->GetX()[0];
  TH1F *hgphiy = new TH1F("hgphiy", ";m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi#pi^{+}#pi^{-}}", gphiy->GetN(), gphiy->GetX()[0] - 0.5 * wid, gphiy->GetX()[gphiy->GetN() - 1] + 0.5 * wid);
  for (int i = 0; i < gphiy->GetN(); ++i)
  {
    hgphiy->SetBinContent(i + 1, gphiy->GetY()[i]);
    hgphiy->SetBinError(i + 1, gphiy->GetEY()[i]);
  }
  hgphiy->SetTitleSize(0.07, "XYZ");
  hgphiy->SetLabelSize(0.07, "XYZ");
  hgphiy->Draw("e");
  outputfig->cd();
  hgphiy->Write(Form("h%s_phiy", name.Data()), TObject::kWriteDelete);
  hgphiy->SetMinimum(-100.);
  // save the histo to root file
  w.factory(Form("Voigtian::sig_y_data(m_y_data[%f, %f],mean_y_data[%f],width_y_data[%f],sigma_y_data[%f])", m2pi2k_min2, m2pi2k_max2, fsb_mc->GetParameter(1), fsb_mc->GetParameter(3), fsb_mc->GetParameter(2)));
  w.factory("Chebychev::bkg_y_data(m_y_data,{c0_y_data[-10,+10], c1_y_data[-10,+10], c2_y_data[-10,+10], c3_y_data[-10,+10]})"); // c0_y_data[-10,+10], c1_y_data[-10,+10], c2_y_data[-10,+10], c3_y_data[-10,+10]
  w.factory("SUM:::model_y_data(nsig_y_data[0,-1000000,+1000000]*sig_y_data, nbkg_y_data[0,+1000000]*bkg_y_data)"); 
  //nsig_y_data[0,-1000000,+1000000]*sig_y_data, nbkg_y_data[0,+1000000]*bkg_y_data
  w.var("m_y_data")->SetTitle("m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2})");
  RooDataHist dh_y_data("dh_y_data", "dh_y_data", *w.var("m_y_data"), Import(*hgphiy));
  RooPlot *fr_y_data = w.var("m_y_data")->frame(Title(" "));
  w.pdf("model_y_data")->fitTo(dh_y_data);

  // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
  dh_y_data.plotOn(fr_y_data, RooFit::Name("ldh_y_data"));
  w.pdf("model_y_data")->plotOn(fr_y_data, Components(*w.pdf("sig_y_data")), LineColor(kRed), RooFit::Name("lsig_y_data"));
  w.pdf("model_y_data")->plotOn(fr_y_data, Components(*w.pdf("bkg_y_data")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_y_data"));
  w.pdf("model_y_data")->plotOn(fr_y_data, RooFit::Name("lmodel_y_data"));
  // w.pdf("model_y_data")->paramOn(fr_y_data, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_y_data"), *w.var("nbkg_y_data")))); //,*w.var("mean_y_data"),*w.var("width_y_data"),*w.var("sigma_y_data"))));
  fr_y_data->Draw();

  double N_y_data = w.var("nsig_y_data")->getVal();
  double dN_y_data = w.var("nsig_y_data")->getError();

  TLatex lat_phiye;
  lat_phiye.SetTextSize(0.05);
  lat_phiye.SetTextAlign(13); //align at top
  lat_phiye.SetNDC();
  lat_phiye.SetTextColor(kBlue);
  // lat_phiye.DrawLatex(0.60, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
  lat_phiye.DrawLatex(0.60, 0.40, Form("N_{sig} = %0.2f#pm%0.2f", N_y_data, dN_y_data));
  lat_phiye.DrawLatex(0.60, 0.34, Form("#mu = %0.3f#pm%0.3f", w.var("mean_y_data")->getVal(), w.var("mean_y_data")->getError()));
  lat_phiye.DrawLatex(0.60, 0.28, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_y_data")->getVal(), w.var("sigma_y_data")->getError()));
  lat_phiye.DrawLatex(0.60, 0.22, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_y_data")->getVal(), w.var("width_y_data")->getError()));

  // +++++++++++++++++++ Systematic Errors +++++++++++++++++++

  // if(name == "16")
  double arr_sig_mean_16[3] = {0.26, 0.27, 0.23}; //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
  double arr_sig_width_16[3] = {0.26, 0.23, 0.29}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
  double arr_bkg_poln_16[3] = {0.26, 0.16, 4.57}; //pol:4,3,5
  double arr_fit_range_16[3] = {0.26, 2.87, 0.16}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
  double arr_binning_16[3] = {0.26, 0.36, 0.28}; //binning: 100,90,110
  double arr_beambun_16[3] = {0.26, 0.38, 0.29}; //bunches: w8, w2, w4
  double arr_p_dttof_16[3] = {0.29, 0.33, 0.39}; // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_16[3] = {0.29, 0.19, 0.12}; // kin_chisq<55, 45, 65
  double arr_mm2_16[3] = {0.29, 0.15, 0.31}; // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "17")
  double arr_sig_mean_17[3] = {0.27, 0.30, 0.22}; //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
  double arr_sig_width_17[3] = {0.27, 0.24, 0.30}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
  double arr_bkg_poln_17[3] = {0.27, 0.21, 0.30}; //pol:4,3,5
  double arr_fit_range_17[3] = {0.27, 0.33, 0.27}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
  double arr_binning_17[3] = {0.27, 0.24, 0.32}; //binning: 100,90,110
  double arr_beambun_17[3] = {0.27, 0.27, 0.31}; //bunches: w8, w2, w4
  double arr_p_dttof_17[3] = {0.33, 0.33, 0.35}; // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_17[3] = {0.33, 0.33, 0.36}; // kin_chisq<55, 45, 65
  double arr_mm2_17[3] = {0.33, 0.27, 0.31}; // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18")
  double arr_sig_mean_18[3] = {0.07, 0.15, -0.01}; //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
  double arr_sig_width_18[3] = {0.07, 0.05, 0.09}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
  double arr_bkg_poln_18[3] = {0.07, -0.01, -0.14};  //pol:4,3,5
  double arr_fit_range_18[3] = {0.07, 0.00, 0.04}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
  double arr_binning_18[3] = {0.07, 0.10, 0.21}; //binning: 100,90,110
  double arr_beambun_18[3] = {0.07, 0.10, 0.08}; //bunches: w8, w2, w4
  double arr_p_dttof_18[3] = {0.12, 0.11, 0.09}; // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18[3] = {0.12, 0.17, 0.15}; // kin_chisq<55, 45, 65
  double arr_mm2_18[3] = {0.12, 0.11, 0.16}; // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18l")
  double arr_sig_mean_18l[3] = {0.11, 0.06, 0.14}; //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
  double arr_sig_width_18l[3] = {0.11, 0.10, 0.12}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
  double arr_bkg_poln_18l[3] = {0.11, 0.02, 0.30}; //pol:4,3,5
  double arr_fit_range_18l[3] = {0.11, 0.32, 0.04}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
  double arr_binning_18l[3] = {0.11, 0.10, 0.12}; //binning: 100,90,110
  double arr_beambun_18l[3] = {0.11, 0.15, 0.18}; //bunches: w8, w2, w4
  double arr_p_dttof_18l[3] = {0.08, 0.09, 0.07}; // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18l[3] = {0.08, 0.11, 0.03}; // kin_chisq<55, 45, 65
  double arr_mm2_18l[3] = {0.08, 0.09, 0.07}; // abs(mm2)<0.035, 0.025, 0.045

  double arr_sig_mean[3], arr_sig_width[3], arr_bkg_poln[3], arr_fit_range[3], arr_binning[3], arr_beambun[3], arr_p_dttof[3], arr_kin_chisq[3], arr_mm2[3];

  if(name == "16")
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
  if(name == "17")
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
  if(name == "18")
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
  if(name == "18l")
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
 
  double err_sig_mean = TMath::RMS(3, arr_sig_mean) / arr_sig_mean[0];
  double err_sig_width = TMath::RMS(3, arr_sig_width) / arr_sig_width[0];
  double err_bkg_poln = TMath::RMS(3, arr_bkg_poln) / arr_bkg_poln[0];
  double err_fit_range = TMath::RMS(3, arr_fit_range) / arr_fit_range[0];
  double err_binning = TMath::RMS(3, arr_binning) / arr_binning[0];
  double err_beambun = TMath::RMS(3, arr_beambun) / arr_beambun[0];
  double err_p_dttof = TMath::RMS(3, arr_p_dttof) / arr_p_dttof[0];
  double err_kin_chisq = TMath::RMS(3, arr_kin_chisq) / arr_kin_chisq[0];
  double err_mm2 = TMath::RMS(3, arr_mm2) / arr_mm2[0];

  double err_sys = TMath::Sqrt(err_sig_mean*err_sig_mean + err_sig_width*err_sig_width + err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range + err_binning*err_binning + err_beambun*err_beambun + err_p_dttof*err_p_dttof + err_kin_chisq*err_kin_chisq + err_mm2*err_mm2);

  // ++++++++++++++++++++++++++++ Cross-section  +++++++++++++++++++++++
  double N_y_gen = h_beame_tru->Integral(1, 10);
  double eff_y = N_y_mc / N_y_gen; // Efficiency = N_observed/N_generated
  double deff_y = eff_y * (dN_y_mc / N_y_mc);
  double lumi_y = h_tagged_flux->Integral(1, 10); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.273 barns^-1
  double T = 1.273;                               //Target thickness [b^{-1}]
  double br_phi = 0.492;
  double dbr_phi = 0.005;
  double xsec_y = 1e9 * N_y_data / (eff_y * lumi_y * T * br_phi);
  double dxsec_y_stat = xsec_y * (dN_y_data / N_y_data);
  double dxsec_y_sys = xsec_y * TMath::Sqrt((err_sys*err_sys) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
  // double dxsec_y_sys = xsec_y * TMath::Sqrt((err_sys / N_y_data) * (err_sys / N_y_data) + (deff_y / eff_y) * (deff_y / eff_y) + (dbr_phi / br_phi) * (dbr_phi / br_phi));
  double dxsec_y_tot = TMath::Sqrt(dxsec_y_stat*dxsec_y_stat + dxsec_y_sys*dxsec_y_sys);

  test << std::setprecision(2) << std::fixed;
  test << "xsec_y = " << xsec_y << " | sys: " << dxsec_y_sys << endl;

  // ++++++++++++++++++++++++++++ Upper Limit  90% CL  +++++++++++++++++

  // ++++++++++ gaussian method
  // double ul = mu + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, sig, mu), sig);
  // double ul = xsec_y + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, dxsec_y_tot, xsec_y), dxsec_y_tot);

  // ++++++++++ profile likelihood method
  // double dN_y_tot = sqrt(dN_y_data*dN_y_data + err_sys*err_sys);
  double dN_y_tot = dN_y_data;
  TGraph *gprof = profileLikelihoodScan(dh_y_data, *(w.pdf("model_y_data")), w.var("nsig_y_data"), N_y_data - 5 * dN_y_tot, N_y_data + 5 * dN_y_tot); //N_y_data - 5 * dN_y_tot, N_y_data + 5 * dN_y_tot
  TCanvas *cprof = new TCanvas("cprof", "", 900, 600);
  cprof->cd();
  cprof->SetGrid();
  gprof->SetTitle("; Yield_{#phi#pi^{+}#pi^{-}}; -Log(L)");
  gprof->SetLineWidth(2);
  gprof->SetMarkerStyle(20);
  gprof->SetMarkerSize(1.5);
  gprof->Draw("AP");
  gprof->Fit("pol2", "R", "", gprof->GetHistogram()->GetXaxis()->GetXmin(), gprof->GetHistogram()->GetXaxis()->GetXmax());

  TCanvas *cprofxsec = new TCanvas("cprofxsec", "", 900, 600);
  TGraph *gprofxsec = new TGraph();
  for (int i = 0; i < gprof->GetN(); i++)
  {
    double x_gprof, y_gprof;
    gprof->GetPoint(i, x_gprof, y_gprof);
    double xsec_gprof = 1e9 * x_gprof / (eff_y * lumi_y * T * br_phi);
    cout << "i = " << i << " | x_gprof = " << x_gprof << " | y_gprof = " << y_gprof << " | xsec_gprof = " << xsec_gprof << endl;
    gprofxsec->SetPoint(i, xsec_gprof, exp(-y_gprof));
  }
  cprofxsec->cd();
  gprofxsec->SetTitle("; #sigma [nb]; Likelihood");
  gprofxsec->SetLineWidth(2);
  gprofxsec->SetMarkerStyle(20);
  gprofxsec->SetMarkerSize(1.5);
  gprofxsec->Draw("AP");
  TF1 *func_profxsec = new TF1("func_profxsec", "gaus", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());
  gprofxsec->Fit(func_profxsec, "R", "", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());

  //creat a convolution between the profile likelihood and gaus with same mean and a s.t.d as sys err
  TCanvas *c_twogaus_conv = new TCanvas("c_twogaus_conv", "", 900, 600);
  c_twogaus_conv->cd();
  TF1 *func_twogaus_conv = new TF1("func_twogaus_conv", "gaus", -1.5, 1.5);
  func_twogaus_conv->SetParameters(func_profxsec->GetParameter(0), func_profxsec->GetParameter(1), dxsec_y_tot);
  func_twogaus_conv->SetTitle(";#sigma (nb);Likelihood");
  func_twogaus_conv->Draw();

  double ul = func_twogaus_conv->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, func_twogaus_conv->GetParameter(2), func_twogaus_conv->GetParameter(1)), func_twogaus_conv->GetParameter(2)); //bayesian
  // double ul = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, func_profxsec->GetParameter(2), func_profxsec->GetParameter(1)), func_profxsec->GetParameter(2));//bayesian
  double ul2 = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9, func_profxsec->GetParameter(2)); //frequentist
  TLine *l_mean = new TLine(func_profxsec->GetParameter(1), 0, func_profxsec->GetParameter(1), 1);
  l_mean->SetLineWidth(2);
  l_mean->SetLineStyle(10);
  // l_mean->Draw("same");
  TLine *l_sigma = new TLine(func_profxsec->GetParameter(1), 0.5, func_profxsec->GetParameter(1) + func_profxsec->GetParameter(2), 0.5);
  l_sigma->SetLineWidth(2);
  l_sigma->SetLineStyle(10);
  // l_sigma->Draw("same");
  TLine *l_ul = new TLine(ul, 0, ul, 0.45);
  l_ul->SetLineWidth(2);
  l_ul->SetLineColor(kBlue);
  l_ul->Draw("same");
  // TLine *l_ul2 = new TLine(ul2, 0, ul2, 0.45);
  // l_ul2->SetLineWidth(2);
  // l_ul2->SetLineColor(28);
  // l_ul2->Draw("same");
  //  +++++++++++++++++++

  fprintf(table_ul_yphi2pi, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), N_y_gen, N_y_mc, dN_y_mc, N_y_data, dN_y_data, eff_y * 100, deff_y * 100, xsec_y, dxsec_y_stat, dxsec_y_sys, ul);

  fprintf(table_ul_yphi2pi_sys, "%0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", err_sig_mean * 100, err_sig_width * 100, err_bkg_poln * 100, err_fit_range * 100, err_binning * 100, err_beambun * 100, err_p_dttof * 100, err_kin_chisq * 100, err_mm2 * 100);

  cprof->Print(Form("output/fig_yphi2pi/c_prof_%s.root", name.Data()), "root");
  cprof->Print(Form("output/fig_yphi2pi/c_prof_%s.eps", name.Data()), "eps");
  cprof->Print(Form("output/fig_yphi2pi/c_prof_%s.png", name.Data()), "png");
  cprofxsec->Print(Form("output/fig_yphi2pi/c%s_profxsec.root", name.Data()), "root");
  cprofxsec->Print(Form("output/fig_yphi2pi/c%s_profxsec.eps", name.Data()), "eps");
  cprofxsec->Print(Form("output/fig_yphi2pi/c%s_profxsec.png", name.Data()), "png");
  c_twogaus_conv->Print(Form("output/fig_yphi2pi/c%s_twogaus_conv.root", name.Data()), "root");
  c_twogaus_conv->Print(Form("output/fig_yphi2pi/c%s_twogaus_conv.eps", name.Data()), "eps");
  c_twogaus_conv->Print(Form("output/fig_yphi2pi/c%s_twogaus_conv.png", name.Data()), "png");

  cgphiy->Print(Form("output/fig_yphi2pi/c_gphiy_%s.root", name.Data()), "root");
  cgphiy->Print(Form("output/fig_yphi2pi/c_gphiy_%s.eps", name.Data()), "eps");
  cgphiy->Print(Form("output/fig_yphi2pi/c_gphiy_%s.png", name.Data()), "png");
  c_tagged_flux->Print(Form("output/fig_yphi2pi/c%s_tagged_flux.root", name.Data()), "root");
  c_tagged_flux->Print(Form("output/fig_yphi2pi/c%s_tagged_flux.eps", name.Data()), "eps");
  c_tagged_flux->Print(Form("output/fig_yphi2pi/c%s_tagged_flux.png", name.Data()), "png");
  c_beame_tru->Print(Form("output/fig_yphi2pi/c%s_beame_tru.root", name.Data()), "root");
  c_beame_tru->Print(Form("output/fig_yphi2pi/c%s_beame_tru.eps", name.Data()), "eps");
  c_beame_tru->Print(Form("output/fig_yphi2pi/c%s_beame_tru.png", name.Data()), "png");

  outputfig->Print();

  fprintf(table_ul_yphi2pi, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
  fclose(table_ul_yphi2pi);
  gSystem->Exec(Form("pdflatex table_ul_yphi2pi_%s.tex", name.Data()));

  fprintf(table_ul_yphi2pi_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
  fclose(table_ul_yphi2pi_sys);
  gSystem->Exec(Form("pdflatex table_ul_yphi2pi_sys_%s.tex", name.Data()));
}