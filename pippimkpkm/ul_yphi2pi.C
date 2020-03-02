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

  TFile *fdata = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_%s_flat.root", name.Data()));
  TFile *fmc = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_yphi2pi_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/thrown_yphi2pi_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_11366_11555.root");
  if(name == "17") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  if(name == "18") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_40856_42577.root");
 

  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/ul_yphi2pi.root","UPDATE");

  FILE *table_ul_yphi2pi = fopen(Form("table_ul_yphi2pi_%s.tex", name.Data()),"w");
  fprintf(table_ul_yphi2pi,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\usepackage{pbox}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections and upper limit}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & \\pbox{10cm}{$N_{generated}$\\\\(MC)} & \\pbox{10cm}{$N_{measured}$\\\\(MC)} & \\pbox{10cm}{$N_{measured}$\\\\(Data)} & $\\epsilon$ [$\\%$] & $\\sigma$ [nb] & $90\\%$ CL limit [nb] \\\\ \n \\hline\n");
  //  \\pbox{10cm}{$90\\%$ CL limit [nb]\\\\(Bayesian)} & 

  FILE *table_ul_yphi2pi_sys = fopen(Form("table_ul_yphi2pi_sys_%s.tex", name.Data()),"w");
  fprintf(table_ul_yphi2pi_sys,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors}\n \\begin{tabular}{|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ $\\%$ ] & Fit Range [$\\%$] & $\\phi$-mass bins [$\\%$]  & Y Mean [$\\%$] & Y width [$\\%$] \\\\ \n \\hline\n");

  RooWorkspace w("w", kTRUE);

  double mkk_min    = 0.99, mkk_max    = 1.2;
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2;

  // gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // ++++++++++++++++++++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.0);
  h_tagged_flux->Draw("e");

  // ++++++++++++++++++++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.0);
  h_beame_tru->Draw("e");  

  // // *********************** Phi(1020) MC *********************************
  // TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  // TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", "MC signal; m_{K^{+}K^{-}} [GeV/c^{2}]; Counts", 200, 1.005, 1.035);
  // tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");//+cutlist+" && kin_chisq<25)"
  // c_PhiMass_postcuts->cd();
  // h_PhiMass_postcuts->Draw("e");

  // // w.factory("BreitWigner::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004])");
  // // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  // w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.001,0.01])"); //sigma_PhiMass_mc[0.001,0.01], mean_PhiMass_mc[1.015,1.022]
  // w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-10,10], c1_PhiMass_mc[-10,10]})");
  // w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000000]*sig_PhiMass_mc, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)");//, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); //nsig[0,100000000]*sig2,
  // w.var("m_PhiMass_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
  // RooDataHist dh_PhiMass_mc("dh_PhiMass_mc", "dh_PhiMass_mc", *w.var("m_PhiMass_mc"), Import(*h_PhiMass_postcuts));
  // RooPlot *fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title("K^{+}K^{-}"));
  // // fr_PhiMass_mc->SetTitleOffset(0.90, "X");
  // // fr_PhiMass_mc->SetTitleSize(0.06, "XYZ");
  // // fr_PhiMass_mc->SetLabelSize(0.06, "xyz");
  // w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc);

  // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
  // dh_PhiMass_mc.plotOn(fr_PhiMass_mc, RooFit::Name("ldh_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig_PhiMass_mc")), LineColor(kRed), RooFit::Name("lsig_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("bkg_PhiMass_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, RooFit::Name("lmodel_PhiMass_mc"));
  // // w.pdf("model_PhiMass_mc")->paramOn(fr_PhiMass_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_PhiMass_mc"), *w.var("nbkg_PhiMass_mc")))); //,*w.var("mean_PhiMass_mc"),*w.var("width_PhiMass_mc"),*w.var("sigma_PhiMass_mc"))));
  // fr_PhiMass_mc->Draw();

  // TLegend *l_phi_mc = new TLegend(0.2, 0.65, 0.4, 0.85);
  // l_phi_mc->SetFillColor(kWhite);
  // l_phi_mc->SetLineColor(kWhite);
  // // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lmodel_PhiMass_mc"), "total", "l");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lsig_PhiMass_mc"), "Voigtian", "l");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lbkg_PhiMass_mc"), "pol 2nd", "l");
  // l_phi_mc->Draw();

  // TLatex lat_PhiMass_mc;
  // lat_PhiMass_mc.SetTextSize(0.05);
  // lat_PhiMass_mc.SetTextAlign(13); //align at top
  // lat_PhiMass_mc.SetNDC();
  // lat_PhiMass_mc.SetTextColor(kBlue);
  // lat_PhiMass_mc.DrawLatex(0.62, 0.87, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_PhiMass_mc")->getVal(), w.var("nsig_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.78, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_PhiMass_mc")->getVal(), w.var("nbkg_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.68, Form("#mu = %0.3f#pm%0.3f", w.var("mean_PhiMass_mc")->getVal(), w.var("mean_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.58, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_PhiMass_mc")->getVal(), w.var("width_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.48, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_PhiMass_mc")->getVal(), w.var("sigma_PhiMass_mc")->getError()));

  // // TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 1.005, 1.035);
  // // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  // // fsb->SetLineColor(2);
  // // fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1, 1, 1);
  // // fsb->SetParLimits(0, 0, 10000);
  // // fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
  // // fsb->FixParameter(2, 0.004);        //fsb->SetParLimits(2, 0.008, 0.010);
  // // fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

  // // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 1.005, 1.035);
  // // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  // // fs->SetLineColor(4);

  // // TF1 *fb = new TF1("fb", "pol4(4)", 1.005, 1.035); //pol2(3)
  // // fb->SetLineColor(28);
  // // fb->SetLineStyle(2);

  // // h_PhiMass_postcuts->Fit("fsb", "", "", 1.005, 1.035);
  // // double par[8]; //6
  // // fsb->GetParameters(&par[0]);
  // // fs->SetParameters(&par[0]);
  // // fb->SetParameters(&par[4]); //3

  // // fs->Draw("same");
  // // fb->Draw("same");

  // c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_PhiMass_postcuts_fitted.root", "root");
  // c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_PhiMass_postcuts_fitted.eps", "eps");
  // c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_PhiMass_postcuts_fitted.png", "png");

  // *********************** Y(2175) MC *********************************
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts", "c_YMass_postcuts", 900, 600);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_postcuts", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 2, 2.5); //[2, 2.5]
  tmc->Project("h_YMass_postcuts", "kpkmpippim_mf", "w8*(kpkmpippim_uni && is_truecombo && kpkm_mf>1.005 && kpkm_mf<1.035)");
  // TH1F *h_YMass_postcuts = new TH1F("h_YMass_beambunchcut", ";m_{#phif_{0}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);
  // t->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && abs(pippim_mf-0.99)<0.1)");
  c_YMass_postcuts->cd();
  h_YMass_postcuts->Draw("e"); //"hist"

  // // w.factory("BreitWigner::sig_YMass_mc(m_YMass_mc[1.7, 3.2],mean_YMass_mc[2.0,2.5],width_YMass_mc[0.079])");
  // // w.factory("Voigtian::sig_YMass_mc(m_YMass_mc[1.7, 3.2],mean_YMass_mc[2.0,2.5],width_YMass_mc[0.03,0.1],sigma_YMass_mc[0.001,0.01])"); //width:0.079 sigma_YMass_mc[0.001,0.01], mean_YMass_mc[1.015,1.022]
  // w.factory("Gaussian::sig_YMass_mc1(m_YMass_mc[1.7, 3.2],mean_YMass_mc1[2.0,2.2],sigma_YMass_mc1[0.01,0.1])"); //width:0.079 sigma_YMass_mc[0.001,0.01], mean_YMass_mc[1.015,1.022]
  // w.factory("Gaussian::sig_YMass_mc2(m_YMass_mc,mean_YMass_mc2[2.2,2.5],sigma_YMass_mc2[0.01,0.1])"); //width:0.079 sigma_YMass_mc[0.001,0.01], mean_YMass_mc[1.015,1.022]
  // w.factory("Chebychev::bkg_YMass_mc(m_YMass_mc,{c0_YMass_mc[-10,10], c1_YMass_mc[-10,10]})");
  // // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  // w.factory("SUM:::model_YMass_mc(nsig_YMass_mc[0,100000000]*sig_YMass_mc1,nsig_YMass_mc[0,100000000]*sig_YMass_mc2, nbkg_YMass_mc[0,100000000]*bkg_YMass_mc)");//, nbkg_YMass_mc[0,100000000]*bkg_YMass_mc)"); //nsig[0,100000000]*sig2,
  // w.var("m_YMass_mc")->SetTitle("m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]");
  // RooDataHist dh_YMass_mc("dh_YMass_mc", "dh_YMass_mc", *w.var("m_YMass_mc"), Import(*h_YMass_postcuts));
  // RooPlot *fr_YMass_mc = w.var("m_YMass_mc")->frame(Title("K^{+}K^{-}"));
  // // fr_YMass_mc->SetTitleOffset(0.90, "X");
  // // fr_YMass_mc->SetTitleSize(0.06, "XYZ");
  // // fr_YMass_mc->SetLabelSize(0.06, "xyz");
  // w.pdf("model_YMass_mc")->fitTo(dh_YMass_mc);

  // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
  // dh_YMass_mc.plotOn(fr_YMass_mc, RooFit::Name("ldh_YMass_mc"));
  // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig_YMass_mc1")), LineColor(kRed), RooFit::Name("lsig_YMass_mc1"));
  // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig_YMass_mc2")), LineColor(kRed), RooFit::Name("lsig_YMass_mc2"));
  // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("bkg_YMass_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_YMass_mc"));
  // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, RooFit::Name("lmodel_YMass_mc"));
  // // w.pdf("model_YMass_mc")->paramOn(fr_YMass_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_YMass_mc"), *w.var("nbkg_YMass_mc")))); //,*w.var("mean_YMass_mc"),*w.var("width_YMass_mc"),*w.var("sigma_YMass_mc"))));
  // fr_YMass_mc->Draw();

  // TLatex lat_YMass_mc;
  // lat_YMass_mc.SetTextSize(0.05);
  // lat_YMass_mc.SetTextAlign(13); //align at top
  // lat_YMass_mc.SetNDC();
  // lat_YMass_mc.SetTextColor(kBlue);
  // lat_YMass_mc.DrawLatex(0.65, 0.88, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_YMass_mc")->getVal(), w.var("nsig_YMass_mc")->getError()));
  // lat_YMass_mc.DrawLatex(0.65, 0.78, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_YMass_mc")->getVal(), w.var("nbkg_YMass_mc")->getError()));
  // lat_YMass_mc.DrawLatex(0.65, 0.68, Form("#mu = %0.3f#pm%0.3f", w.var("mean_YMass_mc1")->getVal(), w.var("mean_YMass_mc1")->getError()));
  // // lat_YMass_mc.DrawLatex(0.65, 0.58, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_YMass_mc")->getVal(), w.var("width_YMass_mc")->getError()));
  // lat_YMass_mc.DrawLatex(0.65, 0.48, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_YMass_mc1")->getVal(), w.var("sigma_YMass_mc1")->getError()));

  TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 2, 2.5); // + pol2(4) , [2, 2.5]
  // TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  fsb_mc->SetLineColor(2);
  fsb_mc->SetParameters(700, 2.175, 0.004, 0.002, 1, 1,1,1,1);//, 1, 1,1
  fsb_mc->SetParLimits(0, 0, 10000);
  fsb_mc->SetParLimits(1, 2.0, 2.2); // 1.018, 1.021
  fsb_mc->SetParLimits(2, 0.01, 0.1);
  // fsb_mc->FixParameter(2, 0.079);
  // fsb_mc->SetParLimits(3, 0.079-0.14, 0.079+0.14);
  fsb_mc->FixParameter(3, 0.079);

  // TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", 2, 2.5);
  fs_mc->SetLineColor(4);

  TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", 2, 2.5); //pol2(3)
  fb_mc->SetLineColor(28);
  fb_mc->SetLineStyle(2);

  h_YMass_postcuts->Fit("fsb_mc", "", "", 2, 2.5);
  double par_mc[8]; //6
  fsb_mc->GetParameters(&par_mc[0]);
  fs_mc->SetParameters(&par_mc[0]);
  fb_mc->SetParameters(&par_mc[4]); //3

  fs_mc->Draw("same");
  fb_mc->Draw("same");
  fsb_mc->Draw("same");

  double N_y_mc = fs_mc->Integral(2, 2.5) / h_YMass_postcuts->GetBinWidth(1);
  double dN_y_mc = N_y_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

  TLegend *l_y_mc = new TLegend(0.2, 0.65, 0.35, 0.85);
  l_y_mc->SetTextSize(0.05);
  l_y_mc->SetFillColor(kWhite);
  l_y_mc->SetLineColor(kWhite);
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  l_y_mc->AddEntry("fsb_mc", "total", "l");
  l_y_mc->AddEntry("fs_mc", "Voigtian", "l");
  l_y_mc->AddEntry("fb_mc", "pol 4^{th}", "l");
  l_y_mc->Draw();

  TLatex lat_mc;
  lat_mc.SetTextSize(0.04);
  lat_mc.SetTextAlign(13); //align at top
  lat_mc.SetNDC();
  lat_mc.SetTextColor(kBlue);
  lat_mc.DrawLatex(0.6, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
  lat_mc.DrawLatex(0.6, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_y_mc, dN_y_mc));
  lat_mc.DrawLatex(0.6, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
  lat_mc.DrawLatex(0.6, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
  lat_mc.DrawLatex(0.6, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));
  // lat_mc.Draw("same");

  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_YMass_postcuts_fitted.root", "root");
  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_YMass_postcuts_fitted.eps", "eps");
  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/cmc_YMass_postcuts_fitted.png", "png");

  // cout<<"=============================  no problem up to here ! ========================"<<endl;

  // ======================================== Y vs. Phi(1020) ===============================================
  // root -l 'ul_yphi2pi.C+(100,50,1,1)'

  // TCanvas *cphiy=new TCanvas("cphiy","cphiy",900, 600);//1500, 800

  // TCanvas *cphiy1 = new TCanvas("cphiy1","cphiy1", 1500, 800);//1500, 800
  // cphiy1->Divide(5, 5);
  // TCanvas *cphiy2 = new TCanvas("cphiy2", "cphiy2", 1500, 800);
  // cphiy2->Divide(5, 5);

  TCanvas *cgphiy=new TCanvas("cgphiy","cgphiy",900, 600);//1500, 800
  TGraphErrors *gphiy = scanphi("y", "yphi2pi", Form("%s",name.Data()), n2k, n2pi2k);
  
  // gphiy = new TGraphErrors(); //n2pi2k
  // gphiy->SetMarkerStyle(20);
  // gphiy->SetMarkerSize(1.0);
  // gphiy->SetMarkerColor(1);
  // gphiy->SetTitle(Form("(%s);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}",name.Data()));

  // // TCanvas *cgnophiy=new TCanvas("cgnophiy","cgnophiy",1500, 800);
  // // TGraphErrors *gnophiy= new TGraphErrors;//(n2pi2k)
  // // gnophiy->SetMarkerStyle(20);
  // // gnophiy->SetMarkerSize(1.0);
  // // gnophiy->SetMarkerColor(2);

  // // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  // cphiy->cd();
  // TH2F *h1d2 = new TH2F("h1d2",Form("(%s);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]",name.Data()), n2pi2k, m2pi2k_min, m2pi2k_max, n2k, mkk_min, mkk_max);
  // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*(kpkm_uni || kpkmpippim_uni)");

  // h1d2->Draw("colz");


  // // gnophiy->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1, Eg2));

  // for (int i = 1; i <= n2pi2k; ++i) //n2pi2k
  // {
  //   // cout << i << " " << flush;

  //   if (i < 26)
  //     cphiy1->cd(i); //i
  //   if (i > 25 && i < 51)
  //     cphiy2->cd(i - 25);
  //   // if(i > 50 && i<76) cphiy3->cd(i-50);
  //   // if(i > 75) cphiy4->cd(i-75);

  //   TH1D *hphiy_py = h1d2->ProjectionY(Form("_hphiy_py_%d", i), i, i);

  //   // hphiy_py->Draw();
  //   // if (hphiy_py->Integral(1, 50) < 100)
  //   //   continue;

  //   // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  //   TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 0.98, 1.2);
  //   fsb->SetLineColor(2);
  //   fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1,1,1,1,1);
  //   fsb->SetParLimits(0, 0, 10000);
  //   fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
  //   fsb->FixParameter(2, w.var("sigma_PhiMass_mc")->getVal());   //fsb->SetParLimits(2, 0.008, 0.010); // 0.0042
  //   fsb->FixParameter(3, w.var("width_PhiMass_mc")->getVal());        //fsb->SetParLimits(3, 0.001,0.01);// 0.002

  //   // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  //   TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
  //   fs->SetLineColor(4);

  //   TF1 *fb = new TF1("fb", "pol4(4)", mkk_min, mkk_max); //3
  //   fb->SetLineColor(28);
  //   fb->SetLineStyle(2);

  //   hphiy_py->Fit("fsb", "", "", mkk_min, mkk_max);
  //   double par[10]; //6
  //   fsb->GetParameters(&par[0]);
  //   fs->SetParameters(&par[0]);
  //   fb->SetParameters(&par[4]); //3

  //   fs->Draw("same");
  //   fb->Draw("same");

  //   double N1 = fs->Integral(mkk_min, mkk_max) / hphiy_py->GetBinWidth(1);
  //   double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

  //   if (N1 <= 0)
  //     gphiy->RemovePoint(i - 1);

  //   gphiy->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
  //   gphiy->SetPointError(i - 1, 0, dN1);

  //   // gnophiy->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
  //   // gnophiy->SetPointError(i - 1, 0, dNbkg);
  //   // ofs_ul_yphi2pi << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

  //   TLatex lat_phiy;
  //   lat_phiy.SetTextSize(0.09);
  //   lat_phiy.SetTextAlign(13); //align at top
  //   lat_phiy.SetNDC();
  //   lat_phiy.SetTextColor(kBlue);
  //   lat_phiy.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
  //   lat_phiy.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
  //   lat_phiy.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
  //   lat_phiy.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
  //   lat_phiy.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

  //   // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

  //   hphiy_py->Write();
  //   cgphiy->Update();
  //   // c2->Update();
  //   //sleep(1);
  //   }

  //   // cout << endl;

    cgphiy->cd();
    gphiy->Draw("ALP");
  //   gphiy->Write();
  //   gphiy->Write("gphiy", TObject::kWriteDelete);
  //   gphiy->SetMinimum(-70.);
    // cgnophiy->cd();//j);
    // gnophiy->Draw("AP");

    // TH1 *hphiy = gphiy->GetHistogram()->Draw();

    double wid = gphiy->GetX()[1] - gphiy->GetX()[0];
    TH1F *hgphiy = new TH1F("hgphiy", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", gphiy->GetN(), gphiy->GetX()[0] - 0.5 * wid, gphiy->GetX()[gphiy->GetN() - 1] + 0.5 * wid);
    for (int i = 0; i < gphiy->GetN(); ++i)
    {
      hgphiy->SetBinContent(i + 1, gphiy->GetY()[i]);
      hgphiy->SetBinError(i + 1, gphiy->GetEY()[i]);
    }
    hgphiy->Draw("e");
    hgphiy->SetMinimum(-100.);
    // save the histo to root file
    w.factory(Form("Voigtian::sig_y_data(m_y_data[2, 3],mean_y_data[%f],width_y_data[%f],sigma_y_data[%f])", fsb_mc->GetParameter(1),fsb_mc->GetParameter(3), fsb_mc->GetParameter(2))); //  m2pi2k_min, m2pi2k_max,  //sigma_y_data[0.0001,0.01], mean_y_data[1.011,1.030], m_y_data[1.93,3]
    // w.factory(Form("BreitWigner::sig_y_data(m_y_data[%f,%f],mean_y_data[1.018,1.021],width_y_data[0.004])", mkk_min, mkk_max));
    w.factory("Chebychev::bkg_y_data(m_y_data,{c0_y_data[-100000,+100000], c1_y_data[-100000,+100000], c2_y_data[-100000,+100000]})");//, c2_y_data[-100000,+100000]
    w.factory("SUM:::model_y_data(nsig_y_data[-1000,+100000]*sig_y_data, nbkg_y_data[0,+100000]*bkg_y_data)"); //nsig[0,100000000]*sig2,
    w.var("m_y_data")->SetTitle("#phi#pi^{+}#pi^{-} [GeV/c^{2}]");
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

    // TLegend *l_y_data = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_y_data->SetFillColor(kWhite);
    // l_y_data->AddEntry(fr_y_data->findObject("lsig_y_data"), Form("N_{Sig} = %.2f", w.var("nsig_y_data")->getVal()), "l");
    // l_y_data->AddEntry(fr_y_data->findObject("lbkg_y_data"), Form("N_{Bkg} = %.2f", w.var("nbkg_y_data")->getVal()), "l");
    // l_y_data->Draw();

    double N_y_data = w.var("nsig_y_data")->getVal();
    double dN_y_data = w.var("nsig_y_data")->getError();

    TLatex lat_phiye;
    lat_phiye.SetTextSize(0.04);
    lat_phiye.SetTextAlign(13); //align at top
    lat_phiye.SetNDC();
    lat_phiye.SetTextColor(kBlue);
    // lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phiye.DrawLatex(0.60, 0.40, Form("N_{sig} = %0.2f#pm%0.2f", N_y_data, dN_y_data));
    lat_phiye.DrawLatex(0.60, 0.34, Form("#mu = %0.3f#pm%0.3f", w.var("mean_y_data")->getVal(), w.var("mean_y_data")->getError()));
    lat_phiye.DrawLatex(0.60, 0.28, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_y_data")->getVal(), w.var("sigma_y_data")->getError()));
    lat_phiye.DrawLatex(0.60, 0.22, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_y_data")->getVal(), w.var("width_y_data")->getError()));

    // +++++++++++++++++++ Systematic Errors
    double bkg_poln[3] = {-274, -286, -261}; // Bkg pol degree: 4, 3, 5
    double err_bkg_poln = TMath::RMS(3, bkg_poln);
    double fit_range[3] = {-274, -320, -225}; // range phi2pi:[2, 3], [2, 2.8], [2, 3.2]
    double err_fit_range = TMath::RMS(3, fit_range);
    double slice_numb[3] = {-274, -339, -325}; // K+K- bining: 50, 40, 60
    double err_slice_numb = TMath::RMS(3, slice_numb);
    double mean_value[3] = {-274, -278, -269}; // Y Mean value : 2.194, 2.183, 2.205  
    double err_mean_value = TMath::RMS(3, mean_value);
    double width_value[3] = {-274, -230, -319}; // Y Width value : 0.079, 0.065, 0.093  
    double err_width_value = TMath::RMS(3, width_value);
    
    double err_sys = TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range + err_slice_numb*err_slice_numb + err_mean_value*err_mean_value + err_width_value*err_width_value);

  // ++++++++++++++++++++++++++++ Cross-section  +++++++++++++++++++++++
    double N_y_gen = h_beame_tru->Integral(1,10);
    double eff_y = N_y_mc / N_y_gen; // Efficiency = N_observed/N_generated
    double deff_y = eff_y * (dN_y_mc / N_y_mc);
    double lumi_y = h_tagged_flux->Integral(1,10); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_y = 1e9 * N_y_data / (eff_y * lumi_y * T * br_phi);
    double dxsec_y_stat = xsec_y * (dN_y_data / N_y_data);
    double dxsec_y_sys = xsec_y * TMath::Sqrt((err_sys / N_y_data) * (err_sys / N_y_data) + (deff_y / eff_y) * (deff_y / eff_y) + (dbr_phi / br_phi) * (dbr_phi / br_phi));
    double dxsec_y_tot = TMath::Sqrt(dxsec_y_stat*dxsec_y_stat + dxsec_y_sys*dxsec_y_sys);

    // ++++++++++++++++++++++++++++ Upper Limit  90% CL  +++++++++++++++++

    // ++++++++++ gaussian method
    // double ul = mu + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, sig, mu), sig);
    // double ul = xsec_y + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, dxsec_y_tot, xsec_y), dxsec_y_tot);

    // ++++++++++ profile likelihood method
    double dN_y_tot = sqrt(dN_y_data*dN_y_data + err_sys*err_sys);
    TGraph *gprof = profileLikelihoodScan(dh_y_data, *(w.pdf("model_y_data")), w.var("nsig_y_data"), N_y_data - 5 * dN_y_tot, N_y_data + 5 * dN_y_tot);
    TCanvas *cprof = new TCanvas("cprof","",900,600);
    cprof->cd();
    cprof->SetGrid();
    gprof->SetTitle("; Yield_{#phi#pi^{+}#pi^{-}}; -Log(L)");
    gprof->SetLineWidth(2);
    gprof->SetMarkerStyle(20);
    gprof->SetMarkerSize(1.0);
    gprof->Draw("AP");
    gprof->Fit("pol2", "R", "", gprof->GetHistogram()->GetXaxis()->GetXmin(), gprof->GetHistogram()->GetXaxis()->GetXmax());

    TCanvas *cprofxsec = new TCanvas("cprofxsec","",900,600);
    TGraph *gprofxsec = new TGraph();
    for (int i = 0; i < gprof->GetN(); i++)
    {
      double x_gprof, y_gprof;
      gprof->GetPoint(i,x_gprof,y_gprof);
      double xsec_gprof = 1e9 * x_gprof / (eff_y * lumi_y * T * br_phi);
      cout << "i = " << i << " | x_gprof = "<<x_gprof<<" | y_gprof = " <<y_gprof<<" | xsec_gprof = " <<xsec_gprof<<endl;
      gprofxsec->SetPoint(i, xsec_gprof, exp(-y_gprof));
    }
    cprofxsec->cd();
    gprofxsec->SetTitle("; #sigma [nb]; Likelihood");
    gprofxsec->SetLineWidth(2);
    gprofxsec->SetMarkerStyle(20);
    gprofxsec->SetMarkerSize(1.0);
    gprofxsec->Draw("AP");
    TF1 *func_profxsec = new TF1("func_profxsec", "gaus", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());
    gprofxsec->Fit(func_profxsec, "R", "", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());
    double ul = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, func_profxsec->GetParameter(2), func_profxsec->GetParameter(1)), func_profxsec->GetParameter(2));
    double ul2 = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9,func_profxsec->GetParameter(2));
    TLine *l_mean = new TLine(func_profxsec->GetParameter(1), 0, func_profxsec->GetParameter(1), 1);
    l_mean->SetLineWidth(2);
    l_mean->SetLineStyle(10);
    l_mean->Draw("same");
    TLine *l_sigma = new TLine(func_profxsec->GetParameter(1), 0.5, func_profxsec->GetParameter(1)+func_profxsec->GetParameter(2), 0.5);
    l_sigma->SetLineWidth(2);
    l_sigma->SetLineStyle(10);
    l_sigma->Draw("same");
    TLine *l_ul = new TLine(ul, 0, ul, 0.45);
    l_ul->SetLineWidth(2);
    l_ul->SetLineColor(kBlue);
    l_ul->Draw("same");
    // TLine *l_ul2 = new TLine(ul2, 0, ul2, 0.45);
    // l_ul2->SetLineWidth(2);
    // l_ul2->SetLineColor(28);
    // l_ul2->Draw("same");
    //  +++++++++++++++++++

    fprintf(table_ul_yphi2pi, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.3f & %0.3f $\\pm$ %0.3f $\\pm$ %0.3f & %0.3f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), N_y_gen, N_y_mc, dN_y_mc, N_y_data, dN_y_data, eff_y * 100, deff_y * 100, xsec_y, dxsec_y_stat, dxsec_y_sys, ul);    

    fprintf(table_ul_yphi2pi_sys, "%0.2f - %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), abs(err_bkg_poln*100/bkg_poln[0]), abs(err_fit_range*100/fit_range[0]), abs(err_slice_numb*100/slice_numb[0]), abs(err_mean_value*100/mean_value[0]), abs(err_width_value*100/width_value[0]));    

    // int j =1;
    // gphiy->Write(Form("grphiy_%d", j), TObject::kWriteDelete);

    // cphiy1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c1_phiy_%s.root", name.Data()), "root");
    // cphiy1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c1_phiy_%s.eps", name.Data()), "eps");
    // cphiy1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c1_phiy_%s.png", name.Data()), "png");
    // cphiy2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c2_phiy_%s.root", name.Data()), "root");
    // cphiy2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c2_phiy_%s.eps", name.Data()), "eps");
    // cphiy2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c2_phiy_%s.png", name.Data()), "png");
    // cphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_phiy_%s.root", name.Data()), "root");
    // cphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_phiy_%s.eps", name.Data()), "eps");
    // cphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_phiy_%s.png", name.Data()), "png");
    // cgnophiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gnophiy.root", "root");
    // cgnophiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gnophiy.eps", "eps");

    cgphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gphiy_%s.root", name.Data()), "root");
    cgphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gphiy_%s.eps", name.Data()), "eps");
    cgphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gphiy_%s.png", name.Data()), "png");
    cgphiy->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_gphiy_%s.pdf", name.Data()), "pdf");
    cprof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_prof_%s.root", name.Data()), "root");
    cprof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_prof_%s.eps", name.Data()), "eps");
    cprof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_prof_%s.png", name.Data()), "png");
    cprof->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c_prof_%s.pdf", name.Data()), "pdf");
    cprofxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_profxsec.root", name.Data()), "root");
    cprofxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_profxsec.eps", name.Data()), "eps");
    cprofxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_profxsec.png", name.Data()), "png");
    cprofxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_profxsec.pdf", name.Data()), "pdf");


    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_tagged_flux.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_tagged_flux.png", name.Data()), "png");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_beame_tru.root", name.Data()), "root");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_beame_tru.eps", name.Data()), "eps");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphi2pi/c%s_beame_tru.png", name.Data()), "png");

    outputfig->Print();

    fprintf(table_ul_yphi2pi, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_ul_yphi2pi);
    gSystem->Exec(Form("pdflatex table_ul_yphi2pi_%s.tex",name.Data()));

    fprintf(table_ul_yphi2pi_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_ul_yphi2pi_sys);
    gSystem->Exec(Form("pdflatex table_ul_yphi2pi_sys_%s.tex",name.Data()));
}


    // // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", 1.71, 2.06);
    // // //  TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.98, 1.2);
    // // fsb->SetLineColor(2);
    // // fsb->SetParameters(68, 1.8, 0.09, 1, 1);
    // // fsb->SetParLimits(0, 0, 10000);
    // // fsb->SetParLimits(1, 1.75, 1.85); // 1.018, 1.021
    // // fsb->SetParLimits(2, 0.03, 0.3); //fsb->SetParLimits(2, 0.008, 0.010);
    // //                                  //  fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    // // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", 1.71, 2.06);
    // // //  TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
    // // fs->SetLineColor(4);

    // // TF1 *fb = new TF1("fb", "pol2(3)", 1.71, 2.06); //3
    // // fb->SetLineColor(28);
    // // fb->SetLineStyle(2);

    // // gphiy->Fit("fsb", "", "", 1.71, 2.06);
    // // double par[7]; //6
    // // fsb->GetParameters(&par[0]);
    // // fs->SetParameters(&par[0]);
    // // fb->SetParameters(&par[3]); //3

    // // fs->Draw("same");
    // // fb->Draw("same");

    // // TLatex lat_phiye;
    // // lat_phiye.SetTextSize(0.04);
    // // lat_phiye.SetTextAlign(13); //align at top
    // // lat_phiye.SetNDC();
    // // lat_phiye.SetTextColor(kBlue);
    // // lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    // // lat_phiye.DrawLatex(0.68, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", fs->Integral(1.71, 2.06), fs->Integral(1.71, 2.06) * fsb->GetParError(0) / fsb->GetParameter(0)));
    // // lat_phiye.DrawLatex(0.68, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    // // lat_phiye.DrawLatex(0.68, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
    // //  lat_phiye.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    // // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol4(3)", 1.71, 2.06);
    // TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)",  2, 3);
    // fsb_data->SetLineColor(2);
    // fsb_data->SetParameters(800, 2.194, 0.016, 0.093, 1,1,1,1,1);
    // // fsb_data->SetParLimits(0, 0, 10000);
    // // fsb_data->SetParLimits(1, 2.0, 2.5); // 1.018, 1.021
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1));
    // fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); // 0.016, 
    // fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); // 0.079, 0.065, 0.093

    // // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::BreitWigner(x,[1],[2])", 1.7, 3.2);
    // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", 2, 3);// 2.0, 2.5
    // fs_data->SetLineColor(4);

    // TF1 *fb_data = new TF1("fb_data", "pol4(4)",  2, 3); //3
    // fb_data->SetLineColor(28);
    // fb_data->SetLineStyle(2);

    // gphiy->Fit("fsb_data", "", "",  2, 3);
    // double par_data[10]; //6
    // fsb_data->GetParameters(&par_data[0]);
    // fs_data->SetParameters(&par_data[0]);
    // fb_data->SetParameters(&par_data[4]); //3

    // double N_y_data = fs_data->Integral(2, 3)*n2pi2k/(m2pi2k_max-m2pi2k_min);
    // double dN_y_data = N_y_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    // fs_data->Draw("same");
    // fb_data->Draw("same");

    // TLatex lat_phiye;
    // lat_phiye.SetTextSize(0.04);
    // lat_phiye.SetTextAlign(13); //align at top
    // lat_phiye.SetNDC();
    // lat_phiye.SetTextColor(kBlue);
    // lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    // lat_phiye.DrawLatex(0.68, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_y_data, dN_y_data));
    // lat_phiye.DrawLatex(0.68, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    // lat_phiye.DrawLatex(0.68, 0.67, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    // lat_phiye.DrawLatex(0.68, 0.60, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));