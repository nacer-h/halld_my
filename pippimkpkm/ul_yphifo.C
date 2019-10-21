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

void ul_yphifo(TString name, int nkk=100, int n2pi2k=50, int ne=1, int nt=1) // , TString cut="&& kin_chisq<30 && abs(mm2)<0.015"
{
   TFile *fdata = NULL;
  if(name == "data_16") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_16_flat.root");
  if(name == "data_17") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_17_flat.root");
  if(name == "data_18") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_18_flat.root");
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TFile *fps = NULL;
  if(name == "data_16") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_11366_11555.root");
  if(name == "data_17") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  if(name == "data_18") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_40856_42577.root");
  TFile *fmc = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_phifo_genr8_17v3_flat.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TFile *ftru = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_thrown_phifo_17v3.root");
  TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/ul_yphifo.root","UPDATE");

  FILE *table_ul_yphifo = fopen("table_ul_yphifo.tex","w");
  fprintf(table_ul_yphifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Total cross-sections and upper limit}\n \\begin{tabularx}{\\textwidth}{|c|X|X|X|X|X|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [$\\%$] & $\\sigma$ [nb] & $90\\%$ CL limit [nb] \\\\ \n \\hline\n");

  FILE *table_ul_yphifo_sys = fopen("table_ul_yphifo_sys.tex","w");
  fprintf(table_ul_yphifo_sys,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Systematic errors}\n \\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ $\\%$ ] & Fit Range [$\\%$] & $\\phi$-mass bins [$\\%$]  & Y Mean [$\\%$] & Y width [$\\%$] \\\\ \n \\hline\n");

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
  cout << "h_tagged_flux" << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.0);
  h_tagged_flux->Draw("e");

  // ++++++++++++++++++++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru" << h_beame_tru << endl;
  // h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.0);
  h_beame_tru->Draw("e");  

  // *********************** Phi(1020) MC *********************************
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", "MC signal; m_{K^{+}K^{-}} [GeV/c^{2}]; Counts", 200, 1.005, 1.035);
  tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");//+cutlist+" && kin_chisq<25)"
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");

  // w.factory("BreitWigner::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004])");
  // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.001,0.01])"); //sigma_PhiMass_mc[0.001,0.01], mean_PhiMass_mc[1.015,1.022]
  w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-10,10], c1_PhiMass_mc[-10,10]})");
  w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000000]*sig_PhiMass_mc, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)");//, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); //nsig[0,100000000]*sig2,
  w.var("m_PhiMass_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
  RooDataHist dh_PhiMass_mc("dh_PhiMass_mc", "dh_PhiMass_mc", *w.var("m_PhiMass_mc"), Import(*h_PhiMass_postcuts));
  RooPlot *fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title("K^{+}K^{-}"));
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
  lat_PhiMass_mc.DrawLatex(0.62, 0.87, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_PhiMass_mc")->getVal(), w.var("nsig_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.78, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_PhiMass_mc")->getVal(), w.var("nbkg_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.68, Form("#mu = %0.3f#pm%0.3f", w.var("mean_PhiMass_mc")->getVal(), w.var("mean_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.58, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_PhiMass_mc")->getVal(), w.var("width_PhiMass_mc")->getError()));
  lat_PhiMass_mc.DrawLatex(0.62, 0.48, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_PhiMass_mc")->getVal(), w.var("sigma_PhiMass_mc")->getError()));

  // TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 1.005, 1.035);
  // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  // fsb->SetLineColor(2);
  // fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1, 1, 1);
  // fsb->SetParLimits(0, 0, 10000);
  // fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
  // fsb->FixParameter(2, 0.004);        //fsb->SetParLimits(2, 0.008, 0.010);
  // fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

  // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 1.005, 1.035);
  // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  // fs->SetLineColor(4);

  // TF1 *fb = new TF1("fb", "pol4(4)", 1.005, 1.035); //pol2(3)
  // fb->SetLineColor(28);
  // fb->SetLineStyle(2);

  // h_PhiMass_postcuts->Fit("fsb", "", "", 1.005, 1.035);
  // double par[8]; //6
  // fsb->GetParameters(&par[0]);
  // fs->SetParameters(&par[0]);
  // fb->SetParameters(&par[4]); //3

  // fs->Draw("same");
  // fb->Draw("same");

  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted.root", "root");
  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted.eps", "eps");
  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted.png", "png");

  // *********************** Y(2175) MC *********************************
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts", "c_YMass_postcuts", 900, 600);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_beambunchcut", "(MC);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 2, 2.5); //[2, 2.5]
  tmc->Project("h_YMass_beambunchcut", "kpkmpippim_mf", "w8*(kpkmpippim_uni && is_truecombo && kpkm_mf>1.005 && kpkm_mf<1.035)");
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

  double N_y_mc = fs_mc->Integral(2, 2.5) / h_YMass_postcuts->GetBinWidth(1);
  double dN_y_mc = N_y_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

  TLegend *l_y_mc = new TLegend(0.2, 0.65, 0.35, 0.85);
  l_y_mc->SetTextSize(0.04);
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
  lat_mc.DrawLatex(0.6, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
  lat_mc.DrawLatex(0.6, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N_y_mc, dN_y_mc));
  lat_mc.DrawLatex(0.6, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
  lat_mc.DrawLatex(0.6, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
  lat_mc.DrawLatex(0.6, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));
  // lat_mc.Draw("same");

  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_YMass_postcuts_fitted.root", "root");
  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_YMass_postcuts_fitted.eps", "eps");
  c_YMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/cmc_YMass_postcuts_fitted.png", "png");

  // cout<<"=============================  no problem up to here ! ========================"<<endl;

  // ======================================== Y vs. Phi(1020) ===============================================
  // root -l 'ul_yphifo.C+(100,50,1,1)'

  TCanvas *cphiyall=new TCanvas("cphiyall","cphiyall",900, 600);//1500, 800

  TCanvas *cphiyall1 = new TCanvas("cphiyall1","cphiyall1", 1500, 800);//1500, 800
  cphiyall1->Divide(5, 5);
  TCanvas *cphiyall2 = new TCanvas("cphiyall2", "cphiyall2", 1500, 800);
  cphiyall2->Divide(5, 5);

  TCanvas *cgphiyall=new TCanvas("cgphiyall","cgphiyall",900, 600);//1500, 800
  TGraphErrors *gphiyall;

  // TCanvas *cgnophiyall=new TCanvas("cgnophiyall","cgnophiyall",1500, 800);
  // TGraphErrors *gnophiyall= new TGraphErrors;//(n2pi2k)
  // gnophiyall->SetMarkerStyle(20);
  // gnophiyall->SetMarkerSize(1.0);
  // gnophiyall->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  cphiyall->cd();
  TH2F *h1d2 = new TH2F("h1d2",Form("(%s);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]",name.Data()), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
  tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*(kpkm_uni || kpkmpippim_uni)");

  h1d2->Draw("colz");

  gphiyall = new TGraphErrors(); //n2pi2k
  gphiyall->SetMarkerStyle(20);
  gphiyall->SetMarkerSize(1.0);
  gphiyall->SetMarkerColor(1);
  gphiyall->SetTitle(Form("(%s);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}",name.Data()));

  // gnophiyall->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1, Eg2));

  for (int i = 1; i <= n2pi2k; ++i) //n2pi2k
  {
    // cout << i << " " << flush;

    if (i < 26)
      cphiyall1->cd(i); //i
    if (i > 25 && i < 51)
      cphiyall2->cd(i - 25);
    // if(i > 50 && i<76) cphiyall3->cd(i-50);
    // if(i > 75) cphiyall4->cd(i-75);

    TH1D *hphiyall_py = h1d2->ProjectionY(Form("_hphiyall_py_%d", i), i, i);

    // hphiyall_py->Draw();
    // if (hphiyall_py->Integral(1, 50) < 100)
    //   continue;

    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
    TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 0.98, 1.2);
    fsb->SetLineColor(2);
    fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1,1,1,1,1);
    fsb->SetParLimits(0, 0, 10000);
    fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb->FixParameter(2, w.var("sigma_PhiMass_mc")->getVal());   //fsb->SetParLimits(2, 0.008, 0.010); // 0.0042
    fsb->FixParameter(3, w.var("width_PhiMass_mc")->getVal());        //fsb->SetParLimits(3, 0.001,0.01);// 0.002

    // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
    fs->SetLineColor(4);

    TF1 *fb = new TF1("fb", "pol4(4)", mkk_min, mkk_max); //3
    fb->SetLineColor(28);
    fb->SetLineStyle(2);

    hphiyall_py->Fit("fsb", "", "", mkk_min, mkk_max);
    double par[10]; //6
    fsb->GetParameters(&par[0]);
    fs->SetParameters(&par[0]);
    fb->SetParameters(&par[4]); //3

    fs->Draw("same");
    fb->Draw("same");

    double N1 = fs->Integral(mkk_min, mkk_max) / hphiyall_py->GetBinWidth(1);
    double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

    if (N1 <= 0)
      gphiyall->RemovePoint(i - 1);

    gphiyall->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
    gphiyall->SetPointError(i - 1, 0, dN1);

    // gnophiyall->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
    // gnophiyall->SetPointError(i - 1, 0, dNbkg);
    // ofs_ul_yphifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

    TLatex lat_phiyall;
    lat_phiyall.SetTextSize(0.09);
    lat_phiyall.SetTextAlign(13); //align at top
    lat_phiyall.SetNDC();
    lat_phiyall.SetTextColor(kBlue);
    lat_phiyall.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    lat_phiyall.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
    lat_phiyall.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    lat_phiyall.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
    lat_phiyall.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

    hphiyall_py->Write();
    cgphiyall->Update();
    // c2->Update();
    //sleep(1);
    }

    // cout << endl;

    cgphiyall->cd();
    gphiyall->Draw("AP");
    gphiyall->SetMinimum(-70.);
    // cgnophiyall->cd();//j);
    // gnophiyall->Draw("AP");

    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", 1.71, 2.06);
    // //  TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.98, 1.2);
    // fsb->SetLineColor(2);
    // fsb->SetParameters(68, 1.8, 0.09, 1, 1);
    // fsb->SetParLimits(0, 0, 10000);
    // fsb->SetParLimits(1, 1.75, 1.85); // 1.018, 1.021
    // fsb->SetParLimits(2, 0.03, 0.3); //fsb->SetParLimits(2, 0.008, 0.010);
    //                                  //  fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", 1.71, 2.06);
    // //  TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
    // fs->SetLineColor(4);

    // TF1 *fb = new TF1("fb", "pol2(3)", 1.71, 2.06); //3
    // fb->SetLineColor(28);
    // fb->SetLineStyle(2);

    // gphiyall->Fit("fsb", "", "", 1.71, 2.06);
    // double par[7]; //6
    // fsb->GetParameters(&par[0]);
    // fs->SetParameters(&par[0]);
    // fb->SetParameters(&par[3]); //3

    // fs->Draw("same");
    // fb->Draw("same");

    // TLatex lat_phiye;
    // lat_phiye.SetTextSize(0.04);
    // lat_phiye.SetTextAlign(13); //align at top
    // lat_phiye.SetNDC();
    // lat_phiye.SetTextColor(kBlue);
    // lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    // lat_phiye.DrawLatex(0.68, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", fs->Integral(1.71, 2.06), fs->Integral(1.71, 2.06) * fsb->GetParError(0) / fsb->GetParameter(0)));
    // lat_phiye.DrawLatex(0.68, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    // lat_phiye.DrawLatex(0.68, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
    //  lat_phiye.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol4(3)", 1.71, 2.06);
    TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)",  2, 3);
    fsb_data->SetLineColor(2);
    fsb_data->SetParameters(800, 2.194, 0.016, 0.093, 1,1,1,1,1);
    // fsb_data->SetParLimits(0, 0, 10000);
    // fsb_data->SetParLimits(1, 2.0, 2.5); // 1.018, 1.021
    fsb_data->FixParameter(1, fsb_mc->GetParameter(1));
    fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); // 0.016, 
    fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); // 0.079, 0.065, 0.093

    // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::BreitWigner(x,[1],[2])", 1.7, 3.2);
    TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", 2, 3);// 2.0, 2.5
    fs_data->SetLineColor(4);

    TF1 *fb_data = new TF1("fb_data", "pol4(4)",  2, 3); //3
    fb_data->SetLineColor(28);
    fb_data->SetLineStyle(2);

    gphiyall->Fit("fsb_data", "", "",  2, 3);
    double par_data[10]; //6
    fsb_data->GetParameters(&par_data[0]);
    fs_data->SetParameters(&par_data[0]);
    fb_data->SetParameters(&par_data[4]); //3

    double N_y_data = fs_data->Integral(2, 3)*n2pi2k/(m2pi2k_max-m2pi2k_min);
    double dN_y_data = N_y_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    fs_data->Draw("same");
    fb_data->Draw("same");

    TLatex lat_phiye;
    lat_phiye.SetTextSize(0.04);
    lat_phiye.SetTextAlign(13); //align at top
    lat_phiye.SetNDC();
    lat_phiye.SetTextColor(kBlue);
    lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phiye.DrawLatex(0.68, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_y_data, dN_y_data));
    lat_phiye.DrawLatex(0.68, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phiye.DrawLatex(0.68, 0.67, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phiye.DrawLatex(0.68, 0.60, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

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

    // ++++++++++++++++++++++++++++ Upper Limit  90% CL  +++++++++++++++++
    // double ul = mu + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, sig, mu), sig);
    double dxsec_y_tot = TMath::Sqrt(dxsec_y_stat*dxsec_y_stat + dxsec_y_sys*dxsec_y_sys);
    double ul = xsec_y + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, dxsec_y_tot, xsec_y), dxsec_y_tot);

    fprintf(table_ul_yphifo, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), N_y_gen, N_y_mc, dN_y_mc, N_y_data, dN_y_data, eff_y * 100, deff_y * 100, xsec_y, dxsec_y_stat, dxsec_y_sys, ul);    

    fprintf(table_ul_yphifo_sys, "%0.2f - %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), abs(err_bkg_poln*100/bkg_poln[0]), abs(err_fit_range*100/fit_range[0]), abs(err_slice_numb*100/slice_numb[0]), abs(err_mean_value*100/mean_value[0]), abs(err_width_value*100/width_value[0]));    

    // int j =1;
    // gphiyall->Write(Form("grphiyall_%d", j), TObject::kWriteDelete);

    cphiyall1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiyall_%s.root", name.Data()), "root");
    cphiyall1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiyall_%s.eps", name.Data()), "eps");
    cphiyall1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiyall_%s.png", name.Data()), "png");
    cphiyall2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiyall_%s.root", name.Data()), "root");
    cphiyall2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiyall_%s.eps", name.Data()), "eps");
    cphiyall2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiyall_%s.png", name.Data()), "png");
    cphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyall_%s.root", name.Data()), "root");
    cphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyall_%s.eps", name.Data()), "eps");
    cphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyall_%s.png", name.Data()), "png");
    cgphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyall_%s.root", name.Data()), "root");
    cgphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyall_%s.eps", name.Data()), "eps");
    cgphiyall->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyall_%s.png", name.Data()), "png");
    // cgnophiyall->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gnophiyall.root", "root");
    // cgnophiyall->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gnophiyall.eps", "eps");

/*
  // ======================================== Y vs. Eg ===============================================
  // root -l 'ul_yphifo.C+(50,50,4,4)'
  // Eg= 3, 7.9, 8.56, 9.66, 12
  TCanvas *cphiye=new TCanvas("cphiye","cphiye",1500, 800);//1500, 800
  cphiye->Divide(3,2);

  TCanvas *cphiye1[ne];
  TCanvas *cphiye2[ne];
  // TCanvas *cphiye3[ne];
  // TCanvas *cphiye4[ne];

  TCanvas *cgphiye=new TCanvas("cgphiye","cgphiye",1500, 800);//1500, 800
  cgphiye->Divide(3,2);
  TGraphErrors *gphiye[ne];

  // TCanvas *cgnophiye=new TCanvas("cgnophiye","cgnophiye",1500, 800);
  // TGraphErrors *gnophiye= new TGraphErrors;//(n2pi2k)
  // gnophiye->SetMarkerStyle(20);
  // gnophiye->SetMarkerSize(1.0);
  // gnophiye->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  double Egmin = 3; //6.3 = hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 12; //11.6 = hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax-Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    cphiye->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<E_{#gamma}<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // TH2F *h1d2 = new TH2F("h1d2","(Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    // 4: 6, 8.1, 8.7, 9.3, 12
    // 3: 6, 8.3, 9.6, 12
    TH2F *h1d2;

    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "3<E_{#gamma}<7.6 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<7.6)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "7.6<E_{#gamma}<8.2 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>7.6 && beam_p4_kin.E()<8.2)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "8.2<E_{#gamma}<8.5 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.2 && beam_p4_kin.E()<8.5)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "8.5<E_{#gamma}<8.85 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.5 && beam_p4_kin.E()<8.85)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "8.85<E_{#gamma}<10 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.85 && beam_p4_kin.E()<10)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "10<E_{#gamma}<12 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>10 && beam_p4_kin.E()<12)");
    }

    h1d2->Draw("colz");

    cphiye1[j] = new TCanvas(Form("cphiye1_%d", j), Form("cphiye1_%d", j), 1500, 800);//1500, 800
    cphiye1[j]->Divide(5, 5);
    cphiye2[j] = new TCanvas(Form("cphiye2_%d", j), Form("cphiye2_%d", j), 1500, 800);
    cphiye2[j]->Divide(5, 5);
    // cphiye3[j] = new TCanvas(Form("cphiye3_%d", j), Form("cphiye3_%d", j), 1500, 800);
    // cphiye3[j]->Divide(5, 5);
    // cphiye4[j] = new TCanvas(Form("cphiye4_%d", j), Form("cphiye4_%d", j), 1500, 800);
    // cphiye4[j]->Divide(5, 5);

    gphiye[j] = new TGraphErrors();//n2pi2k
    gphiye[j]->SetMarkerStyle(20);
    gphiye[j]->SetMarkerSize(1.0);
    gphiye[j]->SetMarkerColor(1);
    // gphiye[j]->SetTitle("(Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    // gphiye[j]->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", Eg1[j], Eg2[j]));
    if(j==1) gphiye[j]->SetTitle("3<E_{#gamma}<7.6 (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiye[j]->SetTitle("7.6<E_{#gamma}<8.2 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiye[j]->SetTitle("8.2<E_{#gamma}<8.5 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiye[j]->SetTitle("8.5<E_{#gamma}<8.85 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiye[j]->SetTitle("8.85<E_{#gamma}<10 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiye[j]->SetTitle("10<E_{#gamma}<12 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    // gnophiye->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1[j], Eg2[j]));


    for (int i = 1; i <=n2pi2k; ++i)//n2pi2k
    {
      // cout << i << " " << flush;

      if(i < 26) cphiye1[j]->cd(i);//i
      if(i > 25 && i<51) cphiye2[j]->cd(i-25);
      // if(i > 50 && i<76) cphiye3[j]->cd(i-50);
      // if(i > 75) cphiye4[j]->cd(i-75);

      TH1D *hphiye_py = h1d2->ProjectionY(Form("_hphiye_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiye_py->Draw();
      if (hphiye_py->Integral(1, 50) < 100)
        continue;
      // fb2->Draw("same");

      // // w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.022],width_phiye[0.004],sigma_phiye[0.00001,0.01])", mkk_min, mkk_max)); //1.018,1.022
      // w.factory(Form("BreitWigner::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.021],width_phiye[0.008, 0.010])", mkk_min, mkk_max));
      // // w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.022],width_phiye[0.004],sigma_phiye[0.0001,0.01])", mkk_min, mkk_max)); //50:mean_phiye[1.01,1.025],sigma_phiye[0.0001,0.01],100:mean_phiye[1.018,1.022]
      // w.factory("Chebychev::bkg_phiye(m_phiye,{c0_phiye[-10,10], c1_phiye[-10,10], c2_phiye[-10,10]})");//
      // w.factory("SUM:::model_phiye(nsig_phiye[0,100000000]*sig_phiye, nbkg_phiye[0,100000000]*bkg_phiye)"); //nsig[0,100000000]*sig2,
      // w.var("m_phiye")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phiye("dh_phiye", "dh_phiye", *w.var("m_phiye"), Import(*hphiye_py));
      // RooPlot *fr_phiye = w.var("m_phiye")->frame(Title(" "));
      // // fr_phiye->SetTitleOffset(0.90, "X");
      // // fr_phiye->SetTitleSize(0.06, "XYZ");
      // // fr_phiye->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phiye")->fitTo(dh_phiye);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phiye.plotOn(fr_phiye, RooFit::Name("ldh_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("sig_phiye")), LineColor(kRed), RooFit::Name("lsig_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("bkg_phiye")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, RooFit::Name("lmodel_phiye"));
      // // w.pdf("model_phiye")->paramOn(fr_phiye, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiye"), *w.var("nbkg_phiye")))); //,*w.var("mean_phiye"),*w.var("width_phiye"),*w.var("sigma_phiye"))));
      // fr_phiye->Draw();

      // TLegend *l_phiye = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phiye->SetFillColor(kWhite);
      // l_phiye->SetLineColor(kWhite);
      // // l_phiye->AddEntry(fr_phiye->findObject("ldh_phiye"), "Data", "p");
      // // l_phiye->AddEntry(fr_phiye->findObject("lmodel_phiye"), "total", "l");
      // l_phiye->AddEntry(fr_phiye->findObject("lsig_phiye"), Form("N_{Sig} = %.2f", w.var("nsig_phiye")->getVal()), "l");
      // l_phiye->AddEntry(fr_phiye->findObject("lbkg_phiye"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiye")->getVal()), "l");
      // l_phiye->Draw();

      // // double N1  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // // double dN1 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      // double N1 = w.var("nsig_phiye")->getVal();
      // double dN1 = w.var("nsig_phiye")->getError();

      // // // Integrate normalized pdf over subrange
      // // w.var("m_phiye")->setRange("sigregion",1.005,1.035);
      // // // RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
      // // RooAbsReal* fbkg = w.pdf("bkg_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      // // RooAbsReal* fsig = w.pdf("sig_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      // // // double Nbkg = -fbkg->getVal()*(w.var("nsig_phiye")->getVal()+w.var("nbkg_phiye")->getVal())+fsig->getVal()*w.var("nsig_phiye")->getVal();
      // // double Nbkg = fbkg->getVal()*w.var("nbkg_phiye")->getVal();
      // // double dNbkg = 0;//fbkg->getError(); // w.var("nbkg_phiye")->getError();

      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.98, 1.2);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//3
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiye_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiye_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiye[j]->RemovePoint(i-1);

      gphiye[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiye[j]->SetPointError(i - 1, 0, dN1);

      // gnophiye->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophiye->SetPointError(i - 1, 0, dNbkg);
      // ofs_ul_yphifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

      TLatex lat_phiye;
      lat_phiye.SetTextSize(0.09);
      lat_phiye.SetTextAlign(13); //align at top
      lat_phiye.SetNDC();
      lat_phiye.SetTextColor(kBlue);
      lat_phiye.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phiye.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiye.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phiye.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phiye.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

      hphiye_py->Write();
      cgphiye->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphiye->cd(j);
    gphiye[j]->Draw("AP");

    // cgnophiye->cd();//j);
    // gnophiye->Draw("AP");

    // int j =1;
    // gphiye->Write(Form("grphiye_%d", j), TObject::kWriteDelete);

    cphiye1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiye_%d.root", j), "root");
    cphiye1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiye_%d.eps", j), "eps");
    cphiye1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiye_%d.png", j), "png");
    cphiye2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiye_%d.root", j), "root");
    cphiye2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiye_%d.eps", j), "eps");
    cphiye2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiye_%d.png", j), "png");
    // cphiye3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c3_phiye_%d.root", j), "root");
    // cphiye3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c3_phiye_%d.eps", j), "eps");
    // cphiye4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c4_phiye_%d.root", j), "root");
    // cphiye4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c4_phiye_%d.eps", j), "eps");

  }

  cphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiye.root", "root");
  cphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiye.eps", "eps");
  cphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiye.png", "png");
  cgphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiye.root", "root");
  cgphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiye.eps", "eps");
  cgphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiye.png", "png");
  // cgnophiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gnophiye.root", "root");
  // cgnophiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gnophiye.eps", "eps");
*/
    /*
  // ======================================== Y vs. -t ===============================================

  TCanvas *cphiyt=new TCanvas("cphiyt","cphiyt", 10, 10, 1500, 800);
  cphiyt->Divide(3,2);

  TCanvas *cphiyt1[nt];
  TCanvas *cphiyt2[nt];

  TCanvas *cgphiyt=new TCanvas("cgphiyt","cgphiyt", 10, 10, 1500, 800);
  cgphiyt->Divide(3,2);
  TGraphErrors *gphiyt[nt];

  double tmin = 0; // 0.1 = hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double tmax = 10;   // 2 = hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double tstep = (tmax-tmin) / nt;
  double t1[nt];
  double t2[nt];

  for (int j = 1; j <= nt; ++j)
  {
    t1[j] = tmin + ((j - 1) * tstep); //i * tstep;
    t2[j] = tmin + (j * tstep);
    cout << "########  j = " << j << " | t1[" << j << "] = " << t1[j] << " | t2[" << j << "] = " << t2[j] << endl;

    cphiyt->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<-t<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    TH2F *h1d2;
    // 3: 0, 0.7, 1.5, 4
    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "0<-t<0.45 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0 && -t_kin<0.45)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "0.45<-t<0.73 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0.45 && -t_kin<0.73)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "0.73<-t<1.1 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0.73 && -t_kin<1.1)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "1.1<-t<1.63 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>1.1 && -t_kin<1.63)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "1.63<-t<2.63 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>1.63 && -t_kin<2.63)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "2.63<-t<10 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>2.63 && -t_kin<10)");
    }

    h1d2->Draw("colz");

    cphiyt1[j] = new TCanvas(Form("cphiyt1_%d", j), Form("cphiyt1_%d", j), 1500, 800);
    cphiyt1[j]->Divide(5, 5);
    cphiyt2[j] = new TCanvas(Form("cphiyt2_%d", j), Form("cphiyt2_%d", j), 1500, 800);
    cphiyt2[j]->Divide(5, 5);

    gphiyt[j] = new TGraphErrors();//n2pi2k
    gphiyt[j]->SetMarkerStyle(20);
    gphiyt[j]->SetMarkerSize(1.0);
    gphiyt[j]->SetMarkerColor(1);
    // gphiyt[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphiyt[j]->SetTitle("0<-t<0.45 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiyt[j]->SetTitle("0.45<-t<0.73 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiyt[j]->SetTitle("0.73<-t<1.1 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiyt[j]->SetTitle("1.1<-t<1.63 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiyt[j]->SetTitle("1.63<-t<2.63 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiyt[j]->SetTitle("2.63<-t<10 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi2k; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphiyt1[j]->cd(i);
      if(i > 25 && i<51) cphiyt2[j]->cd(i-25);

      TH1D *hphiyt_py = h1d2->ProjectionY(Form("_hphiyt_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiyt_py->Draw();
      if (hphiyt_py->Integral(1, 50) < 100)
        continue;
      // fb2->Draw("same");

      // w.factory(Form("Voigtian::sig_phiyt(m_phiyt[%f,%f],mean_phiyt[1.018,1.022],width_phiyt[0.004],sigma_phiyt[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phiy[0.001,0.01], mean_phiy[1.016,1.022]
      // w.factory("Chebychev::bkg_phiyt(m_phiyt,{c0_phiyt[-1,1], c1_phiyt[-1,1]})");//, c2_phiyt[-1,1]
      // w.factory("SUM:::model_phiyt(nsig_phiyt[0,100000000]*sig_phiyt, nbkg_phiyt[0,100000000]*bkg_phiyt)"); //nsig[0,100000000]*sig2,
      // w.var("m_phiyt")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phiyt("dh_phiyt", "dh_phiyt", *w.var("m_phiyt"), Import(*hphiyt_py));
      // RooPlot *fr_phiyt = w.var("m_phiyt")->frame(Title(" "));
      // // fr_phiyt->SetTitleOffset(0.90, "X");
      // // fr_phiyt->SetTitleSize(0.06, "XYZ");
      // // fr_phiyt->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phiyt")->fitTo(dh_phiyt);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phiyt.plotOn(fr_phiyt, RooFit::Name("ldh_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("sig_phiyt")), LineColor(kRed), RooFit::Name("lsig_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("bkg_phiyt")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, RooFit::Name("lmodel_phiyt"));
      // // w.pdf("model_phiyt")->paramOn(fr_phiyt, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiyt"), *w.var("nbkg_phiyt")))); //,*w.var("mean_phiyt"),*w.var("width_phiyt"),*w.var("sigma_phiyt"))));
      // fr_phiyt->Draw();

      // TLegend *l_phiyt = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phiyt->SetFillColor(kWhite);
      // l_phiyt->SetLineColor(kWhite);
      // // l_phiyt->AddEntry(fr_phiyt->findObject("ldh_phiyt"), "Data", "p");
      // // l_phiyt->AddEntry(fr_phiyt->findObject("lmodel_phiyt"), "total", "l");
      // l_phiyt->AddEntry(fr_phiyt->findObject("lsig_phiyt"), Form("N_{Sig} = %.2f", w.var("nsig_phiyt")->getVal()), "l");
      // l_phiyt->AddEntry(fr_phiyt->findObject("lbkg_phiyt"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiyt")->getVal()), "l");
      // l_phiyt->Draw();

      // // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      // double N2 = w.var("nsig_phiyt")->getVal();
      // double dN2 = w.var("nsig_phiyt")->getError();

      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//4
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiyt_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiyt_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiyt[j]->RemovePoint(i-1);

      gphiyt[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiyt[j]->SetPointError(i - 1, 0, dN1);

      TLatex lat_phiyt;
      lat_phiyt.SetTextSize(0.09);
      lat_phiyt.SetTextAlign(13); //align at top
      lat_phiyt.SetNDC();
      lat_phiyt.SetTextColor(kBlue);
      lat_phiyt.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phiyt.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiyt.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phiyt.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phiyt.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      //cgphiyt->Update();
      // c2->Update();
      //sleep(1);


    }

    // cout << endl;

    cgphiyt->cd(j);
    gphiyt[j]->Draw("AP");

    // int j =1;
    // gphiyt->Write(Form("grphiyt_%d", j), TObject::kWriteDelete);

    cphiyt1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiyt_%d.root", j), "root");
    cphiyt1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiyt_%d.eps", j), "eps");
    cphiyt2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiyt_%d.root", j), "root");
    cphiyt2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiyt_%d.eps", j), "eps");
  }

  cphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyt.root", "root");
  cphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyt.eps", "eps");
  cphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiyt.png", "png");
  cgphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyt.root", "root");
  cgphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyt.eps", "eps");
  cgphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiyt.png", "png");
*/
    /*
  // ======================================== Y vs. (Eg,-t) ===============================================

  TCanvas *cphiy = new TCanvas("cphiy", "cphiy", 10, 10, 1500, 800);
  cphiy->Divide(3,2);

  TCanvas *cphiy1[ne*nt];
  TCanvas *cphiy2[ne*nt];

  TCanvas *cgphiy=new TCanvas("cgphiy","cgphiy", 10, 10, 1500, 800);
  cgphiy->Divide(3,2);
  TGraphErrors *gphiy[ne*nt];


  for (int j = 1; j <= ne*nt; ++j)
  {
    cout << "########  j = " << j <<endl;

    cphiy->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<-t<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    TH2F *h1d2;
    // 3: 0, 0.7, 1.5, 4
    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>1.13 && -t_kin<10)");
    }

    h1d2->Draw("colz");

    cphiy1[j] = new TCanvas(Form("cphiy1_%d", j), Form("cphiy1_%d", j), 1500, 800);
    cphiy1[j]->Divide(5, 5);
    cphiy2[j] = new TCanvas(Form("cphiy2_%d", j), Form("cphiy2_%d", j), 1500, 800);
    cphiy2[j]->Divide(5, 5);

    gphiy[j] = new TGraphErrors();//n2pi2k
    gphiy[j]->SetMarkerStyle(20);
    gphiy[j]->SetMarkerSize(1.0);
    gphiy[j]->SetMarkerColor(1);
    // gphiy[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphiy[j]->SetTitle("(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiy[j]->SetTitle("(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiy[j]->SetTitle("(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiy[j]->SetTitle("(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiy[j]->SetTitle("(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiy[j]->SetTitle("(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi2k; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphiy1[j]->cd(i);
      if(i > 25 && i<51) cphiy2[j]->cd(i-25);

      TH1D *hphiy_py = h1d2->ProjectionY(Form("_hphiy_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiy_py->Draw();
      if(hphiy_py->Integral(1,100)<100)
        continue;
      // fb2->Draw("same");

      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 100000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//4
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiy_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiy_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiy[j]->RemovePoint(i-1);

      gphiy[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiy[j]->SetPointError(i - 1, 0, dN1);

      TLatex lat_phiy;
      lat_phiy.SetTextSize(0.09);
      lat_phiy.SetTextAlign(13); //align at top
      lat_phiy.SetNDC();
      lat_phiy.SetTextColor(kBlue);
      lat_phiy.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
      lat_phiy.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiy.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
      lat_phiy.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
      lat_phiy.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    }

    // cout << endl;

    cgphiy->cd(j);
    gphiy[j]->Draw("AP");

    // int j =1;
    // gphiy->Write(Form("grphiy_%d", j), TObject::kWriteDelete);

    cphiy1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiy_%d.root", j), "root");
    cphiy1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiy_%d.eps", j), "eps");
    cphiy1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c1_phiy_%d.png", j), "png");
    cphiy2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiy_%d.root", j), "root");
    cphiy2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiy_%d.eps", j), "eps");
    cphiy2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c2_phiy_%d.png", j), "png");
  }

  cphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiy.root", "root");
  cphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiy.eps", "eps");
  cphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_phiy.png", "png");
  cgphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiy.root", "root");
  cgphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiy.eps", "eps");
  cgphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_gphiy.png", "png");
*/

  c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c%s_tagged_flux.root", name.Data()), "root");
  c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c%s_tagged_flux.eps", name.Data()), "eps");
  c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c%s_tagged_flux.png", name.Data()), "png");
  c_beame_tru->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_beame_tru.root", "root");
  c_beame_tru->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_beame_tru.eps", "eps");
  c_beame_tru->Print("/data.local/nacer/halld_my/pippimkpkm/fig_ul_yphifo/c_beame_tru.png", "png");

  outputfig->Print();

  fprintf(table_ul_yphifo,"\\hline\n \\end{tabularx}\n \\end{center}\n \\end{minipage}\n \\end{table}\n \\end{document}\n");
  fclose(table_ul_yphifo);
  gSystem->Exec("pdflatex table_ul_yphifo.tex");
 
  fprintf(table_ul_yphifo_sys,"\\hline\n \\end{tabularx}\n \\end{center}\n \\end{minipage}\n \\end{table}\n \\end{document}\n");
  fclose(table_ul_yphifo_sys);
  gSystem->Exec("pdflatex table_ul_yphifo_sys.tex");
 
}
