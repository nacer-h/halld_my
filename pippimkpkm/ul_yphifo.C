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
#include "TFrame.h"
#include "TLegendEntry.h"
#include <fstream>
#include "TStyle.h"
#include <TMultiGraph.h>
#include <iostream>
#include "TSystem.h"
#include "TLine.h"
#include "scanphi.C"
#include "showpull.C"
#include "profileLikelihoodScan.C"
#include "TROOT.h"


#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooStats/NumberCountingUtils.h"
#endif
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooHist.h"
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

void ul_yphifo(TString name, int n2k=100, int n2pi2k=50, int ne=1, int nt=1) // , TString cut="&& kin_chisq<30 && abs(mm2)<0.015"
{
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", name.Data()));
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_yphifo_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/thrown_yphifo_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_11366_11555.root");
  if(name == "17") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_30274_31057.root");
  if(name == "18") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_40856_42577.root");
 if(name == "18l") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_50677_51768.root");

  TFile *outputfig = new TFile("output/fig_ul_yphifo/ul_yphifo.root","UPDATE");
  if(outputfig->IsOpen()) printf("File opened successfully\n");

  FILE *table_ul_yphifo = fopen(Form("table_ul_yphifo_%s.tex", name.Data()),"w");
  fprintf(table_ul_yphifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [$\\%$] & $\\sigma$ [nb] & Significance \\\\ \n \\hline\n");

  // FILE *table_ul_yphifo_sys = fopen(Form("table_ul_yphifo_sys_%s.tex", name.Data()),"w");
  // fprintf(table_ul_yphifo_sys,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors}\n \\begin{tabular}{|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ $\\%$ ] & Fit Range [$\\%$] & $\\phi$-mass bins [$\\%$]  & Y Mean [$\\%$] & Y width [$\\%$] \\\\ \n \\hline\n");

  ofstream ofs_yphifo("yphifo_tab.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.99, mkk_max = 1.2; //0.99
  double m2pi_min = 0.85, m2pi_max = 1.06; // [0.85, 1.05]
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2; //  [2.0, 3.0]
  double m2pi2k_min2 = 2.0,  m2pi2k_max2 = 3.0; //  [2.0, 3.0]

  // gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // ++++++++++++++++++++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux" << h_tagged_flux << endl;
  // h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.0);
  h_tagged_flux->Draw("e");
  // h_tagged_flux->Write();
  // ++++++++++++++++++++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  // h_beame_tru->Rebin(100);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");
  // h_beame_tru->Write();

  // outputfig->Print();

  // *********************** Phi(1020) MC *********************************
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 1.0, 1.04);
  tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni && is_truecombo)");//+cutlist+" && kin_chisq<25)"
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");
  // w.factory("BreitWigner::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004])");
  // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  w.factory("Voigtian::sig_phi_mc(m_phi_mc[1.0, 1.04],mean_phi_mc[1.016,1.022],width_phi_mc[0.004],sigma_phi_mc[0.0001,0.1])"); //m_phi_mc[1.005, 1.035],mean_phi_mc[1.018,1.021],width_phi_mc[0.004],sigma_phi_mc[0.0001,0.01]
  w.factory("Chebychev::bkg_phi_mc(m_phi_mc,{c0_phi_mc[-100000,100000], c1_phi_mc[-100000,100000], c2_phi_mc[-100000,100000]})");
  w.factory("SUM:::model_phi_mc(nsig_phi_mc[0,100000]*sig_phi_mc, nbkg_phi_mc[0,100000]*bkg_phi_mc)");//, nbkg_phi_mc[0,100000000]*bkg_phi_mc)"); //nsig[0,100000000]*sig2,
  w.var("m_phi_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
  RooDataHist dh_phi_mc("dh_phi_mc", "dh_phi_mc", *w.var("m_phi_mc"), Import(*h_PhiMass_postcuts));
  RooPlot *fr_phi_mc = w.var("m_phi_mc")->frame(Title(" "));
  // fr_phi_mc->SetTitleOffset(0.90, "X");
  // fr_phi_mc->SetTitleSize(0.06, "XYZ");
  // fr_phi_mc->SetLabelSize(0.06, "xyz");
  w.pdf("model_phi_mc")->fitTo(dh_phi_mc);

  // //result = w.pdf("model")->fitTo(dh_phi,Extended(kTRUE),Save());
  dh_phi_mc.plotOn(fr_phi_mc, RooFit::Name("ldh_phi_mc"));
  w.pdf("model_phi_mc")->plotOn(fr_phi_mc, Components(*w.pdf("sig_phi_mc")), LineColor(kRed), RooFit::Name("lsig_phi_mc"));
  w.pdf("model_phi_mc")->plotOn(fr_phi_mc, Components(*w.pdf("bkg_phi_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phi_mc"));
  w.pdf("model_phi_mc")->plotOn(fr_phi_mc, RooFit::Name("lmodel_phi_mc"));
  // w.pdf("model_phi_mc")->paramOn(fr_phi_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_phi_mc"), *w.var("nbkg_phi_mc")))); //,*w.var("mean_phi_mc"),*w.var("width_phi_mc"),*w.var("sigma_phi_mc"))));
  fr_phi_mc->Draw();

  TLegend *l_phi_mc = new TLegend(0.2, 0.65, 0.4, 0.85);
  l_phi_mc->SetFillColor(kWhite);
  l_phi_mc->SetLineColor(kWhite);
  // l_phi_mc->AddEntry(fr_phi_mc->findObject("ldh_phi_mc"), "Data", "p");
  l_phi_mc->AddEntry(fr_phi_mc->findObject("lmodel_phi_mc"), "total", "l");
  l_phi_mc->AddEntry(fr_phi_mc->findObject("lsig_phi_mc"), "Voigtian", "l");
  l_phi_mc->AddEntry(fr_phi_mc->findObject("lbkg_phi_mc"), "pol 2nd", "l");
  l_phi_mc->Draw();

  TLatex lat_phi_mc;
  lat_phi_mc.SetTextSize(0.05);
  lat_phi_mc.SetTextAlign(13); //align at top
  lat_phi_mc.SetNDC();
  lat_phi_mc.SetTextColor(kBlue);
  lat_phi_mc.DrawLatex(0.62, 0.86, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_phi_mc")->getVal(), w.var("nsig_phi_mc")->getError()));
  lat_phi_mc.DrawLatex(0.62, 0.80, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_phi_mc")->getVal(), w.var("nbkg_phi_mc")->getError()));
  lat_phi_mc.DrawLatex(0.62, 0.74, Form("#mu = %0.3f#pm%0.3f", w.var("mean_phi_mc")->getVal(), w.var("mean_phi_mc")->getError()));
  lat_phi_mc.DrawLatex(0.62, 0.68, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_phi_mc")->getVal(), w.var("width_phi_mc")->getError()));
  lat_phi_mc.DrawLatex(0.62, 0.62, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_phi_mc")->getVal(), w.var("sigma_phi_mc")->getError()));

  // TF1 *fsb_phi_mc = new TF1("fsb_phi_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol3(4)", 1.005, 1.035);
  // // TF1 *fsb_phi_mc = new TF1("fsb_phi_mc", "gaus(0) + pol3(3)", 1.005, 1.035);
  // // TF1 *fsb_phi_mc = new TF1("fsb_phi_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  // fsb_phi_mc->SetLineColor(2);
  // fsb_phi_mc->SetParameters(4000, 1.019, 0.004, 0.003, 1,1,1,1);
  // // fsb_phi_mc->SetParLimits(0, 0, 10000);
  // fsb_phi_mc->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
  // fsb_phi_mc->SetParLimits(2, 0.008, 0.010);
  // fsb_phi_mc->SetParLimits(3, 0.001,0.01);// 0.001,0.01

  // TF1 *fs_phi_mc = new TF1("fs_phi_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", 1.005, 1.035);
  // // TF1 *fs_phi_mc = new TF1("fs_phi_mc", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  // fs_phi_mc->SetLineColor(4);

  // TF1 *fb_phi_mc = new TF1("fb_phi_mc", "pol3(4)", 1.005, 1.035); //pol2(3)
  // fb_phi_mc->SetLineColor(28);
  // fb_phi_mc->SetLineStyle(2);

  // h_phi_postcuts->Fit("fsb_phi_mc", "", "", 1.005, 1.035);
  // double par_phi_mc[10]; //6
  // fsb_phi_mc->GetParameters(&par_phi_mc[0]);
  // fs_phi_mc->SetParameters(&par_phi_mc[0]);
  // fb_phi_mc->SetParameters(&par_phi_mc[4]); //3

  // fs_phi_mc->Draw("same");
  // fb_phi_mc->Draw("same");

  c_PhiMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted_%s.root",name.Data()), "root");
  c_PhiMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted_%s.eps",name.Data()), "eps");
  c_PhiMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_PhiMass_postcuts_fitted_%s.png",name.Data()), "png");

  // *********************** Y(2175) MC *********************************
  TCanvas *c_YMass_postcuts = new TCanvas("c_YMass_postcuts", "c_YMass_postcuts", 900, 600);
  TH1F *h_YMass_postcuts = new TH1F("h_YMass_postcuts", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 2, 2.5); //[2, 2.5]
  tmc->Project("h_YMass_postcuts", "kpkmpippim_mf", "w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && pippim_mf>0.9 && pippim_mf<1.0)"); // is_truecombo && 
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

  TF1 *fsb_y_mc = new TF1("fsb_y_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 2, 2.5); // + pol2(4) , [2, 2.5]
  // TF1 *fsb_y_mc = new TF1("fsb_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
  fsb_y_mc->SetLineColor(2);
  fsb_y_mc->SetParameters(700, 2.175, 0.004, 0.002, 1, 1,1,1,1);//, 1, 1,1
  fsb_y_mc->SetParLimits(0, 0, 10000);
  fsb_y_mc->SetParLimits(1, 2.0, 2.2); // 1.018, 1.021
  fsb_y_mc->SetParLimits(2, 0.01, 0.1);
  // fsb_y_mc->FixParameter(2, 0.079);
  // fsb_y_mc->SetParLimits(3, 0.079-0.14, 0.079+0.14);
  fsb_y_mc->FixParameter(3, 0.079);

  // TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
  TF1 *fs_y_mc = new TF1("fs_y_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", 2, 2.5);
  fs_y_mc->SetLineColor(4);

  TF1 *fb_y_mc = new TF1("fb_y_mc", "pol4(4)", 2, 2.5); //pol2(3)
  fb_y_mc->SetLineColor(28);
  fb_y_mc->SetLineStyle(2);

  h_YMass_postcuts->Fit("fsb_y_mc", "", "", 2, 2.5);
  double par_y_mc[8]; //6
  fsb_y_mc->GetParameters(&par_y_mc[0]);
  fs_y_mc->SetParameters(&par_y_mc[0]);
  fb_y_mc->SetParameters(&par_y_mc[4]); //3

  fs_y_mc->Draw("same");
  fb_y_mc->Draw("same");
  fsb_y_mc->Draw("same");

  double N_y_mc = fs_y_mc->Integral(2, 2.5) / h_YMass_postcuts->GetBinWidth(1);
  double dN_y_mc = N_y_mc * fsb_y_mc->GetParError(0) / fsb_y_mc->GetParameter(0);

  TLegend *l_y_mc = new TLegend(0.2, 0.65, 0.35, 0.85);
  l_y_mc->SetTextSize(0.05);
  l_y_mc->SetFillColor(kWhite);
  l_y_mc->SetLineColor(kWhite);
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  l_y_mc->AddEntry("fsb_y_mc", "total", "l");
  l_y_mc->AddEntry("fs_y_mc", "Voigtian", "l");
  l_y_mc->AddEntry("fb_y_mc", "pol 4^{th}", "l");
  l_y_mc->Draw();

  TLatex lat_mc;
  lat_mc.SetTextSize(0.05);
  lat_mc.SetTextAlign(13); //align at top
  lat_mc.SetNDC();
  lat_mc.SetTextColor(kBlue);
  lat_mc.DrawLatex(0.6, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_y_mc->GetChisquare() / fsb_y_mc->GetNDF()));
  lat_mc.DrawLatex(0.6, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_y_mc, dN_y_mc));
  lat_mc.DrawLatex(0.6, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_y_mc->GetParameter(1), fsb_y_mc->GetParError(1)));
  lat_mc.DrawLatex(0.6, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_y_mc->GetParameter(2), fsb_y_mc->GetParError(2)));
  lat_mc.DrawLatex(0.6, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_y_mc->GetParameter(3), fsb_y_mc->GetParError(3)));
  // lat_mc.Draw("same");

  c_YMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_YMass_postcuts_fitted_%s.root",name.Data()), "root");
  c_YMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_YMass_postcuts_fitted_%s.eps",name.Data()), "eps");
  c_YMass_postcuts->Print(Form("output/fig_ul_yphifo/cmc_YMass_postcuts_fitted_%s.png",name.Data()), "png");

  // cout<<"##############  no problem up to here ! ##################"<<endl;

  // ================================== fo vs. Phi(1020) ===========================
  TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 900, 600); //1500, 800
  TGraphErrors *gphifo = scanphi("fo", "yphifo", Form("%s",name.Data()), n2k, n2pi2k);
  cgphifo->cd();
  gphifo->Draw("AP");

  TF1 *fsb_fo_data = new TF1("fsb_fo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", m2pi_min, m2pi_max);
  fsb_fo_data->SetLineColor(2);
  fsb_fo_data->SetParameters(30, 0.980, 0.01, 1, 1, 1, 1);//30, 0.980, 0.01, 1, 1, 1, 1
  // fsb->SetParLimits(0, 0, 10000);
  fsb_fo_data->SetParLimits(1, 0.95, 0.99);//1, 0.95, 0.99
  fsb_fo_data->SetParLimits(2, 0.01, 0.1);//2, 0.01, 0.1
  // fsb_fo_data->FixParameter(1, fsb_fo_mc->GetParameter(1));
  // fsb_fo_data->FixParameter(2, fsb_fo_mc->GetParameter(2));

  // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
  TF1 *fs_fo_data = new TF1("fs_fo_data", "[0]*TMath::BreitWigner(x,[1],[2])", m2pi_min, m2pi_max);
  fs_fo_data->SetLineColor(4);

  TF1 *fb_fo_data = new TF1("fb_fo_data", "pol3", m2pi_min, m2pi_max); //pol3(3)
  fb_fo_data->SetLineColor(28);
  fb_fo_data->SetLineStyle(2);

  gphifo->Fit("fsb_fo_data", "", "", m2pi_min, m2pi_max);
  double par_fo_data[10]; //6
  fsb_fo_data->GetParameters(&par_fo_data[0]);
  fs_fo_data->SetParameters(&par_fo_data[0]);
  fb_fo_data->SetParameters(&par_fo_data[3]); //3

  double N_fo_data = fs_fo_data->Integral(m2pi_min, m2pi_max) * n2pi2k / (m2pi_max - m2pi_min);
  double dN_fo_data = N_fo_data * fsb_fo_data->GetParError(0) / fsb_fo_data->GetParameter(0);

  fs_fo_data->Draw("same");
  fb_fo_data->Draw("same");
  // gPad->Modified();
  // gPad->Update(); 

  TLatex lat_phifo;
  lat_phifo.SetTextSize(0.05);
  lat_phifo.SetTextAlign(13); //align at top
  lat_phifo.SetNDC();
  lat_phifo.SetTextColor(kBlue);
  lat_phifo.DrawLatex(0.6, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_fo_data->GetChisquare() / fsb_fo_data->GetNDF()));
  lat_phifo.DrawLatex(0.6, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_fo_data, dN_fo_data));
  lat_phifo.DrawLatex(0.6, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(1), fsb_fo_data->GetParError(1)));
  lat_phifo.DrawLatex(0.6, 0.68, Form("#Gamma = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(2), fsb_fo_data->GetParError(2)));

  double bandin = 0.08, bandout = 0.14;
  //data_17: (1sigma: 0.08, 0.11), (3sigma: 0.05, 0.08), (6sigma: 0.1, 0.15)
  //data_18: (1sigma: 0.08, 0.12), (3sigma: 0.09, 0.13), (3sigma: 0.1, 0.15)

  TLine *l_phifo_1 = new TLine(fsb_fo_data->GetParameter(1)-bandin, 0, fsb_fo_data->GetParameter(1)-bandin, gphifo->Eval(fsb_fo_data->GetParameter(1)-bandin));
  TLine *l_phifo_2 = new TLine(fsb_fo_data->GetParameter(1)+bandin, 0, fsb_fo_data->GetParameter(1)+bandin, gphifo->Eval(fsb_fo_data->GetParameter(1)+bandin));
  TLine *l_phifo_3 = new TLine(fsb_fo_data->GetParameter(1)-bandout, 0, fsb_fo_data->GetParameter(1)-bandout, gphifo->Eval(fsb_fo_data->GetParameter(1)-bandout));
  TLine *l_phifo_4 = new TLine(fsb_fo_data->GetParameter(1)+bandout, 0, fsb_fo_data->GetParameter(1)+bandout, gphifo->Eval(fsb_fo_data->GetParameter(1)+bandout));

  l_phifo_1->SetLineColor(kRed);
  l_phifo_1->SetLineWidth(2);
  // l_phifo_1->SetLineStyle(10);
  l_phifo_1->Draw("same");
  l_phifo_2->SetLineColor(kRed);
  l_phifo_2->SetLineWidth(2);
  // l_phifo_2->SetLineStyle(10);
  l_phifo_2->Draw("same");
  l_phifo_3->SetLineColor(kMagenta);
  l_phifo_3->SetLineWidth(2);
  // l_phifo_3->SetLineStyle(10);
  l_phifo_3->Draw("same");
  l_phifo_4->SetLineColor(kMagenta);
  l_phifo_4->SetLineWidth(2);
  // l_phifo_4->SetLineStyle(10);
  l_phifo_4->Draw("same");

    // ++++++ convert to histogram
    TCanvas *chgphifo = new TCanvas("chgphifo", "chgphifo", 900, 600);
    double widphifo = gphifo->GetX()[1] - gphifo->GetX()[0];
    TH1F *hgphifo = new TH1F("hgphifo", ";m_{#phif_{0}} [GeV/c^{2}];N_{#phi}", gphifo->GetN(), gphifo->GetX()[0] - 0.5 * widphifo, gphifo->GetX()[gphifo->GetN() - 1] + 0.5 * widphifo);
    for (int i = 0; i < gphifo->GetN(); ++i)
    {
      hgphifo->SetBinContent(i + 1, gphifo->GetY()[i]);
      hgphifo->SetBinError(i + 1, gphifo->GetEY()[i]);
    }
    chgphifo->cd();
    hgphifo->Draw("e");
    outputfig->cd();
    hgphifo->Write(Form("h%s_phifo", name.Data()), TObject::kWriteDelete);
    hgphifo->Sumw2();

  // ======================================== Y vs. Phi(1020) with fo signal===================================
  // root -l 'ul_yphifo.C+(100,50,1,1)'

  TCanvas *cphiy=new TCanvas("cphiy","cphiy",900, 600);//1500, 800

  TCanvas *cphiy1 = new TCanvas("cphiy1","cphiy1", 1500, 800);//1500, 800
  cphiy1->Divide(5, 5);
  TCanvas *cphiy2 = new TCanvas("cphiy2", "cphiy2", 1500, 800);
  cphiy2->Divide(5, 5);

  TCanvas *cgphiy=new TCanvas("cgphiy","cgphiy",900, 600);//1500, 800
  TGraphErrors *gphiy;

  // TCanvas *cgnophiy=new TCanvas("cgnophiy","cgnophiy",1500, 800);
  // TGraphErrors *gnophiy= new TGraphErrors;//(n2pi2k)
  // gnophiy->SetMarkerStyle(20);
  // gnophiy->SetMarkerSize(1.0);
  // gnophiy->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  cphiy->cd();
  TH2D *h1d2 = new TH2D("h1d2",";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, n2k, mkk_min, mkk_max);
  tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni) && abs(pippim_mf-%f)<%f)",fsb_fo_data->GetParameter(1), bandin));//<0.1, 0.08, 0.05

  h1d2->Draw("colz");

  gphiy = new TGraphErrors(); //n2pi2k
  gphiy->SetMarkerStyle(20);
  gphiy->SetMarkerSize(1.0);
  gphiy->SetMarkerColor(1);
  gphiy->SetTitle(";m_{#phif_{0}} [GeV/c^{2}];N_{#phi}");

  // gnophiy->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1, Eg2));

  for (int i = 1; i <= n2pi2k; ++i) //n2pi2k
  {
    // cout << i << " " << flush;

    if (i < 26)
      cphiy1->cd(i); //i
    if (i > 25 && i < 51)
      cphiy2->cd(i - 25);
    // if(i > 50 && i<76) cphiy3->cd(i-50);
    // if(i > 75) cphiy4->cd(i-75);

    TH1D *hphiy_py = h1d2->ProjectionY(Form("_hphiy_py_%d", i), i, i);

    // hphiy_py->Draw();
    // if (hphiy_py->Integral(1, 50) < 100)
    //   continue;

    // TF1 *fsb_phiy_data = new TF1("fsb_phiy_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
    TF1 *fsb_phiy_data = new TF1("fsb_phiy_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 0.98, 1.2);
    fsb_phiy_data->SetLineColor(2);
    fsb_phiy_data->SetParameters(1433, 1.019, 0.004, 0.002, 1,1,1,1,1);
    fsb_phiy_data->SetParLimits(0, 0, 10000);
    fsb_phiy_data->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb_phiy_data->FixParameter(2, w.var("sigma_phi_mc")->getVal());   //fsb_phiy_data->SetParLimits(2, 0.008, 0.010); // 0.0042
    fsb_phiy_data->FixParameter(3, w.var("width_phi_mc")->getVal());        //fsb_phiy_data->SetParLimits(3, 0.001,0.01);// 0.002

    // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    TF1 *fs_phiy_data = new TF1("fs_phiy_data", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
    fs_phiy_data->SetLineColor(4);

    TF1 *fb_phiy_data = new TF1("fb_phiy_data", "pol4(4)", mkk_min, mkk_max); //3
    fb_phiy_data->SetLineColor(28);
    fb_phiy_data->SetLineStyle(2);

    hphiy_py->Fit("fsb_phiy_data", "", "", mkk_min, mkk_max);
    double par[10]; //6
    fsb_phiy_data->GetParameters(&par[0]);
    fs_phiy_data->SetParameters(&par[0]);
    fb_phiy_data->SetParameters(&par[4]); //3

    fs_phiy_data->Draw("same");
    fb_phiy_data->Draw("same");

    double N1 = fs_phiy_data->Integral(mkk_min, mkk_max) / hphiy_py->GetBinWidth(1);
    double dN1 = N1 * fsb_phiy_data->GetParError(0) / fsb_phiy_data->GetParameter(0);

    // if (N1 <= 0)
    //   gphiy->RemovePoint(i - 1);

    gphiy->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
    gphiy->SetPointError(i - 1, 0, dN1);

    // gnophiy->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
    // gnophiy->SetPointError(i - 1, 0, dNbkg);
    // ofs_ul_yphifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

    TLatex lat_phiy;
    lat_phiy.SetTextSize(0.05);
    lat_phiy.SetTextAlign(13); //align at top
    lat_phiy.SetNDC();
    lat_phiy.SetTextColor(kBlue);
    lat_phiy.DrawLatex(0.45, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_phiy_data->GetChisquare() / fsb_phiy_data->GetNDF()));
    lat_phiy.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
    lat_phiy.DrawLatex(0.45, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_phiy_data->GetParameter(1), fsb_phiy_data->GetParError(1)));
    lat_phiy.DrawLatex(0.45, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_phiy_data->GetParameter(2), fsb_phiy_data->GetParError(2)));
    lat_phiy.DrawLatex(0.45, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_phiy_data->GetParameter(3), fsb_phiy_data->GetParError(3)));

    // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

    // hphiy_py->Write();
    // cgphiy->Update();
    // c2->Update();
    //sleep(1);
    }

    // cout << endl;

    cgphiy->cd();
    gphiy->Draw("AP");
    // gphiy->SetMinimum(-70.);

    // ++++++ convert to histogram
    TCanvas *chgphiy = new TCanvas("chgphiy", "chgphiy", 900, 600);
    double widphiy = gphiy->GetX()[1] - gphiy->GetX()[0];
    TH1F *hgphiy = new TH1F("hgphiy", ";m_{#phif_{0}} [GeV/c^{2}];N_{#phi}", gphiy->GetN(), gphiy->GetX()[0] - 0.5 * widphiy, gphiy->GetX()[gphiy->GetN() - 1] + 0.5 * widphiy);
    for (int i = 0; i < gphiy->GetN(); ++i)
    {
      hgphiy->SetBinContent(i + 1, gphiy->GetY()[i]);
      hgphiy->SetBinError(i + 1, gphiy->GetEY()[i]);
    }
    chgphiy->cd();
    hgphiy->Draw("e");
    outputfig->cd();
    hgphiy->Write(Form("h%s_phiy", name.Data()), TObject::kWriteDelete);
    hgphiy->Sumw2();

    // ======================================== Y vs. Phi(1020) with fo side band ================================
    // root -l 'ul_yphifo.C+(100,50,1,1)'

    TCanvas *cphiy_nofo = new TCanvas("cphiy_nofo", "cphiy_nofo", 900, 600); //1500, 800

    TCanvas *cphiy_nofo1 = new TCanvas("cphiy_nofo1", "cphiy_nofo1", 1500, 800); //1500, 800
    cphiy_nofo1->Divide(5, 5);
    TCanvas *cphiy_nofo2 = new TCanvas("cphiy_nofo2", "cphiy_nofo2", 1500, 800);
    cphiy_nofo2->Divide(5, 5);

    TCanvas *cgphiy_nofo = new TCanvas("cgphiy_nofo", "cgphiy_nofo", 900, 600); //1500, 800
    TGraphErrors *gphiy_nofo;

    // TCanvas *cgnophiy_nofo=new TCanvas("cgnophiy_nofo","cgnophiy_nofo",1500, 800);
    // TGraphErrors *gnophiy_nofo= new TGraphErrors;//(n2pi2k)
    // gnophiy_nofo->SetMarkerStyle(20);
    // gnophiy_nofo->SetMarkerSize(1.0);
    // gnophiy_nofo->SetMarkerColor(2);

    // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

    cphiy_nofo->cd();
    TH2D *h1d2_nofo = new TH2D("h1d2_nofo", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, n2k, mkk_min, mkk_max);
    tdata->Project("h1d2_nofo", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni) && abs(pippim_mf-%f)>%f & abs(pippim_mf-%f)<%f)",fsb_fo_data->GetParameter(1), bandin, fsb_fo_data->GetParameter(1), bandout));  // abs(pippim_mf-%f)>0.1 & abs(pippim_mf-%f)<0.15, abs(pippim_mf-%f)>0.08 & abs(pippim_mf-%f)<0.15, 

    h1d2_nofo->Draw("colz");

    gphiy_nofo = new TGraphErrors(); //n2pi2k
    gphiy_nofo->SetMarkerStyle(20);
    gphiy_nofo->SetMarkerSize(1.0);
    gphiy_nofo->SetMarkerColor(1);
    gphiy_nofo->SetTitle(";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    // gnophiy_nofo->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1, Eg2));

    for (int i = 1; i <= n2pi2k; ++i) //n2pi2k
    {
      // cout << i << " " << flush;

      if (i < 26)
        cphiy_nofo1->cd(i); //i
      if (i > 25 && i < 51)
        cphiy_nofo2->cd(i - 25);
      // if(i > 50 && i<76) cphiy_nofo3->cd(i-50);
      // if(i > 75) cphiy_nofo4->cd(i-75);

      TH1D *hphiy_nofo_py = h1d2_nofo->ProjectionY(Form("_hphiy_nofo_py_%d", i), i, i);

      // hphiy_nofo_py->Draw();
      // if (hphiy_nofo_py->Integral(1, 50) < 100)
      //   continue;

      // TF1 *fsb_phiy_nofo_data = new TF1("fsb_phiy_nofo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb_phiy_nofo_data = new TF1("fsb_phiy_nofo_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 0.98, 1.2);
      fsb_phiy_nofo_data->SetLineColor(2);
      fsb_phiy_nofo_data->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1, 1, 1);
      fsb_phiy_nofo_data->SetParLimits(0, 0, 10000);
      fsb_phiy_nofo_data->SetParLimits(1, 1.015, 1.022);                        // 1.018, 1.021
      fsb_phiy_nofo_data->FixParameter(2, w.var("sigma_phi_mc")->getVal()); //fsb_phiy_nofo_data->SetParLimits(2, 0.008, 0.010); // 0.0042
      fsb_phiy_nofo_data->FixParameter(3, w.var("width_phi_mc")->getVal()); //fsb_phiy_nofo_data->SetParLimits(3, 0.001,0.01);// 0.002

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs_phiy_nofo_data = new TF1("fs_phiy_nofo_data", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
      fs_phiy_nofo_data->SetLineColor(4);

      TF1 *fb_phiy_nofo_data = new TF1("fb_phiy_nofo_data", "pol4(4)", mkk_min, mkk_max); //3
      fb_phiy_nofo_data->SetLineColor(28);
      fb_phiy_nofo_data->SetLineStyle(2);

      hphiy_nofo_py->Fit("fsb_phiy_nofo_data", "", "", mkk_min, mkk_max);
      double par[10]; //6
      fsb_phiy_nofo_data->GetParameters(&par[0]);
      fs_phiy_nofo_data->SetParameters(&par[0]);
      fb_phiy_nofo_data->SetParameters(&par[4]); //3

      fs_phiy_nofo_data->Draw("same");
      fb_phiy_nofo_data->Draw("same");

      double N1 = fs_phiy_nofo_data->Integral(mkk_min, mkk_max) / hphiy_nofo_py->GetBinWidth(1);
      double dN1 = N1 * fsb_phiy_nofo_data->GetParError(0) / fsb_phiy_nofo_data->GetParameter(0);

      // if (N1 <= 0)
      //   gphiy_nofo->RemovePoint(i - 1);

      gphiy_nofo->SetPoint(i - 1, h1d2_nofo->GetXaxis()->GetBinCenter(i), N1);
      gphiy_nofo->SetPointError(i - 1, 0, dN1);

      // gnophiy_nofo->SetPoint(i - 1, h1d2_nofo->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophiy_nofo->SetPointError(i - 1, 0, dNbkg);
      // ofs_ul_yphifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2_nofo->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2_nofo->GetYaxis()->GetBinCenter(i)<<endl;

      TLatex lat_phiy_nofo;
      lat_phiy_nofo.SetTextSize(0.05);
      lat_phiy_nofo.SetTextAlign(13); //align at top
      lat_phiy_nofo.SetNDC();
      lat_phiy_nofo.SetTextColor(kBlue);
      lat_phiy_nofo.DrawLatex(0.45, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_phiy_nofo_data->GetChisquare() / fsb_phiy_nofo_data->GetNDF()));
      lat_phiy_nofo.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiy_nofo.DrawLatex(0.45, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_phiy_nofo_data->GetParameter(1), fsb_phiy_nofo_data->GetParError(1)));
      lat_phiy_nofo.DrawLatex(0.45, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_phiy_nofo_data->GetParameter(2), fsb_phiy_nofo_data->GetParError(2)));
      lat_phiy_nofo.DrawLatex(0.45, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_phiy_nofo_data->GetParameter(3), fsb_phiy_nofo_data->GetParError(3)));

      // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

      // hphiy_nofo_py->Write();
      // cgphiy_nofo->Upda
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphiy_nofo->cd();
    gphiy_nofo->Draw("AP");
    // gphiy_nofo->SetMinimum(-70.);

  // ++++++ convert to histogram
    TCanvas *chgphiy_nofo = new TCanvas("chgphiy_nofo", "chgphiy_nofo", 900, 600);
    double wid_nofo = gphiy_nofo->GetX()[1] - gphiy_nofo->GetX()[0];
    TH1F *hgphiy_nofo = new TH1F("hgphiy_nofo", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", gphiy_nofo->GetN(), gphiy_nofo->GetX()[0] - 0.5 * wid_nofo, gphiy_nofo->GetX()[gphiy_nofo->GetN() - 1] + 0.5 * wid_nofo);
    for (int i = 0; i < gphiy_nofo->GetN(); ++i)
    {
      hgphiy_nofo->SetBinContent(i + 1, gphiy_nofo->GetY()[i]);
      hgphiy_nofo->SetBinError(i + 1, gphiy_nofo->GetEY()[i]);
    }
    chgphiy->cd(); //chgphiy_nofo
    hgphiy_nofo->Draw("same");
    outputfig->cd();
    hgphiy_nofo->Write(Form("h%s_phiy_nofo", name.Data()), TObject::kWriteDelete);
    hgphiy_nofo->Sumw2();

    // =============== Subtracted
    TCanvas *chgphiy_sub = new TCanvas("chgphiy_sub", "chgphiy_sub", 900, 600);
    chgphiy_sub->Divide(2);
    TH1F *hgphiy_sub =  new TH1F("hgphiy_sub", ";m_{#phif_{0}} [GeV/c^{2}];N_{#phi}", gphiy->GetN(), gphiy->GetX()[0] - 0.5 * widphiy, gphiy->GetX()[gphiy->GetN() - 1] + 0.5 * widphiy);
    hgphiy_sub->Add(hgphiy, hgphiy_nofo, 1., -1.);
    cout << "hgphiy_sub = " << hgphiy_sub << endl;
    chgphiy_sub->cd();
    hgphiy_sub->Draw("e");
    outputfig->cd();
    hgphiy_sub->Write(Form("h%s_phiy_sub", name.Data()), TObject::kWriteDelete);

    // =========== Y(2175) with fo signal region fitted with Likelihood ===========

    w.factory(Form("Voigtian::sig_y_data(m_y_data[%f, %f],mean_y_data[%f],width_y_data[%f],sigma_y_data[%f])", m2pi2k_min2, m2pi2k_max2, fsb_y_mc->GetParameter(1),fsb_y_mc->GetParameter(3), fsb_y_mc->GetParameter(2))); // m_y_data[1.93,3], m2pi2k_min2, m2pi2k_max2,  //sigma_y_data[0.0001,0.01], mean_y_data[1.011,1.030]
    // w.factory(Form("BreitWigner::sig_y_data(m_y_data[%f,%f],mean_y_data[1.018,1.021],width_y_data[0.004])", mkk_min, mkk_max));
    w.factory("Chebychev::bkg_y_data(m_y_data,{c0_y_data[-INF,+INF], c1_y_data[-INF,+INF], c2_y_data[-INF,+INF]})");// c0_y_data[-100000,+100000], c1_y_data[-100000,+100000], c2_y_data[-100000,+100000]
    w.factory("Chebychev::bkg_y_data2(m_y_data,{c0_y_data[-INF,+INF], c1_y_data[-INF,+INF], c2_y_data[-INF,+INF]})");// c0_y_data[-100000,+100000], c1_y_data[-100000,+100000], c2_y_data[-100000,+100000]
    w.factory("SUM:::model_y_data(nsig_y_data[0,-1000,+1000]*sig_y_data, nbkg_y_data[0,+100000]*bkg_y_data)"); //nsig_y_data[-100000,+100000]*sig_y_data, nbkg_y_data[0,+100000]*bkg_y_data
    w.factory("SUM:::model_y_data2(nbkg_y_data[0,+100000]*bkg_y_data2)"); //nbkg_y_data[0,+100000]*bkg_y_data2
    w.var("m_y_data")->SetTitle("m_{#phif_{0}} [GeV/c^{2}]");
    RooDataHist dh_y_data("dh_y_data", "dh_y_data", *w.var("m_y_data"), Import(*hgphiy));//hgphiy_sub
    RooPlot *fr_y_data = w.var("m_y_data")->frame(Title(" "));

    RooAbsPdf *model1 = w.pdf("model_y_data");
    RooAbsPdf *model2 = w.pdf("model_y_data2");
    cout <<model1<<"  "<<model2<<" | minll_sig_plus_bg = "<<endl; 
    
    RooFitResult *fitFull = model1->fitTo(dh_y_data, Save(true));
    double minll_sig_plus_bg = fitFull->minNll();

    RooFitResult *fitbg = model2->fitTo(dh_y_data, Save(true));
    double minll_bg = fitbg->minNll();

    double Z;
    // if(name == "17") Z = sqrt(-2 * (minll_sig_plus_bg - minll_bg)); 
    // if(name == "18") Z = sqrt(2 * (minll_sig_plus_bg - minll_bg));
    Z = sqrt(-2 * (minll_sig_plus_bg - minll_bg));

    double N_y_data = w.var("nsig_y_data")->getVal();
    double dN_y_data = w.var("nsig_y_data")->getError();

    ofs_yphifo << "model1 = " << model1 <<" | minll_sig_plus_bg = "<<minll_sig_plus_bg<< " | model2 = " << model2<<" | minll_bg = "<<minll_bg<<" | Z = "<<Z<<endl; // " | minll_bg = "<<minll_bg<<" | Z = "<<Z<<

    // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    dh_y_data.plotOn(fr_y_data, RooFit::Name("ldh_y_data"));
    w.pdf("model_y_data")->plotOn(fr_y_data, Components(*w.pdf("sig_y_data")), LineColor(kRed), RooFit::Name("lsig_y_data"));
    w.pdf("model_y_data")->plotOn(fr_y_data, Components(*w.pdf("bkg_y_data")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_y_data"));
    w.pdf("model_y_data")->plotOn(fr_y_data, RooFit::Name("lmodel_y_data"));
    RooHist* hpull = fr_y_data->pullHist();
    w.pdf("model_y_data2")->plotOn(fr_y_data, LineColor(kMagenta));
   // w.pdf("model_y_data")->paramOn(fr_y_data, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_y_data"), *w.var("nbkg_y_data")))); //,*w.var("mean_y_data"),*w.var("width_y_data"),*w.var("sigma_y_data"))));
    fr_y_data->Draw();

    // TLegend *l_y_data = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_y_data->SetFillColor(kWhite);
    // l_y_data->AddEntry(fr_y_data->findObject("lsig_y_data"), Form("N_{Sig} = %.2f", w.var("nsig_y_data")->getVal()), "l");
    // l_y_data->AddEntry(fr_y_data->findObject("lbkg_y_data"), Form("N_{Bkg} = %.2f", w.var("nbkg_y_data")->getVal()), "l");
    // l_y_data->Draw();

    // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
    // -------------------------------------------------------

    // // Construct a histogram with the residuals of the data w.r.t. the curve
    // RooHist* hresid = fr_y_data->residHist();

    // // Construct a histogram with the pulls of the data w.r.t the curve
    // RooHist* hpull = fr_y_data->pullHist();

    // // Create a new frame to draw the residual distribution and add the distribution to the frame
    // RooPlot* fr_resi = w.var("m_y_data")->frame(Title("Residual Distribution"));
    // fr_resi->addPlotable(hresid, "P");

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot *fr_pull = w.var("m_y_data")->frame(Title(" "));
    fr_pull->addPlotable(hpull, "P");

    fr_pull->SetStats(0);
    fr_pull->GetXaxis()->SetLabelSize(0.18);
    fr_pull->GetXaxis()->SetTitleSize(0.18);
    fr_pull->GetXaxis()->SetTitleOffset(0.95);
    fr_pull->GetXaxis()->SetTickLength(0.05);

    fr_pull->GetYaxis()->SetLabelSize(0.15);
    fr_pull->GetYaxis()->SetTitleSize(0.15);
    fr_pull->GetYaxis()->SetTitleOffset(0.32);
    fr_pull->GetYaxis()->SetNdivisions(5);

    // fr_pull->SetXTitle(fr_y_data->GetXaxis()->GetTitle());
    fr_pull->SetYTitle("pull [#sigma] ");

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);//0, 0.25, 1, 1
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.251);//0, 0, 1, 0.251

    pad1->SetBottomMargin(0.0);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.4);//0.4
    pad2->SetBorderMode(0);

    pad1->Draw();
    pad2->Draw();

    // TCanvas* c_phiy_sub_resipull = new TCanvas("c_phiy_sub_resipull", "c_phiy_sub_resipull", 900, 300);
    // // c_phiy_sub_resipull->Divide(3);
    // c_phiy_sub_resipull->cd();
    // gPad->SetLeftMargin(0.15);
    // fr_y_data->GetYaxis()->SetTitleOffset(1.6);
    pad1->cd();
    fr_y_data->Draw();

    TLatex lat_phiye;
    lat_phiye.SetTextSize(0.05);
    lat_phiye.SetTextAlign(13); //align at top
    lat_phiye.SetNDC();
    lat_phiye.SetTextColor(kBlue);
    // lat_phiye.DrawLatex(0.68, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phiye.DrawLatex(0.68, 0.86, Form("N_{sig} = %0.2f#pm%0.2f", N_y_data, dN_y_data));
    lat_phiye.DrawLatex(0.68, 0.80, Form("#mu = %0.3f#pm%0.3f", w.var("mean_y_data")->getVal(), w.var("mean_y_data")->getError()));
    lat_phiye.DrawLatex(0.68, 0.74, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_y_data")->getVal(), w.var("sigma_y_data")->getError()));
    lat_phiye.DrawLatex(0.68, 0.68, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_y_data")->getVal(), w.var("width_y_data")->getError()));
    lat_phiye.DrawLatex(0.68, 0.62, Form("Z = %0.3f", Z));

    // c_phiy_sub_resipull->cd(2);
    // gPad->SetLeftMargin(0.15);
    // fr_resi->GetYaxis()->SetTitleOffset(1.6);
    // fr_resi->Draw();
    // c_phiy_sub_resipull->cd(3);
    // gPad->SetLeftMargin(0.15);
    // fr_pull->GetYaxis()->SetTitleOffset(1.6);
    pad2->cd();
    pad2->SetGridy();
    fr_pull->Draw();

    // Enhance line at 0
    TLine l;
    l.SetLineColor(4);
    l.DrawLine(fr_pull->GetXaxis()->GetXmin(), 0, fr_pull->GetXaxis()->GetXmax(), 0);

    // =========== Y(2175) with fo side bands subtracted fitted with Chi2 ===========

    TCanvas *chgphiy_sub2 = new TCanvas("chgphiy_sub2", "chgphiy_sub2", 900, 600);
    TH1F *hgphiy_sub2 =  new TH1F("hgphiy_sub2", ";m_{#phif_{0}} [GeV/c^{2}];N_{#phi}", gphiy->GetN(), gphiy->GetX()[0] - 0.5 * widphiy, gphiy->GetX()[gphiy->GetN() - 1] + 0.5 * widphiy);
    hgphiy_sub2->Add(hgphiy, hgphiy_nofo, 1., -1.);
    cout << "hgphiy_sub2 = " << hgphiy_sub2 << endl;
    chgphiy_sub2->cd();
    // hgphiy_sub2->Draw("e");
    // hgphiy_sub2->Write();

    TF1 *fsb_phiy_sub2 = new TF1("fsb_phiy_sub2", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol3(4)", m2pi2k_min2, m2pi2k_max2);
    // TF1 *fsb_phiy_sub2 = new TF1("fsb_phiy_sub2", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", m2pi2k_min2, m2pi2k_max2);
    // TF1 *fsb_phiy_sub2 = new TF1("fsb_phiy_sub2", "gaus(0) + pol3(3)", m2pi2k_min2, m2pi2k_max2);
    fsb_phiy_sub2->SetLineColor(2);
    fsb_phiy_sub2->SetParameters(200, 2.08, 0.02, 0.02, 1, 1, 1, 1);
    fsb_phiy_sub2->SetNpx(500);
    // fsb_phiy_sub2->SetParLimits(0, 0, 10000);
    fsb_phiy_sub2->SetParLimits(1, 2.0, 2.18);//2.2                    // 1.018, 1.021
    fsb_phiy_sub2->SetParLimits(2, 0.001, 0.1); //fsb_phiy_sub2->SetParLimits(2, 0.008, 0.010); // 0.0042
    fsb_phiy_sub2->SetParLimits(3, 0.001, 0.1); //fsb_phiy_sub2->SetParLimits(2, 0.008, 0.010); // 0.0042

    TF1 *fs_phiy_sub2 = new TF1("fs_phiy_sub2", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    // TF1 *fs_phiy_sub2 = new TF1("fs_phiy_sub2", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    // TF1 *fs_phiy_sub2 = new TF1("fs_phiy_sub2", "gaus(0)", m2pi2k_min2, m2pi2k_max2);//gaus(0)
    fs_phiy_sub2->SetLineColor(4);

    TF1 *fb_phiy_sub2 = new TF1("fb_phiy_sub2", "pol3(4)", m2pi2k_min2, m2pi2k_max2); //3
    fb_phiy_sub2->SetLineColor(28);
    fb_phiy_sub2->SetLineStyle(2);

    hgphiy_sub2->Fit("fsb_phiy_sub2", "", "", m2pi2k_min2, m2pi2k_max2);
    double par[10]; //6
    fsb_phiy_sub2->GetParameters(&par[0]);
    fs_phiy_sub2->SetParameters(&par[0]);
    fb_phiy_sub2->SetParameters(&par[4]); //3

    // TCanvas *chgphiy_sub2_pull = new TCanvas("chgphiy_sub2_pull", "chgphiy_sub2_pull", 900, 600);
    TH1F *hgphiy_sub2_pull = showpull(hgphiy_sub2, fsb_phiy_sub2);
    // chgphiy_sub2_pull->cd();
    hgphiy_sub2_pull->Draw("e");
    hgphiy_sub2->Draw("e");
    fs_phiy_sub2->Draw("same");
    fb_phiy_sub2->Draw("same");


    double N1 = fs_phiy_sub2->Integral(m2pi2k_min2, m2pi2k_max2) / hgphiy_sub2->GetBinWidth(1);
    double dN1 = N1 * fsb_phiy_sub2->GetParError(0) / fsb_phiy_sub2->GetParameter(0);
    double Z2 = N1 / dN1;
    //sqrt(N1 + fb_phiy_sub2->Integral(2.05, 2.15) / hgphiy_sub2->GetBinWidth(1));

    // if (N1 <= 0)
    //   gphiy_nofo->RemovePoint(i - 1);

    // gnophiy_nofo->SetPoint(i - 1, h1d2_nofo->GetXaxis()->GetBinCenter(i), Nbkg);
    // gnophiy_nofo->SetPointError(i - 1, 0, dNbkg);
    // ofs_ul_yphifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2_nofo->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2_nofo->GetYaxis()->GetBinCenter(i)<<endl;

    TLatex lat_phiy_sub2;
    lat_phiy_sub2.SetTextSize(0.05);
    lat_phiy_sub2.SetTextAlign(13); //align at top
    lat_phiy_sub2.SetNDC();
    lat_phiy_sub2.SetTextColor(kBlue);
    lat_phiy_sub2.DrawLatex(0.62, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_phiy_sub2->GetChisquare() / fsb_phiy_sub2->GetNDF()));
    lat_phiy_sub2.DrawLatex(0.62, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
    lat_phiy_sub2.DrawLatex(0.62, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_phiy_sub2->GetParameter(1), fsb_phiy_sub2->GetParError(1)));
    lat_phiy_sub2.DrawLatex(0.62, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_phiy_sub2->GetParameter(2), fsb_phiy_sub2->GetParError(2)));
    lat_phiy_sub2.DrawLatex(0.62, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_phiy_sub2->GetParameter(3), fsb_phiy_sub2->GetParError(3)));
    lat_phiy_sub2.DrawLatex(0.62, 0.56, Form("Z = %0.3f", Z2));

    // // +++++++++++++++++++ Systematic Errors
    // double bkg_poln[3] = {-274, -286, -261}; // Bkg pol degree: 4, 3, 5
    // double err_bkg_poln = TMath::RMS(3, bkg_poln);
    // double fit_range[3] = {-274, -320, -225}; // range phi2pi:[2, 3], [2, 2.8], [2, 3.2]
    // double err_fit_range = TMath::RMS(3, fit_range);
    // double slice_numb[3] = {-274, -339, -325}; // K+K- bining: 50, 40, 60
    // double err_slice_numb = TMath::RMS(3, slice_numb);
    // double mean_value[3] = {-274, -278, -269}; // Y Mean value : 2.194, 2.183, 2.205
    // double err_mean_value = TMath::RMS(3, mean_value);
    // double width_value[3] = {-274, -230, -319}; // Y Width value : 0.079, 0.065, 0.093
    // double err_width_value = TMath::RMS(3, width_value);

    // double err_sys = TMath::Sqrt(err_bkg_poln * err_bkg_poln + err_fit_range * err_fit_range + err_slice_numb * err_slice_numb + err_mean_value * err_mean_value + err_width_value * err_width_value);

  // ++++++++++++++++++++++++++++ Cross-section  +++++++++++++++++++++++
    double N_y_gen = h_beame_tru->Integral(1,100);
    double eff_y = N_y_mc / N_y_gen; // Efficiency = N_observed/N_generated
    double deff_y = eff_y * (dN_y_mc / N_y_mc);
    double lumi_y = h_tagged_flux->Integral(1,100); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_y = 1e9 * N_y_data / (eff_y * lumi_y * T * br_phi);
    double dxsec_y_stat = xsec_y * (dN_y_data / N_y_data);
    double dxsec_y_sys = xsec_y * TMath::Sqrt((deff_y / eff_y) * (deff_y / eff_y) + (dbr_phi / br_phi) * (dbr_phi / br_phi)); // (err_sys / N_y_data) * (err_sys / N_y_data) + 
    double dxsec_y_tot = TMath::Sqrt(dxsec_y_stat*dxsec_y_stat + dxsec_y_sys*dxsec_y_sys);

    // // ++++++++++++++++++++++++++++ Upper Limit  90% CL  +++++++++++++++++

    // // ++++++++++ gaussian method
    // // double ul = mu + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, sig, mu), sig);
    // // double ul = xsec_y + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, dxsec_y_tot, xsec_y), dxsec_y_tot);

    // // ******** hide messages from roofit
    // gROOT->ProcessLine( "gErrorIgnoreLevel = 4001;");   
    // gROOT->ProcessLine( "gPrintViaErrorHandler = kTRUE;");   
    
    // if (RooMsgService::instance().numStreams()>0)
    // {
    //     RooMsgService::instance().deleteStream(0);
    //     RooMsgService::instance().deleteStream(1);
    // } 

    // // ++++++++++ profile likelihood method
    // double dN_y_tot = sqrt(dN_y_data*dN_y_data); // + err_sys*err_sys
    // TGraph *gprof = profileLikelihoodScan(dh_y_data, *(w.pdf("model_y_data")), w.var("nsig_y_data"), N_y_data - 4 * dN_y_tot, N_y_data + 4 * dN_y_tot);//N_y_data - 4 * dN_y_tot, N_y_data + 4 * dN_y_tot

    // TCanvas *cprof = new TCanvas("cprof","",900,600);
    // cprof->cd();
    // cprof->SetGrid();
    // gprof->SetTitle("; Yield_{#phi#pi^{+}#pi^{-}}; -Log(L)");
    // gprof->SetLineWidth(2);
    // gprof->SetMarkerStyle(20);
    // gprof->SetMarkerSize(1.0);
    // gprof->Draw("AP");
    // // gprof->Fit("pol2", "R", "", gprof->GetHistogram()->GetXaxis()->GetXmin(), gprof->GetHistogram()->GetXaxis()->GetXmax());

    // TCanvas *cprofxsec = new TCanvas("cprofxsec","",900,600);
    // TGraph *gprofxsec = new TGraph();
    // for (int i = 0; i < gprof->GetN(); i++)
    // {
    //   double x_gprof, y_gprof;
    //   gprof->GetPoint(i,x_gprof,y_gprof);
    //   double xsec_gprof = 1e9 * x_gprof / (eff_y * lumi_y * T * br_phi);
    //   cout << "i = " << i << " | x_gprof = "<<x_gprof<<" | y_gprof = " <<y_gprof<<" | xsec_gprof = " <<xsec_gprof<<endl;
    //   gprofxsec->SetPoint(i, xsec_gprof, exp(-y_gprof));
    // }
    // cprofxsec->cd();
    // gprofxsec->SetTitle("; #sigma [nb]; Likelihood");
    // gprofxsec->SetLineWidth(2);
    // gprofxsec->SetMarkerStyle(20);
    // gprofxsec->SetMarkerSize(1.0);
    // gprofxsec->Draw("AP");
    // TF1 *func_profxsec = new TF1("func_profxsec", "gaus", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());
    // gprofxsec->Fit(func_profxsec, "R", "", gprofxsec->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec->GetHistogram()->GetXaxis()->GetXmax());
    // double ul = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, func_profxsec->GetParameter(2), func_profxsec->GetParameter(1)), func_profxsec->GetParameter(2));
    // double ul2 = func_profxsec->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9,func_profxsec->GetParameter(2));
    // TLine *l_mean = new TLine(func_profxsec->GetParameter(1), 0, func_profxsec->GetParameter(1), 1);
    // l_mean->SetLineWidth(2);
    // l_mean->SetLineStyle(10);
    // l_mean->Draw("same");
    // TLine *l_sigma = new TLine(func_profxsec->GetParameter(1), 0.5, func_profxsec->GetParameter(1)+func_profxsec->GetParameter(2), 0.5);
    // l_sigma->SetLineWidth(2);
    // l_sigma->SetLineStyle(10);
    // l_sigma->Draw("same");
    // TLine *l_ul = new TLine(ul, 0, ul, 0.45);
    // l_ul->SetLineWidth(2);
    // l_ul->SetLineColor(kBlue);
    // l_ul->Draw("same");
    // TLine *l_ul2 = new TLine(ul2, 0, ul2, 0.45);
    // l_ul2->SetLineWidth(2);
    // l_ul2->SetLineColor(28);
    // l_ul2->Draw("same");
    //  +++++++++++++++++++

    fprintf(table_ul_yphifo, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f & %0.2f $\\pm$ %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), N_y_gen, N_y_mc, dN_y_mc, N_y_data, dN_y_data, eff_y * 100, xsec_y, dxsec_y_stat, Z);    

    // fprintf(table_ul_yphi2pi_sys, "%0.2f - %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", h_beame_tru->GetXaxis()->GetXmin(), h_beame_tru->GetXaxis()->GetXmax(), abs(err_bkg_poln*100/bkg_poln[0]), abs(err_fit_range*100/fit_range[0]), abs(err_slice_numb*100/slice_numb[0]), abs(err_mean_value*100/mean_value[0]), abs(err_width_value*100/width_value[0]));    

    // int j =1;
    // gphiy->Write(Form("grphiy_%d", j), TObject::kWriteDelete);

    // cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.root", name.Data()), "root");
    // cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.eps", name.Data()), "eps");
    // cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.png", name.Data()), "png");
    // cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.root", name.Data()), "root");
    // cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.eps", name.Data()), "eps");
    // cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.png", name.Data()), "png");
    // cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.root", name.Data()), "root");
    // cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.eps", name.Data()), "eps");
    // cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.png", name.Data()), "png");
    // cgnophiy->Print("output/fig_ul_yphifo/c_gnophiy.root", "root");
    // cgnophiy->Print("output/fig_ul_yphifo/c_gnophiy.eps", "eps");

    cgphiy->Print(Form("output/fig_ul_yphifo/c_gphiy_%s.root", name.Data()), "root");
    cgphiy->Print(Form("output/fig_ul_yphifo/c_gphiy_%s.eps", name.Data()), "eps");
    cgphiy->Print(Form("output/fig_ul_yphifo/c_gphiy_%s.png", name.Data()), "png");  
    chgphiy_sub->Print(Form("output/fig_ul_yphifo/c_hgphiy_sub_%s.root", name.Data()), "root");
    chgphiy_sub->Print(Form("output/fig_ul_yphifo/c_hgphiy_sub_%s.eps", name.Data()), "eps");
    chgphiy_sub->Print(Form("output/fig_ul_yphifo/c_hgphiy_sub_%s.png", name.Data()), "png");//cgphifo
    cgphifo->Print(Form("output/fig_ul_yphifo/c_gphifo_sideband_%s.root", name.Data()), "root");
    cgphifo->Print(Form("output/fig_ul_yphifo/c_gphifo_sideband_%s.eps", name.Data()), "eps");
    cgphifo->Print(Form("output/fig_ul_yphifo/c_gphifo_sideband_%s.png", name.Data()), "png");//
    // cprof->Print(Form("output/fig_ul_yphifo/c_prof_%s.root", name.Data()), "root");
    // cprof->Print(Form("output/fig_ul_yphifo/c_prof_%s.eps", name.Data()), "eps");
    // cprof->Print(Form("output/fig_ul_yphifo/c_prof_%s.png", name.Data()), "png");
    // cprof->Print(Form("output/fig_ul_yphifo/c_prof_%s.pdf", name.Data()), "pdf");
    // cprofxsec->Print(Form("output/fig_ul_yphifo/c_profxsec_%s.root", name.Data()), "root");
    // cprofxsec->Print(Form("output/fig_ul_yphifo/c_profxsec_%s.eps", name.Data()), "eps");
    // cprofxsec->Print(Form("output/fig_ul_yphifo/c_profxsec_%s.png", name.Data()), "png");
    // cprofxsec->Print(Form("output/fig_ul_yphifo/c_profxsec_%s.pdf", name.Data()), "pdf");

    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.png", name.Data()), "png");

    cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.root", name.Data()), "root");
    cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.eps", name.Data()), "eps");
    cphiy1->Print(Form("output/fig_ul_yphifo/c1_phiy_%s.png", name.Data()), "png");
    cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.root", name.Data()), "root");
    cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.eps", name.Data()), "eps");
    cphiy2->Print(Form("output/fig_ul_yphifo/c2_phiy_%s.png", name.Data()), "png");
    cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.root", name.Data()), "root");
    cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.eps", name.Data()), "eps");
    cphiy->Print(Form("output/fig_ul_yphifo/c_phiy_%s.png", name.Data()), "png");
    cphiy_nofo->Print(Form("output/fig_ul_yphifo/c_phiy_nofo_%s.root", name.Data()), "root");
    cphiy_nofo->Print(Form("output/fig_ul_yphifo/c_phiy_nofo_%s.eps", name.Data()), "eps");
    cphiy_nofo->Print(Form("output/fig_ul_yphifo/c_phiy_nofo_%s.png", name.Data()), "png");
    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_ul_yphifo/c%s_tagged_flux.png", name.Data()), "png");
    c_beame_tru->Print(Form("output/fig_ul_yphifo/c_beame_tru_%s.root", name.Data()), "root");
    c_beame_tru->Print(Form("output/fig_ul_yphifo/c_beame_tru_%s.eps", name.Data()), "eps");
    c_beame_tru->Print(Form("output/fig_ul_yphifo/c_beame_tru_%s.png", name.Data()), "png");
    chgphiy_sub2->Print(Form("output/fig_ul_yphifo/chgphiy_sub2_%s.root", name.Data()), "root");
    chgphiy_sub2->Print(Form("output/fig_ul_yphifo/chgphiy_sub2_%s.eps", name.Data()), "eps");
    chgphiy_sub2->Print(Form("output/fig_ul_yphifo/chgphiy_sub2_%s.png", name.Data()), "png");
    chgphiy->Print(Form("output/fig_ul_yphifo/chgphiy_%s.root", name.Data()), "root");
    chgphiy->Print(Form("output/fig_ul_yphifo/chgphiy_%s.eps", name.Data()), "eps");
    chgphiy->Print(Form("output/fig_ul_yphifo/chgphiy_%s.png", name.Data()), "png");
    chgphiy_nofo->Print(Form("output/fig_ul_yphifo/chgphiy_nofo_%s.root", name.Data()), "root");
    chgphiy_nofo->Print(Form("output/fig_ul_yphifo/chgphiy_nofo_%s.eps", name.Data()), "eps");
    chgphiy_nofo->Print(Form("output/fig_ul_yphifo/chgphiy_nofo_%s.png", name.Data()), "png");
    chgphifo->Print(Form("output/fig_ul_yphifo/chgphifo_%s.root", name.Data()), "root");
    chgphifo->Print(Form("output/fig_ul_yphifo/chgphifo_%s.eps", name.Data()), "eps");
    chgphifo->Print(Form("output/fig_ul_yphifo/chgphifo_%s.png", name.Data()), "png");

    // outputfig->Write();
    outputfig->Print();
    // outputfig->Close();

    fprintf(table_ul_yphifo, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_ul_yphifo);
    gSystem->Exec(Form("pdflatex table_ul_yphifo_%s.tex", name.Data()));

    // fprintf(table_ul_yphifo_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    // fclose(table_ul_yphifo_sys);
    // gSystem->Exec(Form("pdflatex table_ul_yphifo_sys_%s.tex", name.Data()));
  
}
