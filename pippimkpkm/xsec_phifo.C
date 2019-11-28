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

void xsec_phifo(TString name, int n2k=100, int n2pi=50, int ne=1, int nt=1)//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fdata = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_%s_flat.root", name.Data()));
  TFile *fmc = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_phifo_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/thrown_phifo_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_11366_11555.root");
  if(name == "17") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  if(name == "18") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_40856_42577.root");

  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/xsec_phifo.root", "UPDATE");
  
  // ofstream table_phi("table_phi.tex", ofstream::out);
  // table_phi <<"\\documentclass[8pt]{extarticle}" << endl;
  // table_phi <<"\\usepackage[margin=0.1in]{geometry}" << endl; 
  // table_phi <<"\\usepackage{tabularx}" << endl;
  // table_phi <<"\\begin{document}" << endl;
  // table_phi <<"\\begin{table}[!htbp]" << endl;
  // table_phi <<"\\begin{minipage}{\\textwidth}" << endl;
  // table_phi <<"\\begin{center}" << endl;
  // table_phi <<"\\caption{Total cross-sections in photon energy bins} "<< endl;
  // // table_phi <<"\\label{tab : table1}" << endl;
  // table_phi <<"\\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}" << endl;
  // table_phi <<"\\hline" << endl;
  // table_phi <<"$E_{\\gamma}[ GeV^{2}]$ & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [\\%]& $\\sigma [nb]$ \\\\" << endl;
  // // table_phi << "\\alpha & \\beta & \\gamma \\" << endl;
  // table_phi <<"\\hline" << endl;

  ofstream ofs_phifo("phifo_tab.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.98, mkk_max = 1.2; //mkk_min = 0.98, mkk_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // ##############################################  fo ALL ##############################################

  FILE *table_xsec_phifo = fopen("table_xsec_phifo.tex","w");
  fprintf(table_xsec_phifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Total cross-sections}\n \\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma \\times BR_{f_{0}\\rightarrow\\pi^{+}\\pi^{-}}$ [nb] \\\\ \n \\hline\n");

  // FILE *table_xsec_phifo_sys = fopen("table_xsec_phifo_sys.tex","w");
  // fprintf(table_xsec_phifo_sys,"\\documentclass[12pt]{extarticle}\n \\usepackage{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabularx}{\\textwidth}{|X|X|X|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ \\% ] & Fit range [ \\% ]  \\\\ \n \\hline\n");  

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  // h_tagged_flux->Rebin(100);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
   TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  // h_beame_tru->Rebin(100);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
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

  // *********************** f0(980) MC *********************************
  TCanvas *c_foMass_postcuts = new TCanvas("c_foMass_postcuts", "c_foMass_postcuts", 1000, 600);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", "MC signal; m_{#pi^{+}#pi^{-}} [GeV/c^{2}]; Counts", 200, 0.85, 1.05);
  tmc->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni)");//+cutlist+" && kin_chisq<25)"
  c_foMass_postcuts->cd();
  h_foMass_postcuts->Draw("e");

  // // w.factory("BreitWigner::sig_foMass_mc(m_foMass_mc[1.005, 1.035],mean_foMass_mc[1.015,1.022],width_foMass_mc[0.004])");
  // // w.factory("ArgusBG::bkg_fo_mc(m_fo_mc, 1.04, c0_fo_mc[-50,-10])");
  // w.factory("Voigtian::sig_foMass_mc(m_foMass_mc[1.005, 1.035],mean_foMass_mc[1.015,1.022],width_foMass_mc[0.004],sigma_foMass_mc[0.001,0.01])"); //sigma_foMass_mc[0.001,0.01], mean_foMass_mc[1.015,1.022]
  // w.factory("Chebychev::bkg_foMass_mc(m_foMass_mc,{c0_foMass_mc[-10,10], c1_foMass_mc[-10,10]})");
  // w.factory("SUM:::model_foMass_mc(nsig_foMass_mc[0,100000000]*sig_foMass_mc, nbkg_foMass_mc[0,100000000]*bkg_foMass_mc)");//, nbkg_foMass_mc[0,100000000]*bkg_foMass_mc)"); //nsig[0,100000000]*sig2,
  // w.var("m_foMass_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
  // RooDataHist dh_foMass_mc("dh_foMass_mc", "dh_foMass_mc", *w.var("m_foMass_mc"), Import(*h_foMass_postcuts));
  // RooPlot *fr_foMass_mc = w.var("m_foMass_mc")->frame(Title("K^{+}K^{-}"));
  // // fr_foMass_mc->SetTitleOffset(0.90, "X");
  // // fr_foMass_mc->SetTitleSize(0.06, "XYZ");
  // // fr_foMass_mc->SetLabelSize(0.06, "xyz");
  // w.pdf("model_foMass_mc")->fitTo(dh_foMass_mc);

  // // //result = w.pdf("model")->fitTo(dh_foMass,Extended(kTRUE),Save());
  // dh_foMass_mc.plotOn(fr_foMass_mc, RooFit::Name("ldh_foMass_mc"));
  // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("sig_foMass_mc")), LineColor(kRed), RooFit::Name("lsig_foMass_mc"));
  // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("bkg_foMass_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_foMass_mc"));
  // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, RooFit::Name("lmodel_foMass_mc"));
  // // w.pdf("model_foMass_mc")->paramOn(fr_foMass_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_foMass_mc"), *w.var("nbkg_foMass_mc")))); //,*w.var("mean_foMass_mc"),*w.var("width_foMass_mc"),*w.var("sigma_foMass_mc"))));
  // fr_foMass_mc->Draw();

  // TLegend *l_fo_mc = new TLegend(0.2, 0.65, 0.4, 0.85);
  // l_fo_mc->SetFillColor(kWhite);
  // l_fo_mc->SetLineColor(kWhite);
  // // l_fo_mc->AddEntry(fr_foMass_mc->findObject("ldh_foMass_mc"), "Data", "p");
  // l_fo_mc->AddEntry(fr_foMass_mc->findObject("lmodel_foMass_mc"), "total", "l");
  // l_fo_mc->AddEntry(fr_foMass_mc->findObject("lsig_foMass_mc"), "Voigtian", "l");
  // l_fo_mc->AddEntry(fr_foMass_mc->findObject("lbkg_foMass_mc"), "pol 2nd", "l");
  // l_fo_mc->Draw();

  // TLatex lat_foMass_mc;
  // lat_foMass_mc.SetTextSize(0.05);
  // lat_foMass_mc.SetTextAlign(13); //align at top
  // lat_foMass_mc.SetNDC();
  // lat_foMass_mc.SetTextColor(kBlue);
  // lat_foMass_mc.DrawLatex(0.62, 0.87, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_foMass_mc")->getVal(), w.var("nsig_foMass_mc")->getError()));
  // lat_foMass_mc.DrawLatex(0.62, 0.78, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_foMass_mc")->getVal(), w.var("nbkg_foMass_mc")->getError()));
  // lat_foMass_mc.DrawLatex(0.62, 0.68, Form("#mu = %0.3f#pm%0.3f", w.var("mean_foMass_mc")->getVal(), w.var("mean_foMass_mc")->getError()));
  // lat_foMass_mc.DrawLatex(0.62, 0.58, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_foMass_mc")->getVal(), w.var("width_foMass_mc")->getError()));
  // lat_foMass_mc.DrawLatex(0.62, 0.48, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_foMass_mc")->getVal(), w.var("sigma_foMass_mc")->getError()));

   TF1 *fsb_fo_mc = new TF1("fsb_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", 0.85, 1.05);
   fsb_fo_mc->SetLineColor(2);
   fsb_fo_mc->SetParameters(1800, 0.988, 0.06, 1,1,1,1); //1800, 0.988, 0.06, 1,1,1,1
  //  fsb_fo_mc->SetParLimits(0, 0, 10000);
   fsb_fo_mc->SetParLimits(1, 0.95, 0.1); // 1.018, 1.021
   fsb_fo_mc->SetParLimits(2, 0.001, 0.6); //fsb->SetParLimits(2, 0.008, 0.010);
   // fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

   // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
   TF1 *fs_fo_mc = new TF1("fs_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2])", 0.85, 1.05);
   fs_fo_mc->SetLineColor(4);

   TF1 *fb_fo_mc = new TF1("fb_fo_mc", "pol3(3)", 0.85, 1.05); //pol2(3)
   fb_fo_mc->SetLineColor(28);
   fb_fo_mc->SetLineStyle(2);

   h_foMass_postcuts->Fit("fsb_fo_mc", "", "", 0.85, 1.05);
   double par_fo_mc[10]; //6
   fsb_fo_mc->GetParameters(&par_fo_mc[0]);
   fs_fo_mc->SetParameters(&par_fo_mc[0]);
   fb_fo_mc->SetParameters(&par_fo_mc[3]); //3

   double N_fo_mc = fs_fo_mc->Integral(0.85, 1.05) * n2pi / (m2pi_max - m2pi_min);
   double dN_fo_mc = N_fo_mc * fsb_fo_mc->GetParError(0) / fsb_fo_mc->GetParameter(0);

   fs_fo_mc->Draw("same");
   fb_fo_mc->Draw("same");

   TLatex lat_fo_mc;
   lat_fo_mc.SetTextSize(0.05);
   lat_fo_mc.SetTextAlign(13); //align at top
   lat_fo_mc.SetNDC();
   lat_fo_mc.SetTextColor(kBlue);
   lat_fo_mc.DrawLatex(0.2, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_fo_mc->GetChisquare() / fsb_fo_mc->GetNDF()));
   lat_fo_mc.DrawLatex(0.2, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N_fo_mc, dN_fo_mc));
   lat_fo_mc.DrawLatex(0.2, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_fo_mc->GetParameter(1), fsb_fo_mc->GetParError(1)));
   lat_fo_mc.DrawLatex(0.2, 0.58, Form("#Gamma = %0.3f#pm%0.3f", fsb_fo_mc->GetParameter(2), fsb_fo_mc->GetParError(2)));

   // ======================================== fo vs. Phi(1020) ===============================================
   // root -l 'xsec_phifo.C+(100,50,1,1)'

   TCanvas *cphifo = new TCanvas("cphifo", "cphifo", 900, 600); // 900, 600
   TCanvas *cphifo1 = new TCanvas("cphifo1", "cphifo1", 1500, 800);
   cphifo1->Divide(5, 5);
   TCanvas *cphifo2 = new TCanvas("cphifo2", "cphifo2", 1500, 800);
   cphifo2->Divide(5, 5);
   TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 900, 600);
   TGraphErrors *gphifo;

   // TCanvas *cgphifo_width = new TCanvas("cgphifo_width", "cgphifo_width", 1500, 400);//
   // cgphifo_width->Divide(3,1);
   // TGraphErrors *gphifo_width[ne];

   // TCanvas *cgphifo_mean = new TCanvas("cgphifo_mean", "cgphifo_mean", 1500, 400);//
   // cgphifo_mean->Divide(3,1);
   // TGraphErrors *gphifo_mean[ne];

   // TCanvas *cgnophifo=new TCanvas("cgnophifo","cgnophifo",1500, 800);
   // TGraphErrors *gnophifo = new TGraphErrors(n2pi);
   // gnophifo->SetMarkerStyle(20);
   // gnophifo->SetMarkerSize(1.0);
   // gnophifo->SetMarkerColor(2);

   // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

   cphifo->cd();
   TH2F *h2d2 = new TH2F("h2d2", "(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, n2k, mkk_min, mkk_max);
   tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>5.9 && beam_p4_kin.E()<11.9)");
   h2d2->Draw("colz");

   gphifo = new TGraphErrors(); //n2pi
   gphifo->SetMarkerStyle(20);
   gphifo->SetMarkerSize(1.0);
   gphifo->SetMarkerColor(1);
   gphifo->SetMinimum(0.);
   gphifo->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

   // gphifo_width = new TGraphErrors(n2pi);
   // gphifo_width->SetMarkerStyle(20);
   // gphifo_width->SetMarkerSize(1.0);
   // gphifo_width->SetMarkerColor(1);
   // gphifo_width->SetMinimum(0.);
   // gphifo_width->SetMaximum(0.01);
   // gphifo->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

   // gphifo_mean = new TGraphErrors(n2pi);
   // gphifo_mean->SetMarkerStyle(20);
   // gphifo_mean->SetMarkerSize(1.0);
   // gphifo_mean->SetMarkerColor(1);
   // gphifo_mean->SetMinimum(0.);
   // gphifo_mean->SetMaximum(0.01);
   // gphifo->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

   // gnophifo->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{bkg}", Eg1[j], Eg2[j]));

   for (int i = 1; i <= n2pi; ++i) //n2pi
   {
     // cout << i << " " << flush;

     if (i < 26)
       cphifo1->cd(i); // if (i < 26) cphifo1->cd(i);
     if (i > 25)
       cphifo2->cd(i - 25); // if (i > 25 && i < 51) cphifo2->cd(i - 25);

     TH1D *hphifo_py = h2d2->ProjectionY(Form("_hphifo_py_%d", i), i, i);

     // hslice2->Fit("fsb","q","",0.99,1.08);
     // hslice2->Fit("fsb","qm","",0.99,1.12);
     // fs->SetParameters(fsb->GetParameters());
     // fb2->SetParameters(fsb->GetParameters());
     hphifo_py->Draw();
     // fb2->Draw("same");

     TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
     // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
     fsb->SetLineColor(2);
     fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
     fsb->SetParLimits(0, 0, 10000);
     fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
     fsb->FixParameter(2, 0.004);        //fsb->SetParLimits(2, 0.008, 0.010);
     fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

     TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
     // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
     fs->SetLineColor(4);

     TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max); //pol2(3)
     fb->SetLineColor(28);
     fb->SetLineStyle(2);

     hphifo_py->Fit("fsb", "", "", mkk_min, mkk_max);
     double par[7]; //6
     fsb->GetParameters(&par[0]);
     fs->SetParameters(&par[0]);
     fb->SetParameters(&par[4]); //3

     fs->Draw("same");
     fb->Draw("same");

     double N2 = fs->Integral(mkk_min, mkk_max) / hphifo_py->GetBinWidth(1);
     double dN2 = N2 * fsb->GetParError(0) / fsb->GetParameter(0);

     gphifo->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
     gphifo->SetPointError(i - 1, 0, dN2);

     // gphifo_width[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(2));
     // gphifo_width[j]->SetPointError(i - 1, 0, fsb->GetParError(2));

     // gphifo_mean[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(3));
     // gphifo_mean[j]->SetPointError(i - 1, 0, fsb->GetParError(3));

     TLatex lat_phifo;
     lat_phifo.SetTextSize(0.09);
     lat_phifo.SetTextAlign(13); //align at top
     lat_phifo.SetNDC();
     lat_phifo.SetTextColor(kBlue);
     lat_phifo.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
     lat_phifo.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
     lat_phifo.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
     lat_phifo.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
     lat_phifo.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

     // gnophifo->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), Nbkg);
     // gnophifo->SetPointError(i - 1, 0, dNbkg);
     // ofs_xsec_phifo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h2d2->GetYaxis()->GetBinCenter(" << i << ") = " << h2d2->GetYaxis()->GetBinCenter(i)<<endl;
     // ofs_xsec_phifo << " i = " << i << " |par0 = " << fsb->GetParameter(0) << " | parerr0 = " << fsb->GetParError(0) << " | par1 " << fsb->GetParameter(1) << " | par2 " << fsb->GetParameter(2) <<endl;
     hphifo_py->Write();
     cgphifo->Update();
     // c2->Update();
     //sleep(1);
    }

    // cout << endl;

    cgphifo->cd();
    gphifo->Draw("AP");
    // TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
    TF1 *fsb_fo_data = new TF1("fsb_fo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", 0.85, 1.05);
    fsb_fo_data->SetLineColor(2);
    fsb_fo_data->SetParameters(1433, 0.980, 0.1, 1,1,1,1);
    // fsb->SetParLimits(0, 0, 10000);
    fsb_fo_data->SetParLimits(1, 0.95, 0.98); 
    fsb_fo_data->SetParLimits(2, 0.001, 0.6); 
    // fsb_fo_data->FixParameter(1, fsb_fo_mc->GetParameter(1)); 
    // fsb_fo_data->FixParameter(2, fsb_fo_mc->GetParameter(2));

    // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    TF1 *fs_fo_data = new TF1("fs_fo_data", "[0]*TMath::BreitWigner(x,[1],[2])", 0.85, 1.05);
    fs_fo_data->SetLineColor(4);

    TF1 *fb_fo_data = new TF1("fb", "pol3(3)", 0.85, 1.05); //pol2(3)
    fb_fo_data->SetLineColor(28);
    fb_fo_data->SetLineStyle(2);

    gphifo->Fit("fsb_fo_data", "", "", 0.85, 1.05);
    double par_fo_data[10]; //6
    fsb_fo_data->GetParameters(&par_fo_data[0]);
    fs_fo_data->SetParameters(&par_fo_data[0]);
    fb_fo_data->SetParameters(&par_fo_data[3]); //3

    double N_fo_data = fs_fo_data->Integral(0.85, 1.05)*n2pi/(m2pi_max-m2pi_min);
    double dN_fo_data = N_fo_data * fsb_fo_data->GetParError(0) / fsb_fo_data->GetParameter(0); 

    fs_fo_data->Draw("same");
    fb_fo_data->Draw("same");

    TLatex lat_phifo;
    lat_phifo.SetTextSize(0.05);
    lat_phifo.SetTextAlign(13); //align at top
    lat_phifo.SetNDC();
    lat_phifo.SetTextColor(kBlue);
    lat_phifo.DrawLatex(0.6, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_fo_data->GetChisquare() / fsb_fo_data->GetNDF()));
    lat_phifo.DrawLatex(0.6, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N_fo_data, dN_fo_data));
    lat_phifo.DrawLatex(0.6, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(1), fsb_fo_data->GetParError(1)));
    lat_phifo.DrawLatex(0.6, 0.58, Form("#Gamma = %0.3f#pm%0.3f", fsb_fo_data->GetParameter(2), fsb_fo_data->GetParError(2)));

   // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phifo = N_fo_mc / h_beame_tru->Integral(0,100); // Efficiency = N_observed/N_generated
    double deff_phifo = eff_phifo * (dN_fo_data / N_fo_mc);

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phifo = h_tagged_flux->Integral(0,100); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    // double br_fo = ?;
    // double dbr_fo = ?;
    double xsec_phifo = 1e9 * N_fo_data / (eff_phifo * lumi_phifo * T * br_phi);
    double dxsec_phifo_stat = xsec_phifo * (dN_fo_data / N_fo_data);
    double dxsec_phifo_sys = xsec_phifo * TMath::Sqrt((deff_phifo / eff_phifo) * (deff_phifo / eff_phifo) + (dbr_phi / br_phi) * (dbr_phi / br_phi)); //(err_sys / N_fo_data) * (err_sys / N_fo_data) + 

    fprintf(table_xsec_phifo, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h_tagged_flux->GetXaxis()->GetXmin(), h_tagged_flux->GetXaxis()->GetXmax(), h_beame_tru->Integral(0,100), N_fo_mc, dN_fo_mc, N_fo_data, dN_fo_data, eff_phifo * 100, deff_phifo * 100, xsec_phifo, dxsec_phifo_stat, dxsec_phifo_sys);

    // fprintf(table_xsec_phifo_sys, "%0.2f & %0.2f & %0.2f \\\\ \n", h2phie_mc->GetXaxis()->GetBinCenter(i), err_bkg_poln * 100 / N_bkg_poln_1[i - 1], err_fit_range * 100 / fit_range_1[i - 1]);
    // table_phi << std::setprecision(2) << std::fixed;
    // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_fo_mc<<"$\\pm$"<<dN_fo_mc<<"&"<<N_fo_data<<"$\\pm$"<<dN_fo_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phifo<<"$\\pm$"<<dxsec_phifo_stat<<"$\\pm$"<<dxsec_phifo_sys<<" \\\\"<< endl;  


   c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.root", name.Data()), "root");
   c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.eps", name.Data()), "eps");
   c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.png", name.Data()), "png");
   c_foMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.root", name.Data()), "root");
   c_foMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.eps", name.Data()), "eps");
   c_foMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.png", name.Data()), "png");
   cphifo1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c1_phifo_%s.root", name.Data()), "root");
   cphifo1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c1_phifo_%s.eps", name.Data()), "eps");
   cphifo1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c1_phifo_%s.png", name.Data()), "png");
   cphifo2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c2_phifo_%s.root", name.Data()), "root");
   cphifo2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c2_phifo_%s.eps", name.Data()), "eps");
   cphifo2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c2_phifo_%s.png", name.Data()), "png");
   cphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_phifo_%s.root", name.Data()), "root");
   cphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_phifo_%s.eps", name.Data()), "eps");
   cphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_phifo_%s.png", name.Data()), "png");
   cgphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_%s.root", name.Data()), "root");
   cgphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_%s.eps", name.Data()), "eps");
   cgphifo->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_%s.png", name.Data()), "png");

   c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_tagged_flux.root", name.Data()), "root");
   c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_tagged_flux.eps", name.Data()), "eps");
   c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_tagged_flux.png", name.Data()), "png");
   c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_beame_tru_%s.root", name.Data()), "root");
   c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_beame_tru_%s.eps", name.Data()), "eps");
   c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_beame_tru_%s.png", name.Data()), "png");

   outputfig->Print();

   fprintf(table_xsec_phifo, "\\hline\n \\end{tabularx}\n \\end{center}\n \\end{minipage}\n \\end{table}\n \\end{document}\n");
   fclose(table_xsec_phifo);
   gSystem->Exec("pdflatex table_xsec_phifo.tex");

  //  fprintf(table_xsec_phifo_sys, "\\hline\n \\end{tabularx}\n \\end{center}\n \\end{minipage}\n \\end{table}\n \\end{document}\n");
  //  fclose(table_xsec_phifo_sys);
  //  gSystem->Exec("pdflatex table_xsec_phifo_sys.tex");

   // cgphifo_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_width.root", "root");
   // cgphifo_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_width.eps", "eps");
   // cgphifo_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_mean.root", "root");
   // cgphifo_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gphifo_mean.eps", "eps");

   // cgnophifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gnophifo.root", "root");
   // cgnophifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c_gnophifo.eps", "eps");

   /*
  // ========================= Phi vs. Eg ======================

  //root -l 'xsec_phifo.C+("data_17",100,10,1)'

  FILE *table_xsec_phifo = fopen("table_xsec_phifo.tex","w");
  fprintf(table_xsec_phifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Total cross-sections in photon energy bins}\n \\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma$ [nb] \\\\ \n \\hline\n");

  FILE *table_xsec_phifo_sys = fopen("table_xsec_phifo_sys.tex","w");
  fprintf(table_xsec_phifo_sys,"\\documentclass[12pt]{extarticle}\n \\usepackage{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabularx}{\\textwidth}{|X|X|X|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ \\% ] & Fit range [ \\% ]  \\\\ \n \\hline\n");  

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(100);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
   TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown");// new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  // ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  h_beame_tru->Rebin(100);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");

  // Data --- Phi(1020)
  TCanvas *cphie = new TCanvas("cphie", "cphie", 900, 600);
  cphie->cd();
  TH2D *h2phie = new TH2D("h2phie", "Data; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, 5.9, 11.9, n2k, mkk_min, mkk_max);
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");

  h2phie->Draw("colz");

  TCanvas *cphie1 = new TCanvas("cphie1", "cphie1", 1500, 600);
  cphie1->Divide(3, 1);
  TCanvas *cgphie = new TCanvas("cgphie", "cgphie", 900, 600);
  TGraphErrors *gphie = new TGraphErrors();
  gphie->SetMarkerStyle(20);
  gphie->SetMarkerSize(1.5);
  gphie->SetMarkerColor(1);
  gphie->SetMinimum(0.0);
  gphie->SetTitle("#phi(1020) Yield vs. E_{#gamma} (data); E_{#gamma} [GeV]; N_{#phi}");

  // Monte Carlo --- Phi(1020)
  TCanvas *cphie_mc = new TCanvas("cphie_mc", "cphie_mc", 900, 600);
  cphie_mc->cd();
  TH2D *h2phie_mc = new TH2D("h2phie_mc", "MC; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, 5.9, 11.9, n2k, mkk_min, mkk_max);//0.98, 1.2
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");
  h2phie_mc->Draw("colz");

  TCanvas *cphie1_mc = new TCanvas("cphie1_mc", "cphie1_mc", 1500, 600);
  cphie1_mc->Divide(3, 1);
  TCanvas *cgphie_mc = new TCanvas("cgphie_mc", "cgphie_mc", 900, 600);
  TGraphErrors *gphie_mc = new TGraphErrors();
  gphie_mc->SetMarkerStyle(20);
  gphie_mc->SetMarkerSize(1.5);
  gphie_mc->SetMarkerColor(1);
  gphie_mc->SetMinimum(0.0);
  gphie_mc->SetTitle("#phi(1020) Yield vs. E_{#gamma} (MC); E_{#gamma} [GeV]; N_{#phi}");

  // Monte Carlo --- fo(980)
  TCanvas *cfoe_mc = new TCanvas("cfoe_mc", "cfoe_mc", 900, 600);
  cfoe_mc->cd();
  TH2D *h2foe_mc = new TH2D("h2foe_mc", "MC; E_{#gamma} [GeV]; m_{#pi^{+}#pi^{-}} [GeV/c^{2}]", ne, 5.9, 11.9, n2k, m2pi_min, m2pi_max);//0.98, 1.2
  tmc->Project("h2foe_mc", "pippim_mf:beam_p4_kin.E()", "w8*(pippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035)");
  h2foe_mc->Draw("colz");

  TCanvas *cfoe1_mc = new TCanvas("cfoe1_mc", "cfoe1_mc", 1500, 600);
  cfoe1_mc->Divide(3, 1);
  TCanvas *cgfoe_mc = new TCanvas("cgfoe_mc", "cgfoe_mc", 900, 600);
  TGraphErrors *gfoe_mc = new TGraphErrors();
  gfoe_mc->SetMarkerStyle(20);
  gfoe_mc->SetMarkerSize(1.5);
  gfoe_mc->SetMarkerColor(1);
  gfoe_mc->SetMinimum(0.0);
  gfoe_mc->SetTitle("f_{0}(980) Yield vs. E_{#gamma} (MC); E_{#gamma} [GeV]; N_{f_{0}}");

  // ===================== fo vs. Phi(1020) ========================
  TCanvas *cphifo = new TCanvas("cphifo", "cphifo", 1500, 600); // 900, 600
  cphifo->Divide(3, 1);
  TCanvas *cphifo1 = new TCanvas("cphifo1", "cphifo1", 1500, 800);
  cphifo1->Divide(5, 5);
  TCanvas *cphifo2 = new TCanvas("cphifo2", "cphifo2", 1500, 800);
  cphifo2->Divide(5, 5);
  TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 1500, 600);
  cgphifo->Divide(3, 1);
  TGraphErrors *gphifo;

  // Efficiency
  TCanvas *cgphifoeeff = new TCanvas("cgphifoeeff", "cgphifoeeff", 900, 600);
  TGraphErrors *gphifoeeff = new TGraphErrors();
  gphifoeeff->SetMarkerStyle(100);
  gphifoeeff->SetMarkerSize(1.5);
  gphifoeeff->SetMarkerColor(kBlue);
  gphifoeeff->SetMinimum(0.0);
  gphifoeeff->SetTitle("#phi f_{0} Effeciency vs. E_{#gamma}; E_{#gamma} [GeV]; #epsilon [%]");

  // Cross-section
  TCanvas *cgphifoexsec = new TCanvas("cgphifoexsec", "cgphifoexsec", 900, 600);
  TGraphErrors *gphifoexsec = new TGraphErrors();
  gphifoexsec->SetMarkerStyle(100);
  gphifoexsec->SetMarkerSize(1.5);
  gphifoexsec->SetMarkerColor(1);
  gphifoexsec->SetMinimum(0.0);
  gphifoexsec->SetTitle("#phi f_{0} total Cross-section vs. E_{#gamma}; E_{#gamma} [GeV]; #sigma [nb]"); //#phi(1020) flux normalized yield, Yield_{#phi}

  double Egmin = 5.9;  //6.3;//= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 11.9; //11.6;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax - Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];


  for (int i = 1; i <= ne; ++i)
  {
    cout << i << " " << flush;
    Eg1[i] = Egmin + ((i - 1) * Egstep); //i * Egstep;
    Eg2[i] = Egmin + (i * Egstep);

    // ++++++++++++++++++++++++++++ mc --- Phi(1020) +++++++++++++++++++++++
    cphie1_mc->cd(i);
    TH1D *hphie_mc_py = h2phie_mc->ProjectionY(Form("hphie_mc_py_%d", i), i, i);
    hphie_mc_py->SetTitle(Form("E_{#gamma} = %.2f [GeV]", h2phie_mc->GetXaxis()->GetBinCenter(i)));
    hphie_mc_py->Draw("e");

    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1, 1, 1, 1, 1); // voigt
    fsb_mc->SetParLimits(1, 1.015, 1.022);                            //fsb_mc->FixParameter(1, 1.019);
    // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_mc->FixParameter(3, 0.0042); //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max);
    fb_mc->SetLineColor(28);
    fb_mc->SetLineStyle(2);

    hphie_mc_py->Fit("fsb_mc", "", "", mkk_min, mkk_max);
    double par_mc[10]; //6
    fsb_mc->GetParameters(&par_mc[0]);
    fs_mc->SetParameters(&par_mc[0]);
    fb_mc->SetParameters(&par_mc[4]); //4

    fs_mc->Draw("same");
    fb_mc->Draw("same");

    double N_phie_mc = fs_mc->Integral(mkk_min, mkk_max) / hphie_mc_py->GetBinWidth(1);
    double dN_phie_mc = N_phie_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

    TLatex lat_phie_mc;
    lat_phie_mc.SetTextSize(0.06);
    lat_phie_mc.SetTextAlign(13); //align at top
    lat_phie_mc.SetNDC();
    lat_phie_mc.SetTextColor(kBlue);
    lat_phie_mc.DrawLatex(0.45, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
    lat_phie_mc.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phie_mc, dN_phie_mc));
    lat_phie_mc.DrawLatex(0.45, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
    lat_phie_mc.DrawLatex(0.45, 0.66, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phie_mc.DrawLatex(0.45, 0.59, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

    gphie_mc->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), N_phie_mc);
    gphie_mc->SetPointError(i - 1, 0, dN_phie_mc);

    // ++++++++++++++++++++++++++++ mc --- fo(980) +++++++++++++++++++++++
    cfoe1_mc->cd(i);
    TH1D *hfoe_mc_py = h2foe_mc->ProjectionY(Form("hfoe_mc_py_%d", i), i, i);
    hfoe_mc_py->SetTitle(Form("E_{#gamma} = %.2f [GeV]", h2foe_mc->GetXaxis()->GetBinCenter(i)));
    hfoe_mc_py->Draw("e");

    TF1 *fsb_foe_mc = new TF1("fsb_foe_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.8, 1.15); // + pol2(4) , [2, 2.5]
    // TF1 *fsb_foe_mc = new TF1("fsb_foe_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", 0.8, 1.15);
    fsb_foe_mc->SetLineColor(2);
    fsb_foe_mc->SetParameters(500, 0.98, 0.01, 0.06, 1,1,1); // sigma=0.01,
    // // fsb_foe_mc->SetParLimits(0, 0, 10000);
    // fsb_foe_mc->SetParLimits(1, 0.9, 1.0); // 1.018, 1.021
    // fsb_foe_mc->SetParLimits(2, 0.008, 0.05);
    // // fsb_foe_mc->FixParameter(2, 0.079);
    // fsb_foe_mc->SetParLimits(3, 0.03, 0.08);
    // // fsb_foe_mc->FixParameter(3, 0.05);

    // TF1 *fs_foe_mc = new TF1("fs_foe_mc", "[0]*TMath::BreitWigner(x,[1],[2])", 0.8, 1.15);
    TF1 *fs_foe_mc = new TF1("fs_foe_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.8, 1.15);
    fs_foe_mc->SetLineColor(4);

    TF1 *fb_foe_mc = new TF1("fb_foe_mc", "pol2(4)", 0.8, 1.15); //pol2(3)
    fb_foe_mc->SetLineColor(28);
    fb_foe_mc->SetLineStyle(2);

    hfoe_mc_py->Fit("fsb_foe_mc", "", "", 0.8, 1.15);
    double par_foe_mc[8]; //6
    fsb_foe_mc->GetParameters(&par_foe_mc[0]);
    fs_foe_mc->SetParameters(&par_foe_mc[0]);
    fb_foe_mc->SetParameters(&par_foe_mc[4]); //3

    fs_foe_mc->Draw("same");
    fb_foe_mc->Draw("same");

    double N_foe_mc = fs_foe_mc->Integral(0.8, 1.15) / hfoe_mc_py->GetBinWidth(1);
    double dN_foe_mc = N_foe_mc * fsb_foe_mc->GetParError(0) / fsb_foe_mc->GetParameter(0);

    TLegend *l_foe_mc = new TLegend(0.2, 0.65, 0.35, 0.85);
    l_foe_mc->SetTextSize(0.04);
    l_foe_mc->SetFillColor(kWhite);
    l_foe_mc->SetLineColor(kWhite);
    // l_phi_foe_mc->AddEntry(fr_PhiMass_foe_mc->findObject("ldh_PhiMass_foe_mc"), "Data", "p");
    l_foe_mc->AddEntry("fsb_foe_mc", "total", "l");
    l_foe_mc->AddEntry("fs_foe_mc", "Voigtian", "l");
    l_foe_mc->AddEntry("fb_foe_mc", "pol 4^{th}", "l");
    l_foe_mc->Draw();

    TLatex lat_foe_mc;
    lat_foe_mc.SetTextSize(0.04);
    lat_foe_mc.SetTextAlign(13); //align at top
    lat_foe_mc.SetNDC();
    lat_foe_mc.SetTextColor(kBlue);
    lat_foe_mc.DrawLatex(0.6, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_foe_mc->GetChisquare() / fsb_foe_mc->GetNDF()));
    lat_foe_mc.DrawLatex(0.6, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N_foe_mc, dN_foe_mc));
    lat_foe_mc.DrawLatex(0.6, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_foe_mc->GetParameter(1), fsb_foe_mc->GetParError(1)));
    lat_foe_mc.DrawLatex(0.6, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb_foe_mc->GetParameter(2), fsb_foe_mc->GetParError(2)));
    lat_foe_mc.DrawLatex(0.6, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_foe_mc->GetParameter(3), fsb_foe_mc->GetParError(3)));

    // ========================= fo vs. Phi(1020) ===========================

    cphifo->cd(i);
    TH2F *h2d2 = new TH2F("h2d2", "(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, n2k, mkk_min, mkk_max);
    tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[i], Eg2[i]));
    h2d2->Draw("colz");

    ofs_phifo << "########  i = " << i << " | Eg1[" << i << "] = " << Eg1[i] << " | Eg2[" << i << "] = " << Eg2[i] << endl;

    gphifo = new TGraphErrors(); //n2pi
    gphifo->SetMarkerStyle(20);
    gphifo->SetMarkerSize(1.0);
    gphifo->SetMarkerColor(1);
    gphifo->SetMinimum(0.);
    gphifo->SetTitle(Form("E_{#gamma} = %.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", h2foe_mc->GetXaxis()->GetBinCenter(i)));

    // gnophifo->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{bkg}", Eg1[j], Eg2[j]));


    for (int j = 1; j <= n2pi; ++j) //n2pi
    {
      // cout << j << " " << flush;

      if (j < 26)
        cphifo1->cd(j); // if (j < 26) cphifo1->cd(j);
      if (j > 25)
        cphifo2->cd(j - 25); // if (j > 25 && j < 51) cphifo2->cd(j - 25);

      TH1D *hphifo_py = h2d2->ProjectionY(Form("_hphifo_py_%d", j), j, j);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifo_py->Draw();
      // fb2->Draw("same");

      TF1 *fsb_phifo_data = new TF1("fsb_phifo_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      // TF1 *fsb_phifo_data = new TF1("fsb_phifo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      fsb_phifo_data->SetLineColor(2);
      fsb_phifo_data->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb_phifo_data->SetParLimits(0, 0, 10000);
      fsb_phifo_data->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
      fsb_phifo_data->FixParameter(2, 0.004);        //fsb_phifo_data->SetParLimits(2, 0.008, 0.010);
      fsb_phifo_data->FixParameter(3, 0.002);        //fsb_phifo_data->SetParLimits(3, 0.001,0.01);// 0.001,0.01

      TF1 *fs_phifo_data = new TF1("fs_phifo_data", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      fs_phifo_data->SetLineColor(4);

      TF1 *fb_phifo_data = new TF1("fb_phifo_data", "pol2(4)", mkk_min, mkk_max); //pol2(3)
      fb_phifo_data->SetLineColor(28);
      fb_phifo_data->SetLineStyle(2);

      hphifo_py->Fit("fsb_phifo_data", "", "", mkk_min, mkk_max);
      double par[7]; //6
      fsb_phifo_data->GetParameters(&par[0]);
      fs_phifo_data->SetParameters(&par[0]);
      fb_phifo_data->SetParameters(&par[4]); //3

      fs_phifo_data->Draw("same");
      fb_phifo_data->Draw("same");

      double N2 = fs_phifo_data->Integral(mkk_min, mkk_max) / hphifo_py->GetBinWidth(1);
      double dN2 = N2 * fsb_phifo_data->GetParError(0) / fsb_phifo_data->GetParameter(0);

      gphifo->SetPoint(j - 1, h2d2->GetXaxis()->GetBinCenter(j), N2);
      gphifo->SetPointError(j - 1, 0, dN2);

      // gphifo_width[j]->SetPoint(j - 1, h2d2->GetXaxis()->GetBinCenter(j), fsb_phifo_data->GetParameter(2));
      // gphifo_width[j]->SetPointError(j - 1, 0, fsb_phifo_data->GetParError(2));

      // gphifo_mean[j]->SetPoint(j - 1, h2d2->GetXaxis()->GetBinCenter(j), fsb_phifo_data->GetParameter(3));
      // gphifo_mean[j]->SetPointError(j - 1, 0, fsb_phifo_data->GetParError(3));

      TLatex lat_phifo;
      lat_phifo.SetTextSize(0.09);
      lat_phifo.SetTextAlign(13); //align at top
      lat_phifo.SetNDC();
      lat_phifo.SetTextColor(kBlue);
      lat_phifo.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb_phifo_data->GetChisquare() / fsb_phifo_data->GetNDF()));
      lat_phifo.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
      lat_phifo.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_phifo_data->GetParameter(1), fsb_phifo_data->GetParError(1)));
      lat_phifo.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb_phifo_data->GetParameter(2), fsb_phifo_data->GetParError(2)));
      lat_phifo.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_phifo_data->GetParameter(3), fsb_phifo_data->GetParError(3)));

      hphifo_py->Write();
      cgphifo->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphifo->cd(i);
    gphifo->Draw("AP");
    TF1 *fsb_foe_data = new TF1("fsb_foe_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", 0.85, 1.15);
    // TF1 *fsb_foe_data = new TF1("fsb_foe_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)", 0.85, 1.05);
    fsb_foe_data->SetLineColor(2);
    fsb_foe_data->SetParameters(1433, 0.980, 0.01, 0.06, 1,1,1,1,1); //0.01, 
    // fsb_foe_data->SetParLimits(0, 0, 10000);
    fsb_foe_data->FixParameter(1, fsb_foe_mc->GetParameter(1));//SetParLimits(1, 0.95, 0.98); // 1.018, 1.021
    fsb_foe_data->FixParameter(2, fsb_foe_mc->GetParameter(2)); //SetParLimits(2, 0.001, 0.6); //fsb_foe_data->SetParLimits(2, 0.008, 0.010);
    fsb_foe_data->FixParameter(3, fsb_foe_mc->GetParameter(3)); //SetParLimits(2, 0.001, 0.6); //fsb_foe_data->SetParLimits(2, 0.008, 0.010);
    // fsb_foe_data->FixParameter(3, fsb_foe_mc->GetParameter(3));        //fsb_foe_data->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_foe_data = new TF1("fs_foe_data", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.85, 1.15);
    // TF1 *fs_foe_data = new TF1("fs_foe_data", "[0]*TMath::BreitWigner(x,[1],[2])", 0.85, 1.05);
    fs_foe_data->SetLineColor(4);

    TF1 *fb_foe_data = new TF1("fb_foe_data", "pol4(4)", 0.85, 1.15); //pol2(3)
    fb_foe_data->SetLineColor(28);
    fb_foe_data->SetLineStyle(2);

    gphifo->Fit("fsb_foe_data", "", "", 0.85, 1.15);
    double par[10]; //6
    fsb_foe_data->GetParameters(&par[0]);
    fs_foe_data->SetParameters(&par[0]);
    fb_foe_data->SetParameters(&par[4]); //3

    double N_foe_data = fs_foe_data->Integral(0.85, 1.15) * n2pi / (m2pi_max - m2pi_min);
    double dN_foe_data = N_foe_data * fsb_foe_data->GetParError(0) / fsb_foe_data->GetParameter(0);

    fs_foe_data->Draw("same");
    fb_foe_data->Draw("same");

    TLatex lat_phifo;
    lat_phifo.SetTextSize(0.05);
    lat_phifo.SetTextAlign(13); //align at top
    lat_phifo.SetNDC();
    lat_phifo.SetTextColor(kBlue);
    lat_phifo.DrawLatex(0.6, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_foe_data->GetChisquare() / fsb_foe_data->GetNDF()));
    lat_phifo.DrawLatex(0.6, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N_foe_data, dN_foe_data));
    lat_phifo.DrawLatex(0.6, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb_foe_data->GetParameter(1), fsb_foe_data->GetParError(1)));
    lat_phifo.DrawLatex(0.6, 0.58, Form("#Gamma = %0.3f#pm%0.3f", fsb_foe_data->GetParameter(2), fsb_foe_data->GetParError(2)));

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phifoe = N_foe_mc / h_beame_tru->GetBinContent(i); // Efficiency = N_observed/N_generated
    double deff_phifoe = eff_phifoe * (dN_foe_data / N_foe_mc);
    gphifoeeff->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phifoe * 100);
    gphifoeeff->SetPointError(i - 1, 0, deff_phifoe * 100); //->GetBinError(i)

    // +++++++++++++++++++ Systematic Errors
    // double width_value[3] = {-169.74, -141.06, -197.77}; // Y Width value : 0.079, 0.065, 0.093
    // double err_width_value = TMath::RMS(3, width_value);

    double N_bkg_poln_1[10] = {1689, 1753, 2465, 6581, 8566, 2607, 3110, 2557, 2450, 160}; // Bkg pol degree: 4
    double N_bkg_poln_2[10] = {1685, 1749, 2459, 6567, 8548, 2677, 3103, 2552, 2444, 182}; // Bkg pol degree: 3
    double N_bkg_poln_3[10] = {1694, 1756, 2471, 6595, 8583, 2611, 3118, 2562, 2456, 161}; // Bkg pol degree: 5
    double bkg_poln[3] = {N_bkg_poln_1[i - 1], N_bkg_poln_2[i - 1], N_bkg_poln_3[i - 1]};
    double err_bkg_poln = TMath::RMS(3, bkg_poln);
    double fit_range_1[10] = {1689, 1753, 2465, 6581, 8566, 2607, 3110, 2557, 2450, 160}; // Fit range: [0.98, 1.2]
    double fit_range_2[10] = {1685, 1767, 2477, 6598, 8551, 2600, 3094, 2532, 2439, 182}; // Fit range: [0.98, 1.18]
    double fit_range_3[10] = {1693, 1736, 2458, 6609, 8522, 2582, 3104, 2540, 2433, 163}; // Fit range: [0.98, 1.22]
    double fit_range[3] = {fit_range_1[i - 1], fit_range_2[i - 1], fit_range_3[i - 1]};
    double err_fit_range = TMath::RMS(3, fit_range);
    // double fit_func_1[10] = {1689.35, 1752.50, 2465.09, 6580.90, 8565.59, 2606.73, 3110.43, 2557.01, 2449.74, 160.43}; // Fit function: Voigt
    // double fit_func_2[10] = {2157.41, 2302.83, 3158.83, 7835.29, 10889.65, 3077.05, 3893.53, 3240.87, 3128.20, 217.33}; // Fit function: Breit-wigner
    // double fit_func_3[10] = {1608.70, 1734.22, 2367.27, 6358.65, 8255.21, 2452.65, 2940.46, 2476.66, 2322.49, 203.03}; // Fit function: Double Gaus
    // double fit_func[3] = {fit_func_1[i-1], fit_func_2[i-1], fit_func_3[i-1]};
    // double err_fit_func = TMath::RMS(3, fit_func);

    double err_sys = TMath::Sqrt(err_bkg_poln * err_bkg_poln + err_fit_range * err_fit_range);

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phifoe = h_tagged_flux->GetBinContent(i); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    if (lumi_phifoe <= 0)
    {
      //gphiexsec->RemovePoint(i - 1);
      continue;
    }
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_phifoe = 1e9 * N_foe_data / (eff_phifoe * lumi_phifoe * T * br_phi);
    double dxsec_phifoe_stat = xsec_phifoe * (dN_foe_data / N_foe_data);
    double dxsec_phifoe_sys = xsec_phifoe * TMath::Sqrt((err_sys / N_foe_data) * (err_sys / N_foe_data) + (deff_phifoe / eff_phifoe) * (deff_phifoe / eff_phifoe) + (dbr_phi / br_phi) * (dbr_phi / br_phi));
    // double xsec_phi = N_foe / lumi_phi;
    // double dxsec_phi = dN_foe / lumi_phi;
    gphifoexsec->SetPoint(gphifoexsec->GetN(), h2phie->GetXaxis()->GetBinCenter(i), xsec_phifoe);
    gphifoexsec->SetPointError(gphifoexsec->GetN() - 1, 0, dxsec_phifoe_stat);

    fprintf(table_xsec_phifo, "%0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h2phie_mc->GetXaxis()->GetBinCenter(i), h_beame_tru->GetBinContent(i), N_foe_mc, dN_foe_mc, N_foe_data, dN_foe_data, eff_phifoe * 100, deff_phifoe * 100, xsec_phifoe, dxsec_phifoe_stat, dxsec_phifoe_sys);

    fprintf(table_xsec_phifo_sys, "%0.2f & %0.2f & %0.2f \\\\ \n", h2phie_mc->GetXaxis()->GetBinCenter(i), err_bkg_poln * 100 / N_bkg_poln_1[i - 1], err_fit_range * 100 / fit_range_1[i - 1]);
    // table_phi << std::setprecision(2) << std::fixed;
    // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_foe_mc<<"$\\pm$"<<dN_foe_mc<<"&"<<N_foe_data<<"$\\pm$"<<dN_foe_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phifo<<"$\\pm$"<<dxsec_phifo_stat<<"$\\pm$"<<dxsec_phifo_sys<<" \\\\"<< endl;
  }

  cgphie->cd();
  gphie->Draw("AP");
  gphie->Write(Form("h%s_gphie", name.Data()), TObject::kWriteDelete);
  cgphie_mc->cd();
  gphie_mc->Draw("AP");
  cgphifoeeff->cd();
  gphifoeeff->Draw("AP");
  cgphifoexsec->cd();
  gphifoexsec->Draw("AP");
  gphifoexsec->Write(Form("h%s_gphifoexsec", name.Data()), TObject::kWriteDelete);
  // int j =1;
  // gphie->Write(Form("grphie_%d", j), TObject::kWriteDelete);

  // // total Cross-sections
  // TMultiGraph *mgxsec = new TMultiGraph();
  // TCanvas *cmgxsec = new TCanvas("cmgxsec", "cmgxsec", 900, 600);
  // TLegend *lmgxsec = new TLegend(0.86, 0.86, 0.98, 0.98);
  // lmgxsec->SetTextSize(0.04);
  // lmgxsec->SetBorderSize(0);
  // cmgxsec->cd();
  // cmgxsec->SetGrid();
  // TGraphErrors *hdata_16_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_16_gphiexsec");
  // cout << " ***** hdata_16_gphiexsec = " << hdata_16_gphiexsec << endl;
  // hdata_16_gphiexsec->SetMarkerColor(1);
  // hdata_16_gphiexsec->SetMarkerSize(1.5);
  // hdata_16_gphiexsec->SetLineColor(1);
  // hdata_16_gphiexsec->SetMarkerStyle(20);
  // mgxsec->Add(hdata_16_gphiexsec);
  // TGraphErrors *hdata_17_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_17_gphiexsec");
  // cout << " ***** hdata_17_gphiexsec = " << hdata_17_gphiexsec << endl;
  // hdata_17_gphiexsec->SetMarkerColor(kBlue);
  // hdata_17_gphiexsec->SetMarkerSize(1.5);
  // hdata_17_gphiexsec->SetLineColor(kBlue);
  // hdata_17_gphiexsec->SetMarkerStyle(20);
  // mgxsec->Add(hdata_17_gphiexsec);
  // TGraphErrors *hdata_18_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_18_gphiexsec");
  // cout << " ***** hdata_18_gphiexsec = " << hdata_18_gphiexsec << endl;
  // hdata_18_gphiexsec->SetMarkerColor(kRed);
  // hdata_18_gphiexsec->SetMarkerSize(1.5);
  // hdata_18_gphiexsec->SetLineColor(kRed);
  // hdata_18_gphiexsec->SetMarkerStyle(20);
  // mgxsec->Add(hdata_18_gphiexsec);
  // mgxsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. Beam Energy; E_{#gamma} [GeV]; #sigma [nb]");
  // mgxsec->Draw("AP");
  // mgxsec->SetMinimum(0.);
  // cmgxsec->Modified();
  // lmgxsec->AddEntry(hdata_16_gphiexsec, "2016", "p");
  // lmgxsec->AddEntry(hdata_17_gphiexsec, "2017", "p");
  // lmgxsec->AddEntry(hdata_18_gphiexsec, "2018", "p");
  // lmgxsec->Draw();
  // cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.eps", "eps");
  // cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.png", "png");
  // cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.root", "root");

  // // total Yields
  // TMultiGraph *mgphie = new TMultiGraph();
  // TCanvas *cmgphie = new TCanvas("cmgphie", "cmgphie", 900, 600);
  // TLegend *lmgphie = new TLegend(0.86, 0.86, 0.98, 0.98);
  // lmgphie->SetTextSize(0.04);
  // lmgphie->SetBorderSize(0);
  // cmgphie->cd();
  // cmgphie->SetGrid();
  // TGraphErrors *hdata_16_gphie = (TGraphErrors *)outputfig->Get("hdata_16_gphie");
  // cout << " ***** hdata_16_gphie = " << hdata_16_gphie << endl;
  // hdata_16_gphie->SetMarkerColor(1);
  // hdata_16_gphie->SetMarkerSize(1.5);
  // hdata_16_gphie->SetLineColor(1);
  // hdata_16_gphie->SetMarkerStyle(20);
  // mgphie->Add(hdata_16_gphie);
  // TGraphErrors *hdata_17_gphie = (TGraphErrors *)outputfig->Get("hdata_17_gphie");
  // cout << " ***** hdata_17_gphie = " << hdata_17_gphie << endl;
  // hdata_17_gphie->SetMarkerColor(kBlue);
  // hdata_17_gphie->SetMarkerSize(1.5);
  // hdata_17_gphie->SetLineColor(kBlue);
  // hdata_17_gphie->SetMarkerStyle(20);
  // mgphie->Add(hdata_17_gphie);
  // TGraphErrors *hdata_18_gphie = (TGraphErrors *)outputfig->Get("hdata_18_gphie");
  // cout << " ***** hdata_18_gphie = " << hdata_18_gphie << endl;
  // hdata_18_gphie->SetMarkerColor(kRed);
  // hdata_18_gphie->SetMarkerSize(1.5);
  // hdata_18_gphie->SetLineColor(kRed);
  // hdata_18_gphie->SetMarkerStyle(20);
  // mgphie->Add(hdata_18_gphie);
  // mgphie->SetTitle("#phi(1020) Yield vs. Beam energy (data); E_{#gamma} [GeV]; N_{#phi}");
  // mgphie->Draw("AP");
  // mgphie->SetMinimum(0.);
  // cmgphie->Modified();
  // lmgphie->AddEntry(hdata_16_gphie, "2016", "p");
  // lmgphie->AddEntry(hdata_17_gphie, "2017", "p");
  // lmgphie->AddEntry(hdata_18_gphie, "2018", "p");
  // lmgphie->Draw();
  // cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphie.eps", "eps");
  // cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphie.png", "png");
  // cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphie.root", "root");

  cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie.root", name.Data()), "root");
  cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie.eps", name.Data()), "eps");
  cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie.png", name.Data()), "png");
  cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie1.root", name.Data()), "root");
  cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie1.eps", name.Data()), "eps");
  cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phie1.png", name.Data()), "png");
  cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphie.root", name.Data()), "root");
  cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphie.eps", name.Data()), "eps");
  cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphie.png", name.Data()), "png");
  cphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie_mc.root", "root");
  cphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie_mc.eps", "eps");
  cphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie_mc.png", "png");
  cphie1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie1_mc.root", "root");
  cphie1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie1_mc.eps", "eps");
  cphie1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphie1_mc.png", "png");
  cgphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphie_mc.root", "root");
  cgphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphie_mc.eps", "eps");
  cgphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphie_mc.png", "png");
  cgphifoeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphifoeeff.root", "root");
  cgphifoeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphifoeeff.eps", "eps");
  cgphifoeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphifoeeff.png", "png");
  cgphifoexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphifoexsec.root", name.Data()), "root");
  cgphifoexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphifoexsec.eps", name.Data()), "eps");
  cgphifoexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphifoexsec.png", name.Data()), "png");
*/
   /*
  // ======================================== Phi vs. -t ===============================================

  //root -l 'xsec_phifo.C+("data_17",100,1,10)'

  FILE *table_xsec_phifo = fopen("table_xsec_phifo.tex","w");
  fprintf(table_xsec_phifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Total cross-sections in momentum transfer bins}\n \\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}\n \\hline\n -t $[GeV^{2}]$ & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma$ [nb] \\\\ \n \\hline\n");

  FILE *table_xsec_phifo_sys = fopen("table_xsec_phifo_sys.tex","w");
  fprintf(table_xsec_phifo_sys,"\\documentclass[12pt]{extarticle}\n \\usepackage{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\begin{minipage}{\\textwidth}\n \\begin{center}\n \\caption{Systematic errors in momentum transfer bins}\n \\begin{tabularx}{\\textwidth}{|X|X|X|}\n \\hline\n -t $[GeV^{2}]$ & polynomial degrees [ \\% ] & Fit range [ \\% ]  \\\\ \n \\hline\n");  
  
  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux" << h_tagged_flux << endl;
  h_tagged_flux->Rebin(100);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", 1, 6, 12);
  ttru->Project("h_beame_tru", "ThrownBeam__P4.E()", "ThrownBeam__P4.E()>6");
  cout << "h_beame_tru" << h_beame_tru << endl;
  // h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");

  // root -l 'xsec_phifo.C+(100,20,20)'
  // Data
  TCanvas *cphit = new TCanvas("cphit", "cphit", 900, 600);
  cphit->cd();
  TH2D *h2phit = new TH2D("h2phit", "Data; -t [GeV^{2}]; m_{K^{+}K^{-}} [GeV/c^{2}]", nt, 0, 4, n2k, mkk_min, mkk_max);
  tdata->Project("h2phit", "kpkm_mf:-t_kin", "w8*(kpkm_uni && beam_p4_kin.E()>6)");

  h2phit->Draw("colz");

  TCanvas *cphit1 = new TCanvas("cphit1", "cphit1", 1500, 600);
  cphit1->Divide(5, 2);
  TCanvas *cgphit = new TCanvas("cgphit", "cgphit", 900, 600);
  TGraphErrors *gphit = new TGraphErrors();
  gphit->SetMarkerStyle(20);
  gphit->SetMarkerSize(1.5);
  gphit->SetMarkerColor(1);
  gphit->SetMinimum(0.0);
  gphit->SetTitle("#phi(1020) Yield vs. Momentum transfer (data); -t [GeV^{2}]; N_{#phi}");

  // Monte Carlo
  TCanvas *cphit_mc = new TCanvas("cphit_mc", "cphit_mc", 900, 600);
  cphit_mc->cd();
  TH2D *h2phit_mc = new TH2D("h2phit_mc", "MC; -t [GeV^{2}]; m_{K^{+}K^{-}} [GeV/c^{2}]", nt, 0, 4, n2k, mkk_min, mkk_max);//0.98, 1.2
  tmc->Project("h2phit_mc", "kpkm_mf:-t_kin", "w8*(kpkm_uni && beam_p4_kin.E()>6)");
  h2phit_mc->Draw("colz");

  TCanvas *cphit1_mc = new TCanvas("cphit1_mc", "cphit1_mc", 1500, 600);
  cphit1_mc->Divide(5, 2);
  TCanvas *cgphit_mc = new TCanvas("cgphit_mc", "cgphit_mc", 900, 600);
  TGraphErrors *gphit_mc = new TGraphErrors();
  gphit_mc->SetMarkerStyle(20);
  gphit_mc->SetMarkerSize(1.5);
  gphit_mc->SetMarkerColor(1);
  gphit_mc->SetMinimum(0.0);
  gphit_mc->SetTitle("#phi(1020) Yield vs. Momentum transfer (MC); -t [GeV^{2}]; N_{#phi}");

  // Efficiency
  TCanvas *cgphiteff = new TCanvas("cgphiteff", "cgphiteff", 900, 600);
  TGraphErrors *gphiteff = new TGraphErrors();
  gphiteff->SetMarkerStyle(20);
  gphiteff->SetMarkerSize(1.5);
  gphiteff->SetMarkerColor(kBlue);
  gphiteff->SetMinimum(0.0);
  gphiteff->SetTitle("#phi(1020) Effeciency vs. Momentum transfer; -t [GeV^{2}]; #epsilon [%]");

  // Cross-section
  TCanvas *cgphitxsec = new TCanvas("cgphitxsec", "cgphitxsec", 900, 600);
  TGraphErrors *gphitxsec = new TGraphErrors();
  gphitxsec->SetMarkerStyle(20);
  gphitxsec->SetMarkerSize(1.5);
  gphitxsec->SetMarkerColor(1);
  gphitxsec->SetMinimum(0.0);
  gphitxsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. Momentum transfer; -t [GeV^{2}]; #sigma [nb]"); //#phi(1020) flux normalized yield, Yield_{#phi}

  for (int i = 1; i <= nt; ++i)
  {
    cout << i << " " << flush;

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphit1_mc->cd(i);
    TH1D *hphit_mc_py = h2phit_mc->ProjectionY(Form("hphit_mc_py_%d", i), i, i);
    hphit_mc_py->SetTitle(Form("-t = %.2f [GeV^{2}] (MC)",h2phit_mc->GetXaxis()->GetBinCenter(i)));
    hphit_mc_py->Draw("e");

    // w.factory(Form("Voigtian::sig_phit_mc(m_phit_mc[%f,%f],mean_phit_mc[1.015,1.022],width_phit_mc[0.004],sigma_phit_mc[0.002])", mkk_min, mkk_max)); //sigma_phit_mc[0.0001,0.01], mean_phit_mc[1.011,1.030]
    // // w.factory(Form("BreitWigner::sig_phit_mc(m_phit_mc[%f,%f],mean_phit_mc[1.018,1.021],width_phit_mc[0.004])", mkk_min, mkk_max));
    // w.factory("Chebychev::bkg_phit_mc(m_phit_mc,{c0_phit_mc[-10,10], c1_phit_mc[-10,10], c2_phit_mc[-10,10]})");
    // w.factory("SUM:::model_phit_mc(nsig_phit_mc[0,100000000]*sig_phit_mc, nbkg_phit_mc[0,100000000]*bkg_phit_mc)"); //nsig[0,100000000]*sig2,
    // w.var("m_phit_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    // RooDataHist dh_phit_mc("dh_phit_mc", "dh_phit_mc", *w.var("m_phit_mc"), Import(*hphit_mc_py));
    // RooPlot *fr_phit_mc = w.var("m_phit_mc")->frame(Title(Form("-t = %f",h2phit_mc->GetXaxis()->GetBinCenter(i))));
    // w.pdf("model_phit_mc")->fitTo(dh_phit_mc);

    // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    // dh_phit_mc.plotOn(fr_phit_mc, RooFit::Name("ldh_phit_mc"));
    // w.pdf("model_phit_mc")->plotOn(fr_phit_mc, Components(*w.pdf("sig_phit_mc")), LineColor(kRed), RooFit::Name("lsig_phit_mc"));
    // w.pdf("model_phit_mc")->plotOn(fr_phit_mc, Components(*w.pdf("bkg_phit_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phit_mc"));
    // w.pdf("model_phit_mc")->plotOn(fr_phit_mc, RooFit::Name("lmodel_phit_mc"));
    // // w.pdf("model_phit_mc")->paramOn(fr_phit_mc, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phit_mc"), *w.var("nbkg_phit_mc")))); //,*w.var("mean_phit_mc"),*w.var("width_phit_mc"),*w.var("sigma_phit_mc"))));
    // fr_phit_mc->Draw();

    // TLegend *l_phit_mc = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_phit_mc->SetFillColor(kWhite);
    // l_phit_mc->AddEntry(fr_phit_mc->findObject("lsig_phit_mc"), Form("N_{Sig} = %.2f", w.var("nsig_phit_mc")->getVal()), "l");
    // l_phit_mc->AddEntry(fr_phit_mc->findObject("lbkg_phit_mc"), Form("N_{Bkg} = %.2f", w.var("nbkg_phit_mc")->getVal()), "l");
    // l_phit_mc->Draw();

    // double N_phit_mc = w.var("nsig_phit_mc")->getVal();
    // double dN_phit_mc = w.var("nsig_phit_mc")->getError();

    // TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol4(3)", mkk_min, mkk_max);
    // fsb_mc->SetParameters(1433, 1.019, 0.005, 1,1,1,1,1); // breitwigner
    // TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    // TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5]) + pol4(6)", mkk_min, mkk_max);
    // fsb_mc->SetParameters(1433, 1.019, 0.003, 1433, 1.023, 0.005, 1,1,1,1,1); // double gaus
    // TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])", mkk_min, mkk_max);
   
    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1); // voigt
    fsb_mc->SetParLimits(1, 1.015, 1.022); //fsb_mc->FixParameter(1, 1.019);
    // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_mc->FixParameter(3, 0.0042);   //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max);
    fb_mc->SetLineColor(28);
    fb_mc->SetLineStyle(2);

    hphit_mc_py->Fit("fsb_mc", "", "", mkk_min, mkk_max);
    double par_mc[10]; //6
    fsb_mc->GetParameters(&par_mc[0]);
    fs_mc->SetParameters(&par_mc[0]);
    fb_mc->SetParameters(&par_mc[4]); //4

    fs_mc->Draw("same");
    fb_mc->Draw("same");

    double N_phit_mc = fs_mc->Integral(mkk_min, mkk_max) / hphit_mc_py->GetBinWidth(1);
    double dN_phit_mc = N_phit_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

    TLatex lat_phit_mc;
    lat_phit_mc.SetTextSize(0.06);
    lat_phit_mc.SetTextAlign(13); //align at top
    lat_phit_mc.SetNDC();
    lat_phit_mc.SetTextColor(kBlue);
    lat_phit_mc.DrawLatex(0.45, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
    lat_phit_mc.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_mc, dN_phit_mc));
    lat_phit_mc.DrawLatex(0.45, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
    lat_phit_mc.DrawLatex(0.45, 0.66, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phit_mc.DrawLatex(0.45, 0.59, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

    gphit_mc->SetPoint(i - 1, h2phit_mc->GetXaxis()->GetBinCenter(i), N_phit_mc);
    gphit_mc->SetPointError(i - 1, 0, dN_phit_mc);

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphit1->cd(i);
    TH1D *hphit_py = h2phit->ProjectionY(Form("hphit_py_%d", i), i, i);
    hphit_py->SetTitle(Form("-t = %.2f [GeV^{2}] (Data)",h2phit->GetXaxis()->GetBinCenter(i)));
    hphit_py->Draw("e");

    // w.factory(Form("Voigtian::sig_phit(m_phit[%f,%f],mean_phit[1.015, 1.022],width_phit[0.004],sigma_phit[0.002])", mkk_min, mkk_max)); //sigma_phit[0.0001,0.01], mean_phit[1.016,1.022 or 1.010,1.030]
    // // w.factory(Form("BreitWigner::sig_phit(m_phit[%f,%f],mean_phit[1.018,1.021],width_phit[0.004])", mkk_min, mkk_max));
    // w.factory("Chebychev::bkg_phit(m_phit,{c0_phit[-10,10], c1_phit[-10,10], c2_phit[-10,10]})");
    // w.factory("SUM:::model_phit(nsig_phit[0,100000000]*sig_phit, nbkg_phit[0,100000000]*bkg_phit)"); //nsig[0,100000000]*sig2,
    // w.var("m_phit")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    // RooDataHist dh_phit("dh_phit", "dh_phit", *w.var("m_phit"), Import(*hphit_py));
    // RooPlot *fr_phit = w.var("m_phit")->frame(Title(Form("-t = %f",h2phit->GetXaxis()->GetBinCenter(i))));
    // w.pdf("model_phit")->fitTo(dh_phit);
    // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    // dh_phit.plotOn(fr_phit, RooFit::Name("ldh_phit"));
    // w.pdf("model_phit")->plotOn(fr_phit, Components(*w.pdf("sig_phit")), LineColor(kRed), RooFit::Name("lsig_phit"));
    // w.pdf("model_phit")->plotOn(fr_phit, Components(*w.pdf("bkg_phit")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phit"));
    // w.pdf("model_phit")->plotOn(fr_phit, RooFit::Name("lmodel_phit"));
    // // w.pdf("model_phit")->paramOn(fr_phit, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phit"), *w.var("nbkg_phit")))); //,*w.var("mean_phit"),*w.var("width_phit"),*w.var("sigma_phit"))));
    // fr_phit->Draw();

    // TLegend *l_phit = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_phit->SetFillColor(kWhite);
    // l_phit->SetLineColor(kWhite);
    // // l_phit->AddEntry(fr_phit->findObject("ldh_phit"), "Data", "p");
    // // l_phit->AddEntry(fr_phit->findObject("lmodel_phit"), "total", "l");
    // l_phit->AddEntry(fr_phit->findObject("lsig_phit"), Form("N_{Sig} = %.2f", w.var("nsig_phit")->getVal()), "l");
    // l_phit->AddEntry(fr_phit->findObject("lbkg_phit"), Form("N_{Bkg} = %.2f", w.var("nbkg_phit")->getVal()), "l");
    // l_phit->Draw();

    // double N_phit = w.var("nsig_phit")->getVal();
    // double dN_phit = w.var("nsig_phit")->getError();

    // //+++++++++ BreitWigner + pol4
    // TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol4(3)", mkk_min, mkk_max);
    // fsb_data->SetParameters(1433, 1.019, 0.005, 1,1,1,1,1);  // breiwigner
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1));
    // // fsb_data->FixParameter(2, fsb_mc->GetParameter(3));
    
    // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    // fs_data->SetLineColor(4);

    // TF1 *fb_data = new TF1("fb_data", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    // fb_data->SetLineColor(28);
    // fb_data->SetLineStyle(2);

    // hphit_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    // double par_data[10]; //6
    // fsb_data->GetParameters(&par_data[0]);
    // fs_data->SetParameters(&par_data[0]);
    // fb_data->SetParameters(&par_data[3]);

    // //+++++++++ Double Gaus + pol4
    // TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5]) + pol4(6)", mkk_min, mkk_max);
    // fsb_data->SetParameters(1433, 1.019, 0.003, 1433, 1.023, 0.005, 1,1,1,1,1);  // double gaus
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1));
    // // fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); // double gaus
    // fsb_data->SetParLimits(4, 1.19, 1.023);
    // // fsb_data->FixParameter(5, fsb_mc->GetParameter(3)); // double gaus
    // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[4],[5])", mkk_min, mkk_max);
    // fs_data->SetLineColor(4);

    // TF1 *fb_data = new TF1("fb_data", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    // fb_data->SetLineColor(28);
    // fb_data->SetLineStyle(2);

    // hphit_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    // double par_data[10]; //6
    // fsb_data->GetParameters(&par_data[0]);
    // fs_data->SetParameters(&par_data[0]);
    // fb_data->SetParameters(&par_data[6]);

    // //+++++++++ Voigt + |x-a|^b*|x-c|^d
    // TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + [4]*(x>[5])*(x<[6])*abs(x-[5])^[7]*abs([6]-x)^[8]", mkk_min, mkk_max);
    // fsb_data->SetLineColor(2);
    // fsb_data->SetParameters(1433, 1.019, 0.003, 0.0042, 1,0.98,1.2,0.5,4);  // voigt 1,0.98,1.2,0.5,4
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); //fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    // fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); //fsb->SetParLimits(2, 0.008, 0.010);
    // fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01
    // // fsb_data->SetParLimits(5, 0.97, 0.99);
    // fsb_data->FixParameter(4, 19);
    // fsb_data->FixParameter(5, 0.98);
    // fsb_data->FixParameter(6, 2.82); //2.82
    // fsb_data->FixParameter(7, 0.67);
    // fsb_data->FixParameter(8, 6.82);

    // TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    // fs_data->SetLineColor(4);

    // TF1 *fb_data = new TF1("fb_data", "[0]*(x>[1])*(x<[2])*abs(x-[1])^[3]*abs([2]-x)^[4]", mkk_min, mkk_max); //
    // fb_data->SetLineColor(28);
    // fb_data->SetLineStyle(2);

    // hphit_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    // double par_data[10]; //6
    // fsb_data->GetParameters(&par_data[0]);
    // fs_data->SetParameters(&par_data[0]);
    // fb_data->SetParameters(&par_data[4]); //4

    // +++++++++ Voigt + pol4
    TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    fsb_data->SetLineColor(2);
    fsb_data->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1);  // voigt
    fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); //fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_data->SetLineColor(4);

    TF1 *fb_data = new TF1("fb_data", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    fb_data->SetLineColor(28);
    fb_data->SetLineStyle(2);

    hphit_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    double par_data[10]; //6
    fsb_data->GetParameters(&par_data[0]);
    fs_data->SetParameters(&par_data[0]);
    fb_data->SetParameters(&par_data[4]); //4

    fs_data->Draw("same");
    fb_data->Draw("same");

    double N_phit_data = fs_data->Integral(mkk_min, mkk_max) / hphit_py->GetBinWidth(1);
    double dN_phit_data = N_phit_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    TLatex lat_phit_data;
    lat_phit_data.SetTextSize(0.06);
    lat_phit_data.SetTextAlign(13); //align at top
    lat_phit_data.SetNDC();
    lat_phit_data.SetTextColor(kBlue);
    lat_phit_data.DrawLatex(0.45, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phit_data.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_data, dN_phit_data));
    lat_phit_data.DrawLatex(0.45, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phit_data.DrawLatex(0.45, 0.66, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phit_data.DrawLatex(0.45, 0.59, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

    gphit->SetPoint(i - 1, h2phit->GetXaxis()->GetBinCenter(i), N_phit_data);
    gphit->SetPointError(i - 1, 0, dN_phit_data);    

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phi = N_phit_mc / h_beame_tru->GetBinContent(1); // Efficiency = N_observed/N_generated
    double deff_phi = eff_phi * (dN_phit_mc / N_phit_mc);
    gphiteff->SetPoint(i - 1, h2phit_mc->GetXaxis()->GetBinCenter(i), eff_phi*100);
    gphiteff->SetPointError(i - 1, 0, deff_phi*100); //->GetBinError(i)

    // +++++++++++++++++++ Systematic Errors 
    // double width_value[3] = {-169.74, -141.06, -197.77}; // Y Width value : 0.079, 0.065, 0.093  
    // double err_width_value = TMath::RMS(3, width_value);

    double N_bkg_poln_1[10] = {4274, 8971, 6453, 4383, 2550, 1544, 1056, 652, 441, 319}; // Bkg pol degree: 4
    double N_bkg_poln_2[10] = {4259, 8953, 6438, 4376, 2546, 1586, 1093, 673, 440, 338}; // Bkg pol degree: 3
    double N_bkg_poln_3[10] = {4288, 8989, 6468, 4390, 2555, 1547, 1057, 653, 442, 320}; // Bkg pol degree: 5
    double bkg_poln[3] = {N_bkg_poln_1[i-1], N_bkg_poln_2[i-1], N_bkg_poln_3[i-1]};
    double err_bkg_poln = TMath::RMS(3, bkg_poln);
    double fit_range_1[10] = {4274, 8971, 6453, 4383, 2550, 1544, 1056, 652, 441, 319}; // Fit range: [0.98, 1.2]
    double fit_range_2[10] = {4255, 8964, 6468, 4363, 2540, 1537, 1067, 649, 462, 321}; // Fit range: [0.98, 1.18]
    double fit_range_3[10] = {4241, 8894, 6522, 4350, 2549, 1545, 1042, 647, 465, 317}; // Fit range: [0.98, 1.22]
    double fit_range[3] = {fit_range_1[i-1], fit_range_2[i-1], fit_range_3[i-1]};
    double err_fit_range = TMath::RMS(3, fit_range);
    // double fit_func_1[10] = {1689.35, 1752.50, 2465.09, 6580.90, 8565.59, 2606.73, 3110.43, 2557.01, 2449.74, 160.43}; // Fit function: Voigt
    // double fit_func_2[10] = {2157.41, 2302.83, 3158.83, 7835.29, 10889.65, 3077.05, 3893.53, 3240.87, 3128.20, 217.33}; // Fit function: Breit-wigner
    // double fit_func_3[10] = {1608.70, 1734.22, 2367.27, 6358.65, 8255.21, 2452.65, 2940.46, 2476.66, 2322.49, 203.03}; // Fit function: Double Gaus
    // double fit_func[3] = {fit_func_1[i-1], fit_func_2[i-1], fit_func_3[i-1]};
    // double err_fit_func = TMath::RMS(3, fit_func);

    double err_sys = TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range);

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phi = h_tagged_flux->GetBinContent(1); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    if (lumi_phi <= 0)
    {
      //gphitxsec->RemovePoint(i - 1);
      continue;
    }
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_phi2pi = 1e9 * N_phit_data / (eff_phi * lumi_phi * T * br_phi);
    double dxsec_phi2pi_stat = xsec_phi2pi * (dN_phit_data / N_phit_data);
    double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys / N_phit_data)*(err_sys / N_phit_data) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double xsec_phi = N_phit / lumi_phi;
    // double dxsec_phi = dN_phit / lumi_phi;  
    gphitxsec->SetPoint(gphitxsec->GetN(), h2phit->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    gphitxsec->SetPointError(gphitxsec->GetN() - 1, 0, dxsec_phi2pi_stat);

    fprintf(table_xsec_phifo, "%0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h2phit_mc->GetXaxis()->GetBinCenter(i), h_beame_tru->GetBinContent(1), N_phit_mc, dN_phit_mc, N_phit_data, dN_phit_data, eff_phi*100, deff_phi*100, xsec_phi2pi, dxsec_phi2pi_stat, dxsec_phi2pi_sys);
    
    fprintf(table_xsec_phifo_sys, "%0.2f & %0.2f & %0.2f \\\\ \n", h2phit_mc->GetXaxis()->GetBinCenter(i), err_bkg_poln*100/N_bkg_poln_1[i-1], err_fit_range*100/fit_range_1[i-1]);
    // table_phi << std::setprecision(2) << std::fixed;
    // table_phi <<h2phit_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_phit_mc<<"$\\pm$"<<dN_phit_mc<<"&"<<N_phit_data<<"$\\pm$"<<dN_phit_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi2pi<<"$\\pm$"<<dxsec_phi2pi_stat<<"$\\pm$"<<dxsec_phi2pi_sys<<" \\\\"<< endl;
    }

    cgphit->cd();
    gphit->Draw("AP");
    gphit->Write(Form("h%s_gphit", name.Data()), TObject::kWriteDelete);
    cgphit_mc->cd();
    gphit_mc->Draw("AP");
    cgphiteff->cd();
    gphiteff->Draw("AP");
    cgphitxsec->cd();
    gphitxsec->Draw("AP");
    gphitxsec->Write(Form("h%s_gphitxsec", name.Data()), TObject::kWriteDelete);
    // int j =1;
    // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);

    // total Cross-sections
    TMultiGraph *mgxsec = new TMultiGraph();
    TCanvas *cmgxsec = new TCanvas("cmgxsec", "cmgxsec", 900, 600);
    TLegend *lmgxsec = new TLegend(0.86, 0.86, 0.98, 0.98);
    lmgxsec->SetTextSize(0.04);
    lmgxsec->SetBorderSize(0);
    cmgxsec->cd();
    cmgxsec->SetGrid();
    TGraphErrors *hdata_16_gphitxsec = (TGraphErrors *)outputfig->Get("hdata_16_gphitxsec");
    cout << " ***** hdata_16_gphitxsec = " << hdata_16_gphitxsec << endl;
    hdata_16_gphitxsec->SetMarkerColor(1);
    hdata_16_gphitxsec->SetMarkerSize(1.5);
    hdata_16_gphitxsec->SetLineColor(1);
    hdata_16_gphitxsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_16_gphitxsec);
    TGraphErrors *hdata_17_gphitxsec = (TGraphErrors *)outputfig->Get("hdata_17_gphitxsec");
    cout << " ***** hdata_17_gphitxsec = " << hdata_17_gphitxsec << endl;
    hdata_17_gphitxsec->SetMarkerColor(kBlue);
    hdata_17_gphitxsec->SetMarkerSize(1.5);
    hdata_17_gphitxsec->SetLineColor(kBlue);
    hdata_17_gphitxsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_17_gphitxsec);
    TGraphErrors *hdata_18_gphitxsec = (TGraphErrors *)outputfig->Get("hdata_18_gphitxsec");
    cout << " ***** hdata_18_gphitxsec = " << hdata_18_gphitxsec << endl;
    hdata_18_gphitxsec->SetMarkerColor(kRed);
    hdata_18_gphitxsec->SetMarkerSize(1.5);
    hdata_18_gphitxsec->SetLineColor(kRed);
    hdata_18_gphitxsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_18_gphitxsec);
    mgxsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. Momentum transfer; -t [GeV^{2}]; #sigma [nb]");
    mgxsec->Draw("AP");
    mgxsec->SetMinimum(0.);
    cmgxsec->Modified();
    lmgxsec->AddEntry(hdata_16_gphitxsec, "2016", "p");
    lmgxsec->AddEntry(hdata_17_gphitxsec, "2017", "p");
    lmgxsec->AddEntry(hdata_18_gphitxsec, "2018", "p");
    lmgxsec->Draw();
    cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.eps", "eps");
    cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.png", "png");
    cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgxsec.root", "root");

    // total Yields
    TMultiGraph *mgphit = new TMultiGraph();
    TCanvas *cmgphit = new TCanvas("cmgphit", "cmgphit", 900, 600);
    TLegend *lmgphit = new TLegend(0.86, 0.86, 0.98, 0.98);
    lmgphit->SetTextSize(0.04);
    lmgphit->SetBorderSize(0);
    cmgphit->cd();
    cmgphit->SetGrid();
    TGraphErrors *hdata_16_gphit = (TGraphErrors *)outputfig->Get("hdata_16_gphit");
    cout << " ***** hdata_16_gphit = " << hdata_16_gphit << endl;
    hdata_16_gphit->SetMarkerColor(1);
    hdata_16_gphit->SetMarkerSize(1.5);
    hdata_16_gphit->SetLineColor(1);
    hdata_16_gphit->SetMarkerStyle(20);
    mgphit->Add(hdata_16_gphit);
    TGraphErrors *hdata_17_gphit = (TGraphErrors *)outputfig->Get("hdata_17_gphit");
    cout << " ***** hdata_17_gphit = " << hdata_17_gphit << endl;
    hdata_17_gphit->SetMarkerColor(kBlue);
    hdata_17_gphit->SetMarkerSize(1.5);
    hdata_17_gphit->SetLineColor(kBlue);
    hdata_17_gphit->SetMarkerStyle(20);
    mgphit->Add(hdata_17_gphit);
    TGraphErrors *hdata_18_gphit = (TGraphErrors *)outputfig->Get("hdata_18_gphit");
    cout << " ***** hdata_18_gphit = " << hdata_18_gphit << endl;
    hdata_18_gphit->SetMarkerColor(kRed);
    hdata_18_gphit->SetMarkerSize(1.5);
    hdata_18_gphit->SetLineColor(kRed);
    hdata_18_gphit->SetMarkerStyle(20);
    mgphit->Add(hdata_18_gphit);
    mgphit->SetTitle("#phi(1020) Yield vs. Momentum transfer (data); -t [GeV^{2}]; N_{#phi}");
    mgphit->Draw("AP");
    mgphit->SetMinimum(0.);
    cmgphit->Modified();
    lmgphit->AddEntry(hdata_16_gphit, "2016", "p");
    lmgphit->AddEntry(hdata_17_gphit, "2017", "p");
    lmgphit->AddEntry(hdata_18_gphit, "2018", "p");
    lmgphit->Draw();
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphit.eps", "eps");
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphit.png", "png");
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cmgphit.root", "root");

    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit.root", name.Data()), "root");
    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit.eps", name.Data()), "eps");
    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit.png", name.Data()), "png");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit1.root", name.Data()), "root");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit1.eps", name.Data()), "eps");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_phit1.png", name.Data()), "png");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphit.root", name.Data()), "root");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphit.eps", name.Data()), "eps");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphit.png", name.Data()), "png");
    cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.root", "root");
    cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.eps", "eps");
    cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.png", "png");
    cphit1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_mc.root", "root");
    cphit1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_mc.eps", "eps");
    cphit1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_mc.png", "png");
    cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.root", "root");
    cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.eps", "eps");
    cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.png", "png");
    cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.root", "root");
    cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.eps", "eps");
    cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.png", "png");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphitxsec.root", name.Data()), "root");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphitxsec.eps", name.Data()), "eps");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/c%s_gphitxsec.png", name.Data()), "png");
*/
   /*
  // ======================================== Phi vs. (-t,E) ==================================================
  // root -l 'xsec_phifo.C+(100,6,20)'
  //++++++  Data
  TCanvas *cbeamet = new TCanvas("cbeamet", "cbeamet", 600, 400);
  cbeamet->cd();
  TH2D *h2beamet = new TH2D("h2beamet", "Data ;E_{#gamma} [GeV/c^{2}];-t [GeV/c]", 100, 6.3, 11.6, 100, 0.1, 2);
  tdata->Project("h2beamet", "-t_kin:beam_p4_kin.E()", "w8*("+cut+")");
  h2beamet->Draw("colz");

  TCanvas *cphit = new TCanvas("cphit", "cphit", 1500, 800);
  cphit->Divide(2, 3);

  TCanvas *cphit1[nt];

  TCanvas *cgphit = new TCanvas("cgphit", "cgphit", 1500, 800);
  cgphit->Divide(2, 3);
  TGraphErrors *gphit[nt];

  //++++++  mc
  TCanvas *cbeamet_mc = new TCanvas("cbeamet_mc", "cbeamet_mc", 600, 400);
  cbeamet_mc->cd();
  TH2D *h2beamet_mc = new TH2D("h2beamet_mc", "MC ;E_{#gamma} [GeV/c^{2}];-t [GeV/c]", 100, 6.3, 11.6, 100, 0.1, 2);
  tmc->Project("h2beamet_mc", "-t_kin:beam_p4_kin.E()", "w8*("+cut+")");
  h2beamet_mc->Draw("colz");

  TCanvas *cphit_mc = new TCanvas("cphit_mc", "cphit_mc", 1500, 800);
  cphit_mc->Divide(2, 3);

  TCanvas *cphit1_mc[nt];

  TCanvas *cgphit_mc = new TCanvas("cgphit_mc", "cgphit_mc", 1500, 800);
  cgphit_mc->Divide(2, 3);
  TGraphErrors *gphit_mc[nt];

  //++++++ Efficiency
  TCanvas *cgphiteff = new TCanvas("cgphiteff", "cgphiteff", 1500, 800);
  cgphiteff->Divide(2, 3);
  TGraphErrors *gphiteff[nt];

  //+++++++ Differential Cross-section
  TCanvas *cgphitxsec = new TCanvas("cgphitxsec", "cgphitxsec", 1500, 800);
  cgphitxsec->Divide(2, 3);
  TGraphErrors *gphitxsec[nt];

  double Egmin = 6.3;  //= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 11.6; //= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax - Egmin) / nt;
  double Eg1[nt];
  double Eg2[nt];

  for (int j = 1; j <= nt; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    //++++++  Data
    cphit->cd(j);
    TH2D *h2phit = new TH2D("h2phit", Form("%f<E_{#gamma}<%f (Data);-t [GeV/c];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), nt, 0.1, 2, n2k, 0.98, 1.2);
    tdata->Project("h2phit", "kpkm_mf:-t_kin", Form("w8*(kpkm_uni &&" + cut + "&& beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    h2phit->Draw("colz");

    cphit1[j] = new TCanvas(Form("cphit1_%d", j), Form("cphit1_%d", j), 1500, 800);
    cphit1[j]->Divide(5, 4);

    gphit[j] = new TGraphErrors(nt);
    gphit[j]->SetMarkerStyle(20);
    gphit[j]->SetMarkerSize(1.0);
    gphit[j]->SetMarkerColor(1);
    gphit[j]->SetMinimum(0.0);
    gphit[j]->SetTitle(Form("%f<E_{#gamma}<%f (Data); -t [GeV/c]; N_{#phi}", Eg1[j], Eg2[j]));

    //++++++  mc
    cphit_mc->cd(j);
    TH2D *h2phit_mc = new TH2D("h2phit_mc", Form("%f<E_{#gamma}<%f (MC);-t [GeV/c];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), nt, 0.1, 2, n2k, 0.98, 1.2);
    tmc->Project("h2phit_mc", "kpkm_mf:-t_kin", Form("w8*(kpkm_uni &&" + cut + "&& beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    h2phit_mc->Draw("colz");

    cphit1_mc[j] = new TCanvas(Form("cphit1_mc_%d", j), Form("cphit1_mc_%d", j), 1500, 800);
    cphit1_mc[j]->Divide(5, 4);

    gphit_mc[j] = new TGraphErrors(nt);
    gphit_mc[j]->SetMarkerStyle(20);
    gphit_mc[j]->SetMarkerSize(1.0);
    gphit_mc[j]->SetMarkerColor(1);
    gphit_mc[j]->SetMinimum(0.0);
    gphit_mc[j]->SetTitle(Form("%f<E_{#gamma}<%f (MC); -t [GeV/c]; N_{#phi}", Eg1[j], Eg2[j]));

    //++++++ Efficiency
    gphiteff[j] = new TGraphErrors(nt);
    gphiteff[j]->SetMarkerStyle(20);
    gphiteff[j]->SetMarkerSize(1.0);
    gphiteff[j]->SetMarkerColor(1);
    gphiteff[j]->SetMinimum(0.0);
    gphiteff[j]->SetTitle(Form("%f<E_{#gamma}<%f; -t [GeV/c]; #epsilon_{#phi}", Eg1[j], Eg2[j]));

    //+++++++ Differential Cross-section
    gphitxsec[j] = new TGraphErrors(nt);
    gphitxsec[j]->SetMarkerStyle(20);
    gphitxsec[j]->SetMarkerSize(1.0);
    gphitxsec[j]->SetMarkerColor(1);
    gphitxsec[j]->SetMinimum(0.0);
    gphitxsec[j]->SetTitle(Form("%f<E_{#gamma}<%f (flux normalized yield); -t [GeV/c]; Yield_{#phi}", Eg1[j], Eg2[j]));// d#sigma_{#phi}/dt [nb/s]

    for (int i = 1; i <= nt; ++i)
    {
      cout << i << " " << flush;

      // +++++++++++++++++++++++++ data  ++++++++++++++++++++
      cphit1[j]->cd(i);
      TH1D *hphit_py = h2phit->ProjectionY(Form("hphit_py_%d_%d", j, i), i, i);
      hphit_py->Draw("e");

      w.factory(Form("Voigtian::sig_phit(m_phit[%f,%f],mean_phit[1.016,1.024],width_phit[0.004],sigma_phit[0.001,0.1])", mkk_min, mkk_max)); //sigma_phit[0.001,0.01], mean_phit[1.016,1.022]
      w.factory("Chebychev::bkg_phit(m_phit,{c0_phit[-1,1], c1_phit[-1,1]})");//, c2_phit[-1,1]
      w.factory("SUM:::model_phit(nsig_phit[0,100000000]*sig_phit, nbkg_phit[0,100000000]*bkg_phit)"); //nsig[0,100000000]*sig2,
      w.var("m_phit")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phit("dh_phit", "dh_phit", *w.var("m_phit"), Import(*hphit_py));
      RooPlot *fr_phit = w.var("m_phit")->frame(Title(Form("-t = %f", h2phit->GetXaxis()->GetBinCenter(i))));
      w.pdf("model_phit")->fitTo(dh_phit);
      // result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phit.plotOn(fr_phit, RooFit::Name("ldh_phit"));
      w.pdf("model_phit")->plotOn(fr_phit, Components(*w.pdf("sig_phit")), LineColor(kRed), RooFit::Name("lsig_phit"));
      w.pdf("model_phit")->plotOn(fr_phit, Components(*w.pdf("bkg_phit")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phit"));
      w.pdf("model_phit")->plotOn(fr_phit, RooFit::Name("lmodel_phit"));
      fr_phit->Draw();
      TLegend *l_phit = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phit->SetFillColor(kWhite);
      l_phit->SetLineColor(kWhite);
      l_phit->AddEntry(fr_phit->findObject("lsig_phit"), Form("N_{Sig} = %.2f", w.var("nsig_phit")->getVal()), "l");
      l_phit->AddEntry(fr_phit->findObject("lbkg_phit"), Form("N_{Bkg} = %.2f", w.var("nbkg_phit")->getVal()), "l");
      l_phit->Draw();

      double N_phit = w.var("nsig_phit")->getVal();
      double dN_phit = w.var("nsig_phit")->getError();

      gphit[j]->SetPoint(i - 1, h2phit->GetXaxis()->GetBinCenter(i), N_phit);
      gphit[j]->SetPointError(i - 1, 0, dN_phit);

      // +++++++++++++++++++++++++ MC  ++++++++++++++++++++

      cphit1_mc[j]->cd(i);
      TH1D *hphit_mc_py = h2phit_mc->ProjectionY(Form("hphit_mc_py_%d_%d", j, i), i, i);
      hphit_mc_py->Draw("e");

      w.factory(Form("Voigtian::sig_phit_mc(m_phit_mc[%f,%f],mean_phit_mc[1.005,1.03],width_phit_mc[0.004],sigma_phit_mc[0.0001,0.8])", mkk_min, mkk_max)); //sigma_phit_mc[0.001,0.01], mean_phit_mc[1.016,1.022]
      w.factory("Chebychev::bkg_phit_mc(m_phit_mc,{c0_phit_mc[-1,1], c1_phit_mc[-1,1]})"); //, c2_phit_mc[-1,1]
      w.factory("SUM:::model_phit_mc(nsig_phit_mc[0,100000000]*sig_phit_mc, nbkg_phit_mc[0,100000000]*bkg_phit_mc)"); //nsig[0,100000000]*sig2,
      w.var("m_phit_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phit_mc("dh_phit_mc", "dh_phit_mc", *w.var("m_phit_mc"), Import(*hphit_mc_py));
      RooPlot *fr_phit_mc = w.var("m_phit_mc")->frame(Title(Form("-t = %f", h2phit_mc->GetXaxis()->GetBinCenter(i))));
      w.pdf("model_phit_mc")->fitTo(dh_phit_mc);
      // result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phit_mc.plotOn(fr_phit_mc, RooFit::Name("ldh_phit_mc"));
      w.pdf("model_phit_mc")->plotOn(fr_phit_mc, Components(*w.pdf("sig_phit_mc")), LineColor(kRed), RooFit::Name("lsig_phit_mc"));
      w.pdf("model_phit_mc")->plotOn(fr_phit_mc, Components(*w.pdf("bkg_phit_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phit_mc"));
      w.pdf("model_phit_mc")->plotOn(fr_phit_mc, RooFit::Name("lmodel_phit_mc"));
      fr_phit_mc->Draw();
      TLegend *l_phit_mc = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phit_mc->SetFillColor(kWhite);
      l_phit_mc->SetLineColor(kWhite);
      l_phit_mc->AddEntry(fr_phit_mc->findObject("lsig_phit_mc"), Form("N_{Sig} = %.2f", w.var("nsig_phit_mc")->getVal()), "l");
      l_phit_mc->AddEntry(fr_phit_mc->findObject("lbkg_phit_mc"), Form("N_{Bkg} = %.2f", w.var("nbkg_phit_mc")->getVal()), "l");
      l_phit_mc->Draw();

      double N_phit_mc = w.var("nsig_phit_mc")->getVal();
      double dN_phit_mc = w.var("nsig_phit_mc")->getError();

      gphit_mc[j]->SetPoint(i - 1, h2phit_mc->GetXaxis()->GetBinCenter(i), N_phit_mc);
      gphit_mc[j]->SetPointError(i - 1, 0, dN_phit_mc);

      // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
      double eff_phi = N_phit_mc / h_beame_tru->GetBinContent(j); // Efficiency = N_observed/N_generated
      gphiteff[j]->SetPoint(i - 1, h2phit_mc->GetXaxis()->GetBinCenter(i), eff_phi);
      gphiteff[j]->SetPointError(i - 1, 0, dN_phit_mc / h_beame_tru->GetBinContent(j));

      // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
      Int_t bin1 = h_tagged_flux->FindBin(Eg1[j]);
      Int_t bin2 = h_tagged_flux->FindBin(Eg2[j])-1;
      double lumi_phi = h_tagged_flux->Integral(bin1, bin2); // * 1.273; // Luminosity = N_gama * T ,  T = 1.26 barns^-1
      double br_phi = 0.489;
      // double xsec_phi = 1e9 * N_phit / (eff_phi * lumi_phi * br_phi);
      // double dxsec_phi = 1e9 * dN_phit / (eff_phi * lumi_phi * br_phi);
      double xsec_phi = N_phit / lumi_phi;
      double dxsec_phi = dN_phit / lumi_phi;     
      gphitxsec[j]->SetPoint(i - 1, h2phit->GetXaxis()->GetBinCenter(i), xsec_phi);
      gphitxsec[j]->SetPointError(i - 1, 0, dxsec_phi);
    
      ofs_xsec_phifo << " j = " << j << " i = " << i << " | N_phit = " << N_phit << " | N_phit_mc = " << N_phit_mc << " | h_beame_tru->GetBinContent(i) = " << h_beame_tru->GetBinContent(j) << " | h_tagged_flux->Integral(bin1, bin2) = " << h_tagged_flux->Integral(bin1, bin2) << " | eff_phi = " << eff_phi << " | xsec_phi = " << xsec_phi << endl;

    }

    cout << endl;

    cgphit->cd(j);
    gphit[j]->Draw("AP");
    cgphit_mc->cd(j);
    gphit_mc[j]->Draw("AP");
    cgphiteff->cd(j);
    gphiteff[j]->Draw("AP");
    cgphitxsec->cd(j);
    gphitxsec[j]->Draw("AP");

    cphit1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_%d.root", j), "root");
    cphit1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_%d.eps", j), "eps");
    cphit1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_%d.png", j), "png");
    cphit1_mc[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_mc_%d.root", j), "root");
    cphit1_mc[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit1_mc_%d.png", j), "png");
  }
  // int j =1;
  // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);

 cbeamet->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cbeamet.root", "root");
  cbeamet->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cbeamet.eps", "eps");
  cbeamet->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cbeamet.png", "png");
  cphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit.root", "root");
  cphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit.eps", "eps");
  cphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit.png", "png");
  cgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit.root", "root");
  cgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit.eps", "eps");
  cgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit.png", "png");
  cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.root", "root");
  cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.eps", "eps");
  cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cphit_mc.png", "png");
  cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.root", "root");
  cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.eps", "eps");
  cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphit_mc.png", "png");
  cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.root", "root");
  cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.eps", "eps");
  cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphiteff.png", "png");
  cgphitxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphitxsec.root", "root");
  cgphitxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphitxsec.eps", "eps");
  cgphitxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phifo/cgphitxsec.png", "png");
*/


   // table_phi << "\\hline" << endl;
   // table_phi << "\\end{tabularx}" << endl;
   // table_phi << "\\end{center}" << endl;
   // table_phi << "\\end{minipage}" << endl;
   // table_phi << "\\end{table}" << endl;
   // table_phi << "\\end{document}" << endl;
   // table_phi.close();
   // gSystem->Exec("pdflatex table_phi.tex");
  
}
