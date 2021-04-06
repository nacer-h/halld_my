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
#include "TGaxis.h"
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

void xsec_phi2pi(TString name, int n2k=100, int ne=10, int nt=10)//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  // TFile *fdata = NULL;
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", name.Data())); // _chi100 ,20200506/
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_phi2pi_%s_flat.root", name.Data())); // _chi100  ver03_11/
  TFile *ftru = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/thrown_phi2pi_%s.root", name.Data()));
  // TFile *fthr = new TFile(Form("output/input/tree_thrown_phi2pi_genr8_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TTree *tthr = (TTree *)fthr->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_11366_11555.root");
  if(name == "17") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_30274_31057.root");
  if(name == "18") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_40856_42577.root");
  if(name == "18l") fps = new TFile("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/flux_50677_51768.root");

  TFile *outputfig = new TFile("output/fig_phi2pi/xsec_phi2pi.root", "UPDATE");

  ofstream test;
  test.open ("test.txt");
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
  // table_phi <<"$E_{\\gamma}[ GeV^{2}]$ & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\varepsilon$ [\\%]& $\\sigma (nb)$ \\\\" << endl;
  // // table_phi << "\\alpha & \\beta & \\gamma \\" << endl;
  // table_phi <<"\\hline" << endl;

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.99, mkk_max = 1.2; // mkk_min = 0.99, mkk_max = 1.2, 1.15, 1.25
  double Eg_min = 6.5, Eg_max = 11.6;
  double Eg_step = (Eg_max-Eg_min)/ne;
  double t_min = 0.0, t_max = 4.0;
  double t_step = (t_max-t_min)/nt;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetLabelSize(0.07, "XYZ");
  // gStyle->SetTitleAlign(13);
/*
  // ======================================== Phi vs. Eg ===============================================

  //root -l 'xsec_phi2pi.C+("data_17",100,10,1)'

  FILE *table_xsec_phi2pi = fopen(Form("table_xsec_phi2pi_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{siunitx} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\small\n \\caption{Total cross-sections in photon energy bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ (GeV) & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\varepsilon$ [\\SI{70}{\\percent}] & $\\sigma$ (nb) \\\\ \n \\hline\n");

  FILE *table_xsec_phi2pi_sys = fopen(Form("table_xsec_phi2pi_sys_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi_sys,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{makecell} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ (GeV) & \\thead{Bkg deg\\\\(3,4,5)}& \\thead{Fit range\\\\(0.99; 1,15, 1.2, 1.25)} & \\thead{bining\\\\(90, 100, 110)} & \\thead{Accidentals\\\\(2, 4, 8)} & \\thead{$\\Delta T_{RF-TOF}$ proton\\\\($\\pm$0.2, $\\pm$0.3, $\\pm$0.4)} & \\thead{$\\chi^{2}$ KinFit\\\\(45, 55, 65)} & \\thead{$MM^{2}$\\\\($\\pm$0.025, $\\pm$0.035, $\\pm$0.045)} \\\\ \n \\hline\n");

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");
  h_tagged_flux->Write(Form("h%s_tagged_flux", name.Data()), TObject::kWriteDelete);

  // +++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown"); // new TH1F("h_beame_tru", "MC truth; E_{#gamma} (GeV); Counts", ne, Eg_min, Eg_max)
  // tthr->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");
  h_beame_tru->Write(Form("h%s_beame_tru", name.Data()), TObject::kWriteDelete);

  // +++++++++++++++++++++ Data  +++++++++++++++++++++ 
  TCanvas *cphie = new TCanvas("cphie", "cphie", 900, 600);
  cphie->cd();
  TH2D *h2phie = new TH2D("h2phie", ";E_{#gamma} (GeV);m_{K^{+}K^{-}} (GeV/c^{2});Counts", ne, Eg_min, Eg_max, n2k, mkk_min, mkk_max);
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");// && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(kpkm_uni && abs(rf_dt)<1.5*4.008), w4*(kpkm_uni && abs(rf_dt)<2.5*4.008), w8*(kpkm_uni)
  // TH2D *h2phie = (TH2D *)fdata->Get("h2_PhiMassVsBeamE_KinFit"); 
  // h2phie->RebinX(10);
  cout << " ***** h2phie = " << h2phie << endl;
  h2phie->Draw("colz");

  TCanvas *cphie1 = new TCanvas("cphie1", "cphie1", 900, 1200);
  cphie1->Divide(2, 5);
  TCanvas *cgphie = new TCanvas("cgphie", "cgphie", 900, 600);
  TGraphErrors *gphie = new TGraphErrors();
  gphie->SetMarkerStyle(20);
  gphie->SetMarkerSize(1.5);
  gphie->SetMinimum(0.0);
  gphie->SetTitle("; E_{#gamma} (GeV); N_{#phi #pi^{+} #pi^{-}}");

  //  +++++++++++++++++++++ Monte Carlo +++++++++++++++++++++ 
  TCanvas *cphie_mc = new TCanvas("cphie_mc", "cphie_mc", 900, 600);
  cphie_mc->cd();
  TH2D *h2phie_mc = new TH2D("h2phie_mc", ";E_{#gamma} (GeV);m_{K^{+}K^{-}} (GeV/c^{2});Counts", ne, Eg_min, Eg_max, n2k, mkk_min, mkk_max);//0.98, 1.2
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");
  // TH2D *h2phie_mc = (TH2D *)fmc->Get("h2_PhiMassVsBeamE_KinFit");
  // h2phie_mc->RebinX(10);
  cout << " ***** h2phie_mc = " << h2phie_mc << endl;
  h2phie_mc->Draw("colz");

  TCanvas *cphie1_mc = new TCanvas("cphie1_mc", "cphie1_mc", 900, 1200);
  cphie1_mc->Divide(2, 5);
  TCanvas *cgphie_mc = new TCanvas("cgphie_mc", "cgphie_mc", 900, 600);
  TGraphErrors *gphie_mc = new TGraphErrors();
  gphie_mc->SetMarkerStyle(20);
  gphie_mc->SetMarkerSize(1.5);
  gphie_mc->SetMinimum(0.0);
  gphie_mc->SetTitle("; E_{#gamma} (GeV); N_{#phi #pi^{+} #pi^{-}}");

  // Efficiency
  TCanvas *cgphieeff = new TCanvas("cgphieeff", "cgphieeff", 900, 600);
  TGraphErrors *gphieeff = new TGraphErrors();
  gphieeff->SetMarkerStyle(20);
  gphieeff->SetMarkerSize(1.5);
  gphieeff->SetMarkerColor(kBlue);
  gphieeff->SetMinimum(0.0);
  gphieeff->SetTitle("; E_{#gamma} (GeV); #varepsilon (%)");

  // Cross-section
  TCanvas *cgphiexsec = new TCanvas("cgphiexsec", "cgphiexsec", 900, 600);
  TGraphErrors *gphiexsec = new TGraphErrors();
  gphiexsec->SetMarkerStyle(20);
  gphiexsec->SetMarkerSize(1.5);
  gphiexsec->SetMinimum(0.0);
  gphiexsec->SetTitle("; E_{#gamma} (GeV); #sigma (nb)"); //#phi(1020) flux normalized yield, Yield_{#phi}
  // double Eg = 0;
  double testxsec[ne];
  for (int i = 1; i <= ne; ++i)
  {
    // cout << i << " " << flush;
    double Eg1 = h2phie_mc->GetXaxis()->GetBinLowEdge(i);
    double Eg2 = h2phie_mc->GetXaxis()->GetBinLowEdge(i)+Eg_step;
    double Eg = h2phie_mc->GetXaxis()->GetBinCenter(i);
    // if(i == 1) Eg = Eg1;
    // else Eg = Eg2;

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphie1_mc->cd(i);
    TH1D *hphie_mc_py = h2phie_mc->ProjectionY(Form("hphie_mc_py_%d", i), i, i);
    hphie_mc_py->SetTitle(Form("%.2f<E_{#gamma}<%.2f (GeV);m_{K^{+}K^{-}} (GeV/c^{2});Counts",Eg1,Eg2));
    // hphie_mc_py->SetLabelSize(.06, "xyz");
    hphie_mc_py->Draw("e");
   
    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);//pol4(4)
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1); // 1,1,1,1,1
    fsb_mc->SetParLimits(1, 1.015, 1.022); //fsb_mc->FixParameter(1, 1.019);
    // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.0001, 0.01);
    fsb_mc->FixParameter(3, 0.0042);   //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    fb_mc->SetLineColor(28);
    fb_mc->SetLineStyle(2);

    hphie_mc_py->Fit("fsb_mc", "", "", mkk_min, mkk_max);
    double par_mc[fsb_mc->GetNpar()]; //6
    fsb_mc->GetParameters(&par_mc[0]);
    fs_mc->SetParameters(&par_mc[0]);
    fb_mc->SetParameters(&par_mc[4]); //4

    fs_mc->Draw("same");
    fb_mc->Draw("same");

    double N_phie_mc = fs_mc->Integral(mkk_min, mkk_max) / hphie_mc_py->GetBinWidth(1);
    double dN_phie_mc = N_phie_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

    TLatex lat_phie_mc;
    lat_phie_mc.SetTextSize(0.07);
    lat_phie_mc.SetTextAlign(13); //align at top
    lat_phie_mc.SetNDC();
    lat_phie_mc.SetTextColor(kBlue);
    lat_phie_mc.DrawLatex(0.5, 0.80, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
    lat_phie_mc.DrawLatex(0.5, 0.72, Form("N_{sig} = %0.2f#pm%0.2f", N_phie_mc, dN_phie_mc));
    lat_phie_mc.DrawLatex(0.5, 0.64, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
    lat_phie_mc.DrawLatex(0.5, 0.56, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phie_mc.DrawLatex(0.5, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

    gphie_mc->SetPoint(i - 1, Eg, N_phie_mc);
    gphie_mc->SetPointError(i - 1, 0, dN_phie_mc);

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphie1->cd(i);
    TH1D *hphie_py = h2phie->ProjectionY(Form("hphie_py_%d", i), i, i);
    hphie_py->SetTitle(Form("%.2f<E_{#gamma}<%.2f (GeV);m_{K^{+}K^{-}} (GeV/c^{2});Counts",Eg1, Eg2));
    // hphie_py->SetLabelSize(.06, "xyz");
    hphie_py->Draw("e");

    // +++++++++ Voigt + pol4
    TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);//pol4(4)
    fsb_data->SetLineColor(2);
    fsb_data->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1);  // 1,1,1,1,1
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); // 1.018, 1.021
    fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); // mean+0.25*sigma
    fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); //
    fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); // 0.001,0.01

    TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_data->SetLineColor(4);

    TF1 *fb_data = new TF1("fb_data", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    fb_data->SetLineColor(28);
    fb_data->SetLineStyle(2);

    hphie_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    double par_data[fsb_data->GetNpar()]; //6
    fsb_data->GetParameters(&par_data[0]);
    fs_data->SetParameters(&par_data[0]);
    fb_data->SetParameters(&par_data[4]); //4

    fs_data->Draw("same");
    fb_data->Draw("same");

    double N_phie_data = fs_data->Integral(mkk_min, mkk_max) / hphie_py->GetBinWidth(1);
    double dN_phie_data = N_phie_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    TLatex lat_phie_data;
    lat_phie_data.SetTextSize(0.07);
    lat_phie_data.SetTextAlign(13); //align at top
    lat_phie_data.SetNDC();
    lat_phie_data.SetTextColor(kBlue);
    lat_phie_data.DrawLatex(0.5, 0.80, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phie_data.DrawLatex(0.5, 0.72, Form("N_{sig} = %0.2f#pm%0.2f", N_phie_data, dN_phie_data));
    lat_phie_data.DrawLatex(0.5, 0.64, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phie_data.DrawLatex(0.5, 0.30, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phie_data.DrawLatex(0.5, 0.23, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

    gphie->SetPoint(i - 1, Eg, N_phie_data);
    gphie->SetPointError(i - 1, 0, dN_phie_data);

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    // if (N_phie_mc <= 0)
    // {
    //   //gphiexsec->RemovePoint(i - 1);
    //   continue;
    // }
    double eff_phi = N_phie_mc / h_beame_tru->GetBinContent(i); // Efficiency = N_observed/N_generated
    double deff_phi = eff_phi * (dN_phie_mc / N_phie_mc);
    gphieeff->SetPoint(i-1, Eg, eff_phi*100);
    gphieeff->SetPointError(i-1, 0, deff_phi*100); //->GetBinError(i)
    // gphieeff->SetPoint(gphieeff->GetN(), h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phi*100);
    // gphieeff->SetPointError(gphieeff->GetN()-1, 0, deff_phi*100); //->GetBinError(i)

    // ============================ Systematic Errors ============================

    // ------------------------ 2016

    double optimal_val_16[10] = {80.72, 71.29, 78.90, 64.05, 59.32, 61.04, 58.58, 58.17, 60.35, 46.74};
    // Bkg pol order
    double bkg_poln_16_1[10] = {80.38, 74.24, 81.46, 66.16, 59.50, 63.05, 61.68, 60.56, 62.79, 46.60}; // pol3
    double bkg_poln_16_2[10] = {79.58, 71.40, 79.17, 64.12, 59.68, 61.12, 58.68, 58.26, 60.72, 47.08}; // pol5
    // fit range
    double fit_range_16_1[10] = {81.57, 73.92, 79.91, 65.61, 61.27, 61.39, 60.70, 59.48, 62.35, 49.39}; // [0.99, 1.15]
    double fit_range_16_2[10] = {80.62, 71.12, 78.81, 64.76, 60.05, 61.09, 59.83, 59.54, 60.62, 49.31}; // [0.99, 1.25]
    // Binning
    double binning_16_1[10] = {81.58, 71.52, 78.81, 63.92, 59.25, 59.71, 58.44, 57.73, 60.04, 46.98}; // 90
    double binning_16_2[10] = {82.06, 70.17, 78.89, 64.25, 59.18, 61.04, 59.01, 57.85, 60.37, 46.33}; //110
    // beam bunches
    double beambun_16_1[10] = {80.91, 69.83, 76.24, 63.25, 59.07, 60.00, 58.93, 58.70, 59.20, 46.42}; // w2
    double beambun_16_2[10] = {79.56, 70.67, 77.11, 64.51, 59.44, 60.40, 59.38, 58.52, 60.09, 46.78}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_16_opt[10] = {80.76, 71.24, 79.41, 64.09, 59.67, 61.00, 58.82, 58.61, 60.34, 46.91};
    // abs(p_dttof)<0.3
    double p_dttof_16_1[10] = {79.51, 70.93, 78.37, 64.14, 58.90, 60.61, 58.03, 58.72, 60.13, 46.89}; // 0.2
    double p_dttof_16_2[10] = {81.62, 72.81, 80.99, 64.84, 60.37, 61.07, 59.68, 58.94, 60.46, 47.22}; // 0.4
    // kin_chisq<55
    double kin_chisq_16_1[10] = {76.14, 68.10, 76.67, 63.84, 57.73, 60.84, 58.10, 58.73, 59.51, 46.47}; // 45
    double kin_chisq_16_2[10] = {81.22, 76.21, 81.14, 65.14, 61.45, 60.55, 58.11, 58.68, 62.52, 47.94}; // 65
    // abs(mm2)<0.035
    double mm2_16_1[10] = {80.06, 71.31, 79.21, 63.03, 59.34, 59.61, 59.01, 58.39, 59.14, 46.87}; // 0.025
    double mm2_16_2[10] = {81.29, 70.02, 79.65, 63.80, 59.70, 60.75, 58.99, 58.73, 59.72, 46.94}; // 0.045

    // ------------------------ 2017

    double optimal_val_17[10] = {65.47, 60.77, 56.99, 54.84, 52.70, 54.45, 48.81, 48.64, 45.55, 41.70};
    // Bkg pol order
    double bkg_poln_17_1[10] = {65.34, 60.69, 57.04, 54.78, 52.64, 54.39, 48.73, 48.58, 45.49, 41.64}; // pol3
    double bkg_poln_17_2[10] = {65.51, 60.86, 57.17, 54.92, 52.77, 54.63, 48.89, 48.71, 45.62, 41.77}; // pol5
    // fit range
    double fit_range_17_1[10] = {66.59, 62.23, 58.39, 55.98, 53.79, 55.62, 50.25, 50.07, 47.25, 43.59}; // [0.99, 1.15]
    double fit_range_17_2[10] = {65.74, 61.29, 57.57, 55.20, 52.96, 54.90, 49.74, 49.35, 46.38, 42.44}; // [0.99, 1.25]
    // Binning
    double binning_17_1[10] = {65.11, 60.80, 57.11, 54.88, 52.62, 54.54, 48.81, 48.65, 45.53, 41.93}; // 90
    double binning_17_2[10] = {65.28, 60.85, 57.18, 54.81, 52.55, 54.44, 48.77, 48.77, 45.54, 41.94}; // 110
    // beam bunches
    double beambun_17_1[10] = {65.12, 60.31, 56.96, 54.35, 52.32, 54.80, 48.55, 48.60, 45.77, 41.54}; // w2
    double beambun_17_2[10] = {65.41, 60.39, 56.86, 54.80, 52.67, 54.49, 48.71, 48.50, 45.68, 41.59}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_17_opt[10] = {65.46, 60.94, 57.03, 54.96, 52.75, 54.53, 48.91, 48.61, 45.67, 41.93};
    // abs(p_dttof)<0.3
    double p_dttof_17_1[10] = {64.36, 60.18, 56.21, 54.07, 52.00, 53.97, 48.30, 47.97, 45.16, 41.23}; // 0.2
    double p_dttof_17_2[10] = {66.24, 61.63, 57.69, 55.41, 53.11, 54.77, 49.19, 48.81, 45.98, 42.54}; // 0.4
    // kin_chisq<55
    double kin_chisq_17_1[10] = {64.14, 59.26, 56.53, 54.13, 51.75, 53.96, 47.87, 48.06, 43.98, 40.57}; // 45
    double kin_chisq_17_2[10] = {66.44, 61.99, 57.88, 55.73, 53.40, 55.13, 49.49, 48.73, 46.20, 42.62}; // 65
    // abs(mm2)<0.035
    double mm2_17_1[10] = {65.58, 60.88, 56.69, 54.61, 52.46, 54.14, 48.35, 47.95, 45.01, 41.61}; // 0.025
    double mm2_17_2[10] = {65.79, 61.12, 57.02, 55.17, 52.69, 54.80, 49.14, 48.74, 45.94, 42.16}; // 0.045

    // ------------------------ 2018 Spring

    double optimal_val_18[10] = {64.18, 59.59, 59.34, 54.55, 53.26, 52.55, 49.31, 46.60, 44.85, 43.89};
    // Bkg pol order
    double bkg_poln_18_1[10] = {64.07, 59.48, 59.23, 54.44, 53.16, 52.45, 49.19, 46.49, 44.75, 43.66}; // pol3
    double bkg_poln_18_2[10] = {64.31, 59.82, 59.45, 54.67, 53.36, 52.65, 49.43, 46.72, 44.96, 44.02}; // pol5
    // fit range
    double fit_range_18_1[10] = {63.98, 59.68, 59.72, 54.49, 53.65, 53.06, 50.27, 47.52, 45.87, 45.80}; // [0.99, 1.15]
    double fit_range_18_2[10] = {64.59, 60.39, 60.07, 55.39, 53.76, 53.34, 50.55, 47.72, 46.00, 44.95};  // [0.99, 1.25]
    // Binning
    double binning_18_1[10] = {64.30, 59.94, 59.21, 54.42, 53.12, 52.65, 49.42, 46.52, 44.90, 43.89}; // 90
    double binning_18_2[10] = {63.99, 59.76, 59.33, 54.44, 53.17, 52.53, 49.36, 46.58, 44.80, 43.98}; // 110
    // beam bunches
    double beambun_18_1[10] = {64.28, 59.15, 59.56, 54.08, 52.63, 52.69, 49.48, 46.39, 44.72, 43.83}; // w2
    double beambun_18_2[10] = {64.20, 59.44, 59.28, 54.33, 52.94, 52.56, 49.39, 46.71, 44.86, 43.77}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18_opt[10] = {63.44, 58.94, 58.31, 53.79, 52.41, 51.57, 48.69, 45.88, 44.22, 43.25};
    // abs(p_dttof)<0.3
    double p_dttof_18_1[10] = {62.71, 58.33, 57.62, 53.26, 51.82, 51.07, 48.27, 45.51, 43.96, 42.93}; // 0.2
    double p_dttof_18_2[10] = {63.84, 59.17, 58.65, 54.10, 52.59, 51.88, 48.90, 46.04, 44.50, 43.37}; // 0.4
    // kin_chisq<55
    double kin_chisq_18_1[10] = {60.11, 56.59, 56.14, 51.30, 49.88, 49.34, 46.54, 43.28, 41.19, 40.89}; // 45
    double kin_chisq_18_2[10] = {64.71, 60.37, 59.84, 55.61, 54.38, 53.75, 50.88, 48.20, 46.43, 44.97}; // 65
    // abs(mm2)<0.035
    double mm2_18_1[10] = {63.49, 58.76, 58.46, 53.78, 52.36, 51.34, 48.94, 45.72, 43.95, 43.15}; // 0.025
    double mm2_18_2[10] = {63.56, 58.99, 58.34, 53.87, 52.38, 51.66, 48.93, 45.95, 44.32, 43.36}; // 0.045

    // ------------------------ 2018 Fall

    double optimal_val_18l[10] = {63.89, 64.25, 63.18, 60.87, 60.32, 61.48, 58.54, 57.41, 54.69, 53.00};
    // Bkg pol order
    double bkg_poln_18l_1[10] = {63.80, 64.30, 63.10, 60.78, 60.25, 61.40, 58.46, 57.12, 54.60, 52.91}; // pol3
    double bkg_poln_18l_2[10] = {63.98, 64.45, 63.26, 60.96, 60.41, 61.57, 58.63, 57.50, 54.78, 53.11}; // pol5
    // fit range
    double fit_range_18l_1[10] = {63.41, 63.88, 63.58, 61.31, 61.07, 63.06, 60.20, 59.25, 56.49, 55.42}; // [0.99, 1.15]
    double fit_range_18l_2[10] = {64.26, 64.59, 63.60, 61.37, 60.76, 61.97, 59.02, 58.36, 55.58, 54.21}; // [0.99, 1.25]
    // Binning
    double binning_18l_1[10] = {63.90, 64.24, 63.29, 60.88, 60.31, 61.56, 58.46, 57.41, 54.76, 53.13}; // 90
    double binning_18l_2[10] = {63.69, 64.34, 63.23, 60.77, 60.27, 61.38, 58.47, 57.35, 54.74, 52.91}; // 110
    // beam bunches
    double beambun_18l_1[10] = {64.28, 64.25, 62.83, 60.05, 59.66, 61.17, 58.22, 57.27, 54.66, 52.67}; // w2
    double beambun_18l_2[10] = {63.89, 64.17, 62.92, 60.55, 60.14, 61.38, 58.37, 57.28, 54.73, 52.90}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18l_opt[10] = {63.62, 64.12, 63.12, 60.85, 60.27, 61.52, 58.58, 57.46, 54.77, 53.09};
    // abs(p_dttof)<0.3
    double p_dttof_18l_1[10] = {62.33, 63.28, 62.27, 60.05, 59.63, 60.89, 57.87, 56.89, 54.29, 52.57}; // 0.2
    double p_dttof_18l_2[10] = {63.99, 64.61, 63.51, 61.23, 60.59, 61.90, 58.95, 57.81, 55.16, 53.33}; // 0.4
    // kin_chisq<55
    double kin_chisq_18l_1[10] = {63.16, 63.81, 62.54, 59.88, 59.35, 60.56, 57.49, 56.37, 53.51, 51.83}; // 45
    double kin_chisq_18l_2[10] = {64.01, 64.39, 63.38, 61.46, 60.92, 62.30, 59.69, 58.11, 55.78, 53.80}; // 65
    // abs(mm2)<0.035
    double mm2_18l_1[10] = {63.61, 64.21, 62.97, 60.60, 60.02, 61.39, 58.32, 56.95, 54.33, 52.39}; // 0.025
    double mm2_18l_2[10] = {63.67, 64.17, 63.24, 60.98, 60.36, 61.83, 58.68, 57.68, 55.00, 53.32}; // 0.045

    // --------------- total  

    double optimal_val[10], bkg_poln_1[10], bkg_poln_2[10], fit_range_1[10], fit_range_2[10], binning_1[10], binning_2[10], beambun_1[10], beambun_2[10], cutsvar_opt[10], p_dttof_1[10], p_dttof_2[10], kin_chisq_1[10], kin_chisq_2[10], mm2_1[10], mm2_2[10];

    if (name == "16")
    {
      optimal_val[i-1] = optimal_val_16[i-1];    
      bkg_poln_1[i-1] = bkg_poln_16_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_16_2[i-1];
      fit_range_1[i-1] = fit_range_16_1[i-1];
      fit_range_2[i-1] = fit_range_16_2[i-1];
      binning_1[i-1] = binning_16_1[i-1];
      binning_2[i-1] = binning_16_2[i-1];
      beambun_1[i-1] = beambun_16_1[i-1];
      beambun_2[i-1] = beambun_16_2[i-1];
      cutsvar_opt[i-1] = cutsvar_16_opt[i-1];
      p_dttof_1[i-1] = p_dttof_16_1[i-1];
      p_dttof_2[i-1] = p_dttof_16_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_16_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_16_2[i-1];
      mm2_1[i-1] = mm2_16_1[i-1];
      mm2_2[i-1] = mm2_16_2[i-1];
    }
    if (name == "17")
    {
      optimal_val[i-1] = optimal_val_17[i-1];    
      bkg_poln_1[i-1] = bkg_poln_17_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_17_2[i-1];
      fit_range_1[i-1] = fit_range_17_1[i-1];
      fit_range_2[i-1] = fit_range_17_2[i-1];
      binning_1[i-1] = binning_17_1[i-1];
      binning_2[i-1] = binning_17_2[i-1];
      beambun_1[i-1] = beambun_17_1[i-1];
      beambun_2[i-1] = beambun_17_2[i-1];
      cutsvar_opt[i-1] = cutsvar_17_opt[i-1];
      p_dttof_1[i-1] = p_dttof_17_1[i-1];
      p_dttof_2[i-1] = p_dttof_17_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_17_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_17_2[i-1];
      mm2_1[i-1] = mm2_17_1[i-1];
      mm2_2[i-1] = mm2_17_2[i-1];
    }
    if (name == "18")
    {
      optimal_val[i-1] = optimal_val_18[i-1];    
      bkg_poln_1[i-1] = bkg_poln_18_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_18_2[i-1];
      fit_range_1[i-1] = fit_range_18_1[i-1];
      fit_range_2[i-1] = fit_range_18_2[i-1];
      binning_1[i-1] = binning_18_1[i-1];
      binning_2[i-1] = binning_18_2[i-1];
      beambun_1[i-1] = beambun_18_1[i-1];
      beambun_2[i-1] = beambun_18_2[i-1];
      cutsvar_opt[i-1] = cutsvar_18_opt[i-1];
      p_dttof_1[i-1] = p_dttof_18_1[i-1];
      p_dttof_2[i-1] = p_dttof_18_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_18_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_18_2[i-1];
      mm2_1[i-1] = mm2_18_1[i-1];
      mm2_2[i-1] = mm2_18_2[i-1];
    }
    if (name == "18l")
    {
      optimal_val[i-1] = optimal_val_18l[i-1];    
      bkg_poln_1[i-1] = bkg_poln_18l_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_18l_2[i-1];
      fit_range_1[i-1] = fit_range_18l_1[i-1];
      fit_range_2[i-1] = fit_range_18l_2[i-1];
      binning_1[i-1] = binning_18l_1[i-1];
      binning_2[i-1] = binning_18l_2[i-1];
      beambun_1[i-1] = beambun_18l_1[i-1];
      beambun_2[i-1] = beambun_18l_2[i-1];
      cutsvar_opt[i-1] = cutsvar_18l_opt[i-1];
      p_dttof_1[i-1] = p_dttof_18l_1[i-1];
      p_dttof_2[i-1] = p_dttof_18l_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_18l_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_18l_2[i-1];
      mm2_1[i-1] = mm2_18l_1[i-1];
      mm2_2[i-1] = mm2_18l_2[i-1];
    }

    // total
    double arr_bkg_poln[3] = {optimal_val[i-1], bkg_poln_1[i-1], bkg_poln_2[i-1]};
    double err_bkg_poln = TMath::RMS(3, arr_bkg_poln)/optimal_val[i-1];
    double arr_fit_range[3] = {optimal_val[i-1], fit_range_1[i-1], fit_range_2[i-1]};
    double err_fit_range = TMath::RMS(3, arr_fit_range)/optimal_val[i-1];
    double arr_binning[3] = {optimal_val[i-1], binning_1[i-1], binning_2[i-1]};
    double err_binning = TMath::RMS(3, arr_binning)/optimal_val[i-1];
    double arr_beambun[3] = {optimal_val[i-1], beambun_1[i-1], beambun_2[i-1]};
    double err_beambun = TMath::RMS(3, arr_beambun)/optimal_val[i-1];
    double arr_p_dttof[3] = {cutsvar_opt[i-1], p_dttof_1[i-1], p_dttof_2[i-1]};
    double err_p_dttof = TMath::RMS(3, arr_p_dttof)/optimal_val[i-1];
    double arr_kin_chisq[3] = {cutsvar_opt[i-1], kin_chisq_1[i-1], kin_chisq_2[i-1]};
    double err_kin_chisq = TMath::RMS(3, arr_kin_chisq)/optimal_val[i-1];
    double arr_mm2[3] = {cutsvar_opt[i-1], mm2_1[i-1], mm2_2[i-1]};
    double err_mm2 = TMath::RMS(3, arr_mm2)/optimal_val[i-1];

    double err_sys = TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range + err_binning*err_binning + err_beambun*err_beambun + err_p_dttof*err_p_dttof + err_kin_chisq*err_kin_chisq + err_mm2*err_mm2);

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phi = h_tagged_flux->GetBinContent(i); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    // if (lumi_phi <= 0)
    // {
    //   //gphiexsec->RemovePoint(i - 1);
    //   continue;
    // }
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_phi2pi = 1e9 * N_phie_data / (eff_phi * lumi_phi * T * br_phi);
    double dxsec_phi2pi_stat = xsec_phi2pi * (dN_phie_data / N_phie_data);
    double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys*err_sys) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys*err_sys) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((dN_phie_data_sys / N_phie_data)*(dN_phie_data_sys / N_phie_data) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    double dxsec_phi2pi_tot = TMath::Sqrt(dxsec_phi2pi_stat*dxsec_phi2pi_stat + dxsec_phi2pi_sys*dxsec_phi2pi_sys);
    // double xsec_phi = N_phie / lumi_phi;
    // double dxsec_phi = dN_phie / lumi_phi;  
    gphiexsec->SetPoint(i-1, Eg, xsec_phi2pi);
    gphiexsec->SetPointError(i-1, 0, dxsec_phi2pi_tot);
    // gphiexsec->SetPoint(gphiexsec->GetN(), h2phie->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    // gphiexsec->SetPointError(gphiexsec->GetN() - 1, 0, dxsec_phi2pi_stat);

    fprintf(table_xsec_phi2pi, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", Eg1, Eg2, h_beame_tru->GetBinContent(i), N_phie_mc, dN_phie_mc, N_phie_data, dN_phie_data, eff_phi*100, deff_phi*100, xsec_phi2pi, dxsec_phi2pi_stat, dxsec_phi2pi_sys);

    fprintf(table_xsec_phi2pi_sys, "%0.2f - %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", Eg1, Eg2, err_bkg_poln * 100, err_fit_range * 100, err_binning * 100, err_beambun * 100, err_p_dttof * 100, err_kin_chisq * 100, err_mm2 * 100);

    //err_bkg_poln * 100 / optimal_val[i - 1]
    testxsec[i-1] = xsec_phi2pi;
    // test << std::setprecision(2) << std::fixed;
    // test << "Eg1 - Eg2 = " << Eg1 << "-" << Eg2 << " | xsec_phi2pi = " << xsec_phi2pi << endl;

    // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_phie_mc<<"$\\pm$"<<dN_phie_mc<<"&"<<N_phie_data<<"$\\pm$"<<dN_phie_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi2pi<<"$\\pm$"<<dxsec_phi2pi_stat<<"$\\pm$"<<dxsec_phi2pi_sys<<" \\\\"<< endl;
  }

  test << std::setprecision(2) << std::fixed;
  test << testxsec[0] << ", " << testxsec[1] << ", " << testxsec[2] << ", " << testxsec[3] << ", " << testxsec[4] << ", " << testxsec[5] << ", " << testxsec[6] << ", " << testxsec[7] << ", " << testxsec[8] << ", " << testxsec[9] << endl;

  cgphie->cd();
  gphie->Draw("AP");
  gphie->Write(Form("h%s_gphie", name.Data()), TObject::kWriteDelete);
  cgphie_mc->cd();
  gphie_mc->Draw("AP");
  gphie_mc->Write(Form("h%s_gphie_mc", name.Data()), TObject::kWriteDelete);
  cgphieeff->cd();
  gphieeff->Draw("AP");
  gphieeff->Write(Form("h%s_gphieeff", name.Data()), TObject::kWriteDelete);
  cgphiexsec->cd();
  gphiexsec->Draw("AP");
  gphiexsec->Write(Form("h%s_gphiexsec", name.Data()), TObject::kWriteDelete);
  // int j =1;
  // gphie->Write(Form("grphie_%d", j), TObject::kWriteDelete);

  // *********** total Beam flux
  TCanvas *ctot_tagged_flux = new TCanvas("ctot_tagged_flux", "ctot_tagged_flux", 900, 600);
  TLegend *ltot_tagged_flux = new TLegend(0.84, 0.80, 0.98, 0.98);
  ltot_tagged_flux->SetTextSize(0.06);
  ltot_tagged_flux->SetBorderSize(0);
  ctot_tagged_flux->cd();
  ctot_tagged_flux->SetGrid();
  TH1F *h16_tagged_flux = (TH1F *)outputfig->Get("h16_tagged_flux");
  cout << " ***** h16_tagged_flux = " << h16_tagged_flux << endl;
  // h16_tagged_flux->SetMarkerColor(kBlue);
  h16_tagged_flux->SetMarkerSize(1.5);
  h16_tagged_flux->SetMarkerStyle(kFullSquare);
  TH1F *h17_tagged_flux = (TH1F *)outputfig->Get("h17_tagged_flux");
  cout << " ***** h17_tagged_flux = " << h17_tagged_flux << endl;
  h17_tagged_flux->SetMarkerSize(1.5);
  h17_tagged_flux->SetMarkerColor(kBlue);
  h17_tagged_flux->SetMarkerStyle(kFullCircle);
  TH1F *h18_tagged_flux = (TH1F *)outputfig->Get("h18_tagged_flux");
  cout << " ***** h18_tagged_flux = " << h18_tagged_flux << endl;
  h18_tagged_flux->SetMarkerColor(kRed);
  h18_tagged_flux->SetMarkerSize(1.5);
  h18_tagged_flux->SetMarkerStyle(kFullTriangleUp);
  TH1F *h18l_tagged_flux = (TH1F *)outputfig->Get("h18l_tagged_flux");
  cout << " ***** h18l_tagged_flux = " << h18l_tagged_flux << endl;
  h18l_tagged_flux->SetMarkerColor(kMagenta);
  h18l_tagged_flux->SetMarkerSize(1.5);
  h18l_tagged_flux->SetMarkerStyle(kOpenTriangleUp);
  h18_tagged_flux->SetMinimum(0.);
  h18_tagged_flux->Draw("e");
  h18l_tagged_flux->Draw("esame");
  h17_tagged_flux->Draw("esame");
  // h16_tagged_flux->Draw("esame");
  ctot_tagged_flux->Modified();
  // ltot_tagged_flux->AddEntry(h16_tagged_flux, "2016", "p");
  ltot_tagged_flux->AddEntry(h17_tagged_flux, "2017", "p");
  ltot_tagged_flux->AddEntry(h18_tagged_flux, "2018 S", "p");
  ltot_tagged_flux->AddEntry(h18l_tagged_flux, "2018 F", "p");
  ltot_tagged_flux->Draw();
  ctot_tagged_flux->Print("output/fig_phi2pi/ctot_tagged_flux.eps", "eps");
  ctot_tagged_flux->Print("output/fig_phi2pi/ctot_tagged_flux.png", "png");
  ctot_tagged_flux->Print("output/fig_phi2pi/ctot_tagged_flux.root", "root");

  // *********** total Thrown Beam energy
  TCanvas *ctot_beame_tru = new TCanvas("ctot_beame_tru", "ctot_beame_tru", 900, 600);
  TLegend *ltot_beame_tru = new TLegend(0.84, 0.80, 0.98, 0.98);
  ltot_beame_tru->SetTextSize(0.06);
  ltot_beame_tru->SetBorderSize(0);
  ctot_beame_tru->cd();
  ctot_beame_tru->SetGrid();
  TH1F *h16_beame_tru = (TH1F *)outputfig->Get("h16_beame_tru");
  cout << " ***** h16_beame_tru = " << h16_beame_tru << endl;
  // h16_beame_tru->SetMarkerColor(kBlue);
  h16_beame_tru->SetMarkerSize(1.5);
  h16_beame_tru->SetMarkerStyle(kFullSquare);
  TH1F *h17_beame_tru = (TH1F *)outputfig->Get("h17_beame_tru");
  cout << " ***** h17_beame_tru = " << h17_beame_tru << endl;
  h17_beame_tru->SetMarkerSize(1.5);
  h17_beame_tru->SetMarkerColor(kBlue);
  h17_beame_tru->SetMarkerStyle(kFullCircle);
  TH1F *h18_beame_tru = (TH1F *)outputfig->Get("h18_beame_tru");
  cout << " ***** h18_beame_tru = " << h18_beame_tru << endl;
  h18_beame_tru->SetMarkerColor(kRed);
  h18_beame_tru->SetMarkerSize(1.5);
  h18_beame_tru->SetMarkerStyle(kFullTriangleUp);
  TH1F *h18l_beame_tru = (TH1F *)outputfig->Get("h18l_beame_tru");
  cout << " ***** h18l_beame_tru = " << h18l_beame_tru << endl;
  h18l_beame_tru->SetMarkerColor(kMagenta);
  h18l_beame_tru->SetMarkerSize(1.5);
  h18l_beame_tru->SetMarkerStyle(kOpenTriangleUp);
  h17_beame_tru->SetMinimum(0.);
  h18_beame_tru->Draw("e");
  h18l_beame_tru->Draw("esame");
  h17_beame_tru->Draw("esame");
  // h16_beame_tru->Draw("esame");
  ctot_beame_tru->Modified();
  // ltot_beame_tru->AddEntry(h16_beame_tru, "2016", "p");
  ltot_beame_tru->AddEntry(h17_beame_tru, "2017", "p");
  ltot_beame_tru->AddEntry(h18_beame_tru, "2018 S", "p");
  ltot_beame_tru->AddEntry(h18l_beame_tru, "2018 F", "p");
  ltot_beame_tru->Draw();
  ctot_beame_tru->Print("output/fig_phi2pi/ctot_beame_tru.eps", "eps");
  ctot_beame_tru->Print("output/fig_phi2pi/ctot_beame_tru.png", "png");
  ctot_beame_tru->Print("output/fig_phi2pi/ctot_beame_tru.root", "root");

  // *********** total Yields MC
  TMultiGraph *mgphie_mc = new TMultiGraph();
  TCanvas *cmgphie_mc = new TCanvas("cmgphie_mc", "cmgphie_mc", 900, 600);
  TLegend *lmgphie_mc = new TLegend(0.84, 0.80, 0.98, 0.98);
  lmgphie_mc->SetTextSize(0.06);
  lmgphie_mc->SetBorderSize(0);
  cmgphie_mc->cd();
  cmgphie_mc->SetGrid();
  TGraphErrors *h16_gphie_mc = (TGraphErrors *)outputfig->Get("h16_gphie_mc");
  cout << " ***** h16_gphie_mc = " << h16_gphie_mc << endl;
  // h16_gphie_mc->SetMarkerColor(kBlue);
  h16_gphie_mc->SetMarkerSize(1.5);
  h16_gphie_mc->SetMarkerStyle(kFullSquare);
  // mgphie_mc->Add(h16_gphie_mc);
  TGraphErrors *h17_gphie_mc = (TGraphErrors *)outputfig->Get("h17_gphie_mc");
  cout << " ***** h17_gphie_mc = " << h17_gphie_mc << endl;
  h17_gphie_mc->SetMarkerSize(1.5);
  h17_gphie_mc->SetMarkerColor(kBlue);
  h17_gphie_mc->SetMarkerStyle(kFullCircle);
  mgphie_mc->Add(h17_gphie_mc);
  TGraphErrors *h18_gphie_mc = (TGraphErrors *)outputfig->Get("h18_gphie_mc");
  cout << " ***** h18_gphie_mc = " << h18_gphie_mc << endl;
  h18_gphie_mc->SetMarkerColor(kRed);
  h18_gphie_mc->SetMarkerSize(1.5);
  h18_gphie_mc->SetMarkerStyle(kFullTriangleUp);
  mgphie_mc->Add(h18_gphie_mc);
  TGraphErrors *h18l_gphie_mc = (TGraphErrors *)outputfig->Get("h18l_gphie_mc");
  cout << " ***** h18l_gphie_mc = " << h18l_gphie_mc << endl;
  h18l_gphie_mc->SetMarkerColor(kMagenta);
  h18l_gphie_mc->SetMarkerSize(1.5);
  h18l_gphie_mc->SetMarkerStyle(kOpenTriangleUp);
  mgphie_mc->Add(h18l_gphie_mc);
  mgphie_mc->SetTitle("; E_{#gamma} (GeV); N_{#phi #pi^{+} #pi^{-}}");
  mgphie_mc->Draw("AP");
  mgphie_mc->SetMinimum(0.);
  cmgphie_mc->Modified();
  // lmgphie_mc->AddEntry(h16_gphie_mc, "2016", "p");
  lmgphie_mc->AddEntry(h17_gphie_mc, "2017", "p");
  lmgphie_mc->AddEntry(h18_gphie_mc, "2018 S", "p");
  lmgphie_mc->AddEntry(h18l_gphie_mc, "2018 F", "p");
  lmgphie_mc->Draw();
  cmgphie_mc->Print("output/fig_phi2pi/cmgphie_mc.eps", "eps");
  cmgphie_mc->Print("output/fig_phi2pi/cmgphie_mc.png", "png");
  cmgphie_mc->Print("output/fig_phi2pi/cmgphie_mc.root", "root");

  // *********** total Yields Data
  TMultiGraph *mgphie = new TMultiGraph();
  TCanvas *cmgphie = new TCanvas("cmgphie", "cmgphie", 900, 600);
  TLegend *lmgphie = new TLegend(0.84, 0.80, 0.98, 0.98);
  lmgphie->SetTextSize(0.06);
  lmgphie->SetBorderSize(0);
  cmgphie->cd();
  cmgphie->SetGrid();
  TGraphErrors *h16_gphie = (TGraphErrors *)outputfig->Get("h16_gphie");
  cout << " ***** h16_gphie = " << h16_gphie << endl;
  // h16_gphie->SetMarkerColor(kBlue);
  h16_gphie->SetMarkerSize(1.5);
  h16_gphie->SetMarkerStyle(kFullSquare);
  // mgphie->Add(h16_gphie);
  TGraphErrors *h17_gphie = (TGraphErrors *)outputfig->Get("h17_gphie");
  cout << " ***** h17_gphie = " << h17_gphie << endl;
  h17_gphie->SetMarkerSize(1.5);
  h17_gphie->SetMarkerColor(kBlue);
  h17_gphie->SetMarkerStyle(kFullCircle);
  mgphie->Add(h17_gphie);
  TGraphErrors *h18_gphie = (TGraphErrors *)outputfig->Get("h18_gphie");
  cout << " ***** h18_gphie = " << h18_gphie << endl;
  h18_gphie->SetMarkerColor(kRed);
  h18_gphie->SetMarkerSize(1.5);
  h18_gphie->SetMarkerStyle(kFullTriangleUp);
  mgphie->Add(h18_gphie);
  TGraphErrors *h18l_gphie = (TGraphErrors *)outputfig->Get("h18l_gphie");
  cout << " ***** h18l_gphie = " << h18l_gphie << endl;
  h18l_gphie->SetMarkerColor(kMagenta);
  h18l_gphie->SetMarkerSize(1.5);
  h18l_gphie->SetMarkerStyle(kOpenTriangleUp);
  mgphie->Add(h18l_gphie);
  mgphie->SetTitle("; E_{#gamma} (GeV); N_{#phi #pi^{+} #pi^{-}}");
  mgphie->Draw("AP");
  mgphie->SetMinimum(0.);
  cmgphie->Modified();
  // lmgphie->AddEntry(h16_gphie, "2016", "p");
  lmgphie->AddEntry(h17_gphie, "2017", "p");
  lmgphie->AddEntry(h18_gphie, "2018 S", "p");
  lmgphie->AddEntry(h18l_gphie, "2018 F", "p");
  lmgphie->Draw();
  cmgphie->Print("output/fig_phi2pi/cmgphie.eps", "eps");
  cmgphie->Print("output/fig_phi2pi/cmgphie.png", "png");
  cmgphie->Print("output/fig_phi2pi/cmgphie.root", "root");

  // *********** total Efficiency
  TCanvas *cmgeeff = new TCanvas("cmgeeff", "cmgeeff", 800, 800);
  cmgeeff->cd();
  cmgeeff->SetGrid();

  float aeeff = 0.3;
  float beeff = 0.02;

  TPad *padeeff1 = new TPad("padeeff1", "padeeff1", 0, aeeff - beeff, 1, 1);
  padeeff1->SetBottomMargin(beeff);
  padeeff1->SetGridx();
  padeeff1->SetGridy();
  cmgeeff->cd();
  padeeff1->Draw();

  TPad *padeeff2 = new TPad("padeeff2", "padeeff2", 0, 0, 1, aeeff * (1 - beeff));
  padeeff2->SetTopMargin(0);
  padeeff2->SetBottomMargin(0.4);
  padeeff2->SetFillColor(0);
  padeeff2->SetFillStyle(0);
  padeeff2->SetGridy();
  cmgeeff->cd();
  padeeff2->Draw();

  padeeff1->cd();
  // TMultiGraph *mgeeff = new TMultiGraph();
  TLegend *lmgeeff = new TLegend(0.85, 0.78, 0.98, 0.98);
  lmgeeff->SetTextSize(0.06);
  lmgeeff->SetBorderSize(0);
  TGraphErrors *h16_gphieeff = (TGraphErrors *)outputfig->Get("h16_gphieeff");
  cout << " ***** h16_gphieeff = " << h16_gphieeff << endl;
  // h16_gphieeff->SetMarkerColor(kBlue);
  // h16_gphieeff->SetMarkerSize(1.5);
  h16_gphieeff->SetMarkerStyle(kFullSquare);
  h16_gphieeff->SetTitle("; E_{#gamma} (GeV); #varepsilon (%)");
  h16_gphieeff->SetMinimum(0.);
  // mgeeff->Add(h16_gphieeff);
  TGraphErrors *h17_gphieeff = (TGraphErrors *)outputfig->Get("h17_gphieeff");
  cout << " ***** h17_gphieeff = " << h17_gphieeff << endl;
  h17_gphieeff->SetMarkerColor(kBlue);
  h17_gphieeff->SetMarkerSize(1.5);
  h17_gphieeff->SetMarkerStyle(kFullCircle);
  h17_gphieeff->SetTitle("; E_{#gamma} (GeV); #varepsilon (%)");
  h17_gphieeff->SetMinimum(0.);
  // mgeeff->Add(h17_gphieeff);
  TGraphErrors *h18_gphieeff = (TGraphErrors *)outputfig->Get("h18_gphieeff");
  cout << " ***** h18_gphieeff = " << h18_gphieeff << endl;
  h18_gphieeff->SetMarkerColor(kRed);
  h18_gphieeff->SetMarkerSize(1.5);
  h18_gphieeff->SetMarkerStyle(kFullTriangleUp);
  // mgeeff->Add(h18_gphieeff);
  TGraphErrors *h18l_gphieeff = (TGraphErrors *)outputfig->Get("h18l_gphieeff");
  cout << " ***** h18l_gphieeff = " << h18l_gphieeff << endl;
  h18l_gphieeff->SetMarkerColor(kMagenta);
  h18l_gphieeff->SetMarkerSize(1.5);
  h18l_gphieeff->SetMarkerStyle(kOpenTriangleUp);
  // mgeeff->Add(h18l_gphieeff);
  // mgeeff->SetTitle("; E_{#gamma} (GeV); #varepsilon (%)");
  // mgeeff->Draw("AP");
  // mgeeff->SetMinimum(0.);
  // h16_gphieeff->Draw("AP");
  h18l_gphieeff->Draw("AP");
  h17_gphieeff->Draw("PSAME");
  h18_gphieeff->Draw("PSAME");

  gPad->Modified();
  gPad->Update();
  cmgeeff->Modified();
  // lmgeeff->AddEntry(h16_gphieeff, "2016", "p");
  lmgeeff->AddEntry(h17_gphieeff, "2017", "p");
  lmgeeff->AddEntry(h18_gphieeff, "2018 S", "p");
  lmgeeff->AddEntry(h18l_gphieeff, "2018 F", "p");
  lmgeeff->Draw();

  padeeff2->cd(); // pad2 becomes the current pad
  // TGraphErrors *gr_efferate_17 = new TGraphErrors();
  TGraphErrors *gr_efferate_18 = new TGraphErrors();
  TGraphErrors *gr_efferate_18l = new TGraphErrors();
  // gr_efferate_17->SetMarkerColor(kBlue);
  // gr_efferate_17->SetMarkerStyle(2);
  // gr_efferate_17->SetMarkerSize(2);
  gr_efferate_18->SetMarkerColor(kRed);
  gr_efferate_18->SetMarkerStyle(2);
  gr_efferate_18->SetMarkerSize(2);
  gr_efferate_18l->SetMarkerColor(kMagenta);
  gr_efferate_18l->SetMarkerStyle(2);
  gr_efferate_18l->SetMarkerSize(2);

  // Eg = 0;
  for (int i = 1; i <= ne; ++i)
  {
    double Eg1 = h2phie_mc->GetXaxis()->GetBinLowEdge(i);
    double Eg2 = h2phie_mc->GetXaxis()->GetBinLowEdge(i) + Eg_step;
    double Eg = h2phie_mc->GetXaxis()->GetBinCenter(i);
    // if (i == 1)
    //   Eg = Eg1;
    // else
    //   Eg = Eg2;

    // double efferate_17 = (h16_gphieeff->GetY()[i - 1] - h17_gphieeff->GetY()[i - 1]) / (h16_gphieeff->GetY()[i - 1]);
    double efferate_18 = (h17_gphieeff->GetY()[i - 1] - h18_gphieeff->GetY()[i - 1]) / (h17_gphieeff->GetY()[i - 1]);
    double efferate_18l = (h17_gphieeff->GetY()[i - 1] - h18l_gphieeff->GetY()[i - 1]) / (h17_gphieeff->GetY()[i - 1]);
    // cout << efferate_18 = " << efferate_18 << " | efferate_18l = " << efferate_18l << endl;
    // gr_efferate_17->SetPoint(i-1, Eg, efferate_17);//i-1, h2phie_mc->GetXaxis()->GetBinCenter(i), efferate_17
    // gr_efferate_17->SetPointError(i-1, 0, 0);
    gr_efferate_18->SetPoint(i - 1, Eg, efferate_18);
    gr_efferate_18->SetPointError(i - 1, 0, 0);
    gr_efferate_18l->SetPoint(i - 1, Eg, efferate_18l);
    gr_efferate_18l->SetPointError(i - 1, 0, 0);
    }
    
    // gr_efferate_17->Draw("AP");
    gr_efferate_18->Draw("AP");
    gr_efferate_18l->Draw("PSAME");
    TLine leffe;
    leffe.SetLineStyle(kDashed);
    leffe.DrawLine(Eg_min, 0, Eg_max, 0);

    // Y axis ratio plot settings
    gr_efferate_18->GetHistogram()->SetStats(0);
    gr_efferate_18->SetTitle(";E_{#gamma} (GeV);Ratio");
    // gr_efferate_18->GetYaxis()->SetNdivisions(506);
    gr_efferate_18->GetYaxis()->SetTitleSize(35);
    gr_efferate_18->GetYaxis()->SetTitleFont(43);
    gr_efferate_18->GetYaxis()->SetTitleOffset(1.55);
    gr_efferate_18->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efferate_18->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_efferate_18->GetXaxis()->SetTitleSize(38);
    gr_efferate_18->GetXaxis()->SetTitleFont(43);
    gr_efferate_18->GetXaxis()->SetTitleOffset(3.);
    gr_efferate_18->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efferate_18->GetXaxis()->SetLabelSize(38);
    gr_efferate_18->GetHistogram()->SetMaximum(1.0);   // along          
    gr_efferate_18->GetHistogram()->SetMinimum(-1.0);  //   Y 
    gPad->Update();
    gPad->Modified();
    cmgeeff->Print("output/fig_phi2pi/cmgeeff.eps", "eps");
    cmgeeff->Print("output/fig_phi2pi/cmgeeff.png", "png");
    cmgeeff->Print("output/fig_phi2pi/cmgeeff.root", "root");

    // ************ total Cross-sections

    TCanvas *cmgexsec = new TCanvas("cmgexsec", "cmgexsec", 800, 800);
    // cmgexsec->Divide(2);
    cmgexsec->cd();
    cmgexsec->SetGrid();

    float aexsec = 0.3;
    float bexsec = 0.02;

    TPad *padexsec1 = new TPad("padexsec1", "padexsec1", 0, aexsec - bexsec, 1, 1);
    padexsec1->SetBottomMargin(bexsec);
    padexsec1->SetGridx();
    padexsec1->SetGridy();
    cmgexsec->cd();
    padexsec1->Draw();

    TPad *padexsec2 = new TPad("padexsec2", "padexsec2", 0, 0, 1, aexsec * (1 - bexsec));
    padexsec2->SetTopMargin(0);
    padexsec2->SetBottomMargin(0.4);
    padexsec2->SetFillColor(0);
    padexsec2->SetFillStyle(0);
    padexsec2->SetGridy();
    cmgexsec->cd();
    padexsec2->Draw();

    padexsec1->cd();
    // TMultiGraph *mgexsec = new TMultiGraph();
    TLegend *lmgexsec = new TLegend(0.8, 0.78, 0.93, 0.98);
    lmgexsec->SetTextSize(0.06);
    lmgexsec->SetBorderSize(0);
    TGraphErrors *h16_gphiexsec = (TGraphErrors *)outputfig->Get("h16_gphiexsec");
    cout << " ***** h16_gphiexsec = " << h16_gphiexsec << endl;
    h16_gphiexsec->SetMarkerStyle(kFullSquare);
    h16_gphiexsec->SetMarkerSize(1.5);
    h16_gphiexsec->SetTitle("; E_{#gamma} (GeV); #sigma (nb)");
    h16_gphiexsec->SetMinimum(0.);
    // h16_gphiexsec->SetMarkerColor(kBlue);
    // mgexsec->Add(h16_gphiexsec);
    TGraphErrors *h17_gphiexsec = (TGraphErrors *)outputfig->Get("h17_gphiexsec");
    cout << " ***** h17_gphiexsec = " << h17_gphiexsec << endl;
    h17_gphiexsec->SetMarkerStyle(kFullCircle);
    h17_gphiexsec->SetMarkerSize(1.5);
    h17_gphiexsec->SetMarkerColor(kBlue);
    h17_gphiexsec->SetLineColor(kBlue);
    h17_gphiexsec->SetTitle("; E_{#gamma} (GeV); #sigma (nb)");
    h17_gphiexsec->SetMinimum(0.);
    // mgexsec->Add(h17_gphiexsec);
    TGraphErrors *h18_gphiexsec = (TGraphErrors *)outputfig->Get("h18_gphiexsec");
    cout << " ***** h18_gphiexsec = " << h18_gphiexsec << endl;
    h18_gphiexsec->SetMarkerStyle(kFullTriangleUp);//SetMarkerColor(kRed);
    h18_gphiexsec->SetMarkerSize(1.5);
    h18_gphiexsec->SetMarkerColor(kRed);
    h18_gphiexsec->SetLineColor(kRed);
    // mgexsec->Add(h18_gphiexsec);
    TGraphErrors *h18l_gphiexsec = (TGraphErrors *)outputfig->Get("h18l_gphiexsec");
    cout << " ***** h18l_gphiexsec = " << h18l_gphiexsec << endl;
    h18l_gphiexsec->SetMarkerStyle(kOpenTriangleUp);//SetMarkerColor(kMagenta);
    h18l_gphiexsec->SetMarkerSize(1.5);
    h18l_gphiexsec->SetMarkerColor(kMagenta);
    h18l_gphiexsec->SetLineColor(kMagenta);
    // mgexsec->Add(h18l_gphiexsec);
    // mgexsec->SetTitle(";E_{#gamma} (GeV);#sigma (nb)");
    // mgexsec->Draw("AP");
    // mgexsec->SetMinimum(0.);
    // h16_gphiexsec->Draw("AP");
    h17_gphiexsec->Draw("AP");
    h18_gphiexsec->Draw("PSAME");
    h18l_gphiexsec->Draw("PSAME");
    gPad->Modified();
    gPad->Update();
    // lmgexsec->AddEntry(h16_gphiexsec, "2016", "p");
    lmgexsec->AddEntry(h17_gphiexsec, "2017", "p");
    lmgexsec->AddEntry(h18_gphiexsec, "2018 S", "p");
    lmgexsec->AddEntry(h18l_gphiexsec, "2018 F", "p");
    lmgexsec->Draw();
 
    padexsec2->cd(); // pad2 becomes the current pad
    // TGraphErrors *gr_exsecrate_17 = new TGraphErrors();
    TGraphErrors *gr_exsecrate_18 = new TGraphErrors();
    TGraphErrors *gr_exsecrate_18l = new TGraphErrors();
    // gr_exsecrate_17->SetMarkerColor(kBlue);
    // gr_exsecrate_17->SetMarkerStyle(2);//20
    // gr_exsecrate_17->SetMarkerSize(2);
    gr_exsecrate_18->SetMarkerColor(kRed);
    gr_exsecrate_18->SetMarkerStyle(2);
    gr_exsecrate_18->SetMarkerSize(2);
    gr_exsecrate_18l->SetMarkerColor(kMagenta);
    gr_exsecrate_18l->SetMarkerStyle(2);
    gr_exsecrate_18l->SetMarkerSize(2);

    // Eg = 0;
    for (int i = 1; i <= ne; ++i)
    {
      double Eg1 = h2phie_mc->GetXaxis()->GetBinLowEdge(i);
      double Eg2 = h2phie_mc->GetXaxis()->GetBinLowEdge(i) + Eg_step;
      double Eg = h2phie_mc->GetXaxis()->GetBinCenter(i);
      // if (i == 1)
      //   Eg = Eg1;
      // else
      //   Eg = Eg2;

      // double exsecrate_17 = (h16_gphiexsec->GetY()[i-1] - h17_gphiexsec->GetY()[i-1])/(h16_gphiexsec->GetY()[i-1]);
      double exsecrate_18 = (h17_gphiexsec->GetY()[i-1] - h18_gphiexsec->GetY()[i-1])/(h17_gphiexsec->GetY()[i-1]);
      double exsecrate_18l = (h17_gphiexsec->GetY()[i-1] - h18l_gphiexsec->GetY()[i-1])/(h17_gphiexsec->GetY()[i-1]);
      // cout<<"exsecrate_17 = "<<exsecrate_17<<" | exsecrate_18 = "<<exsecrate_18<<" | exsecrate_18l = "<<exsecrate_18l<<endl;
      // gr_exsecrate_17->SetPoint(i-1, Eg, exsecrate_17);//h2phie->GetXaxis()->GetBinCenter(i)
      // gr_exsecrate_17->SetPointError(i-1, 0, 0);
      gr_exsecrate_18->SetPoint(i-1, Eg, exsecrate_18);
      gr_exsecrate_18->SetPointError(i-1, 0, 0);
      gr_exsecrate_18l->SetPoint(i-1, Eg, exsecrate_18l);
      gr_exsecrate_18l->SetPointError(i-1, 0, 0);
    }

    // gr_exsecrate_17->Draw("AP");
    gr_exsecrate_18->Draw("AP");
    gr_exsecrate_18l->Draw("PSAME");
    TLine lexsec;
    lexsec.SetLineStyle(kDashed);
    lexsec.DrawLine(Eg_min, 0, Eg_max, 0);

    // Y axis ratio plot settings
    gr_exsecrate_18->GetHistogram()->SetStats(0);
    gr_exsecrate_18->SetTitle(";E_{#gamma} (GeV);Ratio");
    // gr_exsecrate_18->GetYaxis()->SetNdivisions(506);
    gr_exsecrate_18->GetYaxis()->SetTitleSize(35);
    gr_exsecrate_18->GetYaxis()->SetTitleFont(43);
    gr_exsecrate_18->GetYaxis()->SetTitleOffset(1.55);
    gr_exsecrate_18->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_exsecrate_18->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_exsecrate_18->GetXaxis()->SetTitleSize(38);
    gr_exsecrate_18->GetXaxis()->SetTitleFont(43);
    gr_exsecrate_18->GetXaxis()->SetTitleOffset(3.);
    gr_exsecrate_18->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_exsecrate_18->GetXaxis()->SetLabelSize(38);

    gr_exsecrate_18->GetHistogram()->SetMaximum(0.5);   // along          
    gr_exsecrate_18->GetHistogram()->SetMinimum(-0.5);  //   Y  
    gPad->Update();
    gPad->Modified();
    cmgexsec->Print("output/fig_phi2pi/cmgexsec.eps", "eps");
    cmgexsec->Print("output/fig_phi2pi/cmgexsec.png", "png");
    cmgexsec->Print("output/fig_phi2pi/cmgexsec.root", "root");

    // +++++++++++++++++++ PDG cross section average
    TCanvas *cgexsec_avg = new TCanvas("cgexsec_avg", "cgexsec_avg", 900, 600);
    cgexsec_avg->cd();
    cgexsec_avg->SetGrid();
    TLegend *lgexsec_avg = new TLegend(0.67, 0.78, 0.93, 0.98);
    lgexsec_avg->SetTextSize(0.06);
    lgexsec_avg->SetBorderSize(0);
    TGraphErrors *gphiexsec_avg = new TGraphErrors();
    gphiexsec_avg->SetMarkerStyle(20);
    gphiexsec_avg->SetMarkerSize(1.5);
    gphiexsec_avg->SetMarkerColor(28);
    gphiexsec_avg->SetLineColor(28);
    gphiexsec_avg->SetMinimum(0.0);
    gphiexsec_avg->SetTitle("; E_{#gamma} (GeV); #sigma (nb)");
    // TGraphErrors *gphiexsec_17corr = new TGraphErrors();
    // gphiexsec_17corr->SetMarkerStyle(kOpenCircle);
    // gphiexsec_17corr->SetMarkerSize(1.5);
    // gphiexsec_17corr->SetMarkerColor(kBlue);
    // gphiexsec_17corr->SetLineColor(kBlue);
    // gphiexsec_17corr->SetMinimum(0.0);
    // gphiexsec_17corr->SetTitle("; E_{#gamma} (GeV); #sigma (nb)");

    // vector<double> xi = {1.0,0.8,1.1}, dxi = {0.1,0.2,0.05};
    // pair<double, double> x = pdgavg(xi, dxi);
    // Eg = 0;
    vector<double> xsec(3), dxsec(3);
    pair<double, double> xsec_pair;

    for (int i = 1; i <= ne; ++i)
    {
      xsec = {h17_gphiexsec->GetY()[i - 1], h18_gphiexsec->GetY()[i - 1], h18l_gphiexsec->GetY()[i - 1]};
      dxsec = {h17_gphiexsec->GetErrorY(i - 1), h18_gphiexsec->GetErrorY(i - 1), h18l_gphiexsec->GetErrorY(i - 1)};

      double Eg1 = h2phie_mc->GetXaxis()->GetBinLowEdge(i);
      double Eg2 = h2phie_mc->GetXaxis()->GetBinLowEdge(i) + Eg_step;
      double Eg = h2phie_mc->GetXaxis()->GetBinCenter(i);
      // if (i == 1)
      //   Eg = Eg1;
      // else
      //   Eg = Eg2;

      xsec_pair = pdgavg(xsec, dxsec);  
      
      gphiexsec_avg->SetPoint(i-1, Eg, xsec_pair.first);
      gphiexsec_avg->SetPointError(i-1, 0, xsec_pair.second);

      // double w_17 = 0.76; //0.76
      
      // gphiexsec_17corr->SetPoint(i-1, Eg, w_17*h17_gphiexsec->GetY()[i - 1]);
      // gphiexsec_17corr->SetPointError(i-1, 0, w_17*h17_gphiexsec->GetErrorY(i - 1));

      // test << std::setprecision(2) << std::fixed;
      // test << "i-1 = " << i - 1 << " | xsec[0] = " << xsec[0] << " | dxsec[0] = " << dxsec[0] << " | xsec[1] = " << xsec[1] << " | dxsec[1] = " << dxsec[1] <<  " | xsec[2] = " << xsec[2] << " | dxsec[2] = " << dxsec[2] << " | xsec_pair.first = " << xsec_pair.first << " | xsec_pair.second = " << xsec_pair.second <<  " | w_17 = " << w_17 << endl;
    }

    // h17_gphiexsec->Draw("AP");
    gphiexsec_avg->Draw("AP");
    // gphiexsec_17corr->Draw("PSAME");
    // lgexsec_avg->AddEntry(h17_gphiexsec, "2017", "p");
    lgexsec_avg->AddEntry(gphiexsec_avg, "#bar{#sigma}_{(2017, 2018 S, 2018 F)}", "p");
    // lgexsec_avg->AddEntry(gphiexsec_17corr, "2017 corrected", "p");
    lgexsec_avg->Draw();

    cgexsec_avg->Print("output/fig_phi2pi/cgexsec_avg.eps", "eps");
    cgexsec_avg->Print("output/fig_phi2pi/cgexsec_avg.png", "png");
    cgexsec_avg->Print("output/fig_phi2pi/cgexsec_avg.root", "root");

    cphie->Print(Form("output/fig_phi2pi/c2_phie_%s.root", name.Data()), "root");
    cphie->Print(Form("output/fig_phi2pi/c2_phie_%s.eps", name.Data()), "eps");
    cphie->Print(Form("output/fig_phi2pi/c2_phie_%s.png", name.Data()), "png");
    cphie1->Print(Form("output/fig_phi2pi/c_phie1_%s.root", name.Data()), "root");
    cphie1->Print(Form("output/fig_phi2pi/c_phie1_%s.eps", name.Data()), "eps");
    cphie1->Print(Form("output/fig_phi2pi/c_phie1_%s.png", name.Data()), "png");
    cgphie->Print(Form("output/fig_phi2pi/c_gphie_%s.root", name.Data()), "root");
    cgphie->Print(Form("output/fig_phi2pi/c_gphie_%s.eps", name.Data()), "eps");
    cgphie->Print(Form("output/fig_phi2pi/c_gphie_%s.png", name.Data()), "png");
    cphie_mc->Print(Form("output/fig_phi2pi/c2_phie_mc_%s.root", name.Data()), "root");
    cphie_mc->Print(Form("output/fig_phi2pi/c2_phie_mc_%s.eps", name.Data()), "eps");
    cphie_mc->Print(Form("output/fig_phi2pi/c2_phie_mc_%s.png", name.Data()),"png");
    cphie1_mc->Print(Form("output/fig_phi2pi/c_phie1_mc_%s.root", name.Data()),"root");
    cphie1_mc->Print(Form("output/fig_phi2pi/c_phie1_mc_%s.eps", name.Data()),"eps");
    cphie1_mc->Print(Form("output/fig_phi2pi/c_phie1_mc_%s.png", name.Data()),"png");
    cgphie_mc->Print(Form("output/fig_phi2pi/c_gphie_mc_%s.root", name.Data()),"root");
    cgphie_mc->Print(Form("output/fig_phi2pi/c_gphie_mc_%s.eps", name.Data()),"eps");
    cgphie_mc->Print(Form("output/fig_phi2pi/c_gphie_mc_%s.png", name.Data()),"png");
    cgphieeff->Print(Form("output/fig_phi2pi/c_gphieeff_%s.root", name.Data()),"root");
    cgphieeff->Print(Form("output/fig_phi2pi/c_gphieeff_%s.eps", name.Data()),"eps");
    cgphieeff->Print(Form("output/fig_phi2pi/c_gphieeff_%s.png", name.Data()),"png");
    cgphiexsec->Print(Form("output/fig_phi2pi/c_gphiexsec_%s.root", name.Data()), "root");
    cgphiexsec->Print(Form("output/fig_phi2pi/c_gphiexsec_%s.eps", name.Data()), "eps");
    cgphiexsec->Print(Form("output/fig_phi2pi/c_gphiexsec_%s.png", name.Data()), "png");
    c_beame_tru->Print(Form("output/fig_phi2pi/c_beame_tru_%s.root", name.Data()),"root");
    c_beame_tru->Print(Form("output/fig_phi2pi/c_beame_tru_%s.eps", name.Data()),"eps");
    c_beame_tru->Print(Form("output/fig_phi2pi/c_beame_tru_%s.png", name.Data()),"png");
    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.png", name.Data()), "png");
    
    outputfig->Print();

    fprintf(table_xsec_phi2pi, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_xsec_phi2pi);
    gSystem->Exec(Form("pdflatex table_xsec_phi2pi_%s.tex", name.Data()));

    fprintf(table_xsec_phi2pi_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_xsec_phi2pi_sys);
    gSystem->Exec(Form("pdflatex table_xsec_phi2pi_sys_%s.tex", name.Data()));

    // table_phi << "\\hline" << endl;
    // table_phi << "\\end{tabularx}" << endl;
    // table_phi << "\\end{center}" << endl;
    // table_phi << "\\end{minipage}" << endl;
    // table_phi << "\\end{table}" << endl;
    // table_phi << "\\end{document}" << endl;
    // table_phi.close();
    // gSystem->Exec("pdflatex table_phi.tex");
*/
  // ======================================== Phi vs. -t ===============================================

  //root -l 'xsec_phi2pi.C+("data_17",100,1,10)'
  // latex: \\SI{70}{\\percent} -> \\%

  FILE *table_xsec_phi2pi = fopen(Form("table_xsec_phi2pi_t_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{siunitx} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\small\n \\caption{Total cross-sections in momentum transfer bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|}\n \\hline\n -t $(GeV/c)^{2}$ & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\varepsilon$ [\\%] & $\\sigma$ (nb) \\\\ \n \\hline\n");

  FILE *table_xsec_phi2pi_sys = fopen(Form("table_xsec_phi2pi_sys_t_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi_sys,"\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{makecell} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors in momentum transfer bins}\n \\begin{tabular}{|c|c|c|c|c|c|c|c|}\n \\hline\n -t $(GeV/c)^{2}$ & \\thead{Bkg deg\\\\(3,4,5)}& \\thead{Fit range\\\\(0.99; 1,15, 1.2, 1.25)} & \\thead{bining\\\\(90, 100, 110)} & \\thead{Accidentals\\\\(2, 4, 8)} & \\thead{$\\Delta T_{RF-TOF}$ proton\\\\($\\pm$0.2, $\\pm$0.3, $\\pm$0.4)} & \\thead{$\\chi^{2}$ KinFit\\\\(45, 55, 65)} & \\thead{$MM^{2}$\\\\($\\pm$0.025, $\\pm$0.035, $\\pm$0.045)} \\\\ \n \\hline\n");
  
  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(100);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Momentum Transfer
  TCanvas *c_h_t_Thrown = new TCanvas("c_h_t_Thrown", "c_h_t_Thrown", 900, 600);
  c_h_t_Thrown->cd();
  TH1F *h_t_Thrown = (TH1F *)ftru->Get("h_t_Thrown");
  cout << "h_t_Thrown = " << h_t_Thrown << endl;
  h_t_Thrown->Rebin(10);
  h_t_Thrown->SetMarkerStyle(20);
  h_t_Thrown->SetMarkerSize(1.5);
  h_t_Thrown->Draw("e");
  h_t_Thrown->Write(Form("h%s_t_tru", name.Data()), TObject::kWriteDelete);

  // ++++++++++++++++++++++++++++ Data  ++++++++++++++++++++++++++++ 
  TCanvas *cphit = new TCanvas("cphit", "cphit", 900, 600);
  cphit->cd();
  TH2D *h2phit = new TH2D("h2phit", ";-t (GeV/c)^{2};m_{K^{+}K^{-}} (GeV/c^{2});Counts", nt, 0, 4, n2k, mkk_min, mkk_max);
  tdata->Project("h2phit", "kpkm_mf:-t_kin", "w8*(kpkm_uni)");// && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(kpkm_uni && abs(rf_dt)<1.5*4.008), w4*(kpkm_uni && abs(rf_dt)<2.5*4.008), w8*(kpkm_uni)

  h2phit->Draw("colz");

  TCanvas *cphit1 = new TCanvas("cphit1", "cphit1", 900, 1200);
  cphit1->Divide(2, 5);
  TCanvas *cgphit = new TCanvas("cgphit", "cgphit", 900, 600);
  TGraphErrors *gphit = new TGraphErrors();
  gphit->SetMarkerStyle(20);
  gphit->SetMarkerSize(1.5);
  gphit->SetMarkerColor(1);
  gphit->SetMinimum(0.0);
  gphit->SetTitle(";-t (GeV/c)^{2}; N_{#phi #pi^{+} #pi^{-}}");

  //  ++++++++++++++++++++++++++++ Monte Carlo  ++++++++++++++++++++++++++++ 
  TCanvas *cphit_mc = new TCanvas("cphit_mc", "cphit_mc", 900, 600);
  cphit_mc->cd();
  TH2D *h2phit_mc = new TH2D("h2phit_mc", ";-t (GeV/c)^{2};m_{K^{+}K^{-}} (GeV/c^{2});Counts", nt, 0, 4, n2k, mkk_min, mkk_max);//0.98, 1.2
  tmc->Project("h2phit_mc", "kpkm_mf:-t_kin", "w8*(kpkm_uni)");
  h2phit_mc->Draw("colz");

  TCanvas *cphit1_mc = new TCanvas("cphit1_mc", "cphit1_mc", 900, 1200);
  cphit1_mc->Divide(2, 5);
  TCanvas *cgphit_mc = new TCanvas("cgphit_mc", "cgphit_mc", 900, 600);
  TGraphErrors *gphit_mc = new TGraphErrors();
  gphit_mc->SetMarkerStyle(20);
  gphit_mc->SetMarkerSize(1.5);
  gphit_mc->SetMarkerColor(1);
  gphit_mc->SetMinimum(0.0);
  gphit_mc->SetTitle("; -t (GeV/c)^{2}; N_{#phi #pi^{+} #pi^{-}}");

  // Efficiency
  TCanvas *cgphiteff = new TCanvas("cgphiteff", "cgphiteff", 900, 600);
  TGraphErrors *gphiteff = new TGraphErrors();
  gphiteff->SetMarkerStyle(20);
  gphiteff->SetMarkerSize(1.5);
  gphiteff->SetMarkerColor(kBlue);
  gphiteff->SetMinimum(0.0);
  gphiteff->SetTitle("; -t (GeV/c)^{2}; #varepsilon (%)");

  // Cross-section
  TCanvas *cgphitxsec = new TCanvas("cgphitxsec", "cgphitxsec", 900, 600);
  TGraphErrors *gphitxsec = new TGraphErrors();
  gphitxsec->SetMarkerStyle(20);
  gphitxsec->SetMarkerSize(1.5);
  gphitxsec->SetMarkerColor(1);
  gphitxsec->SetMinimum(0.0);
  gphitxsec->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)"); //#phi(1020) flux normalized yield, Yield_{#phi}

  double testxsec[nt];

  for (int i = 1; i <= nt; ++i)
  {
    cout << i << " " << flush;
    double t1 = h2phit_mc->GetXaxis()->GetBinLowEdge(i);
    double t2 = h2phit_mc->GetXaxis()->GetBinLowEdge(i) + t_step;
    double t = h2phit_mc->GetXaxis()->GetBinCenter(i);

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphit1_mc->cd(i);
    TH1D *hphit_mc_py = h2phit_mc->ProjectionY(Form("hphit_mc_py_%d", i), i, i);
    hphit_mc_py->SetTitle(Form("%.2f<-t<%.2f (GeV/c)^{2};m_{K^{+}K^{-}} (GeV/c^{2});Counts",t1,t2));
    hphit_mc_py->Draw("e");
   
    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max); // pol4
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1); // 1,1,1,1,1
    fsb_mc->SetParLimits(1, 1.015, 1.022); //fsb_mc->FixParameter(1, 1.019);
    // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_mc->FixParameter(3, 0.0042);   //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max); //pol4
    fb_mc->SetLineColor(28);
    fb_mc->SetLineStyle(2);

    hphit_mc_py->Fit("fsb_mc", "", "", mkk_min, mkk_max);
    double par_mc[fsb_mc->GetNpar()];
    fsb_mc->GetParameters(&par_mc[0]);
    fs_mc->SetParameters(&par_mc[0]);
    fb_mc->SetParameters(&par_mc[4]); //4

    fs_mc->Draw("same");
    fb_mc->Draw("same");

    double N_phit_mc = fs_mc->Integral(mkk_min, mkk_max) / hphit_mc_py->GetBinWidth(1);
    double dN_phit_mc = N_phit_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

    TLatex lat_phit_mc;
    lat_phit_mc.SetTextSize(0.07);
    lat_phit_mc.SetTextAlign(13); //align at top
    lat_phit_mc.SetNDC();
    lat_phit_mc.SetTextColor(kBlue);
    lat_phit_mc.DrawLatex(0.5, 0.80, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
    lat_phit_mc.DrawLatex(0.5, 0.72, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_mc, dN_phit_mc));
    lat_phit_mc.DrawLatex(0.5, 0.64, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
    lat_phit_mc.DrawLatex(0.5, 0.56, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phit_mc.DrawLatex(0.5, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

    gphit_mc->SetPoint(i - 1, t, N_phit_mc); // h2phit_mc->GetXaxis()->GetBinCenter(i)
    gphit_mc->SetPointError(i - 1, 0, dN_phit_mc);

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphit1->cd(i);
    TH1D *hphit_py = h2phit->ProjectionY(Form("hphit_py_%d", i), i, i);
    hphit_py->SetTitle(Form("%.2f<-t<%.2f (GeV/c)^{2};m_{K^{+}K^{-}} (GeV/c^{2});Counts",t1,t2));
    hphit_py->Draw("e");

    // +++++++++ Voigt + pol4
    TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max); //pol4
    fsb_data->SetLineColor(2);
    fsb_data->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1);  // 1,1,1,1,1
    fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); //fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_data = new TF1("fs_data", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_data->SetLineColor(4);

    TF1 *fb_data = new TF1("fb_data", "pol4(4)", mkk_min, mkk_max); //pol4(4)
    fb_data->SetLineColor(28);
    fb_data->SetLineStyle(2);

    hphit_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    double par_data[fsb_data->GetNpar()];
    fsb_data->GetParameters(&par_data[0]);
    fs_data->SetParameters(&par_data[0]);
    fb_data->SetParameters(&par_data[4]); //4

    fs_data->Draw("same");
    fb_data->Draw("same");

    double N_phit_data = fs_data->Integral(mkk_min, mkk_max) / hphit_py->GetBinWidth(1);
    double dN_phit_data = N_phit_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    TLatex lat_phit_data;
    lat_phit_data.SetTextSize(0.07);
    lat_phit_data.SetTextAlign(13); //align at top
    lat_phit_data.SetNDC();
    lat_phit_data.SetTextColor(kBlue);
    lat_phit_data.DrawLatex(0.5, 0.80, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phit_data.DrawLatex(0.5, 0.72, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_data, dN_phit_data));
    lat_phit_data.DrawLatex(0.5, 0.64, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phit_data.DrawLatex(0.5, 0.30, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phit_data.DrawLatex(0.5, 0.22, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

    gphit->SetPoint(i - 1, t, N_phit_data);
    gphit->SetPointError(i - 1, 0, dN_phit_data);

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phi = N_phit_mc / h_t_Thrown->GetBinContent(i); // Efficiency = N_observed/N_generated
    double deff_phi = eff_phi * (dN_phit_mc / N_phit_mc);
    gphiteff->SetPoint(i - 1, t, eff_phi*100);
    gphiteff->SetPointError(i - 1, 0, deff_phi*100); //->GetBinError(i)

    // +++++++++++++++++++ Systematic Errors 

    // ------------------------ 2016

    double optimal_val_16[10] = {16.27, 18.99, 12.77, 8.97, 5.70, 3.48, 2.29, 1.57, 0.93, 0.70};
    // Bkg pol order
    double bkg_poln_16_1[10] = {16.24, 18.96, 12.76, 8.99, 5.78, 3.62, 2.40, 1.56, 1.03, 0.78}; // pol3
    double bkg_poln_16_2[10] = {16.35, 19.01, 12.79, 9.02, 5.68, 3.49, 2.29, 1.52, 0.94, 0.71}; // pol5
    // fit range
    double fit_range_16_1[10] = {16.97, 19.40, 13.11, 9.25, 5.80, 3.60, 2.37, 1.56, 0.99, 0.77}; //[0.99,1.15]
    double fit_range_16_2[10] = {16.45, 19.18, 13.02, 9.20, 5.71, 3.52, 2.34, 1.55, 0.93, 0.74}; //[0.99,1.25]
    // Binning
    double binning_16_1[10] = {16.24, 18.96, 12.75, 9.00, 5.72, 3.48, 2.30, 1.57, 0.95, 0.71}; // 90
    double binning_16_2[10] = {16.21, 18.88, 12.78, 9.03, 5.72, 3.49, 2.33, 1.59, 0.96, 0.74}; //110
    // beam bunches
    double beambun_16_1[10] = {15.85, 18.83, 12.78, 8.96, 5.72, 3.51, 2.28, 1.50, 0.88, 0.74}; // w2
    double beambun_16_2[10] = {16.12, 19.03, 12.77, 9.01, 5.67, 3.54, 2.29, 1.53, 0.92, 0.70}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_16_opt[10] = {16.34, 19.00, 12.83, 9.01, 5.72, 3.50, 2.33, 1.61, 0.94, 0.73};
    // abs(p_dttof)<0.3
    double p_dttof_16_1[10] = {16.36, 18.88, 12.61, 8.91, 5.64, 3.41, 2.30, 1.59, 0.96, 0.70}; // 0.2
    double p_dttof_16_2[10] = {16.39, 19.13, 12.95, 9.08, 5.72, 3.56, 2.36, 1.66, 1.00, 0.77}; // 0.4
    // kin_chisq<55
    double kin_chisq_16_1[10] = {15.93, 18.74, 12.37, 8.98, 5.55, 3.47, 2.32, 1.51, 1.00, 0.88}; // 45
    double kin_chisq_16_2[10] = {16.60, 19.26, 13.06, 9.14, 5.80, 3.65, 2.49, 1.59, 0.97, 0.89}; // 65
    // abs(mm2)<0.035
    double mm2_16_1[10] = {16.14, 18.99, 12.77, 8.91, 5.71, 3.52, 2.36, 1.57, 0.96, 0.74}; // 0.025
    double mm2_16_2[10] = {16.25, 19.05, 12.82, 8.98, 5.74, 3.49, 2.29, 1.61, 0.93, 0.73}; // 0.045

    // ------------------------ 2017

    double optimal_val_17[10] = {12.67, 14.42, 10.27, 6.90, 4.33, 2.86, 1.96, 1.41, 1.00, 0.78};
    // Bkg pol order
    double bkg_poln_17_1[10] = {12.64, 14.40, 10.26, 6.89, 4.31, 2.87, 1.99, 1.45, 1.02, 0.81}; // pol3
    double bkg_poln_17_2[10] = {12.69, 14.44, 10.28, 6.90, 4.33, 2.88, 1.96, 1.42, 0.99, 0.79}; // pol5
    // fit range
    double fit_range_17_1[10] = {13.16, 14.75, 10.31, 7.00, 4.42, 2.91, 1.98, 1.46, 1.00, 0.80}; // [0.99,1.15]
    double fit_range_17_2[10] = {12.81, 14.58, 10.37, 6.95, 4.37, 2.90, 1.95, 1.42, 1.00, 0.79}; // [0.99,1.25]
    // Binning
    double binning_17_1[10] = {12.70, 14.41, 10.26, 6.91, 4.33, 2.86, 1.94, 1.42, 1.00, 0.78}; // 90
    double binning_17_2[10] = {12.68, 14.41, 10.24, 6.90, 4.33, 2.86, 1.95, 1.41, 0.99, 0.78}; // 110
    // beam bunches
    double beambun_17_1[10] = {12.61, 14.34, 10.23, 6.86, 4.31, 2.88, 1.96, 1.40, 0.99, 0.79}; // w2
    double beambun_17_2[10] = {12.65, 14.38, 10.26, 6.88, 4.33, 2.85, 1.96, 1.41, 1.00, 0.79}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_17_opt[10] = {12.63, 14.42, 10.28, 6.92, 4.35, 2.90, 1.98, 1.43, 1.01, 0.80};
    // abs(p_dttof)<0.3
    double p_dttof_17_1[10] = {12.60, 14.27, 10.12, 6.81, 4.27, 2.80, 1.90, 1.40, 0.95, 0.77}; // 0.2
    double p_dttof_17_2[10] = {12.65, 14.55, 10.38, 6.98, 4.39, 2.93, 2.00, 1.46, 1.04, 0.82}; // 0.4
    // kin_chisq<55
    double kin_chisq_17_1[10] = {12.23, 14.15, 10.10, 6.82, 4.33, 2.87, 1.93, 1.42, 1.02, 0.79}; // 45
    double kin_chisq_17_2[10] = {12.94, 14.64, 10.38, 6.96, 4.39, 2.87, 2.00, 1.42, 1.02, 0.83}; // 65
    // abs(mm2)<0.035
    double mm2_17_1[10] = {12.48, 14.36, 10.23, 6.89, 4.33, 2.86, 1.99, 1.44, 1.02, 0.81}; // 0.025
    double mm2_17_2[10] = {12.67, 14.42, 10.30, 6.93, 4.35, 2.89, 1.97, 1.42, 1.02, 0.81}; // 0.045

    // ------------------------ 2018 Spring

    double optimal_val_18[10] = {13.88, 14.39, 9.53, 6.33, 4.13, 2.74, 1.73, 1.27, 0.85, 0.57};
    // Bkg pol order
    double bkg_poln_18_1[10] = {13.82, 14.36, 9.50, 6.32, 4.12, 2.74, 1.73, 1.26, 0.85, 0.57}; // pol3
    double bkg_poln_18_2[10] = {13.92, 14.41, 9.54, 6.34, 4.13, 2.75, 1.74, 1.27, 0.85, 0.57}; // pol5
    // fit range
    double fit_range_18_1[10] = {14.29, 14.39, 9.52, 6.29, 4.21, 2.74, 1.78, 1.32, 0.89, 0.62}; //[0.99,1.15]
    double fit_range_18_2[10] = {14.17, 14.59, 9.65, 6.43, 4.18, 2.76, 1.75, 1.29, 0.87, 0.59}; //[0.99,1.25]
    // Binning
    double binning_18_1[10] = {13.88, 14.38, 9.53, 6.33, 4.12, 2.74, 1.73, 1.27, 0.85, 0.58}; // 90
    double binning_18_2[10] = {13.87, 14.37, 9.53, 6.33, 4.13, 2.74, 1.73, 1.27, 0.85, 0.57}; // 110
    // beam bunches
    double beambun_18_1[10] = {13.78, 14.36, 9.49, 6.28, 4.11, 2.73, 1.69, 1.25, 0.86, 0.56}; // w2
    double beambun_18_2[10] = {13.83, 14.39, 9.50, 6.32, 4.11, 2.73, 1.72, 1.27, 0.85, 0.56}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18_opt[10] = {13.66, 14.18, 9.39, 6.24, 4.06, 2.72, 1.73, 1.26, 0.85, 0.56};
    // abs(p_dttof)<0.3
    double p_dttof_18_1[10] = {13.64, 14.11, 9.32, 6.17, 4.00, 2.66, 1.69, 1.22, 0.83, 0.54}; // 0.2
    double p_dttof_18_2[10] = {13.67, 14.25, 9.45, 6.28, 4.09, 2.73, 1.74, 1.28, 0.86, 0.57}; // 0.4
    // kin_chisq<55
    double kin_chisq_18_1[10] = {12.95, 13.48, 8.98, 5.94, 3.87, 2.54, 1.63, 1.20, 0.86, 0.54}; // 45
    double kin_chisq_18_2[10] = {14.15, 14.74, 9.72, 6.46, 4.23, 2.80, 1.78, 1.28, 0.89, 0.61}; // 65
    // abs(mm2)<0.035
    double mm2_18_1[10] = {13.64, 14.19, 9.41, 6.25, 4.08, 2.74, 1.73, 1.26, 0.86, 0.57}; // 0.025
    double mm2_18_2[10] = {13.67, 14.19, 9.38, 6.24, 4.06, 2.72, 1.72, 1.28, 0.87, 0.56}; // 0.045

    // ------------------------ 2018 Fall

    double optimal_val_18l[10] = {14.59, 15.86, 10.79, 7.12, 4.55, 2.87, 1.90, 1.32, 0.90, 0.68};
    // Bkg pol order
    double bkg_poln_18l_1[10] = {14.57, 15.84, 10.78, 7.09, 4.55, 2.87, 1.90, 1.32, 0.91, 0.68}; // pol3
    double bkg_poln_18l_2[10] = {14.62, 15.88, 10.81, 7.13, 4.56, 2.87, 1.90, 1.32, 0.91, 0.68}; // pol5
    // fit range
    double fit_range_18l_1[10] = {15.14, 15.96, 10.81, 7.12, 4.54, 2.91, 1.89, 1.36, 0.93, 0.71}; // [0.99,1.15]
    double fit_range_18l_2[10] = {14.78, 16.00, 10.88, 7.19, 4.59, 2.88, 1.91, 1.34, 0.91, 0.69}; // [0.99,1.25]
    // Binning
    double binning_18l_1[10] = {14.60, 15.86, 10.80, 7.12, 4.56, 2.87, 1.90, 1.32, 0.91, 0.68}; // 90
    double binning_18l_2[10] = {14.58, 15.85, 10.79, 7.12, 4.55, 2.86, 1.90, 1.31, 0.91, 0.68}; // 110
    // beam bunches
    double beambun_18l_1[10] = {14.49, 15.76, 10.72, 7.09, 4.53, 2.87, 1.90, 1.32, 0.92, 0.70}; // w2
    double beambun_18l_2[10] = {14.58, 15.85, 10.75, 7.11, 4.53, 2.86, 1.89, 1.32, 0.91, 0.69}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18l_opt[10] = {14.55, 15.82, 10.78, 7.13, 4.58, 2.89, 1.91, 1.34, 0.92, 0.69};
    // abs(p_dttof)<0.3
    double p_dttof_18l_1[10] = {14.52, 15.71, 10.67, 7.03, 4.50, 2.84, 1.88, 1.30, 0.88, 0.65}; // 0.2
    double p_dttof_18l_2[10] = {14.58, 15.91, 10.84, 7.17, 4.60, 2.92, 1.93, 1.35, 0.92, 0.69}; // 0.4
    // kin_chisq<55
    double kin_chisq_18l_1[10] = {14.16, 15.55, 10.68, 7.06, 4.51, 2.86, 1.90, 1.33, 0.92, 0.67}; // 45
    double kin_chisq_18l_2[10] = {14.87, 16.03, 10.86, 7.21, 4.61, 2.90, 1.94, 1.32, 0.91, 0.69}; // 65
    // abs(mm2)<0.035
    double mm2_18l_1[10] = {14.45, 15.78, 10.75, 7.11, 4.56, 2.89, 1.92, 1.34, 0.92, 0.68}; // 0.025
    double mm2_18l_2[10] = {14.58, 15.85, 10.80, 7.13, 4.58, 2.90, 1.93, 1.34, 0.92, 0.69}; // 0.045

    // --------------- total  

    double optimal_val[10], bkg_poln_1[10], bkg_poln_2[10], fit_range_1[10], fit_range_2[10], binning_1[10], binning_2[10], beambun_1[10], beambun_2[10], cutsvar_opt[10], p_dttof_1[10], p_dttof_2[10], kin_chisq_1[10], kin_chisq_2[10], mm2_1[10], mm2_2[10];

    if (name == "16")
    {
      optimal_val[i-1] = optimal_val_16[i-1];    
      bkg_poln_1[i-1] = bkg_poln_16_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_16_2[i-1];
      fit_range_1[i-1] = fit_range_16_1[i-1];
      fit_range_2[i-1] = fit_range_16_2[i-1];
      binning_1[i-1] = binning_16_1[i-1];
      binning_2[i-1] = binning_16_2[i-1];
      beambun_1[i-1] = beambun_16_1[i-1];
      beambun_2[i-1] = beambun_16_2[i-1];
      cutsvar_opt[i-1] = cutsvar_16_opt[i-1];
      p_dttof_1[i-1] = p_dttof_16_1[i-1];
      p_dttof_2[i-1] = p_dttof_16_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_16_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_16_2[i-1];
      mm2_1[i-1] = mm2_16_1[i-1];
      mm2_2[i-1] = mm2_16_2[i-1];
    }
    if (name == "17")
    {
      optimal_val[i-1] = optimal_val_17[i-1];    
      bkg_poln_1[i-1] = bkg_poln_17_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_17_2[i-1];
      fit_range_1[i-1] = fit_range_17_1[i-1];
      fit_range_2[i-1] = fit_range_17_2[i-1];
      binning_1[i-1] = binning_17_1[i-1];
      binning_2[i-1] = binning_17_2[i-1];
      beambun_1[i-1] = beambun_17_1[i-1];
      beambun_2[i-1] = beambun_17_2[i-1];
      cutsvar_opt[i-1] = cutsvar_17_opt[i-1];
      p_dttof_1[i-1] = p_dttof_17_1[i-1];
      p_dttof_2[i-1] = p_dttof_17_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_17_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_17_2[i-1];
      mm2_1[i-1] = mm2_17_1[i-1];
      mm2_2[i-1] = mm2_17_2[i-1];
    }
    if (name == "18")
    {
      optimal_val[i-1] = optimal_val_18[i-1];    
      bkg_poln_1[i-1] = bkg_poln_18_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_18_2[i-1];
      fit_range_1[i-1] = fit_range_18_1[i-1];
      fit_range_2[i-1] = fit_range_18_2[i-1];
      binning_1[i-1] = binning_18_1[i-1];
      binning_2[i-1] = binning_18_2[i-1];
      beambun_1[i-1] = beambun_18_1[i-1];
      beambun_2[i-1] = beambun_18_2[i-1];
      cutsvar_opt[i-1] = cutsvar_18_opt[i-1];
      p_dttof_1[i-1] = p_dttof_18_1[i-1];
      p_dttof_2[i-1] = p_dttof_18_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_18_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_18_2[i-1];
      mm2_1[i-1] = mm2_18_1[i-1];
      mm2_2[i-1] = mm2_18_2[i-1];
    }
    if (name == "18l")
    {
      optimal_val[i-1] = optimal_val_18l[i-1];    
      bkg_poln_1[i-1] = bkg_poln_18l_1[i-1];
      bkg_poln_2[i-1] = bkg_poln_18l_2[i-1];
      fit_range_1[i-1] = fit_range_18l_1[i-1];
      fit_range_2[i-1] = fit_range_18l_2[i-1];
      binning_1[i-1] = binning_18l_1[i-1];
      binning_2[i-1] = binning_18l_2[i-1];
      beambun_1[i-1] = beambun_18l_1[i-1];
      beambun_2[i-1] = beambun_18l_2[i-1];
      cutsvar_opt[i-1] = cutsvar_18l_opt[i-1];
      p_dttof_1[i-1] = p_dttof_18l_1[i-1];
      p_dttof_2[i-1] = p_dttof_18l_2[i-1];
      kin_chisq_1[i-1] = kin_chisq_18l_1[i-1];
      kin_chisq_2[i-1] = kin_chisq_18l_2[i-1];
      mm2_1[i-1] = mm2_18l_1[i-1];
      mm2_2[i-1] = mm2_18l_2[i-1];
    }

    // total
    double arr_bkg_poln[3] = {optimal_val[i-1], bkg_poln_1[i-1], bkg_poln_2[i-1]};
    double err_bkg_poln = TMath::RMS(3, arr_bkg_poln)/optimal_val[i-1];
    double arr_fit_range[3] = {optimal_val[i-1], fit_range_1[i-1], fit_range_2[i-1]};
    double err_fit_range = TMath::RMS(3, arr_fit_range)/optimal_val[i-1];
    double arr_binning[3] = {optimal_val[i-1], binning_1[i-1], binning_2[i-1]};
    double err_binning = TMath::RMS(3, arr_binning)/optimal_val[i-1];
    double arr_beambun[3] = {optimal_val[i-1], beambun_1[i-1], beambun_2[i-1]};
    double err_beambun = TMath::RMS(3, arr_beambun)/optimal_val[i-1];
    double arr_p_dttof[3] = {cutsvar_opt[i-1], p_dttof_1[i-1], p_dttof_2[i-1]};
    double err_p_dttof = TMath::RMS(3, arr_p_dttof)/optimal_val[i-1];
    double arr_kin_chisq[3] = {cutsvar_opt[i-1], kin_chisq_1[i-1], kin_chisq_2[i-1]};
    double err_kin_chisq = TMath::RMS(3, arr_kin_chisq)/optimal_val[i-1];
    double arr_mm2[3] = {cutsvar_opt[i-1], mm2_1[i-1], mm2_2[i-1]};
    double err_mm2 = TMath::RMS(3, arr_mm2)/optimal_val[i-1];

    double err_sys = TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range + err_binning*err_binning + err_beambun*err_beambun + err_p_dttof*err_p_dttof + err_kin_chisq*err_kin_chisq + err_mm2*err_mm2);

    // double err_sys = 0; //TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range);

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phi = h_tagged_flux->GetBinContent(1); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    // if (lumi_phi <= 0)
    // {
    //   //gphitxsec->RemovePoint(i - 1);
    //   continue;
    // }
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_phi2pi = 1e9 * N_phit_data / (eff_phi * lumi_phi * T * br_phi);
    double dxsec_phi2pi_stat = xsec_phi2pi * (dN_phit_data / N_phit_data);
    double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys*err_sys) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys*err_sys) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((dN_phit_data_sys / N_phit_data)*(dN_phit_data_sys / N_phit_data) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    double dxsec_phi2pi_tot = TMath::Sqrt(dxsec_phi2pi_stat*dxsec_phi2pi_stat + dxsec_phi2pi_sys*dxsec_phi2pi_sys);
    // double xsec_phi = N_phit / lumi_phi;
    // double dxsec_phi = dN_phit / lumi_phi;  
    gphitxsec->SetPoint(i-1, t, xsec_phi2pi);
    gphitxsec->SetPointError(i-1, 0, dxsec_phi2pi_tot);
    // gphitxsec->SetPoint(gphitxsec->GetN(), t, xsec_phi2pi);
    // gphitxsec->SetPointError(gphitxsec->GetN() - 1, 0, dxsec_phi2pi_stat);
    
    fprintf(table_xsec_phi2pi, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", t1, t2, h_t_Thrown->GetBinContent(i), N_phit_mc, dN_phit_mc, N_phit_data, dN_phit_data, eff_phi*100, deff_phi*100, xsec_phi2pi, dxsec_phi2pi_stat, dxsec_phi2pi_sys);

    fprintf(table_xsec_phi2pi_sys, "%0.2f - %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", t1, t2, err_bkg_poln * 100, err_fit_range * 100, err_binning * 100, err_beambun * 100, err_p_dttof * 100, err_kin_chisq * 100, err_mm2 * 100);

    testxsec[i-1] = xsec_phi2pi;

    // test << std::setprecision(3) << std::fixed;
    // test << "i-1 = " << i - 1 << " | xsec_phi2pi = " << xsec_phi2pi << endl;
    // " | dxsec_phi2pi_sys = " << dxsec_phi2pi_sys <<" | dxsec_phi2pi_tot = " << dxsec_phi2pi_tot <<  " | error = " << gphiexsec->GetErrorY(i - 1) << endl;
    }
    
    test << std::setprecision(2) << std::fixed;
    test << testxsec[0] << ", " << testxsec[1] << ", " << testxsec[2] << ", " << testxsec[3] << ", " << testxsec[4] << ", " << testxsec[5] << ", " << testxsec[6] << ", " << testxsec[7] << ", " << testxsec[8] << ", " << testxsec[9] << endl;

    cgphit->cd();
    gphit->Draw("AP");
    gphit->Write(Form("h%s_gphit", name.Data()), TObject::kWriteDelete);
    cgphit_mc->cd();
    gphit_mc->Draw("AP");
    gphit_mc->Write(Form("h%s_gphit_mc", name.Data()), TObject::kWriteDelete);
    cgphiteff->cd();
    gphiteff->Draw("AP");
    gphiteff->Write(Form("h%s_gphiteff", name.Data()), TObject::kWriteDelete);
    cgphitxsec->cd();
    gphitxsec->Draw("AP");
    gphitxsec->Write(Form("h%s_gphitxsec", name.Data()), TObject::kWriteDelete);
    // int j =1;
    // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);

     // *********** total Thrown momentum transfer
    TCanvas *ctot_t_tru = new TCanvas("ctot_t_tru", "ctot_t_tru", 900, 600);
    TLegend *ltot_t_tru = new TLegend(0.84, 0.80, 0.98, 0.98);
    ltot_t_tru->SetTextSize(0.06);
    ltot_t_tru->SetBorderSize(0);
    ctot_t_tru->cd();
    ctot_t_tru->SetGrid();
    // TH1F *h16_t_tru = (TH1F *)outputfig->Get("h16_t_tru");
    // cout << " ***** h16_t_tru = " << h16_t_tru << endl;
    // // h16_t_tru->SetMarkerColor(kBlue);
    // h16_t_tru->SetMarkerSize(1.5);
    // h16_t_tru->SetMarkerStyle(kFullSquare);
    TH1F *h17_t_tru = (TH1F *)outputfig->Get("h17_t_tru");
    cout << " ***** h17_t_tru = " << h17_t_tru << endl;
    h17_t_tru->SetMarkerSize(1.5);
    h17_t_tru->SetMarkerColor(kBlue);
    h17_t_tru->SetMarkerStyle(kFullCircle);
    TH1F *h18_t_tru = (TH1F *)outputfig->Get("h18_t_tru");
    cout << " ***** h18_t_tru = " << h18_t_tru << endl;
    h18_t_tru->SetMarkerColor(kRed);
    h18_t_tru->SetMarkerSize(1.5);
    h18_t_tru->SetMarkerStyle(kFullTriangleUp);
    TH1F *h18l_t_tru = (TH1F *)outputfig->Get("h18l_t_tru");
    cout << " ***** h18l_t_tru = " << h18l_t_tru << endl;
    h18l_t_tru->SetMarkerColor(kMagenta);
    h18l_t_tru->SetMarkerSize(1.5);
    h18l_t_tru->SetMarkerStyle(kOpenTriangleUp);
    h17_t_tru->SetTitle(";E_{#gamma}:Counts");
    h17_t_tru->SetMinimum(0.);
    h18_t_tru->Draw("e");
    h18l_t_tru->Draw("esame");
    h17_t_tru->Draw("esame");
    // h16_t_tru->Draw("esame");
    ctot_t_tru->Modified();
    // ltot_t_tru->AddEntry(h16_t_tru, "2016", "p");
    ltot_t_tru->AddEntry(h17_t_tru, "2017", "p");
    ltot_t_tru->AddEntry(h18_t_tru, "2018 S", "p");
    ltot_t_tru->AddEntry(h18l_t_tru, "2018 F", "p");
    ltot_t_tru->Draw();
    ctot_t_tru->Print("output/fig_phi2pi/ctot_t_tru.eps", "eps");
    ctot_t_tru->Print("output/fig_phi2pi/ctot_t_tru.png", "png");
    ctot_t_tru->Print("output/fig_phi2pi/ctot_t_tru.root", "root");

    // total Yields MC
    TMultiGraph *mgphit_mc = new TMultiGraph();
    TCanvas *cmgphit_mc = new TCanvas("cmgphit_mc", "cmgphit_mc", 900, 600);
    TLegend *lmgphit_mc = new TLegend(0.84, 0.80, 0.98, 0.98);
    lmgphit_mc->SetTextSize(0.06);
    lmgphit_mc->SetBorderSize(0);
    cmgphit_mc->cd();
    cmgphit_mc->SetGrid();
    // TGraphErrors *h16_gphit_mc = (TGraphErrors *)outputfig->Get("h16_gphit_mc");
    // cout << " ***** h16_gphit_mc = " << h16_gphit_mc << endl;
    // // h16_gphit_mc->SetMarkerColor(kBlue);
    // h16_gphit_mc->SetMarkerSize(1.5);
    // h16_gphit_mc->SetMarkerStyle(kFullSquare);
    // mgphit_mc->Add(h16_gphit_mc);
    TGraphErrors *h17_gphit_mc = (TGraphErrors *)outputfig->Get("h17_gphit_mc");
    cout << " ***** h17_gphit_mc = " << h17_gphit_mc << endl;
    h17_gphit_mc->SetMarkerSize(1.5);
    h17_gphit_mc->SetMarkerColor(kBlue);
    h17_gphit_mc->SetMarkerStyle(kFullCircle);
    mgphit_mc->Add(h17_gphit_mc);
    TGraphErrors *h18_gphit_mc = (TGraphErrors *)outputfig->Get("h18_gphit_mc");
    cout << " ***** h18_gphit_mc = " << h18_gphit_mc << endl;
    h18_gphit_mc->SetMarkerColor(kRed);
    h18_gphit_mc->SetMarkerSize(1.5);
    h18_gphit_mc->SetMarkerStyle(kFullTriangleUp);
    mgphit_mc->Add(h18_gphit_mc);
    TGraphErrors *h18l_gphit_mc = (TGraphErrors *)outputfig->Get("h18l_gphit_mc");
    cout << " ***** h18l_gphit_mc = " << h18l_gphit_mc << endl;
    h18l_gphit_mc->SetMarkerColor(kMagenta);
    h18l_gphit_mc->SetMarkerSize(1.5);
    h18l_gphit_mc->SetMarkerStyle(kOpenTriangleUp);
    mgphit_mc->Add(h18l_gphit_mc);
    mgphit_mc->SetTitle(";-t (GeV/c)^{2}; N_{#phi #pi^{+} #pi^{-}}");
    mgphit_mc->Draw("AP");
    mgphit_mc->SetMinimum(0.);
    cmgphit_mc->Modified();
    // lmgphit_mc->AddEntry(h16_gphit_mc, "2016", "p");
    lmgphit_mc->AddEntry(h17_gphit_mc, "2017", "p");
    lmgphit_mc->AddEntry(h18_gphit_mc, "2018 S", "p");
    lmgphit_mc->AddEntry(h18l_gphit_mc, "2018 F", "p");
    lmgphit_mc->Draw();
    cmgphit_mc->Print("output/fig_phi2pi/cmgphit_mc.eps", "eps");
    cmgphit_mc->Print("output/fig_phi2pi/cmgphit_mc.png", "png");
    cmgphit_mc->Print("output/fig_phi2pi/cmgphit_mc.root", "root");

    // total Yields Data
    TMultiGraph *mgphit = new TMultiGraph();
    TCanvas *cmgphit = new TCanvas("cmgphit", "cmgphit", 900, 600);
    TLegend *lmgphit = new TLegend(0.84, 0.80, 0.98, 0.98);
    lmgphit->SetTextSize(0.06);
    lmgphit->SetBorderSize(0);
    cmgphit->cd();
    cmgphit->SetGrid();
    // TGraphErrors *h16_gphit = (TGraphErrors *)outputfig->Get("h16_gphit");
    // cout << " ***** h16_gphit = " << h16_gphit << endl;
    // // h16_gphit->SetMarkerColor(kBlue);
    // h16_gphit->SetMarkerSize(1.5);
    // h16_gphit->SetMarkerStyle(kFullSquare);
    // mgphit->Add(h16_gphit);
    TGraphErrors *h17_gphit = (TGraphErrors *)outputfig->Get("h17_gphit");
    cout << " ***** h17_gphit = " << h17_gphit << endl;
    h17_gphit->SetMarkerSize(1.5);
    h17_gphit->SetMarkerColor(kBlue);
    h17_gphit->SetMarkerStyle(kFullCircle);
    mgphit->Add(h17_gphit);
    TGraphErrors *h18_gphit = (TGraphErrors *)outputfig->Get("h18_gphit");
    cout << " ***** h18_gphit = " << h18_gphit << endl;
    h18_gphit->SetMarkerColor(kRed);
    h18_gphit->SetMarkerSize(1.5);
    h18_gphit->SetMarkerStyle(kFullTriangleUp);
    mgphit->Add(h18_gphit);
    TGraphErrors *h18l_gphit = (TGraphErrors *)outputfig->Get("h18l_gphit");
    cout << " ***** h18l_gphit = " << h18l_gphit << endl;
    h18l_gphit->SetMarkerColor(kMagenta);
    h18l_gphit->SetMarkerSize(1.5);
    h18l_gphit->SetMarkerStyle(kOpenTriangleUp);
    mgphit->Add(h18l_gphit);
    mgphit->SetTitle(";-t (GeV/c)^{2}; N_{#phi #pi^{+} #pi^{-}}");
    mgphit->Draw("AP");
    mgphit->SetMinimum(0.);
    cmgphit->Modified();
    // lmgphit->AddEntry(h16_gphit, "2016", "p");
    lmgphit->AddEntry(h17_gphit, "2017", "p");
    lmgphit->AddEntry(h18_gphit, "2018 S", "p");
    lmgphit->AddEntry(h18l_gphit, "2018 F", "p");
    lmgphit->Draw();
    cmgphit->Print("output/fig_phi2pi/cmgphit.eps", "eps");
    cmgphit->Print("output/fig_phi2pi/cmgphit.png", "png");
    cmgphit->Print("output/fig_phi2pi/cmgphit.root", "root");

    // *********** total Efficiency
    TCanvas *cmgteff = new TCanvas("cmgteff", "cmgteff", 800, 800);
    cmgteff->cd();
    cmgteff->SetGrid();

    float ateff = 0.3;
    float bteff = 0.02;

    TPad *padteff1 = new TPad("padteff1", "padteff1", 0, ateff - bteff, 1, 1);
    padteff1->SetBottomMargin(bteff);
    padteff1->SetGridx();
    padteff1->SetGridy();
    cmgteff->cd();
    padteff1->Draw();

    TPad *padteff2 = new TPad("padteff2", "padteff2", 0, 0, 1, ateff * (1 - bteff));
    padteff2->SetTopMargin(0);
    padteff2->SetBottomMargin(0.4);
    padteff2->SetFillColor(0);
    padteff2->SetFillStyle(0);
    padteff2->SetGridy();
    cmgteff->cd();
    padteff2->Draw();

    padteff1->cd();
    // TMultiGraph *mgteff = new TMultiGraph();
    TLegend *lmgteff = new TLegend(0.8, 0.78, 0.93, 0.98);
    lmgteff->SetTextSize(0.06);
    lmgteff->SetBorderSize(0);
    // TGraphErrors *h16_gphiteff = (TGraphErrors *)outputfig->Get("h16_gphiteff");
    // cout << " ***** h16_gphiteff = " << h16_gphiteff << endl;
    // // h16_gphiteff->SetMarkerColor(kBlue);
    // h16_gphiteff->SetMarkerSize(1.5);
    // h16_gphiteff->SetMarkerStyle(kFullSquare);
    // h16_gphiteff->SetTitle(";-t (GeV/c)^{2}; #varepsilon (%)");
    // h16_gphiteff->SetMinimum(0.);
    TGraphErrors *h17_gphiteff = (TGraphErrors *)outputfig->Get("h17_gphiteff");
    cout << " ***** h17_gphiteff = " << h17_gphiteff << endl;
    h17_gphiteff->SetMarkerColor(kBlue);
    h17_gphiteff->SetMarkerSize(1.5);
    h17_gphiteff->SetMarkerStyle(kFullCircle);
    h17_gphiteff->SetTitle(";-t (GeV/c)^{2}; #varepsilon (%)");
    h17_gphiteff->SetMinimum(0.);
    // mgteff->Add(h17_gphiteff);
    TGraphErrors *h18_gphiteff = (TGraphErrors *)outputfig->Get("h18_gphiteff");
    cout << " ***** h18_gphiteff = " << h18_gphiteff << endl;
    h18_gphiteff->SetMarkerColor(kRed);
    h18_gphiteff->SetMarkerSize(1.5);
    h18_gphiteff->SetMarkerStyle(kFullTriangleUp);
    // mgteff->Add(h18_gphiteff);
    TGraphErrors *h18l_gphiteff = (TGraphErrors *)outputfig->Get("h18l_gphiteff");
    cout << " ***** h18l_gphiteff = " << h18l_gphiteff << endl;
    h18l_gphiteff->SetMarkerColor(kMagenta);
    h18l_gphiteff->SetMarkerSize(1.5);
    h18l_gphiteff->SetMarkerStyle(kOpenTriangleUp);
    // mgteff->Add(h18l_gphiteff);
    // mgteff->SetTitle(";-t (GeV/c)^{2}; #varepsilon (%)");
    // mgteff->Draw("AP");
    // mgteff->SetMinimum(0.);
    // h16_gphiteff->Draw("AP");
    h18l_gphiteff->Draw("AP");
    h18_gphiteff->Draw("PSAME");
    h17_gphiteff->Draw("PSAME");
    cmgteff->Modified();
    // lmgteff->AddEntry(h16_gphiteff, "2016", "p");
    lmgteff->AddEntry(h17_gphiteff, "2017", "p");
    lmgteff->AddEntry(h18_gphiteff, "2018 S", "p");
    lmgteff->AddEntry(h18l_gphiteff, "2018 F", "p");
    lmgteff->Draw();

    padteff2->cd(); // pad2 becomes the current pad
    // TGraphErrors *gr_efftrate_17 = new TGraphErrors();
    TGraphErrors *gr_efftrate_18 = new TGraphErrors();
    TGraphErrors *gr_efftrate_18l = new TGraphErrors();
    // gr_efftrate_17->SetMarkerColor(kBlue);
    // gr_efftrate_17->SetMarkerStyle(2);
    // gr_efftrate_17->SetMarkerSize(2); 
    gr_efftrate_18->SetMarkerColor(kRed);
    gr_efftrate_18->SetMarkerStyle(2);
    gr_efftrate_18->SetMarkerSize(2);
    gr_efftrate_18l->SetMarkerColor(kMagenta);
    gr_efftrate_18l->SetMarkerStyle(2);
    gr_efftrate_18l->SetMarkerSize(2);

    for (int i = 1; i <= nt; ++i)
    {
      double t1 = h2phit_mc->GetXaxis()->GetBinLowEdge(i);
      double t2 = h2phit_mc->GetXaxis()->GetBinLowEdge(i) + t_step;
      double t = h2phit_mc->GetXaxis()->GetBinCenter(i);
      
      // double efftrate_17 = (h16_gphiteff->GetY()[i-1] - h17_gphiteff->GetY()[i-1])/(h16_gphiteff->GetY()[i-1]);
      double efftrate_18 = (h17_gphiteff->GetY()[i-1] - h18_gphiteff->GetY()[i-1])/(h17_gphiteff->GetY()[i-1]);
      double efftrate_18l = (h17_gphiteff->GetY()[i-1] - h18l_gphiteff->GetY()[i-1])/(h17_gphiteff->GetY()[i-1]);
      // cout<<"efftrate_17 = "<<efftrate_17<<" | efftrate_18 = "<<efftrate_18<<" | efftrate_18l = "<<efftrate_18l<<endl;
      // gr_efftrate_17->SetPoint(i-1, t, abs(efftrate_17));
      // gr_efftrate_17->SetPointError(i-1, 0, 0);
      gr_efftrate_18->SetPoint(i-1, t, abs(efftrate_18));
      gr_efftrate_18->SetPointError(i-1, 0, 0);
      gr_efftrate_18l->SetPoint(i-1, t, abs(efftrate_18l));
      gr_efftrate_18l->SetPointError(i-1, 0, 0);
    }

    // gr_efftrate_17->Draw("AP");
    gr_efftrate_18->Draw("AP");
    gr_efftrate_18l->Draw("PSAME");
    TLine lefft;
    lefft.SetLineStyle(kDashed);
    lefft.DrawLine(t_min, 0, t_max, 0);

    // Y axis ratio plot settings
    gr_efftrate_18->GetHistogram()->SetStats(0);
    gr_efftrate_18->SetTitle(";-t (GeV/c)^{2};Ratio");
    gr_efftrate_18->GetYaxis()->SetNdivisions(506);
    gr_efftrate_18->GetYaxis()->SetTitleSize(35);
    gr_efftrate_18->GetYaxis()->SetTitleFont(43);
    gr_efftrate_18->GetYaxis()->SetTitleOffset(1.55);
    gr_efftrate_18->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efftrate_18->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_efftrate_18->GetXaxis()->SetTitleSize(38);
    gr_efftrate_18->GetXaxis()->SetTitleFont(43);
    gr_efftrate_18->GetXaxis()->SetTitleOffset(3.);
    gr_efftrate_18->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efftrate_18->GetXaxis()->SetLabelSize(38);
    gr_efftrate_18->GetHistogram()->SetMaximum(1.0);   // along          
    gr_efftrate_18->GetHistogram()->SetMinimum(-1.0);  //   Y  
    gPad->Update();
    gPad->Modified();    
    cmgteff->Print("output/fig_phi2pi/cmgteff.eps", "eps");
    cmgteff->Print("output/fig_phi2pi/cmgteff.png", "png");
    cmgteff->Print("output/fig_phi2pi/cmgteff.root", "root");

    // ************ total Cross-sections

    TCanvas *cmgtxsec = new TCanvas("cmgtxsec", "cmgtxsec", 800, 800);//1000, 600
    // cmgtxsec->Divide(2);
    cmgtxsec->cd();
    cmgtxsec->SetGrid();

    float atxsec = 0.3;
    float btxsec = 0.02;

    TPad *padtxsec1 = new TPad("padtxsec1", "padtxsec1", 0, atxsec - btxsec, 1, 1);
    padtxsec1->SetBottomMargin(btxsec);
    padtxsec1->SetGridx();
    padtxsec1->SetGridy();
    cmgtxsec->cd();
    padtxsec1->Draw();
   
    TPad *padtxsec2 = new TPad("padtxsec2", "padtxsec2", 0, 0, 1, atxsec * (1 - btxsec));
    padtxsec2->SetTopMargin(0);
    padtxsec2->SetBottomMargin(0.4);//0
    padtxsec2->SetFillColor(0);
    padtxsec2->SetFillStyle(0);
    padtxsec2->SetGridy();
    padtxsec2->Draw();

    padtxsec1->cd();
    // TMultiGraph *mgtxsec = new TMultiGraph();
    TLegend *lmgtxsec = new TLegend(0.8, 0.78, 0.93, 0.98);
    lmgtxsec->SetTextSize(0.06);
    lmgtxsec->SetBorderSize(0);
    // TGraphErrors *h16_gphitxsec = (TGraphErrors *)outputfig->Get("h16_gphitxsec");
    // cout << " ***** h16_gphitxsec = " << h16_gphitxsec << endl;
    // h16_gphitxsec->SetMarkerStyle(kFullSquare);
    // h16_gphitxsec->SetMarkerSize(1.5);
    // h16_gphitxsec->SetMinimum(0.);
    // h16_gphitxsec->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");
    // h16_gphitxsec->SetMarkerColor(kBlue);
    // mgtxsec->Add(h16_gphitxsec);
    TGraphErrors *h17_gphitxsec = (TGraphErrors *)outputfig->Get("h17_gphitxsec");
    cout << " ***** h17_gphitxsec = " << h17_gphitxsec << endl;
    h17_gphitxsec->SetMarkerStyle(kFullCircle);
    h17_gphitxsec->SetMarkerSize(1.5);
    h17_gphitxsec->SetMarkerColor(kBlue);
    h17_gphitxsec->SetLineColor(kBlue);
    h17_gphitxsec->SetMinimum(0.);
    h17_gphitxsec->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");
    // mgtxsec->Add(h17_gphitxsec);
    TGraphErrors *h18_gphitxsec = (TGraphErrors *)outputfig->Get("h18_gphitxsec");
    cout << " ***** h18_gphitxsec = " << h18_gphitxsec << endl;
    h18_gphitxsec->SetMarkerStyle(kFullTriangleUp);//SetMarkerColor(kRed);
    h18_gphitxsec->SetMarkerSize(1.5);
    h18_gphitxsec->SetMarkerColor(kRed);
    h18_gphitxsec->SetLineColor(kRed);
    // mgtxsec->Add(h18_gphitxsec);
    TGraphErrors *h18l_gphitxsec = (TGraphErrors *)outputfig->Get("h18l_gphitxsec");
    cout << " ***** h18l_gphitxsec = " << h18l_gphitxsec << endl;
    h18l_gphitxsec->SetMarkerStyle(kOpenTriangleUp);//SetMarkerColor(kMagenta);
    h18l_gphitxsec->SetMarkerSize(1.5);
    h18l_gphitxsec->SetMarkerColor(kMagenta);
    h18l_gphitxsec->SetLineColor(kMagenta);
    // mgtxsec->Add(h18l_gphitxsec);
    // mgtxsec->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");
    // mgtxsec->Draw("AP");
    // mgtxsec->SetMinimum(0.);
    // h16_gphitxsec->Draw("AP");
    h18l_gphitxsec->Draw("AP");
    h18_gphitxsec->Draw("PSAME");
    h17_gphitxsec->Draw("PSAME");
    cmgtxsec->Modified();
    // lmgtxsec->AddEntry(h16_gphitxsec, "2016", "p");
    lmgtxsec->AddEntry(h17_gphitxsec, "2017", "p");
    lmgtxsec->AddEntry(h18_gphitxsec, "2018 S", "p");
    lmgtxsec->AddEntry(h18l_gphitxsec, "2018 F", "p");
    lmgtxsec->Draw();

    padtxsec2->cd(); // pad2 becomes the current pad
    // TGraphErrors *gr_txsecrate_17 = new TGraphErrors();
    TGraphErrors *gr_txsecrate_18 = new TGraphErrors();
    TGraphErrors *gr_txsecrate_18l = new TGraphErrors();
    gr_txsecrate_18->SetMarkerColor(kBlue);
    gr_txsecrate_18->SetMarkerStyle(2);
    gr_txsecrate_18->SetMarkerSize(2); 
    gr_txsecrate_18->SetMarkerColor(kRed);
    gr_txsecrate_18->SetMarkerStyle(2);
    gr_txsecrate_18->SetMarkerSize(2); 
    gr_txsecrate_18l->SetMarkerColor(kMagenta);
    gr_txsecrate_18l->SetMarkerStyle(2);
    gr_txsecrate_18l->SetMarkerSize(2);

    for (int i = 1; i <= nt; ++i)
    {
      double t1 = h2phit->GetXaxis()->GetBinLowEdge(i);
      double t2 = h2phit->GetXaxis()->GetBinLowEdge(i) + t_step;
      double t = h2phit_mc->GetXaxis()->GetBinCenter(i);

      // double txsecrate_17 = (h16_gphitxsec->GetY()[i-1] - h17_gphitxsec->GetY()[i-1])/(h16_gphitxsec->GetY()[i-1]);
      double txsecrate_18 = (h17_gphitxsec->GetY()[i-1] - h18_gphitxsec->GetY()[i-1])/(h17_gphitxsec->GetY()[i-1]);
      double txsecrate_18l = (h17_gphitxsec->GetY()[i-1] - h18l_gphitxsec->GetY()[i-1])/(h17_gphitxsec->GetY()[i-1]);
      // cout<<"txsecrate_17 = "<<txsecrate_17<<" | txsecrate_18 = "<<txsecrate_18<<" | txsecrate_18l = "<<txsecrate_18l<<endl;
      // gr_txsecrate_17->SetPoint(i-1, t, abs(txsecrate_17));
      // gr_txsecrate_17->SetPointError(i-1, 0, 0);
      gr_txsecrate_18->SetPoint(i-1, t, abs(txsecrate_18));
      gr_txsecrate_18->SetPointError(i-1, 0, 0);
      gr_txsecrate_18l->SetPoint(i-1, t, abs(txsecrate_18l));
      gr_txsecrate_18l->SetPointError(i-1, 0, 0);
    }

    // gr_txsecrate_17->Draw("AP");
    gr_txsecrate_18->Draw("AP");
    gr_txsecrate_18l->Draw("PSAME");
    TLine ltxsec;
    ltxsec.SetLineStyle(kDashed);
    ltxsec.DrawLine(t_min, 0, t_max, 0);

    // Y axis ratio plot settings
    gr_txsecrate_18->GetHistogram()->SetStats(0);
    gr_txsecrate_18->SetTitle(";-t (GeV/c)^{2};Ratio");
    gr_txsecrate_18->GetYaxis()->SetNdivisions(506);
    gr_txsecrate_18->GetYaxis()->SetTitleSize(35);//24
    gr_txsecrate_18->GetYaxis()->SetTitleFont(43);
    gr_txsecrate_18->GetYaxis()->SetTitleOffset(1.55);//1.55
    gr_txsecrate_18->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_txsecrate_18->GetYaxis()->SetLabelSize(35);//30
    // X axis ratio plot settings
    gr_txsecrate_18->GetXaxis()->SetTitleSize(38);//35
    gr_txsecrate_18->GetXaxis()->SetTitleFont(43);
    gr_txsecrate_18->GetXaxis()->SetTitleOffset(3.0); //3.
    gr_txsecrate_18->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_txsecrate_18->GetXaxis()->SetLabelSize(38);//35
    
    gr_txsecrate_18->GetHistogram()->SetMaximum(0.5);   // along          
    gr_txsecrate_18->GetHistogram()->SetMinimum(-0.5);  //   Y  
    gPad->Update();
    gPad->Modified();
    cmgtxsec->Update();
    cmgtxsec->Print("output/fig_phi2pi/cmgtxsec.eps", "eps");
    cmgtxsec->Print("output/fig_phi2pi/cmgtxsec.png", "png");
    cmgtxsec->Print("output/fig_phi2pi/cmgtxsec.root", "root");

      // +++++++++++++++++++ PDG cross section average
    TCanvas *cgtxsec_avg = new TCanvas("cgtxsec_avg", "cgtxsec_avg", 900, 600);
    cgtxsec_avg->cd();
    cgtxsec_avg->SetGrid();
    TLegend *lgtxsec_avg = new TLegend(0.67, 0.78, 0.93, 0.98);
    lgtxsec_avg->SetTextSize(0.06);
    lgtxsec_avg->SetBorderSize(0);
    TGraphErrors *gphitxsec_avg = new TGraphErrors();
    gphitxsec_avg->SetMarkerStyle(20);
    gphitxsec_avg->SetMarkerSize(1.5);
    gphitxsec_avg->SetMarkerColor(28);
    gphitxsec_avg->SetLineColor(28);
    gphitxsec_avg->SetMinimum(0.0);
    gphitxsec_avg->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");
    // TGraphErrors *gphitxsec_17corr = new TGraphErrors();
    // gphitxsec_17corr->SetMarkerStyle(kOpenCircle);
    // gphitxsec_17corr->SetMarkerSize(1.5);
    // gphitxsec_17corr->SetMarkerColor(kBlue);
    // gphitxsec_17corr->SetLineColor(kBlue);
    // gphitxsec_17corr->SetMinimum(0.0);
    // gphitxsec_17corr->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");

    // vector<double> xi = {1.0,0.8,1.1}, dxi = {0.1,0.2,0.05};
    // pair<double, double> x = pdgavg(xi, dxi);
    // Eg = 0;
    vector<double> xsec(3), dxsec(3);
    pair<double, double> xsec_pair;

    for (int i = 1; i <= nt; ++i)
    {
      xsec = {h17_gphitxsec->GetY()[i - 1], h18_gphitxsec->GetY()[i - 1], h18l_gphitxsec->GetY()[i - 1]};
      dxsec = {h17_gphitxsec->GetErrorY(i - 1), h18_gphitxsec->GetErrorY(i - 1), h18l_gphitxsec->GetErrorY(i - 1)};

      double t1 = h2phit_mc->GetXaxis()->GetBinLowEdge(i);
      double t2 = h2phit_mc->GetXaxis()->GetBinLowEdge(i) + t_step;
      double t = h2phit_mc->GetXaxis()->GetBinCenter(i);
      // if (i == 1)
      //   t = t1;
      // else
      //   t = t2;

      xsec_pair = pdgavg(xsec, dxsec);  
      
      gphitxsec_avg->SetPoint(i-1, t, xsec_pair.first);
      gphitxsec_avg->SetPointError(i-1, 0, xsec_pair.second);

      double w_17 = 0.63; //0.63
      // if(i ==5)
      // {
      //   w_17 = xsec_pair.first/h17_gphitxsec->GetY()[i-1];
      // } 
      
      // gphitxsec_17corr->SetPoint(i-1, t, w_17*h17_gphitxsec->GetY()[i - 1]);
      // gphitxsec_17corr->SetPointError(i-1, 0, w_17*h17_gphitxsec->GetErrorY(i - 1));

      // test << std::setprecision(2) << std::fixed;
      // test << "i-1 = " << i - 1 << " | xsec[0] = " << xsec[0] << " | dxsec[0] = " << dxsec[0] << " | xsec[1] = " << xsec[1] << " | dxsec[1] = " << dxsec[1] <<  " | xsec[2] = " << xsec[2] << " | dxsec[2] = " << dxsec[2] << " | xsec_pair.first = " << xsec_pair.first << " | xsec_pair.second = " << xsec_pair.second <<  " | w_17 = " << w_17 << endl;
    }

    // h17_gphitxsec->Draw("AP");
    gphitxsec_avg->Draw("AP");
    // gphitxsec_17corr->Draw("PSAME");
    // lgtxsec_avg->AddEntry(h17_gphitxsec, "2017", "p");
    lgtxsec_avg->AddEntry(gphitxsec_avg, "#bar{#sigma}_{(2016, 2018 S, 2018 F)}", "p");
    // lgtxsec_avg->AddEntry(gphitxsec_17corr, "2017 corrected", "p");
    lgtxsec_avg->Draw();

    cgtxsec_avg->Print("output/fig_phi2pi/cgtxsec_avg.eps", "eps");
    cgtxsec_avg->Print("output/fig_phi2pi/cgtxsec_avg.png", "png");
    cgtxsec_avg->Print("output/fig_phi2pi/cgtxsec_avg.root", "root");

    cphit->Print(Form("output/fig_phi2pi/c2_phit_%s.root", name.Data()), "root");
    cphit->Print(Form("output/fig_phi2pi/c2_phit_%s.eps", name.Data()), "eps");
    cphit->Print(Form("output/fig_phi2pi/c2_phit_%s.png", name.Data()), "png");
    cphit1->Print(Form("output/fig_phi2pi/c_phit1_%s.root", name.Data()), "root");
    cphit1->Print(Form("output/fig_phi2pi/c_phit1_%s.eps", name.Data()), "eps");
    cphit1->Print(Form("output/fig_phi2pi/c_phit1_%s.png", name.Data()), "png");
    cgphit->Print(Form("output/fig_phi2pi/c_gphit_%s.root", name.Data()), "root");
    cgphit->Print(Form("output/fig_phi2pi/c_gphit_%s.eps", name.Data()), "eps");
    cgphit->Print(Form("output/fig_phi2pi/c_gphit_%s.png", name.Data()), "png");
    cphit_mc->Print(Form("output/fig_phi2pi/c2_phit_mc_%s.root", name.Data()), "root");
    cphit_mc->Print(Form("output/fig_phi2pi/c2_phit_mc_%s.eps", name.Data()), "eps");
    cphit_mc->Print(Form("output/fig_phi2pi/c2_phit_mc_%s.png", name.Data()), "png");
    cphit1_mc->Print(Form("output/fig_phi2pi/c_phit1_mc_%s.root", name.Data()), "root");
    cphit1_mc->Print(Form("output/fig_phi2pi/c_phit1_mc_%s.eps", name.Data()), "eps");
    cphit1_mc->Print(Form("output/fig_phi2pi/c_phit1_mc_%s.png", name.Data()), "png");
    cgphit_mc->Print(Form("output/fig_phi2pi/c_gphit_mc_%s.root", name.Data()), "root");
    cgphit_mc->Print(Form("output/fig_phi2pi/c_gphit_mc_%s.eps", name.Data()), "eps");
    cgphit_mc->Print(Form("output/fig_phi2pi/c_gphit_mc_%s.png", name.Data()), "png");
    cgphiteff->Print(Form("output/fig_phi2pi/c_gphiteff_%s.root", name.Data()), "root");
    cgphiteff->Print(Form("output/fig_phi2pi/c_gphiteff_%s.eps", name.Data()), "eps");
    cgphiteff->Print(Form("output/fig_phi2pi/c_gphiteff_%s.png", name.Data()), "png");
    cgphitxsec->Print(Form("output/fig_phi2pi/c_gphitxsec_%s.root", name.Data()), "root");
    cgphitxsec->Print(Form("output/fig_phi2pi/c_gphitxsec_%s.eps", name.Data()), "eps");
    cgphitxsec->Print(Form("output/fig_phi2pi/c_gphitxsec_%s.png", name.Data()), "png");

    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_phi2pi/c_tagged_flux_%s.png", name.Data()), "png");
    c_h_t_Thrown->Print(Form("output/fig_phi2pi/c_h_t_Thrown_%s.root", name.Data()), "root");
    c_h_t_Thrown->Print(Form("output/fig_phi2pi/c_h_t_Thrown_%s.eps", name.Data()), "eps");
    c_h_t_Thrown->Print(Form("output/fig_phi2pi/c_h_t_Thrown_%s.png", name.Data()), "png");

    outputfig->Print();

    fprintf(table_xsec_phi2pi, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_xsec_phi2pi);
    gSystem->Exec(Form("pdflatex table_xsec_phi2pi_t_%s.tex", name.Data()));

    fprintf(table_xsec_phi2pi_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_xsec_phi2pi_sys);
    gSystem->Exec(Form("pdflatex table_xsec_phi2pi_sys_t_%s.tex", name.Data()));

    // table_phi << "\\hline" << endl;
    // table_phi << "\\end{tabularx}" << endl;
    // table_phi << "\\end{center}" << endl;
    // table_phi << "\\end{minipage}" << endl;
    // table_phi << "\\end{table}" << endl;
    // table_phi << "\\end{document}" << endl;
    // table_phi.close();
    // gSystem->Exec("pdflatex table_phi.tex");

}
