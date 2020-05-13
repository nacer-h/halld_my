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
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_chi100_flat.root", name.Data())); // _chi100 ,20200506/
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_phi2pi_%s_chi100_flat.root", name.Data())); // _chi100  ver03_11/
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

  TFile *outputfig = new TFile("output/fig_xsec_phi2pi/xsec_phi2pi.root", "UPDATE");

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

  double mkk_min = 0.99, mkk_max = 1.2; // mkk_min = 0.99, mkk_max = 1.2 [0.99, 1.15] [0.99, 1.25]
  double Eg_min = 6.5, Eg_max = 11.6;
  double Eg_step = (Eg_max-Eg_min)/ne;
  double t_min = 0.0, t_max = 4.0;
  double t_step = (t_max-t_min)/nt;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetLabelSize(0.07, "XYZ");
  // gStyle->SetTitleAlign(13);

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
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035)");// && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(kpkm_uni && abs(rf_dt)<1.5*4.008), w4*(kpkm_uni && abs(rf_dt)<2.5*4.008)
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
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035)");
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
   
    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);//4
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1); // voigt
    fsb_mc->SetParLimits(1, 1.015, 1.022); //fsb_mc->FixParameter(1, 1.019);
    // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.0001, 0.01);
    fsb_mc->FixParameter(3, 0.0042);   //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max);
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
    TF1 *fsb_data = new TF1("fsb_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    fsb_data->SetLineColor(2);
    fsb_data->SetParameters(1433, 1.019, 0.003, 0.0042, 1,1,1,1,1);  // voigt
    // fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); //fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb_data->FixParameter(1, fsb_mc->GetParameter(1)); // mean+0.25*sigma
    fsb_data->FixParameter(2, fsb_mc->GetParameter(2)); //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_data->FixParameter(3, fsb_mc->GetParameter(3)); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

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

    // ++++++++++++++++++++++++++++++++++++ Systematic Errors

    // ------------------------ 2016

    double optimal_val_16[10] = {65.632, 60.481, 67.237, 55.121, 51.016, 52.497, 51.333, 51.765, 53.202, 41.924};
    // Bkg pol order
    double bkg_poln_16_1[10] = {65.486, 63.168, 69.457, 56.988, 50.880, 54.286, 51.244, 53.964, 55.461, 41.800}; // pol3
    double bkg_poln_16_2[10] = {64.646, 60.569, 67.323, 55.187, 51.058, 52.571, 51.476, 51.859, 53.282, 42.054}; // pol5
    // fit range
    double fit_range_16_1[10] = {65.952, 62.810, 68.275, 56.425, 52.319, 53.081, 53.042, 53.027, 54.822, 44.000}; // [0.99, 1.15]
    double fit_range_16_2[10] = {65.198, 60.806, 67.513, 55.729, 51.692, 52.899, 52.555, 53.053, 53.517, 44.292}; // [0.99, 1.25]
    // Binning
    double binning_16_1[10] = {66.263, 60.908, 67.190, 54.888, 50.938, 51.760, 51.170, 51.447, 52.509, 42.146}; // 90
    double binning_16_2[10] = {66.318, 59.802, 67.213, 55.245, 50.811, 52.615, 51.405, 51.341, 52.952, 41.534}; //110
    // beam bunches
    double beambun_16_1[10] = {65.145, 59.249, 65.158, 54.479, 50.700, 52.053, 51.083, 51.917, 52.087, 41.724}; // w2
    double beambun_16_2[10] = {64.130, 60.387, 65.834, 55.679, 51.152, 51.893, 51.766, 51.950, 52.997, 42.083}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_16_opt[10] = {65.560, 60.516, 67.353, 55.112, 51.077, 52.439, 51.574, 51.948, 53.267, 41.978};
    // abs(p_dttof)<0.3
    double p_dttof_16_1[10] = {64.445, 60.386, 66.368, 55.239, 50.409, 52.240, 50.913, 52.024, 53.078, 41.877}; // 0.2
    double p_dttof_16_2[10] = {66.350, 61.990, 68.878, 55.696, 51.591, 52.412, 52.226, 52.248, 53.319, 42.076}; // 0.4
    // kin_chisq<55
    double kin_chisq_16_1[10] = {64.205, 58.134, 65.079, 54.818, 49.332, 52.058, 50.769, 52.114, 52.706, 41.166}; // 45
    double kin_chisq_16_2[10] = {66.863, 65.160, 68.865, 56.007, 53.042, 52.294, 51.097, 52.188, 55.383, 42.875}; // 65
    // abs(mm2)<0.035
    double mm2_16_1[10] = {65.012, 60.390, 66.922, 54.181, 50.728, 51.157, 51.508, 51.602, 52.331, 41.908}; // 0.025
    double mm2_16_2[10] = {66.057, 59.691, 67.591, 54.900, 51.139, 52.391, 51.789, 52.324, 52.729, 41.974}; // 0.045

    // ------------------------ 2017

    double optimal_val_17[10] = {96.135, 89.688, 83.519, 79.727, 76.388, 78.125, 70.159, 69.726, 65.149, 59.284};
    // Bkg pol order
    double bkg_poln_17_1[10] = {95.813, 89.577, 83.439, 79.644, 76.322, 78.058, 70.065, 69.669, 65.441, 59.221}; // pol3
    double bkg_poln_17_2[10] = {96.048, 89.813, 83.833, 79.822, 76.464, 78.520, 70.265, 69.797, 65.595, 59.358}; // pol5
    // fit range
    double fit_range_17_1[10] = {97.671, 91.397, 85.678, 81.391, 77.725, 79.967, 72.654, 71.967, 68.072, 61.629}; // [0.99, 1.15]
    double fit_range_17_2[10] = {96.518, 90.369, 84.258, 79.968, 76.473, 78.594, 71.305, 70.412, 66.462, 60.204}; // [0.99, 1.25]
    // Binning
    double binning_17_1[10] = {95.419, 89.667, 83.742, 79.808, 76.324, 78.250, 70.044, 69.705, 65.489, 59.596}; // 90
    double binning_17_2[10] = {95.917, 89.695, 83.791, 79.732, 76.026, 78.152, 70.036, 69.741, 65.523, 59.142}; // 110
    // beam bunches
    double beambun_17_1[10] = {94.688, 89.267, 83.301, 79.000, 75.825, 78.878, 69.469, 69.795, 65.826, 58.672}; // w2
    double beambun_17_2[10] = {95.876, 89.233, 83.310, 79.717, 76.400, 78.138, 70.003, 69.552, 65.707, 58.740}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_17_opt[10] = {96.010, 90.003, 83.619, 79.843, 76.474, 78.137, 70.294, 69.669, 65.669, 59.639};
    // abs(p_dttof)<0.3
    double p_dttof_17_1[10] = {94.516, 88.845, 82.451, 78.544, 75.338, 77.215, 69.382, 68.767, 64.911, 58.808}; // 0.2
    double p_dttof_17_2[10] = {97.070, 91.101, 84.440, 80.499, 77.027, 78.496, 70.471, 70.005, 66.100, 60.447}; // 0.4
    // kin_chisq<55
    double kin_chisq_17_1[10] = {97.571, 91.465, 86.158, 82.221, 77.952, 80.484, 71.565, 71.174, 65.826, 59.888}; // 45 ?
    double kin_chisq_17_2[10] = {94.317, 88.945, 82.149, 78.055, 75.226, 76.392, 68.982, 67.837, 64.864, 59.169}; // 65
    // abs(mm2)<0.035
    double mm2_17_1[10] = {96.254, 90.121, 83.294, 79.408, 76.488, 77.821, 69.541, 69.043, 65.246, 59.511}; // 0.025
    double mm2_17_2[10] = {96.593, 90.391, 83.845, 80.100, 76.263, 78.451, 70.465, 69.811, 65.804, 59.759}; // 0.045

    // ------------------------ 2018 Spring

    double optimal_val_18[10] = {72.196, 65.149, 64.714, 57.988, 55.502, 54.502, 50.931, 47.396, 44.682, 43.801};
    // Bkg pol order
    double bkg_poln_18_1[10] = {72.005, 65.041, 64.599, 57.758, 55.506, 54.402, 50.806, 47.278, 44.714, 43.679}; // pol3
    double bkg_poln_18_2[10] = {72.272, 65.264, 64.838, 58.114, 55.720, 54.610, 51.174, 47.622, 44.926, 44.080}; // pol5
    // fit range
    double fit_range_18_1[10] = {71.940, 65.380, 65.467, 57.894, 55.938, 54.914, 51.951, 48.2895, 45.665, 45.802}; // [0.99, 1.15]
    double fit_range_18_2[10] = {72.798, 66.157, 65.625, 58.914, 56.059, 55.296, 52.258, 48.780, 45.956, 44.999};  // [0.99, 1.25]
    // Binning
    double binning_18_1[10] = {72.407, 65.592, 64.750, 57.822, 55.325, 54.497, 51.002, 47.335, 44.787, 43.814}; // 90
    double binning_18_2[10] = {72.006, 65.307, 64.653, 57.759, 55.332, 54.589, 50.918, 47.438, 44.728, 43.925}; // 110
    // beam bunches
    double beambun_18_1[10] = {72.101, 64.412, 65.002, 57.393, 54.857, 54.423, 51.029, 46.949, 44.727, 43.636}; // w2
    double beambun_18_2[10] = {72.309, 64.823, 64.958, 57.644, 55.339, 54.400, 50.947, 47.389, 44.677, 43.561}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18_opt[10] = {71.491, 64.379, 63.695, 57.062, 54.579, 53.632, 50.281, 46.672, 44.042, 43.370};
    // abs(p_dttof)<0.3
    double p_dttof_18_1[10] = {70.647, 63.698, 62.901, 56.485, 53.968, 53.107, 49.874, 46.289, 43.803, 42.922}; // 0.2
    double p_dttof_18_2[10] = {71.950, 64.656, 64.002, 57.415, 54.784, 53.947, 50.538, 46.947, 44.322, 43.371}; // 0.4
    // kin_chisq<55
    double kin_chisq_18_1[10] = {69.350, 63.756, 63.002, 55.925, 53.362, 52.083, 49.167, 44.966, 41.836, 41.531}; // 45
    double kin_chisq_18_2[10] = {71.617, 65.000, 64.428, 58.068, 55.761, 54.940, 51.634, 48.263, 45.574, 44.172}; // 65
    // abs(mm2)<0.035
    double mm2_18_1[10] = {71.667, 64.178, 63.955, 57.236, 54.678, 53.554, 50.764, 46.800, 44.079, 43.406}; // 0.025
    double mm2_18_2[10] = {71.659, 64.456, 63.551, 57.198, 54.490, 53.525, 50.506, 46.601, 44.055, 43.235}; // 0.045

    // ------------------------ 2018 Fall

    double optimal_val_18l[10] = {70.883, 70.100, 67.510, 64.285, 63.715, 62.693, 59.462, 57.908, 54.524, 51.789};
    // Bkg pol order
    double bkg_poln_18l_1[10] = {70.787, 70.027, 67.370, 64.148, 63.602, 62.612, 59.324, 57.773, 54.339, 51.696}; // pol3
    double bkg_poln_18l_2[10] = {70.985, 70.179, 67.610, 64.395, 63.813, 62.781, 59.499, 58.006, 54.624, 51.931}; // pol5
    // fit range
    double fit_range_18l_1[10] = {70.518, 69.814, 67.940, 64.770, 64.439, 64.344, 60.987, 59.629, 56.179, 54.095}; // [0.99, 1.15]
    double fit_range_18l_2[10] = {71.426, 70.511, 68.049, 64.885, 64.248, 63.349, 59.981, 58.907, 55.463, 53.208}; // [0.99, 1.25]
    // Binning
    double binning_18l_1[10] = {70.889, 70.343, 67.713, 64.300, 63.663, 62.823, 59.369, 57.963, 54.550, 51.965}; // 90
    double binning_18l_2[10] = {70.900, 70.457, 67.565, 64.149, 63.691, 62.625, 59.340, 57.806, 54.590, 51.718}; // 110
    // beam bunches
    double beambun_18l_1[10] = {71.323, 69.612, 67.093, 63.451, 63.068, 62.479, 59.005, 57.739, 54.445, 51.443}; // w2
    double beambun_18l_2[10] = {71.083, 69.692, 67.373, 63.906, 63.460, 62.523, 59.356, 57.715, 54.468, 51.632}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18l_opt[10] = {79.764, 79.120, 73.855, 69.565, 67.630, 65.608, 61.267, 60.338, 55.094, 52.495};
    // abs(p_dttof)<0.3
    double p_dttof_18l_1[10] = {78.132, 78.109, 72.884, 68.628, 66.916, 64.946, 60.674, 59.569, 54.624, 52.020}; // 0.2
    double p_dttof_18l_2[10] = {80.197, 79.698, 74.303, 69.995, 67.967, 65.934, 61.773, 60.681, 55.272, 52.741}; // 0.4
    // kin_chisq<55
    double kin_chisq_18l_1[10] = {81.632, 80.182, 74.952, 70.180, 68.234, 65.768, 61.437, 60.311, 54.751, 52.234}; // 45
    double kin_chisq_18l_2[10] = {78.897, 78.113, 72.903, 69.187, 67.471, 65.440, 61.606, 60.090, 55.555, 52.610}; // 65
    // abs(mm2)<0.035
    double mm2_18l_1[10] = {80.036, 79.208, 73.670, 69.434, 67.580, 65.503, 61.237, 59.864, 54.917, 51.979}; // 0.025
    double mm2_18l_2[10] = {79.681, 79.113, 73.854, 69.663, 67.662, 65.703, 61.461, 60.307, 55.177, 52.623}; // 0.045

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
    // test << std::setprecision(2) << std::fixed;
    // test << "i-1 = " << i - 1 << " | dxsec_phi2pi_stat = " << dxsec_phi2pi_stat << " | dxsec_phi2pi_sys = " << dxsec_phi2pi_sys <<" | dxsec_phi2pi_tot = " << dxsec_phi2pi_tot <<  " | error = " << gphiexsec->GetErrorY(i - 1) << endl;

    // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_phie_mc<<"$\\pm$"<<dN_phie_mc<<"&"<<N_phie_data<<"$\\pm$"<<dN_phie_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi2pi<<"$\\pm$"<<dxsec_phi2pi_stat<<"$\\pm$"<<dxsec_phi2pi_sys<<" \\\\"<< endl;
  }

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
    TLegend *ltot_tagged_flux = new TLegend(0.87, 0.80, 0.98, 0.98);
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
    h16_tagged_flux->Draw("esame");
    ctot_tagged_flux->Modified();
    ltot_tagged_flux->AddEntry(h16_tagged_flux, "2016", "p");
    ltot_tagged_flux->AddEntry(h17_tagged_flux, "2017", "p");
    ltot_tagged_flux->AddEntry(h18_tagged_flux, "2018 S", "p");
    ltot_tagged_flux->AddEntry(h18l_tagged_flux, "2018 F", "p");
    ltot_tagged_flux->Draw();
    ctot_tagged_flux->Print("output/fig_xsec_phi2pi/ctot_tagged_flux.eps", "eps");
    ctot_tagged_flux->Print("output/fig_xsec_phi2pi/ctot_tagged_flux.png", "png");
    ctot_tagged_flux->Print("output/fig_xsec_phi2pi/ctot_tagged_flux.root", "root");

    // *********** total Thrown Beam energy
    TCanvas *ctot_beame_tru = new TCanvas("ctot_beame_tru", "ctot_beame_tru", 900, 600);
    TLegend *ltot_beame_tru = new TLegend(0.87, 0.80, 0.98, 0.98);
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
    h17_beame_tru->Draw("e");
    h18l_beame_tru->Draw("esame");
    h18_beame_tru->Draw("esame");
    h16_beame_tru->Draw("esame");
    ctot_beame_tru->Modified();
    ltot_beame_tru->AddEntry(h16_beame_tru, "2016", "p");
    ltot_beame_tru->AddEntry(h17_beame_tru, "2017", "p");
    ltot_beame_tru->AddEntry(h18_beame_tru, "2018 S", "p");
    ltot_beame_tru->AddEntry(h18l_beame_tru, "2018 F", "p");
    ltot_beame_tru->Draw();
    ctot_beame_tru->Print("output/fig_xsec_phi2pi/ctot_beame_tru.eps", "eps");
    ctot_beame_tru->Print("output/fig_xsec_phi2pi/ctot_beame_tru.png", "png");
    ctot_beame_tru->Print("output/fig_xsec_phi2pi/ctot_beame_tru.root", "root");

    // *********** total Yields MC
    TMultiGraph *mgphie_mc = new TMultiGraph();
    TCanvas *cmgphie_mc = new TCanvas("cmgphie_mc", "cmgphie_mc", 900, 600);
    TLegend *lmgphie_mc = new TLegend(0.87, 0.80, 0.98, 0.98);
    lmgphie_mc->SetTextSize(0.06);
    lmgphie_mc->SetBorderSize(0);
    cmgphie_mc->cd();
    cmgphie_mc->SetGrid();
    TGraphErrors *h16_gphie_mc = (TGraphErrors *)outputfig->Get("h16_gphie_mc");
    cout << " ***** h16_gphie_mc = " << h16_gphie_mc << endl;
    // h16_gphie_mc->SetMarkerColor(kBlue);
    h16_gphie_mc->SetMarkerSize(1.5);
    h16_gphie_mc->SetMarkerStyle(kFullSquare);
    mgphie_mc->Add(h16_gphie_mc);
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
    lmgphie_mc->AddEntry(h16_gphie_mc, "2016", "p");
    lmgphie_mc->AddEntry(h17_gphie_mc, "2017", "p");
    lmgphie_mc->AddEntry(h18_gphie_mc, "2018 S", "p");
    lmgphie_mc->AddEntry(h18l_gphie_mc, "2018 F", "p");
    lmgphie_mc->Draw();
    cmgphie_mc->Print("output/fig_xsec_phi2pi/cmgphie_mc.eps", "eps");
    cmgphie_mc->Print("output/fig_xsec_phi2pi/cmgphie_mc.png", "png");
    cmgphie_mc->Print("output/fig_xsec_phi2pi/cmgphie_mc.root", "root");
    
    // *********** total Yields Data
    TMultiGraph *mgphie = new TMultiGraph();
    TCanvas *cmgphie = new TCanvas("cmgphie", "cmgphie", 900, 600);
    TLegend *lmgphie = new TLegend(0.87, 0.80, 0.98, 0.98);
    lmgphie->SetTextSize(0.06);
    lmgphie->SetBorderSize(0);
    cmgphie->cd();
    cmgphie->SetGrid();
    TGraphErrors *h16_gphie = (TGraphErrors *)outputfig->Get("h16_gphie");
    cout << " ***** h16_gphie = " << h16_gphie << endl;
    // h16_gphie->SetMarkerColor(kBlue);
    h16_gphie->SetMarkerSize(1.5);
    h16_gphie->SetMarkerStyle(kFullSquare);
    mgphie->Add(h16_gphie);
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
    lmgphie->AddEntry(h16_gphie, "2016", "p");
    lmgphie->AddEntry(h17_gphie, "2017", "p");
    lmgphie->AddEntry(h18_gphie, "2018 S", "p");
    lmgphie->AddEntry(h18l_gphie, "2018 F", "p");
    lmgphie->Draw();
    cmgphie->Print("output/fig_xsec_phi2pi/cmgphie.eps", "eps");
    cmgphie->Print("output/fig_xsec_phi2pi/cmgphie.png", "png");
    cmgphie->Print("output/fig_xsec_phi2pi/cmgphie.root", "root");

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
    h16_gphieeff->Draw("AP");
    h17_gphieeff->Draw("PSAME");
    h18_gphieeff->Draw("PSAME");
    h18l_gphieeff->Draw("PSAME");

    gPad->Modified();
    gPad->Update();
    cmgeeff->Modified();
    lmgeeff->AddEntry(h16_gphieeff, "2016", "p");
    lmgeeff->AddEntry(h17_gphieeff, "2017", "p");
    lmgeeff->AddEntry(h18_gphieeff, "2018 S", "p");
    lmgeeff->AddEntry(h18l_gphieeff, "2018 F", "p");
    lmgeeff->Draw();
  
    padeeff2->cd(); // pad2 becomes the current pad
    TGraphErrors *gr_efferate_17 = new TGraphErrors();
    TGraphErrors *gr_efferate_18 = new TGraphErrors();
    TGraphErrors *gr_efferate_18l = new TGraphErrors();
    gr_efferate_17->SetMarkerColor(kBlue);
    gr_efferate_17->SetMarkerStyle(2);
    gr_efferate_17->SetMarkerSize(2);
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

      double efferate_17 = (h16_gphieeff->GetY()[i - 1] - h17_gphieeff->GetY()[i - 1]) / (h16_gphieeff->GetY()[i - 1]);
      double efferate_18 = (h16_gphieeff->GetY()[i - 1] - h18_gphieeff->GetY()[i - 1]) / (h16_gphieeff->GetY()[i - 1]);
      double efferate_18l = (h16_gphieeff->GetY()[i - 1] - h18l_gphieeff->GetY()[i - 1]) / (h16_gphieeff->GetY()[i - 1]);
      cout << "efferate_17 = " << efferate_17 << " | efferate_18 = " << efferate_18 << " | efferate_18l = " << efferate_18l << endl;
      gr_efferate_17->SetPoint(i-1, Eg, efferate_17);//i-1, h2phie_mc->GetXaxis()->GetBinCenter(i), efferate_17
      gr_efferate_17->SetPointError(i-1, 0, 0);
      gr_efferate_18->SetPoint(i-1, Eg, efferate_18);
      gr_efferate_18->SetPointError(i-1, 0, 0);
      gr_efferate_18l->SetPoint(i-1, Eg, efferate_18l);
      gr_efferate_18l->SetPointError(i-1, 0, 0);
    }
    
    gr_efferate_17->Draw("AP");
    gr_efferate_18->Draw("PSAME");
    gr_efferate_18l->Draw("PSAME");
    TLine leffe;
    leffe.SetLineStyle(kDashed);
    leffe.DrawLine(Eg_min, 0, Eg_max, 0);

    // Y axis ratio plot settings
    gr_efferate_17->GetHistogram()->SetStats(0);
    gr_efferate_17->SetTitle(";E_{#gamma} (GeV);Ratio");
    // gr_efferate_17->GetYaxis()->SetNdivisions(506);
    gr_efferate_17->GetYaxis()->SetTitleSize(35);
    gr_efferate_17->GetYaxis()->SetTitleFont(43);
    gr_efferate_17->GetYaxis()->SetTitleOffset(1.55);
    gr_efferate_17->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efferate_17->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_efferate_17->GetXaxis()->SetTitleSize(38);
    gr_efferate_17->GetXaxis()->SetTitleFont(43);
    gr_efferate_17->GetXaxis()->SetTitleOffset(3.);
    gr_efferate_17->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efferate_17->GetXaxis()->SetLabelSize(38);
    gr_efferate_17->GetHistogram()->SetMaximum(1.0);   // along          
    gr_efferate_17->GetHistogram()->SetMinimum(-1.0);  //   Y 
    gPad->Update();
    gPad->Modified();
    cmgeeff->Print("output/fig_xsec_phi2pi/cmgeeff.eps", "eps");
    cmgeeff->Print("output/fig_xsec_phi2pi/cmgeeff.png", "png");
    cmgeeff->Print("output/fig_xsec_phi2pi/cmgeeff.root", "root");

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
    h17_gphiexsec->Draw("AP");
    h16_gphiexsec->Draw("PSAME");
    h18_gphiexsec->Draw("PSAME");
    h18l_gphiexsec->Draw("PSAME");
    gPad->Modified();
    gPad->Update();
    lmgexsec->AddEntry(h16_gphiexsec, "2016", "p");
    lmgexsec->AddEntry(h17_gphiexsec, "2017", "p");
    lmgexsec->AddEntry(h18_gphiexsec, "2018 S", "p");
    lmgexsec->AddEntry(h18l_gphiexsec, "2018 F", "p");
    lmgexsec->Draw();
 
    padexsec2->cd(); // pad2 becomes the current pad
    TGraphErrors *gr_exsecrate_17 = new TGraphErrors();
    TGraphErrors *gr_exsecrate_18 = new TGraphErrors();
    TGraphErrors *gr_exsecrate_18l = new TGraphErrors();
    gr_exsecrate_17->SetMarkerColor(kBlue);
    gr_exsecrate_17->SetMarkerStyle(2);//20
    gr_exsecrate_17->SetMarkerSize(2);
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

      double exsecrate_17 = (h16_gphiexsec->GetY()[i-1] - h17_gphiexsec->GetY()[i-1])/(h16_gphiexsec->GetY()[i-1]);
      double exsecrate_18 = (h16_gphiexsec->GetY()[i-1] - h18_gphiexsec->GetY()[i-1])/(h16_gphiexsec->GetY()[i-1]);
      double exsecrate_18l = (h16_gphiexsec->GetY()[i-1] - h18l_gphiexsec->GetY()[i-1])/(h16_gphiexsec->GetY()[i-1]);
      cout<<"exsecrate_17 = "<<exsecrate_17<<" | exsecrate_18 = "<<exsecrate_18<<" | exsecrate_18l = "<<exsecrate_18l<<endl;
      gr_exsecrate_17->SetPoint(i-1, Eg, exsecrate_17);//h2phie->GetXaxis()->GetBinCenter(i)
      gr_exsecrate_17->SetPointError(i-1, 0, 0);
      gr_exsecrate_18->SetPoint(i-1, Eg, exsecrate_18);
      gr_exsecrate_18->SetPointError(i-1, 0, 0);
      gr_exsecrate_18l->SetPoint(i-1, Eg, exsecrate_18l);
      gr_exsecrate_18l->SetPointError(i-1, 0, 0);
    }

    gr_exsecrate_17->Draw("AP");
    gr_exsecrate_18->Draw("PSAME");
    gr_exsecrate_18l->Draw("PSAME");
    TLine lexsec;
    lexsec.SetLineStyle(kDashed);
    lexsec.DrawLine(Eg_min, 0, Eg_max, 0);

    // Y axis ratio plot settings
    gr_exsecrate_17->GetHistogram()->SetStats(0);
    gr_exsecrate_17->SetTitle(";E_{#gamma} (GeV);Ratio");
    // gr_exsecrate_17->GetYaxis()->SetNdivisions(506);
    gr_exsecrate_17->GetYaxis()->SetTitleSize(35);
    gr_exsecrate_17->GetYaxis()->SetTitleFont(43);
    gr_exsecrate_17->GetYaxis()->SetTitleOffset(1.55);
    gr_exsecrate_17->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_exsecrate_17->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_exsecrate_17->GetXaxis()->SetTitleSize(38);
    gr_exsecrate_17->GetXaxis()->SetTitleFont(43);
    gr_exsecrate_17->GetXaxis()->SetTitleOffset(3.);
    gr_exsecrate_17->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_exsecrate_17->GetXaxis()->SetLabelSize(38);

    gr_exsecrate_17->GetHistogram()->SetMaximum(1.0);   // along          
    gr_exsecrate_17->GetHistogram()->SetMinimum(-1.0);  //   Y  
    gPad->Update();
    gPad->Modified();
    cmgexsec->Print("output/fig_xsec_phi2pi/cmgexsec.eps", "eps");
    cmgexsec->Print("output/fig_xsec_phi2pi/cmgexsec.png", "png");
    cmgexsec->Print("output/fig_xsec_phi2pi/cmgexsec.root", "root");

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
    TGraphErrors *gphiexsec_17corr = new TGraphErrors();
    gphiexsec_17corr->SetMarkerStyle(kOpenCircle);
    gphiexsec_17corr->SetMarkerSize(1.5);
    gphiexsec_17corr->SetMarkerColor(kBlue);
    gphiexsec_17corr->SetLineColor(kBlue);
    gphiexsec_17corr->SetMinimum(0.0);
    gphiexsec_17corr->SetTitle("; E_{#gamma} (GeV); #sigma (nb)");

    // vector<double> xi = {1.0,0.8,1.1}, dxi = {0.1,0.2,0.05};
    // pair<double, double> x = pdgavg(xi, dxi);
    // Eg = 0;
    vector<double> xsec(3), dxsec(3);
    pair<double, double> xsec_pair;

    for (int i = 1; i <= ne; ++i)
    {
      xsec = {h16_gphiexsec->GetY()[i - 1], h18_gphiexsec->GetY()[i - 1], h18l_gphiexsec->GetY()[i - 1]};
      dxsec = {h16_gphiexsec->GetErrorY(i - 1), h18_gphiexsec->GetErrorY(i - 1), h18l_gphiexsec->GetErrorY(i - 1)};

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

      double w_17 = 0.76; //0.76
      // if(i ==5)
      // {
      //   w_17 = xsec_pair.first/h17_gphiexsec->GetY()[i-1];
      // } 
      
      gphiexsec_17corr->SetPoint(i-1, Eg, w_17*h17_gphiexsec->GetY()[i - 1]);
      gphiexsec_17corr->SetPointError(i-1, 0, w_17*h17_gphiexsec->GetErrorY(i - 1));

      test << std::setprecision(2) << std::fixed;
      test << "i-1 = " << i - 1 << " | xsec[0] = " << xsec[0] << " | dxsec[0] = " << dxsec[0] << " | xsec[1] = " << xsec[1] << " | dxsec[1] = " << dxsec[1] <<  " | xsec[2] = " << xsec[2] << " | dxsec[2] = " << dxsec[2] << " | xsec_pair.first = " << xsec_pair.first << " | xsec_pair.second = " << xsec_pair.second <<  " | w_17 = " << w_17 << endl;
    }

    h17_gphiexsec->Draw("AP");
    gphiexsec_avg->Draw("PSAME");
    gphiexsec_17corr->Draw("PSAME");
    lgexsec_avg->AddEntry(h17_gphiexsec, "2017", "p");
    lgexsec_avg->AddEntry(gphiexsec_avg, "#bar{#sigma}_{(2016, 2018 S, 2018 F)}", "p");
    lgexsec_avg->AddEntry(gphiexsec_17corr, "2017 corrected", "p");
    lgexsec_avg->Draw();

    cgexsec_avg->Print("output/fig_xsec_phi2pi/cgexsec_avg.eps", "eps");
    cgexsec_avg->Print("output/fig_xsec_phi2pi/cgexsec_avg.png", "png");
    cgexsec_avg->Print("output/fig_xsec_phi2pi/cgexsec_avg.root", "root");

    cphie->Print(Form("output/fig_xsec_phi2pi/c2_phie_%s.root", name.Data()), "root");
    cphie->Print(Form("output/fig_xsec_phi2pi/c2_phie_%s.eps", name.Data()), "eps");
    cphie->Print(Form("output/fig_xsec_phi2pi/c2_phie_%s.png", name.Data()), "png");
    cphie1->Print(Form("output/fig_xsec_phi2pi/c_phie1_%s.root", name.Data()), "root");
    cphie1->Print(Form("output/fig_xsec_phi2pi/c_phie1_%s.eps", name.Data()), "eps");
    cphie1->Print(Form("output/fig_xsec_phi2pi/c_phie1_%s.png", name.Data()), "png");
    cgphie->Print(Form("output/fig_xsec_phi2pi/c_gphie_%s.root", name.Data()), "root");
    cgphie->Print(Form("output/fig_xsec_phi2pi/c_gphie_%s.eps", name.Data()), "eps");
    cgphie->Print(Form("output/fig_xsec_phi2pi/c_gphie_%s.png", name.Data()), "png");
    cphie_mc->Print(Form("output/fig_xsec_phi2pi/c2_phie_mc_%s.root", name.Data()), "root");
    cphie_mc->Print(Form("output/fig_xsec_phi2pi/c2_phie_mc_%s.eps", name.Data()), "eps");
    cphie_mc->Print(Form("output/fig_xsec_phi2pi/c2_phie_mc_%s.png", name.Data()),"png");
    cphie1_mc->Print(Form("output/fig_xsec_phi2pi/c_phie1_mc_%s.root", name.Data()),"root");
    cphie1_mc->Print(Form("output/fig_xsec_phi2pi/c_phie1_mc_%s.eps", name.Data()),"eps");
    cphie1_mc->Print(Form("output/fig_xsec_phi2pi/c_phie1_mc_%s.png", name.Data()),"png");
    cgphie_mc->Print(Form("output/fig_xsec_phi2pi/c_gphie_mc_%s.root", name.Data()),"root");
    cgphie_mc->Print(Form("output/fig_xsec_phi2pi/c_gphie_mc_%s.eps", name.Data()),"eps");
    cgphie_mc->Print(Form("output/fig_xsec_phi2pi/c_gphie_mc_%s.png", name.Data()),"png");
    cgphieeff->Print(Form("output/fig_xsec_phi2pi/c_gphieeff_%s.root", name.Data()),"root");
    cgphieeff->Print(Form("output/fig_xsec_phi2pi/c_gphieeff_%s.eps", name.Data()),"eps");
    cgphieeff->Print(Form("output/fig_xsec_phi2pi/c_gphieeff_%s.png", name.Data()),"png");
    cgphiexsec->Print(Form("output/fig_xsec_phi2pi/c_gphiexsec_%s.root", name.Data()), "root");
    cgphiexsec->Print(Form("output/fig_xsec_phi2pi/c_gphiexsec_%s.eps", name.Data()), "eps");
    cgphiexsec->Print(Form("output/fig_xsec_phi2pi/c_gphiexsec_%s.png", name.Data()), "png");
    c_beame_tru->Print(Form("output/fig_xsec_phi2pi/c_beame_tru_%s.root", name.Data()),"root");
    c_beame_tru->Print(Form("output/fig_xsec_phi2pi/c_beame_tru_%s.eps", name.Data()),"eps");
    c_beame_tru->Print(Form("output/fig_xsec_phi2pi/c_beame_tru_%s.png", name.Data()),"png");
    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.png", name.Data()), "png");
    
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
/*
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
  cout << "h_tagged_flux" << h_tagged_flux << endl;
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

    double optimal_val_16[10] = {14.775, 17.004, 11.186, 7.667, 4.621, 2.800, 1.779, 1.249, 0.625, 0.557};
    // Bkg pol order
    double bkg_poln_16_1[10] = {14.750, 16.981, 11.153, 7.654, 4.716, 2.912, 1.864, 1.247, 0.691, 0.617}; // pol3
    double bkg_poln_16_2[10] = {14.803, 17.029, 11.204, 7.681, 4.626, 2.810, 1.780, 1.208, 0.630, 0.559}; // pol5
    // fit range
    double fit_range_16_1[10] = {15.373, 17.368, 11.474, 7.904, 4.718, 2.879, 1.818, 1.240, 0.662, 0.622}; //[0.99,1.15]
    double fit_range_16_2[10] = {14.984, 17.259, 11.470, 7.865, 4.656, 2.831, 1.797, 1.223, 0.627, 0.601}; //[0.99,1.25]
    // Binning
    double binning_16_1[10] = {14.757, 16.976, 11.174, 7.677, 4.630, 2.780, 1.771, 1.262, 0.653, 0.562}; // 90
    double binning_16_2[10] = {14.700, 16.985, 11.196, 7.711, 4.623, 2.808, 1.797, 1.252, 0.669, 0.585}; //110
    // beam bunches
    double beambun_16_1[10] = {14.331, 16.917, 11.192, 7.625, 4.646, 2.828, 1.767, 1.201, 0.596, 0.577}; // w2
    double beambun_16_2[10] = {14.611, 17.091, 11.158, 7.704, 4.615, 2.830, 1.785, 1.226, 0.616, 0.564}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_16_opt[10] = {14.819, 17.022, 11.194, 7.676, 4.633, 2.811, 1.801, 1.267, 0.625, 0.575};
    // abs(p_dttof)<0.3
    double p_dttof_16_1[10] = {14.842, 16.910, 11.007, 7.574, 4.569, 2.742, 1.771, 1.245, 0.635, 0.549}; // 0.2
    double p_dttof_16_2[10] = {14.865, 17.117, 11.322, 7.728, 4.638, 2.864, 1.833, 1.297, 0.665, 0.605}; // 0.4
    // kin_chisq<55
    double kin_chisq_16_1[10] = {14.447, 16.742, 10.727, 7.671, 4.453, 2.749, 1.792, 1.205, 0.668, 0.606}; // 45
    double kin_chisq_16_2[10] = {15.106, 17.323, 11.430, 7.792, 4.721, 2.945, 1.942, 1.278, 0.643, 0.707}; // 65
    // abs(mm2)<0.035
    double mm2_16_1[10] = {14.636, 16.993, 11.110, 7.603, 4.606, 2.799, 1.834, 1.249, 0.631, 0.597}; // 0.025
    double mm2_16_2[10] = {14.749, 17.071, 11.193, 7.657, 4.658, 2.800, 1.784, 1.272, 0.608, 0.580}; // 0.045

    // ------------------------ 2017

    double optimal_val_17[10] = {14.935, 19.221, 15.559, 11.793, 7.674, 5.377, 3.823, 2.762, 1.907, 1.509};
    // Bkg pol order
    double bkg_poln_17_1[10] = {14.958, 19.202, 15.545, 11.829, 7.707, 5.404, 3.884, 2.848, 1.944, 1.566}; // pol3
    double bkg_poln_17_2[10] = {15.012, 19.243, 15.573, 11.845, 7.721, 5.410, 3.845, 2.765, 1.908, 1.533}; // pol5
    // fit range
    double fit_range_17_1[10] = {15.536, 19.612, 15.896, 12.055, 7.928, 5.465, 3.868, 2.873, 1.931, 1.545};//[0.99,1.15]
    double fit_range_17_2[10] = {15.116, 19.365, 15.656, 11.888, 7.778, 5.425, 3.790, 2.785, 1.898, 1.531};//[0.99,1.25]
    // Binning
    double binning_17_1[10] = {14.967, 19.218, 15.551, 11.842, 7.679, 5.400, 3.799, 2.789, 1.914, 1.523}; // 90
    double binning_17_2[10] = {14.989, 19.215, 15.528, 11.852, 7.689, 5.377, 3.815, 2.778, 1.900, 1.530}; // 110
    // beam bunches
    double beambun_17_1[10] = {14.905, 19.107, 15.528, 11.789, 7.660, 5.382, 3.826, 2.718, 1.905, 1.531}; // w2
    double beambun_17_2[10] = {14.926, 19.174, 15.562, 11.843, 7.694, 5.338, 3.829, 2.757, 1.932, 1.548}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_17_opt[10] = {14.894, 19.236, 15.582, 11.834, 7.713, 5.441, 3.870, 2.820, 1.935, 1.558};
    // abs(p_dttof)<0.3
    double p_dttof_17_1[10] = {14.852, 19.034, 15.347, 11.696, 7.571, 5.281, 3.735, 2.763, 1.811, 1.517}; // 0.2
    double p_dttof_17_2[10] = {14.919, 19.404, 15.732, 11.944, 7.823, 5.501, 3.942, 2.876, 1.992, 1.596}; // 0.4
    // kin_chisq<55
    double kin_chisq_17_1[10] = {14.816, 19.585, 16.064, 12.337, 8.097, 5.673, 3.963, 2.939, 2.072, 1.674}; // 45
    double kin_chisq_17_2[10] = {15.040, 18.974, 15.227, 11.454, 7.497, 5.153, 3.789, 2.663, 1.856, 1.545}; // 65
    // abs(mm2)<0.035
    double mm2_17_1[10] = {14.738, 19.222, 15.604, 11.892, 7.764, 5.419, 3.944, 2.866, 1.976, 1.609}; // 0.025
    double mm2_17_2[10] = {14.939, 19.223, 15.576, 11.835, 7.696, 5.407, 3.858, 2.781, 1.938, 1.542}; // 0.045

    // ------------------------ 2018 Spring

    double optimal_val_18[10] = {12.969, 14.777, 10.490, 7.294, 4.784, 3.274, 2.106, 1.612, 1.099, 0.746};
    // Bkg pol order
    double bkg_poln_18_1[10] = {12.907, 14.750, 10.454, 7.283, 4.778, 3.270, 2.101, 1.610, 1.083, 0.745}; // pol3
    double bkg_poln_18_2[10] = {13.006, 14.804, 10.507, 7.306, 4.791, 3.279, 2.109, 1.619, 1.102, 0.749}; // pol5
    // fit range
    double fit_range_18_1[10] = {13.325, 14.756, 10.470, 7.250, 4.898, 3.352, 2.159, 1.682, 1.137, 0.818}; //[0.99,1.15]
    double fit_range_18_2[10] = {13.246, 14.973, 10.632, 7.413, 4.859, 3.306, 2.139, 1.662, 1.112, 0.794}; //[0.99,1.25]
    // Binning
    double binning_18_1[10] = {12.954, 14.759, 10.472, 7.287, 4.776, 3.262, 2.092, 1.626, 1.099, 0.765}; // 90
    double binning_18_2[10] = {12.958, 14.749, 10.491, 7.284, 4.782, 3.265, 2.104, 1.623, 1.093, 0.741}; // 110
    // beam bunches
    double beambun_18_1[10] = {12.845, 14.726, 10.397, 7.279, 4.763, 3.250, 2.064, 1.582, 1.111, 0.712}; // w2
    double beambun_18_2[10] = {12.916, 14.751, 10.459, 7.286, 4.768, 3.269, 2.091, 1.611, 1.079, 0.718}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18_opt[10] = {12.737, 14.559, 10.322, 7.196, 4.719, 3.260, 2.094, 1.597, 1.102, 0.741};
    // abs(p_dttof)<0.3
    double p_dttof_18_1[10] = {12.744, 14.483, 10.251, 7.119, 4.652, 3.194, 2.052, 1.547, 1.070, 0.710}; // 0.2
    double p_dttof_18_2[10] = {12.754, 14.633, 10.404, 7.245, 4.754, 3.276, 2.119, 1.610, 1.112, 0.758}; // 0.4
    // kin_chisq<55
    double kin_chisq_18_1[10] = {12.239, 14.167, 10.150, 7.089, 4.656, 3.193, 2.061, 1.576, 1.150, 0.750}; // 45
    double kin_chisq_18_2[10] = {13.115, 14.907, 10.509, 7.282, 4.776, 3.252, 2.103, 1.548, 1.116, 0.742}; // 65
    // abs(mm2)<0.035
    double mm2_18_1[10] = {12.780, 14.606, 10.385, 7.237, 4.763, 3.281, 2.115, 1.603, 1.113, 0.746}; // 0.025
    double mm2_18_2[10] = {12.762, 14.553, 10.301, 7.197, 4.693, 3.251, 2.091, 1.600, 1.109, 0.747}; // 0.045

    // ------------------------ 2018 Fall

    double optimal_val_18l[10] = {13.960, 16.520, 11.536, 7.833, 4.996, 3.106, 2.015, 1.467, 1.014, 0.771};
    // Bkg pol order
    double bkg_poln_18l_1[10] = {13.936, 16.488, 11.525, 7.815, 4.991, 3.104, 2.013, 1.465, 1.011, 0.771}; // pol3
    double bkg_poln_18l_2[10] = {14.006, 16.543, 11.558, 7.832, 5.007, 3.109, 2.018, 1.469, 1.015, 0.773}; // pol5
    // fit range
    double fit_range_18l_1[10] = {14.470, 16.652, 11.543, 7.833, 4.986, 3.173, 2.007, 1.511, 1.048, 0.793};//[0.99,1.15]
    double fit_range_18l_2[10] = {14.189, 16.694, 11.656, 7.913, 5.039, 3.115, 2.027, 1.482, 1.017, 0.769};//[0.99,1.25]
    // Binning
    double binning_18l_1[10] = {13.996, 16.523, 11.545, 7.838, 5.000, 3.108, 2.018, 1.465, 1.014, 0.762}; // 90
    double binning_18l_2[10] = {13.977, 16.499, 11.533, 7.832, 4.998, 3.106, 2.018, 1.463, 1.020, 0.763}; // 110
    // beam bunches
    double beambun_18l_1[10] = {13.895, 16.396, 11.462, 7.813, 4.975, 3.093, 2.019, 1.466, 1.023, 0.792}; // w2
    double beambun_18l_2[10] = {13.952, 16.494, 11.483, 7.827, 4.983, 3.085, 2.014, 1.466, 1.018, 0.777}; // w4
    // cut variations (cut at ttree level not DSelector)
    double cutsvar_18l_opt[10] = {15.266, 17.779, 11.824, 7.737, 4.979, 3.109, 2.137, 1.530, 1.052, 0.792};
    // abs(p_dttof)<0.3
    double p_dttof_18l_1[10] = {15.275, 17.652, 11.713, 7.635, 4.905, 3.049, 2.091, 1.487, 1.010, 0.750}; // 0.2
    double p_dttof_18l_2[10] = {15.304, 17.870, 11.889, 7.785, 5.006, 3.137, 2.151, 1.552, 1.058, 0.796}; // 0.4
    // kin_chisq<55
    double kin_chisq_18l_1[10] = {15.108, 17.808, 11.987, 7.865, 5.043, 3.181, 2.198, 1.571, 1.100, 0.827}; // 45
    double kin_chisq_18l_2[10] = {15.510, 17.801, 11.750, 7.668, 4.904, 3.052, 2.128, 1.485, 1.026, 0.773}; // 65
    // abs(mm2)<0.035
    double mm2_18l_1[10] = {15.193, 17.752, 11.818, 7.742, 4.978, 3.129, 2.153, 1.537, 1.061, 0.785}; // 0.025
    double mm2_18l_2[10] = {15.287, 17.796, 11.836, 7.730, 4.980, 3.113, 2.137, 1.532, 1.052, 0.790}; // 0.045

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

    test << std::setprecision(3) << std::fixed;
    test << "i-1 = " << i - 1 << " | xsec_phi2pi = " << xsec_phi2pi << endl;
    // " | dxsec_phi2pi_sys = " << dxsec_phi2pi_sys <<" | dxsec_phi2pi_tot = " << dxsec_phi2pi_tot <<  " | error = " << gphiexsec->GetErrorY(i - 1) << endl;
    }

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
    TLegend *ltot_t_tru = new TLegend(0.87, 0.80, 0.98, 0.98);
    ltot_t_tru->SetTextSize(0.06);
    ltot_t_tru->SetBorderSize(0);
    ctot_t_tru->cd();
    ctot_t_tru->SetGrid();
    TH1F *h16_t_tru = (TH1F *)outputfig->Get("h16_t_tru");
    cout << " ***** h16_t_tru = " << h16_t_tru << endl;
    // h16_t_tru->SetMarkerColor(kBlue);
    h16_t_tru->SetMarkerSize(1.5);
    h16_t_tru->SetMarkerStyle(kFullSquare);
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
    h16_t_tru->SetTitle(";E_{#gamma}:Counts");
    h16_t_tru->SetMinimum(0.);
    h17_t_tru->Draw("e");
    h18_t_tru->Draw("esame");
    h18l_t_tru->Draw("esame");
    h16_t_tru->Draw("esame");
    ctot_t_tru->Modified();
    ltot_t_tru->AddEntry(h16_t_tru, "2016", "p");
    ltot_t_tru->AddEntry(h17_t_tru, "2017", "p");
    ltot_t_tru->AddEntry(h18_t_tru, "2018 S", "p");
    ltot_t_tru->AddEntry(h18l_t_tru, "2018 F", "p");
    ltot_t_tru->Draw();
    ctot_t_tru->Print("output/fig_xsec_phi2pi/ctot_t_tru.eps", "eps");
    ctot_t_tru->Print("output/fig_xsec_phi2pi/ctot_t_tru.png", "png");
    ctot_t_tru->Print("output/fig_xsec_phi2pi/ctot_t_tru.root", "root");

    // total Yields MC
    TMultiGraph *mgphit_mc = new TMultiGraph();
    TCanvas *cmgphit_mc = new TCanvas("cmgphit_mc", "cmgphit_mc", 900, 600);
    TLegend *lmgphit_mc = new TLegend(0.87, 0.80, 0.98, 0.98);
    lmgphit_mc->SetTextSize(0.06);
    lmgphit_mc->SetBorderSize(0);
    cmgphit_mc->cd();
    cmgphit_mc->SetGrid();
    TGraphErrors *h16_gphit_mc = (TGraphErrors *)outputfig->Get("h16_gphit_mc");
    cout << " ***** h16_gphit_mc = " << h16_gphit_mc << endl;
    // h16_gphit_mc->SetMarkerColor(kBlue);
    h16_gphit_mc->SetMarkerSize(1.5);
    h16_gphit_mc->SetMarkerStyle(kFullSquare);
    mgphit_mc->Add(h16_gphit_mc);
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
    lmgphit_mc->AddEntry(h16_gphit_mc, "2016", "p");
    lmgphit_mc->AddEntry(h17_gphit_mc, "2017", "p");
    lmgphit_mc->AddEntry(h18_gphit_mc, "2018 S", "p");
    lmgphit_mc->AddEntry(h18l_gphit_mc, "2018 F", "p");
    lmgphit_mc->Draw();
    cmgphit_mc->Print("output/fig_xsec_phi2pi/cmgphit_mc.eps", "eps");
    cmgphit_mc->Print("output/fig_xsec_phi2pi/cmgphit_mc.png", "png");
    cmgphit_mc->Print("output/fig_xsec_phi2pi/cmgphit_mc.root", "root");

    // total Yields Data
    TMultiGraph *mgphit = new TMultiGraph();
    TCanvas *cmgphit = new TCanvas("cmgphit", "cmgphit", 900, 600);
    TLegend *lmgphit = new TLegend(0.87, 0.80, 0.98, 0.98);
    lmgphit->SetTextSize(0.06);
    lmgphit->SetBorderSize(0);
    cmgphit->cd();
    cmgphit->SetGrid();
    TGraphErrors *h16_gphit = (TGraphErrors *)outputfig->Get("h16_gphit");
    cout << " ***** h16_gphit = " << h16_gphit << endl;
    // h16_gphit->SetMarkerColor(kBlue);
    h16_gphit->SetMarkerSize(1.5);
    h16_gphit->SetMarkerStyle(kFullSquare);
    mgphit->Add(h16_gphit);
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
    lmgphit->AddEntry(h16_gphit, "2016", "p");
    lmgphit->AddEntry(h17_gphit, "2017", "p");
    lmgphit->AddEntry(h18_gphit, "2018 S", "p");
    lmgphit->AddEntry(h18l_gphit, "2018 F", "p");
    lmgphit->Draw();
    cmgphit->Print("output/fig_xsec_phi2pi/cmgphit.eps", "eps");
    cmgphit->Print("output/fig_xsec_phi2pi/cmgphit.png", "png");
    cmgphit->Print("output/fig_xsec_phi2pi/cmgphit.root", "root");

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
    TGraphErrors *h16_gphiteff = (TGraphErrors *)outputfig->Get("h16_gphiteff");
    cout << " ***** h16_gphiteff = " << h16_gphiteff << endl;
    // h16_gphiteff->SetMarkerColor(kBlue);
    h16_gphiteff->SetMarkerSize(1.5);
    h16_gphiteff->SetMarkerStyle(kFullSquare);
    h16_gphiteff->SetTitle(";-t (GeV/c)^{2}; #varepsilon (%)");
    h16_gphiteff->SetMinimum(0.);
    TGraphErrors *h17_gphiteff = (TGraphErrors *)outputfig->Get("h17_gphiteff");
    cout << " ***** h17_gphiteff = " << h17_gphiteff << endl;
    h17_gphiteff->SetMarkerColor(kBlue);
    h17_gphiteff->SetMarkerSize(1.5);
    h17_gphiteff->SetMarkerStyle(kFullCircle);
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
    h16_gphiteff->Draw("AP");
    h17_gphiteff->Draw("PSAME");
    h18_gphiteff->Draw("PSAME");
    h18l_gphiteff->Draw("PSAME");
    cmgteff->Modified();
    lmgteff->AddEntry(h16_gphiteff, "2016", "p");
    lmgteff->AddEntry(h17_gphiteff, "2017", "p");
    lmgteff->AddEntry(h18_gphiteff, "2018 S", "p");
    lmgteff->AddEntry(h18l_gphiteff, "2018 F", "p");
    lmgteff->Draw();

    padteff2->cd(); // pad2 becomes the current pad
    TGraphErrors *gr_efftrate_17 = new TGraphErrors();
    TGraphErrors *gr_efftrate_18 = new TGraphErrors();
    TGraphErrors *gr_efftrate_18l = new TGraphErrors();
    gr_efftrate_17->SetMarkerColor(kBlue);
    gr_efftrate_17->SetMarkerStyle(2);
    gr_efftrate_17->SetMarkerSize(2); 
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
      
      double efftrate_17 = (h16_gphiteff->GetY()[i-1] - h17_gphiteff->GetY()[i-1])/(h16_gphiteff->GetY()[i-1]);
      double efftrate_18 = (h16_gphiteff->GetY()[i-1] - h18_gphiteff->GetY()[i-1])/(h16_gphiteff->GetY()[i-1]);
      double efftrate_18l = (h16_gphiteff->GetY()[i-1] - h18l_gphiteff->GetY()[i-1])/(h16_gphiteff->GetY()[i-1]);
      cout<<"efftrate_17 = "<<efftrate_17<<" | efftrate_18 = "<<efftrate_18<<" | efftrate_18l = "<<efftrate_18l<<endl;
      gr_efftrate_17->SetPoint(i-1, t, abs(efftrate_17));
      gr_efftrate_17->SetPointError(i-1, 0, 0);
      gr_efftrate_18->SetPoint(i-1, t, abs(efftrate_18));
      gr_efftrate_18->SetPointError(i-1, 0, 0);
      gr_efftrate_18l->SetPoint(i-1, t, abs(efftrate_18l));
      gr_efftrate_18l->SetPointError(i-1, 0, 0);
    }

    gr_efftrate_17->Draw("AP");
    gr_efftrate_18->Draw("PSAME");
    gr_efftrate_18l->Draw("PSAME");
    TLine lefft;
    lefft.SetLineStyle(kDashed);
    lefft.DrawLine(t_min, 0, t_max, 0);

    // Y axis ratio plot settings
    gr_efftrate_17->GetHistogram()->SetStats(0);
    gr_efftrate_17->SetTitle(";-t (GeV/c)^{2};Ratio");
    gr_efftrate_17->GetYaxis()->SetNdivisions(506);
    gr_efftrate_17->GetYaxis()->SetTitleSize(35);
    gr_efftrate_17->GetYaxis()->SetTitleFont(43);
    gr_efftrate_17->GetYaxis()->SetTitleOffset(1.55);
    gr_efftrate_17->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efftrate_17->GetYaxis()->SetLabelSize(35);
    // X axis ratio plot settings
    gr_efftrate_17->GetXaxis()->SetTitleSize(38);
    gr_efftrate_17->GetXaxis()->SetTitleFont(43);
    gr_efftrate_17->GetXaxis()->SetTitleOffset(3.);
    gr_efftrate_17->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_efftrate_17->GetXaxis()->SetLabelSize(38);
    gr_efftrate_17->GetHistogram()->SetMaximum(1.0);   // along          
    gr_efftrate_17->GetHistogram()->SetMinimum(-1.0);  //   Y  
    gPad->Update();
    gPad->Modified();    
    cmgteff->Print("output/fig_xsec_phi2pi/cmgteff.eps", "eps");
    cmgteff->Print("output/fig_xsec_phi2pi/cmgteff.png", "png");
    cmgteff->Print("output/fig_xsec_phi2pi/cmgteff.root", "root");

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
    TGraphErrors *h16_gphitxsec = (TGraphErrors *)outputfig->Get("h16_gphitxsec");
    cout << " ***** h16_gphitxsec = " << h16_gphitxsec << endl;
    h16_gphitxsec->SetMarkerStyle(kFullSquare);
    h16_gphitxsec->SetMarkerSize(1.5);
    h16_gphitxsec->SetMinimum(0.);
    h16_gphitxsec->SetTitle(";-t (GeV/c)^{2}; #sigma (nb)");
    // h16_gphitxsec->SetMarkerColor(kBlue);
    // mgtxsec->Add(h16_gphitxsec);
    TGraphErrors *h17_gphitxsec = (TGraphErrors *)outputfig->Get("h17_gphitxsec");
    cout << " ***** h17_gphitxsec = " << h17_gphitxsec << endl;
    h17_gphitxsec->SetMarkerStyle(kFullCircle);
    h17_gphitxsec->SetMarkerSize(1.5);
    h17_gphitxsec->SetMarkerColor(kBlue);
    h17_gphitxsec->SetLineColor(kBlue);
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
    h17_gphitxsec->Draw("AP");
    h16_gphitxsec->Draw("PSAME");
    h18_gphitxsec->Draw("PSAME");
    h18l_gphitxsec->Draw("PSAME");
    cmgtxsec->Modified();
    lmgtxsec->AddEntry(h16_gphitxsec, "2016", "p");
    lmgtxsec->AddEntry(h17_gphitxsec, "2017", "p");
    lmgtxsec->AddEntry(h18_gphitxsec, "2018 S", "p");
    lmgtxsec->AddEntry(h18l_gphitxsec, "2018 F", "p");
    lmgtxsec->Draw();

    padtxsec2->cd(); // pad2 becomes the current pad
    TGraphErrors *gr_txsecrate_17 = new TGraphErrors();
    TGraphErrors *gr_txsecrate_18 = new TGraphErrors();
    TGraphErrors *gr_txsecrate_18l = new TGraphErrors();
    gr_txsecrate_17->SetMarkerColor(kBlue);
    gr_txsecrate_17->SetMarkerStyle(2);
    gr_txsecrate_17->SetMarkerSize(2); 
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

      double txsecrate_17 = (h16_gphitxsec->GetY()[i-1] - h17_gphitxsec->GetY()[i-1])/(h16_gphitxsec->GetY()[i-1]);
      double txsecrate_18 = (h16_gphitxsec->GetY()[i-1] - h18_gphitxsec->GetY()[i-1])/(h16_gphitxsec->GetY()[i-1]);
      double txsecrate_18l = (h16_gphitxsec->GetY()[i-1] - h18l_gphitxsec->GetY()[i-1])/(h16_gphitxsec->GetY()[i-1]);
      cout<<"txsecrate_17 = "<<txsecrate_17<<" | txsecrate_18 = "<<txsecrate_18<<" | txsecrate_18l = "<<txsecrate_18l<<endl;
      gr_txsecrate_17->SetPoint(i-1, t, abs(txsecrate_17));
      gr_txsecrate_17->SetPointError(i-1, 0, 0);
      gr_txsecrate_18->SetPoint(i-1, t, abs(txsecrate_18));
      gr_txsecrate_18->SetPointError(i-1, 0, 0);
      gr_txsecrate_18l->SetPoint(i-1, t, abs(txsecrate_18l));
      gr_txsecrate_18l->SetPointError(i-1, 0, 0);
    }

    gr_txsecrate_17->Draw("AP");
    gr_txsecrate_18->Draw("PSAME");
    gr_txsecrate_18l->Draw("PSAME");
    TLine ltxsec;
    ltxsec.SetLineStyle(kDashed);
    ltxsec.DrawLine(t_min, 0, t_max, 0);

    // Y axis ratio plot settings
    gr_txsecrate_17->GetHistogram()->SetStats(0);
    gr_txsecrate_17->SetTitle(";-t (GeV/c)^{2};Ratio");
    gr_txsecrate_17->GetYaxis()->SetNdivisions(506);
    gr_txsecrate_17->GetYaxis()->SetTitleSize(35);//24
    gr_txsecrate_17->GetYaxis()->SetTitleFont(43);
    gr_txsecrate_17->GetYaxis()->SetTitleOffset(1.55);//1.55
    gr_txsecrate_17->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_txsecrate_17->GetYaxis()->SetLabelSize(35);//30
    // X axis ratio plot settings
    gr_txsecrate_17->GetXaxis()->SetTitleSize(38);//35
    gr_txsecrate_17->GetXaxis()->SetTitleFont(43);
    gr_txsecrate_17->GetXaxis()->SetTitleOffset(3.0); //3.
    gr_txsecrate_17->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    gr_txsecrate_17->GetXaxis()->SetLabelSize(38);//35
    
    gr_txsecrate_17->GetHistogram()->SetMaximum(1.0);   // along          
    gr_txsecrate_17->GetHistogram()->SetMinimum(-1.0);  //   Y  
    gPad->Update();
    gPad->Modified();
    cmgtxsec->Update();
    cmgtxsec->Print("output/fig_xsec_phi2pi/cmgtxsec.eps", "eps");
    cmgtxsec->Print("output/fig_xsec_phi2pi/cmgtxsec.png", "png");
    cmgtxsec->Print("output/fig_xsec_phi2pi/cmgtxsec.root", "root");

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
    gphitxsec_avg->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. E_{#gamma}; E_{#gamma} (GeV); #sigma (nb)");
    TGraphErrors *gphitxsec_17corr = new TGraphErrors();
    gphitxsec_17corr->SetMarkerStyle(kOpenCircle);
    gphitxsec_17corr->SetMarkerSize(1.5);
    gphitxsec_17corr->SetMarkerColor(kBlue);
    gphitxsec_17corr->SetLineColor(kBlue);
    gphitxsec_17corr->SetMinimum(0.0);
    gphitxsec_17corr->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. E_{#gamma}; E_{#gamma} (GeV); #sigma (nb)");

    // vector<double> xi = {1.0,0.8,1.1}, dxi = {0.1,0.2,0.05};
    // pair<double, double> x = pdgavg(xi, dxi);
    // Eg = 0;
    vector<double> xsec(3), dxsec(3);
    pair<double, double> xsec_pair;

    for (int i = 1; i <= nt; ++i)
    {
      xsec = {h16_gphitxsec->GetY()[i - 1], h18_gphitxsec->GetY()[i - 1], h18l_gphitxsec->GetY()[i - 1]};
      dxsec = {h16_gphitxsec->GetErrorY(i - 1), h18_gphitxsec->GetErrorY(i - 1), h18l_gphitxsec->GetErrorY(i - 1)};

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
      
      gphitxsec_17corr->SetPoint(i-1, t, w_17*h17_gphitxsec->GetY()[i - 1]);
      gphitxsec_17corr->SetPointError(i-1, 0, w_17*h17_gphitxsec->GetErrorY(i - 1));

      // test << std::setprecision(2) << std::fixed;
      // test << "i-1 = " << i - 1 << " | xsec[0] = " << xsec[0] << " | dxsec[0] = " << dxsec[0] << " | xsec[1] = " << xsec[1] << " | dxsec[1] = " << dxsec[1] <<  " | xsec[2] = " << xsec[2] << " | dxsec[2] = " << dxsec[2] << " | xsec_pair.first = " << xsec_pair.first << " | xsec_pair.second = " << xsec_pair.second <<  " | w_17 = " << w_17 << endl;
    }

    h17_gphitxsec->Draw("AP");
    gphitxsec_avg->Draw("PSAME");
    gphitxsec_17corr->Draw("PSAME");
    lgtxsec_avg->AddEntry(h17_gphitxsec, "2017", "p");
    lgtxsec_avg->AddEntry(gphitxsec_avg, "#bar{#sigma}_{(2016, 2018 S, 2018 F)}", "p");
    lgtxsec_avg->AddEntry(gphitxsec_17corr, "2017 corrected", "p");
    lgtxsec_avg->Draw();

    cgtxsec_avg->Print("output/fig_xsec_phi2pi/cgtxsec_avg.eps", "eps");
    cgtxsec_avg->Print("output/fig_xsec_phi2pi/cgtxsec_avg.png", "png");
    cgtxsec_avg->Print("output/fig_xsec_phi2pi/cgtxsec_avg.root", "root");

    cphit->Print(Form("output/fig_xsec_phi2pi/c2_phit_%s.root", name.Data()), "root");
    cphit->Print(Form("output/fig_xsec_phi2pi/c2_phit_%s.eps", name.Data()), "eps");
    cphit->Print(Form("output/fig_xsec_phi2pi/c2_phit_%s.png", name.Data()), "png");
    cphit1->Print(Form("output/fig_xsec_phi2pi/c_phit1_%s.root", name.Data()), "root");
    cphit1->Print(Form("output/fig_xsec_phi2pi/c_phit1_%s.eps", name.Data()), "eps");
    cphit1->Print(Form("output/fig_xsec_phi2pi/c_phit1_%s.png", name.Data()), "png");
    cgphit->Print(Form("output/fig_xsec_phi2pi/c_gphit_%s.root", name.Data()), "root");
    cgphit->Print(Form("output/fig_xsec_phi2pi/c_gphit_%s.eps", name.Data()), "eps");
    cgphit->Print(Form("output/fig_xsec_phi2pi/c_gphit_%s.png", name.Data()), "png");
    cphit_mc->Print(Form("output/fig_xsec_phi2pi/c2_phit_mc_%s.root", name.Data()), "root");
    cphit_mc->Print(Form("output/fig_xsec_phi2pi/c2_phit_mc_%s.eps", name.Data()), "eps");
    cphit_mc->Print(Form("output/fig_xsec_phi2pi/c2_phit_mc_%s.png", name.Data()), "png");
    cphit1_mc->Print(Form("output/fig_xsec_phi2pi/c_phit1_mc_%s.root", name.Data()), "root");
    cphit1_mc->Print(Form("output/fig_xsec_phi2pi/c_phit1_mc_%s.eps", name.Data()), "eps");
    cphit1_mc->Print(Form("output/fig_xsec_phi2pi/c_phit1_mc_%s.png", name.Data()), "png");
    cgphit_mc->Print(Form("output/fig_xsec_phi2pi/c_gphit_mc_%s.root", name.Data()), "root");
    cgphit_mc->Print(Form("output/fig_xsec_phi2pi/c_gphit_mc_%s.eps", name.Data()), "eps");
    cgphit_mc->Print(Form("output/fig_xsec_phi2pi/c_gphit_mc_%s.png", name.Data()), "png");
    cgphiteff->Print(Form("output/fig_xsec_phi2pi/c_gphiteff_%s.root", name.Data()), "root");
    cgphiteff->Print(Form("output/fig_xsec_phi2pi/c_gphiteff_%s.eps", name.Data()), "eps");
    cgphiteff->Print(Form("output/fig_xsec_phi2pi/c_gphiteff_%s.png", name.Data()), "png");
    cgphitxsec->Print(Form("output/fig_xsec_phi2pi/c_gphitxsec_%s.root", name.Data()), "root");
    cgphitxsec->Print(Form("output/fig_xsec_phi2pi/c_gphitxsec_%s.eps", name.Data()), "eps");
    cgphitxsec->Print(Form("output/fig_xsec_phi2pi/c_gphitxsec_%s.png", name.Data()), "png");

    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.root", name.Data()), "root");
    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("output/fig_xsec_phi2pi/c_tagged_flux_%s.png", name.Data()), "png");
    c_h_t_Thrown->Print(Form("output/fig_xsec_phi2pi/c_h_t_Thrown_%s.root", name.Data()), "root");
    c_h_t_Thrown->Print(Form("output/fig_xsec_phi2pi/c_h_t_Thrown_%s.eps", name.Data()), "eps");
    c_h_t_Thrown->Print(Form("output/fig_xsec_phi2pi/c_h_t_Thrown_%s.png", name.Data()), "png");


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
*/
}
