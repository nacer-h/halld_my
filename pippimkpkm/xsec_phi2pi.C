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

void xsec_phi2pi(TString name, int n2k=100, int ne=10, int nt=10)//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  // TFile *fdata = NULL;
  TFile *fdata = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_%s_flat.root", name.Data()));
  TFile *fmc = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_phi2pi_%s_flat.root", name.Data()));
  TFile *ftru = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/thrown_phi2pi_%s.root", name.Data()));
  // TFile *fthr = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/tree_thrown_phi2pi_genr8_%s.root", name.Data()));
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TTree *tthr = (TTree *)fthr->Get("Thrown_Tree");
  TFile *fps = NULL;
  if(name == "16") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_11366_11555.root");
  if(name == "17") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  if(name == "18") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_40856_42577.root");
  if(name == "18l") fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_50677_51768.root");

  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/xsec_phi2pi.root", "UPDATE");
  
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

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.98, mkk_max = 1.2; //mkk_min = 0.98, mkk_max = 1.2;
  double Eg_min = 6.5, Eg_max = 11.6;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  // ======================================== Phi vs. Eg ===============================================

  //root -l 'xsec_phi2pi.C+("data_17",100,10,1)'

  FILE *table_xsec_phi2pi = fopen(Form("table_xsec_phi2pi_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections in photon energy bins}\n \\begin{tabular}{|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma$ [nb] \\\\ \n \\hline\n");

  FILE *table_xsec_phi2pi_sys = fopen(Form("table_xsec_phi2pi_sys_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi_sys,"\\documentclass[12pt]{extarticle}\n \\usepackage{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabular}{|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & polynomial degrees [ \\% ] & Fit range [ \\% ]  \\\\ \n \\hline\n");

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux = " << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.5);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = (TH1F *)ftru->Get("h_beame_Thrown"); // new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, Eg_min, Eg_max)
  // tthr->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru = " << h_beame_tru << endl;
  h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.5);
  h_beame_tru->Draw("e");

  // Data
  TCanvas *cphie = new TCanvas("cphie", "cphie", 900, 600);
  cphie->cd();
  TH2D *h2phie = new TH2D("h2phie", "Data; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, Eg_min, Eg_max, n2k, mkk_min, mkk_max);
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");

  h2phie->Draw("colz");

  TCanvas *cphie1 = new TCanvas("cphie1", "cphie1", 1500, 600);
  cphie1->Divide(5, 2);
  TCanvas *cgphie = new TCanvas("cgphie", "cgphie", 900, 600);
  TGraphErrors *gphie = new TGraphErrors();
  gphie->SetMarkerStyle(20);
  gphie->SetMarkerSize(1.5);
  gphie->SetMarkerColor(1);
  gphie->SetMinimum(0.0);
  gphie->SetTitle("#phi(1020) Yield vs. E_{#gamma} (data); E_{#gamma} [GeV]; N_{#phi}");

  // Monte Carlo
  TCanvas *cphie_mc = new TCanvas("cphie_mc", "cphie_mc", 900, 600);
  cphie_mc->cd();
  TH2D *h2phie_mc = new TH2D("h2phie_mc", "MC; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, Eg_min, Eg_max, n2k, mkk_min, mkk_max);//0.98, 1.2
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");
  h2phie_mc->Draw("colz");

  TCanvas *cphie1_mc = new TCanvas("cphie1_mc", "cphie1_mc", 1500, 600);
  cphie1_mc->Divide(5, 2);
  TCanvas *cgphie_mc = new TCanvas("cgphie_mc", "cgphie_mc", 900, 600);
  TGraphErrors *gphie_mc = new TGraphErrors();
  gphie_mc->SetMarkerStyle(20);
  gphie_mc->SetMarkerSize(1.5);
  gphie_mc->SetMarkerColor(1);
  gphie_mc->SetMinimum(0.0);
  gphie_mc->SetTitle("#phi(1020) Yield vs. E_{#gamma} (MC); E_{#gamma} [GeV]; N_{#phi}");

  // Efficiency
  TCanvas *cgphieeff = new TCanvas("cgphieeff", "cgphieeff", 900, 600);
  TGraphErrors *gphieeff = new TGraphErrors();
  gphieeff->SetMarkerStyle(20);
  gphieeff->SetMarkerSize(1.5);
  gphieeff->SetMarkerColor(kBlue);
  gphieeff->SetMinimum(0.0);
  gphieeff->SetTitle("#phi(1020) Effeciency vs. E_{#gamma}; E_{#gamma} [GeV]; #epsilon [%]");

  // Cross-section
  TCanvas *cgphiexsec = new TCanvas("cgphiexsec", "cgphiexsec", 900, 600);
  TGraphErrors *gphiexsec = new TGraphErrors();
  gphiexsec->SetMarkerStyle(20);
  gphiexsec->SetMarkerSize(1.5);
  gphiexsec->SetMarkerColor(1);
  gphiexsec->SetMinimum(0.0);
  gphiexsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. E_{#gamma}; E_{#gamma} [GeV]; #sigma [nb]"); //#phi(1020) flux normalized yield, Yield_{#phi}

  for (int i = 1; i <= ne; ++i)
  {
    cout << i << " " << flush;

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphie1_mc->cd(i);
    TH1D *hphie_mc_py = h2phie_mc->ProjectionY(Form("hphie_mc_py_%d", i), i, i);
    hphie_mc_py->SetTitle(Form("E_{#gamma} = %.2f [GeV]",h2phie_mc->GetXaxis()->GetBinCenter(i)));
    hphie_mc_py->Draw("e");

    // w.factory(Form("Voigtian::sig_phie_mc(m_phie_mc[%f,%f],mean_phie_mc[1.015,1.022],width_phie_mc[0.004],sigma_phie_mc[0.002])", mkk_min, mkk_max)); //sigma_phie_mc[0.0001,0.01], mean_phie_mc[1.011,1.030]
    // // w.factory(Form("BreitWigner::sig_phie_mc(m_phie_mc[%f,%f],mean_phie_mc[1.018,1.021],width_phie_mc[0.004])", mkk_min, mkk_max));
    // w.factory("Chebychev::bkg_phie_mc(m_phie_mc,{c0_phie_mc[-10,10], c1_phie_mc[-10,10], c2_phie_mc[-10,10]})");
    // w.factory("SUM:::model_phie_mc(nsig_phie_mc[0,100000000]*sig_phie_mc, nbkg_phie_mc[0,100000000]*bkg_phie_mc)"); //nsig[0,100000000]*sig2,
    // w.var("m_phie_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    // RooDataHist dh_phie_mc("dh_phie_mc", "dh_phie_mc", *w.var("m_phie_mc"), Import(*hphie_mc_py));
    // RooPlot *fr_phie_mc = w.var("m_phie_mc")->frame(Title(Form("E_{#gamma} = %f",h2phie_mc->GetXaxis()->GetBinCenter(i))));
    // w.pdf("model_phie_mc")->fitTo(dh_phie_mc);

    // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    // dh_phie_mc.plotOn(fr_phie_mc, RooFit::Name("ldh_phie_mc"));
    // w.pdf("model_phie_mc")->plotOn(fr_phie_mc, Components(*w.pdf("sig_phie_mc")), LineColor(kRed), RooFit::Name("lsig_phie_mc"));
    // w.pdf("model_phie_mc")->plotOn(fr_phie_mc, Components(*w.pdf("bkg_phie_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phie_mc"));
    // w.pdf("model_phie_mc")->plotOn(fr_phie_mc, RooFit::Name("lmodel_phie_mc"));
    // // w.pdf("model_phie_mc")->paramOn(fr_phie_mc, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phie_mc"), *w.var("nbkg_phie_mc")))); //,*w.var("mean_phie_mc"),*w.var("width_phie_mc"),*w.var("sigma_phie_mc"))));
    // fr_phie_mc->Draw();

    // TLegend *l_phie_mc = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_phie_mc->SetFillColor(kWhite);
    // l_phie_mc->AddEntry(fr_phie_mc->findObject("lsig_phie_mc"), Form("N_{Sig} = %.2f", w.var("nsig_phie_mc")->getVal()), "l");
    // l_phie_mc->AddEntry(fr_phie_mc->findObject("lbkg_phie_mc"), Form("N_{Bkg} = %.2f", w.var("nbkg_phie_mc")->getVal()), "l");
    // l_phie_mc->Draw();

    // double N_phie_mc = w.var("nsig_phie_mc")->getVal();
    // double dN_phie_mc = w.var("nsig_phie_mc")->getError();

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
    lat_phie_mc.SetTextSize(0.05);
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

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphie1->cd(i);
    TH1D *hphie_py = h2phie->ProjectionY(Form("hphie_py_%d", i), i, i);
    hphie_py->SetTitle(Form("E_{#gamma} = %.2f [GeV]",h2phie->GetXaxis()->GetBinCenter(i)));
    hphie_py->Draw("e");

    // w.factory(Form("Voigtian::sig_phie(m_phie[%f,%f],mean_phie[1.015, 1.022],width_phie[0.004],sigma_phie[0.002])", mkk_min, mkk_max)); //sigma_phie[0.0001,0.01], mean_phie[1.016,1.022 or 1.010,1.030]
    // // w.factory(Form("BreitWigner::sig_phie(m_phie[%f,%f],mean_phie[1.018,1.021],width_phie[0.004])", mkk_min, mkk_max));
    // w.factory("Chebychev::bkg_phie(m_phie,{c0_phie[-10,10], c1_phie[-10,10], c2_phie[-10,10]})");
    // w.factory("SUM:::model_phie(nsig_phie[0,100000000]*sig_phie, nbkg_phie[0,100000000]*bkg_phie)"); //nsig[0,100000000]*sig2,
    // w.var("m_phie")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    // RooDataHist dh_phie("dh_phie", "dh_phie", *w.var("m_phie"), Import(*hphie_py));
    // RooPlot *fr_phie = w.var("m_phie")->frame(Title(Form("E_{#gamma} = %f",h2phie->GetXaxis()->GetBinCenter(i))));
    // w.pdf("model_phie")->fitTo(dh_phie);
    // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    // dh_phie.plotOn(fr_phie, RooFit::Name("ldh_phie"));
    // w.pdf("model_phie")->plotOn(fr_phie, Components(*w.pdf("sig_phie")), LineColor(kRed), RooFit::Name("lsig_phie"));
    // w.pdf("model_phie")->plotOn(fr_phie, Components(*w.pdf("bkg_phie")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phie"));
    // w.pdf("model_phie")->plotOn(fr_phie, RooFit::Name("lmodel_phie"));
    // // w.pdf("model_phie")->paramOn(fr_phie, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phie"), *w.var("nbkg_phie")))); //,*w.var("mean_phie"),*w.var("width_phie"),*w.var("sigma_phie"))));
    // fr_phie->Draw();

    // TLegend *l_phie = new TLegend(0.5, 0.7, 0.8, 0.9);
    // l_phie->SetFillColor(kWhite);
    // l_phie->SetLineColor(kWhite);
    // // l_phie->AddEntry(fr_phie->findObject("ldh_phie"), "Data", "p");
    // // l_phie->AddEntry(fr_phie->findObject("lmodel_phie"), "total", "l");
    // l_phie->AddEntry(fr_phie->findObject("lsig_phie"), Form("N_{Sig} = %.2f", w.var("nsig_phie")->getVal()), "l");
    // l_phie->AddEntry(fr_phie->findObject("lbkg_phie"), Form("N_{Bkg} = %.2f", w.var("nbkg_phie")->getVal()), "l");
    // l_phie->Draw();

    // double N_phie = w.var("nsig_phie")->getVal();
    // double dN_phie = w.var("nsig_phie")->getError();

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

    // hphie_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
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

    // hphie_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
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

    // hphie_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
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

    hphie_py->Fit("fsb_data", "", "", mkk_min, mkk_max);
    double par_data[10]; //6
    fsb_data->GetParameters(&par_data[0]);
    fs_data->SetParameters(&par_data[0]);
    fb_data->SetParameters(&par_data[4]); //4

    fs_data->Draw("same");
    fb_data->Draw("same");

    double N_phie_data = fs_data->Integral(mkk_min, mkk_max) / hphie_py->GetBinWidth(1);
    double dN_phie_data = N_phie_data * fsb_data->GetParError(0) / fsb_data->GetParameter(0);

    TLatex lat_phie_data;
    lat_phie_data.SetTextSize(0.05);
    lat_phie_data.SetTextAlign(13); //align at top
    lat_phie_data.SetNDC();
    lat_phie_data.SetTextColor(kBlue);
    lat_phie_data.DrawLatex(0.45, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phie_data.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phie_data, dN_phie_data));
    lat_phie_data.DrawLatex(0.45, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phie_data.DrawLatex(0.45, 0.66, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phie_data.DrawLatex(0.45, 0.59, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

    gphie->SetPoint(i - 1, h2phie->GetXaxis()->GetBinCenter(i), N_phie_data);
    gphie->SetPointError(i - 1, 0, dN_phie_data);    

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    // if (N_phie_mc <= 0)
    // {
    //   //gphiexsec->RemovePoint(i - 1);
    //   continue;
    // }
    double eff_phi = N_phie_mc / h_beame_tru->GetBinContent(i); // Efficiency = N_observed/N_generated
    double deff_phi = eff_phi * (dN_phie_mc / N_phie_mc);
    gphieeff->SetPoint(i-1, h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phi*100);
    gphieeff->SetPointError(i-1, 0, deff_phi*100); //->GetBinError(i)
    // gphieeff->SetPoint(gphieeff->GetN(), h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phi*100);
    // gphieeff->SetPointError(gphieeff->GetN()-1, 0, deff_phi*100); //->GetBinError(i)

    // +++++++++++++++++++ Systematic Errors 
    // double width_value[3] = {-169.74, -141.06, -197.77}; // Y Width value : 0.079, 0.065, 0.093  
    // double err_width_value = TMath::RMS(3, width_value);

    double N_bkg_poln_1[10] = {1689, 1753, 2465, 6581, 8566, 2607, 3110, 2557, 2450, 160}; // Bkg pol degree: 4
    double N_bkg_poln_2[10] = {1685, 1749, 2459, 6567, 8548, 2677, 3103, 2552, 2444, 182}; // Bkg pol degree: 3
    double N_bkg_poln_3[10] = {1694, 1756, 2471, 6595, 8583, 2611, 3118, 2562, 2456, 161}; // Bkg pol degree: 5
    double bkg_poln[3] = {N_bkg_poln_1[i-1], N_bkg_poln_2[i-1], N_bkg_poln_3[i-1]};
    double err_bkg_poln = TMath::RMS(3, bkg_poln);
    double fit_range_1[10] = {1689, 1753, 2465, 6581, 8566, 2607, 3110, 2557, 2450, 160}; // Fit range: [0.98, 1.2]
    double fit_range_2[10] = {1685, 1767, 2477, 6598, 8551, 2600, 3094, 2532, 2439, 182}; // Fit range: [0.98, 1.18]
    double fit_range_3[10] = {1693, 1736, 2458, 6609, 8522, 2582, 3104, 2540, 2433, 163}; // Fit range: [0.98, 1.22]
    double fit_range[3] = {fit_range_1[i-1], fit_range_2[i-1], fit_range_3[i-1]};
    double err_fit_range = TMath::RMS(3, fit_range);
    // double fit_func_1[10] = {1689.35, 1752.50, 2465.09, 6580.90, 8565.59, 2606.73, 3110.43, 2557.01, 2449.74, 160.43}; // Fit function: Voigt
    // double fit_func_2[10] = {2157.41, 2302.83, 3158.83, 7835.29, 10889.65, 3077.05, 3893.53, 3240.87, 3128.20, 217.33}; // Fit function: Breit-wigner
    // double fit_func_3[10] = {1608.70, 1734.22, 2367.27, 6358.65, 8255.21, 2452.65, 2940.46, 2476.66, 2322.49, 203.03}; // Fit function: Double Gaus
    // double fit_func[3] = {fit_func_1[i-1], fit_func_2[i-1], fit_func_3[i-1]};
    // double err_fit_func = TMath::RMS(3, fit_func);

    double err_sys = TMath::Sqrt(err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range);

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
    double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys / N_phie_data)*(err_sys / N_phie_data) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double xsec_phi = N_phie / lumi_phi;
    // double dxsec_phi = dN_phie / lumi_phi;  
    gphiexsec->SetPoint(i-1, h2phie->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    gphiexsec->SetPointError(i-1, 0, dxsec_phi2pi_stat);
    // gphiexsec->SetPoint(gphiexsec->GetN(), h2phie->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    // gphiexsec->SetPointError(gphiexsec->GetN() - 1, 0, dxsec_phi2pi_stat);

    fprintf(table_xsec_phi2pi, "%0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h2phie_mc->GetXaxis()->GetBinCenter(i), h_beame_tru->GetBinContent(i), N_phie_mc, dN_phie_mc, N_phie_data, dN_phie_data, eff_phi*100, deff_phi*100, xsec_phi2pi, dxsec_phi2pi_stat, dxsec_phi2pi_sys);
    
    fprintf(table_xsec_phi2pi_sys, "%0.2f & %0.2f & %0.2f \\\\ \n", h2phie_mc->GetXaxis()->GetBinCenter(i), err_bkg_poln*100/N_bkg_poln_1[i-1], err_fit_range*100/fit_range_1[i-1]);
    // table_phi << std::setprecision(2) << std::fixed;
    // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_phie_mc<<"$\\pm$"<<dN_phie_mc<<"&"<<N_phie_data<<"$\\pm$"<<dN_phie_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi2pi<<"$\\pm$"<<dxsec_phi2pi_stat<<"$\\pm$"<<dxsec_phi2pi_sys<<" \\\\"<< endl;
    }

    cgphie->cd();
    gphie->Draw("AP");
    gphie->Write(Form("h%s_gphie", name.Data()), TObject::kWriteDelete);
    cgphie_mc->cd();
    gphie_mc->Draw("AP");
    cgphieeff->cd();
    gphieeff->Draw("AP");
    gphieeff->Write(Form("h%s_gphieeff", name.Data()), TObject::kWriteDelete);
    cgphiexsec->cd();
    gphiexsec->Draw("AP");
    gphiexsec->Write(Form("h%s_gphiexsec", name.Data()), TObject::kWriteDelete);
    // int j =1;
    // gphie->Write(Form("grphie_%d", j), TObject::kWriteDelete);


    // total Yields
    TMultiGraph *mgphie = new TMultiGraph();
    TCanvas *cmgphie = new TCanvas("cmgphie", "cmgphie", 900, 600);
    TLegend *lmgphie = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgphie->SetTextSize(0.05);
    lmgphie->SetBorderSize(0);
    cmgphie->cd();
    cmgphie->SetGrid();
    TGraphErrors *h16_gphie = (TGraphErrors *)outputfig->Get("h16_gphie");
    cout << " ***** h16_gphie = " << h16_gphie << endl;
    h16_gphie->SetMarkerColor(1);
    h16_gphie->SetMarkerSize(1.5);
    h16_gphie->SetLineColor(1);
    h16_gphie->SetMarkerStyle(20);
    mgphie->Add(h16_gphie);
    TGraphErrors *h17_gphie = (TGraphErrors *)outputfig->Get("h17_gphie");
    cout << " ***** h17_gphie = " << h17_gphie << endl;
    h17_gphie->SetMarkerColor(kBlue);
    h17_gphie->SetMarkerSize(1.5);
    h17_gphie->SetLineColor(kBlue);
    h17_gphie->SetMarkerStyle(20);
    mgphie->Add(h17_gphie);
    TGraphErrors *h18_gphie = (TGraphErrors *)outputfig->Get("h18_gphie");
    cout << " ***** h18_gphie = " << h18_gphie << endl;
    h18_gphie->SetMarkerColor(kRed);
    h18_gphie->SetMarkerSize(1.5);
    h18_gphie->SetLineColor(kRed);
    h18_gphie->SetMarkerStyle(20);
    mgphie->Add(h18_gphie);
    TGraphErrors *h18l_gphie = (TGraphErrors *)outputfig->Get("h18l_gphie");
    cout << " ***** h18l_gphie = " << h18l_gphie << endl;
    h18l_gphie->SetMarkerColor(kMagenta);
    h18l_gphie->SetMarkerSize(1.5);
    h18l_gphie->SetLineColor(kMagenta);
    h18l_gphie->SetMarkerStyle(20);
    mgphie->Add(h18l_gphie);
    mgphie->SetTitle("; E_{#gamma} [GeV]; N_{#phi}");
    mgphie->Draw("AP");
    mgphie->SetMinimum(0.);
    cmgphie->Modified();
    lmgphie->AddEntry(h16_gphie, "2016", "p");
    lmgphie->AddEntry(h17_gphie, "2017", "p");
    lmgphie->AddEntry(h18_gphie, "2018 S", "p");
    lmgphie->AddEntry(h18l_gphie, "2018 F", "p");
    lmgphie->Draw();
    cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphie.eps", "eps");
    cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphie.png", "png");
    cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphie.root", "root");

    // total Efficiency
    TMultiGraph *mgeeff = new TMultiGraph();
    TCanvas *cmgeeff = new TCanvas("cmgeeff", "cmgeeff", 900, 600);
    TLegend *lmgeeff = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgeeff->SetTextSize(0.05);
    lmgeeff->SetBorderSize(0);
    cmgeeff->cd();
    cmgeeff->SetGrid();
    TGraphErrors *h16_gphieeff = (TGraphErrors *)outputfig->Get("h16_gphieeff");
    cout << " ***** h16_gphieeff = " << h16_gphieeff << endl;
    h16_gphieeff->SetMarkerColor(1);
    h16_gphieeff->SetMarkerSize(1.5);
    h16_gphieeff->SetLineColor(1);
    h16_gphieeff->SetMarkerStyle(20);
    mgeeff->Add(h16_gphieeff);
    TGraphErrors *h17_gphieeff = (TGraphErrors *)outputfig->Get("h17_gphieeff");
    cout << " ***** h17_gphieeff = " << h17_gphieeff << endl;
    h17_gphieeff->SetMarkerColor(kBlue);
    h17_gphieeff->SetMarkerSize(1.5);
    h17_gphieeff->SetLineColor(kBlue);
    h17_gphieeff->SetMarkerStyle(20);
    mgeeff->Add(h17_gphieeff);
    TGraphErrors *h18_gphieeff = (TGraphErrors *)outputfig->Get("h18_gphieeff");
    cout << " ***** h18_gphieeff = " << h18_gphieeff << endl;
    h18_gphieeff->SetMarkerColor(kRed);
    h18_gphieeff->SetMarkerSize(1.5);
    h18_gphieeff->SetLineColor(kRed);
    h18_gphieeff->SetMarkerStyle(20);
    mgeeff->Add(h18_gphieeff);
    TGraphErrors *h18l_gphieeff = (TGraphErrors *)outputfig->Get("h18l_gphieeff");
    cout << " ***** h18l_gphieeff = " << h18l_gphieeff << endl;
    h18l_gphieeff->SetMarkerColor(kMagenta);
    h18l_gphieeff->SetMarkerSize(1.5);
    h18l_gphieeff->SetLineColor(kMagenta);
    h18l_gphieeff->SetMarkerStyle(20);
    mgeeff->Add(h18l_gphieeff);
    mgeeff->SetTitle("; E_{#gamma} [GeV]; #epsilon [%]");
    mgeeff->Draw("AP");
    mgeeff->SetMinimum(0.);
    cmgeeff->Modified();
    lmgeeff->AddEntry(h16_gphieeff, "2016", "p");
    lmgeeff->AddEntry(h17_gphieeff, "2017", "p");
    lmgeeff->AddEntry(h18_gphieeff, "2018 S", "p");
    lmgeeff->AddEntry(h18l_gphieeff, "2018 F", "p");
    lmgeeff->Draw();
    cmgeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgeeff.eps", "eps");
    cmgeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgeeff.png", "png");
    cmgeeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgeeff.root", "root");

    // total Cross-sections
    TMultiGraph *mgexsec = new TMultiGraph();
    TCanvas *cmgexsec = new TCanvas("cmgexsec", "cmgexsec", 900, 600);
    TLegend *lmgexsec = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgexsec->SetTextSize(0.05);
    lmgexsec->SetBorderSize(0);
    cmgexsec->cd();
    cmgexsec->SetGrid();
    TGraphErrors *h16_gphiexsec = (TGraphErrors *)outputfig->Get("h16_gphiexsec");
    cout << " ***** h16_gphiexsec = " << h16_gphiexsec << endl;
    h16_gphiexsec->SetMarkerColor(1);
    h16_gphiexsec->SetMarkerSize(1.5);
    h16_gphiexsec->SetLineColor(1);
    h16_gphiexsec->SetMarkerStyle(20);
    mgexsec->Add(h16_gphiexsec);
    TGraphErrors *h17_gphiexsec = (TGraphErrors *)outputfig->Get("h17_gphiexsec");
    cout << " ***** h17_gphiexsec = " << h17_gphiexsec << endl;
    h17_gphiexsec->SetMarkerColor(kBlue);
    h17_gphiexsec->SetMarkerSize(1.5);
    h17_gphiexsec->SetLineColor(kBlue);
    h17_gphiexsec->SetMarkerStyle(20);
    mgexsec->Add(h17_gphiexsec);
    TGraphErrors *h18_gphiexsec = (TGraphErrors *)outputfig->Get("h18_gphiexsec");
    cout << " ***** h18_gphiexsec = " << h18_gphiexsec << endl;
    h18_gphiexsec->SetMarkerColor(kRed);
    h18_gphiexsec->SetMarkerSize(1.5);
    h18_gphiexsec->SetLineColor(kRed);
    h18_gphiexsec->SetMarkerStyle(20);
    mgexsec->Add(h18_gphiexsec);
    TGraphErrors *h18l_gphiexsec = (TGraphErrors *)outputfig->Get("h18l_gphiexsec");
    cout << " ***** h18l_gphiexsec = " << h18l_gphiexsec << endl;
    h18l_gphiexsec->SetMarkerColor(kMagenta);
    h18l_gphiexsec->SetMarkerSize(1.5);
    h18l_gphiexsec->SetLineColor(kMagenta);
    h18l_gphiexsec->SetMarkerStyle(20);
    mgexsec->Add(h18l_gphiexsec);
    mgexsec->SetTitle("; E_{#gamma} [GeV]; #sigma [nb]");
    mgexsec->Draw("AP");
    mgexsec->SetMinimum(0.);
    cmgexsec->Modified();
    lmgexsec->AddEntry(h16_gphiexsec, "2016", "p");
    lmgexsec->AddEntry(h17_gphiexsec, "2017", "p");
    lmgexsec->AddEntry(h18_gphiexsec, "2018 S", "p");
    lmgexsec->AddEntry(h18l_gphiexsec, "2018 F", "p");
    lmgexsec->Draw();
    cmgexsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgexsec.eps", "eps");
    cmgexsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgexsec.png", "png");
    cmgexsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgexsec.root", "root");

    cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie.root", name.Data()), "root");
    cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie.eps", name.Data()), "eps");
    cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie.png", name.Data()), "png");
    cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1.root", name.Data()), "root");
    cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1.eps", name.Data()), "eps");
    cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1.png", name.Data()), "png");
    cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie.root", name.Data()), "root");
    cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie.eps", name.Data()), "eps");
    cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie.png", name.Data()), "png");
    cphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie_mc.root", name.Data()), "root");
    cphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie_mc.eps", name.Data()), "eps");
    cphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie_mc.png", name.Data()),"png");
    cphie1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1_mc.root", name.Data()),"root");
    cphie1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1_mc.eps", name.Data()),"eps");
    cphie1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phie1_mc.png", name.Data()),"png");
    cgphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie_mc.root", name.Data()),"root");
    cgphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie_mc.eps", name.Data()),"eps");
    cgphie_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphie_mc.png", name.Data()),"png");
    cgphieeff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphieeff.root", name.Data()),"root");
    cgphieeff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphieeff.eps", name.Data()),"eps");
    cgphieeff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphieeff.png", name.Data()),"png");
    cgphiexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiexsec.root", name.Data()), "root");
    cgphiexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiexsec.eps", name.Data()), "eps");
    cgphiexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiexsec.png", name.Data()), "png");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_beame_tru.root", name.Data()),"root");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_beame_tru.eps", name.Data()),"eps");
    c_beame_tru->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_beame_tru.png", name.Data()),"png");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.png", name.Data()), "png");
    
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

  FILE *table_xsec_phi2pi = fopen(Form("table_xsec_phi2pi_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections in momentum transfer bins}\n \\begin{tabular}{|c|c|c|c|c|c|}\n \\hline\n -t $[GeV^{2}]$ & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma$ [nb] \\\\ \n \\hline\n");

  FILE *table_xsec_phi2pi_sys = fopen(Form("table_xsec_phi2pi_sys_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phi2pi_sys,"\\documentclass[12pt]{extarticle}\n \\usepackage{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Systematic errors in momentum transfer bins}\n \\begin{tabular}{|c|c|c|}\n \\hline\n -t $[GeV^{2}]$ & polynomial degrees [ \\% ] & Fit range [ \\% ]  \\\\ \n \\hline\n");  
  
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
    lat_phit_mc.SetTextSize(0.05);
    lat_phit_mc.SetTextAlign(13); //align at top
    lat_phit_mc.SetNDC();
    lat_phit_mc.SetTextColor(kBlue);
    lat_phit_mc.DrawLatex(0.45, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
    lat_phit_mc.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_mc, dN_phit_mc));
    lat_phit_mc.DrawLatex(0.45, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
    lat_phit_mc.DrawLatex(0.45, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phit_mc.DrawLatex(0.45, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

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
    lat_phit_data.SetTextSize(0.05);
    lat_phit_data.SetTextAlign(13); //align at top
    lat_phit_data.SetNDC();
    lat_phit_data.SetTextColor(kBlue);
    lat_phit_data.DrawLatex(0.45, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb_data->GetChisquare() / fsb_data->GetNDF()));
    lat_phit_data.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phit_data, dN_phit_data));
    lat_phit_data.DrawLatex(0.45, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb_data->GetParameter(1), fsb_data->GetParError(1)));
    lat_phit_data.DrawLatex(0.45, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb_data->GetParameter(2), fsb_data->GetParError(2)));
    lat_phit_data.DrawLatex(0.45, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb_data->GetParameter(3), fsb_data->GetParError(3)));

    gphit->SetPoint(i - 1, h2phit->GetXaxis()->GetBinCenter(i), N_phit_data);
    gphit->SetPointError(i - 1, 0, dN_phit_data);    

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phi = N_phit_mc / h_t_Thrown->GetBinContent(i); // Efficiency = N_observed/N_generated
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
    double dxsec_phi2pi_sys = xsec_phi2pi * TMath::Sqrt((err_sys / N_phit_data)*(err_sys / N_phit_data) + (deff_phi / eff_phi)*(deff_phi / eff_phi) + (dbr_phi / br_phi)*(dbr_phi / br_phi));
    // double xsec_phi = N_phit / lumi_phi;
    // double dxsec_phi = dN_phit / lumi_phi;  
    gphitxsec->SetPoint(i-1, h2phit->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    gphitxsec->SetPointError(i-1, 0, dxsec_phi2pi_stat);
    // gphitxsec->SetPoint(gphitxsec->GetN(), h2phit->GetXaxis()->GetBinCenter(i), xsec_phi2pi);
    // gphitxsec->SetPointError(gphitxsec->GetN() - 1, 0, dxsec_phi2pi_stat);

    fprintf(table_xsec_phi2pi, "%0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h2phit_mc->GetXaxis()->GetBinCenter(i), h_t_Thrown->GetBinContent(i), N_phit_mc, dN_phit_mc, N_phit_data, dN_phit_data, eff_phi*100, deff_phi*100, xsec_phi2pi, dxsec_phi2pi_stat, dxsec_phi2pi_sys);
    
    fprintf(table_xsec_phi2pi_sys, "%0.2f & %0.2f & %0.2f \\\\ \n", h2phit_mc->GetXaxis()->GetBinCenter(i), err_bkg_poln*100/N_bkg_poln_1[i-1], err_fit_range*100/fit_range_1[i-1]);
    // table_phi << std::setprecision(2) << std::fixed;
    // table_phi <<h2phit_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_t_Thrown->GetBinContent(i)<<"&"<<N_phit_mc<<"$\\pm$"<<dN_phit_mc<<"&"<<N_phit_data<<"$\\pm$"<<dN_phit_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi2pi<<"$\\pm$"<<dxsec_phi2pi_stat<<"$\\pm$"<<dxsec_phi2pi_sys<<" \\\\"<< endl;
    }

    cgphit->cd();
    gphit->Draw("AP");
    gphit->Write(Form("h%s_gphit", name.Data()), TObject::kWriteDelete);
    cgphit_mc->cd();
    gphit_mc->Draw("AP");
    cgphiteff->cd();
    gphiteff->Draw("AP");
    gphiteff->Write(Form("h%s_gphiteff", name.Data()), TObject::kWriteDelete);
    cgphitxsec->cd();
    gphitxsec->Draw("AP");
    gphitxsec->Write(Form("h%s_gphitxsec", name.Data()), TObject::kWriteDelete);
    // int j =1;
    // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);


    // total Yields
    TMultiGraph *mgphit = new TMultiGraph();
    TCanvas *cmgphit = new TCanvas("cmgphit", "cmgphit", 900, 600);
    TLegend *lmgphit = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgphit->SetTextSize(0.05);
    lmgphit->SetBorderSize(0);
    cmgphit->cd();
    cmgphit->SetGrid();
    TGraphErrors *h16_gphit = (TGraphErrors *)outputfig->Get("h16_gphit");
    cout << " ***** h16_gphit = " << h16_gphit << endl;
    h16_gphit->SetMarkerColor(1);
    h16_gphit->SetMarkerSize(1.5);
    h16_gphit->SetLineColor(1);
    h16_gphit->SetMarkerStyle(20);
    mgphit->Add(h16_gphit);
    TGraphErrors *h17_gphit = (TGraphErrors *)outputfig->Get("h17_gphit");
    cout << " ***** h17_gphit = " << h17_gphit << endl;
    h17_gphit->SetMarkerColor(kBlue);
    h17_gphit->SetMarkerSize(1.5);
    h17_gphit->SetLineColor(kBlue);
    h17_gphit->SetMarkerStyle(20);
    mgphit->Add(h17_gphit);
    TGraphErrors *h18_gphit = (TGraphErrors *)outputfig->Get("h18_gphit");
    cout << " ***** h18_gphit = " << h18_gphit << endl;
    h18_gphit->SetMarkerColor(kRed);
    h18_gphit->SetMarkerSize(1.5);
    h18_gphit->SetLineColor(kRed);
    h18_gphit->SetMarkerStyle(20);
    mgphit->Add(h18_gphit);
    TGraphErrors *h18l_gphit = (TGraphErrors *)outputfig->Get("h18l_gphit");
    cout << " ***** h18l_gphit = " << h18l_gphit << endl;
    h18l_gphit->SetMarkerColor(kMagenta);
    h18l_gphit->SetMarkerSize(1.5);
    h18l_gphit->SetLineColor(kMagenta);
    h18l_gphit->SetMarkerStyle(20);
    mgphit->Add(h18l_gphit);
    mgphit->SetTitle("; -t [GeV^{2}]; N_{#phi}");
    mgphit->Draw("AP");
    mgphit->SetMinimum(0.);
    cmgphit->Modified();
    lmgphit->AddEntry(h16_gphit, "2016", "p");
    lmgphit->AddEntry(h17_gphit, "2017", "p");
    lmgphit->AddEntry(h18_gphit, "2018 S", "p");
    lmgphit->AddEntry(h18l_gphit, "2018 F", "p");
    lmgphit->Draw();
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphit.eps", "eps");
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphit.png", "png");
    cmgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgphit.root", "root");

    // total Efficiency
    TMultiGraph *mgteff = new TMultiGraph();
    TCanvas *cmgteff = new TCanvas("cmgteff", "cmgteff", 900, 600);
    TLegend *lmgteff = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgteff->SetTextSize(0.05);
    lmgteff->SetBorderSize(0);
    cmgteff->cd();
    cmgteff->SetGrid();
    TGraphErrors *h16_gphiteff = (TGraphErrors *)outputfig->Get("h16_gphiteff");
    cout << " ***** h16_gphiteff = " << h16_gphiteff << endl;
    h16_gphiteff->SetMarkerColor(1);
    h16_gphiteff->SetMarkerSize(1.5);
    h16_gphiteff->SetLineColor(1);
    h16_gphiteff->SetMarkerStyle(20);
    mgteff->Add(h16_gphiteff);
    TGraphErrors *h17_gphiteff = (TGraphErrors *)outputfig->Get("h17_gphiteff");
    cout << " ***** h17_gphiteff = " << h17_gphiteff << endl;
    h17_gphiteff->SetMarkerColor(kBlue);
    h17_gphiteff->SetMarkerSize(1.5);
    h17_gphiteff->SetLineColor(kBlue);
    h17_gphiteff->SetMarkerStyle(20);
    mgteff->Add(h17_gphiteff);
    TGraphErrors *h18_gphiteff = (TGraphErrors *)outputfig->Get("h18_gphiteff");
    cout << " ***** h18_gphiteff = " << h18_gphiteff << endl;
    h18_gphiteff->SetMarkerColor(kRed);
    h18_gphiteff->SetMarkerSize(1.5);
    h18_gphiteff->SetLineColor(kRed);
    h18_gphiteff->SetMarkerStyle(20);
    mgteff->Add(h18_gphiteff);
    TGraphErrors *h18l_gphiteff = (TGraphErrors *)outputfig->Get("h18l_gphiteff");
    cout << " ***** h18l_gphiteff = " << h18l_gphiteff << endl;
    h18l_gphiteff->SetMarkerColor(kMagenta);
    h18l_gphiteff->SetMarkerSize(1.5);
    h18l_gphiteff->SetLineColor(kMagenta);
    h18l_gphiteff->SetMarkerStyle(20);
    mgteff->Add(h18l_gphiteff);
    mgteff->SetTitle("; -t [GeV^{2}]; #epsilon [%]");
    mgteff->Draw("AP");
    mgteff->SetMinimum(0.);
    cmgteff->Modified();
    lmgteff->AddEntry(h16_gphiteff, "2016", "p");
    lmgteff->AddEntry(h17_gphiteff, "2017", "p");
    lmgteff->AddEntry(h18_gphiteff, "2018 S", "p");
    lmgteff->AddEntry(h18l_gphiteff, "2018 F", "p");
    lmgteff->Draw();
    cmgteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgteff.eps", "eps");
    cmgteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgteff.png", "png");
    cmgteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgteff.root", "root");

    // total Cross-sections
    TMultiGraph *mgtxsec = new TMultiGraph();
    TCanvas *cmgtxsec = new TCanvas("cmgtxsec", "cmgtxsec", 900, 600);
    TLegend *lmgtxsec = new TLegend(0.87, 0.82, 0.98, 0.98);
    lmgtxsec->SetTextSize(0.05);
    lmgtxsec->SetBorderSize(0);
    cmgtxsec->cd();
    cmgtxsec->SetGrid();
    TGraphErrors *h16_gphitxsec = (TGraphErrors *)outputfig->Get("h16_gphitxsec");
    cout << " ***** h16_gphitxsec = " << h16_gphitxsec << endl;
    h16_gphitxsec->SetMarkerColor(1);
    h16_gphitxsec->SetMarkerSize(1.5);
    h16_gphitxsec->SetLineColor(1);
    h16_gphitxsec->SetMarkerStyle(20);
    mgtxsec->Add(h16_gphitxsec);
    TGraphErrors *h17_gphitxsec = (TGraphErrors *)outputfig->Get("h17_gphitxsec");
    cout << " ***** h17_gphitxsec = " << h17_gphitxsec << endl;
    h17_gphitxsec->SetMarkerColor(kBlue);
    h17_gphitxsec->SetMarkerSize(1.5);
    h17_gphitxsec->SetLineColor(kBlue);
    h17_gphitxsec->SetMarkerStyle(20);
    mgtxsec->Add(h17_gphitxsec);
    TGraphErrors *h18_gphitxsec = (TGraphErrors *)outputfig->Get("h18_gphitxsec");
    cout << " ***** h18_gphitxsec = " << h18_gphitxsec << endl;
    h18_gphitxsec->SetMarkerColor(kRed);
    h18_gphitxsec->SetMarkerSize(1.5);
    h18_gphitxsec->SetLineColor(kRed);
    h18_gphitxsec->SetMarkerStyle(20);
    mgtxsec->Add(h18_gphitxsec);
    TGraphErrors *h18l_gphitxsec = (TGraphErrors *)outputfig->Get("h18l_gphitxsec");
    cout << " ***** h18l_gphitxsec = " << h18l_gphitxsec << endl;
    h18l_gphitxsec->SetMarkerColor(kMagenta);
    h18l_gphitxsec->SetMarkerSize(1.5);
    h18l_gphitxsec->SetLineColor(kMagenta);
    h18l_gphitxsec->SetMarkerStyle(20);
    mgtxsec->Add(h18l_gphitxsec);
    mgtxsec->SetTitle("; -t [GeV^{2}]; #sigma [nb]");
    mgtxsec->Draw("AP");
    mgtxsec->SetMinimum(0.);
    cmgtxsec->Modified();
    lmgtxsec->AddEntry(h16_gphitxsec, "2016", "p");
    lmgtxsec->AddEntry(h17_gphitxsec, "2017", "p");
    lmgtxsec->AddEntry(h18_gphitxsec, "2018 S", "p");
    lmgtxsec->AddEntry(h18l_gphitxsec, "2018 F", "p");
    lmgtxsec->Draw();
    cmgtxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgtxsec.eps", "eps");
    cmgtxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgtxsec.png", "png");
    cmgtxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/cmgtxsec.root", "root");

    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit.root", name.Data()), "root");
    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit.eps", name.Data()), "eps");
    cphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit.png", name.Data()), "png");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1.root", name.Data()), "root");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1.eps", name.Data()), "eps");
    cphit1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1.png", name.Data()), "png");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit.root", name.Data()), "root");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit.eps", name.Data()), "eps");
    cgphit->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit.png", name.Data()), "png");
    cphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit_mc.root", name.Data()), "root");
    cphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit_mc.eps", name.Data()), "eps");
    cphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit_mc.png", name.Data()), "png");
    cphit1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1_mc.root", name.Data()), "root");
    cphit1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1_mc.eps", name.Data()), "eps");
    cphit1_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_phit1_mc.png", name.Data()), "png");
    cgphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit_mc.root", name.Data()), "root");
    cgphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit_mc.eps", name.Data()), "eps");
    cgphit_mc->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphit_mc.png", name.Data()), "png");
    cgphiteff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiteff.root", name.Data()), "root");
    cgphiteff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiteff.eps", name.Data()), "eps");
    cgphiteff->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphiteff.png", name.Data()), "png");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphitxsec.root", name.Data()), "root");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphitxsec.eps", name.Data()), "eps");
    cgphitxsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_gphitxsec.png", name.Data()), "png");

    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.eps", name.Data()), "eps");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_tagged_flux.png", name.Data()), "png");
    c_h_t_Thrown->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_h_t_Thrown.root", name.Data()), "root");
    c_h_t_Thrown->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_h_t_Thrown.eps", name.Data()), "eps");
    c_h_t_Thrown->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_xsec_phi2pi/c%s_h_t_Thrown.png", name.Data()), "png");

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

}
