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

void scanphi(TString name, int n2k=100, int ne=10, int nt=10)//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
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
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/scanphi.root", "UPDATE");

  ofstream ofs_scanphi("tableau_phi.txt", ofstream::out);
  
  ofstream table_phi("table_phi.tex", ofstream::out);
  table_phi <<"\\documentclass{article}" << endl;
  table_phi <<"\\usepackage[margin=0.5in]{geometry}" << endl; 
  table_phi <<"\\usepackage{tabularx}" << endl;
  table_phi <<"\\begin{document}" << endl;
  table_phi <<"\\begin{table}[!htbp]" << endl;
  table_phi <<"\\begin{center}" << endl;
  table_phi <<"\\caption{Total cross-sections in photon energy bins.} "<< endl;
  // table_phi <<"\\label{tab : table1}" << endl;
  table_phi <<"\\begin{tabularx}{\\textwidth}{|c|X|X|X|X|c|}" << endl;
  table_phi <<"\\hline" << endl;
  table_phi <<"$E_{\\gamma} [GeV]$ & $N_{generated}$ & $N_{observed}$ & $N_{measured}$ & $\\epsilon$ [\\%]& $\\sigma [nb]$ \\\\" << endl;
  // table_phi << "\\alpha & \\beta & \\gamma \\" << endl;
  table_phi <<"\\hline" << endl;

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.98, mkk_max = 1.2;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 900, 600);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux" << h_tagged_flux << endl;
  h_tagged_flux->Rebin(10);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.0);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 900, 600);
  c_beame_tru->cd();
  TH1F *h_beame_tru = new TH1F("h_beame_tru", "MC truth; E_{#gamma} [GeV]; Counts", ne, 5.9, 11.9);
  ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru" << h_beame_tru << endl;
  // h_beame_tru->Rebin(10);
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.0);
  h_beame_tru->Draw("e");

  // ======================================== Phi vs. Eg ===============================================
  // root -l 'scanphi.C+(100,20,20)'
  // Data
  TCanvas *cphie = new TCanvas("cphie", "cphie", 900, 600);
  cphie->cd();
  TH2D *h2phie = new TH2D("h2phie", "Data; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, 5.9, 11.9, n2k, 0.98, 1.2);
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");
  h2phie->Draw("colz");

  TCanvas *cphie1 = new TCanvas("cphie1", "cphie1", 1500, 600);
  cphie1->Divide(5, 2);
  TCanvas *cgphie = new TCanvas("cgphie", "cgphie", 900, 600);
  TGraphErrors *gphie = new TGraphErrors();
  gphie->SetMarkerStyle(20);
  gphie->SetMarkerSize(1.0);
  gphie->SetMarkerColor(1);
  gphie->SetMinimum(0.0);
  gphie->SetTitle("#phi(1020) Yield vs. E_{#gamma} (data); E_{#gamma} [GeV]; N_{#phi}");

  // Monte Carlo
  TCanvas *cphie_mc = new TCanvas("cphie_mc", "cphie_mc", 900, 600);
  cphie_mc->cd();
  TH2D *h2phie_mc = new TH2D("h2phie_mc", "MC; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]", ne, 5.9, 11.9, n2k, 0.98, 1.2);
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni)");
  h2phie_mc->Draw("colz");

  TCanvas *cphie1_mc = new TCanvas("cphie1_mc", "cphie1_mc", 1500, 600);
  cphie1_mc->Divide(5, 2);
  TCanvas *cgphie_mc = new TCanvas("cgphie_mc", "cgphie_mc", 900, 600);
  TGraphErrors *gphie_mc = new TGraphErrors();
  gphie_mc->SetMarkerStyle(20);
  gphie_mc->SetMarkerSize(1.0);
  gphie_mc->SetMarkerColor(1);
  gphie_mc->SetMinimum(0.0);
  gphie_mc->SetTitle("#phi(1020) Yield vs. E_{#gamma} (MC); E_{#gamma} [GeV]; N_{#phi}");

  // Efficiency
  TCanvas *cgphieeff = new TCanvas("cgphieeff", "cgphieeff", 900, 600);
  TGraphErrors *gphieeff = new TGraphErrors();
  gphieeff->SetMarkerStyle(20);
  gphieeff->SetMarkerSize(1.0);
  gphieeff->SetMarkerColor(1);
  gphieeff->SetMinimum(0.0);
  gphieeff->SetTitle("#phi(1020) Effeciency vs. E_{#gamma}; E_{#gamma} [GeV]; #epsilon [%]");

  // Cross-section
  TCanvas *cgphiexsec = new TCanvas("cgphiexsec", "cgphiexsec", 900, 600);
  TGraphErrors *gphiexsec = new TGraphErrors();
  gphiexsec->SetMarkerStyle(20);
  gphiexsec->SetMarkerSize(1.0);
  gphiexsec->SetMarkerColor(1);
  gphiexsec->SetMinimum(0.0);
  gphiexsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. E_{#gamma}; E_{#gamma} [GeV]; #sigma [nb]"); //#phi(1020) flux normalized yield, Yield_{#phi}

  for (int i = 1; i <= ne; ++i)
  {
    cout << i << " " << flush;

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphie1->cd(i);
    TH1D *hphie_py = h2phie->ProjectionY(Form("hphie_py_%d", i), i, i);
    // hphie_py->SetTitle(Form("E_{#gamma} = %f",h2phie->GetXaxis()->GetBinCenter(i)));
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

    TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
    fsb->SetLineColor(2);
    fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1, 1, 1);
    fsb->SetParLimits(0, 0, 10000);
    fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb->FixParameter(2, 0.004);        //fsb->SetParLimits(2, 0.008, 0.010);
    fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    fs->SetLineColor(4);

    TF1 *fb = new TF1("fb", "pol4(4)", mkk_min, mkk_max); //pol2(3)
    fb->SetLineColor(28);
    fb->SetLineStyle(2);

    hphie_py->Fit("fsb", "", "", mkk_min, mkk_max);
    double par[8]; //6
    fsb->GetParameters(&par[0]);
    fs->SetParameters(&par[0]);
    fb->SetParameters(&par[4]); //3

    fs->Draw("same");
    fb->Draw("same");

    double N_phie = fs->Integral(mkk_min, mkk_max) / hphie_py->GetBinWidth(1);
    double dN_phie = N_phie * fsb->GetParError(0) / fsb->GetParameter(0);

    TLatex lat_phie;
    lat_phie.SetTextSize(0.06);
    lat_phie.SetTextAlign(13); //align at top
    lat_phie.SetNDC();
    lat_phie.SetTextColor(kBlue);
    lat_phie.DrawLatex(0.45, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    lat_phie.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N_phie, dN_phie));
    lat_phie.DrawLatex(0.45, 0.73, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    lat_phie.DrawLatex(0.45, 0.66, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
    lat_phie.DrawLatex(0.45, 0.59, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    gphie->SetPoint(i - 1, h2phie->GetXaxis()->GetBinCenter(i), N_phie);
    gphie->SetPointError(i - 1, 0, dN_phie);

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphie1_mc->cd(i);
    TH1D *hphie_mc_py = h2phie_mc->ProjectionY(Form("hphie_mc_py_%d", i), i, i);
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

    TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", mkk_min, mkk_max);
    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
    fsb_mc->SetLineColor(2);
    fsb_mc->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1, 1, 1);
    fsb_mc->SetParLimits(0, 0, 10000);
    fsb_mc->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb_mc->FixParameter(2, 0.004);        //fsb->SetParLimits(2, 0.008, 0.010);
    fsb_mc->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    // TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
    fs_mc->SetLineColor(4);

    TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", mkk_min, mkk_max); //pol2(3)
    fb_mc->SetLineColor(28);
    fb_mc->SetLineStyle(2);

    hphie_mc_py->Fit("fsb_mc", "", "", mkk_min, mkk_max);
    double par_mc[8]; //6
    fsb_mc->GetParameters(&par_mc[0]);
    fs_mc->SetParameters(&par_mc[0]);
    fb_mc->SetParameters(&par_mc[4]); //3

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
    lat_phie_mc.DrawLatex(0.45, 0.66, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
    lat_phie_mc.DrawLatex(0.45, 0.59, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

    gphie_mc->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), N_phie_mc);
    gphie_mc->SetPointError(i - 1, 0, dN_phie_mc);

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phi = N_phie_mc / h_beame_tru->GetBinContent(i); // Efficiency = N_observed/N_generated
    double deff_phi = eff_phi * (dN_phie_mc / N_phie_mc);
    gphieeff->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phi*100);
    gphieeff->SetPointError(i - 1, 0, deff_phi*100); //->GetBinError(i)

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phi = h_tagged_flux->GetBinContent(i); // * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    if (lumi_phi <= 0)
    {
      //gphiexsec->RemovePoint(i - 1);
      continue;
    }
    double T = 1.273; //Target thickness [b^{-1}]
    double br_phi = 0.492;
    double dbr_phi = 0.005;
    double xsec_phi = 1e9 * N_phie / (eff_phi * lumi_phi * T * br_phi);
    double dxsec_phi = xsec_phi * sqrt((dN_phie / N_phie) * (dN_phie / N_phie) + (deff_phi / eff_phi) * (deff_phi / eff_phi) + (dbr_phi / br_phi) * (dbr_phi / br_phi));
    // double xsec_phi = N_phie / lumi_phi;
    // double dxsec_phi = dN_phie / lumi_phi;
    gphiexsec->SetPoint(/*i - 1*/ gphiexsec->GetN(), h2phie->GetXaxis()->GetBinCenter(i), xsec_phi);
    gphiexsec->SetPointError(/*i - 1*/ gphiexsec->GetN() - 1, 0, dxsec_phi);

    // c2->Update();
    //sleep(1);
    ofs_scanphi << " i = " << i << " | N_{#phi} = " << N_phie << "#pm" << dN_phie << " |  N_{#phi MC} = " << N_phie_mc << "#pm" << dN_phie_mc << " | N_{generated} = " << h_beame_tru->GetBinContent(i) << " | N_{observed} = " << h_tagged_flux->GetBinContent(i) << " | #epsilon = " << eff_phi << "#pm" << deff_phi << " | #sigma = " << xsec_phi << "#pm" << dxsec_phi << endl;

    table_phi << std::setprecision(4); // << std::fixed;
    table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_phie_mc<<"$\\pm$"<<dN_phie_mc<<"&"<<N_phie<<"$\\pm$"<<dN_phie<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phi<<"$\\pm$"<<dxsec_phi<<" \\\\"<< endl;
    }

    cgphie->cd();
    gphie->Draw("AP");
    gphie->Write(Form("h%s_gphie", name.Data()), TObject::kWriteDelete);
    cgphie_mc->cd();
    gphie_mc->Draw("AP");
    cgphieeff->cd();
    gphieeff->Draw("AP");
    cgphiexsec->cd();
    gphiexsec->Draw("AP");
    gphiexsec->Write(Form("h%s_gphiexsec", name.Data()), TObject::kWriteDelete);
    // int j =1;
    // gphie->Write(Form("grphie_%d", j), TObject::kWriteDelete);

    // total Cross-sections
    TMultiGraph *mgxsec = new TMultiGraph();
    TCanvas *cmgxsec = new TCanvas("cmgxsec", "cmgxsec", 900, 600);
    TLegend *lmgxsec = new TLegend(0.86, 0.86, 0.98, 0.98);
    lmgxsec->SetTextSize(0.04);
    lmgxsec->SetBorderSize(0);
    cmgxsec->cd();
    cmgxsec->SetGrid();
    TGraphErrors *hdata_16_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_16_gphiexsec");
    cout << " ***** hdata_16_gphiexsec = " << hdata_16_gphiexsec << endl;
    hdata_16_gphiexsec->SetMarkerColor(1);
    hdata_16_gphiexsec->SetMarkerSize(1.5);
    hdata_16_gphiexsec->SetLineColor(1);
    hdata_16_gphiexsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_16_gphiexsec);
    TGraphErrors *hdata_17_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_17_gphiexsec");
    cout << " ***** hdata_17_gphiexsec = " << hdata_17_gphiexsec << endl;
    hdata_17_gphiexsec->SetMarkerColor(kBlue);
    hdata_17_gphiexsec->SetMarkerSize(1.5);
    hdata_17_gphiexsec->SetLineColor(kBlue);
    hdata_17_gphiexsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_17_gphiexsec);
    TGraphErrors *hdata_18_gphiexsec = (TGraphErrors *)outputfig->Get("hdata_18_gphiexsec");
    cout << " ***** hdata_18_gphiexsec = " << hdata_18_gphiexsec << endl;
    hdata_18_gphiexsec->SetMarkerColor(kRed);
    hdata_18_gphiexsec->SetMarkerSize(1.5);
    hdata_18_gphiexsec->SetLineColor(kRed);
    hdata_18_gphiexsec->SetMarkerStyle(20);
    mgxsec->Add(hdata_18_gphiexsec);
    mgxsec->SetTitle("#phi#pi^{+}#pi^{-} total Cross-section vs. E_{#gamma}; E_{#gamma} [GeV]; #sigma [nb]");
    mgxsec->Draw("AP");
    mgxsec->SetMinimum(0.);
    cmgxsec->Modified();
    lmgxsec->AddEntry(hdata_16_gphiexsec, "2016", "p");
    lmgxsec->AddEntry(hdata_17_gphiexsec, "2017", "p");
    lmgxsec->AddEntry(hdata_18_gphiexsec, "2018", "p");
    lmgxsec->Draw();
    cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cmgxsec.eps", "eps");
    cmgxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cmgxsec.root", "root");

    // total Yields
    TMultiGraph *mgphie = new TMultiGraph();
    TCanvas *cmgphie = new TCanvas("cmgphie", "cmgphie", 900, 600);
    TLegend *lmgphie = new TLegend(0.86, 0.86, 0.98, 0.98);
    lmgphie->SetTextSize(0.04);
    lmgphie->SetBorderSize(0);
    cmgphie->cd();
    cmgphie->SetGrid();
    TGraphErrors *hdata_16_gphie = (TGraphErrors *)outputfig->Get("hdata_16_gphie");
    cout << " ***** hdata_16_gphie = " << hdata_16_gphie << endl;
    hdata_16_gphie->SetMarkerColor(1);
    hdata_16_gphie->SetMarkerSize(1.5);
    hdata_16_gphie->SetLineColor(1);
    hdata_16_gphie->SetMarkerStyle(20);
    mgphie->Add(hdata_16_gphie);
    TGraphErrors *hdata_17_gphie = (TGraphErrors *)outputfig->Get("hdata_17_gphie");
    cout << " ***** hdata_17_gphie = " << hdata_17_gphie << endl;
    hdata_17_gphie->SetMarkerColor(kBlue);
    hdata_17_gphie->SetMarkerSize(1.5);
    hdata_17_gphie->SetLineColor(kBlue);
    hdata_17_gphie->SetMarkerStyle(20);
    mgphie->Add(hdata_17_gphie);
    TGraphErrors *hdata_18_gphie = (TGraphErrors *)outputfig->Get("hdata_18_gphie");
    cout << " ***** hdata_18_gphie = " << hdata_18_gphie << endl;
    hdata_18_gphie->SetMarkerColor(kRed);
    hdata_18_gphie->SetMarkerSize(1.5);
    hdata_18_gphie->SetLineColor(kRed);
    hdata_18_gphie->SetMarkerStyle(20);
    mgphie->Add(hdata_18_gphie);
    mgphie->SetTitle("#phi(1020) Yield vs. E_{#gamma} (data); E_{#gamma} [GeV]; N_{#phi}");
    mgphie->Draw("AP");
    mgphie->SetMinimum(0.);
    cmgphie->Modified();
    lmgphie->AddEntry(hdata_16_gphie, "2016", "p");
    lmgphie->AddEntry(hdata_17_gphie, "2017", "p");
    lmgphie->AddEntry(hdata_18_gphie, "2018", "p");
    lmgphie->Draw();
    cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cmgphie.eps", "eps");
    cmgphie->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cmgphie.root", "root");

    cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_phie.root", name.Data()), "root");
    cphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_phie.eps", name.Data()), "eps");
    cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_phie1.root", name.Data()), "root");
    cphie1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_phie1.eps", name.Data()), "eps");
    cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_gphie.root", name.Data()), "root");
    cgphie->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_gphie.eps", name.Data()), "eps");
    cphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphie_mc.root", "root");
    cphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphie_mc.eps", "eps");
    cphie1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1_mc.root", "root");
    cphie1_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1_mc.eps", "eps");
    cgphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie_mc.root", "root");
    cgphie_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie_mc.eps", "eps");
    cgphieeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphieeff.root", "root");
    cgphieeff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphieeff.eps", "eps");
    cgphiexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_gphiexsec.root", name.Data()), "root");
    cgphiexsec->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_gphiexsec.eps", name.Data()), "eps");
    /*
  // ======================================== Phi vs. -t ==================================================
  // root -l 'scanphi.C+(100,6,20)'
  //++++++  Data
  TCanvas *cbeamet = new TCanvas("cbeamet", "cbeamet", 600, 400);
  cbeamet->cd();
  TH2D *h2beamet = new TH2D("h2beamet", "Data ;E_{#gamma} [GeV/c^{2}];-t [GeV/c]", 100, 6.3, 11.6, 100, 0.1, 2);
  tdata->Project("h2beamet", "-t_kin:beam_p4_kin.E()", "w8*("+cut+")");
  h2beamet->Draw("colz");

  TCanvas *cphit = new TCanvas("cphit", "cphit", 1500, 800);
  cphit->Divide(2, 3);

  TCanvas *cphit1[ne];

  TCanvas *cgphit = new TCanvas("cgphit", "cgphit", 1500, 800);
  cgphit->Divide(2, 3);
  TGraphErrors *gphit[ne];

  //++++++  mc
  TCanvas *cbeamet_mc = new TCanvas("cbeamet_mc", "cbeamet_mc", 600, 400);
  cbeamet_mc->cd();
  TH2D *h2beamet_mc = new TH2D("h2beamet_mc", "MC ;E_{#gamma} [GeV/c^{2}];-t [GeV/c]", 100, 6.3, 11.6, 100, 0.1, 2);
  tmc->Project("h2beamet_mc", "-t_kin:beam_p4_kin.E()", "w8*("+cut+")");
  h2beamet_mc->Draw("colz");

  TCanvas *cphit_mc = new TCanvas("cphit_mc", "cphit_mc", 1500, 800);
  cphit_mc->Divide(2, 3);

  TCanvas *cphit1_mc[ne];

  TCanvas *cgphit_mc = new TCanvas("cgphit_mc", "cgphit_mc", 1500, 800);
  cgphit_mc->Divide(2, 3);
  TGraphErrors *gphit_mc[ne];

  //++++++ Efficiency
  TCanvas *cgphiteff = new TCanvas("cgphiteff", "cgphiteff", 1500, 800);
  cgphiteff->Divide(2, 3);
  TGraphErrors *gphiteff[ne];

  //+++++++ Differential Cross-section
  TCanvas *cgphitxsec = new TCanvas("cgphitxsec", "cgphitxsec", 1500, 800);
  cgphitxsec->Divide(2, 3);
  TGraphErrors *gphitxsec[ne];

  double Egmin = 6.3;  //= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 11.6; //= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax - Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
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
    
      ofs_scanphi << " j = " << j << " i = " << i << " | N_phit = " << N_phit << " | N_phit_mc = " << N_phit_mc << " | h_beame_tru->GetBinContent(i) = " << h_beame_tru->GetBinContent(j) << " | h_tagged_flux->Integral(bin1, bin2) = " << h_tagged_flux->Integral(bin1, bin2) << " | eff_phi = " << eff_phi << " | xsec_phi = " << xsec_phi << endl;

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

    cphit1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_%d.root", j), "root");
    cphit1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_%d.eps", j), "eps");
    cphit1_mc[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_mc_%d.root", j), "root");
    cphit1_mc[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_mc_%d.eps", j), "eps");
  }
  // int j =1;
  // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);

 cbeamet->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cbeamet.root", "root");
  cbeamet->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cbeamet.eps", "eps");
  cphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit.root", "root");
  cphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit.eps", "eps");
  cgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit.root", "root");
  cgphit->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit.eps", "eps");
  cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit_mc.root", "root");
  cphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cphit_mc.eps", "eps");
  cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit_mc.root", "root");
  cgphit_mc->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit_mc.eps", "eps");
  cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiteff.root", "root");
  cgphiteff->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiteff.eps", "eps");
  cgphitxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphitxsec.root", "root");
  cgphitxsec->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/cgphitxsec.eps", "eps");
*/
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_tagged_flux.root", name.Data()), "root");
    c_tagged_flux->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c%s_tagged_flux.eps", name.Data()), "eps");

    c_beame_tru->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_beame_tru.root", "root");
    c_beame_tru->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_beame_tru.eps", "eps");

    outputfig->Print();

    table_phi <<"\\hline"<<endl;
    table_phi <<"\\end{tabularx}" << endl;
    table_phi <<"\\end{center}" << endl;
    table_phi <<"\\end{table}" << endl;
    table_phi <<"\\end{document}" << endl;
    table_phi.close();
    gSystem->Exec("pdflatex table_phi.tex");

}
