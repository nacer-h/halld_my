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
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", name.Data()));
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

  TFile *outputfig = new TFile("output/fig_xsec_phifo/xsec_phifo.root", "UPDATE");
  
  // ofstream ofs_phifo("phifo_tab.txt", ofstream::out);

  ofstream test;
  test.open ("test.txt");

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.99, mkk_max = 1.2; //mkk_min = 0.98, mkk_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;
  double m2pi_min2 = 0.83, m2pi_max2 = 1.16; // [0.83, 1.14], m2pi_min2 = 0.83, m2pi_max2 = 1.16

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // ##############################################  fo ALL ##############################################

  FILE *table_xsec_phifo = fopen(Form("table_xsec_phifo_%s.tex", name.Data()),"w");
  fprintf(table_xsec_phifo,"\\documentclass[8pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\captionsetup{labelformat=empty}\n \\begin{document}\n \\begin{table}[!htbp]\n \\centering\n \\caption{Total cross-sections}\n \\begin{tabular}{|c|c|c|c|c|c|}\n \\hline\n $E_{\\gamma}$ [GeV] & $N_{generated}~(MC)$ & $N_{measured}~(MC)$ & $N_{measured}~(Data)$ & $\\epsilon$ [ \\% ] & $\\sigma \\times BR_{f_{0}\\rightarrow\\pi^{+}\\pi^{-}}$ [nb] \\\\ \n \\hline\n");

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

  // *********************** f0(980) MC *********************************
  TCanvas *c_foMass_postcuts = new TCanvas("c_foMass_postcuts", "c_foMass_postcuts", 1000, 600);
  TH1F *h_foMass_postcuts = new TH1F("h_foMass_postcuts", ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 200, m2pi_min2, m2pi_max2);//0.85, 1.05
  tmc->Project("h_foMass_postcuts", "pippim_mf", "w8*(pippim_uni && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)"); // && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(abs(rf_dt)<1.5*4.008), w4*(abs(rf_dt)<2.5*4.008)
  c_foMass_postcuts->cd();
  h_foMass_postcuts->Draw("e");
  
  TF1 *fsb_fo_mc = new TF1("fsb_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2]) + pol1(3)", m2pi_min2, m2pi_max2); //[0]*TMath::BreitWigner(x,[1],[2]) + pol3(3)
   fsb_fo_mc->SetLineColor(2);
   fsb_fo_mc->SetParameters(1800, 0.988, 0.06, 1,1); //1800, 0.988, 0.06, 1,1,1,1
  //  fsb_fo_mc->SetParLimits(0, 0, 10000);
   fsb_fo_mc->SetParLimits(1, 0.95, 0.1); // 1.018, 1.021
   fsb_fo_mc->SetParLimits(2, 0.001, 0.1); //fsb->SetParLimits(2, 0.008, 0.010);
  //  fsb_fo_mc->FixParameter(1, fsb_fo_mc->GetParameter(1)); //0.988
  //  fsb_fo_mc->FixParameter(2, fsb_fo_mc->GetParameter(2)); // 0.059, 0.070, 0.063, 0.062

   // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
   TF1 *fs_fo_mc = new TF1("fs_fo_mc", "[0]*TMath::BreitWigner(x,[1],[2])", m2pi_min2, m2pi_max2);
   fs_fo_mc->SetLineColor(4);

   TF1 *fb_fo_mc = new TF1("fb_fo_mc", "pol1(3)", m2pi_min2, m2pi_max2); //pol2(3)
   fb_fo_mc->SetLineColor(28);
   fb_fo_mc->SetLineStyle(2);

   h_foMass_postcuts->Fit("fsb_fo_mc", "", "", m2pi_min2, m2pi_max2);
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

   TCanvas *cphifo = new TCanvas("cphifo", "cphifo", 900, 600); // 900, 600
   TCanvas *cphifo1 = new TCanvas("cphifo1", "cphifo1", 1500, 800);
   cphifo1->Divide(5, 5);
   TCanvas *cphifo2 = new TCanvas("cphifo2", "cphifo2", 1500, 800);
   cphifo2->Divide(5, 5);
   TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 900, 600);
   TGraphErrors *gphifo;

  //  TGraphErrors *gphifo = scanphi("fo", "phifo", Form("%s",name.Data()), n2k, n2pi);
  //  cout<<"******************** gphifo = "<<gphifo<<endl;
  //  cgphifo->cd();
  //  gphifo->Draw("AP");

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
   TH2D *h2d2 = new TH2D("h2d2", ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});m_{K^{+}K^{-}} (GeV/c^{2});Counts", n2pi, m2pi_min, m2pi_max, n2k, mkk_min, mkk_max);
   tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)"); // && (abs(p_dttof)<0.3 || p_dttof == -999) && kin_chisq<55 && abs(mm2)<0.035, w2*(abs(rf_dt)<1.5*4.008), w4*(abs(rf_dt)<2.5*4.008)
   // TH2D *h2d2 = (TH2D *)fdata->Get("h2_PhiVsfoMass_KinFit");
  //  h2d2->RebinX(2);
   cout << " ***** h2d2 = " << h2d2 << endl;
   h2d2->Draw("colz");

   gphifo = new TGraphErrors(); //n2pi
   gphifo->SetMarkerStyle(20);
   gphifo->SetMarkerSize(1.0);
   gphifo->SetMarkerColor(1);
   gphifo->SetMinimum(0.);
   gphifo->SetTitle(";m_{#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi #pi^{+}#pi^{-}}");

   // gphifo_width = new TGraphErrors(n2pi);
   // gphifo_width->SetMarkerStyle(20);
   // gphifo_width->SetMarkerSize(1.0);
   // gphifo_width->SetMarkerColor(1);
   // gphifo_width->SetMinimum(0.);
   // gphifo_width->SetMaximum(0.01);
   // gphifo->SetTitle("(Data);m_{#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi}");

   // gphifo_mean = new TGraphErrors(n2pi);
   // gphifo_mean->SetMarkerStyle(20);
   // gphifo_mean->SetMarkerSize(1.0);
   // gphifo_mean->SetMarkerColor(1);
   // gphifo_mean->SetMinimum(0.);
   // gphifo_mean->SetMaximum(0.01);
   // gphifo->SetTitle("(Data);m_{#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi}");

   // gnophifo->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} (GeV/c^{2});N_{bkg}", Eg1[j], Eg2[j]));

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
    //  fsb->SetParLimits(0, 0, 10000);
     fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
     fsb->FixParameter(2, w.var("sigma_PhiMass_mc")->getVal());  //fsb->SetParLimits(2, 0.008, 0.010);
     fsb->FixParameter(3, w.var("width_PhiMass_mc")->getVal()); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

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
     lat_phifo.SetTextSize(0.05);
     lat_phifo.SetTextAlign(13); //align at top
     lat_phifo.SetNDC();
     lat_phifo.SetTextColor(kBlue);
    //  lat_phifo.DrawLatex(0.45, 0.86, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
     lat_phifo.DrawLatex(0.45, 0.80, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
     lat_phifo.DrawLatex(0.45, 0.74, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
     lat_phifo.DrawLatex(0.45, 0.68, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
     lat_phifo.DrawLatex(0.45, 0.62, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

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

  // TF1 *fsb_fo_data = new TF1("fsb_fo_data", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", m2pi_min2, m2pi_max2);
  TF1 *fsb_fo_data = new TF1("fsb_fo_data", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", m2pi_min2, m2pi_max2); // [0]*TMath::BreitWigner(x,[1],[2]) + expo(3)+pol0(5)
  fsb_fo_data->SetParameters(80, 0.970, 0.06, 6.7, 2.8, 1,1,1);
  fsb_fo_data->SetLineColor(2);
  fsb_fo_data->SetParLimits(1, 0.96, 0.99); // 1, 0.96, 0.99
  fsb_fo_data->SetParLimits(2, 0.01, 0.10);  // 2, 0.01, 0.10
 
  if(name == "16")
  {
  fsb_fo_data->SetParLimits(2, 0.060, 0.070);; //0.971, 0.960, 0.970, 0.965
  // fsb_fo_data->FixParameter(2, 0.071); //0.071, 0.095, 0.064, 0.074
  }
  // if(name == "17")
  // {
  // fsb_fo_data->FixParameter(1, 0.960); //0.971, 0.960, 0.970, 0.965
  // fsb_fo_data->FixParameter(2, 0.095); // 0.071, 0.095, 0.064, 0.074
  // }
  // if(name == "18")
  // {
  // fsb_fo_data->FixParameter(1, 0.970); //0.971, 0.960, 0.970, 0.965
  // fsb_fo_data->FixParameter(2, 0.064); // 0.071, 0.095, 0.064, 0.074
  // }
  // if(name == "18l")
  // {
  // fsb_fo_data->FixParameter(1, 0.965); //0.971, 0.960, 0.970, 0.965
  // fsb_fo_data->FixParameter(2, 0.074); // 0.071, 0.095, 0.064, 0.074
  // }

  // TF1 *fs_fo_data = new TF1("fs_fo_data", "[0]*TMath::Voigt(x - [1], [2], [3])", m2pi_min2, m2pi_max2);
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
  double arr_sig_mean_16[3] = {28.906, 30.856, 25.109}; //mean: mu,-0.01, +0.01 -> 1% of 0.990.  
  double arr_sig_width_16[3] = {28.906, 27.983, 29.818}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_16[3] = {28.906, 24.677, 25.358};  //pol:0,1,2
  double arr_fit_range_16[3] = {28.906, 30.589, 27.356}; //fit: [0.83, 1.14], [0.83, 1.12], [0.83, 1.16]
  double arr_binning_16[3] = {28.906, 28.563, 28.357};   //binning: 100,90,110
  double arr_beambun_16[3] = {28.906, 28.164, 27.665};   //bunches: w8, w4, w2
  double arr_p_dttof_16[3] = {27.780, 28.974, 26.655};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_16[3] = {27.780, 25.954, 26.344}; // kin_chisq<55, 45, 65
  double arr_mm2_16[3] = {27.780, 25.354, 27.403};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "17")
  double arr_sig_mean_17[3] = {40.711, 43.443, 34.830}; //mean: mu,-0.01, +0.01 -> 1% of 0.990. 
  double arr_sig_width_17[3] = {40.711, 39.688, 41.729}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_17[3] = {40.711, 32.637, 34.809};  //pol:0,1,2
  double arr_fit_range_17[3] = {40.711, 42.637, 38.764}; //fit: [0.83, 1.14], [0.83, 1.12], [0.83, 1.16]
  double arr_binning_17[3] = {40.711, 41.628, 40.957};   //binning: 100,90,110
  double arr_beambun_17[3] = {40.711, 40.145, 39.461};   //bunches: w8, w4, w2
  double arr_p_dttof_17[3] = {41.541, 42.819, 41.183};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_17[3] = {41.541, 46.866, 41.620}; // kin_chisq<55, 45, 65
  double arr_mm2_17[3] = {41.541, 42.382, 43.550};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18")
  double arr_sig_mean_18[3] = {18.241, 17.911, 16.529}; //mean: mu,-0.01, +0.01 -> 1% of 0.990. 
  double arr_sig_width_18[3] = {18.241, 18.148, 19.291}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_18[3] = {18.241, 18.469, 19.046};  //pol:0,1,2
  double arr_fit_range_18[3] = {18.241, 19.399, 18.229}; //fit: [0.83, 1.14], [0.83, 1.12], [0.83, 1.16]
  double arr_binning_18[3] = {18.241, 18.455, 18.553};   //binning: 100,90,110
  double arr_beambun_18[3] = {18.241, 18.453, 18.893};   //bunches: w8, w4, w2
  double arr_p_dttof_18[3] = {18.575, 18.395, 18.250};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18[3] = {18.575, 18.432, 18.337}; // kin_chisq<55, 45, 65
  double arr_mm2_18[3] = {18.575, 18.814, 18.509};       // abs(mm2)<0.035, 0.025, 0.045

  // if(name == "18l")
  double arr_sig_mean_18l[3] = {21.780, 21.932, 18.836}; //mean: mu,-0.01, +0.01 -> 1% of 0.990. 
  double arr_sig_width_18l[3] = {21.780, 21.160, 22.399}; //width: gamma, -0.002, +0.002 -> 5% of 0.048.
  double arr_bkg_poln_18l[3] = {21.780, 17.763, 18.951};  //pol:0,1,2
  double arr_fit_range_18l[3] = {21.780, 22.929, 20.484}; //fit: [0.83, 1.14], [0.83, 1.12], [0.83, 1.16]
  double arr_binning_18l[3] = {21.780, 21.718, 21.778};   //binning: 100,90,110
  double arr_beambun_18l[3] = {21.780, 21.986, 22.376};   //bunches: w8, w4, w2
  double arr_p_dttof_18l[3] = {21.885, 21.307, 22.238};   // abs(p_dttof)<0.3, 0.2, 0.4
  double arr_kin_chisq_18l[3] = {21.885, 22.622, 22.080}; // kin_chisq<55, 45, 65
  double arr_mm2_18l[3] = {21.885, 21.803, 21.968};       // abs(mm2)<0.035, 0.025, 0.045

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
  double dxsec_phifo_sys = xsec_phifo * TMath::Sqrt((err_sys*err_sys) + (dbr_phi / br_phi)*(dbr_phi / br_phi)); //(err_sys / N_fo_data) * (err_sys / N_fo_data) +
  // double dxsec_phifo_sys = xsec_phifo * TMath::Sqrt((deff_phifo / eff_phifo) * (deff_phifo / eff_phifo) + (dbr_phi / br_phi) * (dbr_phi / br_phi)); //(err_sys / N_fo_data) * (err_sys / N_fo_data) +
  double dxsec_phifo_tot = TMath::Sqrt(dxsec_phifo_stat*dxsec_phifo_stat + dxsec_phifo_sys*dxsec_phifo_sys);

  test << std::setprecision(3) << std::fixed;
  test << "xsec_phifo = " << xsec_phifo << endl;

  fprintf(table_xsec_phifo, "%0.2f - %0.2f & %0.f & %0.f $\\pm$ %0.f & %0.f $\\pm$ %0.f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f $\\pm$ %0.2f \\\\ \n", h_tagged_flux->GetXaxis()->GetXmin(), h_tagged_flux->GetXaxis()->GetXmax(), h_beame_tru->Integral(0, 100), N_fo_mc, dN_fo_mc, N_fo_data, dN_fo_data, eff_phifo * 100, deff_phifo * 100, xsec_phifo, dxsec_phifo_stat, dxsec_phifo_sys);

  fprintf(table_xsec_phifo_sys, "%0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", err_bkg_poln * 100, err_fit_range * 100, err_binning * 100, err_beambun * 100, err_p_dttof * 100, err_kin_chisq * 100, err_mm2 * 100);

  // table_phi << std::setprecision(2) << std::fixed;
  // table_phi <<h2phie_mc->GetXaxis()->GetBinCenter(i)<<"&"<<h_beame_tru->GetBinContent(i)<<"&"<<N_fo_mc<<"$\\pm$"<<dN_fo_mc<<"&"<<N_fo_data<<"$\\pm$"<<dN_fo_data<<"&"<<eff_phi*100<<"$\\pm$"<<deff_phi*100<<"&"<<xsec_phifo<<"$\\pm$"<<dxsec_phifo_stat<<"$\\pm$"<<dxsec_phifo_sys<<" \\\\"<< endl;

  // c_PhiMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.root", name.Data()), "root");
  // c_PhiMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.eps", name.Data()), "eps");
  // c_PhiMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_PhiMass_postcuts_fitted_%s.png", name.Data()), "png");
  c_foMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.root", name.Data()), "root");
  c_foMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.eps", name.Data()), "eps");
  c_foMass_postcuts->Print(Form("output/fig_xsec_phifo/cmc_foMass_postcuts_fitted_%s.png", name.Data()), "png");
  cphifo1->Print(Form("output/fig_xsec_phifo/c1_phifo_%s.root", name.Data()), "root");
  cphifo1->Print(Form("output/fig_xsec_phifo/c1_phifo_%s.eps", name.Data()), "eps");
  cphifo1->Print(Form("output/fig_xsec_phifo/c1_phifo_%s.png", name.Data()), "png");
  cphifo2->Print(Form("output/fig_xsec_phifo/c2_phifo_%s.root", name.Data()), "root");
  cphifo2->Print(Form("output/fig_xsec_phifo/c2_phifo_%s.eps", name.Data()), "eps");
  cphifo2->Print(Form("output/fig_xsec_phifo/c2_phifo_%s.png", name.Data()), "png");
  cphifo->Print(Form("output/fig_xsec_phifo/c_phifo_%s.root", name.Data()), "root");
  cphifo->Print(Form("output/fig_xsec_phifo/c_phifo_%s.eps", name.Data()), "eps");
  cphifo->Print(Form("output/fig_xsec_phifo/c_phifo_%s.png", name.Data()), "png");
  cgphifo->Print(Form("output/fig_xsec_phifo/c_gphifo_%s.root", name.Data()), "root");
  cgphifo->Print(Form("output/fig_xsec_phifo/c_gphifo_%s.eps", name.Data()), "eps");
  cgphifo->Print(Form("output/fig_xsec_phifo/c_gphifo_%s.png", name.Data()), "png");

  // c_tagged_flux->Print(Form("output/fig_xsec_phifo/c%s_tagged_flux.root", name.Data()), "root");
  // c_tagged_flux->Print(Form("output/fig_xsec_phifo/c%s_tagged_flux.eps", name.Data()), "eps");
  // c_tagged_flux->Print(Form("output/fig_xsec_phifo/c%s_tagged_flux.png", name.Data()), "png");
  // c_beame_tru->Print(Form("output/fig_xsec_phifo/c_beame_tru_%s.root", name.Data()), "root");
  // c_beame_tru->Print(Form("output/fig_xsec_phifo/c_beame_tru_%s.eps", name.Data()), "eps");
  // c_beame_tru->Print(Form("output/fig_xsec_phifo/c_beame_tru_%s.png", name.Data()), "png");

  outputfig->Print();

  fprintf(table_xsec_phifo, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
  fclose(table_xsec_phifo);
  gSystem->Exec(Form("pdflatex table_xsec_phifo_%s.tex", name.Data()));

   fprintf(table_xsec_phifo_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
   fclose(table_xsec_phifo_sys);
   gSystem->Exec(Form("pdflatex table_xsec_phifo_sys_%s.tex", name.Data()));

   // cgphifo_width->Print("output/fig_xsec_phifo/c_gphifo_width.root", "root");
   // cgphifo_width->Print("output/fig_xsec_phifo/c_gphifo_width.eps", "eps");
   // cgphifo_mean->Print("output/fig_xsec_phifo/c_gphifo_mean.root", "root");
   // cgphifo_mean->Print("output/fig_xsec_phifo/c_gphifo_mean.eps", "eps");

   // cgnophifo->Print("output/fig_xsec_phifo/c_gnophifo.root", "root");
   // cgnophifo->Print("output/fig_xsec_phifo/c_gnophifo.eps", "eps");

}
