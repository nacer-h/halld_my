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

void scanphi(int n2k=100, int ne=20, int nt=20, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fdata = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm_17v21.root");
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TFile *fmc = new TFile("/Users/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v03.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TFile *fps = new TFile("/Users/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  TFile *ftru = new TFile("/Users/nacer/halld_my/pippimkpkm/input/tree_thrown_phifo_genr8_17v3.root");
  TTree *ttru = (TTree *)ftru->Get("Thrown_Tree");  
  TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/scanphi.root","UPDATE");

  ofstream ofs_scanphi("tableau_phi.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.98, mkk_max = 1.2;

  // +++ PS flux
  TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 600, 400);
  c_tagged_flux->cd();
  TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
  cout << "h_tagged_flux" << h_tagged_flux << endl;
  // h_tagged_flux->Rebin(50);
  h_tagged_flux->SetMarkerStyle(20);
  h_tagged_flux->SetMarkerSize(1.0);
  h_tagged_flux->Draw("e");

  // +++ Thrown Beam Energy
  TCanvas *c_beame_tru = new TCanvas("cbeame_tru", "cbeame_tru", 600, 400);
  c_beame_tru->cd();
  TH1F *h_beame_tru = new TH1F("h_beame_tru","; E_{#gamma} [GeV]; Counts", ne, 6.3, 11.6);
  ttru->Project("h_beame_tru", "ThrownBeam__P4.E()");
  cout << "h_beame_tru" << h_beame_tru << endl;
  h_beame_tru->SetMarkerStyle(20);
  h_beame_tru->SetMarkerSize(1.0);
  h_beame_tru->Draw("e");  

  // ======================================== Phi vs. Eg ===============================================
  // root -l 'scanphi.C+(100,20,20)'
  // Data
  TCanvas *cphie = new TCanvas("cphie","cphie",600, 400);
  cphie->cd();
  TH2D *h2phie = new TH2D("h2phie","Data; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]",ne,6.3,11.6, n2k,0.98,1.2);
  tdata->Project("h2phie", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni &&"+cut+")");
  h2phie->Draw("colz");

  TCanvas *cphie1 = new TCanvas("cphie1", "cphie1", 1500, 800);
  cphie1->Divide(5,4);
  TCanvas *cgphie=new TCanvas("cgphie","cgphie",600,400);
  TGraphErrors *gphie = new TGraphErrors(ne);
  gphie->SetMarkerStyle(20);
  gphie->SetMarkerSize(1.0);
  gphie->SetMarkerColor(1);
  gphie->SetMinimum(0.0);
  gphie->SetTitle("#phi(1020) Yield vs. E_{#gamma} (data); E_{#gamma} [GeV]; N_{#phi}");

  // Monte Carlo
  TCanvas *cphie_mc = new TCanvas("cphie_mc","cphie_mc",600, 400);
  cphie_mc->cd();
  TH2D *h2phie_mc = new TH2D("h2phie_mc","MC; E_{#gamma} [GeV]; m_{K^{+}K^{-}} [GeV/c^{2}]",ne,6.3,11.6, n2k,0.98,1.2);
  tmc->Project("h2phie_mc", "kpkm_mf:beam_p4_kin.E()", "w8*(kpkm_uni &&"+cut+")");
  h2phie_mc->Draw("colz");

  TCanvas *cphie1_mc = new TCanvas("cphie1_mc", "cphie1_mc", 1500, 800);
  cphie1_mc->Divide(5,4);
  TCanvas *cgphie_mc=new TCanvas("cgphie_mc","cgphie_mc",600,400);
  TGraphErrors *gphie_mc = new TGraphErrors(ne);
  gphie_mc->SetMarkerStyle(20);
  gphie_mc->SetMarkerSize(1.0);
  gphie_mc->SetMarkerColor(1);
  gphie_mc->SetMinimum(0.0);
  gphie_mc->SetTitle("#phi(1020) Yield vs. E_{#gamma} (MC); E_{#gamma} [GeV]; N_{#phi}");

  // Efficiency
  TCanvas *cgphieeff=new TCanvas("cgphieeff","cgphieeff",600,400);
  TGraphErrors *gphieeff = new TGraphErrors(ne);
  gphieeff->SetMarkerStyle(20);
  gphieeff->SetMarkerSize(1.0);
  gphieeff->SetMarkerColor(1);
  gphieeff->SetMinimum(0.0);
  gphieeff->SetTitle("#phi(1020) Effeciency vs. E_{#gamma}; E_{#gamma} [GeV]; #epsilon_{#phi}");  

  // Cross-section
  TCanvas *cgphiexsec=new TCanvas("cgphiexsec","cgphiexsec",600,400);
  TGraphErrors *gphiexsec = new TGraphErrors(ne);
  gphiexsec->SetMarkerStyle(20);
  gphiexsec->SetMarkerSize(1.0);
  gphiexsec->SetMarkerColor(1);
  gphiexsec->SetMinimum(0.0);
  gphiexsec->SetTitle("#phi(1020) flux normalized yield vs. E_{#gamma}; E_{#gamma} [GeV]; Yield_{#phi}");

  for (int i = 1; i <= ne; ++i)
  {
    cout << i << " " << flush;

    // +++++++++++++++++++++++++ data  ++++++++++++++++++++
    cphie1->cd(i);
    TH1D *hphie_py = h2phie->ProjectionY(Form("hphie_py_%d",i), i, i);
    // hphie_py->SetTitle(Form("E_{#gamma} = %f",h2phie->GetXaxis()->GetBinCenter(i)));
    hphie_py->Draw("e");

    w.factory(Form("Voigtian::sig_phie(m_phie[%f,%f],mean_phie[1.010,1.030],width_phie[0.004],sigma_phie[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phie[0.001,0.01], mean_phie[1.016,1.022]
    w.factory("Chebychev::bkg_phie(m_phie,{c0_phie[-10,10], c1_phie[-10,10], c2_phie[-10,10]})");
    w.factory("SUM:::model_phie(nsig_phie[0,100000000]*sig_phie, nbkg_phie[0,100000000]*bkg_phie)"); //nsig[0,100000000]*sig2,
    w.var("m_phie")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    RooDataHist dh_phie("dh_phie", "dh_phie", *w.var("m_phie"), Import(*hphie_py));
    RooPlot *fr_phie = w.var("m_phie")->frame(Title(Form("E_{#gamma} = %f",h2phie->GetXaxis()->GetBinCenter(i))));
    w.pdf("model_phie")->fitTo(dh_phie);
    // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    dh_phie.plotOn(fr_phie, RooFit::Name("ldh_phie"));
    w.pdf("model_phie")->plotOn(fr_phie, Components(*w.pdf("sig_phie")), LineColor(kRed), RooFit::Name("lsig_phie"));
    w.pdf("model_phie")->plotOn(fr_phie, Components(*w.pdf("bkg_phie")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phie"));
    w.pdf("model_phie")->plotOn(fr_phie, RooFit::Name("lmodel_phie"));
    // w.pdf("model_phie")->paramOn(fr_phie, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phie"), *w.var("nbkg_phie")))); //,*w.var("mean_phie"),*w.var("width_phie"),*w.var("sigma_phie"))));
    fr_phie->Draw();

    TLegend *l_phie = new TLegend(0.5, 0.7, 0.8, 0.9);
    l_phie->SetFillColor(kWhite);
    l_phie->SetLineColor(kWhite);
    // l_phie->AddEntry(fr_phie->findObject("ldh_phie"), "Data", "p");
    // l_phie->AddEntry(fr_phie->findObject("lmodel_phie"), "total", "l");
    l_phie->AddEntry(fr_phie->findObject("lsig_phie"), Form("N_{Sig} = %.2f", w.var("nsig_phie")->getVal()), "l");
    l_phie->AddEntry(fr_phie->findObject("lbkg_phie"), Form("N_{Bkg} = %.2f", w.var("nbkg_phie")->getVal()), "l");
    l_phie->Draw();

    double N_phie = w.var("nsig_phie")->getVal();
    double dN_phie = w.var("nsig_phie")->getError();

    gphie->SetPoint(i - 1, h2phie->GetXaxis()->GetBinCenter(i), N_phie);
    gphie->SetPointError(i - 1, 0, dN_phie);

    // ++++++++++++++++++++++++++++ mc  +++++++++++++++++++++++
    cphie1_mc->cd(i);
    TH1D *hphie_mc_py = h2phie_mc->ProjectionY(Form("hphie_mc_py_%d",i), i, i);
    hphie_mc_py->Draw("e");

    w.factory(Form("Voigtian::sig_phie_mc(m_phie_mc[%f,%f],mean_phie_mc[1.011,1.030],width_phie_mc[0.004],sigma_phie_mc[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phie_mc[0.001,0.01], mean_phie_mc[1.016,1.022]
    w.factory("Chebychev::bkg_phie_mc(m_phie_mc,{c0_phie_mc[-10,10], c1_phie_mc[-10,10], c2_phie_mc[-10,10]})");
    w.factory("SUM:::model_phie_mc(nsig_phie_mc[0,100000000]*sig_phie_mc, nbkg_phie_mc[0,100000000]*bkg_phie_mc)"); //nsig[0,100000000]*sig2,
    w.var("m_phie_mc")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    RooDataHist dh_phie_mc("dh_phie_mc", "dh_phie_mc", *w.var("m_phie_mc"), Import(*hphie_mc_py));
    RooPlot *fr_phie_mc = w.var("m_phie_mc")->frame(Title(Form("E_{#gamma} = %f",h2phie_mc->GetXaxis()->GetBinCenter(i))));
    w.pdf("model_phie_mc")->fitTo(dh_phie_mc);

    // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    dh_phie_mc.plotOn(fr_phie_mc, RooFit::Name("ldh_phie_mc"));
    w.pdf("model_phie_mc")->plotOn(fr_phie_mc, Components(*w.pdf("sig_phie_mc")), LineColor(kRed), RooFit::Name("lsig_phie_mc"));
    w.pdf("model_phie_mc")->plotOn(fr_phie_mc, Components(*w.pdf("bkg_phie_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phie_mc"));
    w.pdf("model_phie_mc")->plotOn(fr_phie_mc, RooFit::Name("lmodel_phie_mc"));
    // w.pdf("model_phie_mc")->paramOn(fr_phie_mc, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phie_mc"), *w.var("nbkg_phie_mc")))); //,*w.var("mean_phie_mc"),*w.var("width_phie_mc"),*w.var("sigma_phie_mc"))));
    fr_phie_mc->Draw();

    TLegend *l_phie_mc = new TLegend(0.5, 0.7, 0.8, 0.9);
    l_phie_mc->SetFillColor(kWhite);
    l_phie_mc->AddEntry(fr_phie_mc->findObject("lsig_phie_mc"), Form("N_{Sig} = %.2f", w.var("nsig_phie_mc")->getVal()), "l");
    l_phie_mc->AddEntry(fr_phie_mc->findObject("lbkg_phie_mc"), Form("N_{Bkg} = %.2f", w.var("nbkg_phie_mc")->getVal()), "l");
    l_phie_mc->Draw();

    double N_phie_mc = w.var("nsig_phie_mc")->getVal();
    double dN_phie_mc = w.var("nsig_phie_mc")->getError();

    gphie_mc->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), N_phie_mc);
    gphie_mc->SetPointError(i - 1, 0, dN_phie_mc);  

    // ++++++++++++++++++++++++++++ efficiency  +++++++++++++++++++++++
    double eff_phi =  N_phie_mc / h_beame_tru->GetBinContent(i); // Efficiency = N_observed/N_generated
    gphieeff->SetPoint(i - 1, h2phie_mc->GetXaxis()->GetBinCenter(i), eff_phi);
    gphieeff->SetPointError(i - 1, 0, dN_phie_mc/h_beame_tru->GetBinContent(i));//->GetBinError(i)

    // ++++++++++++++++++++++++++++ cross-section  +++++++++++++++++++++++
    double lumi_phi = h_tagged_flux->GetBinContent(i);// * 1.273;   // Luminosity = N_gama * T ,  T = 1.26 barns^-1
    // double br_phi = 0.489;
    //double xsec_phi = 1e9 * N_phie/(eff_phi*lumi_phi*br_phi);
    double xsec_phi = N_phie/lumi_phi;
    //double dxsec_phi = 1e9 * dN_phie/(eff_phi*lumi_phi*br_phi);
    double dxsec_phi = dN_phie/lumi_phi;
    gphiexsec->SetPoint(i - 1, h2phie->GetXaxis()->GetBinCenter(i), xsec_phi);
    gphiexsec->SetPointError(i - 1, 0, dxsec_phi);

    // c2->Update();
    //sleep(1);
    ofs_scanphi <<" i = " << i << " | N_phie = " << N_phie << " | N_phie_mc = " << N_phie_mc << " | h_beame_tru->GetBinContent(i) = " << h_beame_tru->GetBinContent(i) <<" | eff_phi = " << eff_phi << " | xsec_phi = " << xsec_phi << endl;
  
  }

  cout << endl;

  cgphie->cd();
  gphie->Draw("AP");
  cgphie_mc->cd();
  gphie_mc->Draw("AP");
  cgphieeff->cd();
  gphieeff->Draw("AP");
  cgphiexsec->cd();
  gphiexsec->Draw("AP");
  // int j =1;
  // gphie->Write(Form("grphie_%d", j), TObject::kWriteDelete);

  cphie->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie.root", "root");
  cphie->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie.eps", "eps");
  cphie1->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1.root", "root");
  cphie1->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1.eps", "eps"); 
  cgphie->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie.root", "root");
  cgphie->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie.eps", "eps"); 
  cphie_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie_mc.root", "root");
  cphie_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie_mc.eps", "eps");
  cphie1_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1_mc.root", "root");
  cphie1_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphie1_mc.eps", "eps"); 
  cgphie_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie_mc.root", "root");
  cgphie_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphie_mc.eps", "eps");   
  cgphieeff->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphieeff.root", "root");
  cgphieeff->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphieeff.eps", "eps");   
  cgphiexsec->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiexsec.root", "root");
  cgphiexsec->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiexsec.eps", "eps");
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

    cphit1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_%d.root", j), "root");
    cphit1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_%d.eps", j), "eps");
    cphit1_mc[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_mc_%d.root", j), "root");
    cphit1_mc[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit1_mc_%d.eps", j), "eps");
  }
  // int j =1;
  // gphit->Write(Form("grphit_%d", j), TObject::kWriteDelete);

 cbeamet->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cbeamet.root", "root");
  cbeamet->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cbeamet.eps", "eps");
  cphit->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit.root", "root");
  cphit->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit.eps", "eps");
  cgphit->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit.root", "root");
  cgphit->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit.eps", "eps");
  cphit_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit_mc.root", "root");
  cphit_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cphit_mc.eps", "eps");
  cgphit_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit_mc.root", "root");
  cgphit_mc->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphit_mc.eps", "eps");
  cgphiteff->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiteff.root", "root");
  cgphiteff->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphiteff.eps", "eps");
  cgphitxsec->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphitxsec.root", "root");
  cgphitxsec->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/cgphitxsec.eps", "eps");
*/
  c_tagged_flux->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/c_tagged_flux.root", "root");
  c_tagged_flux->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/c_tagged_flux.eps", "eps");
 
  outputfig->Print();
}
