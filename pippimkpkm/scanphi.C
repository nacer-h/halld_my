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

TGraphErrors *scanphi(TString name, TString fname_mc="", TString fname_data="", double n2k=100, int n=50)
//, TString cut="kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fmc = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_%s_%s_flat.root", fname_mc.Data(), fname_data.Data()));
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_flat.root", fname_data.Data()));
  TTree *tdata = (TTree *)fdata->Get("ntp");

  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/scanphi.root", "UPDATE");

  RooWorkspace w("w", kTRUE);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  double m2k_min = 0.99, m2k_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2;   

  TCanvas *cphi=new TCanvas("cphi","cphi",900, 600);//1500, 800

  TCanvas *cphi1 = new TCanvas("cphi1","cphi1", 1500, 800);//1500, 800
  cphi1->Divide(5, 5);
  TCanvas *cphi2 = new TCanvas("cphi2", "cphi2", 1500, 800);
  cphi2->Divide(5, 5);

  TCanvas *cgphi=new TCanvas("cgphi","cgphi",900, 600);//1500, 800
  TGraphErrors *gphi;

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  // *********************** Phi(1020) MC *********************************
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", "MC signal; m_{K^{+}K^{-}} [GeV/c^{2}]; Counts", 200, 1.005, 1.035);
  tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni && is_truecombo)");//+cutlist+" && kin_chisq<25)"
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");

  // w.factory("BreitWigner::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.015,1.022],width_PhiMass_mc[0.004])");
  // w.factory("ArgusBG::bkg_Phi_mc(m_Phi_mc, 1.04, c0_Phi_mc[-50,-10])");
  if(name == "y") w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.017,1.021],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.0001,0.01])");
  if(name == "fo") w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005, 1.035],mean_PhiMass_mc[1.017,1.021],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.0001,0.1])");
  w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-100000,100000], c1_PhiMass_mc[-100000,100000]})");
  w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000]*sig_PhiMass_mc, nbkg_PhiMass_mc[0,100000]*bkg_PhiMass_mc)");//, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); //nsig[0,100000000]*sig2,
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

  cphi->cd();
  TH2D *h2d;
  if (name == "y")
  {
  h2d = new TH2D("h2d",";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n, m2pi2k_min, m2pi2k_max, n2k, m2k_min, m2k_max);
  tdata->Project("h2d", "kpkm_mf:kpkmpippim_mf", "w8*(kpkm_uni || kpkmpippim_uni)");
  }
  if(name=="fo")
  {
  h2d = new TH2D("h2d", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n, m2pi_min, m2pi_max, n2k, m2k_min, m2k_max);
  tdata->Project("h2d", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni))");  
  }
  cout << "h2d = " <<h2d<< endl;

  h2d->Draw("colz");

  gphi = new TGraphErrors(); //n2pi2k
  gphi->SetMarkerStyle(20);
  gphi->SetMarkerSize(1.0);
  gphi->SetMarkerColor(1);
  if(name=="y") gphi->SetTitle(";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
  if(name=="fo") gphi->SetTitle(";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

  for (int i = 1; i <= n; ++i) //n2pi2k
  {
    // cout << i << " " << flush;

    if (i < 26)
      cphi1->cd(i); //i
    if (i > 25 && i < 51)
      cphi2->cd(i - 25);
    // if(i > 50 && i<76) cphi3->cd(i-50);
    // if(i > 75) cphi4->cd(i-75);

    TH1D *hphi_py = h2d->ProjectionY(Form("_hphi_py_%d", i), i, i);

    // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", m2k_min, m2k_max);
    TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", m2k_min, m2k_max);
    fsb->SetLineColor(2);
    fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1,1,1,1,1);
    // fsb->SetParLimits(0, 0, 10000);
    fsb->SetParLimits(1, 1.015, 1.022); // 1.018, 1.021
    fsb->FixParameter(2, w.var("sigma_PhiMass_mc")->getVal());   //fsb->SetParLimits(2, 0.008, 0.010); // 0.0042
    fsb->FixParameter(3, w.var("width_PhiMass_mc")->getVal());   //fsb->SetParLimits(3, 0.001,0.01);// 0.002

    // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", m2k_min, m2k_max);
    TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", m2k_min, m2k_max);
    fs->SetLineColor(4);

    TF1 *fb = new TF1("fb", "pol4(4)", m2k_min, m2k_max); //3
    fb->SetLineColor(28);
    fb->SetLineStyle(2);

    hphi_py->Fit("fsb", "", "", m2k_min, m2k_max);
    double par[10]; //6
    fsb->GetParameters(&par[0]);
    fs->SetParameters(&par[0]);
    fb->SetParameters(&par[4]); //3

    fs->Draw("same");
    fb->Draw("same");

    double N1 = fs->Integral(m2k_min, m2k_max) / hphi_py->GetBinWidth(1);
    double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

    // if (N1 <= 0)
    //   gphi->RemovePoint(i - 1);

    gphi->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), N1);
    gphi->SetPointError(i - 1, 0, dN1);

    // gnophiy->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
    // gnophiy->SetPointError(i - 1, 0, dNbkg);
    // ofs_ul_yphi2pi << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

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

    // w.factory(Form("Voigtian::sig_PhiMass(m_PhiMass[%f, %f],mean_PhiMass[1.005,1.035],width_PhiMass[%f],sigma_PhiMass[%f])", m2k_min, m2k_max, w.var("width_PhiMass_mc")->getVal(), w.var("sigma_PhiMass_mc")->getVal())); //sigma_PhiMass[0.001,0.01], mean_PhiMass[1.015,1.022]
    // w.factory("Chebychev::bkg_PhiMass(m_PhiMass,{c0_PhiMass[-10,10], c1_PhiMass[-10,10], c2_PhiMass[-10,10]})");
    // w.factory("SUM:::model_PhiMass(nsig_PhiMass[0,100000000]*sig_PhiMass, nbkg_PhiMass[0,100000000]*bkg_PhiMass)"); //, nbkg_PhiMass[0,100000000]*bkg_PhiMass)"); //nsig[0,100000000]*sig2,
    // w.var("m_PhiMass")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    // RooDataHist dh_PhiMass("dh_PhiMass", "dh_PhiMass", *w.var("m_PhiMass"), Import(*hphi_py));
    // RooPlot *fr_PhiMass = w.var("m_PhiMass")->frame(Title("K^{+}K^{-}"));
    // // fr_PhiMass->SetTitleOffset(0.90, "X");
    // // fr_PhiMass->SetTitleSize(0.06, "XYZ");
    // // fr_PhiMass->SetLabelSize(0.06, "xyz");
    // w.pdf("model_PhiMass")->fitTo(dh_PhiMass);

    // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    // dh_PhiMass.plotOn(fr_PhiMass, RooFit::Name("ldh_PhiMass"));
    // w.pdf("model_PhiMass")->plotOn(fr_PhiMass, Components(*w.pdf("sig_PhiMass")), LineColor(kRed), RooFit::Name("lsig_PhiMass"));
    // w.pdf("model_PhiMass")->plotOn(fr_PhiMass, Components(*w.pdf("bkg_PhiMass")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass"));
    // w.pdf("model_PhiMass")->plotOn(fr_PhiMass, RooFit::Name("lmodel_PhiMass"));
    // // w.pdf("model_PhiMass")->paramOn(fr_PhiMass, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_PhiMass"), *w.var("nbkg_PhiMass")))); //,*w.var("mean_PhiMass"),*w.var("width_PhiMass"),*w.var("sigma_PhiMass"))));
    // fr_PhiMass->Draw();

    // TLegend *l_phi_mc = new TLegend(0.2, 0.65, 0.4, 0.85);
    // l_phi_mc->SetFillColor(kWhite);
    // l_phi_mc->SetLineColor(kWhite);
    // // l_phi_mc->AddEntry(fr_PhiMass->findObject("ldh_PhiMass"), "Data", "p");
    // l_phi_mc->AddEntry(fr_PhiMass->findObject("lmodel_PhiMass"), "total", "l");
    // l_phi_mc->AddEntry(fr_PhiMass->findObject("lsig_PhiMass"), "Voigtian", "l");
    // l_phi_mc->AddEntry(fr_PhiMass->findObject("lbkg_PhiMass"), "pol 2nd", "l");
    // l_phi_mc->Draw();

    // double N1 = w.var("nsig_PhiMass")->getVal();
    // double dN1 = w.var("nsig_PhiMass")->getError();
    // gphi->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), N1);
    // gphi->SetPointError(i - 1, 0, dN1);    

    // TLatex lat_PhiMass;
    // lat_PhiMass.SetTextSize(0.05);
    // lat_PhiMass.SetTextAlign(13); //align at top
    // lat_PhiMass.SetNDC();
    // lat_PhiMass.SetTextColor(kBlue);
    // lat_PhiMass.DrawLatex(0.62, 0.87, Form("N_{Sig} = %0.2f#pm%0.2f", N1, dN1));
    // lat_PhiMass.DrawLatex(0.62, 0.78, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_PhiMass")->getVal(), w.var("nbkg_PhiMass")->getError()));
    // lat_PhiMass.DrawLatex(0.62, 0.68, Form("#mu = %0.3f#pm%0.3f", w.var("mean_PhiMass")->getVal(), w.var("mean_PhiMass")->getError()));
    // lat_PhiMass.DrawLatex(0.62, 0.58, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_PhiMass")->getVal(), w.var("width_PhiMass")->getError()));
    // lat_PhiMass.DrawLatex(0.62, 0.48, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_PhiMass")->getVal(), w.var("sigma_PhiMass")->getError()));

    hphi_py->Write();
    cgphi->Update();
    // c2->Update();
    //sleep(1);
  }

  // cout << endl;

  cgphi->cd();
  gphi->Draw("AP");
  gphi->SetMinimum(0.);

  // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", 1.71, 2.06);
  // //  TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.98, 1.2);
  // fsb->SetLineColor(2);
  // fsb->SetParameters(68, 1.8, 0.09, 1, 1);
  // fsb->SetParLimits(0, 0, 10000);
  // fsb->SetParLimits(1, 1.75, 1.85); // 1.018, 1.021
  // fsb->SetParLimits(2, 0.03, 0.3);  //fsb->SetParLimits(2, 0.008, 0.010);
  //                                   //  fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

  // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", 1.71, 2.06);
  // //  TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
  // fs->SetLineColor(4);

  // TF1 *fb = new TF1("fb", "pol2(3)", 1.71, 2.06); //3
  // fb->SetLineColor(28);
  // fb->SetLineStyle(2);

  // gphi->Fit("fsb", "", "", 1.71, 2.06);
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
  // // lat_phiye.DrawLatex(0.68, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
  // lat_phiye.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));

  c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_PhiMass_postcuts_%s.root", name.Data()), "root");
  c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_PhiMass_postcuts_%s.eps", name.Data()), "eps");
  c_PhiMass_postcuts->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_PhiMass_postcuts_%s.png", name.Data()), "png");
  cphi1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c1_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c1_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c1_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cphi2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c2_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c2_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c2_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cgphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_gphi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cgphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_gphi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cgphi->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanphi/c_gphi_%s_%s.png", name.Data(), fname_data.Data()), "png");

  outputfig->Print();

  return gphi;
}
