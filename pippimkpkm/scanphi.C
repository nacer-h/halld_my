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

  // TFile *outputfig = new TFile("output/fig_scanphi/scanphi.root", "UPDATE");

  RooWorkspace w("w", kTRUE);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetLabelSize(0.06, "XYZ");

  double m2k_min = 0.99, m2k_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2;   

  TCanvas *cphi=new TCanvas("cphi","cphi",900, 600);//1500, 800

  TCanvas *cphi1 = new TCanvas("cphi1","cphi1", 1500, 1200);//1500, 800
  cphi1->Divide(5, 5);
  TCanvas *cphi2 = new TCanvas("cphi2", "cphi2", 1500, 1200);
  cphi2->Divide(5, 5);

  TCanvas *cgphi=new TCanvas("cgphi","cgphi",900, 600);//1500, 800
  TGraphErrors *gphi;

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  // *********************** Phi(1020) MC *********************************
  TCanvas *c_PhiMass_postcuts = new TCanvas("c_PhiMass_postcuts", "c_PhiMass_postcuts", 1000, 600);
  TH1F *h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", "; m_{K^{+}K^{-}} (GeV/c^{2}); Counts", 200, m2k_min, m2k_max);
  tmc->Project("h_PhiMass_postcuts", "kpkm_mf", "w8*(kpkm_uni)");// && is_truecombo +cutlist+" && kin_chisq<25)"
  c_PhiMass_postcuts->cd();
  h_PhiMass_postcuts->Draw("e");

  // w.factory(Form("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[%f, %f],mean_PhiMass_mc[1.018, 1.021],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.0001,0.01])",m2k_min,m2k_max)); //m_PhiMass_mc[1.0, 1.04],mean_PhiMass_mc[1.016,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.0001,0.01]
  // // w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.0, 1.04],mean_PhiMass_mc[1.016,1.022],width_PhiMass_mc[0.004],sigma_PhiMass_mc[0.0001,0.01])");
  // w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-100000,100000], c1_PhiMass_mc[-100000,100000], c2_PhiMass_mc[-100000,100000]})");//[-1000,1000], , c3_PhiMass_mc[-100000,100000]
  // w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000]*sig_PhiMass_mc, nbkg_PhiMass_mc[0,100000]*bkg_PhiMass_mc)");// [0,100000], nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); //nsig[0,100000000]*sig2,
  // w.var("m_PhiMass_mc")->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
  // RooDataHist dh_PhiMass_mc("dh_PhiMass_mc", "dh_PhiMass_mc", *w.var("m_PhiMass_mc"), Import(*h_PhiMass_postcuts));
  // RooPlot *fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title("K^{+}K^{-}"));
  // // fr_PhiMass_mc->SetTitleOffset(0.90, "X");
  // // fr_PhiMass_mc->SetTitleSize(0.06, "XYZ");
  // // fr_PhiMass_mc->SetLabelSize(0.06, "xyz");
  // w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc);

  // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
  // dh_PhiMass_mc.plotOn(fr_PhiMass_mc, RooFit::Name("ldh_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig_PhiMass_mc")), LineColor(kRed), RooFit::Name("lsig_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("bkg_PhiMass_mc")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass_mc"));
  // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, RooFit::Name("lmodel_PhiMass_mc"));
  // // w.pdf("model_PhiMass_mc")->paramOn(fr_PhiMass_mc, Layout(0.5, 0.90, 0.99));//, Parameters(RooArgSet(*w.var("nsig_PhiMass_mc"), *w.var("nbkg_PhiMass_mc")))); //,*w.var("mean_PhiMass_mc"),*w.var("width_PhiMass_mc"),*w.var("sigma_PhiMass_mc"))));
  // fr_PhiMass_mc->Draw();

  // TLegend *l_phi_mc = new TLegend(0.5, 0.65, 0.6, 0.82);
  // l_phi_mc->SetTextSize(0.05);
  // l_phi_mc->SetFillColor(kWhite);
  // l_phi_mc->SetLineColor(kWhite);
  // // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"), "Data", "p");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lmodel_PhiMass_mc"), "total", "l");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lsig_PhiMass_mc"), "Voigtian", "l");
  // l_phi_mc->AddEntry(fr_PhiMass_mc->findObject("lbkg_PhiMass_mc"), "pol 2nd", "l");
  // l_phi_mc->Draw();

  // TLatex lat_PhiMass_mc;
  // lat_PhiMass_mc.SetTextSize(0.05);
  // lat_PhiMass_mc.SetTextAlign(13); //align at top
  // lat_PhiMass_mc.SetNDC();
  // lat_PhiMass_mc.SetTextColor(kBlue);
  // lat_PhiMass_mc.DrawLatex(0.62, 0.80, Form("N_{Sig} = %0.2f#pm%0.2f", w.var("nsig_PhiMass_mc")->getVal(), w.var("nsig_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.73, Form("N_{Bkg} = %0.2f#pm%0.2f", w.var("nbkg_PhiMass_mc")->getVal(), w.var("nbkg_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.66, Form("#mu = %0.3f#pm%0.3f", w.var("mean_PhiMass_mc")->getVal(), w.var("mean_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.59, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_PhiMass_mc")->getVal(), w.var("width_PhiMass_mc")->getError()));
  // lat_PhiMass_mc.DrawLatex(0.62, 0.52, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_PhiMass_mc")->getVal(), w.var("sigma_PhiMass_mc")->getError()));

  TF1 *fsb_mc = new TF1("fsb_mc", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol4(4)", m2k_min, m2k_max); //4
  fsb_mc->SetLineColor(2);
  fsb_mc->SetParameters(1433, 1.019, 0.003, 0.0042, 1, 1, 1, 1, 1); // voigt
  fsb_mc->SetParLimits(1, 1.015, 1.022);                            //fsb_mc->FixParameter(1, 1.019);
  // fsb_mc->FixParameter(2, 0.0028);   //fsb->SetParLimits(2, 0.0001, 0.01);
  fsb_mc->FixParameter(3, 0.0042); //0.0042   //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

  TF1 *fs_mc = new TF1("fs_mc", "[0]*TMath::Voigt(x - [1], [2], [3])", m2k_min, m2k_max);
  fs_mc->SetLineColor(4);

  TF1 *fb_mc = new TF1("fb_mc", "pol4(4)", m2k_min, m2k_max);
  fb_mc->SetLineColor(28);
  fb_mc->SetLineStyle(2);

  h_PhiMass_postcuts->Fit("fsb_mc", "", "", m2k_min, m2k_max);
  double par_mc[fsb_mc->GetNpar()]; //6
  fsb_mc->GetParameters(&par_mc[0]);
  fs_mc->SetParameters(&par_mc[0]);
  fb_mc->SetParameters(&par_mc[4]); //4

  fs_mc->Draw("same");
  fb_mc->Draw("same");

  double N_phie_mc = fs_mc->Integral(m2k_min, m2k_max) / h_PhiMass_postcuts->GetBinWidth(1);
  double dN_phie_mc = N_phie_mc * fsb_mc->GetParError(0) / fsb_mc->GetParameter(0);

  TLatex lat_phie_mc;
  lat_phie_mc.SetTextSize(0.06);
  lat_phie_mc.SetTextAlign(13); //align at top
  lat_phie_mc.SetNDC();
  lat_phie_mc.SetTextColor(kBlue);
  // lat_phie_mc.DrawLatex(0.5, 0.80, Form("#chi^{2}/NDF = %0.2f", fsb_mc->GetChisquare() / fsb_mc->GetNDF()));
  lat_phie_mc.DrawLatex(0.5, 0.72, Form("N_{sig} = %0.2f#pm%0.2f", N_phie_mc, dN_phie_mc));
  lat_phie_mc.DrawLatex(0.5, 0.64, Form("#mu = %0.3f#pm%0.3f", fsb_mc->GetParameter(1), fsb_mc->GetParError(1)));
  lat_phie_mc.DrawLatex(0.5, 0.56, Form("#sigma = %0.3f#pm%0.3f", fsb_mc->GetParameter(2), fsb_mc->GetParError(2)));
  lat_phie_mc.DrawLatex(0.5, 0.48, Form("#Gamma = %0.3f#pm%0.3f", fsb_mc->GetParameter(3), fsb_mc->GetParError(3)));

  cphi->cd();
  TH2D *h2d;
  if (name == "y")
  {
  h2d = new TH2D("h2d",";m_{K^{+}K^{-}#pi^{+}#pi^{-}} (GeV/c^{2});m_{K^{+}K^{-}} (GeV/c^{2});Counts", n, m2pi2k_min, m2pi2k_max, n2k, m2k_min, m2k_max);
  tdata->Project("h2d", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)");
  // h2d = (TH2D *)fdata->Get("h2_PhiVsYMass_KinFit");
  // cout << " ***** h2d = " << h2d << endl;
  // h2d->RebinX(2);
  }
  if(name=="fo")
  {
  h2d = new TH2D("h2d", ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});m_{K^{+}K^{-}} (GeV/c^{2});Counts", n, m2pi_min, m2pi_max, n2k, m2k_min, m2k_max);
  tdata->Project("h2d", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>6.5 && beam_p4_kin.E()<11.6)");
  // h2d = (TH2D *)fdata->Get("h2_PhiVsfoMass_KinFit");
  // cout << " ***** h2d = " << h2d << endl;
  // h2d->RebinX(2);  
  }
  cout << "h2d = " <<h2d<< endl;

  h2d->Draw("colz");

  gphi = new TGraphErrors(); //n2pi2k
  gphi->SetMarkerStyle(20);
  gphi->SetMarkerSize(1.5);
  gphi->SetMarkerColor(1);
  if(name=="y") gphi->SetTitle(";m_{K^{+}K^{-}#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi#pi^{+}#pi^{-}}");
  if(name=="fo") gphi->SetTitle(";m_{#pi^{+}#pi^{-}} (GeV/c^{2});N_{#phi#pi^{+}#pi^{-}}");

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
    fsb->SetParLimits(1, 1.016, 1.022); // 1.018, 1.021
    fsb->FixParameter(2, fsb_mc->GetParameter(2));   // w.var("sigma_PhiMass_mc")->getVal()
    fsb->FixParameter(3, fsb_mc->GetParameter(3));   // w.var("width_PhiMass_mc")->getVal()

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

    TLatex lat_phiy;
    lat_phiy.SetTextSize(0.06);
    lat_phiy.SetTextAlign(13); //align at top
    lat_phiy.SetNDC();
    lat_phiy.SetTextColor(kBlue);
    lat_phiy.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    lat_phiy.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
    lat_phiy.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    lat_phiy.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
    lat_phiy.DrawLatex(0.45, 0.43, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));

    // hphi_py->Write();
    cgphi->Update();
    // c2->Update();
    //sleep(1);
  }

  // cout << endl;

  cgphi->cd();
  gphi->Draw("AP");
  gphi->SetMinimum(0.);

  c_PhiMass_postcuts->Print(Form("output/fig_scanphi/c_PhiMass_postcuts_%s_%s.root", fname_mc.Data(), fname_data.Data()), "root");
  c_PhiMass_postcuts->Print(Form("output/fig_scanphi/c_PhiMass_postcuts_%s_%s.eps", fname_mc.Data(), fname_data.Data()), "eps");
  c_PhiMass_postcuts->Print(Form("output/fig_scanphi/c_PhiMass_postcuts_%s_%s.png", fname_mc.Data(), fname_data.Data()), "png");
  cphi1->Print(Form("output/fig_scanphi/c1_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi1->Print(Form("output/fig_scanphi/c1_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi1->Print(Form("output/fig_scanphi/c1_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cphi2->Print(Form("output/fig_scanphi/c2_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi2->Print(Form("output/fig_scanphi/c2_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi2->Print(Form("output/fig_scanphi/c2_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cphi->Print(Form("output/fig_scanphi/c_phi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cphi->Print(Form("output/fig_scanphi/c_phi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cphi->Print(Form("output/fig_scanphi/c_phi_%s_%s.png", name.Data(), fname_data.Data()), "png");
  cgphi->Print(Form("output/fig_scanphi/c_gphi_%s_%s.root", name.Data(), fname_data.Data()), "root");
  cgphi->Print(Form("output/fig_scanphi/c_gphi_%s_%s.eps", name.Data(), fname_data.Data()), "eps");
  cgphi->Print(Form("output/fig_scanphi/c_gphi_%s_%s.png", name.Data(), fname_data.Data()), "png");

  // outputfig->Print();

  return gphi;
}
