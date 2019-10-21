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

void scanfo(TString name, int nkk = 100, int n2pi = 100, int ne = 1, int nt = 1) // TString cut="&& kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
   TFile *fdata = NULL;
  if(name == "data_16") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_16_flat.root");
  if(name == "data_17") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_17_flat.root");
  if(name == "data_18") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_18_flat.root");  
  if(name == "data_all") fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_pippimkpkm_all_flat.root");
  TTree *tdata = (TTree *)fdata->Get("ntp");
  TFile *fmc = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/tree_phifo_genr8_17v3_flat.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TFile *fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/scanfo.root", "UPDATE");

  ofstream ofs_scanfo("tableau_fo.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min = 0.99, mkk_max = 1.2;
  double m2pi_min = 0.3, m2pi_max = 1.2;

  // gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


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

  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/cmc_PhiMass_postcuts_fitted.root", "root");
  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/cmc_PhiMass_postcuts_fitted.eps", "eps");
  c_PhiMass_postcuts->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/cmc_PhiMass_postcuts_fitted.png", "png");

  // cout<<"=============================  no problem up to here ! ========================"<<endl;
/*

 // ======================================== fo vs. Phi(1020) ===============================================
 // root -l 'scanfo.C+(100,50,1,1)'
 
  TCanvas *cphifo_all = new TCanvas("cphifo_all", "cphifo_all", 900, 600); // 900, 600
  TCanvas *cphifo_all1 = new TCanvas("cphifo_all1", "cphifo_all1", 1500, 800);
  cphifo_all1->Divide(5, 5);
  TCanvas *cphifo_all2 = new TCanvas("cphifo_all2", "cphifo_all2", 1500, 800);
  cphifo_all2->Divide(5, 5);
  TCanvas *cgphifo_all = new TCanvas("cgphifo_all", "cgphifo_all", 900, 600);
  TGraphErrors *gphifo_all;

  // TCanvas *cgphifo_all_width = new TCanvas("cgphifo_all_width", "cgphifo_all_width", 1500, 400);//
  // cgphifo_all_width->Divide(3,1);
  // TGraphErrors *gphifo_all_width[ne];

  // TCanvas *cgphifo_all_mean = new TCanvas("cgphifo_all_mean", "cgphifo_all_mean", 1500, 400);//
  // cgphifo_all_mean->Divide(3,1);
  // TGraphErrors *gphifo_all_mean[ne];

  // TCanvas *cgnophifo_all=new TCanvas("cgnophifo_all","cgnophifo_all",1500, 800);
  // TGraphErrors *gnophifo_all = new TGraphErrors(n2pi);
  // gnophifo_all->SetMarkerStyle(20);
  // gnophifo_all->SetMarkerSize(1.0);
  // gnophifo_all->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

    cphifo_all->cd();
    TH2F *h2d2 = new TH2F("h2d2", "(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni))");
    h2d2->Draw("colz");

    gphifo_all = new TGraphErrors();//n2pi
    gphifo_all->SetMarkerStyle(20);
    gphifo_all->SetMarkerSize(1.0);
    gphifo_all->SetMarkerColor(1);
    gphifo_all->SetMinimum(0.);
    gphifo_all->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
   
    // gphifo_all_width = new TGraphErrors(n2pi);
    // gphifo_all_width->SetMarkerStyle(20);
    // gphifo_all_width->SetMarkerSize(1.0);
    // gphifo_all_width->SetMarkerColor(1);
    // gphifo_all_width->SetMinimum(0.);
    // gphifo_all_width->SetMaximum(0.01);
    // gphifo_all->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    
    // gphifo_all_mean = new TGraphErrors(n2pi);
    // gphifo_all_mean->SetMarkerStyle(20);
    // gphifo_all_mean->SetMarkerSize(1.0);
    // gphifo_all_mean->SetMarkerColor(1);
    // gphifo_all_mean->SetMinimum(0.);
    // gphifo_all_mean->SetMaximum(0.01);
    // gphifo_all->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    
    // gnophifo_all->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{bkg}", Eg1[j], Eg2[j]));

    for (int i = 1; i <= n2pi; ++i) //n2pi
    {
      // cout << i << " " << flush;

      if (i < 26) cphifo_all1->cd(i); // if (i < 26) cphifo_all1->cd(i);
      if (i > 25) cphifo_all2->cd(i - 25);// if (i > 25 && i < 51) cphifo_all2->cd(i - 25);
     
      TH1D *hphifo_all_py = h2d2->ProjectionY(Form("_hphifo_all_py_%d", i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifo_all_py->Draw();
      // fb2->Draw("same");
 
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);      
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01
      
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max); //pol2(3)
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphifo_all_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N2 = fs->Integral(mkk_min, mkk_max) / hphifo_all_py->GetBinWidth(1);
      double dN2 = N2 * fsb->GetParError(0) / fsb->GetParameter(0);

      gphifo_all->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifo_all->SetPointError(i - 1, 0, dN2);

      // gphifo_all_width[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(2));
      // gphifo_all_width[j]->SetPointError(i - 1, 0, fsb->GetParError(2));

      // gphifo_all_mean[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(3));
      // gphifo_all_mean[j]->SetPointError(i - 1, 0, fsb->GetParError(3));

      TLatex lat_phifo_all;
      lat_phifo_all.SetTextSize(0.09);
      lat_phifo_all.SetTextAlign(13); //align at top
      lat_phifo_all.SetNDC();
      lat_phifo_all.SetTextColor(kBlue);
      lat_phifo_all.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phifo_all.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
      lat_phifo_all.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phifo_all.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phifo_all.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // gnophifo_all->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophifo_all->SetPointError(i - 1, 0, dNbkg);
      // ofs_scanfo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h2d2->GetYaxis()->GetBinCenter(" << i << ") = " << h2d2->GetYaxis()->GetBinCenter(i)<<endl;
      // ofs_scanfo << " i = " << i << " |par0 = " << fsb->GetParameter(0) << " | parerr0 = " << fsb->GetParError(0) << " | par1 " << fsb->GetParameter(1) << " | par2 " << fsb->GetParameter(2) <<endl;
      hphifo_all_py->Write();
      cgphifo_all->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphifo_all->cd();
    gphifo_all->Draw("AP");
    // TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
    TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", 0.85, 1.05);
    fsb->SetLineColor(2);
    fsb->SetParameters(1433, 0.980, 0.1, 1, 1);
    fsb->SetParLimits(0, 0, 10000);
    fsb->SetParLimits(1, 0.95, 0.98); // 1.018, 1.021
    fsb->SetParLimits(2, 0.001, 0.6);        //fsb->SetParLimits(2, 0.008, 0.010);
    // fsb->FixParameter(3, 0.002);        //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01

    // TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
    TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", 0.85, 1.05);
    fs->SetLineColor(4);

    TF1 *fb = new TF1("fb", "pol2(3)", 0.85, 1.05); //pol2(3)
    fb->SetLineColor(28);
    fb->SetLineStyle(2);

    gphifo_all->Fit("fsb", "", "", 0.85, 1.05);
    double par[7]; //6
    fsb->GetParameters(&par[0]);
    fs->SetParameters(&par[0]);
    fb->SetParameters(&par[3]); //3

    double Nfo_data = fs->Integral(0.85, 1.05)*n2pi/(m2pi_max-m2pi_min);
    double dNfo_data = Nfo_data * fsb->GetParError(0) / fsb->GetParameter(0); 

    fs->Draw("same");
    fb->Draw("same");

    TLatex lat_phifo_all;
    lat_phifo_all.SetTextSize(0.05);
    lat_phifo_all.SetTextAlign(13); //align at top
    lat_phifo_all.SetNDC();
    lat_phifo_all.SetTextColor(kBlue);
    lat_phifo_all.DrawLatex(0.6, 0.87, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
    lat_phifo_all.DrawLatex(0.6, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", Nfo_data, dNfo_data));
    lat_phifo_all.DrawLatex(0.6, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
    lat_phifo_all.DrawLatex(0.6, 0.58, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));

    // cgphifo_all_width->cd(j);
    // gphifo_all_width[j]->Draw("AP");

    // cgphifo_all_mean->cd(j);
    // gphifo_all_mean[j]->Draw("AP");


    // cgnophifo_all->cd();//j);
    // gnophifo_all->Draw("Psame");
    // int j =1;
    // gphifo_all->Write(Form("grphifo_all_%d", j), TObject::kWriteDelete);

    cphifo_all1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_all_%s.root", name.Data()), "root");
    cphifo_all1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_all_%s.eps", name.Data()), "eps");
    cphifo_all1->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_all_%s.png", name.Data()), "png");
    cphifo_all2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_all_%s.root", name.Data()), "root");
    cphifo_all2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_all_%s.eps", name.Data()), "eps");
    cphifo_all2->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_all_%s.png", name.Data()), "png");
    cphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo_all_%s.root", name.Data()), "root");
    cphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo_all_%s.eps", name.Data()), "eps");
    cphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo_all_%s.png", name.Data()), "png");
    cgphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_%s.root", name.Data()), "root");
    cgphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_%s.eps", name.Data()), "eps");
    cgphifo_all->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_%s.png", name.Data()), "png");
    // cgphifo_all_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_width.root", "root");
    // cgphifo_all_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_width.eps", "eps");
    // cgphifo_all_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_mean.root", "root");
    // cgphifo_all_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo_all_mean.eps", "eps");

    // cgnophifo_all->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifo_all.root", "root");
    // cgnophifo_all->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifo_all.eps", "eps");
*/
/*
  // ======================================== fo vs. Eg ===============================================
  // root -l 'scanfo.C+(50,50,4,4)'
  // Eg= 3, 7.9, 8.56, 9.66, 12

  TCanvas *cphifoe = new TCanvas("cphifoe", "cphifoe", 1500, 800); // 900, 600
  cphifoe->Divide(3, 2);

  TCanvas *cphifoe1[ne];
  TCanvas *cphifoe2[ne];
  // TCanvas *cphifoe3[ne];
  // TCanvas *cphifoe4[ne];

  TCanvas *cgphifoe = new TCanvas("cgphifoe", "cgphifoe", 1500, 800);//
  cgphifoe->Divide(3,2);
  TGraphErrors *gphifoe[ne];

  // TCanvas *cgphifoe_width = new TCanvas("cgphifoe_width", "cgphifoe_width", 1500, 400);//
  // cgphifoe_width->Divide(3,1);
  // TGraphErrors *gphifoe_width[ne];

  // TCanvas *cgphifoe_mean = new TCanvas("cgphifoe_mean", "cgphifoe_mean", 1500, 400);//
  // cgphifoe_mean->Divide(3,1);
  // TGraphErrors *gphifoe_mean[ne];

  // TCanvas *cgnophifoe=new TCanvas("cgnophifoe","cgnophifoe",1500, 800);
  // TGraphErrors *gnophifoe = new TGraphErrors(n2pi);
  // gnophifoe->SetMarkerStyle(20);
  // gnophifoe->SetMarkerSize(1.0);
  // gnophifoe->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  double Egmin = 6;  //6.3;//= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 12; //11.6;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax - Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    cphifoe->cd(j);
    // TH2F *h2d2 = new TH2F("h2d2", Form("%.2f<E_{#gamma}<%.2f GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    // TH2F *h2d2 = new TH2F("h2d2", "(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    // 4: 6, 8.1, 8.7, 9.3, 12
    // 3: 6, 8.3, 9.6, 12
    TH2F *h2d2;
    if (j == 1)
    {
      h2d2 = new TH2F("h2d2", "3<E_{#gamma}<7.6 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<7.6)");
    }
    if (j == 2)
    {
      h2d2 = new TH2F("h2d2", "7.6<E_{#gamma}<8.2 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>7.6 && beam_p4_kin.E()<8.2)");
    }
    if (j == 3)
    {
      h2d2 = new TH2F("h2d2", "8.2<E_{#gamma}<8.5 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>8.2 && beam_p4_kin.E()<8.5)");
    }
    if (j == 4)
    {
      h2d2 = new TH2F("h2d2", "8.5<E_{#gamma}<8.85 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>8.5 && beam_p4_kin.E()<8.85)");
    }
    if (j == 5)
    {
      h2d2 = new TH2F("h2d2", "8.85<E_{#gamma}<10 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>8.85 && beam_p4_kin.E()<10)");
    }
    if (j == 6)
    {
      h2d2 = new TH2F("h2d2", "10<E_{#gamma}<12 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>10 && beam_p4_kin.E()<12)");
    }
  
    h2d2->Draw("colz");

    cphifoe1[j] = new TCanvas(Form("cphifoe1_%d", j), Form("cphifoe1_%d", j), 1500, 800);
    cphifoe1[j]->Divide(5, 5);
    cphifoe2[j] = new TCanvas(Form("cphifoe2_%d", j), Form("cphifoe2_%d", j), 1500, 800);
    cphifoe2[j]->Divide(5, 5);
    // cphifoe3[j] = new TCanvas(Form("cphifoe3_%d", j), Form("cphifoe3_%d", j), 1500, 800);
    // cphifoe3[j]->Divide(5, 5);
    // cphifoe4[j] = new TCanvas(Form("cphifoe4_%d", j), Form("cphifoe4_%d", j), 1500, 800);
    // cphifoe4[j]->Divide(5, 5);

    gphifoe[j] = new TGraphErrors();//n2pi
    gphifoe[j]->SetMarkerStyle(20);
    gphifoe[j]->SetMarkerSize(1.0);
    gphifoe[j]->SetMarkerColor(1);
    gphifoe[j]->SetMinimum(0.);
    // gphifoe[j]->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==1) gphifoe[j]->SetTitle("3<E_{#gamma}<7.6 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphifoe[j]->SetTitle("7.6<E_{#gamma}<8.2 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphifoe[j]->SetTitle("8.2<E_{#gamma}<8.5 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphifoe[j]->SetTitle("8.5<E_{#gamma}<8.85 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphifoe[j]->SetTitle("8.85<E_{#gamma}<10 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphifoe[j]->SetTitle("10<E_{#gamma}<12 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    // gphifoe_width[j] = new TGraphErrors(n2pi);
    // gphifoe_width[j]->SetMarkerStyle(20);
    // gphifoe_width[j]->SetMarkerSize(1.0);
    // gphifoe_width[j]->SetMarkerColor(1);
    // gphifoe_width[j]->SetMinimum(0.);
    // gphifoe_width[j]->SetMaximum(0.01);
    // // gphifoe[j]->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    // if(j==1) gphifoe_width[j]->SetTitle("6<E_{#gamma}<8.3 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#Gamma [GeV/c^{2}]");
    // if(j==2) gphifoe_width[j]->SetTitle("8.3<E_{#gamma}<9.3 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#Gamma [GeV/c^{2}]");
    // if(j==3) gphifoe_width[j]->SetTitle("9.3<E_{#gamma}<12 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#Gamma [GeV/c^{2}]");

    // gphifoe_mean[j] = new TGraphErrors(n2pi);
    // gphifoe_mean[j]->SetMarkerStyle(20);
    // gphifoe_mean[j]->SetMarkerSize(1.0);
    // gphifoe_mean[j]->SetMarkerColor(1);
    // gphifoe_mean[j]->SetMinimum(0.);
    // gphifoe_mean[j]->SetMaximum(0.01);
    // // gphifoe[j]->SetTitle("(Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    // if(j==1) gphifoe_mean[j]->SetTitle("6<E_{#gamma}<8.3 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#sigma [GeV/c^{2}]");
    // if(j==2) gphifoe_mean[j]->SetTitle("8.3<E_{#gamma}<9.3 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#sigma [GeV/c^{2}]");
    // if(j==3) gphifoe_mean[j]->SetTitle("9.3<E_{#gamma}<12 GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];#sigma [GeV/c^{2}]");

    // gnophifoe->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{bkg}", Eg1[j], Eg2[j]));

    for (int i = 1; i <= n2pi; ++i) //n2pi
    {
      // cout << i << " " << flush;

      if (i < 26) cphifoe1[j]->cd(i);
      if (i > 25 && i < 51) cphifoe2[j]->cd(i - 25);
      // if(i > 50 && i<76) cphifoe3[j]->cd(i-50);
      // if(i > 75) cphifoe4[j]->cd(i-75);

      TH1D *hphifoe_py = h2d2->ProjectionY(Form("_hphifoe_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifoe_py->Draw();
      // fb2->Draw("same");
 
      // // w.factory(Form("Voigtian::sig_phifoe(m_phifoe[%f,%f],mean_phifoe[1.016,1.022],width_phifoe[0.004],sigma_phifoe[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phifo[0.001,0.01], mean_phifo[1.016,1.022],50:mean_phifoe[1.016,1.032]
      // w.factory(Form("BreitWigner::sig_phifoe(m_phifoe[%f,%f],mean_phifoe[1.018,1.021],width_phifoe[0.008, 0.010])", mkk_min, mkk_max));
      // w.factory("Chebychev::bkg_phifoe(m_phifoe,{c0_phifoe[-10,10], c1_phifoe[-10,10], c2_phifoe[-1,1]})");
      // w.factory("SUM:::model_phifoe(nsig_phifoe[0,100000000]*sig_phifoe, nbkg_phifoe[0,100000000]*bkg_phifoe)"); //nsig[0,100000000]*sig2,
      // w.var("m_phifoe")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phifoe("dh_phifoe", "dh_phifoe", *w.var("m_phifoe"), Import(*hphifoe_py));
      // RooPlot *fr_phifoe = w.var("m_phifoe")->frame(Title(" "));
      // // fr_phifoe->SetTitleOffset(0.90, "X");
      // // fr_phifoe->SetTitleSize(0.06, "XYZ");
      // // fr_phifoe->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phifoe")->fitTo(dh_phifoe);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phifoe.plotOn(fr_phifoe, RooFit::Name("ldh_phifoe"));
      // w.pdf("model_phifoe")->plotOn(fr_phifoe, Components(*w.pdf("sig_phifoe")), LineColor(kRed), RooFit::Name("lsig_phifoe"));
      // w.pdf("model_phifoe")->plotOn(fr_phifoe, Components(*w.pdf("bkg_phifoe")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phifoe"));
      // w.pdf("model_phifoe")->plotOn(fr_phifoe, RooFit::Name("lmodel_phifoe"));
      // // w.pdf("model_phifoe")->paramOn(fr_phifoe, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phifoe"), *w.var("nbkg_phifoe")))); //,*w.var("mean_phifoe"),*w.var("width_phifoe"),*w.var("sigma_phifoe"))));
      // fr_phifoe->Draw();

      // TLegend *l_phifoe = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phifoe->SetFillColor(kWhite);
      // l_phifoe->SetLineColor(kWhite);
      // // l_phifoe->AddEntry(fr_phifoe->findObject("ldh_phifoe"), "Data", "p");
      // // l_phifoe->AddEntry(fr_phifoe->findObject("lmodel_phifoe"), "total", "l");
      // l_phifoe->AddEntry(fr_phifoe->findObject("lsig_phifoe"), Form("N_{Sig} = %.2f", w.var("nsig_phifoe")->getVal()), "l");
      // l_phifoe->AddEntry(fr_phifoe->findObject("lbkg_phifoe"), Form("N_{Bkg} = %.2f", w.var("nbkg_phifoe")->getVal()), "l");
      // l_phifoe->Draw();

      // double N2 = w.var("nsig_phifoe")->getVal();
      // double dN2 = w.var("nsig_phifoe")->getError();

      // // // Integrate normalized pdf over subrange
      // // w.var("m_phifoe")->setRange("sigregion",1.005,1.035);
      // // // RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
      // // RooAbsReal* fbkg = w.pdf("bkg_phifoe")->createIntegral(*w.var("m_phifoe"),NormSet(*w.var("m_phifoe")),Range("sigregion"));
      // // RooAbsReal* fsig = w.pdf("sig_phifoe")->createIntegral(*w.var("m_phifoe"),NormSet(*w.var("m_phifoe")),Range("sigregion"));
      // // // double Nbkg = -fbkg->getVal()*(w.var("nsig_phifoe")->getVal()+w.var("nbkg_phifoe")->getVal())+fsig->getVal()*w.var("nsig_phifoe")->getVal();
      // // double Nbkg = fbkg->getVal()*w.var("nbkg_phifoe")->getVal();
      // // double dNbkg = 0;//fbkg->getError();


      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);      
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01
      
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max); //pol2(3)
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphifoe_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N2 = fs->Integral(mkk_min, mkk_max) / hphifoe_py->GetBinWidth(1);
      double dN2 = N2 * fsb->GetParError(0) / fsb->GetParameter(0);

      gphifoe[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifoe[j]->SetPointError(i - 1, 0, dN2);

      // gphifoe_width[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(2));
      // gphifoe_width[j]->SetPointError(i - 1, 0, fsb->GetParError(2));

      // gphifoe_mean[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), fsb->GetParameter(3));
      // gphifoe_mean[j]->SetPointError(i - 1, 0, fsb->GetParError(3));

      TLatex lat_phifoe;
      lat_phifoe.SetTextSize(0.09);
      lat_phifoe.SetTextAlign(13); //align at top
      lat_phifoe.SetNDC();
      lat_phifoe.SetTextColor(kBlue);
      lat_phifoe.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phifoe.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
      lat_phifoe.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phifoe.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phifoe.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // gnophifoe->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophifoe->SetPointError(i - 1, 0, dNbkg);
      // ofs_scanfo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h2d2->GetYaxis()->GetBinCenter(" << i << ") = " << h2d2->GetYaxis()->GetBinCenter(i)<<endl;
      // ofs_scanfo << " i = " << i << " |par0 = " << fsb->GetParameter(0) << " | parerr0 = " << fsb->GetParError(0) << " | par1 " << fsb->GetParameter(1) << " | par2 " << fsb->GetParameter(2) <<endl;
      hphifoe_py->Write();
      cgphifoe->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphifoe->cd(j);
    gphifoe[j]->Draw("AP");
    
    // cgphifoe_width->cd(j);
    // gphifoe_width[j]->Draw("AP");

    // cgphifoe_mean->cd(j);
    // gphifoe_mean[j]->Draw("AP");


    // cgnophifoe->cd();//j);
    // gnophifoe->Draw("Psame");
    // int j =1;
    // gphifoe->Write(Form("grphifoe_%d", j), TObject::kWriteDelete);

    cphifoe1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifoe_%d.root", j), "root");
    cphifoe1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifoe_%d.eps", j), "eps");
    cphifoe1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifoe_%d.png", j), "png");
    cphifoe2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifoe_%d.root", j), "root");
    cphifoe2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifoe_%d.eps", j), "eps");
    cphifoe2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifoe_%d.png", j), "png");
    // cphifoe3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c3_phifoe_%d.root", j), "root");
    // cphifoe3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c3_phifoe_%d.eps", j), "eps");
    // cphifoe4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c4_phifoe_%d.root", j), "root");
    // cphifoe4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c4_phifoe_%d.eps", j), "eps");
  }

  cphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifoe.root", "root");
  cphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifoe.eps", "eps");
  cphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifoe.png", "png");
  cgphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe.root", "root");
  cgphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe.eps", "eps");
  cgphifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe.png", "png");
  // cgphifoe_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe_width.root", "root");
  // cgphifoe_width->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe_width.eps", "eps");
  // cgphifoe_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe_mean.root", "root");
  // cgphifoe_mean->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe_mean.eps", "eps");

  // cgnophifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifoe.root", "root");
  // cgnophifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifoe.eps", "eps");
  // cgnophifoe->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifoe.png", "png");
*/
    /*
  // ======================================== fo vs. -t ===============================================

  TCanvas *cphifot = new TCanvas("cphifot", "cphifot", 10, 10, 1500, 800); //1500, 800
  cphifot->Divide(3,2);

  TCanvas *cphifot1[nt];
  TCanvas *cphifot2[nt];

  TCanvas *cgphifot = new TCanvas("cgphifot", "cgphifot", 10, 10, 1500, 800);
  cgphifot->Divide(3,2);
  TGraphErrors *gphifot[nt];

  double tmin = 0; // 0.1= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double tmax = 10;   // 2= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double tstep = (tmax-tmin) / nt;
  double t1[nt];
  double t2[nt];

  for (int j = 1; j <= nt; ++j)
  {
    t1[j] = tmin + ((j - 1) * tstep); //i * tstep;
    t2[j] = tmin + (j * tstep);
    cout << "########  j = " << j << " | t1[" << j << "] = " << t1[j] << " | t2[" << j << "] = " << t2[j] << endl;

    cphifot->cd(j);
    // TH2F *h2d2 = new TH2F("h2d2", Form("%.2f<-t<%.2f GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    // h2d2->Draw("colz");

    TH2F *h2d2;
    // 4: 
    // 3: 0, 0.7, 1.5, 4
    if (j == 1)
    {
      h2d2 = new TH2F("h2d2", "0<-t<0.45 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>0 && -t_kin<0.45)");
    }
    if (j == 2)
    {
      h2d2 = new TH2F("h2d2", "0.45<-t<0.73 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>0.45 && -t_kin<0.73)");
    }
    if (j == 3)
    {
      h2d2 = new TH2F("h2d2", "0.73<-t<1.1 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>0.73 && -t_kin<1.1)");
    }  
    if (j == 4)
    {
      h2d2 = new TH2F("h2d2", "1.1<-t<1.63 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>1.1 && -t_kin<1.63)");
    }
    if (j == 5)
    {
      h2d2 = new TH2F("h2d2", "1.63<-t<2.63 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>1.63 && -t_kin<2.63)");
    }
    if (j == 6)
    {
      h2d2 = new TH2F("h2d2", "2.63<-t<10 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && -t_kin>2.63 && -t_kin<10)");
    }
  
    h2d2->Draw("colz");

    cphifot1[j] = new TCanvas(Form("cphifot1_%d", j), Form("cphifot1_%d", j), 1500, 800);
    cphifot1[j]->Divide(5, 5);
    cphifot2[j] = new TCanvas(Form("cphifot2_%d", j), Form("cphifot2_%d", j), 1500, 800);
    cphifot2[j]->Divide(5, 5);

    gphifot[j] = new TGraphErrors();//n2pi
    gphifot[j]->SetMarkerStyle(20);
    gphifot[j]->SetMarkerSize(1.0);
    gphifot[j]->SetMarkerColor(1);
    gphifot[j]->SetMinimum(0.);
    // gphifot[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphifot[j]->SetTitle("0<-t<0.45 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphifot[j]->SetTitle("0.45<-t<0.73 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphifot[j]->SetTitle("0.73<-t<1.1 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphifot[j]->SetTitle("1.1<-t<1.63 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphifot[j]->SetTitle("1.63<-t<2.63 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphifot[j]->SetTitle("2.63<-t<10 GeV^{2} (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi; ++i)
    {
      // cout << i << " " << flush;

      if (i < 26)
        cphifot1[j]->cd(i);
      if (i > 25)
        cphifot2[j]->cd(i - 25);

      TH1D *hphifot_py = h2d2->ProjectionY(Form("_hphifot_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifot_py->Draw();
      // fb2->Draw("same");

      // w.factory(Form("Voigtian::sig_phifot(m_phifot[%f,%f],mean_phifot[1.018,1.022],width_phifot[0.004],sigma_phifot[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phifo[0.001,0.01], mean_phifo[1.016,1.022]
      // w.factory("Chebychev::bkg_phifot(m_phifot,{c0_phifot[-1,1], c1_phifot[-1,1]})"); //, c2_phifot[-1,1]
      // w.factory("SUM:::model_phifot(nsig_phifot[0,100000000]*sig_phifot, nbkg_phifot[0,100000000]*bkg_phifot)"); //nsig[0,100000000]*sig2,
      // w.var("m_phifot")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phifot("dh_phifot", "dh_phifot", *w.var("m_phifot"), Import(*hphifot_py));
      // RooPlot *fr_phifot = w.var("m_phifot")->frame(Title(" "));
      // // fr_phifot->SetTitleOffset(0.90, "X");
      // // fr_phifot->SetTitleSize(0.06, "XYZ");
      // // fr_phifot->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phifot")->fitTo(dh_phifot);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phifot.plotOn(fr_phifot, RooFit::Name("ldh_phifot"));
      // w.pdf("model_phifot")->plotOn(fr_phifot, Components(*w.pdf("sig_phifot")), LineColor(kRed), RooFit::Name("lsig_phifot"));
      // w.pdf("model_phifot")->plotOn(fr_phifot, Components(*w.pdf("bkg_phifot")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phifot"));
      // w.pdf("model_phifot")->plotOn(fr_phifot, RooFit::Name("lmodel_phifot"));
      // // w.pdf("model_phifot")->paramOn(fr_phifot, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phifot"), *w.var("nbkg_phifot")))); //,*w.var("mean_phifot"),*w.var("width_phifot"),*w.var("sigma_phifot"))));
      // fr_phifot->Draw();

      // TLegend *l_phifot = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phifot->SetFillColor(kWhite);
      // l_phifot->SetLineColor(kWhite);
      // // l_phifot->AddEntry(fr_phifot->findObject("ldh_phifot"), "Data", "p");
      // // l_phifot->AddEntry(fr_phifot->findObject("lmodel_phifot"), "total", "l");
      // l_phifot->AddEntry(fr_phifot->findObject("lsig_phifot"), Form("N_{Sig} = %.2f", w.var("nsig_phifot")->getVal()), "l");
      // l_phifot->AddEntry(fr_phifot->findObject("lbkg_phifot"), Form("N_{Bkg} = %.2f", w.var("nbkg_phifot")->getVal()), "l");
      // l_phifot->Draw();

      // // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      // double N2 = w.var("nsig_phifot")->getVal();
      // double dN2 = w.var("nsig_phifot")->getError();

  
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);//1433, 1.019, 0.004, 1, 1, 1
      fsb->SetParLimits(0, 0, 10000);
      fsb->SetParLimits(1, 1.015, 1.022); // 1, 1.018, 1.021
      fsb->FixParameter(2, 0.004); // fsb->SetParLimits(3, 0.0008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);

      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//pol2(3)
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphifot_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N2 = fs->Integral(mkk_min, mkk_max) / hphifot_py->GetBinWidth(1);
      double dN2 = N2 * fsb->GetParError(0) / fsb->GetParameter(0);

      gphifot[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifot[j]->SetPointError(i - 1, 0, dN2);

      TLatex lat_phifot;
      lat_phifot.SetTextSize(0.09);
      lat_phifot.SetTextAlign(13); //align at top
      lat_phifot.SetNDC();
      lat_phifot.SetTextColor(kBlue);
      lat_phifot.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phifot.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
      lat_phifot.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phifot.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phifot.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // if(i==1 && j==1) gphifot[j]->RemovePoint(1);
      // cgphifot->Update();
      // c2->Update();
      //sleep(1); 
    }

    // cout << endl;

    cgphifot->cd(j);
    gphifot[j]->Draw("AP");

    // int j =1;
    // gphifot->Write(Form("grphifot_%d", j), TObject::kWriteDelete);

    cphifot1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifot_%d.root", j), "root");
    cphifot1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifot_%d.eps", j), "eps");
    cphifot1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifot_%d.png", j), "png");
    cphifot2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifot_%d.root", j), "root");
    cphifot2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifot_%d.eps", j), "eps");
    cphifot2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifot_%d.png", j), "png");
  
  }

  cphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifot.root", "root");
  cphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifot.eps", "eps");
  cphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifot.png", "png");
  cgphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifot.root", "root");
  cgphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifot.eps", "eps");
  cgphifot->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifot.png", "png");
*/
/*
  // ======================================== fo vs. (E,-t) ===============================================
  TCanvas *cphifo = new TCanvas("cphifo", "cphifo", 10, 10, 1500, 800); //1500, 800
  cphifo->Divide(3,2);

  TCanvas *cphifo1[ne*nt];
  TCanvas *cphifo2[ne*nt];

  TCanvas *cgphifo = new TCanvas("cgphifo", "cgphifo", 10, 10, 1500, 800);
  cgphifo->Divide(3,2);
  TGraphErrors *gphifo[ne*nt];

  for (int j = 1; j <= ne*nt; ++j)
  {
    cout << "########  j = " << j <<endl;

    cphifo->cd(j);
    // TH2F *h2d2 = new TH2F("h2d2", Form("%.2f<-t<%.2f GeV (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    // h2d2->Draw("colz");

    TH2F *h2d2;
    
    if (j == 1)
    {
      h2d2 = new TH2F("h2d2", "(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 2)
    {
      h2d2 = new TH2F("h2d2", "(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 3)
    {
      h2d2 = new TH2F("h2d2", "(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 4)
    {
      h2d2 = new TH2F("h2d2", "(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 5)
    {
      h2d2 = new TH2F("h2d2", "(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 6)
    {
      h2d2 = new TH2F("h2d2", "(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
      tdata->Project("h2d2", "kpkm_mf:pippim_mf", "w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>1.13 && -t_kin<10)");
    }
  
    h2d2->Draw("colz");

    cphifo1[j] = new TCanvas(Form("cphifo1_%d", j), Form("cphifo1_%d", j), 1500, 800);
    cphifo1[j]->Divide(5, 5);
    cphifo2[j] = new TCanvas(Form("cphifo2_%d", j), Form("cphifo2_%d", j), 1500, 800);
    cphifo2[j]->Divide(5, 5);

    gphifo[j] = new TGraphErrors();//n2pi
    gphifo[j]->SetMarkerStyle(20);
    gphifo[j]->SetMarkerSize(1.0);
    gphifo[j]->SetMarkerColor(1);
    gphifo[j]->SetMinimum(0.);
    // gphifo[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphifo[j]->SetTitle("(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphifo[j]->SetTitle("(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphifo[j]->SetTitle("(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphifo[j]->SetTitle("(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphifo[j]->SetTitle("(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphifo[j]->SetTitle("(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi; ++i)
    {
      // cout << i << " " << flush;

      if (i < 26)
        cphifo1[j]->cd(i);
      if (i > 25)
        cphifo2[j]->cd(i - 25);

      TH1D *hphifo_py = h2d2->ProjectionY(Form("_hphifo_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifo_py->Draw();
      // fb2->Draw("same");
  
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);//1433, 1.019, 0.004, 1, 1, 1
      fsb->SetParLimits(0, 0, 10000);   
      fsb->SetParLimits(1, 1.015,1.022);// 1, 1.018, 1.021
      fsb->FixParameter(2, 0.004); // fsb->SetParLimits(3, 0.0008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);

      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//pol2(3)
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphifo_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N2 = fs->Integral(mkk_min, mkk_max) / hphifo_py->GetBinWidth(1);
      double dN2 = N2 * fsb->GetParError(0) / fsb->GetParameter(0);

      gphifo[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifo[j]->SetPointError(i - 1, 0, dN2);

      TLatex lat_phifo;
      lat_phifo.SetTextSize(0.09);
      lat_phifo.SetTextAlign(13); //align at top
      lat_phifo.SetNDC();
      lat_phifo.SetTextColor(kBlue);
      lat_phifo.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phifo.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N2, dN2));
      lat_phifo.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phifo.DrawLatex(0.45, 0.58, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phifo.DrawLatex(0.45, 0.48, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // if(i==1 && j==1) gphifo[j]->RemovePoint(1);
      // cgphifo->Update();
      // c2->Update();
      //sleep(1); 
    }

    // cout << endl;

    cgphifo->cd(j);
    gphifo[j]->Draw("AP");

    // int j =1;
    // gphifo->Write(Form("grphifo_%d", j), TObject::kWriteDelete);

    cphifo1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_%d.root", j), "root");
    cphifo1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_%d.eps", j), "eps");
    cphifo1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifo_%d.png", j), "png");
    cphifo2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_%d.root", j), "root");
    cphifo2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_%d.eps", j), "eps");
    cphifo2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifo_%d.png", j), "png");
  
  }

  cphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo.root", "root");
  cphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo.eps", "eps");
  cphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifo.png", "png");
  cgphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo.root", "root");
  cgphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo.eps", "eps");
  cgphifo->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifo.png", "png");
*/
    outputfig->Print();
}
