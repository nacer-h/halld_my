//  File: pimpipkmkp.C # Author: Nacer # Date: 29 January 2017 # Email: a.hamdi@gsi.de # Description: Macro to study Y(2175).

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TMath.h"
#include "TStyle.h"
#include "vector"
#include "TF1.h"
#include "TFormula.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "iostream"
#include "iomanip"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TPaveStats.h"
#include "TSpectrum.h"
using namespace std;

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

void significance(TString cut)
{
  // TFile *fsim = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/pippimkpkm__B4_sim_recon17_1v03_%scut.root",cut.Data()));
  TFile *fsim = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_sig_17v03_%scut.root",cut.Data()));
  TFile *fdata = new TFile(Form("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm__B4_17v21_%scut.root",cut.Data()));
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_significance/significance.root","UPDATE");

  // gStyle->SetLabelSize(0.07,"xyz"); // size of axis values
  // gStyle->SetTitleSize(0.07,"XYZ"); // size of axis titles
  // gStyle->SetTitleOffset(0.90,"X");

  
  // RooWorkspace w("w",kTRUE);
  // RooFitResult* result = NULL;
  double cuts[15] = {0.06, 0.05, 0.04, 0.03, 0.02, 0.018, 0.016, 0.014, 0.012, 0.010, 0.008, 0.006, 0.004, 0.002, 0.0008};
  // double cuts_chi2[15] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 4, 2, 1, 0.5};

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  // ******* Selection variable
  TCanvas *cdata_precut = new TCanvas("cdata_precut","cdata_precut",600,400);
  cdata_precut->cd();
  TH1F *hdata_precut = (TH1F*) fdata->Get(Form("h_%s_precut",cut.Data()));
  cout<<"hdata_precut = "<<hdata_precut<<endl;
  hdata_precut->Draw();
  cdata_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_%s_precut.root",cut.Data()), "root");
  cdata_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_%s_precut.eps",cut.Data()), "eps");

  TCanvas *cdata_postcut = new TCanvas("cdata_postcut","cdata_postcut",600,400);
  cdata_postcut->cd();
  TH1F *hdata_postcut = (TH1F*) fdata->Get(Form("h_%s_postcut",cut.Data()));
  cout<<"hdata_postcut = "<<hdata_postcut<<endl;
  hdata_postcut->Rebin(5);
  hdata_postcut->Draw();
  cdata_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_%s_postcut.root",cut.Data()), "root");
  cdata_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_%s_postcut.eps",cut.Data()), "eps");

  // ********************************* Phi(1020) *********************************
    
  TCanvas *cdata_PhiMass = new TCanvas("cdata_PhiMass","cdata_PhiMass",1500,800);
  cdata_PhiMass->Divide(5,3);
  // cdata_PhiMass->Divide(1,3);
  TH1F *hdata_PhiMass;
 
  // TCanvas *cgrdata = new TCanvas("cgrdata","cgrdata",600,400);
  // cgrdata->SetGrid();
  // TGraphErrors *grdata = new TGraphErrors();
  // grdata->SetMarkerStyle(20);
  // grdata->GetXaxis()->SetNdivisions(510, kFALSE);
  
  for(int j=0; j<15; j++)
    {
      cdata_PhiMass->cd(j+1);
      hdata_PhiMass = (TH1F*) fdata->Get(Form("h_PhiMass_%scut_%d",cut.Data(),j));
      hdata_PhiMass->Rebin(5);
      double xmin = hdata_PhiMass->GetXaxis()->GetXmin();
      double xmax = hdata_PhiMass->GetXaxis()->GetXmax();
      hdata_PhiMass->SetMarkerStyle(20);
      hdata_PhiMass->SetMarkerSize(0.7);
      cout<<"hdata_PhiMass = "<<hdata_PhiMass<<endl;
      hdata_PhiMass->Draw("e");
      
      // //w.factory("Voigtian::sig(m[0.99,1.05],mean[1.017,1.021],width[0.001,0.01],sigma[0.001,0.01])");
      // w.factory("BreitWigner::sig(m[0.99,1.05],mean[1.017,1.021],width[0.0001,0.01])");
      // w.factory("Chebychev::bkg(m,{c0[-100,100],c1[-1,1],c2[-1,1]})");
      // w.factory("SUM:::model(nsig[0,100000000]*sig, nbkg[0,100000000]*bkg)"); //nsig[0,100000000]*sig2,
      // w.var("m")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_PhiMass("dh_PhiMass","dh_PhiMass",*w.var("m"),Import(*hdata_PhiMass));
      // RooPlot* fr_PhiMass = w.var("m")->frame(Title("invariant mass K^{+}K^{-}"));
      // fr_PhiMass->SetTitleOffset(0.90,"X");
      // fr_PhiMass->SetTitleSize(0.06,"XYZ");
      // fr_PhiMass->SetLabelSize(0.06,"xyz");
      // w.pdf("model")->fitTo(dh_PhiMass);
      // //cout<<"=============================  no problem up to here ! ========================"<<endl;

      // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_PhiMass.plotOn(fr_PhiMass,RooFit::Name("ldh_PhiMass"));
      // w.pdf("model")->plotOn(fr_PhiMass, Components(*w.pdf("sig")),LineColor(kRed),RooFit::Name("lsig"));
      // w.pdf("model")->plotOn(fr_PhiMass, Components(*w.pdf("bkg")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg"));
      // w.pdf("model")->plotOn(fr_PhiMass,RooFit::Name("lmodel"));
      // w.pdf("model")->paramOn(fr_PhiMass,Layout(0.55,0.95,0.6));//,Parameters(RooArgSet(*w.var("nsig"),*w.var("nbkg"),*w.var("mean"),*w.var("width"),*w.var("sigma"))));
      // fr_PhiMass->Draw();
      
      // TLegend *l_PhiMass = new TLegend(0.13,0.6,0.3,0.85);
      // l_PhiMass->SetFillColor(kWhite);
      // l_PhiMass->SetLineColor(kWhite);
      // l_PhiMass->AddEntry(fr_PhiMass->findObject("ldh_PhiMass"),"Data", "p");
      // l_PhiMass->AddEntry(fr_PhiMass->findObject("lmodel"),"total","l");
      // l_PhiMass->AddEntry(fr_PhiMass->findObject("lbkg"),"Background", "l");
      // l_PhiMass->AddEntry(fr_PhiMass->findObject("lsig"),"Signal", "l");
      // l_PhiMass->Draw();
      
      // double ndsig = w.var("nsig")->getVal();
      // double ndsigerr = w.var("nsig")->getError();
      // double ndbkg = w.var("nbkg")->getVal();
      // double ndbkgerr = w.var("nbkg")->getError();
      // double zd = ndsig/sqrt(ndsig+ndbkg);
      // double zderr = sqrt(((ndsig+(2*ndbkg))*(ndsig+(2*ndbkg))*(ndsigerr*ndsigerr) + (ndsig*ndsig)+(ndbkgerr*ndbkgerr))/(4*(ndsig+ndbkg)*(ndsig+ndbkg)*(ndsig+ndbkg)));
      // cout<<"ndsig = "<<ndsig<<"| ndbkg = "<<ndbkg<<"| zd = "<<zd<<endl;
      // // hdata_PhiMass->Draw("e");
      // hdata_PhiMass->Write();	       
      
      // grdata->SetPoint(j,cuts[j],zd);
      // grdata->SetPointError(j,0,zderr);
    }

  //result->Print();
  // h_PhiMass->SetMarkerStyle(20);
  // h_PhiMass->SetMarkerSize(0.7);
  // h_PhiMass->Draw("e");
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%scut.root",cut.Data()), "root");
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%scut.eps",cut.Data()), "eps");

  // cgrdata->cd();
  // grdata->Draw("AP");
  
  // TF1 *fdatapol6 = new TF1("fdatapol6","pol6");
  // grdata->Fit("fdatapol6","mR","", cuts[12], cuts[4]);
  // double zd_max = fdatapol6->GetMaximumX(cuts[12], cuts[4]);
  // cout<<"zd_max = "<<zd_max<<endl;

  // grdata->Write();
  // grdata->SetTitle("Significance #phi(1020) Vs selection (data); cut; Significance");
  // cgrdata->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_%scut.root",cut.Data()),"root");
  // cgrdata->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_%scut.eps",cut.Data()),"eps");
  

  // ********************************* fo(980) *********************************

  TCanvas *cdata_foMass = new TCanvas("cdata_foMass","cdata_foMass",1500,800);
  cdata_foMass->Divide(5,3);
  // cdata_foMass->Divide(1,3);
  TH1F *hdata_foMass;

  for(int j=0; j<15; j++)
    {
      cdata_foMass->cd(j+1);
      hdata_foMass = (TH1F*) fdata->Get(Form("h_foMass_%scut_%d",cut.Data(),j));
      hdata_foMass->Rebin(5);
      double xmin = hdata_foMass->GetXaxis()->GetXmin();
      double xmax = hdata_foMass->GetXaxis()->GetXmax();
      hdata_foMass->SetMarkerStyle(20);
      hdata_foMass->SetMarkerSize(0.7);
      cout<<"hdata_foMass = "<<hdata_foMass<<endl;
      hdata_foMass->Draw("e");
    }

  cdata_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_foMass_%scut.root",cut.Data()), "root");
  cdata_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_foMass_%scut.eps",cut.Data()), "eps");

  // ********************************* Y(2175) *********************************

  TCanvas *cdata_YMass = new TCanvas("cdata_YMass","cdata_YMass",1500,800);
  cdata_YMass->Divide(5,3);
  // cdata_YMass->Divide(1,3);
  TH1F *hdata_YMass;

  for(int j=0; j<15; j++)
    {
      cdata_YMass->cd(j+1);
      hdata_YMass = (TH1F*) fdata->Get(Form("h_YMass_%scut_%d",cut.Data(),j));
      cout<<"hdata_YMass = "<<hdata_YMass<<endl;
      hdata_YMass->Rebin(5);
      double xmin = hdata_YMass->GetXaxis()->GetXmin();
      double xmax = hdata_YMass->GetXaxis()->GetXmax();
      hdata_YMass->SetMarkerStyle(20);
      hdata_YMass->SetMarkerSize(0.7);
      hdata_YMass->Draw("e");

      // double pardata_y[3];
      // TF1 *fsigdata_y = new TF1("fsigdata_y","gaus"); //[0]*TMath::BreitWigner(x,[1],[2])
      // fsigdata_y->SetLineStyle(9);
      // fsigdata_y->SetLineColor(kRed);
      // fsigdata_y->SetParameters(16000,2.175,0.04);
      // hdata_YMass->Fit(fsigdata_y,"","mR",1.8,2.6);
      // fsig_y->GetParameters(&par_y[0]);
      // TF1 *fbkg_y = (TF1*) gROOT->GetFunction("chebyshev4");
          
      // TF1 *fbkg_phiy = new TF1("fbkg_phiy","pol4");
      // fbkg_phiy->SetLineStyle(9);
      // fbkg_phiy->SetLineColor(4);
      // fbkg_phiy->SetParameters(1.,1.,1.,1.,1.);
      // grphiy->Fit(fbkg_phiy,"","mR",1.46,3.20);
      // fbkg_phiy->GetParameters(&paramphiy[3]);
      // TF1 *fmodel_phiy = new TF1("fmodel_phiy","[0]*TMath::BreitWigner(x,[1],[2])+pol4(3)");
      // fmodel_phiy->SetParameters(&paramphiy[0]);
      // fmodel_phiy->SetParLimits(1,1.60,1.80);
      // // fmodel_phifo->SetParLimits(4,);
      // grphiy->Fit(fmodel_phiy,"","mR",1.46,3.20);
    }

  cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%scut.root",cut.Data()), "root");
  cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%scut.eps",cut.Data()), "eps");

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    MC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  // ******* Selection variable
  TCanvas *csim_precut = new TCanvas("csim_precut","csim_precut",600,400);
  csim_precut->cd();
  TH1F *hsim_precut = (TH1F*) fsim->Get(Form("h_%s_precut",cut.Data()));
  cout<<"hsim_precut = "<<hsim_precut<<endl;
  hsim_precut->Draw();
  csim_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_%s_precut.root",cut.Data()), "root");
  csim_precut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_%s_precut.eps",cut.Data()), "eps");

  TCanvas *csim_postcut = new TCanvas("csim_postcut","csim_postcut",600,400);
  csim_postcut->cd();
  TH1F *hsim_postcut = (TH1F*) fsim->Get(Form("h_%s_postcut",cut.Data()));
  cout<<"hsim_postcut = "<<hsim_postcut<<endl;
  hsim_postcut->Draw();
  csim_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_%s_postcut.root",cut.Data()), "root");
  csim_postcut->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_%s_postcut.eps",cut.Data()), "eps");

 // ********************************* Phi(1020) *********************************

  TCanvas *csim_PhiMass = new TCanvas("csim_PhiMass","csim_PhiMass",1500,800);
  csim_PhiMass->Divide(5,3);
  // csim_PhiMass->Divide(1,3);
  TH1F *hsim_PhiMass;

  // TCanvas *cgrsim_PhiMass = new TCanvas("cgrsim_PhiMass","cgrsim_PhiMass",600,400);
  // cgrsim_PhiMass->SetGrid();
  // TGraphErrors *grsim_PhiMass = new TGraphErrors();
  // grsim_PhiMass->SetMarkerStyle(20);
  // grsim_PhiMass->GetXaxis()->SetNdivisions(510, kFALSE);
    
  for(int j=0; j<15; j++)
    {
      csim_PhiMass->cd(j+1);
      hsim_PhiMass = (TH1F*) fsim->Get(Form("h_PhiMass_%scut_%d",cut.Data(),j));
      hsim_PhiMass->Rebin(5);
      double xmin = hsim_PhiMass->GetXaxis()->GetXmin();
      double xmax = hsim_PhiMass->GetXaxis()->GetXmax();
      hsim_PhiMass->SetMarkerStyle(20);
      hsim_PhiMass->SetMarkerSize(0.7);
      cout<<"hsim_PhiMass = "<<hsim_PhiMass<<endl;
      hsim_PhiMass->Draw("e");
      
      // //w.factory("Voigtian::sig1_mc(m_mc[0.99,1.05],mean1_mc[1.017,1.021],width1_mc[0.0004,0.04],sigma1_mc[0.0004,0.04])");
      // //w.factory("Gaussian::sig1_PhiMass_mc(m_PhiMass_mc[0.99,1.05],mean1_PhiMass_mc[1.019,1.017,1.021],sigma1_PhiMass_mc[0.0001,0.02])");
      // //w.factory("Gaussian::sig2_PhiMass_mc(m_PhiMass_mc[0.99,1.05],mean2_PhiMass_mc[1.019,1.017,1.021],sigma2_PhiMass_mc[0.0001,0.02])");
      // //w.factory("Voigtian::sig1_PhiMass_mc(m_PhiMass_mc[0.99,1.05],mean1_PhiMass_mc[1.017,1.021],width1_PhiMass_mc[0.0001,0.1],sigma1_PhiMass_mc[0.0001,0.1])");
      // w.factory("BreitWigner::sig1_PhiMass_mc(m_PhiMass_mc[0.99,1.05],mean1_PhiMass_mc[1.017,1.021],width1_PhiMass_mc[0.0001,0.1])");//mean1_PhiMass_mc[1.017,1.021],width1_PhiMass_mc[0.0001,0.1]
      // w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-10,10],c1_PhiMass_mc[-10,10],c2_PhiMass_mc[-10,10]})");
      // //w.factory("SUM:::model(nsig[0,100000000]*sig, nbkg[0,100000000]*bkg)"); //nsig[0,100000000]*sig2,
      // w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000000]*sig1_PhiMass_mc, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); // nsig_PhiMass_mc[0,100000000]*sig2_PhiMass_mc,
      // w.var("m_PhiMass_mc")->SetTitle("m_PhiMass_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_PhiMass_mc("dh_PhiMass_mc","dh_PhiMass_mc",*w.var("m_PhiMass_mc"),Import(*hsim_PhiMass));
      // RooPlot* fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
      // fr_PhiMass_mc->SetTitleOffset(0.90,"X");
      // fr_PhiMass_mc->SetTitleSize(0.06,"XYZ");
      // fr_PhiMass_mc->SetLabelSize(0.06,"xyz");
      // w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc);
      // // result = w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc,Extended(kTRUE),Save());
      // dh_PhiMass_mc.plotOn(fr_PhiMass_mc,RooFit::Name("ldh_PhiMass_mc"));
      // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig1_PhiMass_mc")),LineColor(kRed),RooFit::Name("lsig1_PhiMass_mc"));
      // //w.pdf("model_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig2_mc")),LineColor(kMagenta),RooFit::Name("lsig2_mc"));
      // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("bkg_PhiMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_PhiMass_mc"));
      // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc,RooFit::Name("lmodel_PhiMass_mc"));
      // w.pdf("model_PhiMass_mc")->paramOn(fr_PhiMass_mc,Layout(0.55,0.98,0.92));//,Parameters(RooArgSet(*w.var("nsig_mc"),*w.var("nbkg_mc"),*w.var("mean1_mc"),*w.var("mean2_mc"),*w.var("sigma1_mc"),*w.var("sigma2_mc"))));
      // fr_PhiMass_mc->Draw();

      // TLegend *l_PhiMass_mc = new TLegend(0.13,0.6,0.3,0.85);
      // l_PhiMass_mc->SetFillColor(kWhite);
      // l_PhiMass_mc->SetLineColor(kWhite);
      // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"),"Data", "p");
      // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lmodel_PhiMass_mc"),"total","l");
      // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lbkg_PhiMass_mc"),"Background", "l");
      // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lsig1_PhiMass_mc"),"Signal1", "l");
      // // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lsig2_PhiMass_mc"),"Signal2", "l");
      // l_PhiMass_mc->Draw();
      
      // double fscale = 15.31; // fscale = (ndsig/ndbkg)/(nmcsig/nmcbkg) 
      // double nmcsig = w.var("nsig_PhiMass_mc")->getVal();
      // double nmcsigerr = w.var("nsig_PhiMass_mc")->getError();
      // double nmcbkg = w.var("nbkg_PhiMass_mc")->getVal();
      // double nmcbkgerr = w.var("nbkg_PhiMass_mc")->getError();
      // double zmc = nmcsig/sqrt(nmcsig+(fscale*nmcbkg));
      // double zmcerr = sqrt(((nmcsig+(2*fscale*nmcbkg))*(nmcsig+(2*fscale*nmcbkg))*(nmcsigerr*nmcsigerr) + (fscale*fscale)*(nmcsig*nmcsig)+(nmcbkgerr*nmcbkgerr))/(4*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))));
      // cout<<"nmcsig = "<<nmcsig<<"| nmcbkg = "<<nmcbkg<<"| zmc = "<<zmc<<endl;
      // // hsim_PhiMass->Draw("e");
      
      // // double zObs = NumberCountingUtils::BinomialObsZ(nmcsig+nmcbkg, nmcbkg, nmcbkgerr);
      // // cout << " Z value (Gaussian sigma) = "<< zObs << endl;
      
      // grsim_PhiMass->SetPoint(j,cuts[j],zmc);
      // grsim_PhiMass->SetPointError(j,0,zmcerr);
      // hsim_PhiMass->Write();	       
  }

  csim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_PhiMass_%scut.root",cut.Data()), "root");
  csim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_PhiMass_%scut.eps",cut.Data()), "eps");

  // cgrsim_PhiMass->cd();
  // grsim_PhiMass->Draw("AP");

  // TF1 *fsimpol6 = new TF1("fsimpol6","pol6");
  // grsim_PhiMass->Fit("fsimpol6","mR","", cuts[12], cuts[4]);
  // double zmc_max = fsimpol6->GetMaximumX(cuts[12], cuts[4]);
  // cout<<"zmc_max = "<<zmc_max<<endl;

  // //result->Print();
  // grsim_PhiMass->Write();
  // grsim_PhiMass->SetTitle("Significance #phi(1020) Vs selection (MC); cut; Significance");
  // cgrsim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_PhiMass_%scut.root",cut.Data()),"root");
  // cgrsim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_PhiMass_%scut.eps",cut.Data()),"eps");
  
  // for (int i = 0; i < nbins; i++) source[i]=hsim_PhiMass->GetBinContent(i+1);
  // s->Background(source,nbins,70,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
  // hsim_PhiMass_bkg = new TH1F("hsim_PhiMass_bkg","",nbins,xmin,xmax);
  // for (int i = 0; i < nbins; i++) hsim_PhiMass_bkg->SetBinContent(i+1,source[i]);

 // ********************************* fo(980) *********************************

  TCanvas *csim_foMass = new TCanvas("csim_foMass","csim_foMass",1500,800);
  csim_foMass->Divide(5,3);
  //csim_foMass->Divide(1,3);
  TH1F *hsim_foMass;

  // TCanvas *cgrsim_foMass = new TCanvas("cgrsim_foMass","cgrsim_foMass",600,400);
  // cgrsim_foMass->SetGrid();
  // TGraphErrors *grsim_foMass = new TGraphErrors();
  // grsim_foMass->SetMarkerStyle(20);
  // grsim_foMass->GetXaxis()->SetNdivisions(510, kFALSE);
    
  for(int j=0; j<15; j++)
    {
      csim_foMass->cd(j+1);
      hsim_foMass = (TH1F*) fsim->Get(Form("h_foMass_%scut_%d",cut.Data(),j));
      hsim_foMass->Rebin(5);
      double xmin = hsim_foMass->GetXaxis()->GetXmin();
      double xmax = hsim_foMass->GetXaxis()->GetXmax();
      hsim_foMass->SetMarkerStyle(20);
      hsim_foMass->SetMarkerSize(0.7);
      cout<<"hsim_foMass = "<<hsim_foMass<<endl;
      hsim_foMass->Draw("e");

      // double par_fo[3];
      // TF1 *fsig_fo = new TF1("fsig_fo","gaus"); //[0]*TMath::BreitWigner(x,[1],[2])
      // fsig_fo->SetLineStyle(9);
      // // fsig_fo->SetLineColor(kRed);
      // fsig_fo->SetParameters(15000,0.980,0.03);
      // hsim_foMass->Fit(fsig_fo,"","mR",0.6,1.4);

      // //w.factory("BreitWigner::sig1_foMass_mc(m_foMass_mc[0.7,1.2],mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.005,0.08])");//mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.005,0.08]
      // w.factory("Voigtian::sig1_foMass_mc(m_foMass_mc[0.7,1.2],mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.00005,0.08],sigma1_foMass_mc[0.005,0.08])");
      // w.factory("Chebychev::bkg_foMass_mc(m_foMass_mc,{c0_foMass_mc[-10,10],c1_foMass_mc[-10,10],c2_foMass_mc[-10,10]})");
      // //w.factory("SUM:::model(nsig[0,100000000]*sig, nbkg[0,100000000]*bkg)"); //nsig[0,100000000]*sig2,
      // w.factory("SUM:::model_foMass_mc(nsig_foMass_mc[0,100000000]*sig1_foMass_mc, nbkg_foMass_mc[0,100000000]*bkg_foMass_mc)"); //, nsig_mc[0,100000000]*sig2_mc
      // w.var("m_foMass_mc")->SetTitle("m_foMass_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_foMass_mc("dh_foMass_mc","dh_foMass_mc",*w.var("m_foMass_mc"),Import(*hsim_foMass));
      // RooPlot* fr_foMass_mc = w.var("m_foMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
      // fr_foMass_mc->SetTitleOffset(0.90,"X");
      // fr_foMass_mc->SetTitleSize(0.06,"XYZ");
      // fr_foMass_mc->SetLabelSize(0.06,"xyz");
      // w.pdf("model_foMass_mc")->fitTo(dh_foMass_mc);
      // // result = w.pdf("model_foMass_mc")->fitTo(dh_foMass_mc,Extended(kTRUE),Save());
      // dh_foMass_mc.plotOn(fr_foMass_mc,RooFit::Name("ldh_foMass_mc"));
      // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("sig1_foMass_mc")),LineColor(kRed),RooFit::Name("lsig1_foMass_mc"));
      // //w.pdf("model_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("sig2_mc")),LineColor(kMagenta),RooFit::Name("lsig2_mc"));
      // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("bkg_foMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_foMass_mc"));
      // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc,RooFit::Name("lmodel_foMass_mc"));
      // w.pdf("model_foMass_mc")->paramOn(fr_foMass_mc,Layout(0.55,0.98,0.92));//,Parameters(RooArgSet(*w.var("nsig_mc"),*w.var("nbkg_mc"),*w.var("mean1_mc"),*w.var("mean2_mc"),*w.var("sigma1_mc"),*w.var("sigma2_mc"))));
      // fr_foMass_mc->Draw();

      // TLegend *l_foMass_mc = new TLegend(0.13,0.6,0.3,0.85);
      // l_foMass_mc->SetFillColor(kWhite);
      // l_foMass_mc->SetLineColor(kWhite);
      // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("ldh_foMass_mc"),"Data", "p");
      // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lmodel_foMass_mc"),"total","l");
      // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lbkg_foMass_mc"),"Background", "l");
      // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lsig1_foMass_mc"),"Signal1", "l");
      // // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lsig2_foMass_mc"),"Signal2", "l");
      // l_foMass_mc->Draw();
      
      // double fscale = 15.31; // fscale = (ndsig/ndbkg)/(nmcsig/nmcbkg) 
      // double nmcsig = w.var("nsig_foMass_mc")->getVal();
      // double nmcsigerr = w.var("nsig_foMass_mc")->getError();
      // double nmcbkg = w.var("nbkg_foMass_mc")->getVal();
      // double nmcbkgerr = w.var("nbkg_foMass_mc")->getError();
      // double zmc = nmcsig/sqrt(nmcsig+(fscale*nmcbkg));
      // double zmcerr = sqrt(((nmcsig+(2*fscale*nmcbkg))*(nmcsig+(2*fscale*nmcbkg))*(nmcsigerr*nmcsigerr) + (fscale*fscale)*(nmcsig*nmcsig)+(nmcbkgerr*nmcbkgerr))/(4*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))));
      // cout<<"nmcsig = "<<nmcsig<<"| nmcbkg = "<<nmcbkg<<"| zmc = "<<zmc<<endl;
      // // hsim_foMass->Draw("e");
      
      // // double zObs = NumberCountingUtils::BinomialObsZ(nmcsig+nmcbkg, nmcbkg, nmcbkgerr);
      // // cout << " Z value (Gaussian sigma) = "<< zObs << endl;
      
      // grsim_foMass->SetPoint(j,cuts[j],zmc);
      // grsim_foMass->SetPointError(j,0,zmcerr);
      // hsim_foMass->Write();	       
  }

  csim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_foMass_%scut.root",cut.Data()), "root");
  csim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_foMass_%scut.eps",cut.Data()), "eps");

  // cgrsim_foMass->cd();
  // grsim_foMass->Draw("AP");

  // // TF1 *fsimpol6 = new TF1("fsimpol6","pol6");
  // // grsim_foMass->Fit("fsimpol6","mR","", cuts[12], cuts[4]);
  // // double zmc_max = fsimpol6->GetMaximumX(cuts[12], cuts[4]);
  // // cout<<"zmc_max = "<<zmc_max<<endl;

  // // result->Print();
  // grsim_foMass->Write();
  // grsim_foMass->SetTitle("Significance #f_{0}(980) Vs selection (MC); cut; Significance");
  // cgrsim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_foMass_%scut.root",cut.Data()),"root");
  // cgrsim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_foMass_%scut.eps",cut.Data()),"eps");

 // ********************************* Y(2175) *********************************

  TCanvas *csim_YMass = new TCanvas("csim_YMass","csim_YMass",1500,800);
  csim_YMass->Divide(5,3);
  // csim_YMass->Divide(1,3);
  TH1F *hsim_YMass;

  // TCanvas *cgrsim_YMass = new TCanvas("cgrsim_YMass","cgrsim_YMass",600,400);
  // cgrsim_YMass->SetGrid();
  // TGraphErrors *grsim_YMass = new TGraphErrors();
  // grsim_YMass->SetMarkerStyle(20);
  // grsim_YMass->GetXaxis()->SetNdivisions(510, kFALSE);
    
  for(int j=0; j<15; j++)
    {
      csim_YMass->cd(j+1);
      hsim_YMass = (TH1F*) fsim->Get(Form("h_YMass_%scut_%d",cut.Data(),j));
      hsim_YMass->Rebin(5);
      double xmin = hsim_YMass->GetXaxis()->GetXmin();
      double xmax = hsim_YMass->GetXaxis()->GetXmax();
      hsim_YMass->SetMarkerStyle(20);
      hsim_YMass->SetMarkerSize(0.7);
      cout<<"hsim_YMass = "<<hsim_YMass<<endl;
      hsim_YMass->Draw("e");

      // double par_y[4];
      // TF1 *fsig_y = new TF1("fsig_y","[0]*TMath::Voigt(x-[1],[2],[3])"); //[0]*TMath::BreitWigner(x,[1],[2])
      // fsig_y->SetLineStyle(9);
      // // fsig_y->SetLineColor(kRed);
      // fsig_y->SetParameters(15000,2.175,0.04,0.04);
      // hsim_YMass->Fit(fsig_y,"","mR",xmin,xmax);
      // fsig_y->GetParameters(&par_y[0]);
      // TF1 *fbkg_y = (TF1*) gROOT->GetFunction("chebyshev4");
          
      // TF1 *fbkg_phiy = new TF1("fbkg_phiy","pol4");
      // fbkg_phiy->SetLineStyle(9);
      // fbkg_phiy->SetLineColor(4);
      // fbkg_phiy->SetParameters(1.,1.,1.,1.,1.);
      // grphiy->Fit(fbkg_phiy,"","mR",1.46,3.20);
      // fbkg_phiy->GetParameters(&paramphiy[3]);
      // TF1 *fmodel_phiy = new TF1("fmodel_phiy","[0]*TMath::BreitWigner(x,[1],[2])+pol4(3)");
      // fmodel_phiy->SetParameters(&paramphiy[0]);
      // fmodel_phiy->SetParLimits(1,1.60,1.80);
      // // fmodel_phifo->SetParLimits(4,);
      // grphiy->Fit(fmodel_phiy,"","mR",1.46,3.20);
      
      // //w.factory("BreitWigner::sig1_YMass_mc(m_YMass_mc[1.9,2.5],mean1_YMass_mc[2.10,2.20],width1_YMass_mc[0.01,0.1])");
      // w.factory("Voigtian::sig1_YMass_mc(m_YMass_mc[1.9,2.5],mean1_YMass_mc[2.10,2.20],width1_YMass_mc[0.004,0.8],sigma1_YMass_mc[0.004,0.8])");//mean1_YMass_mc[2.10,2.20],width1_YMass_mc[0.004,0.8],sigma1_YMass_mc[0.004,0.8]
      // w.factory("Chebychev::bkg_YMass_mc(m_YMass_mc,{c0_YMass_mc[-10,10],c1_YMass_mc[-10,10],c2_YMass_mc[-10,10]})");
      // //w.factory("SUM:::model(nsig[0,100000000]*sig, nbkg[0,100000000]*bkg)"); //nsig[0,100000000]*sig2,
      // w.factory("SUM:::model_YMass_mc(nsig_YMass_mc[0,100000000]*sig1_YMass_mc, nbkg_YMass_mc[0,100000000]*bkg_YMass_mc)"); //, nsig_mc[0,100000000]*sig2_mc
      // w.var("m_YMass_mc")->SetTitle("m_YMass_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_YMass_mc("dh_YMass_mc","dh_YMass_mc",*w.var("m_YMass_mc"),Import(*hsim_YMass));
      // RooPlot* fr_YMass_mc = w.var("m_YMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
      // fr_YMass_mc->SetTitleOffset(0.90,"X");
      // fr_YMass_mc->SetTitleSize(0.06,"XYZ");
      // fr_YMass_mc->SetLabelSize(0.06,"xyz");
      // w.pdf("model_YMass_mc")->fitTo(dh_YMass_mc);
      // // result = w.pdf("model_YMass_mc")->fitTo(dh_YMass_mc,Extended(kTRUE),Save());
      // dh_YMass_mc.plotOn(fr_YMass_mc,RooFit::Name("ldh_YMass_mc"));
      // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig1_YMass_mc")),LineColor(kRed),RooFit::Name("lsig1_YMass_mc"));
      // //w.pdf("model_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig2_mc")),LineColor(kMagenta),RooFit::Name("lsig2_mc"));
      // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("bkg_YMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_YMass_mc"));
      // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc,RooFit::Name("lmodel_YMass_mc"));
      // w.pdf("model_YMass_mc")->paramOn(fr_YMass_mc,Layout(0.55,0.98,0.92));//,Parameters(RooArgSet(*w.var("nsig_mc"),*w.var("nbkg_mc"),*w.var("mean1_mc"),*w.var("mean2_mc"),*w.var("sigma1_mc"),*w.var("sigma2_mc"))));
      // fr_YMass_mc->Draw();

      // TLegend *l_YMass_mc = new TLegend(0.13,0.6,0.3,0.85);
      // l_YMass_mc->SetFillColor(kWhite);
      // l_YMass_mc->SetLineColor(kWhite);
      // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("ldh_YMass_mc"),"Data", "p");
      // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lmodel_YMass_mc"),"total","l");
      // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lbkg_YMass_mc"),"Background", "l");
      // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lsig1_YMass_mc"),"Signal1", "l");
      // // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lsig2_YMass_mc"),"Signal2", "l");
      // l_YMass_mc->Draw();
      
      // double fscale = 15.31; // fscale = (ndsig/ndbkg)/(nmcsig/nmcbkg) 
      // double nmcsig = w.var("nsig_YMass_mc")->getVal();
      // double nmcsigerr = w.var("nsig_YMass_mc")->getError();
      // double nmcbkg = w.var("nbkg_YMass_mc")->getVal();
      // double nmcbkgerr = w.var("nbkg_YMass_mc")->getError();
      // double zmc = nmcsig/sqrt(nmcsig+(fscale*nmcbkg));
      // double zmcerr = sqrt(((nmcsig+(2*fscale*nmcbkg))*(nmcsig+(2*fscale*nmcbkg))*(nmcsigerr*nmcsigerr) + (fscale*fscale)*(nmcsig*nmcsig)+(nmcbkgerr*nmcbkgerr))/(4*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))));
      // cout<<"nmcsig = "<<nmcsig<<"| nmcbkg = "<<nmcbkg<<"| zmc = "<<zmc<<endl;
      // // hsim_YMass->Draw("e");
      
      // // double zObs = NumberCountingUtils::BinomialObsZ(nmcsig+nmcbkg, nmcbkg, nmcbkgerr);
      // // cout << " Z value (Gaussian sigma) = "<< zObs << endl;
      
      // grsim_YMass->SetPoint(j,cuts[j],zmc);
      // grsim_YMass->SetPointError(j,0,zmcerr);
      // hsim_YMass->Write();	       
  }

  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%scut.root",cut.Data()), "root");
  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%scut.eps",cut.Data()), "eps");

  // cgrsim_YMass->cd();
  // grsim_YMass->Draw("AP");

  // // TF1 *fsimpol6 = new TF1("fsimpol6","pol6");
  // // grsim_YMass->Fit("fsimpol6","mR","", cuts[12], cuts[4]);
  // // double zmc_max = fsimpol6->GetMaximumX(cuts[12], cuts[4]);
  // // cout<<"zmc_max = "<<zmc_max<<endl;

  // // result->Print();
  // grsim_YMass->Write();
  // grsim_YMass->SetTitle("Significance Y(2175) Vs selection (MC); cut; Significance");
  // cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%scut.root",cut.Data()),"root");
  // cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%scut.eps",cut.Data()),"eps");

}
