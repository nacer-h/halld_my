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
#include "iostream"
#include "vector"
#include "TF1.h"
#include "TFormula.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "iostream"
#include "iomanip"
#include "vector"
#include "TLorentzVector.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TPaveStats.h"
using namespace std;

#ifndef __CINT__
#include "RooGlobalFunc.h"
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
using namespace RooFit;

RooRealVar rofit(TH1 *hist, RooRealVar x, double meanmin, double meanmax, double widthmin, double widthmax, double sigmamin, double sigmamax, double c0min, double c0max, double c1min, double c1max,double c2min, double c2max,double c3min, double c3max)
{
  // -------- RooFit steps

  // Declare observable x
  double xmin = hist->GetXaxis()->GetXmin();
  double xmax = hist->GetXaxis()->GetXmax();
  // RooRealVar x("x","x",xmin,xmax,"x");
  x.setRange("range",xmin,xmax);

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist dhist("data","dataset with x",x,Import(*hist));

  // Fit a p.d.f to the data

  // Signal pdf
  RooRealVar mean("#mu","mean",meanmin,meanmax,"GeV/c^{2}");
  RooRealVar width("#Gamma","width",widthmin,widthmax,"GeV/c^{2}");
  RooRealVar sigma("#sigma","sigma",sigmamin,sigmamax,"GeV/c^{2}");
  // RooRealVar mean1_fo("#mu_{1}","mean1",0.900,1.010,"GeV/c^{2}");
  // RooRealVar width_fo("#Gamma","width",0.001,0.1,"GeV/c^{2}");
  // RooRealVar mean2_fo("#mu_{2}","mean2",0.800,1.010,"GeV/c^{2}");
  // RooRealVar sigma1_fo("#sigma_{1}","sigma1",0.02,0.08,"GeV/c^{2}");
  // RooRealVar sigma2_fo("#sigma_{2}","sigma2",0.009,0.05,"GeV/c^{2}");//0.01,0.08
  // RooRealVar alpha_fo("#alpha","alpha",0.5,2.5);
  // RooRealVar n_fo("n","n",0.5,2.5);
  // RooVoigtian sig("sig","Voigtian(x,mean,gamma,sigma)",x,mean,width,sigma);
  // RooBreitWigner sig1_fo("sig1","BreitWigner(m_fo,mean1_fo,width_fo)",m_fo,mean1_fo,width_fo);
  // RooVoigtian sig1_fo("sig1","Voigtian(m_fo,mean1_fo,width_fo,sigma_fo)",m_fo,mean1_fo,width_fo,sigma_fo);
  // RooCBShape sig1_fo("sig1","CBShape(m_fo,mean1_fo,sigma_fo,alpha_fo,n_fo)", m_fo, mean1_fo, sigma_fo, alpha_fo, n_fo);
  // RooLandau sig2_fo("sig2","Landau(m_fo,mean2_fo,sigma2_fo)",m_fo,mean2_fo,sigma2_fo);
  RooBreitWigner sig("sig","BreitWigner(x,mean,width)",x,mean,width);
  RooRealVar nsig("N_{SIG}","number of signal events",0,100000000);
  // RooExtendPdf esig("esig","esig",sig,nsig,"range");

  // Background pdf

  RooRealVar c0("c0","coefficient #0",c0min,c0max);
  RooRealVar c1("c1","coefficient #1",c1min,c1max);
  RooRealVar c2("c2","coefficient #2",c2min,c2max);
  RooRealVar c3("c3","coefficient #3",c3min,c3max);
  // RooPolynomial bkg("pol","pol(c0)",x,RooArgList(c0));
  // RooPolynomial bkg("pol","pol(c0,c1,c2,c3)",x,RooArgList(c0,c1,c2,c3));
  RooChebychev bkg("bkg","Chebychev(c0)",x,RooArgList(c0,c1,c2));
  RooRealVar nbkg("N_{BKG}","number of background events",0,100000000);
  // RooExtendPdf ebkg("ebkg","ebkg",bkg,nbkg,"range");

  // RooRealVar sigfrac("sigfrac","signal fraction",0.5,0.,1.);
  RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
  // RooAddPdf model("model","model",RooArgList(esig,ebkg));
  model.fitTo(dhist);
  RooFitResult* result = model.fitTo(dhist,Extended(kTRUE),Save());
  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frhist = x.frame(Title("invariant Mass"));
  frhist->SetTitleOffset(0.90,"X");
  frhist->SetTitleSize(0.05,"XYZ");
  frhist->SetLabelSize(0.05,"xyz");
  dhist.plotOn(frhist,RooFit::Name("ldhist"));
  model.plotOn(frhist, Components(sig),LineColor(kRed),RooFit::Name("lsig"));
  model.plotOn(frhist, Components(bkg),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg"));
  model.plotOn(frhist,RooFit::Name("lmodel"));
  model.paramOn(frhist,Layout(0.55,0.98,0.92));//,Parameters(RooArgSet(nsig,nbkg,mean,width))
  frhist->Draw();
  TLegend *lhist = new TLegend(0.13,0.6,0.3,0.85);
  lhist->SetFillColor(kWhite);
  lhist->SetLineColor(kWhite);
  lhist->AddEntry(frhist->findObject("ldhist"),"Data", "p");
  lhist->AddEntry(frhist->findObject("lmodel"),"total","l");
  lhist->AddEntry(frhist->findObject("lbkg"),"Background", "l");
  lhist->AddEntry(frhist->findObject("lsig"),"Signal", "l");
  lhist->Draw();

  // RooArgList list(mean, width, "list" );

  return nsig;

  // RooRealVar m_phi("mass","m_{K^{+}K^{-}}",h1_PhiMass_kmkpcut->GetXaxis()->GetXmin(),h1_PhiMass_kmkpcut->GetXaxis()->GetXmax(),"GeV/c^{2}");
  // m_phi.setRange("Range",1.01, 1.03);
  // // rofit(h1_PhiMass_kmkpcut, m_phi, 1.018,1.021, 0.0004,0.09, 0.0001,0.009, -2.,2., -2.,2., -1.,1.);
  //
  // RooAbsArg* arg = rofit(h1_PhiMass_kmkpcut, m_phi, 1.018,1.021, 0.0004,0.09, 0.0001,0.009, -2.,2., -2.,2., -1.,1.).find( "mean" );
  // RooRealVar* mean = (RooRealVar*) arg->getPrameters("mean")

  // RooAbsArg* arg = list.find( "mean" );
  // RooAbsArg* arg = list.at( 0 );
  // RooAbsArg& arg = list[ 0 ];
}

void pippimkpkm(TString sim, TString data)
{
  TFile *fsim = new TFile(Form("/data.local/nacer/halld_my/pimpipkmkp/pippimkpkm_%s.root",sim.Data()));
  TFile *fdata = new TFile(Form("/data.local/nacer/halld_my/pimpipkmkp/pippimkpkm_%s.root",data.Data()));
  TFile *outputfig = new TFile(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/y.root",data.Data()),"UPDATE");
  // RooWorkspace *w = new RooWorkspace("w",kTRUE);
  RooWorkspace w("w",kTRUE);


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // *********************************************** Phi(1020)  ***********************************************
  TH1F *h1_PhiMass_kmkpcut = (TH1F*) fsim->Get("PhiMass_kmkpcut");
  h1_PhiMass_kmkpcut->Rebin(6);
  // h1_PhiMass_kmkpcut->SetStats(false);
  h1_PhiMass_kmkpcut->SetLineStyle(1);
  h1_PhiMass_kmkpcut->SetLineWidth(2);
  double xmin_phi = h1_PhiMass_kmkpcut->GetXaxis()->GetXmin();
  double xmax_phi = h1_PhiMass_kmkpcut->GetXaxis()->GetXmax();
  TCanvas *c_PhiMass_kmkpcut = new TCanvas("c_PhiMass_kmkpcut","c_PhiMass_kmkpcut",700,400);
  c_PhiMass_kmkpcut->cd();
  h1_PhiMass_kmkpcut->SetMarkerStyle(20);
  h1_PhiMass_kmkpcut->SetMarkerSize(0.7);
  h1_PhiMass_kmkpcut->Draw("e");

  //+++++++++++++ RooFit  +++++++++++

  // w.factory("BreitWigner::sig1_phi(m_phi[0.99,1.06],mean1_phi[1.016,1.022],width_phi[0.001,0.01])");
  // w.factory("Landau::sig2_phi(m_phi,mean2_phi[1.016,1.022],sigma_phi[0.0008,0.01])");

  w.factory("Voigtian::sig1_phi(m_phi[0.99,1.06],mean1_phi[1.016,1.022],width_phi[0.001,0.02],sigma_phi[0.0008,0.02])");
  // w.var("mean2_phi")->setVal(1.01946);
  // w.var("sigma_phi")->setVal(0.00424);
  w.factory("Chebychev::bkg_phi(m_phi,{c0_phi[-1.,1]})");
  w.factory("SUM:::model_phi(nsig_phi[0,100000000]*sig1_phi,  nbkg_phi[0,100000000]*bkg_phi)"); //nsig_phi[0,100000000]*sig2_phi,
  RooDataHist dh1_PhiMass_kmkpcut("dh1_PhiMass_kmkpcut","dh1_PhiMass_kmkpcut",*w.var("m_phi"),Import(*h1_PhiMass_kmkpcut));
  RooPlot* frh1_PhiMass_kmkpcut = w.var("m_phi")->frame(Title("invariant mass K^{+}K^{-}"));
  frh1_PhiMass_kmkpcut->SetTitleOffset(0.90,"X");
  frh1_PhiMass_kmkpcut->SetTitleSize(0.05,"XYZ");
  frh1_PhiMass_kmkpcut->SetLabelSize(0.05,"xyz");
  w.pdf("model_phi")->fitTo(dh1_PhiMass_kmkpcut);
  RooFitResult* result_phi = w.pdf("model_phi")->fitTo(dh1_PhiMass_kmkpcut,Extended(kTRUE),Save());
  dh1_PhiMass_kmkpcut.plotOn(frh1_PhiMass_kmkpcut,RooFit::Name("ldh1_PhiMass_kmkpcut"));
  w.pdf("model_phi")->plotOn(frh1_PhiMass_kmkpcut, Components(*w.pdf("sig1_phi")),LineColor(kRed),RooFit::Name("lsig1_phi"));
  // w.pdf("model_phi")->plotOn(frh1_PhiMass_kmkpcut, Components(*w.pdf("sig2_phi")),LineColor(kMagenta),RooFit::Name("lsig2_phi"));
  w.pdf("model_phi")->plotOn(frh1_PhiMass_kmkpcut, Components(*w.pdf("bkg_phi")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_phi"));
  w.pdf("model_phi")->plotOn(frh1_PhiMass_kmkpcut,RooFit::Name("lmodel_phi"));
  w.pdf("model_phi")->paramOn(frh1_PhiMass_kmkpcut,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet(nsig,nbkg,mean1,sigma)))
  frh1_PhiMass_kmkpcut->Draw();
  TLegend *lh1_PhiMass_kmkpcut = new TLegend(0.13,0.6,0.3,0.85);
  lh1_PhiMass_kmkpcut->SetFillColor(kWhite);
  lh1_PhiMass_kmkpcut->SetLineColor(kWhite);
  lh1_PhiMass_kmkpcut->AddEntry(frh1_PhiMass_kmkpcut->findObject("ldh1_PhiMass_kmkpcut"),"Data", "p");
  lh1_PhiMass_kmkpcut->AddEntry(frh1_PhiMass_kmkpcut->findObject("lmodel_phi"),"total","l");
  lh1_PhiMass_kmkpcut->AddEntry(frh1_PhiMass_kmkpcut->findObject("lbkg_phi"),"Background", "l");
  lh1_PhiMass_kmkpcut->AddEntry(frh1_PhiMass_kmkpcut->findObject("lsig1_phi"),"BreitWigner", "l");
  // lh1_PhiMass_kmkpcut->AddEntry(frh1_PhiMass_kmkpcut->findObject("lsig2_phi"),"Landau", "l");
  lh1_PhiMass_kmkpcut->Draw();
  result_phi->Print();
  w.Print();
  cout<<" %%%%%%%%%% nsig_phi = "<<w.var("nsig_phi")->getVal()<<endl;
  h1_PhiMass_kmkpcut->Write();
  c_PhiMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_PhiMass_kmkpcut.root", sim.Data(), sim.Data()),"root");
  c_PhiMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_PhiMass_kmkpcut.eps", sim.Data(), sim.Data()),"eps");

  //*************************************** fo(980) ***********************************************
  TH1F *h1_foMass_kmkpcut = (TH1F*) fsim->Get("foMass_kmkpcut");
  h1_foMass_kmkpcut->Rebin(6);
  // h1_foMass_kmkpcut->SetStats(false);
  h1_foMass_kmkpcut->SetLineStyle(1);
  h1_foMass_kmkpcut->SetLineWidth(2);
  TCanvas *c_foMass_kmkpcut = new TCanvas("c_foMass_kmkpcut","c_foMass_kmkpcut",700,400);
  c_foMass_kmkpcut->cd();
  h1_foMass_kmkpcut->SetMarkerStyle(20);
  h1_foMass_kmkpcut->SetMarkerSize(0.7);
  h1_foMass_kmkpcut->Draw("e");

  // w.factory("Gaussian::sig1_fo(m_fo[0.28,1.50],mean1_fo[0.900,1.010],sigma1_fo[0.02,0.08])");
  // w.factory("Gaussian::sig2_fo(m_fo,mean2_fo[0.800,1.010],sigma2_fo[0.009,0.05])");
  // w.factory("Chebychev::bkg_fo(m_fo,{c0_fo[0.,1]})");
  // w.factory("SUM:::model_fo(nsig_fo[0,100000000]*sig1_fo, nsig_fo[0,100000000]*sig2_fo, nbkg_fo[0,100000000]*bkg_fo)");
  w.factory("Voigtian::sig1_fo(m_fo[0.28,1.50],mean1_fo[0.900,1.010],width_fo[0.001,0.02],sigma_fo[0.02,0.08])");
  w.factory("Chebychev::bkg_fo(m_fo,{c0_fo[-1.,1]})");
  w.factory("SUM:::model_fo(nsig_fo[0,100000000]*sig1_fo, nbkg_fo[0,100000000]*bkg_fo)");
  RooDataHist dh1_foMass_kmkpcut("dh1_foMass_kmkpcut","dh1_foMass_kmkpcut",*w.var("m_fo"),Import(*h1_foMass_kmkpcut));
  RooPlot* frh1_foMass_kmkpcut = w.var("m_fo")->frame(Title("invariant mass #pi^{+}#pi^{-}"));
  frh1_foMass_kmkpcut->SetTitleOffset(0.90,"X");
  frh1_foMass_kmkpcut->SetTitleSize(0.05,"XYZ");
  frh1_foMass_kmkpcut->SetLabelSize(0.05,"xyz");
  w.pdf("model_fo")->fitTo(dh1_foMass_kmkpcut);
  RooFitResult* result_fo = w.pdf("model_fo")->fitTo(dh1_foMass_kmkpcut,Extended(kTRUE),Save());
  dh1_foMass_kmkpcut.plotOn(frh1_foMass_kmkpcut,RooFit::Name("ldh1_foMass_kmkpcut"));
  w.pdf("model_fo")->plotOn(frh1_foMass_kmkpcut, Components(*w.pdf("sig1_fo")),LineColor(kRed),RooFit::Name("lsig1_fo"));
  // w.pdf("model_fo")->plotOn(frh1_foMass_kmkpcut, Components(*w.pdf("sig2_fo")),LineColor(kMagenta),RooFit::Name("lsig2_fo"));
  w.pdf("model_fo")->plotOn(frh1_foMass_kmkpcut, Components(*w.pdf("bkg_fo")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_fo"));
  w.pdf("model_fo")->plotOn(frh1_foMass_kmkpcut,RooFit::Name("lmodel_fo"));
  w.pdf("model_fo")->paramOn(frh1_foMass_kmkpcut,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet((*w->var("nsig")),(*w->var("nbkg")),(*w->var("mean2")),(*w->var("sigma2"))))
  frh1_foMass_kmkpcut->Draw();
  TLegend *lh1_foMass_kmkpcut = new TLegend(0.13,0.6,0.3,0.85);
  lh1_foMass_kmkpcut->SetFillColor(kWhite);
  lh1_foMass_kmkpcut->SetLineColor(kWhite);
  lh1_foMass_kmkpcut->AddEntry(frh1_foMass_kmkpcut->findObject("ldh1_foMass_kmkpcut"),"Data", "p");
  lh1_foMass_kmkpcut->AddEntry(frh1_foMass_kmkpcut->findObject("lmodel_fo"),"total","l");
  lh1_foMass_kmkpcut->AddEntry(frh1_foMass_kmkpcut->findObject("lbkg_fo"),"Background", "l");
  lh1_foMass_kmkpcut->AddEntry(frh1_foMass_kmkpcut->findObject("lsig1_fo"),"Gaussian", "l");
  // lh1_foMass_kmkpcut->AddEntry(frh1_foMass_kmkpcut->findObject("lsig2_fo"),"Gaussian", "l");
  lh1_foMass_kmkpcut->Draw();
  // result_fo->Print();
  h1_foMass_kmkpcut->Write();
  c_foMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_foMass_kmkpcut.root", sim.Data(), sim.Data()),"root");
  c_foMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_foMass_kmkpcut.eps", sim.Data(), sim.Data()),"eps");

  // *********************************************** Y(2175)  ***********************************************
  TH1F *h1_YMass_kmkpcut = (TH1F*) fsim->Get("YMass_kmkpcut");
  h1_YMass_kmkpcut->Rebin(6);
  // h1_YMass_kmkpcut->SetStats(false);
  h1_YMass_kmkpcut->SetLineStyle(1);
  h1_YMass_kmkpcut->SetLineWidth(2);

  TCanvas *c_YMass_kmkpcut = new TCanvas("c_YMass_kmkpcut","c_YMass_kmkpcut",700,400);
  c_YMass_kmkpcut->cd();
  h1_YMass_kmkpcut->SetMarkerStyle(20);
  h1_YMass_kmkpcut->SetMarkerSize(0.7);
  h1_YMass_kmkpcut->Draw("e");

  // w.factory("Gaussian::sig1_y(m_y[1.50,3.20],mean1_y[2.100,2.300],sigma1_y[0.01,0.1])");
  // w.factory("Gaussian::sig2_y(m_y,mean2_y[2.100,2.300],sigma2_y[0.01,0.1])");
  // w.factory("Chebychev::bkg_y(m_y,{c0_y[0.,1.]})");
  // w.factory("SUM:::model_y(nsig_y[0,100000000]*sig1_y, nsig_y[0,100000000]*sig2_y, nbkg_y[0,100000000]*bkg_y)");
  w.factory("Voigtian::sig1_y(m_y[1.50,3.20],mean1_y[2.100,2.300],width_y[0.001,0.02],sigma_y[0.01,0.1])");
  w.factory("Chebychev::bkg_y(m_y,{c0_y[-1.,1])}");
  w.factory("SUM:::model_y(nsig_y[0,100000000]*sig1_y, nbkg_y[0,100000000]*bkg_y)");
  RooDataHist dh1_YMass_kmkpcut("dh1_YMass_kmkpcut","dh1_YMass_kmkpcut",*w.var("m_y"),Import(*h1_YMass_kmkpcut));
  RooPlot* frh1_YMass_kmkpcut = w.var("m_y")->frame(Title("invariant mass #pi^{+}#pi^{-}K^{+}K^{-}"));
  frh1_YMass_kmkpcut->SetTitleOffset(0.90,"X");
  frh1_YMass_kmkpcut->SetTitleSize(0.05,"XYZ");
  frh1_YMass_kmkpcut->SetLabelSize(0.05,"xyz");
  w.pdf("model_y")->fitTo(dh1_YMass_kmkpcut);
  RooFitResult* result_y = w.pdf("model_y")->fitTo(dh1_YMass_kmkpcut,Extended(kTRUE),Save());
  dh1_YMass_kmkpcut.plotOn(frh1_YMass_kmkpcut,RooFit::Name("ldh1_YMass_kmkpcut"));
  w.pdf("model_y")->plotOn(frh1_YMass_kmkpcut, Components(*w.pdf("sig1_y")),LineColor(kRed),RooFit::Name("lsig1_y"));
  // w.pdf("model_y")->plotOn(frh1_YMass_kmkpcut, Components(*w.pdf("sig2_y")),LineColor(kMagenta),RooFit::Name("lsig2_y"));
  w.pdf("model_y")->plotOn(frh1_YMass_kmkpcut, Components(*w.pdf("bkg_y")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_y"));
  w.pdf("model_y")->plotOn(frh1_YMass_kmkpcut,RooFit::Name("lmodel_y"));
  w.pdf("model_y")->paramOn(frh1_YMass_kmkpcut,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet((*w->var("nsig")),(*w->var("nbkg")),(*w->var("mean2")),(*w->var("sigma2"))))
  frh1_YMass_kmkpcut->Draw();
  TLegend *lh1_YMass_kmkpcut = new TLegend(0.13,0.6,0.3,0.85);
  lh1_YMass_kmkpcut->SetFillColor(kWhite);
  lh1_YMass_kmkpcut->SetLineColor(kWhite);
  lh1_YMass_kmkpcut->AddEntry(frh1_YMass_kmkpcut->findObject("ldh1_YMass_kmkpcut"),"Data", "p");
  lh1_YMass_kmkpcut->AddEntry(frh1_YMass_kmkpcut->findObject("lmodel_y"),"total","l");
  lh1_YMass_kmkpcut->AddEntry(frh1_YMass_kmkpcut->findObject("lbkg_y"),"Background", "l");
  lh1_YMass_kmkpcut->AddEntry(frh1_YMass_kmkpcut->findObject("lsig1_y"),"Gaussian", "l");
  // lh1_YMass_kmkpcut->AddEntry(frh1_YMass_kmkpcut->findObject("lsig2_y"),"Gaussian", "l");
  lh1_YMass_kmkpcut->Draw();
  result_y->Print();
  h1_YMass_kmkpcut->Write();
  c_YMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_YMass_kmkpcut.root", sim.Data(), sim.Data()),"root");
  c_YMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_YMass_kmkpcut.eps", sim.Data(), sim.Data()),"eps");
return;
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  TH1F *hdata_PhiMass_kmkpcut = (TH1F*) fdata->Get("PhiMass_kmkpcut");
  hdata_PhiMass_kmkpcut->Rebin(6);
  // hdata_PhiMass_kmkpcut->SetStats(false);
  hdata_PhiMass_kmkpcut->SetLineStyle(1);
  hdata_PhiMass_kmkpcut->SetLineWidth(2);

  TCanvas *c_data_PhiMass_kmkpcut = new TCanvas("c_data_PhiMass_kmkpcut","c_data_PhiMass_kmkpcut",700,400);
  c_data_PhiMass_kmkpcut->cd();
  hdata_PhiMass_kmkpcut->SetMarkerStyle(20);
  hdata_PhiMass_kmkpcut->SetMarkerSize(0.7);
  // hdata_PhiMass_kmkpcut->Draw("e");


  // w.factory("BreitWigner::sig1_phi(m_phi[0.99,1.06],mean1_phi[1.016,1.022],width_phi[0.001,0.01])");
  // w.factory("Landau::sig2_phi(m_phi,mean2_phi[1.016,1.022],sigma_phi[0.0008,0.01])");

  w.factory("Voigtian::sig1_phi_data(m_phi_data[0.99,1.06],mean1_phi_data[1.016,1.022],width_phi_data[0.001,0.02],sigma_phi_data[0.0008,0.02])");
  w.factory("Chebychev::bkg_phi_data(m_phi_data,{c0_phi_data[-1.,1],c1_phi_data[-1.,1],c2_phi_data[-1.,1]})");
  w.factory("SUM:::model_phi_data(nsig_phi_data[0,100000000]*sig1_phi_data, nbkg_phi_data[0,100000000]*bkg_phi_data)"); //nsig_phi[0,100000000]*sig2_phi,
  RooDataHist dhdata_PhiMass_kmkpcut("dhdata_PhiMass_kmkpcut","dhdata_PhiMass_kmkpcut",*w.var("m_phi_data"),Import(*hdata_PhiMass_kmkpcut));
  RooPlot* frhdata_PhiMass_kmkpcut = w.var("m_phi_data")->frame(Title("invariant mass K^{+}K^{-}"));
  frhdata_PhiMass_kmkpcut->SetTitleOffset(0.90,"X");
  frhdata_PhiMass_kmkpcut->SetTitleSize(0.05,"XYZ");
  frhdata_PhiMass_kmkpcut->SetLabelSize(0.05,"xyz");
  w.pdf("model_phi_data")->fitTo(dhdata_PhiMass_kmkpcut);
  RooFitResult* result_phi_data = w.pdf("model_phi_data")->fitTo(dhdata_PhiMass_kmkpcut,Extended(kTRUE),Save());
  dhdata_PhiMass_kmkpcut.plotOn(frhdata_PhiMass_kmkpcut,RooFit::Name("ldh1_PhiMass_kmkpcut_data"));
  w.pdf("model_phi_data")->plotOn(frhdata_PhiMass_kmkpcut, Components(*w.pdf("sig1_phi_data")),LineColor(kRed),RooFit::Name("lsig1_phi_data"));
  // w.pdf("model_phi_data")->plotOn(frhdata_PhiMass_kmkpcut, Components(*w.pdf("sig2_phi_data")),LineColor(kMagenta),RooFit::Name("lsig2_phi_data"));
  w.pdf("model_phi_data")->plotOn(frhdata_PhiMass_kmkpcut, Components(*w.pdf("bkg_phi_data")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_phi_data"));
  w.pdf("model_phi_data")->plotOn(frhdata_PhiMass_kmkpcut,RooFit::Name("lmodel_phi_data"));
  w.pdf("model_phi_data")->paramOn(frhdata_PhiMass_kmkpcut,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet(nsig,nbkg,mean1,sigma)))
  frhdata_PhiMass_kmkpcut->Draw();
  TLegend *lhdata_PhiMass_kmkpcut = new TLegend(0.13,0.6,0.3,0.85);
  lhdata_PhiMass_kmkpcut->SetFillColor(kWhite);
  lhdata_PhiMass_kmkpcut->SetLineColor(kWhite);
  lhdata_PhiMass_kmkpcut->AddEntry(frhdata_PhiMass_kmkpcut->findObject("ldh1_PhiMass_kmkpcut_data"),"Data", "p");
  lhdata_PhiMass_kmkpcut->AddEntry(frhdata_PhiMass_kmkpcut->findObject("lmodel_phi_data"),"total","l");
  lhdata_PhiMass_kmkpcut->AddEntry(frhdata_PhiMass_kmkpcut->findObject("lbkg_phi_data"),"Background", "l");
  lhdata_PhiMass_kmkpcut->AddEntry(frhdata_PhiMass_kmkpcut->findObject("lsig1_phi_data"),"BreitWigner", "l");
  // lhdata_PhiMass_kmkpcut->AddEntry(frhdata_PhiMass_kmkpcut->findObject("lsig2_phi_data"),"Landau", "l");
  lhdata_PhiMass_kmkpcut->Draw();
  result_phi_data->Print();

  // RooRealVar m_phi_data("mass","m_{K^{+}K^{-}}",hdata_PhiMass_kmkpcut->GetXaxis()->GetXmin(),hdata_PhiMass_kmkpcut->GetXaxis()->GetXmax(),"GeV/c^{2}");
  // RooRealVar N_phi = rofit(hdata_PhiMass_kmkpcut, m_phi_data, 1.018,1.020, 0.009,0.05, 0.003,0.1, -1.,1., -1.,1.,-1.,1.,-1,1);

  hdata_PhiMass_kmkpcut->Write();
  c_data_PhiMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_PhiMass_kmkpcut.root", data.Data(), data.Data()),"root");
  c_data_PhiMass_kmkpcut->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_PhiMass_kmkpcut.eps", data.Data(), data.Data()),"eps");


  // ---- look for fo(980) in bins of K+K- -----

  TCanvas *ckmkpvspimpip = new TCanvas("ckmkpvspimpip","ckmkpvspimpip",700,400);
  ckmkpvspimpip->cd();
  TH2F *hkmkpvspimpip = (TH2F*) fdata->Get("kmkpvspimpip_kmkpcut");
  hkmkpvspimpip->RebinX(6);
  hkmkpvspimpip->RebinY(6);
  hkmkpvspimpip->Draw("colz");
  hkmkpvspimpip->Write();
  ckmkpvspimpip->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_kmkpvspimpip.root", data.Data(), data.Data()),"root");
  ckmkpvspimpip->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_kmkpvspimpip.eps", data.Data(), data.Data()),"eps");

  TCanvas *cphifo1 = new TCanvas("cphifo1","cphifo1",1500,1500);
  cphifo1->Divide(5,5,0.00001,0.00001);
  TCanvas *cphifo2 = new TCanvas("cphifo2","cphifo2",1500,1500);
  cphifo2->Divide(5,5,0.00001,0.00001);

  TH1D *hprojphifo;
  hkmkpvspimpip->RebinY(2);
  hkmkpvspimpip->RebinX(2);

  TGraphErrors *grphifo = new TGraphErrors();
  grphifo->SetTitle("Yield of #phi(1020) Vs m_{#pi^{+}#pi^{-}}; m_{#pi^{+}#pi^{-}} [GeV/c^{2}]; N_{K^{+}K^{-}}");
  grphifo->SetMarkerStyle(20);
  grphifo->GetXaxis()->SetTitleOffset(0.90);
  grphifo->GetXaxis()->SetTitleSize(0.05);
  grphifo->GetYaxis()->SetTitleSize(0.05);
  grphifo->GetXaxis()->SetLabelSize(0.05);
  grphifo->GetYaxis()->SetLabelSize(0.05);

  TH1D *hbwmeanphifo = new TH1D("hbwmeanphifo","Mean of the Breit-Wigner fit;m [GeV/c^{2}];Counts",50,1.018,1.021);
  TH1D *hbwwidthphifo = new TH1D("hbwwidthphifo","Width of the Breit-Wigner fit;#Gamma [GeV/c^{2}];Counts",50,0.008,0.020);

  double nbinsphifo = hkmkpvspimpip->GetYaxis()->GetNbins();

  for(int i=0; i<nbinsphifo; i++)
  {
    if(i<25) cphifo1->cd(i+1);
    else cphifo2->cd(i-25+1);
    hprojphifo = hkmkpvspimpip->ProjectionX(Form("binphifo%d",i+1),i+1,i+1);
    double xminphifo = hprojphifo->GetXaxis()->GetXmin();
    double xmaxphifo = hprojphifo->GetXaxis()->GetXmax();
    double binwphifo = hprojphifo->GetBinWidth(1);
    hprojphifo->SetMarkerStyle(20);
    hprojphifo->SetMarkerSize(0.7);

    w.factory("Voigtian::sig1_phifo_data(m_phifo_data[0.99,1.06],mean1_phifo_data[1.016,1.022],width_phifo_data[0.001,0.02],sigma_phifo_data[0.0008,0.02])");
    w.factory("Chebychev::bkg_phifo_data(m_phifo_data,{c0_phifo_data[-1.,1],c1_phifo_data[-1.,1],c2_phifo_data[-1.,1]})");
    w.factory("SUM:::model_phifo_data(nsig_phifo_data[0,100000000]*sig1_phifo_data, nbkg_phifo_data[0,100000000]*bkg_phifo_data)");
    RooDataHist dhprojphifo("dhprojphifo","dhprojphifo",*w.var("m_phifo_data"),Import(*hprojphifo));
    RooPlot* frhprojphifo = w.var("m_phifo_data")->frame(Title("invariant mass K^{+}K^{-}"));
    frhprojphifo->SetTitleOffset(0.90,"X");
    frhprojphifo->SetTitleSize(0.05,"XYZ");
    frhprojphifo->SetLabelSize(0.05,"xyz");
    w.pdf("model_phifo_data")->fitTo(dhprojphifo);
    RooFitResult* result_phifo_data = w.pdf("model_phifo_data")->fitTo(dhprojphifo,Extended(kTRUE),Save());
    dhprojphifo.plotOn(frhprojphifo,RooFit::Name("ldh1_PhiMass_kmkpcut_data"));
    w.pdf("model_phifo_data")->plotOn(frhprojphifo, Components(*w.pdf("sig1_phifo_data")),LineColor(kRed),RooFit::Name("lsig1_phifo_data"));
    // w.pdf("model_phifo_data")->plotOn(frhprojphifofo, Components(*w.pdf("sig2_phifo_data")),LineColor(kMagenta),RooFit::Name("lsig2_phifo_data"));
    w.pdf("model_phifo_data")->plotOn(frhprojphifo, Components(*w.pdf("bkg_phifo_data")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_phifo_data"));
    w.pdf("model_phifo_data")->plotOn(frhprojphifo,RooFit::Name("lmodel_phifo_data"));
    w.pdf("model_phifo_data")->paramOn(frhprojphifo,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet(nsig,nbkg,mean1,sigma)))
    frhprojphifo->Draw();
    TLegend *lhprojphifo = new TLegend(0.13,0.6,0.3,0.85);
    lhprojphifo->SetFillColor(kWhite);
    lhprojphifo->SetLineColor(kWhite);
    lhprojphifo->AddEntry(frhprojphifo->findObject("ldh1_phifoMass_kmkpcut_data"),"Data", "p");
    lhprojphifo->AddEntry(frhprojphifo->findObject("lmodel_phi_data"),"total","l");
    lhprojphifo->AddEntry(frhprojphifo->findObject("lbkg_phifo_data"),"Background", "l");
    lhprojphifo->AddEntry(frhprojphifo->findObject("lsig1_phifo_data"),"BreitWigner", "l");
    // lhprojphifofo->AddEntry(frhprojphifofo->findObject("lsig2_phifo_data"),"Landau", "l");
    lhprojphifo->Draw();
    result_phifo_data->Print();

    // w.var("m_phi_data")->setRange("Range",xminphifo,xmaxphifo);
    // // RooRealVar nsig_phifo = rofit(hprojphifo, m_phi_data, 1.018,1.020, 0.009,0.05, 0.003,0.1, -1.,1., -1.,1.,-1.,1.,-1.,1.);
    // RooRealVar nsig_phifo = rofit(hprojphifo, w.var("m_phi_data"), 1.018,1.020, 0.001,0.1, 0.003,0.1, -1.,1., -1.,1.,-1.,1.,-1.,1.);

    double nsigphifo = w.var("nsig_phifo_data")->getVal();
    // double nbkgphifo = nbkg.getVal(nbkg)/binwphifo;
    double pippim_m = hkmkpvspimpip->GetYaxis()->GetBinCenter(i);
    grphifo->SetPoint(i,pippim_m,nsigphifo);
    grphifo->SetPointError(i,0,w.var("nsig_phifo_data")->getError());
    // hbwmeanphifo->Fill(model.getVal(mean));
    // hbwwidthphifo->Fill(model.getVal(width));

    //   cout<<"---------- i = "<<i<<" | pippim_m = "<<pippim_m<<" | nsigphifo = "<<nsigphifo<<" | meanphifo = "<<mean.getVal(mean)<<" | widthphifo = "<<width.getVal(width)<<endl;
  }

  gStyle->SetHistMinimumZero();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);

  cphifo1->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo1.root", data.Data(), data.Data()),"root");
  cphifo1->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo1.eps", data.Data(), data.Data()),"eps");
  cphifo2->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo2.root", data.Data(), data.Data()),"root");
  cphifo2->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo2.eps", data.Data(), data.Data()),"eps");


  TCanvas *cphifo = new TCanvas("cphifo","cphifo",700,400);
  cphifo->cd();
  // double paramphifo[5];
  // TF1 *fsig1_phifo = new TF1("fsig1_phifo","[0]*TMath::BreitWigner(x,[1],[2])");
  // fsig1_phifo->SetLineStyle(9);
  // fsig1_phifo->SetLineColor(kMagenta);
  // fsig1_phifo->SetParameters(900000,0.453,0.267);
  // grphifo->Fit(fsig1_phifo,"","mR",0.25,0.50);
  // fsig1_phifo->GetParameters(&paramphifo[0]);
  // TF1 *fsig2_phifo = new TF1("fsig2_phifo","[0]*TMath::BreitWigner(x,[1],[2])");
  // fsig2_phifo->SetLineStyle(9);
  // fsig2_phifo->SetLineColor(4);
  // fsig2_phifo->SetParameters(9000000,0.642,0.425);
  // grphifo->Fit(fsig2_phifo,"","mR",0.50,1.50);
  // fsig2_phifo->GetParameters(&paramphifo[3]);
  // TF1 *fmodel_phifo = new TF1("fmodel_phifo","[0]*TMath::BreitWigner(x,[1],[2])+[3]*TMath::BreitWigner(x,[4],[5])");
  // fmodel_phifo->SetParameters(&paramphifo[0]);
  // fmodel_phifo->SetParLimits(1,0.25,0.50);
  // fmodel_phifo->SetParLimits(4,0.50,1.50);
  // grphifo->Fit(fmodel_phifo,"","mR",0.25,1.50);

  grphifo->Draw("AP");
  grphifo->Write();
  // fmodel_phifo->DrawCopy("same");
  // fmodel_phifo->Write();
  // fsig1_phifo->DrawCopy("same");
  // fsig1_phifo->Write();
  // fsig2_phifo->DrawCopy("same");
  // fsig2_phifo->Write();
  //
  // TLatex texphifo;
  // texphifo.SetTextAlign(12);
  // texphifo.SetTextSize(0.05);
  // texphifo.SetTextColor(2);
  // texphifo.SetNDC(kTRUE);
  // texphifo.DrawLatex(0.62,0.57,Form("#mu_{?} = %0.3f[GeV/c^{2}]",fmodel_phifo->GetParameter(1)));
  // texphifo.DrawLatex(0.62,0.49,Form("#Gamma_{?} = %0.4f[GeV/c^{2}]",fmodel_phifo->GetParameter(2)));
  // texphifo.DrawLatex(0.62,0.41,Form("#mu_{f_{0}(600)} = %0.3f[GeV/c^{2}]",fmodel_phifo->GetParameter(4)));
  // texphifo.DrawLatex(0.62,0.33,Form("#Gamma_{f_{0}(600)} = %0.4f[GeV/c^{2}]",fmodel_phifo->GetParameter(5)));

  cphifo->Modified();
  cphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo.root", data.Data(), data.Data()),"root");
  cphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phifo.eps", data.Data(), data.Data()),"eps");

  // TCanvas *cbwmeanphifo = new TCanvas("cbwmeanphifo","cbwmeanphifo",700,400);
  // cbwmeanphifo->cd();
  // hbwmeanphifo->SetTitleOffset(0.90,"X");
  // hbwmeanphifo->SetTitleSize(0.05,"XYZ");
  // hbwmeanphifo->SetLabelSize(0.05,"xyz");
  // hbwmeanphifo->Draw();
  // hbwmeanphifo->Write();
  // cbwmeanphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwmeanphifo.root", data.Data(), data.Data()),"root");
  // cbwmeanphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwmeanphifo.eps", data.Data(), data.Data()),"eps");
  //
  // TCanvas *cbwwidthphifo = new TCanvas("cbwwidthphifo","cbwwidthphifo",700,400);
  // cbwwidthphifo->cd();
  // hbwwidthphifo->SetTitleOffset(0.90,"X");
  // hbwwidthphifo->SetTitleSize(0.05,"XYZ");
  // hbwwidthphifo->SetLabelSize(0.05,"xyz");
  // hbwwidthphifo->Draw();
  // hbwwidthphifo->Write();
  // cbwwidthphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwwidthphifo.root", data.Data(), data.Data()),"root");
  // cbwwidthphifo->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwwidthphifo.eps", data.Data(), data.Data()),"eps");

  //----- look for Y(2175) in bins of K+K-  --------

  TCanvas *ckmkpvskmkppimpip = new TCanvas("ckmkpvskmkppimpip","ckmkpvskmkppimpip",700,400);
  ckmkpvskmkppimpip->cd();
  TH2F *hkmkpvskmkppimpip = (TH2F*) fdata->Get("kmkpvskmkppimpip_kmkpcut");
  hkmkpvskmkppimpip->RebinX(6);
  hkmkpvskmkppimpip->RebinY(6);
  hkmkpvskmkppimpip->Draw("colz");
  hkmkpvskmkppimpip->Write();
  ckmkpvskmkppimpip->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_kmkpvskmkppimpip.root", data.Data(), data.Data()),"root");
  ckmkpvskmkppimpip->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_kmkpvskmkppimpip.eps", data.Data(), data.Data()),"eps");

  TCanvas *cphiy1 = new TCanvas("cphiy1","cphiy1",1500,1500);
  cphiy1->Divide(5,5,0.00001,0.00001);
  TCanvas *cphiy2 = new TCanvas("cphiy2","cphiy2",1500,1500);
  cphiy2->Divide(5,5,0.00001,0.00001);

  TH1D *hprojphiy;
  hkmkpvskmkppimpip->RebinY(2);
  hkmkpvskmkppimpip->RebinX(2);

  TGraphErrors *grphiy = new TGraphErrors();
  grphiy->SetTitle("Yield of #phi(1020) Vs m_{#pi^{+}#pi^{-}K^{+}K^{-}}; m_{#pi^{+}#pi^{-}K^{+}K^{-}} [GeV/c^{2}]; N_{K^{+}K^{-}}");
  grphiy->SetMarkerStyle(20);
  grphiy->GetXaxis()->SetTitleOffset(0.90);
  grphiy->GetXaxis()->SetTitleSize(0.05);
  grphiy->GetYaxis()->SetTitleSize(0.05);
  grphiy->GetXaxis()->SetLabelSize(0.05);
  grphiy->GetYaxis()->SetLabelSize(0.05);
  grphiy->SetMarkerSize(1.5);

  TH1D *hbwmeanphiy = new TH1D("hbwmeanphiy","Mean of the Breit-Wigner fit;m [GeV/c^{2}];Counts",50,1.018,1.021);
  TH1D *hbwwidthphiy = new TH1D("hbwwidthphiy","Width of the Breit-Wigner fit;#Gamma [GeV/c^{2}];Counts",50,0.008,0.020);

  TH2F *hphiy  = new TH2F("hphiy","hphiy",50,1.50,3.20,50,0,3000000);

  double nbinsphiy = hkmkpvskmkppimpip->GetYaxis()->GetNbins();

  for(int i=0; i<nbinsphiy; i++)
  {
    if(i<25) cphiy1->cd(i+1);
    else cphiy2->cd(i-25+1);
    hprojphiy = hkmkpvskmkppimpip->ProjectionX(Form("binphiy%d",i+1),i+1,i+1);
    double xminphiy = hprojphiy->GetXaxis()->GetXmin();
    double xmaxphiy = hprojphiy->GetXaxis()->GetXmax();
    double binwphiy = hprojphiy->GetBinWidth(1);
    hprojphiy->SetMarkerStyle(20);
    hprojphiy->SetMarkerSize(0.7);

    w.factory("Voigtian::sig1_phiy_data(m_phiy_data[0.99,1.06],mean1_phiy_data[1.016,1.022],width_phiy_data[0.001,0.02],sigma_phiy_data[0.0008,0.02])");
    w.factory("Chebychev::bkg_phiy_data(m_phiy_data,{c0_phiy_data[-1.,1],c1_phiy_data[-1.,1],c2_phiy_data[-1.,1]})");
    w.factory("SUM:::model_phiy_data(nsig_phiy_data[0,100000000]*sig1_phiy_data, nbkg_phiy_data[0,100000000]*bkg_phiy_data)");
    RooDataHist dhprojphiy("dhprojphiy","dhprojphiy",*w.var("m_phiy_data"),Import(*hprojphiy));
    RooPlot* frhprojphiy = w.var("m_phiy_data")->frame(Title("invariant mass K^{+}K^{-}"));
    frhprojphiy->SetTitleOffset(0.90,"X");
    frhprojphiy->SetTitleSize(0.05,"XYZ");
    frhprojphiy->SetLabelSize(0.05,"xyz");
    w.pdf("model_phiy_data")->fitTo(dhprojphiy);
    RooFitResult* result_phiy_data = w.pdf("model_phiy_data")->fitTo(dhprojphiy,Extended(kTRUE),Save());
    dhprojphiy.plotOn(frhprojphiy,RooFit::Name("ldh1_PhiMass_kmkpcut_data"));
    w.pdf("model_phiy_data")->plotOn(frhprojphiy, Components(*w.pdf("sig1_phiy_data")),LineColor(kRed),RooFit::Name("lsig1_phiy_data"));
    // w.pdf("model_phiy_data")->plotOn(frhprojphiy, Components(*w.pdf("sig2_phiy_data")),LineColor(kMagenta),RooFit::Name("lsig2_phiy_data"));
    w.pdf("model_phiy_data")->plotOn(frhprojphiy, Components(*w.pdf("bkg_phiy_data")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_phiy_data"));
    w.pdf("model_phiy_data")->plotOn(frhprojphiy,RooFit::Name("lmodel_phiy_data"));
    w.pdf("model_phiy_data")->paramOn(frhprojphiy,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet(nsig,nbkg,mean1,sigma)))
    frhprojphiy->Draw();
    TLegend *lhprojphiy = new TLegend(0.13,0.6,0.3,0.85);
    lhprojphiy->SetFillColor(kWhite);
    lhprojphiy->SetLineColor(kWhite);
    lhprojphiy->AddEntry(frhprojphiy->findObject("ldh1_PhiMass_kmkpcut_data"),"Data", "p");
    lhprojphiy->AddEntry(frhprojphiy->findObject("lmodel_phiy_data"),"total","l");
    lhprojphiy->AddEntry(frhprojphiy->findObject("lbkg_phiy_data"),"Background", "l");
    lhprojphiy->AddEntry(frhprojphiy->findObject("lsig1_phiy_data"),"BreitWigner", "l");
    // lhprojphiy->AddEntry(frhprojphiy->findObject("lsig2_phi_data"),"Landau", "l");
    lhprojphiy->Draw();
    result_phiy_data->Print();

    // w.var("m_phi_data")->setRange("Range",xminphiy,xmaxphiy);
    // RooRealVar nsig_phiy = rofit(hprojphiy, w.var("m_phi_data"), 1.018,1.020, 0.009,0.05, 0.003,0.1, -1.,1., -1.,1.,-1.,1.,-1.,1.);

    double nsigphiy = w.var("nsig_phiy_data")->getVal();//binwphiy;
    // double nbkgphiy = model.getVal(nbkg)/binwphiy;
    double pippimkpkm_m = hkmkpvskmkppimpip->GetYaxis()->GetBinCenter(i);
    grphiy->SetPoint(i,pippimkpkm_m,nsigphiy);
    grphiy->SetPointError(i,0,w.var("nsig_phiy_data")->getError());///nsig_phiy.getVal()*nsigphiy);
    // hbwmeanphiy->Fill(mean.getVal(mean));
    // hbwwidthphiy->Fill(width.getVal(width));
    hphiy->Fill(pippimkpkm_m,nsigphiy);

    // cout<<"---------- i = "<<i<<" | pippimkpkm_m = "<<pippimkpkm_m<<" | nsigphiy = "<<nsigphiy<<" | meanphiy = "<<mean.getVal(mean)<<" | widthphiy = "<<width.getVal(width)<<endl;
  }


  cphiy1->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy1.root", data.Data(), data.Data()),"root");
  cphiy1->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy1.eps", data.Data(), data.Data()),"eps");
  cphiy2->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy2.root", data.Data(), data.Data()),"root");
  cphiy2->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy2.eps", data.Data(), data.Data()),"eps");

  TCanvas *cphiy = new TCanvas("cphiy","cphiy",700,400);
  cphiy->cd();
  // double paramphiy[7];
  // TF1 *fsig_phiy = new TF1("fsig_phiy","[0]*TMath::BreitWigner(x,[1],[2])");
  // fsig_phiy->SetLineStyle(9);
  // fsig_phiy->SetLineColor(kMagenta);
  // fsig_phiy->SetParameters(900000,1.73,0.2);
  // grphiy->Fit(fsig_phiy,"","mR",1.60,1.80);
  // fsig_phiy->GetParameters(&paramphiy[0]);
  // // TF1 *fbkg_phiy = (TF1*) gROOT->GetFunction("chebyshev4");
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

  grphiy->Draw("AP");
  grphiy->Write();
  // fmodel_phiy->DrawCopy("same");
  // fmodel_phiy->Write();
  // fsig_phiy->DrawCopy("same");
  // fsig_phiy->Write();
  // fbkg_phiy->DrawCopy("same");
  // fbkg_phiy->Write();
  //
  // TLatex texphiy;
  // texphiy.SetTextAlign(12);
  // texphiy.SetTextSize(0.05);
  // texphiy.SetTextColor(2);
  // texphiy.SetNDC(kTRUE);
  // texphiy.DrawLatex(0.30,0.85,Form("#mu_{Y(1750)} = %0.3f[GeV/c^{2}]",fmodel_phiy->GetParameter(1)));
  // texphiy.DrawLatex(0.30,0.78,Form("#Gamma_{Y(1750)} = %0.4f[GeV/c^{2}]",fmodel_phiy->GetParameter(2)));

  cphiy->Modified();
  cphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy.root", data.Data(), data.Data()),"root");
  cphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_phiy.eps", data.Data(), data.Data()),"eps");

  // TCanvas *c2phiy = new TCanvas("chphiy","c2phiy",700,400);
  // c2phiy->cd();
  // hphiy->SetLineStyle(1);
  // hphiy->SetLineWidth(2);
  // hphiy->SetMarkerStyle(20);
  // hphiy->SetMarkerSize(0.7);
  // hphiy->Draw();

  // w.factory("BreitWigner::sig1_phiy(m_phiy[1.50,3.20],mean1_phiy[1.730,1.50,1.80],width_phiy[0.01,1.])");
  // w.var("m_phiy")->setRange("mphiy_range",1.50,1.80);
  // // w.factory("Landau::model_phiy(m_phiy[1.50,3.20],mean2_phiy[1.50,3.20],sigma_phiy[0.01,1.])");
  // w.factory("Chebychev::bkg_phiy(m_phiy,{c0_phiy[-1.13752e+08],c1_phiy[1.41870e+08],c2_phiy[-3.95781e+07],c3_phiy[5.14407e+06],c4_phiy[-2.56040e+05]})");
  // w.factory("SUM:::model_phiy(nsig_phiy[0,100000000]*sig1_phiy, nbkg_phiy[0,100000000]*bkg_phiy)"); //, nsig_phiy[0,100000000]*sig2_phiy
  // RooDataHist dh1_phiyMass_kmkpcut("dh1_phiyMass_kmkpcut","dh1_phiyMass_kmkpcut",*w.var("m_phiy"),Import(*hphiy));
  // RooPlot* frh1_phiyMass_kmkpcut = w.var("m_phiy")->frame(Title("invariant mass K^{+}K^{-}#pi^{+}#pi^{-}"));
  // frh1_phiyMass_kmkpcut->SetTitleOffset(0.90,"X");
  // frh1_phiyMass_kmkpcut->SetTitleSize(0.05,"XYZ");
  // frh1_phiyMass_kmkpcut->SetLabelSize(0.05,"xyz");
  // w.pdf("model_phiy")->fitTo(dh1_phiyMass_kmkpcut);
  // RooFitResult* result_phiy = w.pdf("model_phiy")->fitTo(dh1_phiyMass_kmkpcut,Extended(kTRUE),Save());
  // dh1_phiyMass_kmkpcut.plotOn(frh1_phiyMass_kmkpcut,RooFit::Name("ldh1_phiyMass_kmkpcut"));
  // w.pdf("model_phiy")->plotOn(frh1_phiyMass_kmkpcut,RooFit::Name("lmodel_phiy"));
  // w.pdf("model_phiy")->plotOn(frh1_phiyMass_kmkpcut, Components(*w.pdf("bkg_phiy")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_phiy"));
  // w.pdf("model_phiy")->plotOn(frh1_phiyMass_kmkpcut, Components(*w.pdf("sig1_phiy")),LineColor(kRed),RooFit::Name("lsig1_phiy"));
  // // w.pdf("model_phiy")->plotOn(frh1_phiyMass_kmkpcut, Components(*w.pdf("sig2_phiy")),LineColor(kMagenta),RooFit::Name("lsig2_phiy"));
  // w.pdf("model_phiy")->paramOn(frh1_phiyMass_kmkpcut,Layout(0.55,0.98,0.92)); //,Parameters(RooArgSet((*w->var("nsig")),(*w->var("nbkg")),(*w->var("mean2")),(*w->var("sigma2"))))
  // frh1_phiyMass_kmkpcut->Draw();
  // TLegend *lh1_phiyMass_kmkpcut = new TLegend(0.13,0.6,0.3,0.85);
  // lh1_phiyMass_kmkpcut->SetFillColor(kWhite);
  // lh1_phiyMass_kmkpcut->SetLineColor(kWhite);
  // lh1_phiyMass_kmkpcut->AddEntry(frh1_phiyMass_kmkpcut->findObject("ldh1_phiyMass_kmkpcut"),"Data", "p");
  // lh1_phiyMass_kmkpcut->AddEntry(frh1_phiyMass_kmkpcut->findObject("lmodel_phiy"),"total","l");
  // lh1_phiyMass_kmkpcut->AddEntry(frh1_phiyMass_kmkpcut->findObject("lbkg_phiy"),"Background", "l");
  // lh1_phiyMass_kmkpcut->AddEntry(frh1_phiyMass_kmkpcut->findObject("lsig1_phiy"),"Gaussian", "l");
  // // lh1_phiyMass_kmkpcut->AddEntry(frh1_phiyMass_kmkpcut->findObject("lsig2_phiy"),"Gaussian", "l");
  // lh1_phiyMass_kmkpcut->Draw();
  // result_phiy->Print();

  //   TCanvas *cbwmeanphiy = new TCanvas("cbwmeanphiy","cbwmeanphiy",700,400);
  //   cbwmeanphiy->cd();
  //   hbwmeanphiy->SetTitleOffset(0.90,"X");
  //   hbwmeanphiy->SetTitleSize(0.05,"XYZ");
  //   hbwmeanphiy->SetLabelSize(0.05,"xyz");
  //   hbwmeanphiy->Draw();
  //   hbwmeanphiy->Write();
  //   cbwmeanphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwmeanphiy.root", data.Data(), data.Data()),"root");
  //   cbwmeanphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwmeanphiy.eps", data.Data(), data.Data()),"eps");
  //
  //   TCanvas *cbwwidthphiy = new TCanvas("cbwwidthphiy","cbwwidthphiy",700,400);
  //   cbwwidthphiy->cd();
  //   hbwwidthphiy->SetTitleOffset(0.90,"X");
  //   hbwwidthphiy->SetTitleSize(0.05,"XYZ");
  //   hbwwidthphiy->SetLabelSize(0.05,"xyz");
  //   hbwwidthphiy->Draw();
  //   hbwwidthphiy->Write();
  //   cbwwidthphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwwidthphiy.root", data.Data(), data.Data()),"root");
  //   cbwwidthphiy->Print(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/c_%s_bwwidthphiy.eps", data.Data(), data.Data()),"eps");

  // ************   Cross Section
  //------ Phi(1020)
  double eff_phi = 100*w.var("nsig_phi")->getVal()/500000; //nsig_phi
  double lumi = 5.6*1e09;
  double br_phi = 0.489;
  double br_fo = 1;
  double br_y = 1;
  double CrossSection_phi =1e06*w.var("nsig_phi_data")->getVal()/(eff_phi*lumi*br_phi*br_fo*br_y);
  // cout<<"#Sigma_{phi(1020)} = "<<CrossSection_phi<<"[b]"<<"N_phi.getVal() = "<<N_phi.getVal()<<endl;
  // TCanvas *c_CrossSection =new TCanvas("c_CrossSection","c_CrossSection",500,600);
  // c_CrossSection->cd();
  // TLatex lat;
  // lat.SetTextAlign(12);
  // lat.SetTextSize(0.04);
  // lat.DrawLatex(0.1,0.8,Form("#sigma_{#phi(1020)} = %0.7f[b]",CrossSection_phi));
  // c_CrossSection->Print("CrossSection.pdf");

  // Other options are CREATE (same as NEW), RECREATE (i.e. replace), UPDATE and READ.
  FILE *tab = fopen(Form("/data.local/nacer/halld_my/pimpipkmkp/fig_%s/tab.tex",data.Data()),"w");
  if (tab!=NULL) printf("File opened successfully\n"); // \\usepackage{epstopdf}\n

  fprintf(tab,"\\documentclass{article}\n \\usepackage{graphicx}\n \\usepackage{epstopdf}\n \\usepackage{amsmath}\n" "\\usepackage{caption}\n" "\\captionsetup{labelformat=empty}\n"
  "\\begin{document}\n"
  "\\begin{table}[[!b]]\n"
  "\\centering\n"
  "\\begin{tabular}{ |c|c|c|c|c|c| }\n"
  "\\hline\n"
  "$\\epsilon [\\%%]$ & $\\mathcal{L}$ & $BR[\\phi{\\rightarrow}K^{+}K^{-}]$ & $BR[f_{0}{\\rightarrow}\\pi^{+}\\pi^{-}]$ & $BR[Y{\\rightarrow}K^{+}K^{-}\\pi^{+}\\pi^{-}]$ & $\\sigma$ [$\\mu$b] \x5c\x5c \n"
  "\\hline\n"
  "%0.f & %0.f & %0.3f & %0.f & %0.f & %0.5f \x5c\x5c \n"
  "\\hline\n"
  "\\end{tabular}\n"
  "\\caption{Cross Section}\n"
  "\\end{table}\n"
  "\\end{document}\n",eff_phi,lumi,br_phi,br_fo,br_y,CrossSection_phi);

  fclose(tab);
  gSystem->Exec(Form("latex -shell-escape /data.local/nacer/halld_my/pimpipkmkp/fig_%s/tab.tex",data.Data()));
  gSystem->Exec(Form("dvips -E -o /data.local/nacer/halld_my/pimpipkmkp/fig_%s/tab.eps tab.dvi",data.Data()));

  outputfig->Print();
  // outputfig->Close();
}
