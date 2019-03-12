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

void scany(int nkk=100, int n2pi2k=50, int ne=1, int nt=1, TString cut="&& kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fdata = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm_17v21.root");
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TFile *fmc = new TFile("/Users/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v03.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TFile *fps = new TFile("/Users/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");  
  TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_scany/scany.root","UPDATE");

  ofstream ofs_scany("tableau_y.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min    = 0.98, mkk_max    = 1.2;
  double m2pi2k_min = 1.6,  m2pi2k_max = 3.2;

  // ======================================== Y vs. Eg ===============================================

  TCanvas *cphiye=new TCanvas("cphiye","cphiye",1500, 800);
  // cphiye->Divide(2,2);

  TCanvas *cphiye1[ne];
  TCanvas *cphiye2[ne];
  TCanvas *cphiye3[ne];
  TCanvas *cphiye4[ne];

  TCanvas *cgphiye=new TCanvas("cgphiye","cgphiye",1500, 800);
  // cgphiye->Divide(2,2);
  TGraphErrors *gphiye[ne];

  TCanvas *cgnophiye=new TCanvas("cgnophiye","cgnophiye",1500, 800);
  TGraphErrors *gnophiye= new TGraphErrors(n2pi2k);
  gnophiye->SetMarkerStyle(20);
  gnophiye->SetMarkerSize(1.0);
  gnophiye->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  double Egmin = 6.3;//= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
	double Egmax = 11.6;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
	double Egstep = (Egmax-Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    cphiye->cd();//j);
    TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<E_{#gamma}<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    h1d2->Draw("colz");

    cphiye1[j] = new TCanvas(Form("cphiye1_%d", j), Form("cphiye1_%d", j), 1500, 800);
    cphiye1[j]->Divide(5, 5);
    cphiye2[j] = new TCanvas(Form("cphiye2_%d", j), Form("cphiye2_%d", j), 1500, 800);
    cphiye2[j]->Divide(5, 5);
    cphiye3[j] = new TCanvas(Form("cphiye3_%d", j), Form("cphiye3_%d", j), 1500, 800);
    cphiye3[j]->Divide(5, 5);
    cphiye4[j] = new TCanvas(Form("cphiye4_%d", j), Form("cphiye4_%d", j), 1500, 800);
    cphiye4[j]->Divide(5, 5);

    gphiye[j] = new TGraphErrors(n2pi2k);
    gphiye[j]->SetMarkerStyle(20);
    gphiye[j]->SetMarkerSize(1.0);
    gphiye[j]->SetMarkerColor(1);
    gphiye[j]->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", Eg1[j], Eg2[j]));

    gnophiye->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1[j], Eg2[j]));


    for (int i = 1; i <= n2pi2k; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphiye1[j]->cd(i);
      if(i > 25 && i<51) cphiye2[j]->cd(i-25);
      if(i > 50 && i<76) cphiye3[j]->cd(i-50);
      if(i > 75) cphiye4[j]->cd(i-75);

      TH1D *hphiye_py = h1d2->ProjectionY(Form("_hphiye_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiye_py->Draw();
      // fb2->Draw("same");

      //w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.016,1.022],width_phiye[0.004],sigma_phiye[0.001,0.6])", mkk_min, mkk_max)); //1.018,1.022
      w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.01,1.025],width_phiye[0.004],sigma_phiye[0.0001,0.01])", mkk_min, mkk_max)); //50:mean_phiye[1.01,1.025],sigma_phiye[0.0001,0.01],100:mean_phiye[1.018,1.022] 
      w.factory("Chebychev::bkg_phiye(m_phiye,{c0_phiye[-1,1], c1_phiye[-1,1], c2_phiye[-1,1]})");//, c2_phiye[-1,1]
      w.factory("SUM:::model_phiye(nsig_phiye[0,100000000]*sig_phiye, nbkg_phiye[0,100000000]*bkg_phiye)"); //nsig[0,100000000]*sig2,
      w.var("m_phiye")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phiye("dh_phiye", "dh_phiye", *w.var("m_phiye"), Import(*hphiye_py));
      RooPlot *fr_phiye = w.var("m_phiye")->frame(Title(" "));
      // fr_phiye->SetTitleOffset(0.90, "X");
      // fr_phiye->SetTitleSize(0.06, "XYZ");
      // fr_phiye->SetLabelSize(0.06, "xyz");
      w.pdf("model_phiye")->fitTo(dh_phiye);
      //cout<<"=========  no problem up to here ! =========="<<endl;

      // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phiye.plotOn(fr_phiye, RooFit::Name("ldh_phiye"));
      w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("sig_phiye")), LineColor(kRed), RooFit::Name("lsig_phiye"));
      w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("bkg_phiye")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiye"));
      w.pdf("model_phiye")->plotOn(fr_phiye, RooFit::Name("lmodel_phiye"));
      // w.pdf("model_phiye")->paramOn(fr_phiye, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiye"), *w.var("nbkg_phiye")))); //,*w.var("mean_phiye"),*w.var("width_phiye"),*w.var("sigma_phiye"))));
      fr_phiye->Draw();

      TLegend *l_phiye = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phiye->SetFillColor(kWhite);
      l_phiye->SetLineColor(kWhite);
      // l_phiye->AddEntry(fr_phiye->findObject("ldh_phiye"), "Data", "p");
      // l_phiye->AddEntry(fr_phiye->findObject("lmodel_phiye"), "total", "l");
      l_phiye->AddEntry(fr_phiye->findObject("lsig_phiye"), Form("N_{Sig} = %.2f", w.var("nsig_phiye")->getVal()), "l");
      l_phiye->AddEntry(fr_phiye->findObject("lbkg_phiye"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiye")->getVal()), "l");
      l_phiye->Draw();

      // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      double N2 = w.var("nsig_phiye")->getVal();
      double dN2 = w.var("nsig_phiye")->getError();      

      gphiye[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N2);
      gphiye[j]->SetPointError(i - 1, 0, dN2);

      // Integrate normalized pdf over subrange
      w.var("m_phiye")->setRange("sigregion",1.005,1.035);
      // RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
      RooAbsReal* fbkg = w.pdf("bkg_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      RooAbsReal* fsig = w.pdf("sig_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      // double Nbkg = -fbkg->getVal()*(w.var("nsig_phiye")->getVal()+w.var("nbkg_phiye")->getVal())+fsig->getVal()*w.var("nsig_phiye")->getVal();       
      double Nbkg = fbkg->getVal()*w.var("nbkg_phiye")->getVal();
      double dNbkg = 0;//fbkg->getError(); // w.var("nbkg_phiye")->getError();
      gnophiye->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
      gnophiye->SetPointError(i - 1, 0, dNbkg);
      ofs_scany << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

      cgphiye->Update();
      // c2->Update();
      //sleep(1);
    }

    cout << endl;

    cgphiye->cd();//j);
    gphiye[j]->Draw("AP");

    // cgnophiye->cd();//j);
    gnophiye->Draw("Psame");

    // int j =1;
    // gphiye->Write(Form("grphiye_%d", j), TObject::kWriteDelete);

    cphiye1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c1_phiye_%d.root", j), "root");
    cphiye1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c1_phiye_%d.eps", j), "eps");
    cphiye2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c2_phiye_%d.root", j), "root");
    cphiye2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c2_phiye_%d.eps", j), "eps");
    cphiye3[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c3_phiye_%d.root", j), "root");
    cphiye3[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c3_phiye_%d.eps", j), "eps");
    cphiye4[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c4_phiye_%d.root", j), "root");
    cphiye4[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c4_phiye_%d.eps", j), "eps");
  
  }

  cphiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_phiye.root", "root");
  cphiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_phiye.eps", "eps");
  cgphiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gphiye.root", "root");
  cgphiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gphiye.eps", "eps");
  cgnophiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gnophiye.root", "root");
  cgnophiye->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gnophiye.eps", "eps");

  // ======================================== Y vs. -t ===============================================
/*
  TCanvas *cphiyt=new TCanvas("cphiyt","cphiyt",1500, 800);
  cphiyt->Divide(3,4);

  TCanvas *cphiyt1[nt];

  TCanvas *cgphiyt=new TCanvas("cgphiyt","cgphiyt",1500, 800);
  cgphiyt->Divide(3,4);
  TGraphErrors *gphiyt[nt];

  double tmin = 0.1;//= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
	double tmax = 2;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
	double tstep = (tmax-tmin) / nt;
  double t1[nt];
  double t2[nt];

  for (int j = 1; j <= nt; ++j)
  {
    t1[j] = tmin + ((j - 1) * tstep); //i * tstep;
    t2[j] = tmin + (j * tstep);
    cout << "########  j = " << j << " | t1[" << j << "] = " << t1[j] << " | t2[" << j << "] = " << t2[j] << endl;

    cphiyt->cd(j);
    TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<-t<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    h1d2->Draw("colz");

    cphiyt1[j] = new TCanvas(Form("cphiyt1_%d", j), Form("cphiyt1_%d", j), 1500, 800);
    cphiyt1[j]->Divide(5, 4);

    gphiyt[j] = new TGraphErrors(n2pi2k);
    gphiyt[j]->SetMarkerStyle(20);
    gphiyt[j]->SetMarkerSize(1.0);
    gphiyt[j]->SetMarkerColor(1);
    gphiyt[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));

    for (int i = 1; i <= n2pi2k; ++i)
    {
      cout << i << " " << flush;

      cphiyt1[j]->cd(i);

      TH1D *hphiyt_py = h1d2->ProjectionY(Form("_hphiyt_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiyt_py->Draw();
      // fb2->Draw("same");

      w.factory(Form("Voigtian::sig_phiyt(m_phiyt[%f,%f],mean_phiyt[1.013,1.03],width_phiyt[0.004],sigma_phiyt[0.001,0.1])", mkk_min, mkk_max)); //sigma_phiy[0.001,0.01], mean_phiy[1.016,1.022]
      w.factory("Chebychev::bkg_phiyt(m_phiyt,{c0_phiyt[-1,1], c1_phiyt[-1,1]})");//, c2_phiyt[-1,1]
      w.factory("SUM:::model_phiyt(nsig_phiyt[0,100000000]*sig_phiyt, nbkg_phiyt[0,100000000]*bkg_phiyt)"); //nsig[0,100000000]*sig2,
      w.var("m_phiyt")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phiyt("dh_phiyt", "dh_phiyt", *w.var("m_phiyt"), Import(*hphiyt_py));
      RooPlot *fr_phiyt = w.var("m_phiyt")->frame(Title(" "));
      // fr_phiyt->SetTitleOffset(0.90, "X");
      // fr_phiyt->SetTitleSize(0.06, "XYZ");
      // fr_phiyt->SetLabelSize(0.06, "xyz");
      w.pdf("model_phiyt")->fitTo(dh_phiyt);
      //cout<<"=========  no problem up to here ! =========="<<endl;

      // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phiyt.plotOn(fr_phiyt, RooFit::Name("ldh_phiyt"));
      w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("sig_phiyt")), LineColor(kRed), RooFit::Name("lsig_phiyt"));
      w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("bkg_phiyt")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiyt"));
      w.pdf("model_phiyt")->plotOn(fr_phiyt, RooFit::Name("lmodel_phiyt"));
      // w.pdf("model_phiyt")->paramOn(fr_phiyt, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiyt"), *w.var("nbkg_phiyt")))); //,*w.var("mean_phiyt"),*w.var("width_phiyt"),*w.var("sigma_phiyt"))));
      fr_phiyt->Draw();

      TLegend *l_phiyt = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phiyt->SetFillColor(kWhite);
      l_phiyt->SetLineColor(kWhite);
      // l_phiyt->AddEntry(fr_phiyt->findObject("ldh_phiyt"), "Data", "p");
      // l_phiyt->AddEntry(fr_phiyt->findObject("lmodel_phiyt"), "total", "l");
      l_phiyt->AddEntry(fr_phiyt->findObject("lsig_phiyt"), Form("N_{Sig} = %.2f", w.var("nsig_phiyt")->getVal()), "l");
      l_phiyt->AddEntry(fr_phiyt->findObject("lbkg_phiyt"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiyt")->getVal()), "l");
      l_phiyt->Draw();

      // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      double N2 = w.var("nsig_phiyt")->getVal();
      double dN2 = 0; //w.var("nsig_phiyt")->getError();

      gphiyt[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N2);
      gphiyt[j]->SetPointError(i - 1, 0, dN2);

      //cgphiyt->Update();
      // c2->Update();
      //sleep(1);

      
    }

    cout << endl;

    cgphiyt->cd(j);
    gphiyt[j]->Draw("AP");

    // int j =1;
    // gphiyt->Write(Form("grphiyt_%d", j), TObject::kWriteDelete);

    cphiyt1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c1_phiyt_%d.root", j), "root");
    cphiyt1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scany/c1_phiyt_%d.eps", j), "eps");
  }

  cphiyt->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_phiyt.root", "root");
  cphiyt->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_phiyt.eps", "eps");
  cgphiyt->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gphiyt.root", "root");
  cgphiyt->Print("/Users/nacer/halld_my/pippimkpkm/fig_scany/c_gphiyt.eps", "eps");
*/
  outputfig->Print();

}
