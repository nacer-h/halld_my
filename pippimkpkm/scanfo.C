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

void scanfo(int nkk=50, int n2pi=100, int ne=1, int nt=1)// TString cut="&& kin_chisq<30 && abs(mm2)<0.015") // && -t_kin<1 && beam_p4_kin.E()>6
{
  TFile *fdata = new TFile("/Users/nacer/halld_my/pippimkpkm/input/pippimkpkm_17v21.root");
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TFile *fmc = new TFile("/Users/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v03.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  TFile *fps = new TFile("/Users/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");  
  TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/scanfo.root","UPDATE");

  ofstream ofs_scanfo("tableau_fo.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min    = 0.98, mkk_max    = 1.2;
  double m2pi_min   = 0.3,  m2pi_max   = 1.2;

  // ======================================== fo vs. Eg ===============================================
  // root -l 'scanfo.C+(100,50,4,4)'
  TCanvas *cphifoe=new TCanvas("cphifoe","cphifoe",1500, 800);
  cphifoe->Divide(2,2);

  TCanvas *cphifoe1[ne];
  TCanvas *cphifoe2[ne];
  // TCanvas *cphifoe3[ne];
  // TCanvas *cphifoe4[ne];

  TCanvas *cgphifoe=new TCanvas("cgphifoe","cgphifoe",1500, 800);
  cgphifoe->Divide(2,2);
  TGraphErrors *gphifoe[ne];

  // TCanvas *cgnophifoe=new TCanvas("cgnophifoe","cgnophifoe",1500, 800);
  // TGraphErrors *gnophifoe = new TGraphErrors(n2pi);
  // gnophifoe->SetMarkerStyle(20);
  // gnophifoe->SetMarkerSize(1.0);
  // gnophifoe->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  double Egmin = 6.3; //6.3;//= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
	double Egmax = 11.6; //11.6;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
	double Egstep = (Egmax-Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    cphifoe->cd(j);
    TH2F *h2d2 = new TH2F("h2d2", Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    h2d2->Draw("colz");

    cphifoe1[j] = new TCanvas(Form("cphifoe1_%d", j), Form("cphifoe1_%d", j), 1500, 800);
    cphifoe1[j]->Divide(5, 5);
    cphifoe2[j] = new TCanvas(Form("cphifoe2_%d", j), Form("cphifoe2_%d", j), 1500, 800);
    cphifoe2[j]->Divide(5, 5);
    // cphifoe3[j] = new TCanvas(Form("cphifoe3_%d", j), Form("cphifoe3_%d", j), 1500, 800);
    // cphifoe3[j]->Divide(5, 5);
    // cphifoe4[j] = new TCanvas(Form("cphifoe4_%d", j), Form("cphifoe4_%d", j), 1500, 800);
    // cphifoe4[j]->Divide(5, 5);

    gphifoe[j] = new TGraphErrors(n2pi);
    gphifoe[j]->SetMarkerStyle(20);
    gphifoe[j]->SetMarkerSize(1.0);
    gphifoe[j]->SetMarkerColor(1);
    gphifoe[j]->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", Eg1[j], Eg2[j]));

    // gnophifoe->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{bkg}", Eg1[j], Eg2[j]));

    for (int i = 1; i <= n2pi; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphifoe1[j]->cd(i);
      if(i > 25 && i<51) cphifoe2[j]->cd(i-25);
      // if(i > 50 && i<76) cphifoe3[j]->cd(i-50);
      // if(i > 75) cphifoe4[j]->cd(i-75);

      TH1D *hphifoe_py = h2d2->ProjectionY(Form("_hphifoe_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifoe_py->Draw();
      // fb2->Draw("same");

      //w.factory(Form("Voigtian::sig_phifoe(m_phifoe[%f,%f],mean_phifoe[1.016,1.021],width_phifoe[0.004],sigma_phifoe[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phifo[0.001,0.01], mean_phifo[1.016,1.022]
      w.factory(Form("Voigtian::sig_phifoe(m_phifoe[%f,%f],mean_phifoe[1.016,1.03],width_phifoe[0.004],sigma_phifoe[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phifo[0.001,0.01], mean_phifo[1.016,1.022],50:mean_phifoe[1.016,1.032]
      w.factory("Chebychev::bkg_phifoe(m_phifoe,{c0_phifoe[-1,1], c1_phifoe[-1,1]})");//, c2_phifoe[-1,1]
      w.factory("SUM:::model_phifoe(nsig_phifoe[0,100000000]*sig_phifoe, nbkg_phifoe[0,100000000]*bkg_phifoe)"); //nsig[0,100000000]*sig2,
      w.var("m_phifoe")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phifoe("dh_phifoe", "dh_phifoe", *w.var("m_phifoe"), Import(*hphifoe_py));
      RooPlot *fr_phifoe = w.var("m_phifoe")->frame(Title(" "));
      // fr_phifoe->SetTitleOffset(0.90, "X");
      // fr_phifoe->SetTitleSize(0.06, "XYZ");
      // fr_phifoe->SetLabelSize(0.06, "xyz");
      w.pdf("model_phifoe")->fitTo(dh_phifoe);
      //cout<<"=========  no problem up to here ! =========="<<endl;

      // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phifoe.plotOn(fr_phifoe, RooFit::Name("ldh_phifoe"));
      w.pdf("model_phifoe")->plotOn(fr_phifoe, Components(*w.pdf("sig_phifoe")), LineColor(kRed), RooFit::Name("lsig_phifoe"));
      w.pdf("model_phifoe")->plotOn(fr_phifoe, Components(*w.pdf("bkg_phifoe")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phifoe"));
      w.pdf("model_phifoe")->plotOn(fr_phifoe, RooFit::Name("lmodel_phifoe"));
      // w.pdf("model_phifoe")->paramOn(fr_phifoe, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phifoe"), *w.var("nbkg_phifoe")))); //,*w.var("mean_phifoe"),*w.var("width_phifoe"),*w.var("sigma_phifoe"))));
      fr_phifoe->Draw();

      TLegend *l_phifoe = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phifoe->SetFillColor(kWhite);
      l_phifoe->SetLineColor(kWhite);
      // l_phifoe->AddEntry(fr_phifoe->findObject("ldh_phifoe"), "Data", "p");
      // l_phifoe->AddEntry(fr_phifoe->findObject("lmodel_phifoe"), "total", "l");
      l_phifoe->AddEntry(fr_phifoe->findObject("lsig_phifoe"), Form("N_{Sig} = %.2f", w.var("nsig_phifoe")->getVal()), "l");
      l_phifoe->AddEntry(fr_phifoe->findObject("lbkg_phifoe"), Form("N_{Bkg} = %.2f", w.var("nbkg_phifoe")->getVal()), "l");
      l_phifoe->Draw();

      // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      double N2 = w.var("nsig_phifoe")->getVal();
      double dN2 = w.var("nsig_phifoe")->getError();

      gphifoe[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifoe[j]->SetPointError(i - 1, 0, dN2);

      // // Integrate normalized pdf over subrange
      // w.var("m_phifoe")->setRange("sigregion",1.005,1.035);
      // // RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
      // RooAbsReal* fbkg = w.pdf("bkg_phifoe")->createIntegral(*w.var("m_phifoe"),NormSet(*w.var("m_phifoe")),Range("sigregion"));
      // RooAbsReal* fsig = w.pdf("sig_phifoe")->createIntegral(*w.var("m_phifoe"),NormSet(*w.var("m_phifoe")),Range("sigregion"));
      // // double Nbkg = -fbkg->getVal()*(w.var("nsig_phifoe")->getVal()+w.var("nbkg_phifoe")->getVal())+fsig->getVal()*w.var("nsig_phifoe")->getVal();       
      // double Nbkg = fbkg->getVal()*w.var("nbkg_phifoe")->getVal();
      // double dNbkg = 0;//fbkg->getError();
      // gnophifoe->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophifoe->SetPointError(i - 1, 0, dNbkg);
      // ofs_scanfo << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h2d2->GetYaxis()->GetBinCenter(" << i << ") = " << h2d2->GetYaxis()->GetBinCenter(i)<<endl;

      cgphifoe->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphifoe->cd(j);
    gphifoe[j]->Draw("AP");

    // cgnophifoe->cd();//j);
    // gnophifoe->Draw("Psame");
    // int j =1;
    // gphifoe->Write(Form("grphifoe_%d", j), TObject::kWriteDelete);

    cphifoe1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifoe_%d.root", j), "root");
    cphifoe1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifoe_%d.eps", j), "eps");
    cphifoe2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifoe_%d.root", j), "root");
    cphifoe2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifoe_%d.eps", j), "eps");
    // cphifoe3[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c3_phifoe_%d.root", j), "root");
    // cphifoe3[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c3_phifoe_%d.eps", j), "eps");
    // cphifoe4[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c4_phifoe_%d.root", j), "root");
    // cphifoe4[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c4_phifoe_%d.eps", j), "eps");
 
  }

  cphifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifoe.root", "root");
  cphifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifoe.eps", "eps");
  cgphifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe.root", "root");
  cgphifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifoe.eps", "eps");
  // cgnophifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifoe.root", "root");
  // cgnophifoe->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gnophifoe.eps", "eps");

  // ======================================== fo vs. -t ===============================================
/*
  TCanvas *cphifot=new TCanvas("cphifot","cphifot",1500, 800);
  cphifot->Divide(2,2);

  TCanvas *cphifot1[nt];
  TCanvas *cphifot2[nt];

  TCanvas *cgphifot=new TCanvas("cgphifot","cgphifot",1500, 800);
  cgphifot->Divide(2,2);
  TGraphErrors *gphifot[nt];

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

    cphifot->cd(j);
    TH2F *h2d2 = new TH2F("h2d2", Form("%.2f<-t<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi, m2pi_min, m2pi_max, nkk, mkk_min, mkk_max);
    tdata->Project("h2d2", "kpkm_mf:pippim_mf", Form("w8*((kpkm_uni || pippim_uni) && beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    h2d2->Draw("colz");

    cphifot1[j] = new TCanvas(Form("cphifot1_%d", j), Form("cphifot1_%d", j), 1500, 800);
    cphifot1[j]->Divide(5, 5);
    cphifot2[j] = new TCanvas(Form("cphifot2_%d", j), Form("cphifot2_%d", j), 1500, 800);
    cphifot2[j]->Divide(5, 5);

    gphifot[j] = new TGraphErrors(n2pi);
    gphifot[j]->SetMarkerStyle(20);
    gphifot[j]->SetMarkerSize(1.0);
    gphifot[j]->SetMarkerColor(1);
    gphifot[j]->SetMinimum(0.);
    gphifot[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));

    for (int i = 1; i <= n2pi; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphifot1[j]->cd(i);
      if(i > 25) cphifot2[j]->cd(i-25);

      TH1D *hphifot_py = h2d2->ProjectionY(Form("_hphifot_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphifot_py->Draw();
      // fb2->Draw("same");

      w.factory(Form("Voigtian::sig_phifot(m_phifot[%f,%f],mean_phifot[1.018,1.022],width_phifot[0.004],sigma_phifot[0.001,0.01])", mkk_min, mkk_max)); //sigma_phifo[0.001,0.01], mean_phifo[1.016,1.022]
      w.factory("Chebychev::bkg_phifot(m_phifot,{c0_phifot[-1,1], c1_phifot[-1,1], c2_phifot[-1,1]})");
      w.factory("SUM:::model_phifot(nsig_phifot[0,100000000]*sig_phifot, nbkg_phifot[0,100000000]*bkg_phifot)"); //nsig[0,100000000]*sig2,
      w.var("m_phifot")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      RooDataHist dh_phifot("dh_phifot", "dh_phifot", *w.var("m_phifot"), Import(*hphifot_py));
      RooPlot *fr_phifot = w.var("m_phifot")->frame(Title(" "));
      // fr_phifot->SetTitleOffset(0.90, "X");
      // fr_phifot->SetTitleSize(0.06, "XYZ");
      // fr_phifot->SetLabelSize(0.06, "xyz");
      w.pdf("model_phifot")->fitTo(dh_phifot);
      //cout<<"=========  no problem up to here ! =========="<<endl;

      // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      dh_phifot.plotOn(fr_phifot, RooFit::Name("ldh_phifot"));
      w.pdf("model_phifot")->plotOn(fr_phifot, Components(*w.pdf("sig_phifot")), LineColor(kRed), RooFit::Name("lsig_phifot"));
      w.pdf("model_phifot")->plotOn(fr_phifot, Components(*w.pdf("bkg_phifot")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phifot"));
      w.pdf("model_phifot")->plotOn(fr_phifot, RooFit::Name("lmodel_phifot"));
      // w.pdf("model_phifot")->paramOn(fr_phifot, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phifot"), *w.var("nbkg_phifot")))); //,*w.var("mean_phifot"),*w.var("width_phifot"),*w.var("sigma_phifot"))));
      fr_phifot->Draw();

      TLegend *l_phifot = new TLegend(0.5, 0.7, 0.8, 0.9);
      l_phifot->SetFillColor(kWhite);
      l_phifot->SetLineColor(kWhite);
      // l_phifot->AddEntry(fr_phifot->findObject("ldh_phifot"), "Data", "p");
      // l_phifot->AddEntry(fr_phifot->findObject("lmodel_phifot"), "total", "l");
      l_phifot->AddEntry(fr_phifot->findObject("lsig_phifot"), Form("N_{Sig} = %.2f", w.var("nsig_phifot")->getVal()), "l");
      l_phifot->AddEntry(fr_phifot->findObject("lbkg_phifot"), Form("N_{Bkg} = %.2f", w.var("nbkg_phifot")->getVal()), "l");
      l_phifot->Draw();

      // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      double N2 = w.var("nsig_phifot")->getVal();
      double dN2 = w.var("nsig_phifot")->getError();

      gphifot[j]->SetPoint(i - 1, h2d2->GetXaxis()->GetBinCenter(i), N2);
      gphifot[j]->SetPointError(i - 1, 0, dN2);

      cgphifot->Update();
      // c2->Update();
      //sleep(1);

      
    }

    // cout << endl;

    cgphifot->cd(j);
    gphifot[j]->Draw("AP");

    // int j =1;
    // gphifot->Write(Form("grphifot_%d", j), TObject::kWriteDelete);

    cphifot1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifot_%d.root", j), "root");
    cphifot1[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c1_phifot_%d.eps", j), "eps");
    cphifot2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifot_%d.root", j), "root");
    cphifot2[j]->Print(Form("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c2_phifot_%d.eps", j), "eps");
  
  }

  cphifot->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifot.root", "root");
  cphifot->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_phifot.eps", "eps");
  cgphifot->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifot.root", "root");
  cgphifot->Print("/Users/nacer/halld_my/pippimkpkm/fig_scanfo/c_gphifot.eps", "eps");
*/
  outputfig->Print();

}
