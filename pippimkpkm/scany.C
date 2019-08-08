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

void scany(int nkk=100, int n2pi2k=100, int ne=1, int nt=1) // , TString cut="&& kin_chisq<30 && abs(mm2)<0.015"
{
  TFile *fdata = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/pippimkpkm_17v21.root");
  TTree *tdata = (TTree*)fdata->Get("ntp");
  TFile *fmc = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v3.root");
  TTree *tmc = (TTree *)fmc->Get("ntp");
  // TFile *fps = new TFile("/data.local/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");  
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_scany/scany.root","UPDATE");

  ofstream ofs_scany("tableau_y.txt", ofstream::out);

  RooWorkspace w("w", kTRUE);

  double mkk_min    = 0.99, mkk_max    = 1.2;
  double m2pi2k_min = 1.7,  m2pi2k_max = 3.2;

  // gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
/*
  // ======================================== Y vs. Eg ===============================================
  // root -l 'scany.C+(50,50,4,4)'
  // Eg= 3, 7.9, 8.56, 9.66, 12
  TCanvas *cphiye=new TCanvas("cphiye","cphiye",1500, 800);//1500, 800
  cphiye->Divide(3,2);

  TCanvas *cphiye1[ne];
  TCanvas *cphiye2[ne];
  // TCanvas *cphiye3[ne];
  // TCanvas *cphiye4[ne];

  TCanvas *cgphiye=new TCanvas("cgphiye","cgphiye",1500, 800);//1500, 800
  cgphiye->Divide(3,2);
  TGraphErrors *gphiye[ne];

  // TCanvas *cgnophiye=new TCanvas("cgnophiye","cgnophiye",1500, 800);
  // TGraphErrors *gnophiye= new TGraphErrors;//(n2pi2k)
  // gnophiye->SetMarkerStyle(20);
  // gnophiye->SetMarkerSize(1.0);
  // gnophiye->SetMarkerColor(2);

  // tdata->SetAlias("w4","((abs(delta_t)<2.004)*1.25-0.25)");

  double Egmin = 3; //6.3 = hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double Egmax = 12; //11.6 = hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double Egstep = (Egmax-Egmin) / ne;
  double Eg1[ne];
  double Eg2[ne];

  for (int j = 1; j <= ne; ++j)
  {
    Eg1[j] = Egmin + ((j - 1) * Egstep); //i * Egstep;
    Eg2[j] = Egmin + (j * Egstep);
    cout << "########  j = " << j << " | Eg1[" << j << "] = " << Eg1[j] << " | Eg2[" << j << "] = " << Eg2[j] << endl;

    cphiye->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<E_{#gamma}<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", Eg1[j], Eg2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // TH2F *h1d2 = new TH2F("h1d2","(Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>%f && beam_p4_kin.E()<%f)", Eg1[j], Eg2[j]));
    // 4: 6, 8.1, 8.7, 9.3, 12
    // 3: 6, 8.3, 9.6, 12  
    TH2F *h1d2;

    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "3<E_{#gamma}<7.6 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<7.6)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "7.6<E_{#gamma}<8.2 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>7.6 && beam_p4_kin.E()<8.2)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "8.2<E_{#gamma}<8.5 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.2 && beam_p4_kin.E()<8.5)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "8.5<E_{#gamma}<8.85 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.5 && beam_p4_kin.E()<8.85)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "8.85<E_{#gamma}<10 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.85 && beam_p4_kin.E()<10)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "10<E_{#gamma}<12 GeV (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>10 && beam_p4_kin.E()<12)");
    }

    h1d2->Draw("colz");

    cphiye1[j] = new TCanvas(Form("cphiye1_%d", j), Form("cphiye1_%d", j), 1500, 800);//1500, 800
    cphiye1[j]->Divide(5, 5);
    cphiye2[j] = new TCanvas(Form("cphiye2_%d", j), Form("cphiye2_%d", j), 1500, 800);
    cphiye2[j]->Divide(5, 5);
    // cphiye3[j] = new TCanvas(Form("cphiye3_%d", j), Form("cphiye3_%d", j), 1500, 800);
    // cphiye3[j]->Divide(5, 5);
    // cphiye4[j] = new TCanvas(Form("cphiye4_%d", j), Form("cphiye4_%d", j), 1500, 800);
    // cphiye4[j]->Divide(5, 5);

    gphiye[j] = new TGraphErrors();//n2pi2k
    gphiye[j]->SetMarkerStyle(20);
    gphiye[j]->SetMarkerSize(1.0);
    gphiye[j]->SetMarkerColor(1);
    // gphiye[j]->SetTitle("(Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    // gphiye[j]->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", Eg1[j], Eg2[j]));
    if(j==1) gphiye[j]->SetTitle("3<E_{#gamma}<7.6 (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiye[j]->SetTitle("7.6<E_{#gamma}<8.2 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiye[j]->SetTitle("8.2<E_{#gamma}<8.5 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiye[j]->SetTitle("8.5<E_{#gamma}<8.85 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiye[j]->SetTitle("8.85<E_{#gamma}<10 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiye[j]->SetTitle("10<E_{#gamma}<12 GeV (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    // gnophiye->SetTitle(Form("%.2f<E_{#gamma}<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{NO #phi}", Eg1[j], Eg2[j]));


    for (int i = 1; i <=n2pi2k; ++i)//n2pi2k
    {
      // cout << i << " " << flush;

      if(i < 26) cphiye1[j]->cd(i);//i
      if(i > 25 && i<51) cphiye2[j]->cd(i-25);
      // if(i > 50 && i<76) cphiye3[j]->cd(i-50);
      // if(i > 75) cphiye4[j]->cd(i-75);

      TH1D *hphiye_py = h1d2->ProjectionY(Form("_hphiye_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiye_py->Draw();
      if (hphiye_py->Integral(1, 50) < 100)
        continue;      
      // fb2->Draw("same");

      // // w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.022],width_phiye[0.004],sigma_phiye[0.00001,0.01])", mkk_min, mkk_max)); //1.018,1.022
      // w.factory(Form("BreitWigner::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.021],width_phiye[0.008, 0.010])", mkk_min, mkk_max));
      // // w.factory(Form("Voigtian::sig_phiye(m_phiye[%f,%f],mean_phiye[1.018,1.022],width_phiye[0.004],sigma_phiye[0.0001,0.01])", mkk_min, mkk_max)); //50:mean_phiye[1.01,1.025],sigma_phiye[0.0001,0.01],100:mean_phiye[1.018,1.022] 
      // w.factory("Chebychev::bkg_phiye(m_phiye,{c0_phiye[-10,10], c1_phiye[-10,10], c2_phiye[-10,10]})");//
      // w.factory("SUM:::model_phiye(nsig_phiye[0,100000000]*sig_phiye, nbkg_phiye[0,100000000]*bkg_phiye)"); //nsig[0,100000000]*sig2,
      // w.var("m_phiye")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phiye("dh_phiye", "dh_phiye", *w.var("m_phiye"), Import(*hphiye_py));
      // RooPlot *fr_phiye = w.var("m_phiye")->frame(Title(" "));
      // // fr_phiye->SetTitleOffset(0.90, "X");
      // // fr_phiye->SetTitleSize(0.06, "XYZ");
      // // fr_phiye->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phiye")->fitTo(dh_phiye);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phiye.plotOn(fr_phiye, RooFit::Name("ldh_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("sig_phiye")), LineColor(kRed), RooFit::Name("lsig_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, Components(*w.pdf("bkg_phiye")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiye"));
      // w.pdf("model_phiye")->plotOn(fr_phiye, RooFit::Name("lmodel_phiye"));
      // // w.pdf("model_phiye")->paramOn(fr_phiye, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiye"), *w.var("nbkg_phiye")))); //,*w.var("mean_phiye"),*w.var("width_phiye"),*w.var("sigma_phiye"))));
      // fr_phiye->Draw();

      // TLegend *l_phiye = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phiye->SetFillColor(kWhite);
      // l_phiye->SetLineColor(kWhite);
      // // l_phiye->AddEntry(fr_phiye->findObject("ldh_phiye"), "Data", "p");
      // // l_phiye->AddEntry(fr_phiye->findObject("lmodel_phiye"), "total", "l");
      // l_phiye->AddEntry(fr_phiye->findObject("lsig_phiye"), Form("N_{Sig} = %.2f", w.var("nsig_phiye")->getVal()), "l");
      // l_phiye->AddEntry(fr_phiye->findObject("lbkg_phiye"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiye")->getVal()), "l");
      // l_phiye->Draw();

      // // double N1  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // // double dN1 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      // double N1 = w.var("nsig_phiye")->getVal();
      // double dN1 = w.var("nsig_phiye")->getError();

      // // // Integrate normalized pdf over subrange
      // // w.var("m_phiye")->setRange("sigregion",1.005,1.035);
      // // // RooAbsReal* igx_sig = gx.createIntegral(x,NormSet(x),Range("signal")) ;
      // // RooAbsReal* fbkg = w.pdf("bkg_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      // // RooAbsReal* fsig = w.pdf("sig_phiye")->createIntegral(*w.var("m_phiye"),NormSet(*w.var("m_phiye")),Range("sigregion"));
      // // // double Nbkg = -fbkg->getVal()*(w.var("nsig_phiye")->getVal()+w.var("nbkg_phiye")->getVal())+fsig->getVal()*w.var("nsig_phiye")->getVal();       
      // // double Nbkg = fbkg->getVal()*w.var("nbkg_phiye")->getVal();
      // // double dNbkg = 0;//fbkg->getError(); // w.var("nbkg_phiye")->getError();
      
      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", 0.98, 1.2);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01      

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", 0.98, 1.2);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//3
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiye_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiye_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiye[j]->RemovePoint(i-1);

      gphiye[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiye[j]->SetPointError(i - 1, 0, dN1);

      // gnophiye->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), Nbkg);
      // gnophiye->SetPointError(i - 1, 0, dNbkg);
      // ofs_scany << " i = " << i << " | Nbkg = " << Nbkg << " | dNbkg = " << dNbkg << " | h1d2->GetYaxis()->GetBinCenter(" << i << ") = " << h1d2->GetYaxis()->GetBinCenter(i)<<endl;

      TLatex lat_phiye;
      lat_phiye.SetTextSize(0.09);
      lat_phiye.SetTextAlign(13); //align at top
      lat_phiye.SetNDC();
      lat_phiye.SetTextColor(kBlue);
      lat_phiye.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phiye.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiye.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phiye.DrawLatex(0.45, 0.58, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phiye.DrawLatex(0.45, 0.48, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      // lat.DrawLatex(0.6, 0.65, Form("N1 = %f",N1));

      hphiye_py->Write();
      cgphiye->Update();
      // c2->Update();
      //sleep(1);
    }

    // cout << endl;

    cgphiye->cd(j);
    gphiye[j]->Draw("AP");

    // cgnophiye->cd();//j);
    // gnophiye->Draw("AP");

    // int j =1;
    // gphiye->Write(Form("grphiye_%d", j), TObject::kWriteDelete);

    cphiye1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiye_%d.root", j), "root");
    cphiye1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiye_%d.eps", j), "eps");
    cphiye2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiye_%d.root", j), "root");
    cphiye2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiye_%d.eps", j), "eps");
    // cphiye3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c3_phiye_%d.root", j), "root");
    // cphiye3[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c3_phiye_%d.eps", j), "eps");
    // cphiye4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c4_phiye_%d.root", j), "root");
    // cphiye4[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c4_phiye_%d.eps", j), "eps");
  
  }

  cphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiye.root", "root");
  cphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiye.eps", "eps");
  cgphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiye.root", "root");
  cgphiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiye.eps", "eps");
  // cgnophiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gnophiye.root", "root");
  // cgnophiye->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gnophiye.eps", "eps");
*/
/*
  // ======================================== Y vs. -t ===============================================

  TCanvas *cphiyt=new TCanvas("cphiyt","cphiyt", 10, 10, 1500, 800);
  cphiyt->Divide(3,2);

  TCanvas *cphiyt1[nt];
  TCanvas *cphiyt2[nt];

  TCanvas *cgphiyt=new TCanvas("cgphiyt","cgphiyt", 10, 10, 1500, 800);
  cgphiyt->Divide(3,2);
  TGraphErrors *gphiyt[nt];

  double tmin = 0; // 0.1 = hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double tmax = 10;   // 2 = hdata_postcut->GetXaxis()->GetBinUpEdge(600);
  double tstep = (tmax-tmin) / nt;
  double t1[nt];
  double t2[nt];

  for (int j = 1; j <= nt; ++j)
  {
    t1[j] = tmin + ((j - 1) * tstep); //i * tstep;
    t2[j] = tmin + (j * tstep);
    cout << "########  j = " << j << " | t1[" << j << "] = " << t1[j] << " | t2[" << j << "] = " << t2[j] << endl;

    cphiyt->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<-t<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    TH2F *h1d2;
    // 3: 0, 0.7, 1.5, 4
    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "0<-t<0.45 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0 && -t_kin<0.45)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "0.45<-t<0.73 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0.45 && -t_kin<0.73)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "0.73<-t<1.1 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>0.73 && -t_kin<1.1)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "1.1<-t<1.63 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>1.1 && -t_kin<1.63)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "1.63<-t<2.63 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>1.63 && -t_kin<2.63)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "2.63<-t<10 GeV^{2} (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && -t_kin>2.63 && -t_kin<10)");
    }

    h1d2->Draw("colz");

    cphiyt1[j] = new TCanvas(Form("cphiyt1_%d", j), Form("cphiyt1_%d", j), 1500, 800);
    cphiyt1[j]->Divide(5, 5);
    cphiyt2[j] = new TCanvas(Form("cphiyt2_%d", j), Form("cphiyt2_%d", j), 1500, 800);
    cphiyt2[j]->Divide(5, 5);

    gphiyt[j] = new TGraphErrors();//n2pi2k
    gphiyt[j]->SetMarkerStyle(20);
    gphiyt[j]->SetMarkerSize(1.0);
    gphiyt[j]->SetMarkerColor(1);
    // gphiyt[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphiyt[j]->SetTitle("0<-t<0.45 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiyt[j]->SetTitle("0.45<-t<0.73 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiyt[j]->SetTitle("0.73<-t<1.1 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiyt[j]->SetTitle("1.1<-t<1.63 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiyt[j]->SetTitle("1.63<-t<2.63 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiyt[j]->SetTitle("2.63<-t<10 GeV^{2} (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi2k; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphiyt1[j]->cd(i);
      if(i > 25 && i<51) cphiyt2[j]->cd(i-25);

      TH1D *hphiyt_py = h1d2->ProjectionY(Form("_hphiyt_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiyt_py->Draw();
      if (hphiyt_py->Integral(1, 50) < 100)
        continue;
      // fb2->Draw("same");

      // w.factory(Form("Voigtian::sig_phiyt(m_phiyt[%f,%f],mean_phiyt[1.018,1.022],width_phiyt[0.004],sigma_phiyt[0.0001,0.01])", mkk_min, mkk_max)); //sigma_phiy[0.001,0.01], mean_phiy[1.016,1.022]
      // w.factory("Chebychev::bkg_phiyt(m_phiyt,{c0_phiyt[-1,1], c1_phiyt[-1,1]})");//, c2_phiyt[-1,1]
      // w.factory("SUM:::model_phiyt(nsig_phiyt[0,100000000]*sig_phiyt, nbkg_phiyt[0,100000000]*bkg_phiyt)"); //nsig[0,100000000]*sig2,
      // w.var("m_phiyt")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
      // RooDataHist dh_phiyt("dh_phiyt", "dh_phiyt", *w.var("m_phiyt"), Import(*hphiyt_py));
      // RooPlot *fr_phiyt = w.var("m_phiyt")->frame(Title(" "));
      // // fr_phiyt->SetTitleOffset(0.90, "X");
      // // fr_phiyt->SetTitleSize(0.06, "XYZ");
      // // fr_phiyt->SetLabelSize(0.06, "xyz");
      // w.pdf("model_phiyt")->fitTo(dh_phiyt);
      // //cout<<"=========  no problem up to here ! =========="<<endl;

      // // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
      // dh_phiyt.plotOn(fr_phiyt, RooFit::Name("ldh_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("sig_phiyt")), LineColor(kRed), RooFit::Name("lsig_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, Components(*w.pdf("bkg_phiyt")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_phiyt"));
      // w.pdf("model_phiyt")->plotOn(fr_phiyt, RooFit::Name("lmodel_phiyt"));
      // // w.pdf("model_phiyt")->paramOn(fr_phiyt, Layout(0.4, 0.90, 0.99), Parameters(RooArgSet(*w.var("nsig_phiyt"), *w.var("nbkg_phiyt")))); //,*w.var("mean_phiyt"),*w.var("width_phiyt"),*w.var("sigma_phiyt"))));
      // fr_phiyt->Draw();

      // TLegend *l_phiyt = new TLegend(0.5, 0.7, 0.8, 0.9);
      // l_phiyt->SetFillColor(kWhite);
      // l_phiyt->SetLineColor(kWhite);
      // // l_phiyt->AddEntry(fr_phiyt->findObject("ldh_phiyt"), "Data", "p");
      // // l_phiyt->AddEntry(fr_phiyt->findObject("lmodel_phiyt"), "total", "l");
      // l_phiyt->AddEntry(fr_phiyt->findObject("lsig_phiyt"), Form("N_{Sig} = %.2f", w.var("nsig_phiyt")->getVal()), "l");
      // l_phiyt->AddEntry(fr_phiyt->findObject("lbkg_phiyt"), Form("N_{Bkg} = %.2f", w.var("nbkg_phiyt")->getVal()), "l");
      // l_phiyt->Draw();

      // // double N2  = fs->Integral(0.99,1.12)/hslice2->GetBinWidth(1);
      // // double dN2 = N2*fsb->GetParError(4)/fsb->GetParameter(4);

      // double N2 = w.var("nsig_phiyt")->getVal();
      // double dN2 = w.var("nsig_phiyt")->getError();

      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 10000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01 

      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//4
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiyt_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiyt_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiyt[j]->RemovePoint(i-1);

      gphiyt[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiyt[j]->SetPointError(i - 1, 0, dN1);

      TLatex lat_phiyt;
      lat_phiyt.SetTextSize(0.09);
      lat_phiyt.SetTextAlign(13); //align at top
      lat_phiyt.SetNDC();
      lat_phiyt.SetTextColor(kBlue);
      lat_phiyt.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f",fsb->GetChisquare()/fsb->GetNDF()));
      lat_phiyt.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiyt.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f",fsb->GetParameter(1),fsb->GetParError(1)));
      lat_phiyt.DrawLatex(0.45, 0.58, Form("#Gamma = %0.3f#pm%0.3f",fsb->GetParameter(2),fsb->GetParError(2)));
      lat_phiyt.DrawLatex(0.45, 0.48, Form("#sigma = %0.3f#pm%0.3f",fsb->GetParameter(3),fsb->GetParError(3)));

      //cgphiyt->Update();
      // c2->Update();
      //sleep(1);

      
    }

    // cout << endl;

    cgphiyt->cd(j);
    gphiyt[j]->Draw("AP");

    // int j =1;
    // gphiyt->Write(Form("grphiyt_%d", j), TObject::kWriteDelete);

    cphiyt1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiyt_%d.root", j), "root");
    cphiyt1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiyt_%d.eps", j), "eps");
    cphiyt2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiyt_%d.root", j), "root");
    cphiyt2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiyt_%d.eps", j), "eps");  
  }

  cphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiyt.root", "root");
  cphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiyt.eps", "eps");
  cgphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiyt.root", "root");
  cgphiyt->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiyt.eps", "eps");
*/

  // ======================================== Y vs. (Eg,-t) ===============================================

  TCanvas *cphiy=new TCanvas("cphiy","cphiy", 10, 10, 1500, 800);
  cphiy->Divide(3,2);

  TCanvas *cphiy1[ne*nt];
  TCanvas *cphiy2[ne*nt];

  TCanvas *cgphiy=new TCanvas("cgphiy","cgphiy", 10, 10, 1500, 800);
  cgphiy->Divide(3,2);
  TGraphErrors *gphiy[ne*nt];


  for (int j = 1; j <= ne*nt; ++j)
  {
    cout << "########  j = " << j <<endl;

    cphiy->cd(j);
    // TH2F *h1d2 = new TH2F("h1d2", Form("%.2f<-t<%.2f (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", t1[j], t2[j]), n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
    // tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", Form("w8*((kpkm_uni || kpkmpippim_uni)" + cut + "&& beam_p4_kin.E()>6.3 && beam_p4_kin.E()<11.6 && -t_kin>%f && -t_kin<%f)", t1[j], t2[j]));
    TH2F *h1d2;
    // 3: 0, 0.7, 1.5, 4
    if (j == 1)
    {
      h1d2 = new TH2F("h1d2", "(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 2)
    {
      h1d2 = new TH2F("h1d2", "(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>3 && beam_p4_kin.E()<8.18 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 3)
    {
      h1d2 = new TH2F("h1d2", "(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 4)
    {
      h1d2 = new TH2F("h1d2", "(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>8.18 && beam_p4_kin.E()<9.15 && -t_kin>1.13 && -t_kin<10)");
    }
    if (j == 5)
    {
      h1d2 = new TH2F("h1d2", "(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>0 && -t_kin<1.13)");
    }
    if (j == 6)
    {
      h1d2 = new TH2F("h1d2", "(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}} [GeV/c^{2}]", n2pi2k, m2pi2k_min, m2pi2k_max, nkk, mkk_min, mkk_max);
      tdata->Project("h1d2", "kpkm_mf:kpkmpippim_mf", "w8*((kpkm_uni || kpkmpippim_uni) && beam_p4_kin.E()>9.15 && beam_p4_kin.E()<12 && -t_kin>1.13 && -t_kin<10)");
    }

    h1d2->Draw("colz");

    cphiy1[j] = new TCanvas(Form("cphiy1_%d", j), Form("cphiy1_%d", j), 1500, 800);
    cphiy1[j]->Divide(5, 5);
    cphiy2[j] = new TCanvas(Form("cphiy2_%d", j), Form("cphiy2_%d", j), 1500, 800);
    cphiy2[j]->Divide(5, 5);

    gphiy[j] = new TGraphErrors();//n2pi2k
    gphiy[j]->SetMarkerStyle(20);
    gphiy[j]->SetMarkerSize(1.0);
    gphiy[j]->SetMarkerColor(1);
    // gphiy[j]->SetTitle(Form("%.2f<-t<%.2f (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}", t1[j], t2[j]));
    if(j==1) gphiy[j]->SetTitle("(3,0)<(E_{#gamma},-t)<(8.18,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==2) gphiy[j]->SetTitle("(3,1.13)<(E_{#gamma},-t)<(8.18,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==3) gphiy[j]->SetTitle("(8.18,0)<(E_{#gamma},-t)<(9.15,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==4) gphiy[j]->SetTitle("(8.18,1.13)<(E_{#gamma},-t)<(9.15,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==5) gphiy[j]->SetTitle("(9.15,0)<(E_{#gamma},-t)<(12,1.13) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    if(j==6) gphiy[j]->SetTitle("(9.15,1.13)<(E_{#gamma},-t)<(12,10) (Data);m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");

    for (int i = 1; i <= n2pi2k; ++i)
    {
      // cout << i << " " << flush;

      if(i < 26) cphiy1[j]->cd(i);
      if(i > 25 && i<51) cphiy2[j]->cd(i-25);

      TH1D *hphiy_py = h1d2->ProjectionY(Form("_hphiy_py_%d_%d", j, i), i, i);

      // hslice2->Fit("fsb","q","",0.99,1.08);
      // hslice2->Fit("fsb","qm","",0.99,1.12);
      // fs->SetParameters(fsb->GetParameters());
      // fb2->SetParameters(fsb->GetParameters());
      hphiy_py->Draw();
      if(hphiy_py->Integral(1,50)<100)
        continue;
      // fb2->Draw("same");

      // TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + pol2(3)", mkk_min, mkk_max);
      TF1 *fsb = new TF1("fsb", "[0]*TMath::Voigt(x - [1], [2], [3]) + pol2(4)", mkk_min, mkk_max);
      fsb->SetLineColor(2);
      fsb->SetParameters(1433, 1.019, 0.004, 0.002, 1, 1, 1);
      fsb->SetParLimits(0, 0, 100000);
      fsb->SetParLimits(1, 1.015,1.022);// 1.018, 1.021
      fsb->FixParameter(2, 0.004); //fsb->SetParLimits(2, 0.008, 0.010);
      fsb->FixParameter(3, 0.002); //fsb->SetParLimits(3, 0.001,0.01);// 0.001,0.01 
      
      // TF1 *fs = new TF1("fs", "[0]*TMath::BreitWigner(x,[1],[2])", mkk_min, mkk_max);
      TF1 *fs = new TF1("fs", "[0]*TMath::Voigt(x - [1], [2], [3])", mkk_min, mkk_max);
      fs->SetLineColor(4);

      TF1 *fb = new TF1("fb", "pol2(4)", mkk_min, mkk_max);//4
      fb->SetLineColor(28);
      fb->SetLineStyle(2);

      hphiy_py->Fit("fsb", "", "", mkk_min, mkk_max);
      double par[7];//6
      fsb->GetParameters(&par[0]);
      fs->SetParameters(&par[0]);
      fb->SetParameters(&par[4]);//3

      fs->Draw("same");
      fb->Draw("same");

      double N1 = fs->Integral(mkk_min, mkk_max) / hphiy_py->GetBinWidth(1);
      double dN1 = N1 * fsb->GetParError(0) / fsb->GetParameter(0);

      if (N1 <= 0)
        gphiy[j]->RemovePoint(i-1);

      gphiy[j]->SetPoint(i - 1, h1d2->GetXaxis()->GetBinCenter(i), N1);
      gphiy[j]->SetPointError(i - 1, 0, dN1);

      TLatex lat_phiy;
      lat_phiy.SetTextSize(0.09);
      lat_phiy.SetTextAlign(13); //align at top
      lat_phiy.SetNDC();
      lat_phiy.SetTextColor(kBlue);
      lat_phiy.DrawLatex(0.45, 0.88, Form("#chi^{2}/NDF = %0.2f", fsb->GetChisquare() / fsb->GetNDF()));
      lat_phiy.DrawLatex(0.45, 0.78, Form("N_{sig} = %0.2f#pm%0.2f", N1, dN1));
      lat_phiy.DrawLatex(0.45, 0.68, Form("#mu = %0.3f#pm%0.3f", fsb->GetParameter(1), fsb->GetParError(1)));
      lat_phiy.DrawLatex(0.45, 0.58, Form("#Gamma = %0.3f#pm%0.3f", fsb->GetParameter(2), fsb->GetParError(2)));
      lat_phiy.DrawLatex(0.45, 0.48, Form("#sigma = %0.3f#pm%0.3f", fsb->GetParameter(3), fsb->GetParError(3)));
      
    }

    // cout << endl;

    cgphiy->cd(j);
    gphiy[j]->Draw("AP");

    // int j =1;
    // gphiy->Write(Form("grphiy_%d", j), TObject::kWriteDelete);

    cphiy1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiy_%d.root", j), "root");
    cphiy1[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c1_phiy_%d.eps", j), "eps");
    cphiy2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiy_%d.root", j), "root");
    cphiy2[j]->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c2_phiy_%d.eps", j), "eps");  
  }

  cphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiy.root", "root");
  cphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_phiy.eps", "eps");
  cgphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiy.root", "root");
  cgphiy->Print("/data.local/nacer/halld_my/pippimkpkm/fig_scany/c_gphiy.eps", "eps");

  outputfig->Print();

}
