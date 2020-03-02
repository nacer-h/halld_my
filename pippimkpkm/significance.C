//  File: pimpipkmkp.C # Author: Nacer # Date: 29 January 2017 # Email: a.hamdi@gsi.de # Description: Macro to study Y(2175).

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

void significance(TString cut, TString name)
{
  // TFile *fsim = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_yphifo_%s_chi100_flat.root", name.Data()));
  // TTree *tsim=(TTree*)fsim->Get("ntp");
  TFile *fdata = new TFile(Form("~/lochebe/ahamdi/gluex_root_analysis/workdir/dataout/data/tree_pippimkpkm_%s_chi100_flat.root", name.Data()));
	TTree *tdata=(TTree*)fdata->Get("ntp");
  TFile *outputfig = new TFile("/data.local/nacer/halld_my/pippimkpkm/fig_significance/significance.root","UPDATE");

  RooWorkspace w("w", kTRUE);
  //RooFitResult* result = NULL;

  // ***************** cuts
  // TString cutlist="&& kpkm_mf>1.005 && kpkm_mf<1.035 && kin_chisq<30 && abs(mm2)<0.015 && (abs(kp_dttof)<0.18 || kp_dttof==-999) && p_x4_kin.Z()<78 && p_x4_kin.Z()>52 && beam_p4_kin.E()>6 &&";
  //  && (abs(kp_dtbcal)<0.75 || kp_dtbcal==-999) && (abs(kp_dtfcal)<0.5 || kp_dtfcal==-999) && (abs(pip_dttof)<0.5 || pip_dttof==-999) && (abs(pip_dtbcal)<1.0 || pip_dtbcal==-999) && (abs(pip_dtfcal)<1.1 || pip_dtfcal==-999) && (abs(p_dttof)<0.6 || p_dttof==-999) && (abs(p_dtbcal)<1.0 || p_dtbcal==-999)
  TString cutlist = "&& kpkm_mf>1.005 && kpkm_mf<1.035 && ((kp_dttof>-0.22 && kp_dttof<0.22) || kp_dttof == -999) && ((kp_dtbcal>-0.65 && kp_dtbcal<0.3) || kp_dtbcal == -999) && ((kp_dtfcal>-2.0 && kp_dtfcal<0.8) || kp_dtfcal == -999) && (abs(pip_dttof)<0.3 || pip_dttof == -999) && (abs(pip_dtbcal)<0.6 || pip_dtbcal == -999) && ((p_dttof>-0.25 && p_dttof<0.25) || p_dttof == -999) && abs(mm2)<0.03 &&";// abs(mm2)<0.03 &&
  
  //  (abs(kp_dttof)<0.135 || kp_dttof == -999) && (abs(kp_dtbcal)<0.4125 || kp_dtbcal == -999) && (abs(kp_dtfcal)<1 || kp_dtfcal == -999) && (abs(pip_dttof)<0.15 || pip_dttof == -999) && (abs(pip_dtbcal)<0.85 || pip_dtbcal == -999) && (abs(p_dttof)<0.27 || p_dttof == -999) && abs(mm2)<0.01 &&>7.95 &&

  double cutsmin = 0; //= hdata_postcut->GetXaxis()->GetBinLowEdge(1);
  double cutsmax = 100;//= hdata_postcut->GetXaxis()->GetBinUpEdge(600);
	double cutsstep = (cutsmax-cutsmin) / 20.;
  double cuts[20];
  for (int i = 0; i < 20; ++i)
  {
    cuts[i] = cutsmax - (i * cutsstep); //i * cutsstep;
    // cuts[i] = cutsmin + (i * cutsstep); //i * cutsstep;
    cout<<"########  i = "<<i<<" | cuts[i] = "<<cuts[i]<<endl;
  }

/*
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Y(2175) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // *********************************************** Data  ***********************************************
  gStyle->SetOptFit(0);

  TCanvas *cdata_YMass = new TCanvas("cdata_YMass", "cdata_YMass", 1500, 800);
  cdata_YMass->Divide(5, 4);
  // cdata_YMass->Divide(1, 1);
  TH1F *hdata_YMass[20];

  double nbkg_y[20], nbkgerr_y[20], test[20];

  for (int j = 0; j < 20; j++)
  {
    cdata_YMass->cd(j+1);
    // hdata_YMass = (TH1F *)fdata->Get(Form("h_YMass_cuts_%d", j));
    hdata_YMass[j] = new TH1F(Form("hdata_YMass_%d",j), ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 1.6, 3.2);//2.5, 3.2
    tdata->Project(Form("hdata_YMass_%d",j), "kpkmpippim_mf", Form("w8*(kpkmpippim_uni"+cutlist+"kin_chisq<%f)",cuts[j])); //cuts[j] -t_kin<1 && beam_p4_kin.E()>6 && pt<0.4
    cout << "hdata_YMass = " << hdata_YMass[j] << endl;
    // hdata_YMass->Rebin(4);
    // double xmin = hdata_YMass->GetXaxis()->GetXmin();
    // double xmax = hdata_YMass->GetXaxis()->GetXmax();
    hdata_YMass[j]->SetMarkerStyle(20);
    hdata_YMass[j]->SetMarkerSize(0.7);
    hdata_YMass[j]->SetTitle(Form("mm^{2} < %.2f",cuts[j]));//#chi^{2}

    TF1 *fbkg_y = new TF1("fbkg_y","pol2");
    fbkg_y->SetParameters(1.,1.);
    hdata_YMass[j]->Fit(fbkg_y,"","mR",2.5, 3.2);
    hdata_YMass[j]->Draw("e");

    // nbkg_y[j] = hdata_YMass[j]->Integral(1, 100);
    nbkg_y[j] = fbkg_y->Integral(2.0, 2.5)/hdata_YMass[j]->GetBinWidth(1);
    nbkgerr_y[j] = sqrt(nbkg_y[j]);

    test[j] = hdata_YMass[j]->Integral(hdata_YMass[j]->FindBin(2.0),hdata_YMass[j]->FindBin(2.5))/hdata_YMass[j]->Integral(hdata_YMass[j]->FindBin(2.5),hdata_YMass[j]->FindBin(3.2)); //  / hdata_YMass[j]->Integral(hdata_YMass[j]->FindBin(2.5),hdata_YMas[j]s->FindBin(3.2)); //fbkg_y->Integral(2.0, 2.5) / fbkg_y->Integral(2.5, 3.2);

    hdata_YMass[j]->Write();
  }

  cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%s_%scut.root", name.Data(), cut.Data()), "root");
  cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%s_%scut.eps", name.Data(), cut.Data()), "eps");

 // *********************************************** Sim  ***********************************************

  TCanvas *csim_YMass = new TCanvas("csim_YMass", "csim_YMass", 1600, 1000);
  csim_YMass->Divide(5, 4);
  // csim_YMass->Divide(1,1);
  TH1F *hsim_YMass[20];

  TCanvas *cgrsim_YMass = new TCanvas("cgrsim_YMass","cgrsim_YMass",900,600);
  cgrsim_YMass->SetGrid();
  TGraphErrors *grsim_YMass = new TGraphErrors();
  grsim_YMass->SetMinimum(0.);
  grsim_YMass->SetMarkerStyle(20);
  grsim_YMass->GetXaxis()->SetNdivisions(510, kFALSE);

  ofstream ofs_YMass("tableau_YMass.txt", ofstream::out);

  double nsig_y[20], nsigerr_y[20];

  for (int j = 0; j < 20; j++)
  {
    csim_YMass->cd(j + 1);
    hsim_YMass[j] = new TH1F(Form("hsim_YMass_%d",j), ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 100, 2.0, 2.5); //2.05, 2.35
    tsim->Project(Form("hsim_YMass_%d",j), "kpkmpippim_mf", Form("w8*(kpkmpippim_uni"+cutlist+"kin_chisq<%f)",cuts[j]));
    // hsim_YMass->Rebin(4);
    hsim_YMass[j]->SetMarkerStyle(20);
    hsim_YMass[j]->SetMarkerSize(0.7);
    cout << "hsim_YMass = " << hsim_YMass[j] << endl;
    hsim_YMass[j]->SetTitle(Form("mm^{2} < %.2f",cuts[j]));//#chi^{2}
    hsim_YMass[j]->Draw("e");

    nsig_y[j] = hsim_YMass[j]->Integral(1, 100);
    nsigerr_y[j] = sqrt(nsig_y[j]);
    double z_y = nsig_y[j] / sqrt(nbkg_y[j]);
    double zerr_y = sqrt(((4 * nbkg_y[j] * nbkg_y[j] * nsigerr_y[j] * nsigerr_y[j]) + (nsig_y[j] * nsig_y[j] * nbkgerr_y[j] * nbkgerr_y[j])) / (4 * nbkg_y[j] * nbkg_y[j] * nbkg_y[j]));

    ofs_YMass << " j = " << j << " | cuts = " << cuts[j] << " | nsig_y = " << nsig_y[j] << " | nsigerr_y = " << nsigerr_y[j] << " | nbkg_y = " << nbkg_y[j] << " | nbkgerr_y = " << nbkgerr_y[j] << " | z_y = " << z_y << " | zerr_y = " << zerr_y <<" | test = " << test[j] << endl;

    grsim_YMass->SetPoint(j, cuts[j], z_y);
    grsim_YMass->SetPointError(j, 0, zerr_y);
    hsim_YMass[j]->Write();
  }

  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%s_%scut.root",name.Data(), cut.Data()), "root");
  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%s_%scut.eps",name.Data(), cut.Data()), "eps");

  cgrsim_YMass->cd();
  grsim_YMass->Draw("AP");

  // // result->Print();
  grsim_YMass->Write();
  grsim_YMass->SetTitle("Significance Y(2175) Vs selection (MC); cut; relative significance");
  cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%s_%scut.root",name.Data(),cut.Data()),"root");
  cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%s_%scut.eps",name.Data(),cut.Data()),"eps");

  ofs_YMass.close();
  outputfig->Print();

*/
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Phi(1020) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  // ************** Data ******************

  TCanvas *cdata_PhiMass = new TCanvas("cdata_PhiMass", "cdata_PhiMass", 1600, 1000);
  cdata_PhiMass->Divide(5, 4);
  // cdata_PhiMass->Divide(1, 1);
  TH1F *hdata_PhiMass[20];// = new TH1F("hdata_PhiMass", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 1.005, 1.035);

  TCanvas *cgrdata_PhiMass = new TCanvas("cgrdata_PhiMass","cgrdata_PhiMass",900,600);
  cgrdata_PhiMass->SetGrid();
  TGraphErrors *grdata_PhiMass = new TGraphErrors();
  grdata_PhiMass->SetMarkerStyle(20);
  grdata_PhiMass->GetXaxis()->SetNdivisions(510, kFALSE);

  ofstream ofs_PhiMass("tableau_PhiMass.txt", ofstream::out);

//   int k = 0;
// for (int i = 0; i < 5; i++)
//   {
   for (int j = 0; j < 20; j++)//j < 20
   {
    
    // cdata_PhiMass->cd(k + 1);
    cdata_PhiMass->cd(j + 1);
    // hdata_PhiMass = (TH1F *)fdata->Get(Form("h_PhiMass_cuts_%d", j));
    hdata_PhiMass[j] = new TH1F(Form("hdata_PhiMass_%d",j), ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 100, 1.005, 1.035);
    tdata->Project(Form("hdata_PhiMass_%d",j), "kpkm_mf", Form("w8*(kpkm_uni"+cutlist+"kin_chisq<%f)",cuts[j]));//kin_chisq
    //(pip_dttof<-990 || abs(pip_dttof)<0.5) && (pip_dtbcal<-990 || abs(pip_dtbcal)<1.0) && (pip_dtfcal<-990 || abs(pip_dtfcal)<0.5) && (kp_dttof<-990 || abs(kp_dttof)<0.3) && (kp_dtbcal<-990 || abs(kp_dtbcal)<0.6) && (kp_dtfcal<-990 || abs(kp_dtfcal)<2.375) && (p_dttof<-990 || abs(p_dttof)<0.48) && (p_dtbcal<-990 || abs(p_dtbcal)<%f) && (p_dtfcal<-990 || abs(p_dtfcal)<1.9))",cuts[j]));
    // hdata_PhiMass->Rebin(4);
    hdata_PhiMass[j]->SetMarkerStyle(20);
    hdata_PhiMass[j]->SetMarkerSize(0.7);
    cout << "hdata_PhiMass = " << hdata_PhiMass[j] << endl;
    hdata_PhiMass[j]->SetTitle(Form("#chi^{2} < %0.f",cuts[j])); //|#Delta T| #chi^{2}
    hdata_PhiMass[j]->Draw("e");

    w.factory("Voigtian::sig_PhiMass_data(m_PhiMass_data[1.005,1.035],mean_PhiMass_data[1.017,1.021],width_PhiMass_data[0.004],sigma_PhiMass_data[0.001,0.01])");//sigma_PhiMass_data[0.001,0.01], mean_PhiMass_data[1.015,1.022]
    // w.factory("BreitWigner::sig(m[0.99,1.05],mean[1.017,1.021],width[0.0001,0.01])");
    // w.factory("ArgusBG::bkg_Phi_data(m_Phi_data, 1.04, c0_Phi_data[-50,-10])");
    w.factory("Chebychev::bkg_PhiMass_data(m_PhiMass_data,{c0_PhiMass_data[-10,10], c1_PhiMass_data[-10,10]})");
    w.factory("SUM:::model_PhiMass_data(nsig_PhiMass_data[0,100000000]*sig_PhiMass_data, nbkg_PhiMass_data[0,100000000]*bkg_PhiMass_data)"); //nsig[0,100000000]*sig2,
    w.var("m_PhiMass_data")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    RooDataHist dh_PhiMass_data("dh_PhiMass_data", "dh_PhiMass_data", *w.var("m_PhiMass_data"), Import(*hdata_PhiMass[j]));
    RooPlot *fr_PhiMass_data = w.var("m_PhiMass_data")->frame(Title("K^{+}K^{-}"));
    // fr_PhiMass_data->SetTitleOffset(0.90, "X");
    // fr_PhiMass_data->SetTitleSize(0.06, "XYZ");
    // fr_PhiMass_data->SetLabelSize(0.06, "xyz");
    w.pdf("model_PhiMass_data")->fitTo(dh_PhiMass_data);
    //cout<<"=============================  no problem up to here ! ========================"<<endl;

    // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    dh_PhiMass_data.plotOn(fr_PhiMass_data, RooFit::Name("ldh_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, Components(*w.pdf("sig_PhiMass_data")), LineColor(kRed), RooFit::Name("lsig_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, Components(*w.pdf("bkg_PhiMass_data")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, RooFit::Name("lmodel_PhiMass_data"));
    // w.pdf("model_PhiMass_data")->paramOn(fr_PhiMass_data, Layout(0.5, 0.90, 0.99) ,Parameters(RooArgSet(*w.var("nsig_PhiMass_data"),*w.var("nbkg_PhiMass_data"))));//,*w.var("mean_PhiMass_data"),*w.var("width_PhiMass_data"),*w.var("sigma_PhiMass_data"))));
    fr_PhiMass_data->Draw();

    TLegend *l_PhiMass_data = new TLegend(0.15, 0.6, 0.35, 0.85);
    l_PhiMass_data->SetTextSize(0.06);
    l_PhiMass_data->SetFillColor(kWhite);
    l_PhiMass_data->SetLineColor(kWhite);
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("ldh_PhiMass_data"), "Data", "p");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lmodel_PhiMass_data"), "Model", "l");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lbkg_PhiMass_data"), "Bkg", "l");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lsig_PhiMass_data"), "Signal", "l");
    l_PhiMass_data->Draw();

    double nsig_phi = w.var("nsig_PhiMass_data")->getVal();
    double nsigerr_phi = w.var("nsig_PhiMass_data")->getError();
    double nbkg_phi = w.var("nbkg_PhiMass_data")->getVal();
    double nbkgerr_phi = w.var("nbkg_PhiMass_data")->getError();
    double z_phi = nsig_phi/sqrt(nsig_phi+nbkg_phi);
    // double zerr_phi = sqrt(((4*nbkg_phi*nbkg_phi*nsigerr_phi*nsigerr_phi)+(nsig_phi*nsig_phi*nbkgerr_phi*nbkgerr_phi))/(4*nbkg_phi*nbkg_phi*nbkg_phi));
    double zerr_phi = sqrt((((nsig_phi+2*nbkg_phi)*(nsig_phi+2*nbkg_phi)*nsigerr_phi*nsigerr_phi)+(nsig_phi*nsig_phi*nbkgerr_phi*nbkgerr_phi))/(4*(nsig_phi+nbkg_phi)*(nsig_phi+nbkg_phi)*(nsig_phi+nbkg_phi)));

    ofs_PhiMass << " j = " << j << " | cuts = " << cuts[j] << " | nsig_phi = " << nsig_phi << " | nsigerr_phi = " << nsigerr_phi << " | nbkg_phi = " << nbkg_phi << " | nbkgerr_phi = " << nbkgerr_phi << " | z_phi = " << z_phi << " | zerr_phi = " << zerr_phi << endl;

    // // hdata_PhiMass->Draw("e");

    grdata_PhiMass->SetPoint(j,cuts[j],z_phi);
    grdata_PhiMass->SetPointError(j,0,zerr_phi);
    // grdata_PhiMass->SetPoint(k,k,z_phi);
    // grdata_PhiMass->SetPointError(k,0,zerr_phi);
    // k++;
    hdata_PhiMass[j]->Write();

    TLatex lat_PhiMass_data;
    lat_PhiMass_data.SetTextSize(0.06);
    lat_PhiMass_data.SetTextAlign(13); //align at top
    lat_PhiMass_data.SetNDC();
    lat_PhiMass_data.SetTextColor(kBlue);
    lat_PhiMass_data.DrawLatex(0.58, 0.85, Form("N_{Sig} = %0.f#pm%0.f", nsig_phi, nsigerr_phi));
    lat_PhiMass_data.DrawLatex(0.58, 0.78, Form("N_{Bkg} = %0.f#pm%0.f", nbkg_phi, nbkgerr_phi));    
  }

  //result->Print();
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%s_%scut.root",name.Data(), cut.Data()), "root");
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%s_%scut.eps",name.Data(), cut.Data()), "eps");

  cgrdata_PhiMass->cd();
  grdata_PhiMass->Draw("AP");

  grdata_PhiMass->Write();
  grdata_PhiMass->SetTitle("#phi(1020) Significance Vs selection (data);cut;S/#sqrt{S+B}");
  cgrdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_PhiMass_%s_%scut.root",name.Data(), cut.Data()),"root");
  cgrdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_PhiMass_%s_%scut.eps",name.Data(), cut.Data()),"eps");

  ofs_PhiMass.close(); 
  outputfig->Print();
  // outputfig->Close();
}

/*
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ********************************* Phi(1020) *********************************

  TCanvas *cdata_PhiMass = new TCanvas("cdata_PhiMass", "cdata_PhiMass", 1500, 800);
  cdata_PhiMass->Divide(5, 4);
  // cdata_PhiMass->Divide(1, 1);
  TH1F *hdata_PhiMass;// = new TH1F("hdata_PhiMass", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 1.005, 1.035);

  TCanvas *cgrdata_PhiMass = new TCanvas("cgrdata_PhiMass","cgrdata_PhiMass",600,400);
  cgrdata_PhiMass->SetGrid();
  TGraphErrors *grdata_PhiMass = new TGraphErrors();
  grdata_PhiMass->SetMarkerStyle(20);
  grdata_PhiMass->GetXaxis()->SetNdivisions(510, kFALSE);

  ofstream ofs_PhiMass("tableau_PhiMass.txt", ofstream::out);

//   int k = 0;
// for (int i = 0; i < 5; i++)
//   {
   for (int j = 0; j < 20; j++)
   {
    
    // cdata_PhiMass->cd(k + 1);
    cdata_PhiMass->cd(j + 1);
    // hdata_PhiMass = (TH1F *)fdata->Get(Form("h_PhiMass_cuts_%d", j));
    hdata_PhiMass = new TH1F("hdata_PhiMass", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 50, 1.005, 1.035);
    tdata->Project("hdata_PhiMass", "kpkm_mf", Form("w8*(kpkm_uni"+cutlist+"kin_chisq<%f)", cuts[j]));
    //(pip_dttof<-990 || abs(pip_dttof)<0.5) && (pip_dtbcal<-990 || abs(pip_dtbcal)<1.0) && (pip_dtfcal<-990 || abs(pip_dtfcal)<0.5) && (kp_dttof<-990 || abs(kp_dttof)<0.3) && (kp_dtbcal<-990 || abs(kp_dtbcal)<0.6) && (kp_dtfcal<-990 || abs(kp_dtfcal)<2.375) && (p_dttof<-990 || abs(p_dttof)<0.48) && (p_dtbcal<-990 || abs(p_dtbcal)<%f) && (p_dtfcal<-990 || abs(p_dtfcal)<1.9))",cuts[j]));
    // hdata_PhiMass->Rebin(4);
    hdata_PhiMass->SetMarkerStyle(20);
    hdata_PhiMass->SetMarkerSize(0.7);
    cout << "hdata_PhiMass = " << hdata_PhiMass << endl;
    hdata_PhiMass->SetTitle(Form("|#chi^2| < %.2f",cuts[j]));//|#Delta T|
    hdata_PhiMass->Draw("e");

    w.factory("Voigtian::sig_PhiMass_data(m_PhiMass_data[1.005,1.035],mean_PhiMass_data[1.015,1.022],width_PhiMass_data[0.004],sigma_PhiMass_data[0.001,0.01])");//sigma_PhiMass_data[0.001,0.01], mean_PhiMass_data[1.016,1.022]
    // w.factory("BreitWigner::sig(m[0.99,1.05],mean[1.017,1.021],width[0.0001,0.01])");
    // w.factory("ArgusBG::bkg_Phi_data(m_Phi_data, 1.04, c0_Phi_data[-50,-10])");
    w.factory("Chebychev::bkg_PhiMass_data(m_PhiMass_data,{c0_PhiMass_data[-10,10], c1_PhiMass_data[-10,10]})");
    w.factory("SUM:::model_PhiMass_data(nsig_PhiMass_data[0,100000000]*sig_PhiMass_data, nbkg_PhiMass_data[0,100000000]*bkg_PhiMass_data)"); //nsig[0,100000000]*sig2,
    w.var("m_PhiMass_data")->SetTitle("m_{K^{+}K^{-}} [GeV/c^{2}]");
    RooDataHist dh_PhiMass_data("dh_PhiMass_data", "dh_PhiMass_data", *w.var("m_PhiMass_data"), Import(*hdata_PhiMass));
    RooPlot *fr_PhiMass_data = w.var("m_PhiMass_data")->frame(Title("K^{+}K^{-}"));
    // fr_PhiMass_data->SetTitleOffset(0.90, "X");
    // fr_PhiMass_data->SetTitleSize(0.06, "XYZ");
    // fr_PhiMass_data->SetLabelSize(0.06, "xyz");
    w.pdf("model_PhiMass_data")->fitTo(dh_PhiMass_data);
    //cout<<"=============================  no problem up to here ! ========================"<<endl;

    // //result = w.pdf("model")->fitTo(dh_PhiMass,Extended(kTRUE),Save());
    dh_PhiMass_data.plotOn(fr_PhiMass_data, RooFit::Name("ldh_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, Components(*w.pdf("sig_PhiMass_data")), LineColor(kRed), RooFit::Name("lsig_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, Components(*w.pdf("bkg_PhiMass_data")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_PhiMass_data"));
    w.pdf("model_PhiMass_data")->plotOn(fr_PhiMass_data, RooFit::Name("lmodel_PhiMass_data"));
    w.pdf("model_PhiMass_data")->paramOn(fr_PhiMass_data, Layout(0.5, 0.90, 0.99) ,Parameters(RooArgSet(*w.var("nsig_PhiMass_data"),*w.var("nbkg_PhiMass_data"))));//,*w.var("mean_PhiMass_data"),*w.var("width_PhiMass_data"),*w.var("sigma_PhiMass_data"))));
    fr_PhiMass_data->Draw();

    TLegend *l_PhiMass_data = new TLegend(0.13, 0.6, 0.3, 0.85);
    l_PhiMass_data->SetFillColor(kWhite);
    l_PhiMass_data->SetLineColor(kWhite);
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("ldh_PhiMass_data"), "Data", "p");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lmodel_PhiMass_data"), "total", "l");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lbkg_PhiMass_data"), "Background", "l");
    l_PhiMass_data->AddEntry(fr_PhiMass_data->findObject("lsig_PhiMass_data"), "Signal", "l");
    l_PhiMass_data->Draw();

    double nsig_phi = w.var("nsig_PhiMass_data")->getVal();
    double nsigerr_phi = w.var("nsig_PhiMass_data")->getError();
    double nbkg_phi = w.var("nbkg_PhiMass_data")->getVal();
    double nbkgerr_phi = w.var("nbkg_PhiMass_data")->getError();
    double z_phi = nsig_phi/sqrt(nsig_phi+nbkg_phi);
    // double zerr_phi = sqrt(((4*nbkg_phi*nbkg_phi*nsigerr_phi*nsigerr_phi)+(nsig_phi*nsig_phi*nbkgerr_phi*nbkgerr_phi))/(4*nbkg_phi*nbkg_phi*nbkg_phi));
    double zerr_phi = sqrt((((nsig_phi+2*nbkg_phi)*(nsig_phi+2*nbkg_phi)*nsigerr_phi*nsigerr_phi)+(nsig_phi*nsig_phi*nbkgerr_phi*nbkgerr_phi))/(4*(nsig_phi+nbkg_phi)*(nsig_phi+nbkg_phi)*(nsig_phi+nbkg_phi)));

    ofs_PhiMass << " j = " << j << " | cuts = " << cuts[j] << " | nsig_phi = " << nsig_phi << " | nsigerr_phi = " << nsigerr_phi << " | nbkg_phi = " << nbkg_phi << " | nbkgerr_phi = " << nbkgerr_phi << " | z_phi = " << z_phi << " | zerr_phi = " << zerr_phi << endl;

    // // hdata_PhiMass->Draw("e");

    grdata_PhiMass->SetPoint(j,cuts[j],z_phi);
    grdata_PhiMass->SetPointError(j,0,zerr_phi);
    // grdata_PhiMass->SetPoint(k,k,z_phi);
    // grdata_PhiMass->SetPointError(k,0,zerr_phi);
    // k++;
    hdata_PhiMass->Write();
  }

  //result->Print();
  // h_PhiMass->SetMarkerStyle(20);
  // h_PhiMass->SetMarkerSize(0.7);
  // h_PhiMass->Draw("e");
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%scut.root", cut.Data()), "root");
  cdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_PhiMass_%scut.eps", cut.Data()), "eps");

  cgrdata_PhiMass->cd();
  grdata_PhiMass->Draw("AP");


  // TF1 *fdatapol6 = new TF1("fdatapol6","pol6");
  // grdata->Fit("fdatapol6","mR","", cuts[12], cuts[4]);
  // double zd_max = fdatapol6->GetMaximumX(cuts[12], cuts[4]);
  // cout<<"zd_max = "<<zd_max<<endl;

  grdata_PhiMass->Write();
  grdata_PhiMass->SetTitle("Significance #phi(1020) Vs selection (data); cut;  S/#sqrt{S+B}");
  cgrdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_PhiMass_%scut.root",cut.Data()),"root");
  cgrdata_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrdata_PhiMass_%scut.eps",cut.Data()),"eps");

  // // ********************************* Y(2175) *********************************

  // TCanvas *cdata_YMass = new TCanvas("cdata_YMass", "cdata_YMass", 1500, 800);
  // cdata_YMass->Divide(5, 4);
  // // cdata_YMass->Divide(1, 1);
  // TH1F *hdata_YMass;

  // double nbkg_y[20], nbkgerr_y[20];

  // for (int j = 0; j < 20; j++)
  // {
  //   cdata_YMass->cd(j + 1);
  //   // hdata_YMass = (TH1F *)fdata->Get(Form("h_YMass_cuts_%d", j));
  //   hdata_YMass = new TH1F("hdata_YMass", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 2.5, 3.2);
  //   tdata->Project("hdata_YMass", "kpkmpippim_mf", Form("w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && kin_chisq<60 && abs(mm2)<0.03 && sqrt((p_x4_kin.X()*p_x4_kin.X())+(p_x4_kin.Y()*p_x4_kin.Y()))<%f)",cuts[j]));
  //   cout << "hdata_YMass = " << hdata_YMass << endl;
  //   hdata_YMass->Rebin(4);
  //   // double xmin = hdata_YMass->GetXaxis()->GetXmin();
  //   // double xmax = hdata_YMass->GetXaxis()->GetXmax();
  //   hdata_YMass->SetMarkerStyle(20);
  //   hdata_YMass->SetMarkerSize(0.7);
  //   hdata_YMass->SetTitle(Form("#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} < %.2f",cuts[j]));
  //   hdata_YMass->Draw("e");

  //   nbkg_y[j] = hdata_YMass->Integral(hdata_YMass->GetXaxis()->GetBinLowEdge(1), hdata_YMass->GetXaxis()->GetBinUpEdge(200));
  //   nbkgerr_y[j] = sqrt(nbkg_y[j]);

  //   // // w.factory("Voigtian::sig_YMass_data(m_YMass_data[1.6,3.2],mean_YMass_data[1.7,2.0],width_YMass_data[0.0001,0.1],sigma_YMass_data[0.0001,0.1])");
  //   // w.factory("Voigtian::sig_YMass_data(m_YMass_data[1.6,3.2],mean_YMass_data[1.7,1.95],width_YMass_data[0.001,0.08],sigma_YMass_data[0.001,0.05])");
  //   // // w.factory("BreitWigner::sig_Y_data(m_Y_data[1.6,3.2],mean_Y_data[1.75,1.85],width_Y_data[0.01,1])");
  //   // // w.factory("ArgusBG::bkg_Y_data(m_Y_data, 3.4, c0_Y_data[-100,-1])");
  //   // w.factory("Chebychev::bkg_YMass_data(m_YMass_data,{c0_YMass_data[-10,10], c1_YMass_data[-10,10] ,c2_YMass_data[-10,10], c3_YMass_data[-10,10]})");
  //   // w.factory("SUM:::model_YMass_data(nsig_YMass_data[0,100000000]*sig_YMass_data, nbkg_YMass_data[0,100000000]*bkg_YMass_data)"); //nsig[0,100000000]*sig2,
  //   // w.var("m_YMass_data")->SetTitle("m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}]");
  //   // RooDataHist dh_YMass_data("dh_YMass_data", "dh_YMass_data", *w.var("m_YMass_data"), Import(*hdata_YMass));
  //   // RooPlot *fr_YMass_data = w.var("m_YMass_data")->frame(Title("invariant mass K^{+}K^{-}"));
  //   // fr_YMass_data->SetTitleOffset(0.90, "X");
  //   // fr_YMass_data->SetTitleSize(0.06, "XYZ");
  //   // fr_YMass_data->SetLabelSize(0.06, "xyz");
  //   // w.pdf("model_YMass_data")->fitTo(dh_YMass_data);
  //   // //cout<<"=============================  no problem up to here ! ========================"<<endl;

  //   // // //result = w.pdf("model")->fitTo(dh_YMass,Extended(kTRUE),Save());
  //   // dh_YMass_data.plotOn(fr_YMass_data, RooFit::Name("ldh_YMass_data"));
  //   // w.pdf("model_YMass_data")->plotOn(fr_YMass_data, Components(*w.pdf("sig_YMass_data")), LineColor(kRed), RooFit::Name("lsig_YMass_data"));
  //   // w.pdf("model_YMass_data")->plotOn(fr_YMass_data, Components(*w.pdf("bkg_YMass_data")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_YMass_data"));
  //   // w.pdf("model_YMass_data")->plotOn(fr_YMass_data, RooFit::Name("lmodel_YMass_data"));
  //   // w.pdf("model_YMass_data")->paramOn(fr_YMass_data, Layout(0.5, 0.90, 0.99),Parameters(RooArgSet(*w.var("nsig_YMass_data"),*w.var("nbkg_YMass_data"),*w.var("mean_YMass_data"),*w.var("width_YMass_data"))));//,*w.var("mean_YMass_data"),*w.var("width_YMass_data"),*w.var("sigma_YMass_data"))));
  //   // fr_YMass_data->Draw();

  //   // TLegend *l_YMass_data = new TLegend(0.35, 0.85, 0.49, 0.99);
  //   // l_YMass_data->SetFillColor(kWhite);
  //   // l_YMass_data->SetLineColor(kWhite);
  //   // l_YMass_data->AddEntry(fr_YMass_data->findObject("ldh_YMass_data"), "Data", "p");
  //   // l_YMass_data->AddEntry(fr_YMass_data->findObject("lmodel_YMass_data"), "total", "l");
  //   // l_YMass_data->AddEntry(fr_YMass_data->findObject("lbkg_YMass_data"), "Background", "l");
  //   // l_YMass_data->AddEntry(fr_YMass_data->findObject("lsig_YMass_data"), "Signal", "l");
  //   // l_YMass_data->Draw();

  //   // nbkg_y[j] = w.var("nbkg_YMass_data")->getVal();
  //   // nbkgerr_y[j] = w.var("nbkg_YMass_data")->getError();

  //   // double pardata_y[3];
  //   // TF1 *fsigdata_y = new TF1("fsigdata_y","gaus"); //[0]*TMath::BreitWigner(x,[1],[2])
  //   // fsigdata_y->SetLineStyle(9);
  //   // fsigdata_y->SetLineColor(kRed);
  //   // fsigdata_y->SetParameters(16000,2.175,0.04);
  //   // hdata_YMass->Fit(fsigdata_y,"","mR",1.8,2.6);
  //   // fsig_y->GetParameters(&par_y[0]);
  //   // TF1 *fbkg_y = (TF1*) gROOT->GetFunction("chebyshev4");

  //   // TF1 *fbkg_phiy = new TF1("fbkg_phiy","pol4");
  //   // fbkg_phiy->SetLineStyle(9);
  //   // fbkg_phiy->SetLineColor(4);
  //   // fbkg_phiy->SetParameters(1.,1.,1.,1.,1.);
  //   // grphiy->Fit(fbkg_phiy,"","mR",1.46,3.20);
  //   // fbkg_phiy->GetParameters(&paramphiy[3]);
  //   // TF1 *fmodel_phiy = new TF1("fmodel_phiy","[0]*TMath::BreitWigner(x,[1],[2])+pol4(3)");
  //   // fmodel_phiy->SetParameters(&paramphiy[0]);
  //   // fmodel_phiy->SetParLimits(1,1.60,1.80);
  //   // // fmodel_phifo->SetParLimits(4,);
  //   // grphiy->Fit(fmodel_phiy,"","mR",1.46,3.20);
  //   hdata_YMass->Write();
  }

  // cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%scut.root", cut.Data()), "root");
  // cdata_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_YMass_%scut.eps", cut.Data()), "eps");

  // // ********************************* fo(980) *********************************

  // TCanvas *cdata_foMass = new TCanvas("cdata_foMass", "cdata_foMass", 1500, 800);
  // cdata_foMass->Divide(5, 4);
  // // cdata_foMass->Divide(1,3);
  // TH1F *hdata_foMass;

  // for (int j = 0; j < 20; j++)
  // {
  //   cdata_foMass->cd(j + 1);
  //   hdata_foMass = (TH1F *)fdata->Get(Form("h_foMass_cuts_%d", j));
  //   hdata_foMass->Rebin(4);
  //   double xmin = hdata_foMass->GetXaxis()->GetXmin();
  //   double xmax = hdata_foMass->GetXaxis()->GetXmax();
  //   hdata_foMass->SetMarkerStyle(20);
  //   hdata_foMass->SetMarkerSize(0.7);
  //   cout << "hdata_foMass = " << hdata_foMass << endl;
  //   hdata_foMass->Draw("e");
  // }

  // cdata_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_foMass_%scut.root", cut.Data()), "root");
  // cdata_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cdata_foMass_%scut.eps", cut.Data()), "eps");

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    MC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ********************************* Phi(1020) *********************************

  TCanvas *csim_PhiMass = new TCanvas("csim_PhiMass", "csim_PhiMass", 1500, 800);
  csim_PhiMass->Divide(5, 4);
  // csim_PhiMass->Divide(1,1);
  TH1F *hsim_PhiMass;

  // TCanvas *cgrsim_PhiMass = new TCanvas("cgrsim_PhiMass","cgrsim_PhiMass",600,400);
  // cgrsim_PhiMass->SetGrid();
  // TGraphErrors *grsim_PhiMass = new TGraphErrors();
  // grsim_PhiMass->SetMarkerStyle(20);
  // grsim_PhiMass->GetXaxis()->SetNdivisions(510, kFALSE);

  for (int j = 0; j < 20; j++)
  {
    csim_PhiMass->cd(j + 1);
    hsim_PhiMass = (TH1F *)fsim->Get(Form("h_PhiMass_cuts_%d", j));
    hsim_PhiMass->Rebin(4);
    double xmin = hsim_PhiMass->GetXaxis()->GetXmin();
    double xmax = hsim_PhiMass->GetXaxis()->GetXmax();
    hsim_PhiMass->SetMarkerStyle(20);
    hsim_PhiMass->SetMarkerSize(0.7);
    cout << "hsim_PhiMass = " << hsim_PhiMass << endl;
    hsim_PhiMass->Draw("e");

    // w.factory("Voigtian::sig_PhiMass_mc(m_PhiMass_mc[1.005,1.035],mean_PhiMass_mc[1.017,1.021],width_PhiMass_mc[0.0004,0.04],sigma_PhiMass_mc[0.0004,0.04])");
    w.factory("Gaussian::sig1_PhiMass_mc(m_PhiMass_mc[1.005,1.035],mean1_PhiMass_mc[1.019,1.017,1.021],sigma1_PhiMass_mc[0.0001,0.02])");
    w.factory("Gaussian::sig2_PhiMass_mc(m_PhiMass_mc[1.005,1.035],mean2_PhiMass_mc[1.019,1.017,1.021],sigma2_PhiMass_mc[0.0001,0.02])");
    // w.factory("BreitWigner::sig1_PhiMass_mc(m_PhiMass_mc[0.99,1.05],mean1_PhiMass_mc[1.017,1.021],width1_PhiMass_mc[0.0001,0.1])");//mean1_PhiMass_mc[1.017,1.021],width1_PhiMass_mc[0.0001,0.1]
    // w.factory("Chebychev::bkg_PhiMass_mc(m_PhiMass_mc,{c0_PhiMass_mc[-10,10],c1_PhiMass_mc[-10,10],c2_PhiMass_mc[-10,10]})");
    w.factory("SUM:::model_PhiMass_mc(nsig_PhiMass_mc[0,100000000]*sig1_PhiMass_mc, nsig_PhiMass_mc[0,100000000]*sig2_PhiMass_mc)");//, nbkg_PhiMass_mc[0,100000000]*bkg_PhiMass_mc)"); // nsig_PhiMass_mc[0,100000000]*sig2_PhiMass_mc,
    w.var("m_PhiMass_mc")->SetTitle("m_PhiMass_{K^{+}K^{-}} [GeV/c^{2}]");
    RooDataHist dh_PhiMass_mc("dh_PhiMass_mc","dh_PhiMass_mc",*w.var("m_PhiMass_mc"),Import(*hsim_PhiMass));
    RooPlot* fr_PhiMass_mc = w.var("m_PhiMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
    fr_PhiMass_mc->SetTitleOffset(0.90,"X");
    fr_PhiMass_mc->SetTitleSize(0.06,"XYZ");
    fr_PhiMass_mc->SetLabelSize(0.06,"xyz");
    w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc);
    // result = w.pdf("model_PhiMass_mc")->fitTo(dh_PhiMass_mc,Extended(kTRUE),Save());
    dh_PhiMass_mc.plotOn(fr_PhiMass_mc,RooFit::Name("ldh_PhiMass_mc"));
    w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig1_PhiMass_mc")),LineColor(kRed),RooFit::Name("lsig1_PhiMass_mc"));
    w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("sig2_PhiMass_mc")),LineColor(kRed),RooFit::Name("lsig2_PhiMass_mc"));
    // w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc, Components(*w.pdf("bkg_PhiMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_PhiMass_mc"));
    w.pdf("model_PhiMass_mc")->plotOn(fr_PhiMass_mc,RooFit::Name("lmodel_PhiMass_mc"));
    w.pdf("model_PhiMass_mc")->paramOn(fr_PhiMass_mc,Layout(0.55,0.98,0.92) ,Parameters(RooArgSet(*w.var("nsig_PhiMass_mc"))));
    fr_PhiMass_mc->Draw();

    TLegend *l_PhiMass_mc = new TLegend(0.13,0.6,0.3,0.85);
    l_PhiMass_mc->SetFillColor(kWhite);
    l_PhiMass_mc->SetLineColor(kWhite);
    l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("ldh_PhiMass_mc"),"Data", "p");
    l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lmodel_PhiMass_mc"),"total","l");
    // l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lbkg_PhiMass_mc"),"Background", "l");
    l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lsig1_PhiMass_mc"),"gaus1", "l");
    l_PhiMass_mc->AddEntry(fr_PhiMass_mc->findObject("lsig2_PhiMass_mc"),"gaus2", "l");
    l_PhiMass_mc->Draw();

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
    hsim_PhiMass->Write();
  }

  csim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_PhiMass_%scut.root", cut.Data()), "root");
  csim_PhiMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_PhiMass_%scut.eps", cut.Data()), "eps");

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

  // // ********************************* fo(980) *********************************

  // TCanvas *csim_foMass = new TCanvas("csim_foMass", "csim_foMass", 1500, 800);
  // csim_foMass->Divide(5, 4);
  // //csim_foMass->Divide(1,3);
  // TH1F *hsim_foMass;

  // // TCanvas *cgrsim_foMass = new TCanvas("cgrsim_foMass","cgrsim_foMass",600,400);
  // // cgrsim_foMass->SetGrid();
  // // TGraphErrors *grsim_foMass = new TGraphErrors();
  // // grsim_foMass->SetMarkerStyle(20);
  // // grsim_foMass->GetXaxis()->SetNdivisions(510, kFALSE);

  // for (int j = 0; j < 20; j++)
  // {
  //   csim_foMass->cd(j + 1);
  //   hsim_foMass = (TH1F *)fsim->Get(Form("h_foMass_cuts_%d", j));
  //   hsim_foMass->Rebin(4);
  //   double xmin = hsim_foMass->GetXaxis()->GetXmin();
  //   double xmax = hsim_foMass->GetXaxis()->GetXmax();
  //   hsim_foMass->SetMarkerStyle(20);
  //   hsim_foMass->SetMarkerSize(0.7);
  //   cout << "hsim_foMass = " << hsim_foMass << endl;
  //   hsim_foMass->Draw("e");

  //   // double par_fo[3];
  //   // TF1 *fsig_fo = new TF1("fsig_fo","gaus"); //[0]*TMath::BreitWigner(x,[1],[2])
  //   // fsig_fo->SetLineStyle(9);
  //   // // fsig_fo->SetLineColor(kRed);
  //   // fsig_fo->SetParameters(15000,0.980,0.03);
  //   // hsim_foMass->Fit(fsig_fo,"","mR",0.6,1.4);

  //   // //w.factory("BreitWigner::sig1_foMass_mc(m_foMass_mc[0.7,1.2],mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.005,0.08])");//mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.005,0.08]
  //   // w.factory("Voigtian::sig1_foMass_mc(m_foMass_mc[0.7,1.2],mean1_foMass_mc[0.95,1.0],width1_foMass_mc[0.00005,0.08],sigma1_foMass_mc[0.005,0.08])");
  //   // w.factory("Chebychev::bkg_foMass_mc(m_foMass_mc,{c0_foMass_mc[-10,10],c1_foMass_mc[-10,10],c2_foMass_mc[-10,10]})");
  //   // //w.factory("SUM:::model(nsig[0,100000000]*sig, nbkg[0,100000000]*bkg)"); //nsig[0,100000000]*sig2,
  //   // w.factory("SUM:::model_foMass_mc(nsig_foMass_mc[0,100000000]*sig1_foMass_mc, nbkg_foMass_mc[0,100000000]*bkg_foMass_mc)"); //, nsig_mc[0,100000000]*sig2_mc
  //   // w.var("m_foMass_mc")->SetTitle("m_foMass_{K^{+}K^{-}} [GeV/c^{2}]");
  //   // RooDataHist dh_foMass_mc("dh_foMass_mc","dh_foMass_mc",*w.var("m_foMass_mc"),Import(*hsim_foMass));
  //   // RooPlot* fr_foMass_mc = w.var("m_foMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
  //   // fr_foMass_mc->SetTitleOffset(0.90,"X");
  //   // fr_foMass_mc->SetTitleSize(0.06,"XYZ");
  //   // fr_foMass_mc->SetLabelSize(0.06,"xyz");
  //   // w.pdf("model_foMass_mc")->fitTo(dh_foMass_mc);
  //   // // result = w.pdf("model_foMass_mc")->fitTo(dh_foMass_mc,Extended(kTRUE),Save());
  //   // dh_foMass_mc.plotOn(fr_foMass_mc,RooFit::Name("ldh_foMass_mc"));
  //   // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("sig1_foMass_mc")),LineColor(kRed),RooFit::Name("lsig1_foMass_mc"));
  //   // //w.pdf("model_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("sig2_mc")),LineColor(kMagenta),RooFit::Name("lsig2_mc"));
  //   // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc, Components(*w.pdf("bkg_foMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_foMass_mc"));
  //   // w.pdf("model_foMass_mc")->plotOn(fr_foMass_mc,RooFit::Name("lmodel_foMass_mc"));
  //   // w.pdf("model_foMass_mc")->paramOn(fr_foMass_mc,Layout(0.55,0.98,0.92));//,Parameters(RooArgSet(*w.var("nsig_mc"),*w.var("nbkg_mc"),*w.var("mean1_mc"),*w.var("mean2_mc"),*w.var("sigma1_mc"),*w.var("sigma2_mc"))));
  //   // fr_foMass_mc->Draw();

  //   // TLegend *l_foMass_mc = new TLegend(0.13,0.6,0.3,0.85);
  //   // l_foMass_mc->SetFillColor(kWhite);
  //   // l_foMass_mc->SetLineColor(kWhite);
  //   // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("ldh_foMass_mc"),"Data", "p");
  //   // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lmodel_foMass_mc"),"total","l");
  //   // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lbkg_foMass_mc"),"Background", "l");
  //   // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lsig1_foMass_mc"),"Signal1", "l");
  //   // // l_foMass_mc->AddEntry(fr_foMass_mc->findObject("lsig2_foMass_mc"),"Signal2", "l");
  //   // l_foMass_mc->Draw();

  //   // double fscale = 15.31; // fscale = (ndsig/ndbkg)/(nmcsig/nmcbkg)
  //   // double nmcsig = w.var("nsig_foMass_mc")->getVal();
  //   // double nmcsigerr = w.var("nsig_foMass_mc")->getError();
  //   // double nmcbkg = w.var("nbkg_foMass_mc")->getVal();
  //   // double nmcbkgerr = w.var("nbkg_foMass_mc")->getError();
  //   // double zmc = nmcsig/sqrt(nmcsig+(fscale*nmcbkg));
  //   // double zmcerr = sqrt(((nmcsig+(2*fscale*nmcbkg))*(nmcsig+(2*fscale*nmcbkg))*(nmcsigerr*nmcsigerr) + (fscale*fscale)*(nmcsig*nmcsig)+(nmcbkgerr*nmcbkgerr))/(4*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))*(nmcsig+(fscale*nmcbkg))));
  //   // cout<<"nmcsig = "<<nmcsig<<"| nmcbkg = "<<nmcbkg<<"| zmc = "<<zmc<<endl;
  //   // // hsim_foMass->Draw("e");

  //   // // double zObs = NumberCountingUtils::BinomialObsZ(nmcsig+nmcbkg, nmcbkg, nmcbkgerr);
  //   // // cout << " Z value (Gaussian sigma) = "<< zObs << endl;

  //   // grsim_foMass->SetPoint(j,cuts[j],zmc);
  //   // grsim_foMass->SetPointError(j,0,zmcerr);
  //   // hsim_foMass->Write();
  // }

  // csim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_foMass_%scut.root", cut.Data()), "root");
  // csim_foMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_foMass_%scut.eps", cut.Data()), "eps");

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

  TCanvas *csim_YMass = new TCanvas("csim_YMass", "csim_YMass", 1500, 800);
  csim_YMass->Divide(5, 4);
  // csim_YMass->Divide(1,1);
  TH1F *hsim_YMass;

  TCanvas *cgrsim_YMass = new TCanvas("cgrsim_YMass","cgrsim_YMass",600,400);
  cgrsim_YMass->SetGrid();
  TGraphErrors *grsim_YMass = new TGraphErrors();
  grsim_YMass->SetMarkerStyle(20);
  grsim_YMass->GetXaxis()->SetNdivisions(510, kFALSE);

    // Other options are CREATE (same as NEW), RECREATE (i.e. replace), UPDATE and READ.
  // FILE *tab = fopen("/data.local/nacer/halld_my/pimpipkmkp/fig_significance/tab.txt","w+");
  // if (tab!=NULL) printf("File opened successfully\n");

  ofstream ofs_YMass("tableau_YMass.txt", ofstream::out);

  double nsig_y[20], nsigerr_y[20];

  for (int j = 0; j < 20; j++)
  {
    csim_YMass->cd(j + 1);
    // hsim_YMass = (TH1F *)fsim->Get(Form("h_YMass_cuts_%d", j));
    hsim_YMass = new TH1F("hsim_YMass", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 2.05, 2.35);
    tsim->Project("hsim_YMass", "kpkmpippim_mf", Form("w8*(kpkmpippim_uni && kpkm_mf>1.005 && kpkm_mf<1.035 && kin_chisq<60 && abs(mm2)<0.03 && sqrt((p_x4_kin.X()*p_x4_kin.X())+(p_x4_kin.Y()*p_x4_kin.Y()))<%f)",cuts[j]));
    hsim_YMass->Rebin(4);
    // double xmin = hsim_YMass->GetXaxis()->GetXmin();
    // double xmax = hsim_YMass->GetXaxis()->GetXmax();
    hsim_YMass->SetMarkerStyle(20);
    hsim_YMass->SetMarkerSize(0.7);
    cout << "hsim_YMass = " << hsim_YMass << endl;
    hsim_YMass->SetTitle(Form("#sqrt{x^{2}_{vertex}+y^{2}_{vertex}} < %.2f",cuts[j]));
    hsim_YMass->Draw("e");

    // // w.factory("BreitWigner::sig1_YMass_mc(m_YMass_mc[1.6,3.2],mean1_YMass_mc[2.10,2.30],width1_YMass_mc[0.05,1.0])");
    // w.factory("Voigtian::sig1_YMass_mc(m_YMass_mc[1.6,3.2],mean1_YMass_mc[1.9,2.3],width1_YMass_mc[0.001,0.1],sigma1_YMass_mc[0.001,0.1])");
    // // w.factory("Landau::sig2_YMass_mc(m_YMass_mc[1.6,3.2],mean2_YMass_mc[2.0,2.4],sigma1_YMass_mc[0.001,0.1])");
    // // w.factory("CBShape::sig1_YMass_mc(m_YMass_mc[1.6,3.2],mean_YMass_mc[2.0,2.2],sigma_YMass_mc[0.001,0.1],alpha_YMass_mc[-10,1.0],n_YMass_mc[0,100])");
    // // w.factory("Gaussian::sig1_YMass_mc(m_YMass_mc[1.6,3.2],mean1_YMass_mc[2.0,2.3],width1_YMass_mc[0.01,0.1])");
    // // w.factory("Gaussian:sig2_YMass_mc(m_YMass_mc[1.6,3.2],mean2_YMass_mc[2.0,2.3],width2_YMass_mc[0.01,0.1])");
    // w.factory("Chebychev::bkg_YMass_mc(m_YMass_mc,{c0_YMass_mc[-10,10]})");//,c1_YMass_mc[-10,10]})");
    // w.factory("SUM:::model_YMass_mc(nsig_YMass_mc[0,100000000]*sig1_YMass_mc, nbkg_YMass_mc[0,100000000]*bkg_YMass_mc)");//
    // w.var("m_YMass_mc")->SetTitle("m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}]");
    // w.var("nsig_YMass_mc")->SetTitle("N_{sig}");
    // RooDataHist dh_YMass_mc("dh_YMass_mc","dh_YMass_mc",*w.var("m_YMass_mc"),Import(*hsim_YMass));
    // RooPlot* fr_YMass_mc = w.var("m_YMass_mc")->frame(Title("invariant mass K^{+}K^{-}"));
    // fr_YMass_mc->SetTitleOffset(0.90,"X");
    // fr_YMass_mc->SetTitleSize(0.06,"XYZ");
    // fr_YMass_mc->SetLabelSize(0.06,"xyz");
    // w.pdf("model_YMass_mc")->fitTo(dh_YMass_mc);
    // // result = w.pdf("model_YMass_mc")->fitTo(dh_YMass_mc,Extended(kTRUE),Save());
    // dh_YMass_mc.plotOn(fr_YMass_mc,RooFit::Name("ldh_YMass_mc"));
    // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig1_YMass_mc")),LineColor(kRed),RooFit::Name("lsig1_YMass_mc"));
    // // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("sig2_YMass_mc")),LineColor(kRed),RooFit::Name("lsig2_YMass_mc"));
    // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc, Components(*w.pdf("bkg_YMass_mc")),LineStyle(kDashed),LineColor(28),RooFit::Name("lbkg_YMass_mc"));
    // w.pdf("model_YMass_mc")->plotOn(fr_YMass_mc,RooFit::Name("lmodel_YMass_mc"));
    // w.pdf("model_YMass_mc")->paramOn(fr_YMass_mc,Layout(0.5,0.90,0.90) ,Parameters(RooArgSet(*w.var("nsig_YMass_mc"))));
    // fr_YMass_mc->Draw();

    // TLegend *l_YMass_mc = new TLegend(0.13,0.6,0.3,0.85);
    // l_YMass_mc->SetFillColor(kWhite);
    // l_YMass_mc->SetLineColor(kWhite);
    // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("ldh_YMass_mc"),"Data", "p");
    // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lmodel_YMass_mc"),"total","l");
    // // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lbkg_YMass_mc"),"Background", "l");
    // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lsig1_YMass_mc"),"gaus1", "l");
    // // l_YMass_mc->AddEntry(fr_YMass_mc->findObject("lsig2_YMass_mc"),"gaus2", "l");
    // l_YMass_mc->Draw();

    // nsig_y[j] = w.var("nsig_YMass_mc")->getVal();
    // nsigerr_y[j] = w.var("nsig_YMass_mc")->getError();
    nsig_y[j] = hsim_YMass->Integral(hsim_YMass->GetXaxis()->GetBinLowEdge(1), hsim_YMass->GetXaxis()->GetBinUpEdge(200));
    nsigerr_y[j] = sqrt(nsig_y[j]);
    double z_y = nsig_y[j] / sqrt(nbkg_y[j]);
    double zerr_y = sqrt(((4 * nbkg_y[j] * nbkg_y[j] * nsigerr_y[j] * nsigerr_y[j]) + (nsig_y[j] * nsig_y[j] * nbkgerr_y[j] * nbkgerr_y[j])) / (4 * nbkg_y[j] * nbkg_y[j] * nbkg_y[j]));

    ofs_YMass << " j = " << j << " | cuts = " << cuts[j] << " | nsig_y = " << nsig_y[j] << " | nsigerr_y = " << nsigerr_y[j] << " | nbkg_y = " << nbkg_y[j] << " | nbkgerr_y = " << nbkgerr_y[j] << " | z_y = " << z_y << " | zerr_y = " << zerr_y << endl;

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

    grsim_YMass->SetPoint(j, cuts[j], z_y);
    grsim_YMass->SetPointError(j, 0, zerr_y);
    hsim_YMass->Write();
  }

  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%scut.root", cut.Data()), "root");
  csim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/csim_YMass_%scut.eps", cut.Data()), "eps");

  cgrsim_YMass->cd();
  grsim_YMass->Draw("AP");

  // // TF1 *fsimpol6 = new TF1("fsimpol6","pol6");
  // // grsim_YMass->Fit("fsimpol6","mR","", cuts[12], cuts[4]);
  // // double zmc_max = fsimpol6->GetMaximumX(cuts[12], cuts[4]);
  // // cout<<"zmc_max = "<<zmc_max<<endl;

  // // result->Print();
  grsim_YMass->Write();
  grsim_YMass->SetTitle("Significance Y(2175) Vs selection (MC); cut; S/#sqrt{B}");
  cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%scut.root",cut.Data()),"root");
  cgrsim_YMass->Print(Form("/data.local/nacer/halld_my/pippimkpkm/fig_significance/cgrsim_YMass_%scut.eps",cut.Data()),"eps");
  ofs_PhiMass.close();
  ofs_YMass.close();
  outputfig->Print();
  // outputfig->Close();
}

*/