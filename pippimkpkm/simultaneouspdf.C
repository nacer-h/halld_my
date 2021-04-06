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
#include "TLine.h"
#include "profileLikelihoodScan.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
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
#include "RooDataHist.h"
#include "RooCmdArg.h"

using namespace RooFit;

void simultaneouspdf(TString name)
{
    TFile *finput = new TFile(Form("output/fig_%s/ul_%s.root", name.Data(), name.Data()));

    FILE *table_simpdf_sys = fopen(Form("table_simpdf_sys_%s.tex", name.Data()), "w");
    fprintf(table_simpdf_sys, "\\documentclass[11pt]{extarticle}\n \\usepackage[margin=0.1in]{geometry}\n \\usepackage{tabularx}\n \\usepackage{caption} \n \\usepackage{makecell} \n \\captionsetup{labelformat=empty}\n \\usepackage{siunitx}\n \\begin{document}\n \\begin{table}[!htbp]\n \\small\n \\centering\n \\caption{Systematic errors in beam energy bins}\n \\begin{tabular}{|c|c|c|c|c|}\n \\hline\n \\thead{Mean$\\pm 0.01$ \\\\$\\pm \\delta\\mu_{PDG}$} & \\thead{Width$\\pm 0.012$\\\\$\\pm \\delta\\Gamma_{PDG}$} & \\thead{Bkg deg\\\\(3,4,5)}& \\thead{Fit range\\\\(2; 3, 2.9, 3.1)} & $\\epsilon \\pm \\delta \\epsilon$ \\\\ \n \\hline\n");

    ofstream test;
    test.open("test.txt");

    // get the the histograms of the different datasets
    //-------------------------------------------------
    TCanvas *c_hist = new TCanvas("c_hist", "c_hist", 1400, 300);
    c_hist->Divide(3);
    TH1F *h_17 = (TH1F *)finput->Get("h17_phiy");
    h_17->SetTitle("2017");
    cout << "h_17 = " << h_17 << endl;
    c_hist->cd(1);
    h_17->Draw("e");
    TH1F *h_18 = (TH1F *)finput->Get("h18_phiy");
    h_18->SetTitle("Spring 2018");
    cout << "h_18 = " << h_18 << endl;
    c_hist->cd(2);
    h_18->Draw("e");
    TH1F *h_18l = (TH1F *)finput->Get("h18l_phiy");
    h_18l->SetTitle("Fall 2018");
    cout << "h_18l = " << h_18l << endl;
    c_hist->cd(3);
    h_18l->Draw("e");

    c_hist->Print(Form("output/fig_%s/chist_%s.root", name.Data(), name.Data()), "root");
    c_hist->Print(Form("output/fig_%s/chist_%s.pdf", name.Data(), name.Data()), "pdf");
    c_hist->Print(Form("output/fig_%s/chist_%s.eps", name.Data(), name.Data()), "eps");

    // Create models for different datasets with same 'width'
    // ------------------------------------------------------
    RooWorkspace w("w", kTRUE);

    //--------- define the xsec as common parameter for the simultaneous fit
    // --- xsec_y = 1e9 * N_y_data / (eff_y * lumi_y * T * br_phi) ---
    double xmin = 2.0, xmax = 3.0; // 3.0, 2.9, 3.1
    double y_mean = 2.192;         // +-0.01
    double y_width = 0.083;        // +-0.012
    double eff_17 = 0.076;         // +-0.0003
    double eff_18 = 0.048;         // +-0.0003
    double eff_18l = 0.099;        // +-0.0003

    if (name == "yphifo")
    {
        y_mean = 2.194; // +-0.01
        eff_17 = 0.084;
        eff_18 = 0.053;
        eff_18l = 0.12;
    }

    w.factory("xsec[-0.03,-1000,+1000]"); // free parameter let to float and use to extract -Log(L)
    w.factory(Form("prod::nsig_17(xsec, %f, 5.89545e+13, 1.273, 0.492, 1e-9)", eff_17));
    w.factory(Form("prod::nsig_18(xsec, %f, 1.69415e+14, 1.273, 0.492, 1e-9)", eff_18));
    w.factory(Form("prod::nsig_18l(xsec, %f, 1.09564e+14, 1.273, 0.492, 1e-9)", eff_18l));

    w.factory(Form("Voigtian::sig_17(m[%f, %f],mean_17[%f],width_17[%f],sigma_17[0.025])", xmin, xmax, y_mean, y_width));
    if (name == "yphi2pi") w.factory("Chebychev::bkg_17(m,{c0_17[-10,+10], c1_17[-10,+10], c2_17[-10,+10], c3_17[-10,+10]})"); //c0_17[-10,+10], c1_17[-10,+10], c2_17[-10,+10], c3_17[-10,+10]
    if (name == "yphifo") w.factory("Chebychev::bkg_17(m,{c0_17[-10,+10], c1_17[-10,+10], c2_17[-10,+10]})"); //c0_17[-10,+10], c1_17[-10,+10], c2_17[-10,+10]
    w.factory("SUM::model_17(nsig_17*sig_17, nbkg_17[0,+1000000]*bkg_17)"); //[0,-1000000,+1000000]
    RooDataHist dh_17("dh_17", "dh_17", *w.var("m"), Import(*h_17));
    
    w.factory(Form("Voigtian::sig_18(m[%f, %f],mean_18[%f],width_18[%f],sigma_18[0.026])", xmin, xmax, y_mean, y_width));
    if (name == "yphi2pi") w.factory("Chebychev::bkg_18(m,{c0_18[-10,+10], c1_18[-10,+10], c2_18[-10,+10], c3_18[-10,+10]})");
    if (name == "yphifo") w.factory("Chebychev::bkg_18(m,{c0_18[-10,+10], c1_18[-10,+10], c2_18[-10,+10]})");
    w.factory("SUM::model_18(nsig_18*sig_18, nbkg_18[0,+1000000]*bkg_18)");
    RooDataHist dh_18("dh_18", "dh_18", *w.var("m"), Import(*h_18));
    
    w.factory(Form("Voigtian::sig_18l(m[%f, %f],mean_18l[%f],width_18l[%f],sigma_18l[0.024])", xmin, xmax, y_mean, y_width));
    if (name == "yphi2pi") w.factory("Chebychev::bkg_18l(m,{c0_18l[-10,+10], c1_18l[-10,+10], c2_18l[-10,+10], c3_18l[-10,+10]})");
    if (name == "yphifo") w.factory("Chebychev::bkg_18l(m,{c0_18l[-10,+10], c1_18l[-10,+10], c2_18l[-10,+10]})");
    w.factory("SUM::model_18l(nsig_18l*sig_18l, nbkg_18l[0,+1000000]*bkg_18l)");
    RooDataHist dh_18l("dh_18l", "dh_18l", *w.var("m"), Import(*h_18l));

    // C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
    // ---------------------------------------------------------------------------
    // Define category to distinguish different datasets
    RooCategory sample("sample", "sample");
    sample.defineType("17");
    sample.defineType("18");
    sample.defineType("18l");

    // Construct combined dataset in (x,sample)
    RooDataHist combData("combData", "combined data", *w.var("m"), Index(sample), Import("17", dh_17), Import("18", dh_18), Import("18l", dh_18l));

    // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
    // -----------------------------------------------------------------------------------
    // Construct a simultaneous pdf using category sample as index
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);
    // Associate models with corresponding datasets
    simPdf.addPdf(*w.pdf("model_17"), "17");
    simPdf.addPdf(*w.pdf("model_18"), "18");
    simPdf.addPdf(*w.pdf("model_18l"), "18l");

    // P e r f o r m   a   s i m u l t a n e o u s   f i t
    // ---------------------------------------------------
    // Perform simultaneous fit of model to the corresponding dataset
    // simPdf.fitTo(combData);
    RooFitResult *fitres = simPdf.fitTo(combData,  Save(true), InitialHesse(true), Extended(true));

    // P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s
    // ----------------------------------------------------------------

    TCanvas *c = new TCanvas("c_simpdf", "c_simpdf", 1400, 300);
    c->Divide(3);

    // RooAbsPdf *pdf_17 = w.pdf("model_17");
    // RooAbsPdf *pdf_18 = w.pdf("model_18");
    // RooAbsPdf *pdf_18l  = w.pdf("model_18l");

    RooPlot *fr_17 = w.var("m")->frame(Title("2017;m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2});"));
    // w.pdf("model_17")->fitTo(dh_17);
    dh_17.plotOn(fr_17, RooFit::Name("ldh_17"));
    w.pdf("model_17")->plotOn(fr_17, Components(*w.pdf("sig_17")), LineColor(kRed), RooFit::Name("lsig_17"));
    w.pdf("model_17")->plotOn(fr_17, Components(*w.pdf("bkg_17")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_17"));
    w.pdf("model_17")->plotOn(fr_17, RooFit::Name("lmodel_17"));
    gPad->SetLeftMargin(0.15);
    c->cd(1);
    fr_17->Draw();
    TLatex lat_17;
    lat_17.SetTextSize(0.05);
    lat_17.SetTextAlign(13); //align at top
    lat_17.SetNDC();
    lat_17.SetTextColor(kBlue);
    lat_17.DrawLatex(0.53, 0.46, Form("xsec = %0.2f#pm%0.2f", w.var("xsec")->getVal(), w.var("xsec")->getError()));
    lat_17.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f#pm%0.2f", w.function("nsig_17")->getVal(), w.function("nsig_17")->getPropagatedError(*fitres)));
    // lat_17.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f", w.function("nsig_17")->getVal()));
    lat_17.DrawLatex(0.53, 0.34, Form("#mu = %0.3f#pm%0.3f", w.var("mean_17")->getVal(), w.var("mean_17")->getError()));
    lat_17.DrawLatex(0.53, 0.28, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_17")->getVal(), w.var("sigma_17")->getError()));
    lat_17.DrawLatex(0.53, 0.22, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_17")->getVal(), w.var("width_17")->getError()));

    RooPlot *fr_18 = w.var("m")->frame(Title("Spring 2018;m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2});"));
    // w.pdf("model_18")->fitTo(dh_18);
    dh_18.plotOn(fr_18, RooFit::Name("ldh_18"));
    w.pdf("model_18")->plotOn(fr_18, Components(*w.pdf("sig_18")), LineColor(kRed), RooFit::Name("lsig_18"));
    w.pdf("model_18")->plotOn(fr_18, Components(*w.pdf("bkg_18")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_18"));
    w.pdf("model_18")->plotOn(fr_18, RooFit::Name("lmodel_18"));
    c->cd(2);
    fr_18->Draw();
    TLatex lat_18;
    lat_18.SetTextSize(0.05);
    lat_18.SetTextAlign(13); //align at top
    lat_18.SetNDC();
    lat_18.SetTextColor(kBlue);
    lat_18.DrawLatex(0.53, 0.46, Form("xsec = %0.2f#pm%0.2f", w.var("xsec")->getVal(), w.var("xsec")->getError()));
    lat_18.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f#pm%0.2f", w.function("nsig_18")->getVal(), w.function("nsig_18")->getPropagatedError(*fitres)));
    // lat_18.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f", w.function("nsig_18")->getVal()));
    lat_18.DrawLatex(0.53, 0.34, Form("#mu = %0.3f#pm%0.3f", w.var("mean_18")->getVal(), w.var("mean_18")->getError()));
    lat_18.DrawLatex(0.53, 0.28, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_18")->getVal(), w.var("sigma_18")->getError()));
    lat_18.DrawLatex(0.53, 0.22, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_18")->getVal(), w.var("width_18")->getError()));

    RooPlot *fr_18l = w.var("m")->frame(Title("Fall 2018;m_{#phi#pi^{+}#pi^{-}} (GeV/c^{2});"));
    // w.pdf("model_18l")->fitTo(dh_18l);
    dh_18l.plotOn(fr_18l, RooFit::Name("ldh_18l"));
    w.pdf("model_18l")->plotOn(fr_18l, Components(*w.pdf("sig_18l")), LineColor(kRed), RooFit::Name("lsig_18l"));
    w.pdf("model_18l")->plotOn(fr_18l, Components(*w.pdf("bkg_18l")), LineStyle(kDashed), LineColor(28), RooFit::Name("lbkg_18l"));
    w.pdf("model_18l")->plotOn(fr_18l, RooFit::Name("lmodel_18l"));
    c->cd(3);
    fr_18l->Draw();
    TLatex lat_18l;
    lat_18l.SetTextSize(0.05);
    lat_18l.SetTextAlign(13); //align at top
    lat_18l.SetNDC();
    lat_18l.SetTextColor(kBlue);
    lat_18l.DrawLatex(0.53, 0.46, Form("xsec = %0.2f#pm%0.2f", w.var("xsec")->getVal(), w.var("xsec")->getError()));
    lat_18l.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f#pm%0.2f", w.function("nsig_18l")->getVal(), w.function("nsig_18l")->getPropagatedError(*fitres)));
    // lat_18l.DrawLatex(0.53, 0.40, Form("N_{sig} = %0.2f", w.function("nsig_18l")->getVal()));
    lat_18l.DrawLatex(0.53, 0.34, Form("#mu = %0.3f#pm%0.3f", w.var("mean_18l")->getVal(), w.var("mean_18l")->getError()));
    lat_18l.DrawLatex(0.53, 0.28, Form("#sigma = %0.3f#pm%0.3f", w.var("sigma_18l")->getVal(), w.var("sigma_18l")->getError()));
    lat_18l.DrawLatex(0.53, 0.22, Form("#Gamma = %0.3f#pm%0.3f", w.var("width_18l")->getVal(), w.var("width_18l")->getError()));

    c->Print(Form("output/fig_%s/csimpdf_%s.root", name.Data(), name.Data()), "root");
    c->Print(Form("output/fig_%s/csimpdf_%s.pdf", name.Data(), name.Data()), "pdf");
    c->Print(Form("output/fig_%s/csimpdf_%s.eps", name.Data(), name.Data()), "eps");

    // +++++++++++++++++++ Systematic Errors +++++++++++++++++++
    // vary poln (all 2 and 4), all together fit range (same time), eff+-deff (same - & +), vary mean and width of Y(2175) (same time).

    // if(name == "yphi2pi")
    double arr_sig_mean_yphi2pi[3] = {0.06, -0.19, 0.34};  //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
    double arr_sig_width_yphi2pi[3] = {0.06, 0.01, 0.10}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
    double arr_bkg_poln_yphi2pi[3] = {0.06, 0.32, 0.04};  //pol:4,3,5
    double arr_fit_range_yphi2pi[3] = {0.06, 0.12, -0.23}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
    double arr_eff_yphi2pi[3] = {0.06, 0.06, 0.06}; //eff: +-0.0003 -- 3 sigma
    
    // if(name == "yphifo")
    double arr_sig_mean_yphifo[3] = {0.08, 0.10, 0.07};  //mean: mu, -0.01, +0.01 -- 2.188+/-0.01
    double arr_sig_width_yphifo[3] = {0.08, 0.06, 0.10}; //width: gamma, -0.012, +0.012 -- 0.083+/-0.012
    double arr_bkg_poln_yphifo[3] = {0.08, 0.25, -0.20};  //pol:4,3,5
    double arr_fit_range_yphifo[3] = {0.08, 0.04, 0.08}; //fit: [2.0, 3.0], [2.0, 2.9], [2.0, 3.1]
    double arr_eff_yphifo[3] = {0.08, 0.08, 0.08}; //eff: +-0.0003 -- 3 sigma
    
    double arr_sig_mean[3], arr_sig_width[3], arr_bkg_poln[3], arr_fit_range[3], arr_eff[3];

    if (name == "yphi2pi")
    {
        memcpy(arr_sig_mean, arr_sig_mean_yphi2pi, sizeof(arr_sig_mean));
        memcpy(arr_sig_width, arr_sig_width_yphi2pi, sizeof(arr_sig_width));
        memcpy(arr_bkg_poln, arr_bkg_poln_yphi2pi, sizeof(arr_bkg_poln));
        memcpy(arr_fit_range, arr_fit_range_yphi2pi, sizeof(arr_fit_range));
        memcpy(arr_eff, arr_eff_yphi2pi, sizeof(arr_eff));
    }

    if (name == "yphifo")
    {
        memcpy(arr_sig_mean, arr_sig_mean_yphifo, sizeof(arr_sig_mean));
        memcpy(arr_sig_width, arr_sig_width_yphifo, sizeof(arr_sig_width));
        memcpy(arr_bkg_poln, arr_bkg_poln_yphifo, sizeof(arr_bkg_poln));
        memcpy(arr_fit_range, arr_fit_range_yphifo, sizeof(arr_fit_range));
        memcpy(arr_eff, arr_eff_yphifo, sizeof(arr_eff));
    }

    double err_sig_mean = TMath::RMS(3, arr_sig_mean) / arr_sig_mean[0];
    double err_sig_width = TMath::RMS(3, arr_sig_width) / arr_sig_width[0];
    double err_bkg_poln = TMath::RMS(3, arr_bkg_poln) / arr_bkg_poln[0];
    double err_fit_range = TMath::RMS(3, arr_fit_range) / arr_fit_range[0];
    double err_eff = TMath::RMS(3, arr_eff) / arr_eff[0];

    double err_sys = TMath::Sqrt(err_sig_mean*err_sig_mean + err_sig_width*err_sig_width + err_bkg_poln*err_bkg_poln + err_fit_range*err_fit_range + err_eff*err_eff);

    double dxsec_stat = w.var("xsec")->getError();
    double dxsec_sys = w.var("xsec")->getVal() * err_sys;
    double dxsec_tot = TMath::Sqrt(dxsec_stat*dxsec_stat + dxsec_sys*dxsec_sys);

    // ++++++++++++++++++ Upper Limit 90% CL  +++++++++++++++++
    // ++++ profile likelihood method
    TGraph *gprof_simpdf = profileLikelihoodScan(combData, simPdf, w.var("xsec"), w.var("xsec")->getVal() - 5 * w.var("xsec")->getError(), w.var("xsec")->getVal() + 5 * w.var("xsec")->getError());

    TCanvas *cprof_simpdf = new TCanvas("cprof_simpdf", "", 600, 400);
    cprof_simpdf->cd();
    cprof_simpdf->SetGrid();
    gprof_simpdf->SetTitle("; #sigma_{#phi#pi^{+}#pi^{-}} [nb]; -Log(L)");
    gprof_simpdf->SetLineWidth(2);
    gprof_simpdf->SetMarkerStyle(20);
    gprof_simpdf->SetMarkerSize(1.5);
    gprof_simpdf->Draw("AP");
    gprof_simpdf->Fit("pol2", "R", "", gprof_simpdf->GetHistogram()->GetXaxis()->GetXmin(), gprof_simpdf->GetHistogram()->GetXaxis()->GetXmax());
    cprof_simpdf->Print(Form("output/fig_%s/cprof_simpdf_%s.root", name.Data(), name.Data()), "root");
    cprof_simpdf->Print(Form("output/fig_%s/cprof_simpdf_%s.pdf", name.Data(), name.Data()), "pdf");
    cprof_simpdf->Print(Form("output/fig_%s/cprof_simpdf_%s.eps", name.Data(), name.Data()), "eps");

    TCanvas *cprofxsec_simpdf = new TCanvas("cprofxsec_simpdf", "", 600, 400);
    TGraph *gprofxsec_simpdf = new TGraph();
    for (int i = 0; i < gprof_simpdf->GetN(); i++)
    {
        double x_gprof_simpdf, y_gprof_simpdf;
        gprof_simpdf->GetPoint(i, x_gprof_simpdf, y_gprof_simpdf);
        gprofxsec_simpdf->SetPoint(i, x_gprof_simpdf, exp(-y_gprof_simpdf));
    }
    cprofxsec_simpdf->cd();
    gprofxsec_simpdf->SetTitle("; #sigma [nb]; Likelihood");
    gprofxsec_simpdf->SetLineWidth(2);
    gprofxsec_simpdf->SetMarkerStyle(20);
    gprofxsec_simpdf->SetMarkerSize(1.5);
    gprofxsec_simpdf->Draw("AP");
    TF1 *func_profxsec_simpdf = new TF1("func_profxsec_simpdf", "gaus", gprofxsec_simpdf->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec_simpdf->GetHistogram()->GetXaxis()->GetXmax());
    gprofxsec_simpdf->Fit(func_profxsec_simpdf, "R", "", gprofxsec_simpdf->GetHistogram()->GetXaxis()->GetXmin(), gprofxsec_simpdf->GetHistogram()->GetXaxis()->GetXmax());

    //creat a convolution between the profile likelihood and gaus with same mean and a s.t.d as sys err
    TCanvas *c_twogaus_conv = new TCanvas("c_twogaus_conv", "", 600, 400);
    c_twogaus_conv->cd();
    TF1 *func_twogaus_conv = new TF1("func_twogaus_conv", "gaus", -1.5, 1.5);
    func_twogaus_conv->SetParameters(func_profxsec_simpdf->GetParameter(0), func_profxsec_simpdf->GetParameter(1), dxsec_tot);
    func_twogaus_conv->SetTitle(";#sigma (nb);Likelihood");
    func_twogaus_conv->Draw();

    double ul = func_twogaus_conv->GetParameter(1) + ROOT::Math::gaussian_quantile(0.9 + 0.1 * ROOT::Math::gaussian_cdf(0, func_twogaus_conv->GetParameter(2), func_twogaus_conv->GetParameter(1)), func_twogaus_conv->GetParameter(2)); //bayesian

    TLine *l_ul = new TLine(ul, 0, ul,func_twogaus_conv->Eval(ul));
    l_ul->SetLineWidth(2);
    l_ul->SetLineColor(kBlue);
    l_ul->Draw("same");

    cprofxsec_simpdf->Print(Form("output/fig_%s/cprofxsec_simpdf_%s.root", name.Data(), name.Data()), "root");
    cprofxsec_simpdf->Print(Form("output/fig_%s/cprofxsec_simpdf_%s.pdf", name.Data(), name.Data()), "pdf");
    cprofxsec_simpdf->Print(Form("output/fig_%s/cprofxsec_simpdf_%s.eps", name.Data(), name.Data()), "eps");
    
    c_twogaus_conv->Print(Form("output/fig_%s/c_twogaus_conv_simpdf_%s.root", name.Data(), name.Data()), "root");
    c_twogaus_conv->Print(Form("output/fig_%s/c_twogaus_conv_simpdf_%s.pdf", name.Data(), name.Data()), "pdf");
    c_twogaus_conv->Print(Form("output/fig_%s/c_twogaus_conv_simpdf_%s.eps", name.Data(), name.Data()), "eps");

    fprintf(table_simpdf_sys, "%0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\ \n", err_sig_mean * 100, err_sig_width * 100, err_bkg_poln * 100, err_fit_range * 100, err_eff * 100);

    test << std::setprecision(2) << std::fixed;
    test << "\u03C3 = " << w.var("xsec")->getVal() << "\u00b1" << dxsec_stat <<  "\u00b1" << dxsec_sys <<  "\u00b1" << dxsec_tot << endl;
    test << "UL = " << ul << endl;

    fprintf(table_simpdf_sys, "\\hline\n \\end{tabular}\n \\end{table}\n \\end{document}\n");
    fclose(table_simpdf_sys);
    gSystem->Exec(Form("pdflatex table_simpdf_sys_%s.tex", name.Data()));

}