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

void xsection()
{
    TFile *fsim = new TFile("/Users/nacer/halld_my/pippimkpkm/input/phifo_genr8_17v03.root");
    TTree *tsim = (TTree *)fsim->Get("ntp");
    TFile *fps = new TFile("/Users/nacer/halld_my/pippimkpkm/input/flux_30274_31057.root");
    TFile *f = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_scanphi/scanphi.root");
    TFile *outputfig = new TFile("/Users/nacer/halld_my/pippimkpkm/fig_xsection/xsection.root", "UPDATE");
    ofstream ofs_xsection("xsection.txt", ofstream::out);

    TString cut = "&& kpkm_mf>1.005 && kpkm_mf<1.035 && kin_chisq<30 && abs(mm2)<0.015 && -t_kin<1 && beam_p4_kin.E()>6 && pt<0.4";

    double mkk_min = 0.98, mkk_max = 1.2;
    double m2pi_min = 0.3, m2pi_max = 1.3;
    double m2pi2k_min = 1.6, m2pi2k_max = 3.2;

    // +++ Phi(1020) Mass in MC
    TCanvas *csim_PhiMass = new TCanvas("csim_PhiMass", "csim_PhiMass", 600, 400);
    csim_PhiMass->cd();
    TH1F *hsim_PhiMass = new TH1F("hsim_PhiMass", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 50, 0.98, 1.2);
    tsim->Project("hsim_PhiMass", "kpkm_mf", "w8*(kpkm_uni"+cut+")");
    hsim_PhiMass->Draw("hist");
    
    // +++ Y(2175) Mass in MC
    TCanvas *csim_YMass = new TCanvas("csim_YMass", "csim_YMass", 600, 400);
    csim_YMass->cd();
    TH1F *hsim_YMass = new TH1F("hsim_YMass", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 50, 1.6, 3.2);
    tsim->Project("hsim_YMass", "kpkmpippim_mf", "w8*(kpkmpippim_uni"+cut+")");
    hsim_YMass->Draw("hist");    

    // +++ PS flux
    TCanvas *c_tagged_flux = new TCanvas("c_tagged_flux", "c_tagged_flux", 600, 400);
    c_tagged_flux->cd();
    TH1F *h_tagged_flux = (TH1F *)fps->Get("tagged_flux");
    cout<<"h_tagged_flux"<<h_tagged_flux<<endl;
    h_tagged_flux->Rebin(20);
    h_tagged_flux->Draw("hist");

    // +++ Phi(1020) Pi+ Pi- Mass in Data
    TCanvas *cphiy = new TCanvas("cphiy", "cphiy", 600, 400);
    cphiy->cd();
    TGraphErrors *grphiy = (TGraphErrors *)f->Get("grphiy_1"); //new TGraphErrors(n2pi2k);
    cout<<"grphiy"<<grphiy<<endl;
    grphiy->SetMarkerStyle(20);
    grphiy->SetMarkerSize(1.0);
    grphiy->SetMarkerColor(1);
    grphiy->GetHistogram()->SetTitle(";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];N_{#phi}");
    grphiy->Draw("AP");

    TF1 *fsb = new TF1("fsb", "[0]*TMath::BreitWigner(x,[1],[2]) + [3]*TMath::BreitWigner(x,[4],[5]) + pol2(6)", 1.6, 3.2); // + [3]*TMath::BreitWigner(x,[4],[5])
    fsb->SetLineColor(1);
    fsb->SetParameters(1000,1.8,0.8, 1000,2.4,0.8, 1,1,1);
    fsb->SetParLimits(0, 1, 10000);
    fsb->SetParLimits(1, 1.7, 1.9);
    fsb->SetParLimits(2, 0.01, 1.0);
    fsb->SetParLimits(3, 1, 10000);
    fsb->SetParLimits(4, 2.1, 2.5);
    fsb->SetParLimits(5, 0.01, 1.0); 

    TF1 *fs1 = new TF1("fs1", "[0]*TMath::BreitWigner(x,[1],[2])", 1.6, 3.2);
    fs1->SetLineColor(4);
    TF1 *fs2 = new TF1("fs2", "[3]*TMath::BreitWigner(x,[4],[5])", 1.6, 3.2);
    fs2->SetLineColor(4);    
    TF1 *fb1 = new TF1("fb1", "pol2(6)", 1.6, 3.2);
    fb1->SetLineColor(2);
    // TF1 *fs2 = new TF1("fs2", "[3]*TMath::BreitWigner(x,[4],[5])", 1.6, 3.2);
    // fs2->SetLineColor(4);

    grphiy->Fit("fsb", "", "", 1.6, 3.2);
    double par[6];
    fsb->GetParameters(&par[0]);
    fs1->SetParameters(&par[0]);
    fs2->SetParameters(&par[3]);
    fb1->SetParameters(&par[6]);
    // fs2->SetParameters(&par[3]);

    // grphiy->Draw("AP");
    fsb->DrawCopy("same");
    fb1->DrawCopy("same");
    fs1->DrawCopy("same");
    fs2->DrawCopy("same");

    /*
    // +++ xsection
    TCanvas *c_xsection = new TCanvas("c_xsection", "c_xsection", 600, 400);
    TGraphErrors *gxsection = new TGraphErrors(n2pi2k);
    gxsection->SetMarkerStyle(20);
    gxsection->SetMarkerSize(1.0);
    gxsection->SetMarkerColor(1);
    gxsection->GetHistogram()->SetTitle(";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];BR #sigma_{Y -> #phi f_{0}} [nb]");

    double eff_phi = hsim_YMass->Integral(1, 50) / 10000000; // Effeciency = N_observed/N_generated
    double lumi_phi = h_tagged_flux->Integral(1, 50) * 1.26;   // Luminosity = N_gama * T ,  T = 1.26 barns
    double br_phi = 0.489;
    double br_y = 1;

    double xsection_phiy = 1e9 * w.var("nsig_phiy")->getVal() / (eff_phi * lumi_phi * br_phi * br_y);
    ofs_xsection << " i = " << i << " | eff_phi = " << eff_phi << " | lumi_phi = " << lumi_phi << " | xsection_phiy = " << xsection_phiy << endl;
    gxsection->SetPoint(i - 1, h2d1->GetXaxis()->GetBinCenter(i), xsection_phiy);
    gxsection->SetPointError(i - 1, 0, 0); //dN1

    c_xsection->cd();
    gxsection->Draw("AP");
    c_xsection->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_xsection.root", "root");
    c_xsection->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_xsection.eps", "eps");
*/
    c_tagged_flux->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_tagged_flux.root", "root");
    c_tagged_flux->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_tagged_flux.eps", "eps");
    csim_PhiMass->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/csim_PhiMass.root", "root");
    csim_PhiMass->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/csim_PhiMass.eps", "eps");
    csim_YMass->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/csim_YMass.root", "root");
    csim_YMass->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/csim_YMass.eps", "eps");
    cphiy->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_grphiy.root", "root");
    cphiy->Print("/Users/nacer/halld_my/pippimkpkm/fig_xsection/c_grphiy.eps", "eps");          

}