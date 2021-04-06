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

TGraph* profileLikelihoodScan(RooDataHist &combData, RooAbsPdf &simPdf, RooRealVar *xsec, double xmin = -10, double xmax = 90, double steps=100, double preminll=0)
{
	bool verbose = true;
	TGraph *res = new TGraph(steps);

	double curx = xmin;
	double dx   = (xmax-xmin)/(steps-1);
	
	double oldxs = xsec->getVal();
	
	double minll  = 1e8;
	double minllx = (xmax-xmin)/2.; 
	int minllbin  = -1;
	
	int cnt=0;
	
	 if (verbose) cout <<"LL scan: "<<flush;
	while (cnt<steps)
	{
		xsec->setVal(curx);
		xsec->setConstant(true);
		RooFitResult *fitres = 0;
		int fitcnt = 0;

		// and fit it to the combined data set for BG only
		while (++fitcnt<10 && (fitres==0 || fitres->covQual()!=3 ) ) 
		{
			if (fitres) delete fitres;
			fitres = simPdf.fitTo(combData, Save(true), InitialHesse(true), Extended(true)/*, Range("xR1")*/);
		}
		
		double nll = 0;
		if (fitres && fitres->covQual()==3)
		{
			nll = fitres->minNll();
			
			if (minll > nll) { minll = nll; minllx = curx; minllbin = cnt; }
			
			if (cnt%10==0) 
			{
				if (verbose) cout<<cnt/10<<flush; 
			}
			else 
			{
				if (verbose) cout <<"("<<nll-preminll<<") "<<flush;
			}
			delete fitres;
		}
		else if (verbose) cout <<"x"<<flush;
		
		res->SetPoint(cnt++,curx,nll);
		curx += dx;
	}
	
	for (int i=0;i<steps;++i)
	{
		double x,y, xl, yl, xh, yh;
		res->GetPoint(i,x,y);
		
		if (i>0 && i<steps-1)
		{
			res->GetPoint(i,xl,yl);
			res->GetPoint(i,xh,yh);
			
			if (fabs(y-0.5*(yl+yh))>2*fabs(yl-yh)) {cout <<"("<<i<<":"<<y<<" -> "<<(yl+yh)*0.5<<")"<<endl;y = (yl+yh)*0.5;}
		}
		
		res->SetPoint(i,x,(y-minll)); //exp(-(y-minll)
	}
	
	//res->Print();
	
	if (verbose) cout <<"\nMinimum = "<<minll<<" at position "<<minllx<<"  bin = "<<minllbin<<endl;
		
	res->SetName("nllgr"); 
	res->SetLineWidth(2);
	
	int leftbin = minllbin, rightbin = minllbin;
	while (res->GetY()[leftbin]<0.5 && leftbin>0) leftbin--;
	while (res->GetY()[rightbin]<0.5 && rightbin<steps) rightbin++;

	double left = res->GetX()[leftbin], right = res->GetX()[rightbin];
	
	//res->Draw("APL");
	res->Fit("pol2","0q","", left, right);
	
	double fitmin  = res->GetFunction("pol2")->GetMinimum(left, right);
	double fitminx = res->GetFunction("pol2")->GetMinimumX(left, right);
	
	if (verbose) cout <<"Minimum = "<<minll<<" at position "<<minllx<<", fit range = "<<left<<" .. "<<right<<", fitted minimum = "<<fitmin<<" at "<<fitminx<<endl;
	
	xsec->setVal(fitminx);
	xsec->setConstant(true);
	RooFitResult *fitres = 0;
	int fitcnt = 0;

	// and fit it to the combined data set for BG only
	while (++fitcnt<10 && (fitres==0 || fitres->covQual()!=3 ) ) 
	{
		if (fitres) delete fitres;
		fitres = simPdf.fitTo(combData, Save(true), InitialHesse(true), Extended(true)/*, Range("xR1")*/);
	}
		
	
	xsec->setConstant(false);
	xsec->setVal(oldxs);
	
	return res;
}	  	
