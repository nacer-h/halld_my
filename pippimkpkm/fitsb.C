#include "TStyle.h"
#include "TString.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TPad.h"
#include "TList.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "Math/QuantFuncMathCore.h"

const int maxdeg = 6;

// --------------------------------------------------------------------
// find center and difference to half maximum height in both directions
double getCenterInterval(TF1 *f, double &whmlo, double &whmhi)
{
	double xmin, xmax;
	f->GetRange(xmin, xmax);
	
	double fmaxy = f->GetMaximum(xmin,xmax);
	double fmaxx = f->GetMaximumX(xmin,xmax);
	
	whmlo = fmaxx - f->GetX(0.5*fmaxy,xmin, fmaxx);
	whmhi = f->GetX(0.5*fmaxy, fmaxx, xmax) - fmaxx;
	
	return fmaxx;
}

// --------------------------------------------------------------------

double intHist(TH1 *h, double xmin, double xmax)
{
	int bmin = h->FindBin(xmin), bmax = h->FindBin(xmax);
	
	double integ = h->Integral(bmin+1, bmax-1);
	
	integ += fabs(h->GetBinLowEdge(bmin+1)-xmin)/h->GetBinWidth(bmin) * h->GetBinContent(bmin);
	integ += fabs(xmax - h->GetBinLowEdge(bmax))/h->GetBinWidth(bmax) * h->GetBinContent(bmax);
	
	return integ;
}

// -------------------------------------------
// Splits a TString into vector of strings; separation character contained in delim 
std::vector<TString> fitsbSplitString(TString s, TString delim)
{
	std::vector<TString> toks;
	if (s=="" && delim=="") return toks;
		
	while(s.Contains(delim))
	{
		int pos = s.Index(delim);
		TString token = TString(s(0, pos)).Strip(TString::kBoth);
		toks.push_back(token);
		s = s(pos + delim.Length(), s.Length());
	}
	toks.push_back(s.Strip(TString::kBoth));
	
	return toks;
}

// --------------------------------------------------------------------

void savePad(TString name, TString exts="C,png")
{
	std::vector<TString> vext = fitsbSplitString(exts,",");
	for (int i=0;i<(int)vext.size();++i)
		gPad->SaveAs(name+"."+vext[i]); 
}

// --------------------------------------------------------------------

int fitsb(TH1 *h, TString type="", TString parms="", double A=0., double rmin=-999., double rmax = -999., double Qfac=-999)
{
	// ############ INIT AND USAGE #################
	// ############ INIT AND USAGE #################
	// ############ INIT AND USAGE #################
	
	gStyle->SetOptFit(112);
	gStyle->SetOptStat(10);
	int i;
	if (h==0) 
	{
		cout <<"\n### fitsb.C ###: Fits a distribution consisting of signal peak and background; prints out the integral of the signal in the fit region\n\n";
		cout <<"USAGE:\nfitsb(<TH1*/TString h>, [TString type], [TString parms], [double N], [double min], [double max], [double qfac])\n"; 
		cout <<"  <h>     : histogram pointer or name of histogram\n";
		cout <<"  [type]  : default = \"gaus2\"; fit function (dgaus,dgaus,egaus,degaus,bigaus,bw,vgt,novo) and bkg (0=none, 1..6=polynom order, 7=exp, 8=argus, 9=phsp like fcn) E.g.: \"bw3\" = breit wigner + pol3\n";
		cout <<"            gaus = Gaussian; dgaus = Double Gaussian (common mean); bigaus = Bifurcated Gaussian; (d)egaus = Gauss with (2) exp tails; bw = Breit-Wigner; vgt = Voigt-Fcn; novo = Novosibirsk Fcn\n";
		cout <<"  [parms] : string with initial fit parameters w/o amplitude param. E.g.: \"3.5,0.2,1,1,1\" -> mean=3.5, sigma=0.2, a0=1, a1=1, a2=1\n";
		cout <<"            a leading 'f' fixes the parameter. E.g.: \"f3.5,0.2\" -> mean=3.5 (fixed), sigma=0.2\n";
		cout <<"  [A]     : amplitude parameter; default = h.Entries()/3. Negative value fixes the parameter at -A.\n";
		cout <<"  [min]   : minimum of fit range; default = h min x\n";
		cout <<"  [max]   : maximum of fit range; default = h max x\n";
		cout <<"  [qfac]  : defines the width of count/integral window in FWHM around peak; default = 3 (gaus) / 5 (non-gaus)\n\n";
		return 0;
	} 

	//h->Draw();
		
	bool showpull = type.EndsWith("P");
	type.ReplaceAll("P","");

	if (type=="") type="gaus2";
		
	// bin width for amplitude param
	double w = h->GetBinWidth(1);
	
	// fit range; if min=max, fit range over the whole non-zero histogram
	int imin=1, imax=h->GetNbinsX();
	if (rmin==-999.)
	{
		while (h->GetBinContent(imin)==0 && h->GetBinError(imin)==0) imin++; 
		rmin = h->GetBinCenter(imin); //h->GetXaxis()->GetXmin(); 
	}
	
	if (rmax==-999.)
	{
		while (h->GetBinContent(imax)==0 && h->GetBinError(imax)==0) imax--;
		rmax = h->GetBinCenter(imax);//h->GetXaxis()->GetXmax();
	}

	if (rmin>rmax) {double tmp=rmin; rmin=rmax;rmax=tmp;}
	
	//cout <<"rmin="<<rmin<<"  rmax="<<rmax<<endl;
	
	TString tmpl;
	int sigmode = 0;
	int paroff  = 3;  // first parameter offset for the bkg polynomial
	
	
	// ############ SIGNAL FUNCTION #################
	// ############ SIGNAL FUNCTION #################
	// ############ SIGNAL FUNCTION #################

	// ### set signal type
	if      (type.BeginsWith("gaus"))   sigmode = 0;  // gaussian
	else if (type.BeginsWith("bw"))     sigmode = 1;  // breit wigner
	else if (type.BeginsWith("vgt"))    sigmode = 2;  // voigtian
	else if (type.BeginsWith("dgaus"))  sigmode = 3;  // double gaussian
	else if (type.BeginsWith("bigaus")) sigmode = 4;  // bifurcated gaussian
	else if (type.BeginsWith("novo"))   sigmode = 5;  // novosibirsk function
	else if (type.BeginsWith("egaus"))  sigmode = 6;  // exp gaus (similar to crystal ball)
	else if (type.BeginsWith("degaus")) sigmode = 7;  // 2-tailed exp gaus 
	else if (type.BeginsWith("dvgt"))   sigmode = 8;  // double vgt = breit-wigner with double-gauss (same mean) 
	//else if (type.BeginsWith("cball"))  sigmode = 8;  // crystall ball  
	
	if (Qfac==-999 && (sigmode==0 || sigmode==4)) Qfac = 3;
	if (Qfac==-999) Qfac = 5;
	
	// ### set fit function template 
	if      (sigmode==0) tmpl = "%f*[0]/[2]/sqrt(2.0*pi)*exp(-0.5*((x-[1])/[2])^2)";  // Gaussian
//	else if (sigmode==1) tmpl = "%f*[0]*TMath::Voigt(x-[1],0,[2],4.)";                // Breit-Wigner via Voigt
	else if (sigmode==1) tmpl = "%f*[0]*[2]/(2.0*pi*((x-[1])^2+([2]/2)^2))";          // Breit-Wigner
	else if (sigmode==2) tmpl = "%f*[0]*TMath::Voigt(x-[1],[2],[3],4.)";              // Voigt
	else if (sigmode==3) tmpl = "%f*[0]*([4]/([4]+1.0)/[2]*exp(-0.5*((x-[1])/[2])^2)+1.0/([4]+1.0)/[3]*exp(-0.5*((x-[1])/[3])^2))/sqrt(2.0*pi)";  // 2 Gaussians with common mean (double gauss)
	else if (sigmode==4) tmpl = "%f*[0]*2./([2]+[3])/sqrt(2.0*pi)*exp(-0.5*(x-[1])^2/([2]*(x<[1])+[3]*(x>=[1]))^2)";                              // bifurcated Gaussian
	else if (sigmode==5) tmpl = "[0]*(([3]>0)*((x-[1])*[3]/[2]<1)*exp( -0.5*(log(1-((x-[1])*[3]/[2]<1)*(x-[1])*[3]/[2]))^2/(2./2.3548*sinh(0.5*[3]*2.3548))^2 - 0.5*[2]^2)+([3]<=0)*(1./[2]/sqrt(2.0*pi)*exp(-0.5*((x-[1])/[2])^2)))";            // Novisibirsk function
	else if (sigmode==6) tmpl = "[0]*(((x-[1])/[2]>=-[3])*TMath::Gaus(x,[1],[2])  +  ((x-[1])/[2]<-[3])*(exp(0.5*[3]^2 + [3]*(x-[1])/[2])))";  // exp gaus 
	else if (sigmode==7) tmpl = "[0]*((((x-[1])/[2])<=-[3])*(exp(0.5*[3]^2+[3]*((x-[1])/[2]))) + (((x-[1])/[2])>-[3] && ((x-[1])/[2])<=[4])*(exp(-0.5*(x-[1])^2/[2]^2)) + (((x-[1])/[2])>[4])*(exp(0.5*[4]^2-[4]*((x-[1])/[2]))))";  // 2-tailed exp gaus 
	else if (sigmode==8) tmpl = "%f*[0]*(([5]/([5]+1.0)*TMath::Voigt(x-[1],[2],[4],4.)+(1.0/([5]+1.0)*TMath::Voigt(x-[1],[3],[4],4.))";              // Double Voigt
	//else if (sigmode==8) tmpl = "[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])"; // crystal ball fcn, pars: 1 = alpha, 2 = n, 3 = sigma, 4 = mu
		
    if (sigmode==2 || sigmode==4 || sigmode==5 || sigmode==6) paroff = 4;
	if (sigmode==3|| sigmode==7 ) paroff = 5;
	if (sigmode==8) paroff = 6;

	// ############ BACKGROUND FUNCTION #################
	// ############ BACKGROUND FUNCTION #################
	// ############ BACKGROUND FUNCTION #################

	// ### degree of bkg polynomial; 8+9 are argus fcn and truncated polynom
	int deg = -1; // no background fcn
	TString degs = type(type.Length()-1,1);  // last character of fcn string could be the digit for bkg polynomial
	TString nums="0123456789";
	if (nums.Contains(degs)) deg = degs.Atoi();
	if (deg>maxdeg && deg<7) deg=maxdeg;
	
	// ### total number of params
	int partot = paroff+deg+1;
	
	//cout <<partot<<"  "<<paroff<<"  "<<deg<<endl;

	if (deg==7) partot-=6;    // correction for number of parameters for abusing deg==7 (->8 parms) for expo with only 2 parameters
	if (deg==8) partot-=5;    // correction for number of parameters for abusing deg==8 (->9 parms) for argus fcn with only 4 parameters
	if (deg==9) partot-=5;    // correction for number of parameters for abusing deg==9 (->10 parms) for f(x) = A * (x>a) * (x<b) * |x-a|^c * |b-x|^d with only 5 parameters
	
	// ### set fit function TFormula
	TString fnc  = TString::Format(tmpl.Data(), w); // will get total fit function
	TString fncs = TString::Format(tmpl.Data(), w); // only signal part
	
	if (deg>=0&&deg<7) fnc+=TString::Format("+pol%d(%d)",deg, paroff);
	else if (deg==7) // ### exponetial
	{
		fnc+=TString::Format("+expo(%d)", paroff);
	}
	else if (deg==8) // ### ARGUS fcn bkg
	{
		fnc+=TString::Format("+(x<[%d])*[%d]*x*(abs(1-(x/[%d])^2))^[%d]*exp([%d]*(1-(x/[%d])^2))", paroff+1, paroff, paroff+1, paroff+2, paroff+3, paroff+1);
	}
	else if (deg==9) // ### f(x) = A * (x>a) * (x<b) * |x-a|^c * |b-x|^d
	{
		fnc+=TString::Format("+[%d]*(x>[%d])*(x<[%d])*abs(x-[%d])^[%d]*abs([%d]-x)^[%d]", paroff, paroff+1, paroff+2, paroff+1, paroff+3, paroff+2, paroff+4);
	}

	cout<<"\nFit fcn: " <<fnc<<endl<<endl;


	// ############ INITIAL PARAMETER SETTINGS #################
	// ############ INITIAL PARAMETER SETTINGS #################
	// ############ INITIAL PARAMETER SETTINGS #################

	TSpectrum spec(1);
	spec.Search(h, 2, "goff", 0.8);

	// ### set the parameter list from par string
	double prm[20]={0};       // the parameters
	double prmr[20]={0};      // possible second par values defining a range
	bool   fixparm[20]={0};   // flags for fixed parameter
	int pcnt = 0;
	if (parms!="") parms=parms+",";
	
	// ### set some defaults
	prm[0] = spec.GetNPeaks()>0 ? spec.GetPositionY()[0] : h->GetEntries()/3.;                   // N
	prm[1] = spec.GetNPeaks()>0 ? spec.GetPositionX()[0] : h->GetBinCenter(h->GetMaximumBin());  // mean estimate
	//prm[0] = h->GetEntries()/3.;                   // N
	//prm[1] = h->GetBinCenter(h->GetMaximumBin());  // mean estimate
	prm[2] = h->GetRMS()/2;                                                                      // sigma/Gamma estimate
	
	// ### perform simple test fit gaus + pol1 to improve parameters
	TF1 ftmp("ftmp","gaus(0)+pol2(3)");
	ftmp.SetParameters(prm[0], prm[1], prm[2],1,1,1);
	h->Fit("ftmp","","qn", prm[1]-2*prm[2], prm[1]+2*prm[2]);
	prm[2] = ftmp.GetParameter(2); 
	
	// ### set estimates for other parameters
	if (sigmode>1)  prm[3] = prm[2]*2.;            // 2nd sigma/Gamma estimate
	if (sigmode==3) prm[4] = 1.0;                  // ratio sigma1/sigma2	
	if (sigmode==6) prm[3] = 1.0;                  // k_L (tail parameter)
	if (sigmode==7) { prm[3] = 1.0; prm[4] = 1.0;} // k_H, k_L (tail parameters)
	if (sigmode==8) prm[4] = prm[3];               // gamma for double vgt
	
	if (deg>=0 && deg<maxdeg+1) 
	{
		for (i=paroff;i<partot;++i) prm[i] = 1.0;      // bkg parms
	}
	else if (deg==7) //expo
	{
		prm[paroff]   = 1;//h->GetBinContent(1)+h->GetBinContent(h->GetNbinsX());
		prm[paroff+1] = (h->GetBinContent(h->GetNbinsX())-h->GetBinContent(1))/fabs((h->GetBinContent(h->GetNbinsX())-h->GetBinContent(1)));
	}
	else if (deg==8)  //argus
	{
		prm[paroff]   = h->GetMaximum()/5.;
		prm[paroff+1] = rmax;
		prm[paroff+2] = prm[paroff+3] = 0.5;
	}
	else if (deg==9)  //phsp bg
	{
		prm[paroff]   = h->GetMaximum()/5.;
		prm[paroff+1] = rmin;
		prm[paroff+2] = rmax;
		prm[paroff+3] = prm[paroff+4] = 0.5;
	}
	
	// set initial amplitude parameter
	if (A>0) prm[0] = A;
	// if A<0 the parameter value will be fixed to -A
	if (A<0) {prm[0] = -A; fixparm[0]=true;}
		
		
	// ############ EVALUATE PRESET PARAMETERS #################
	// ############ EVALUATE PRESET PARAMETERS #################
	// ############ EVALUATE PRESET PARAMETERS #################
	
	// ### analyse parameter string
	while (parms!="")
	{
		TString sparm = parms(0,parms.Index(","));
		if (sparm!="") 
		{
			// fix this parameter?
			if (sparm.BeginsWith("f")) 
			{
				sparm=sparm(1,1000); // cut away the 'f'
				fixparm[1+pcnt] = true;
			}
			// ranged parameter in form 'r10.2|10.5'?
			if (sparm.BeginsWith("r"))
			{
				TString sparm2 = sparm(sparm.Index("|")+1,1000);
				sparm = sparm(1,sparm.Index("|"));
				prmr[1+pcnt]=sparm2.Atof();
			}
			prm[1+pcnt++] = sparm.Atof();
			parms = parms(parms.Index(",")+1,1000);
		}
	} 
	
	
	// ############ SETUP TOTAL FIT FUNCTION #################
	// ############ SETUP TOTAL FIT FUNCTION #################
	// ############ SETUP TOTAL FIT FUNCTION #################

	// setup the total fit fcn (signal + background)
	TF1 ff1("ff1", fnc);
	ff1.SetLineColor(4); ff1.SetNpx(500);  // some style setting

	// set the parameter names
	if (sigmode==0)      ff1.SetParNames("N","#mu","#sigma",                                "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==1) ff1.SetParNames("N","#mu","#Gamma",                                "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==2) ff1.SetParNames("N","#mu","#sigma","#Gamma",                       "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==3) ff1.SetParNames("N","#mu","#sigma_{1}","#sigma_{2}","R",           "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==4) ff1.SetParNames("N","#mu","#sigma_{1}","#sigma_{2}",               "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==5) ff1.SetParNames("A","#mu","#sigma","#tau",                         "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==6) ff1.SetParNames("A","#mu","#sigma","k",                            "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==7) ff1.SetParNames("A","#mu","#sigma","k_{L}","k_{H}",                "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}","a_{5}");
	else if (sigmode==8) {ff1.SetParNames("A","#mu","#sigma_{1}","#sigma_{2}","#Gamma", "R", "a_{0}","a_{1}","a_{2}","a_{3}","a_{4}"); ff1.SetParName(10,"a_{5}");}
	
	if (deg==6) ff1.SetParName(paroff+6, "a_{6}");
	
	if (deg==7) 
	{
		ff1.SetParName(paroff  , "a_{exp}");
		ff1.SetParName(paroff+1, "#lambda");
	}
	
	if (deg==8) 
	{
		ff1.SetParName(paroff  , "a");
		ff1.SetParName(paroff+1, "m_{0}");
		ff1.SetParName(paroff+2, "p");
		ff1.SetParName(paroff+3, "c");
	}
	if (deg==9) 
	{
		ff1.SetParName(paroff  , "A");
		ff1.SetParName(paroff+1, "a");
		ff1.SetParName(paroff+2, "b");
		ff1.SetParName(paroff+3, "c");
		ff1.SetParName(paroff+4, "d");
	}
			
	// ### fix the parameters which were requested to be fixed
	for (i=0;i<partot;++i)
	{
		if (fixparm[i])
			ff1.FixParameter(i,prm[i]);
		else if (prmr[i]!=0)
			ff1.SetParLimits(i,prm[i],prmr[i]);
		else
			ff1.SetParameter(i,prm[i]);
	}
	
	// ############ FIT #################
	// ############ FIT #################
	// ############ FIT #################
	
	h->Fit("ff1","q","goff",rmin, rmax);
	h->Fit("ff1","m","goff",rmin, rmax);
		
	// ############ PLOT #################
	// ############ PLOT #################
	// ############ PLOT #################
	
	// ### Do we plot pull window underneath ? 
	if (showpull)
	{
		TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
		TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.251);
		
		pad1->SetBottomMargin(0.0);
		pad1->SetBorderMode(0);
		pad2->SetTopMargin(0.0);
		pad2->SetBottomMargin(0.4);
		pad2->SetBorderMode(0);

		pad1->Draw();
		pad2->Draw();
		
		TH1F *hres = new TH1F("hres","", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
		hres->SetStats(0);
		hres->GetXaxis()->SetLabelSize(0.14);
		hres->GetXaxis()->SetTitleSize(0.14);
		hres->GetXaxis()->SetTitleOffset(1.1);
		hres->GetXaxis()->SetTickLength(0.05);

		hres->GetYaxis()->SetLabelSize(0.08);
		hres->GetYaxis()->SetTitleSize(0.08);
		hres->GetYaxis()->SetTitleOffset(0.35);
		hres->GetYaxis()->SetNdivisions(5);
		
		hres->SetXTitle(h->GetXaxis()->GetTitle());
		hres->SetYTitle("pull [#sigma] ");
		
		for (int i=1;i<=h->GetNbinsX();++i)
		{
			if (h->GetBinError(i)!=0)
			{
				hres->SetBinContent(i, (h->GetBinContent(i) - ff1.Eval(h->GetBinCenter(i)))/h->GetBinError(i));
				hres->SetBinError(i, 1.0);
			}
		}
		hres->GetYaxis()->SetRangeUser(-5,5);
		hres->GetYaxis()->SetNdivisions(505);

		pad2->cd();
		hres->Draw("e");
		pad2->SetGridy();
		
		TLine l;
		l.SetLineColor(4);
		l.DrawLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
		
		pad1->cd();
		
		// pull up the minimum so that a minimum 0 on the y-axis is not cut
				
		double ymax = 0;
		// this is the maximum for (bin content + bin error)
		for (int i=1;i<=h->GetNbinsX(); ++i) if ((h->GetBinContent(i)+h->GetBinError(i)) > ymax) ymax = (h->GetBinContent(i)+h->GetBinError(i));
		double ymin = -h->GetYaxis()->GetLabelSize()*0.5*ymax;
		
		if (gPad->GetFrame()->GetY1()>ymin) 
		{
			//h->SetMinimum(-0.05*(h->GetYaxis()->GetXmax()-h->GetXaxis()->GetXmin()));
			h->GetYaxis()->SetRangeUser(ymin, ymax*1.05);
		}
	}
	
	// ### Draw main histogram
	
	h->Draw("e");	
	

	// ############ SETUP PURE BACKGROUND FUNCTION #################
	// ############ SETUP PURE BACKGROUND FUNCTION #################
	// ############ SETUP PURE BACKGROUND FUNCTION #################
	
	// ### compute integral of total fit fcn in interval min ... max
	
	TF1 *ff2=0;
	
	// ### do we have a bkg polynomial? 
	if (deg>=0 && deg<=maxdeg)
	{
		// create pure bkg function
		ff2=new TF1("ff2",TString::Format("pol%d",deg));		
		
		// copy parameters from full fcn
		for (i=0;i<deg+1;++i) ff2->SetParameter(i,ff1.GetParameter(paroff+i));
	}
	
	// ### expo bkg
	if (deg==7)
	{
		// create pure bkg function
		ff2=new TF1("ff2","expo(0)");		
		
		// copy parameters from full fcn
		for (i=0;i<2;++i) ff2->SetParameter(i,ff1.GetParameter(paroff+i));
	 }
	// ### argus bkg
	if (deg==8)
	{
		// create pure bkg function
		ff2=new TF1("ff2","(x<[1])*[0]*x*(abs(1-(x/[1])^2))^[2]*exp([3]*(1-(x/[1])^2))");
		
		// copy parameters from full fcn
		for (i=0;i<4;++i) ff2->SetParameter(i,ff1.GetParameter(paroff+i));		
	}
	// ### phsp background 
	if (deg==9)
	{
		// create pure bkg function
		ff2=new TF1("ff2","[0]*(x>[1])*(x<[2])*abs(x-[1])^[3]*abs([2]-x)^[4]");
		
		// copy parameters from full fcn
		for (i=0;i<5;++i) ff2->SetParameter(i,ff1.GetParameter(paroff+i));
	}

	// ### common to all
	// ### add it to the histogram (so that it is shown when drawing the histogram)
	if (ff2) 
	{
		ff2->SetRange(rmin,rmax);
		ff2->SetLineColor(kRed+1); ff2->SetLineStyle(2); ff2->SetNpx(500);  // some style setting
		ff2->SetRange(rmin,rmax);
		h->GetListOfFunctions()->Add(ff2);
	}

	// ############ INTEGRALS / S:N / SIGNIFICANCES #################
	// ############ INTEGRALS / S:N / SIGNIFICANCES #################
	// ############ INTEGRALS / S:N / SIGNIFICANCES #################

	// ### compute bkg intergral
	int integral = ff1.Integral(rmin,rmax)/w;
	int bkg = ff2 ? ff2->Integral(rmin,rmax)/w : 0;
	
	// #signals = #total - #bkg 
	if (bkg>=0) integral -= bkg;
	
	// ### restyle the stats box to adapt different number of parameters
	TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
	if (s)
	{
		s->SetX1NDC(0.77);
		s->SetY1NDC(gStyle->GetStatX()-0.028*(partot+3));
		s->SetTextSize(0.025); 
	}
	
	
	// ### compute S/N and significance
	double StoN = (double)integral/(double)bkg;
	double Sign = (double)integral/sqrt((double)bkg+(double)integral);
	
	// ### print out the var list string for easier copy paste
	cout <<endl<<"Var-list: \"";
	for (int i=1;i<partot;++i)
	{
		if (fixparm[i]) cout<<"f";
		cout << ff1.GetParameter(i);
		if (i<partot-1) cout<<",";
	}
	cout <<"\"";
	
	// ### print some info about S and B
	printf("\n\n[%.3f ... %.3f] : S_int = %7.0f    B_int = %7.0f    S/B = %6.2f    S/sqrt(S+B) = %7.3f (full window)\n\n", rmin, rmax, (double)integral, (double)bkg, StoN, Sign); 


	// ############ QUALITY MEASURES IN ROI #################
	// ############ QUALITY MEASURES IN ROI #################
	// ############ QUALITY MEASURES IN ROI #################

	// ### compute some other intervals	
	TF1 ff1s("ff1s", fncs);
	ff1s.SetParameters(ff1.GetParameters());
	ff1s.SetRange(rmin,rmax);

	// ### get limits and centre of interval for signal FWHM
	double wlo, whi;
	double cntr = getCenterInterval(&ff1s, wlo, whi);

	// ### define ROI = Qfac * FWHM	
	double wloq = cntr-Qfac*wlo, whiq = cntr+Qfac*whi;
	
	// ### integral of signal in ROI
	double signalq = ff1s.Integral(wloq,whiq)/w;
	
	// ### if used, integral of background in ROI
	double bkgq=0, StoNq=0; 
	if (ff2)
	{
		bkgq = ff2->Integral(wloq,whiq)/w;
		
		// ### signal to noise in ROI
		StoNq = signalq/bkgq; 
	}
	
	// ### significance for integrated signal
	double Signq = signalq/sqrt(signalq+bkgq);
	
	// ### counted signal = total counts in ROI minus integrated BG
	double signalqcnt = intHist(h,wloq,whiq)-bkgq;
	
	// ### print quality for ROI S_integrated and S_counted
	printf("[%.3f ... %.3f] : S_int = %7.0f    B_int = %7.0f    S/B = %6.2f    S/sqrt(S+B) = %7.3f (%.1f FWHM window, dotted lines)\n",   wloq, whiq, signalq, bkgq, StoNq, Signq, Qfac); 
	printf("[%.3f ... %.3f] : S_cnt = %7.0f    B_int = %7.0f    S/B = %6.2f    S/sqrt(S+B) = %7.3f (%.1f FWHM window, dotted lines, counted)\n\n", wloq, whiq, signalqcnt, bkgq,signalqcnt/bkgq , Signq = signalqcnt/sqrt(signalqcnt+bkgq), Qfac); 
	
	
	// ############ DRAW LINES AND PRINT NUMBERS IN PLOT #################
	// ############ DRAW LINES AND PRINT NUMBERS IN PLOT #################
	// ############ DRAW LINES AND PRINT NUMBERS IN PLOT #################
	
	TLine ll;
	ll.SetLineColor(14);
	ll.SetLineStyle(2);
	
	// ### draw ROI
	ll.DrawLine(wloq, 0, wloq, ff1.Eval(wloq));
	ll.DrawLine(whiq, 0, whiq, ff1.Eval(whiq));
	
	// ### write numbers of S and B to the histogram	
	TLatex ls,lb;
	ls.SetTextFont(42);
	lb.SetTextFont(42);
	ls.SetNDC(); lb.SetNDC();

	// ### margin infos for writing numbers on pad
	double tm = 1.0-gStyle->GetPadTopMargin();
	double lm = gStyle->GetPadLeftMargin();
	
	double xmin = h->GetXaxis()->GetXmin(), xmax = h->GetXaxis()->GetXmax();
	
	ls.SetTextSize(0.045); ls.SetTextColor(4);
	ls.DrawLatex(lm+0.02,tm-0.06,TString::Format("S_{int} = %d",(int)signalq));

	ls.SetTextColor(14);
	ls.DrawLatex(lm+0.02, tm-0.12,TString::Format("S_{cnt} = %d",(int)(signalqcnt)));

	if (deg>=0)
	{
		lb.SetTextSize(0.045); lb.SetTextColor(kRed+1);
		lb.DrawLatex(lm+0.02, tm-0.18,TString::Format("B_{int} = %d",(int)bkgq));
	}

	gPad->GetCanvas()->cd();

	return integral;
}

// --------------------------------------------------------------------
// accepting a histogram name instead of pointer to TH1F
// --------------------------------------------------------------------
int fitsb(TString hname="", TString type="", TString parms="", double A=0., double rmin=-999., double rmax = -999., double Qfac=-999)
{
	TH1 *h=0;
	if (hname!="") h = (TH1*)gDirectory->Get(hname);
	
	return fitsb(h, type, parms, A, rmin, rmax, Qfac);
}

