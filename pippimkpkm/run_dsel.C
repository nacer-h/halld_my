#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TChain.h"

#include <iostream>

TString datapath = "/lustre/nyx/panda/kgoetzen/gluex/data/";
TString outpath  = "/lustre/hebe/panda/ahamdi//gluex_root_analysis/workdir/dataout/";

void run_dsel(TString selname="", TString datasubdir="", TString treename="", TString filename="", TString outprefix="", int nev=1e9, int firstevt=0) 
{
	if (selname=="" || (datasubdir=="" && treename=="")  || (treename=="" && filename==""))
	{
		cout <<"USAGE:\nroot -l -b -q run_dsel.C+(TString sel, TString sub, TString tree, TString file, TString out, <int N=1e9>, <int N0=0>)\n\n"; 
		cout <<"  sel    : Name of DSelector w/o '.C'.\n";
		cout <<"  sub    : sub directory relative to "<<datapath<<"; can be blank if same as given tree name\n";
		cout <<"  tree   : Name of tree w/o '_Tree'\n";
		cout <<"  file   : input file; either a single file name or just one or more run numbers as string (e.g. \"031002\", or \"031002:031006:031018\")\n";
		cout <<"  out    : output file prefix (file name prefix to be stored in output directory specified in variable 'outpath')\n";
		cout <<"  <N>    : optional number of events \n";
		cout <<"  <N0>   : optional number of first event \n\n";
		return;
	}
	
	TStopwatch timer;
	timer.Start();
	
	gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
	
	TChain *ntp=0;
	
	// is filename a (collection of) run number(s)?
	if (filename.Atoi()>0 || filename.Contains(":"))
	{
		// single run? append ':'
		if (!filename.EndsWith(":")) filename=filename+":";
		
		// since we only got run numbers, treename must have been specified
		ntp = new TChain(treename+"_Tree");
		if (datasubdir=="") datasubdir=treename;
		
		while (filename.Contains(":"))
		{
			int curr_run = (TString(filename(0,filename.Index(":"))).Atoi());
			if (curr_run)
			{
				TString ff = Form("%s/%s/tree_%s_%06d.root",datapath.Data(),datasubdir.Data(),treename.Data(),curr_run);
				int ok = ntp->Add(ff,-1);
				if (ok) cout <<"added "<<ff<<endl;
				
			}
			filename = filename(filename.Index(":")+1,1000);
		}
	}
	// or do we habe explicit filename?
	else
	{
		TString ff = datapath+datasubdir+"/"+filename;
		
		if (treename=="")
		{
			// open file
			TFile *f = new TFile(ff);
			// retrieve tree name
			treename = gDirectory->GetListOfKeys()->At(0)->GetName();
			cout <<"Found tree '"<<treename<<"'"<<endl;
			treename.ReplaceAll("_Tree","");	
			f->Close();
		}
		// create TChain and append file
		ntp = new TChain(treename+"_Tree");
		
		int ok = ntp->Add(ff,-1);
		if (ok) cout <<"added "<<ff<<" with "<<nev<<" events"<<endl;
	}
		
	if (outprefix=="") outprefix = treename;
	
	ntp->Process(selname+".C+",outpath+outprefix,nev,firstevt);
	
	timer.Stop();
	cout << "RealTime: " << timer.RealTime() << " seconds" << endl;
	cout << "CpuTime: " << timer.CpuTime() << " seconds" << endl;
}
