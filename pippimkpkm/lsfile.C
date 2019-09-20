#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <iostream>

unsigned long lsfile(TString fname="")
{
	if (fname=="") { cout <<"USAGE:\nroot -l -b -q -n 'lsfile.C(filename)'\n"; return 0;}
	TFile f(fname);
	TString treename = gDirectory->GetListOfKeys()->At(0)->GetName();
	
	TTree *t=(TTree*)f.Get(treename);
	unsigned long N = t->GetEntries();
	
	printf("%s:%lu\n",treename.Data(),N);
	return(N);
}
