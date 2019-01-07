//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  2 10:11:41 2017 by ROOT version 5.34/32
// from TTree pid/pid
// found on file: treepid.root
//////////////////////////////////////////////////////////

#ifndef piddedx_h
#define piddedx_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class piddedx : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           evnum;
   Int_t           runnum;
   Int_t           tracknum;
   Int_t           nhits;
   Int_t           ntracks;
   Double_t        de[200];   //[nhits]
   Double_t        dx[200];   //[nhits]
   Double_t        dedx[200];   //[nhits]
   Double_t        dedxav;
   Double_t        p;
   Double_t        theta;
   Double_t        costheta;
   Double_t        phi;

   // List of branches
   TBranch        *b_evnum;   //!
   TBranch        *b_runnum;   //!
   TBranch        *b_tracknum;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_de;   //!
   TBranch        *b_dx;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_dedxav;   //!
   TBranch        *b_p;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_costheta;   //!
   TBranch        *b_phi;   //!

   piddedx(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~piddedx() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(piddedx,0);
};

#endif

#ifdef piddedx_cxx
void piddedx::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evnum", &evnum, &b_evnum);
   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("tracknum", &tracknum, &b_tracknum);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("de", de, &b_de);
   fChain->SetBranchAddress("dx", dx, &b_dx);
   fChain->SetBranchAddress("dedx", dedx, &b_dedx);
   fChain->SetBranchAddress("dedxav", &dedxav, &b_dedxav);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("costheta", &costheta, &b_costheta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
}

Bool_t piddedx::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef piddedx_cxx
