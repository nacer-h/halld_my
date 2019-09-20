#ifndef DSelector_pippimkpkm_h
#define DSelector_pippimkpkm_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"
#include "DSelector/DComboTreeHelper.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TString.h"
#include "TObjString.h"

class DSelector_pippimkpkm : public DSelector
{
	public:

		DSelector_pippimkpkm(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_pippimkpkm(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dPiPlusWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dKMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

		/* // Monte Carlo Truth masses: Phi(1020), fo(980), Y(2175) */
		/* TH1F* h_PhiMass_Thrown; */
		/* TH1F* h_foMass_Thrown; */
 		/* TH1F* h_YMass_Thrown; */

		// Measured masses: Phi(1020), fo(980), Y(2175)
		TH1F* h_PhiMass_Measured;
		TH1F* h_foMass_Measured;
 		TH1F* h_YMass_Measured;		
	       		
		// P4 & Vertex kinematic fit masses: Phi(1020), fo(980), Y(2175)
		TH1F* h_PhiMass_KinFit;
		TH1F* h_foMass_KinFit;
 		TH1F* h_YMass_KinFit;

		// 2D
		TH2D* h2_PhiVsfoMass_KinFit;
		TH2D *h2_PhiVsYMass_KinFit;
		TH2D *h2_foVsYMass_KinFit;

		// Beam Accidentals subtracted: Phi(1020), fo(980), Y(2175)
		TH1F* h_PhiMass_beambunchcut;
		TH1F* h_foMass_beambunchcut;
		TH1F* h_YMass_beambunchcut;

		// cuts
		TH1F* h_TaggerAccidentals;
		TH1F* h_TaggerAccidentals_postcut;
		TH1F *h_mm2;
		TH1F *h_chi2;
		TH1F *h_KinFitCL;

		// TH1F* h_KinFitCL;
		// double chi2ndf_min = 1000;

		// int ncuts;
		// double cutsmin, cutsmax, cutsstep, cutscut;

		TLorentzVector locP4_p;
		TLorentzVector locX4_p;
		TLorentzVector locP4_pip;
		TLorentzVector locX4_pip;
		TLorentzVector locP4_kp;
		TLorentzVector locX4_kp;				

		// ######## Phi(1020) #############
		TH1F *h_PhiMass_postcuts;
		// TH1F *h_PhiMass_cuts[20];

		// ######## fo(980) #############
		TH1F *h_foMass_postcuts;
		// TH1F *h_foMass_cuts[20];

		// ######## Y(2175) #############
		TH1F *h_YMass_postcuts;
		// TH1F *h_YMass_cuts[20];

        // *************** bggen ***************
		TH1F* dHistThrownTopologies;
		map<TString, TH1I*> dHistInvariantMass_ThrownTopology;

		// TOOL FOR FLAT TREE OUTPUT
		DComboTreeHelper *dComboTreeHelper;

	ClassDef(DSelector_pippimkpkm, 0);
};

void DSelector_pippimkpkm::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dPiMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(4));
}

#endif // DSelector_pippimkpkm_h
