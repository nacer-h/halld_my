#ifndef DSelector_pippimkpkm_h
#define DSelector_pippimkpkm_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"
#include "DSelector/DComboTreeHelper.h"

#include "TH1I.h"
#include "TH2I.h"

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

		// Beam Accidentals subtracted: Phi(1020), fo(980), Y(2175)
		TH1F* h_PhiMass_beambunchcut;
		TH1F* h_foMass_beambunchcut;
		TH1F* h_YMass_beambunchcut;

		// 2D
		// Measured
		TH2D* h2_PhiVsfoMass_Measured;
		TH2D* h2_PhiVsYMass_Measured;
		TH2D* h2_foVsYMass_Measured;
		// KinFit
		TH2D* h2_PhiVsfoMass_KinFit;
		TH2D* h2_PhiVsYMass_KinFit;
		TH2D* h2_foVsYMass_KinFit;
		// // after cuts
		// TH2D* h2_PhiVsfoMass_postcut;
		// TH2D* h2_PhiVsYMass_postcut;
		// TH2D* h2_foVsYMass_postcut;

		// cuts
		TH1F* h_TaggerAccidentals;
		TH1F* h_TaggerAccidentals_postcut;
		TH1F* h_mm2_precut;		
		TH1F* h_mm2_postcut;
		TH1F* h_chi2_precut;
		TH1F* h_chi2_postcut;
		TH1F* h_KinFitCL;
		TH2D* h2_DeltaTp;
		TH2D* h2_mm2Vschi2;

		// TH1F* h_KinFitCL; 
		double chi2ndf_min = 1000;

		int ncuts = 20;
		double cutsmin, cutsmax, cutsstep, cutscut;

		// ######## Phi(1020) #############
		TH1F* h_PhiMass_cuts[20];

		// ######## fo(980) #############
		TH1F* h_foMass_cuts[20];

		// ######## Y(2175) #############
		TH1F* h_YMass_cuts[20];

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
