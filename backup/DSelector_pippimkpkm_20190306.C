#include "DSelector_pippimkpkm.h"
#include "TLorentzRotation.h"

void DSelector_pippimkpkm::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
    TString option = GetOption();
    if (option=="") option = "pippimkpkm";

	dOutputFileName = ""; //option + ".root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = option + "_flat.root"; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = "ntp"; //if blank, default name will be chosen

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	// PIDs
	deque<Particle_t> locPhiPIDs; locPhiPIDs.push_back(KMinus);	locPhiPIDs.push_back(KPlus);
	deque<Particle_t> locf0PIDs; locf0PIDs.push_back(PiMinus); locf0PIDs.push_back(PiPlus);
	deque<Particle_t> locYPIDs; locYPIDs.push_back(KMinus); locYPIDs.push_back(KPlus); locYPIDs.push_back(PiMinus); locYPIDs.push_back(PiPlus);

	//PID
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	// dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, "meas_precut"));
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "precut"));

	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiPlus, SYS_BCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiPlus, SYS_TOF));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiPlus, SYS_FCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, PiMinus, SYS_BCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiMinus, SYS_TOF));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, PiMinus, SYS_FCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.6, KPlus, SYS_BCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, KPlus, SYS_TOF));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.375, KPlus, SYS_FCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.6, KMinus, SYS_BCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.3, KMinus, SYS_TOF));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 2.375, KMinus, SYS_FCAL));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.0, Proton, SYS_BCAL));
	dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.48, Proton, SYS_TOF));
	// dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 1.9, Proton, SYS_FCAL));
	// dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, "meas_postcut"));     

	// Invariant mass plots
	//Measured
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, locPhiPIDs, 600, 1, 1.2, "Phi_meas_PreCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, locf0PIDs, 600, 0.3, 1.2, "fo_meas_PreCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, locYPIDs, 600, 1.6, 3.2, "y_meas_PreCut"));
	//KinFit
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locPhiPIDs, 600, 1, 1.2, "Phi_PreCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locf0PIDs, 600, 0.3, 1.2, "fo_PreCut"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, locYPIDs, 600, 1.6, 3.2, "y_PreCut"));

	// MISSING MASS
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 600, -0.1, 0.1, "mm2_meas"));

	// Cuts
	// dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 60, "CLcut")
    // dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.03));
	// dAnalysisActions.push_back(new DCutAction_InvariantMass(dComboWrapper, true, locPhiPIDs, 1.005, 1.035, "kpkmcut")

	// PID after cut
	// dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, true, "postcut"));


	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper, "KinFitResults"));
	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions(dAnalysisActions, dComboWrapper, true, 0, locPhiPIDs, 600, 1, 1.2, "CutActionEffect");

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	// // Thrown  masses: Phi(1020), fo(980), Y(2175)
	// h_PhiMass_Thrown = new TH1F("PhiMass_Thrown", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 1, 1.2);
	// h_foMass_Thrown = new TH1F("foMass_Thrown", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
	// h_YMass_Thrown = new TH1F("YMass_Thrown", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.6, 3.2);

	// Measured masses: Phi(1020), fo(980), Y(2175)
	h_PhiMass_Measured = new TH1F("PhiMass_Measured", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 1, 1.2);
	h_foMass_Measured = new TH1F("foMass_Measured", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
	h_YMass_Measured = new TH1F("YMass_Measured", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.6, 3.2);

	// P4 & Vertex kinematic fit masses: Phi(1020), fo(980), Y(2175)
	h_PhiMass_KinFit = new TH1F("PhiMass_KinFit", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 1, 1.2);
	h_foMass_KinFit = new TH1F("foMass_KinFit", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
	h_YMass_KinFit = new TH1F("YMass_KinFit", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.6, 3.2);

	// Beam Accidentals subtracted: Phi(1020), fo(980), Y(2175)
	h_PhiMass_beambunchcut = new TH1F("PhiMass_beambunchcut", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 1, 1.2);
	h_foMass_beambunchcut = new TH1F("foMass_beambunchcut", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
	h_YMass_beambunchcut = new TH1F("YMass_beambunchcut", ";m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.6, 3.2);

     // post Phi cut: Phi(1020), fo(980), Y(2175)
	 h_PhiMass_postcuts = new TH1F("h_PhiMass_postcuts", ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 1.005, 1.035);
     h_foMass_postcuts = new TH1F("h_foMass_postcuts", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
     h_YMass_postcuts = new TH1F("h_YMass_postcuts", ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.6, 3.2); //2.5, 3.2, 2.05, 2.35

	// 2D
	// KinFit
	h2_PhiVsfoMass_KinFit = new TH2D("h2_PhiVsfoMass_KinFit", ";m_{K^{+}K^{-}} [GeV/c^{2}];m_{#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 1, 1.2, 600, 0.3, 1.2);
	h2_PhiVsYMass_KinFit = new TH2D("h2_PhiVsYMass_KinFit", ";m_{K^{+}K^{-}} [GeV/c^{2}];m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 1, 1.2, 600, 1.6, 3.2);
	h2_foVsYMass_KinFit = new TH2D("h2_foVsYMass_KinFit", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 0.3, 1.2, 600, 1.6, 3.2);
	// //after cuts
	// h2_PhiVsfoMass_postcut = new TH2D("h2_PhiVsfoMass_postcut", ";m_{K^{+}K^{-}} [GeV/c^{2}];m_{#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 1, 1.04, 600, 0.3, 1.2);
	// h2_PhiVsYMass_postcut = new TH2D("h2_PhiVsYMass_postcut", ";m_{K^{+}K^{-}} [GeV/c^{2}];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 1, 1.04, 600, 1.6, 3.2);
	// h2_foVsYMass_postcut = new TH2D("h2_foVsYMass_postcut", ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}]", 600, 0.3, 1.04, 600, 1.6, 3.2);

	// cuts
	h_TaggerAccidentals = new TH1F("h_TaggerAccidentals", "Vertex time - RF (ns);#Deltat_{Beam-RF}(ns)", 300, -18, 18);
	h_TaggerAccidentals_postcut = new TH1F("h_TaggerAccidentals_postcut", "Vertex time - RF (ns);#Deltat_{Beam-RF}(ns)", 300, -18, 18);
	h_mm2_precut = new TH1F("h_mm2_precut", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.1, 0.1);
	h_mm2_postcut = new TH1F("h_mm2_postcut", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.03, 0.03);
	h_chi2_precut = new TH1F("h_chi2_precut", "#chi^{2}/ndf of the Kinematic Fit (P4 + Vertex) ;#chi^{2};Counts", 600, 0, 60);
	h_chi2_postcut = new TH1F("h_chi2_postcut", "#chi^{2}/ndf of the Kinematic Fit (P4 + Vertex) ;#chi^{2};Counts", 600, 0, 60);
	h_KinFitCL = new TH1F("KinFitCL", ";Kinematic Fit Confidence Level", 100, 0., 1.);
	h2_mm2Vschi2 = new TH2D("h2_mm2Vschi2", ";#chi^{2} of the Kinematic Fit (P4 + Vertex);Missing Mass Squared (GeV/c^{2})^{2}", 600, 0, 60, 600, -0.03, 0.03);
	// h2_DeltaT_p_tof = new TH2D("h2_DeltaT_p_tof", ";p (GeV/c); TOF #Delta T (ns) (proton)", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_p_bcal = new TH2D("h2_DeltaT_p_bcal", ";p (GeV/c); BCAL #Delta T (ns) (proton)", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_p_fcal = new TH2D("h2_DeltaT_p_fcal", ";p (GeV/c); FCAL #Delta T (ns) (proton)", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_pip_tof = new TH2D("h2_DeltaT_pip_tof", ";p (GeV/c); TOF #Delta T (ns) (#pi^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_pip_bcal = new TH2D("h2_DeltaT_pip_bcal", ";p (GeV/c); BCAL #Delta T (ns) (#pi^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_pip_fcal = new TH2D("h2_DeltaT_pip_fcal", ";p (GeV/c); FCAL #Delta T (ns) (#pi^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_kp_tof = new TH2D("h2_DeltaT_kp_tof", ";p (GeV/c); TOF #Delta T (ns) (K^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_kp_bcal = new TH2D("h2_DeltaT_kp_bcal", ";p (GeV/c); BCAL #Delta T (ns) (K^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);
	// h2_DeltaT_kp_fcal = new TH2D("h2_DeltaT_kp_fcal", ";p (GeV/c); FCAL #Delta T (ns) (K^{+})", 250, 0.0, 10.0, 500, -10.0, 10.0);

	// cutsmin = h_mm2_postcut->GetXaxis()->GetBinLowEdge(1);
	// cutsmax = h_mm2_postcut->GetXaxis()->GetBinUpEdge(600);
	// ncuts = 20;
	// cutsstep = cutsmax / ncuts; //cutsmin

	// ######## Phi(1020) #############
	// for (int icuts = 0; icuts < ncuts; ++icuts)
	// {
	// 	TString hname_PhiMass_cuts = Form("h_PhiMass_cuts_%d", icuts);
	// 	h_PhiMass_cuts[icuts] = new TH1F(hname_PhiMass_cuts, ";m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 200, 1.005, 1.035);
	// }

	// ######## fo(980) #############
	
	// for (int icuts = 0; icuts < ncuts; ++icuts)
	// {
	// 	TString hname_foMass_cuts = Form("h_foMass_cuts_%d", icuts);
	// 	h_foMass_cuts[icuts] = new TH1F(hname_foMass_cuts, ";m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 0.3, 1.2);
	// }

	// ######## Y(2175) #############
	// for (int icuts = 0; icuts < ncuts; ++icuts)
	// {
	// 	TString hname_YMass_cuts = Form("h_YMass_cuts_%d", icuts);
	// 	h_YMass_cuts[icuts] = new TH1F(hname_YMass_cuts, ";m_{#phi#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 200, 2.05, 2.35); //2.5, 3.2, 2.05, 2.35
	// }

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t_kin");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t_meas");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("t_p2p4");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mm2_piask");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("me");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pt");

	//CREATE HELPER AND INITIALIZE WITH DESIRED COMBINATIONS TO BE STORED
	dComboTreeHelper = new DComboTreeHelper(dTreeInterface, dComboWrapper, dFlatTreeInterface, "K+ K-; K+ pi-; K- pi+; pi- p; K- p; pi+ pi-; K+ K- p; K+ K- pi+; K+ K- pi-; K+ pi+ pi-; K- pi+ pi-; K+ pi- p; K- pi+ p; K+ K- pi+ pi-", dThrownWrapper ,"acc:pid");

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
}

Bool_t DSelector_pippimkpkm::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	if (locEntry%1000==0) cout <<"ENTRY "<<locEntry<<",  RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE
	// masses: Phi(1020), fo(980), Y(2175).
	set<map<unsigned int, set<Int_t>>> locUsedSoFar_PhiMass;
	set<map<unsigned int, set<Int_t>>> locUsedSoFar_foMass;
	set<map<unsigned int, set<Int_t>>> locUsedSoFar_YMass;

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID = dKMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlusP4_Measured + locPiMinusP4_Measured + locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured;

        /******************************** Angular distribution  ***************************************/
		// Boost in the Pip_pim_cm
		// TLorentzVector loc2piP4 = locPiPlusP4 + locPiMinusP4;
		// TLorentzVector loc2piP4_cm(loc2piP4);
		// TLorentzVector locPiPlusP4_cm(locPiPlusP4);
		// TLorentzVector locPiMinusP4_cm(locPiMinusP4);
		// // boost cms to its rest frame
		// loc2piP4_cm.Boost( -loc2piP4.BoostVector() );
		// locPiPlusP4_cm.Boost( -loc2piP4.BoostVector() );
		// locPiMinusP4_cm.Boost( -loc2piP4.BoostVector() );

		// cout <<"sum("<<loc2piP4_cm.Px()<<","<<loc2piP4_cm.Py()<<","<<loc2piP4_cm.Pz()<<","<<loc2piP4_cm.E()<<")"<<"theta = "<<loc2piP4_cm.Theta()<<endl;
		// cout <<"pip("<< locPiPlusP4_cm.Px() <<","<<locPiPlusP4_cm.Py()<<","<<locPiPlusP4_cm.Pz()<<","<<locPiPlusP4_cm.E()<<")"<<"theta = "<<locPiPlusP4_cm.Theta()<<endl;
		// cout <<"pim("<< locPiMinusP4_cm.Px() <<","<<locPiMinusP4_cm.Py()<<","<<locPiMinusP4_cm.Pz()<<","<<locPiMinusP4_cm.E()<<")"<<"theta = "<<locPiMinusP4_cm.Theta()<<endl;

		// TLorentzVector loc2piP4 = p1 + p2;
		// 1st: transform from Lab (locProtonP4.P()=0) to CM (-locBeamP4.P().BoostVector()) frame.
		TLorentzVector loc2piP4 = locPiPlusP4 + locPiMinusP4;
		TLorentzRotation resRestBoost(-loc2piP4.BoostVector());
		TLorentzVector beam_res = resRestBoost * locBeamP4;
		TLorentzVector recoil_res = resRestBoost * (locKPlusP4 + locKMinusP4); // Phi(1020) is the recoil
		TLorentzVector locPiPlusP4_res = resRestBoost * locPiPlusP4;

		// helicity frame: z-axis is propagation of resonance X => opposite recoil Phi(1020) in X rest frame
		TVector3 z = -1. * recoil_res.Vect().Unit();

		// y axis perpendicular to production plane
		TVector3 y = beam_res.Vect().Cross(z).Unit();

		// x axis is given to produce a right-handed system
		TVector3 x = y.Cross(z).Unit();

		TVector3 angles((locPiPlusP4_res.Vect()).Dot(x),
						(locPiPlusP4_res.Vect()).Dot(y),
						(locPiPlusP4_res.Vect()).Dot(z));

		Double_t cosTheta = angles.CosTheta();
		Double_t phi = angles.Phi();

        // ********************************         Mandelstam variable         **************************************
		// compute invariant t for Y(2175) with proton as a recoil
		Double_t t_kin = -2 * locProtonP4.M() * (locProtonP4.E() - locProtonP4.M());
		Double_t t_meas = -2 * locProtonP4_Measured.M() * (locProtonP4_Measured.E() - locProtonP4_Measured.M());
		Double_t t_p2p4 = (locBeamP4 - locPiPlusP4 - locPiMinusP4 - locKPlusP4 - locKMinusP4).M2();

        // ********************************         total Pt        **************************************
		Double_t pt = (locPiPlusP4_Measured + locPiMinusP4_Measured + locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured).Pt();

        // ********************************      Vertex cut      **************************************
		//Step 0
		TLorentzVector locPiPlusX4 = dPiPlusWrapper->Get_X4();
		TLorentzVector locPiMinusX4 = dPiMinusWrapper->Get_X4();
		TLorentzVector locKPlusX4 = dKPlusWrapper->Get_X4();
		TLorentzVector locKMinusX4 = dKMinusWrapper->Get_X4();
		TLorentzVector locProtonX4 = dProtonWrapper->Get_X4();
		if((locProtonX4.Z()<52 || locProtonX4.Z()>78) || (locPiPlusX4.Z()<52 || locPiPlusX4.Z()>78) || (locPiMinusX4.Z()<52 || locPiMinusX4.Z()>78)  || (locKPlusX4.Z()<52 || locKPlusX4.Z()>78)  || (locKMinusX4.Z()<52 || locKMinusX4.Z()>78))	
		{
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
        // ******************************** Missing Mass with diferent hypothesis  **************************************
		TLorentzVector locPipAsKp = locKPlusP4_Measured;
		TLorentzVector locPimAsKm = locKMinusP4_Measured;
		locPipAsKp.SetVectM(locKPlusP4_Measured.Vect(),0.139570);
		locPimAsKm.SetVectM(locKMinusP4_Measured.Vect(),0.139570);
		TLorentzVector locMissingP4_PiAsK = locBeamP4_Measured + dTargetP4;
		locMissingP4_PiAsK -= locPipAsKp + locPimAsKm + locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured;

		// ********************************************        Chi2            *******************************************
		double locKinFitChiSq = dComboWrapper->Get_ChiSq_KinFit("");

		// if(dComboWrapper->Get_ConfidenceLevel_KinFit() <= 5.73303e-7) //5.73303e-7
		//   {
		//     dComboWrapper->Set_IsComboCut(true);
		//     continue;
		//   }

		if (locKinFitChiSq > 60)
		{
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		// ******************************************* ACCIDENTAL SUBRACTION cut *******************************************

		// measured tagger time for combo
		TLorentzVector locBeam_X4_KinFit = dComboBeamWrapper->Get_X4();

		// measured RF time for combo
		double locRFTime = dComboWrapper->Get_RFTime();

		double dTargetCenterZ = dComboWrapper->Get_TargetCenter().Z();

		// time difference between tagger and RF (corrected for production vertex position relative to target center)
		double locBeamDeltaT = locBeam_X4_KinFit.T() - (locRFTime + (locBeam_X4_KinFit.Z() - dTargetCenter.Z()) / 29.9792458);

		h_TaggerAccidentals->Fill(locBeamDeltaT);

		// calculate accidental subtraction weight based on time difference
		double AccWeight = 0.; // weight to accidentally subtracted histgorams

		if (fabs(locBeamDeltaT) < 0.5 * 4.008)
		{ // prompt signal recieves a weight of 1
			AccWeight = 1.;
		}
		else
		{ // accidentals recieve a weight of 1/# RF bunches included in TTree (8 in this case)
			AccWeight = -1. / 8.;
		}

		// if(fabs(locBeamDeltaT) > 0.5*4.008)
		// {
		// 	dComboWrapper->Set_IsComboCut(true);
		// 	continue;
		// }

		// // ******************************************** PID Delta T cut ****************************************************
		
		// double locRFTime_meas = dComboWrapper->Get_RFTime_Measured();
		// // Proton
        // double locDeltaT_p_tof = -999;
		// double locDeltaT_p_bcal = -999;
		// double locDeltaT_p_fcal = -999;
		// locP4_p = dProtonWrapper->Get_P4_Measured();
		// locX4_p = dProtonWrapper->Get_X4_Measured();
		// double locPropagatedRFTime_p = locRFTime_meas + (locX4_p.Z() - dTargetCenterZ) / 29.9792458;
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_TOF) locDeltaT_p_tof = locX4_p.T() - locPropagatedRFTime_p;
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_BCAL) locDeltaT_p_bcal = locX4_p.T() - locPropagatedRFTime_p;
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_FCAL) locDeltaT_p_fcal = locX4_p.T() - locPropagatedRFTime_p;

		// // pip
		// double locDeltaT_pip_tof = -999;
		// double locDeltaT_pip_bcal = -999;
		// double locDeltaT_pip_fcal = -999;
		// locP4_pip = dPiPlusWrapper->Get_P4_Measured();
		// locX4_pip = dPiPlusWrapper->Get_X4_Measured();
		// double locPropagatedRFTime_pip = locRFTime_meas + (locX4_pip.Z() - dTargetCenterZ) / 29.9792458;
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_TOF) locDeltaT_pip_tof = locX4_pip.T() - locPropagatedRFTime_pip;
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_BCAL) locDeltaT_pip_bcal = locX4_pip.T() - locPropagatedRFTime_pip;
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_FCAL) locDeltaT_pip_fcal = locX4_pip.T() - locPropagatedRFTime_pip;

		// // kp
		// double locDeltaT_kp_tof = -999;
		// double locDeltaT_kp_bcal = -999;
		// double locDeltaT_kp_fcal = -999;
		// locP4_kp = dKPlusWrapper->Get_P4_Measured();
		// locX4_kp = dKPlusWrapper->Get_X4_Measured();
		// double locPropagatedRFTime_kp = locRFTime_meas + (locX4_kp.Z() - dTargetCenterZ) / 29.9792458;
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_TOF) locDeltaT_kp_tof = locX4_kp.T() - locPropagatedRFTime_kp;
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_BCAL) locDeltaT_kp_bcal = locX4_kp.T() - locPropagatedRFTime_kp;
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_FCAL) locDeltaT_kp_fcal = locX4_kp.T() - locPropagatedRFTime_kp;		

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/

		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MissingMass[KPlus].insert(locKPlusTrackID);
		locUsedThisCombo_MissingMass[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			h_mm2_precut->Fill(locMissingMassSquared, AccWeight);
			h2_mm2Vschi2->Fill(locKinFitChiSq, locMissingMassSquared, AccWeight);			
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		if ((locMissingMassSquared < -0.03) || (locMissingMassSquared > 0.03))
		{
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		// **************************** Phi, fo and Y invariant masses *********************************************

		// Measured
		double PhiMass_Measured = (locKPlusP4_Measured + locKMinusP4_Measured).M();
		double foMass_Measured = (locPiPlusP4_Measured + locPiMinusP4_Measured).M();
		double YMass_Measured = (locKPlusP4_Measured + locKMinusP4_Measured + locPiPlusP4_Measured + locPiMinusP4_Measured).M();

		// P4 & Vertex kinematic fit
		double PhiMass_KinFit = (locKPlusP4 + locKMinusP4).M();
		double foMass_KinFit = (locPiPlusP4 + locPiMinusP4).M();
		double YMass_KinFit = (locKPlusP4 + locKMinusP4 + locPiPlusP4 + locPiMinusP4).M();

		map<unsigned int, set<Int_t>> locUsedThisCombo_PhiMass;
		locUsedThisCombo_PhiMass[abs(PDGtype(KPlus))].insert(locKPlusTrackID);
		locUsedThisCombo_PhiMass[abs(PDGtype(KMinus))].insert(locKMinusTrackID);
		locUsedThisCombo_PhiMass[Unknown].insert(locBeamID); //beam
		map<unsigned int, set<Int_t>> locUsedThisCombo_foMass;
		locUsedThisCombo_foMass[abs(PDGtype(PiPlus))].insert(locPiPlusTrackID);
		locUsedThisCombo_foMass[abs(PDGtype(PiMinus))].insert(locPiMinusTrackID);
		locUsedThisCombo_foMass[Unknown].insert(locBeamID); //beam
		map<unsigned int, set<Int_t>> locUsedThisCombo_YMass;
		locUsedThisCombo_YMass[abs(PDGtype(KPlus))].insert(locKPlusTrackID);
		locUsedThisCombo_YMass[abs(PDGtype(KMinus))].insert(locKMinusTrackID);
		locUsedThisCombo_YMass[abs(PDGtype(PiPlus))].insert(locPiPlusTrackID);
		locUsedThisCombo_YMass[abs(PDGtype(PiMinus))].insert(locPiMinusTrackID);
		locUsedThisCombo_YMass[Unknown].insert(locBeamID); //beam

		// ********* cuts
		h_TaggerAccidentals_postcut->Fill(locBeamDeltaT, AccWeight);
		h_chi2_precut->Fill(locKinFitChiSq, AccWeight);
	    h_chi2_postcut->Fill(locKinFitChiSq, AccWeight);
		h_KinFitCL->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(), AccWeight);
		h_mm2_postcut->Fill(locMissingMassSquared, AccWeight);
        // // dtcut
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_TOF) h2_DeltaT_p_tof->Fill(locP4_p.P(), locDeltaT_p_tof);//, AccWeight);
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_BCAL) h2_DeltaT_p_bcal->Fill(locP4_p.P(), locDeltaT_p_bcal);//, AccWeight);
		// if(dProtonWrapper->Get_Detector_System_Timing() == SYS_FCAL) h2_DeltaT_p_fcal->Fill(locP4_p.P(), locDeltaT_p_fcal);//, AccWeight);
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_TOF) h2_DeltaT_pip_tof->Fill(locP4_pip.P(), locDeltaT_pip_tof);//, AccWeight);
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_BCAL) h2_DeltaT_pip_bcal->Fill(locP4_pip.P(), locDeltaT_pip_bcal);//, AccWeight);
		// if(dPiPlusWrapper->Get_Detector_System_Timing() == SYS_FCAL) h2_DeltaT_pip_fcal->Fill(locP4_pip.P(), locDeltaT_pip_fcal);//, AccWeight);
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_TOF) h2_DeltaT_kp_tof->Fill(locP4_kp.P(), locDeltaT_kp_tof);//, AccWeight);
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_BCAL) h2_DeltaT_kp_bcal->Fill(locP4_kp.P(), locDeltaT_kp_bcal);//, AccWeight);
		// if(dKPlusWrapper->Get_Detector_System_Timing() == SYS_FCAL) h2_DeltaT_kp_fcal->Fill(locP4_kp.P(), locDeltaT_kp_fcal);//, AccWeight);

		// ######## Phi(1020) #############
		if (locUsedSoFar_PhiMass.find(locUsedThisCombo_PhiMass) == locUsedSoFar_PhiMass.end())
		{

			if (fabs(locBeamDeltaT) < 0.5 * 4.008)
			{
				h_PhiMass_Measured->Fill(PhiMass_Measured, AccWeight);
				h_PhiMass_KinFit->Fill(PhiMass_KinFit, AccWeight);
			}

			h_PhiMass_beambunchcut->Fill(PhiMass_KinFit, AccWeight);

			if (PhiMass_KinFit > 1.005 && PhiMass_KinFit < 1.035)
			{
				h_PhiMass_postcuts->Fill(PhiMass_KinFit, AccWeight);
				// for (int icuts = 0; icuts < ncuts; ++icuts)
				// {
				// 	cutscut = cutsmax - icuts * cutsstep;
				// 	if (locMissingMassSquared < cutscut && locMissingMassSquared > -1*cutscut)
				// 	{
				// 		h_PhiMass_cuts[icuts]->Fill(PhiMass_KinFit, AccWeight);
				// 	}
				// }
			}
			locUsedSoFar_PhiMass.insert(locUsedThisCombo_PhiMass);
		}

		// ######## fo(980) #############
		if (locUsedSoFar_foMass.find(locUsedThisCombo_foMass) == locUsedSoFar_foMass.end())
		{

			if (fabs(locBeamDeltaT) < 0.5 * 4.008)
			{
				h_foMass_Measured->Fill(foMass_Measured, AccWeight);
				h_foMass_KinFit->Fill(foMass_KinFit, AccWeight);
			}

			h_foMass_beambunchcut->Fill(foMass_KinFit, AccWeight);

			h2_PhiVsfoMass_KinFit->Fill(PhiMass_KinFit, foMass_KinFit, AccWeight);

			if (PhiMass_KinFit > 1.005 && PhiMass_KinFit < 1.035)
			{
				h_foMass_postcuts->Fill(foMass_KinFit, AccWeight);
				// for (int icuts = 0; icuts < ncuts; ++icuts)
				// {
				// 	cutscut = cutsmax - icuts * cutsstep;
				// 	if (locMissingMassSquared < cutscut && locMissingMassSquared > -1*cutscut)
				// 	{
				// h_foMass_cuts[icuts]->Fill(foMass_KinFit, AccWeight);
				// 	}
				// }				
			}
			locUsedSoFar_foMass.insert(locUsedThisCombo_foMass);
		}

		// ######## Y(2175) #############

		if (locUsedSoFar_YMass.find(locUsedThisCombo_YMass) == locUsedSoFar_YMass.end())
		{

			if (fabs(locBeamDeltaT) < 0.5 * 4.008)
			{
				h_YMass_Measured->Fill(YMass_Measured, AccWeight);
				h_YMass_KinFit->Fill(YMass_KinFit, AccWeight);
			}

			h_YMass_beambunchcut->Fill(YMass_KinFit, AccWeight);

			h2_PhiVsYMass_KinFit->Fill(PhiMass_KinFit, YMass_KinFit, AccWeight);
			h2_foVsYMass_KinFit->Fill(foMass_KinFit, YMass_KinFit, AccWeight);

			if (PhiMass_KinFit > 1.005 && PhiMass_KinFit < 1.035)
			{
				h_YMass_postcuts->Fill(YMass_KinFit, AccWeight);
				// for (int icuts = 0; icuts < ncuts; ++icuts)
				// {
				// 	cutscut = cutsmax - icuts * cutsstep;
				// 	if (locMissingMassSquared < cutscut && locMissingMassSquared > -1*cutscut)
				// 	{
				// h_YMass_cuts[icuts]->Fill(YMass_KinFit, AccWeight);
				// 	}
				// }
			}
			locUsedSoFar_YMass.insert(locUsedThisCombo_YMass);
		}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/
	    dFlatTreeInterface->Fill_Fundamental<Double_t>("t_kin", t_kin);
	    dFlatTreeInterface->Fill_Fundamental<Double_t>("t_meas", t_meas);
	    dFlatTreeInterface->Fill_Fundamental<Double_t>("t_p2p4", t_p2p4);
	    dFlatTreeInterface->Fill_Fundamental<Double_t>("mm2_piask", locMissingP4_PiAsK.M2());
		dFlatTreeInterface->Fill_Fundamental<Double_t>("me", locMissingP4_Measured.E());
		dFlatTreeInterface->Fill_Fundamental<Double_t>("pt", pt);

		//FILL FLAT TREE
		dComboTreeHelper->Fill(locEntry); //fills branches for sub-combinations
		Fill_FlatTree(); //for the active combo
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
	// if(dThrownBeam != NULL)
	// 	double locEnergy = dThrownBeam->Get_P4().E();

	// //Loop over throwns
	// TLorentzVector locPiPlusP4_Thrown;
	// TLorentzVector locPiMinusP4_Thrown;
	// TLorentzVector locKPlusP4_Thrown;
	// TLorentzVector locKMinusP4_Thrown;
	// TLorentzVector locProtonP4_Thrown;

	// //Loop over throwns
	// for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	// {
	// 	//Set branch array indices corresponding to this particle
	// 	dThrownWrapper->Set_ArrayIndex(loc_i);

	// 	//Do stuff with the wrapper here ...

	// 	Particle_t thrown_pid = dThrownWrapper->Get_PID();
	// 	// cout << " loc_i=" << loc_i << " thrown_pid=" << thrown_pid << endl;
	//         TLorentzVector locP4_Thrown = dThrownWrapper->Get_P4();

	// 	if (loc_i == 2) locPiPlusP4_Thrown = locP4_Thrown; // assume order of particles as PID is zero at the moment
	// 	if (loc_i == 3) locPiMinusP4_Thrown = locP4_Thrown;
	// 	if (loc_i == 0) locKPlusP4_Thrown = locP4_Thrown;
	// 	if (loc_i == 1) locKMinusP4_Thrown = locP4_Thrown;
	// 	if (loc_i == 4) locProtonP4_Thrown = locP4_Thrown;

	// 	// **************************** Phi, fo and Y invariant masses *********************************************

	// 	double PhiMass_Thrown = (locKPlusP4_Thrown + locKMinusP4_Thrown).M();
	// 	double foMass_Thrown = (locPiPlusP4_Thrown + locPiMinusP4_Thrown).M();
	// 	double YMass_Thrown = (locKPlusP4_Thrown + locKMinusP4_Thrown + locPiPlusP4_Thrown + locPiMinusP4_Thrown).M();

	// 	h_PhiMass_Thrown->Fill(PhiMass_Thrown);
	// 	h_foMass_Thrown->Fill(foMass_Thrown);
	// 	h_YMass_Thrown->Fill(YMass_Thrown);

	// }
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_pippimkpkm::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
