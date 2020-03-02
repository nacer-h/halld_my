#include "DSelector_p2k2pi_thr.h"

void DSelector_p2k2pi_thr::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	// dOutputFileName = "p2k2pi_thr.root"; //"" for none

    TString option = GetOption();
    if (option=="") option = "pippimkpkm";

    dOutputFileName = option + ".root";//option + ".root"; //"" for none
	// dOutputTreeFileName = option +  ".root"; //"" for none
	// dFlatTreeFileName = option + ".root"; //output flat tree (one combo per tree entry), "" for none
	// dFlatTreeName = "ntp"; //if blank, default name will be chosen

	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools
	// dOutputTreeFileNameMap["Bin1"] = "mcgen_bin1.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["Bin1"] = "mcgen_bin1.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["Bin2"] = "mcgen_bin2.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["Bin3"] = "mcgen_bin3.root"; //key is user-defined, value is output file name

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	dPreviousRunNumber = 0;

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/
	// Thrown  masses: Phi(1020), fo(980), Y(2175)
	h_PhiMass_Thrown = new TH1F("PhiMass_Thrown", "Generated ;m_{K^{+}K^{-}} [GeV/c^{2}];Counts", 600, 0.98, 1.2);
	h_foMass_Thrown = new TH1F("foMass_Thrown", "Generated ;m_{#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 0.3, 1.2);
	h_YMass_Thrown = new TH1F("YMass_Thrown", "Generated ;m_{K^{+}K^{-}#pi^{+}#pi^{-}} [GeV/c^{2}];Counts", 600, 1.7, 3.2);
	h_beame_Thrown = new TH1F("h_beame_Thrown", "Generated ;E_{#gamma} [GeV];Counts", 100, 6.5, 11.6);
	h_t_Thrown = new TH1F("h_t_Thrown", "Generated ;-t [GeV^{2}];Counts", 100, 0, 4);
	h2_beamevst_Thrown = new TH2D("h2_beamevst_Thrown", "Generated; -t [GeV^{2}]; E_{#gamma} [GeV]", 100, 0, 4, 100, 6.5, 11.6);
	
	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
	// dFlatTreeInterface->Register_GetEntryBranch("t_thr");
}

Bool_t DSelector_p2k2pi_thr::Process(Long64_t locEntry)
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
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	if (locEntry%100000==0) cout <<"ENTRY "<<locEntry<<",  RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
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

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/******************************************* LOOP OVER THROWN DATA ***************************************/

	//Thrown beam: just use directly
	double locBeamEnergyUsedForBinning = 0.0;
	if(dThrownBeam != NULL)
		locBeamEnergyUsedForBinning = dThrownBeam->Get_P4().E();

	// throwns P4 for final state particles
	TLorentzVector locPiPlusP4_Thrown;
	TLorentzVector locPiMinusP4_Thrown;
	TLorentzVector locKPlusP4_Thrown;
	TLorentzVector locKMinusP4_Thrown;
	TLorentzVector locProtonP4_Thrown;
	
	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		TLorentzVector locThrownP4 = dThrownWrapper->Get_P4();
		//cout << "Thrown " << loc_i << ": " << locPID << ", " << locThrownP4.Px() << ", " << locThrownP4.Py() << ", " << locThrownP4.Pz() << ", " << locThrownP4.E() << endl;

		// loop over throwns P4 for final state particles
		if (loc_i == 0)
		  locKPlusP4_Thrown = locThrownP4; // assume order of particles as PID is zero at the moment
		if (loc_i == 1)
		  locKMinusP4_Thrown = locThrownP4;
		if (loc_i == 2)
		  locPiPlusP4_Thrown = locThrownP4;
		if (loc_i == 3)
		  locPiMinusP4_Thrown = locThrownP4;
		if (loc_i == 4)
		  locProtonP4_Thrown = locThrownP4;
	}

	// **************************** Phi, fo and Y invariant masses *********************************************
	double PhiMass_Thrown = (locKPlusP4_Thrown + locKMinusP4_Thrown).M();
	double foMass_Thrown = (locPiPlusP4_Thrown + locPiMinusP4_Thrown).M();
	double YMass_Thrown = (locKPlusP4_Thrown + locKMinusP4_Thrown + locPiPlusP4_Thrown + locPiMinusP4_Thrown).M();
	// compute invariant -t for Y(2175) with proton as a recoil
	double t_thr = -2 * locProtonP4_Thrown.M() * (locProtonP4_Thrown.E() - locProtonP4_Thrown.M());
	
	h_PhiMass_Thrown->Fill(PhiMass_Thrown);
	h_foMass_Thrown->Fill(foMass_Thrown);
	h_YMass_Thrown->Fill(YMass_Thrown);
	h_beame_Thrown->Fill(locBeamEnergyUsedForBinning);
	h_t_Thrown->Fill(-t_thr);
	h2_beamevst_Thrown->Fill(-t_thr,locBeamEnergyUsedForBinning);

	// if (locEntry%100000==0) cout <<"locBeamEnergyUsedForBinning = "<<locBeamEnergyUsedForBinning<<" | t_thr = " <<-t_thr<< endl;
	
	//OR Manually:
	//BEWARE: Do not expect the particles to be at the same array indices from one event to the next!!!!
	//Why? Because while your channel may be the same, the pions/kaons/etc. will decay differently each event.

	//BRANCHES: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat#TTree_Format:_Simulated_Data
	TClonesArray** locP4Array = dTreeInterface->Get_PointerToPointerTo_TClonesArray("Thrown__P4");
	TBranch* locPIDBranch = dTreeInterface->Get_Branch("Thrown__PID");

	// Particle_t locKPlusPID_Thrown = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	// locKPlusP4_Thrown = *((TLorentzVector*)(*locP4Array)->At(0));
	// Particle_t locKMinusPID_Thrown = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	// locKMinusP4_Thrown = *((TLorentzVector*)(*locP4Array)->At(1));
	// Particle_t locPiPlusPID_Thrown = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[2]);
	// locPiPlusP4_Thrown = *((TLorentzVector*)(*locP4Array)->At(2));
	// Particle_t locPiMinusPID_Thrown = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[3]);
	// locPiMinusP4_Thrown = *((TLorentzVector*)(*locP4Array)->At(3));
	// Particle_t locProtonPID_Thrown = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[4]);
	// locProtonP4_Thrown = *((TLorentzVector*)(*locP4Array)->At(4));

    // dFlatTreeInterface->Fill_Fundamental<Double_t>("t_thr", t_thr);
/*
	Particle_t locThrown1PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	TLorentzVector locThrown1P4 = *((TLorentzVector*)(*locP4Array)->At(0));
	cout << "Particle 1: " << locThrown1PID << ", " << locThrown1P4.Px() << ", " << locThrown1P4.Py() << ", " << locThrown1P4.Pz() << ", " << locThrown1P4.E() << endl;
	Particle_t locThrown2PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	TLorentzVector locThrown2P4 = *((TLorentzVector*)(*locP4Array)->At(1));
	cout << "Particle 2: " << locThrown2PID << ", " << locThrown2P4.Px() << ", " << locThrown2P4.Py() << ", " << locThrown2P4.Pz() << ", " << locThrown2P4.E() << endl;
*/


	/******************************************* BIN THROWN DATA INTO SEPARATE TREES FOR AMPTOOLS ***************************************/
/*
	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION (along with the output file names)
	if((locBeamEnergyUsedForBinning >= 8.0) && (locBeamEnergyUsedForBinning < 9.0))
		Fill_OutputTree("Bin1"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 9.0) && (locBeamEnergyUsedForBinning < 10.0))
		Fill_OutputTree("Bin2"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 10.0) && (locBeamEnergyUsedForBinning < 11.0))
		Fill_OutputTree("Bin3"); //your user-defined key
*/

	return kTRUE;
}

void DSelector_p2k2pi_thr::Finalize(void)
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
