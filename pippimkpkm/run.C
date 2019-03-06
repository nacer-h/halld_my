void run()//TString runnum) 
{
  TStopwatch timer;
  timer.Start();
  gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");

  // gROOT->ProcessLine(Form("DPROOFLiteManager::Process_Tree(\"/cache/halld/RunPeriod-2016-02/analysis/ver05/tree_pimpipkmkp/%s/*.root\", \"pimpipkmkp_Tree\", \"DSelector_pimpipkmkp.C+\", 8, \"pimpipkmkp_%s.root\")",runnum.Data(),runnum.Data()));

  // gROOT->ProcessLine(Form("DPROOFLiteManager::Process_Tree(\"/cache/halld/RunPeriod-2017-01/analysis/ver21/tree_pippimkpkm__B4/merged/tree_pippimkpkm__B4_%s.root\", \"pippimkpkm__B4_Tree\", \"DSelector_pippimkpkm.C+\", 4, \"output/17v21/pippimkpkm__B4_17v21_mm2vschi2cut_%s.root\")",runnum.Data(),runnum.Data()));

  // DPROOFLiteManager::Process_Tree("/cache/halld/RunPeriod-2017-01/analysis/ver06/tree_pippimkpkm__U1/merged/tree_pippimkpkm__U1_0303*.root", "pippimkpkm__U1_Tree", "DSelector_pippimkpkm.C+", 8, "pippimkpkm__U1_17v6l_test.root")

  // // ********* Data *******
  // gROOT->ProcessLine("DPROOFLiteManager::Process_Tree(\"/cache/halld/RunPeriod-2017-01/analysis/ver21/tree_pippimkpkm__B4/merged/tree_pippimkpkm__B4_031*.root\", \"pippimkpkm__B4_Tree\", \"DSelector_pippimkpkm.C+\", 4, \"output/pippimkpkm__B4_17v21_chi2cut.root\")");

  // ********* Sim *******
  gROOT->ProcessLine("DPROOFLiteManager::Process_Tree(\"simulation/tree_pippimkpkm__B4_genr8_17v03.root\", \"pippimkpkm__B4_Tree\", \"DSelector_pippimkpkm.C+\", 8, \"dataout/mc/pippimkpkm_genr8_17v03_test.root\")");

  timer.Stop();
  cout << "RealTime: " << timer.RealTime() << " seconds" << endl;
  cout << "CpuTime: " << timer.CpuTime() << " seconds" << endl;
  // gROOT->ProcessLine(".q");
}
