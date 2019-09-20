void plotTopology()
{

  gStyle->SetOptStat("");
  gStyle->SetPaintTextFormat("4.3f"); // set precision for histogram text
  gStyle->SetLegendBorderSize(0);

  TFile *f = TFile::Open("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/sim/tree_pippimkpkm_bggen_17v3_postcut_flat.root");

  // plot percentages for different thrown topologies
  TCanvas *ctopocent = new TCanvas("ctopocent","ctopocent",900,600);
  TH1F *hThrownTopologies = (TH1F*)f->Get("hThrownTopologies");
  hThrownTopologies->GetXaxis()->LabelsOption(">"); // order by most common topology
  hThrownTopologies->GetXaxis()->SetRangeUser(0, 10); // only plot first 20 topologies
  hThrownTopologies->Scale(100./hThrownTopologies->GetEntries()); // turn histogram into percentage
  hThrownTopologies->GetYaxis()->SetTitle("Thrown Topology %");
  hThrownTopologies->Draw("htext");

  // draw invariant mass histograms for different thrown topologies
  TCanvas *ctopomass = new TCanvas("ctopomass","ctopomass",900,600);
  TLegend *leg = new TLegend(0.15,0.6,0.4,0.85);
  leg->SetNColumns(2);
  THStack *hStack = new THStack("hstack","");

  int locNumTopologies = 10;
  for(int i=0; i<locNumTopologies; i++)
    {
      TH1I *hMassThrownTopology = (TH1I *)f->Get(Form("hInvariantMass_ThrownTopology_%d", i));
      if(!hMassThrownTopology) break;
      hMassThrownTopology->SetLineColor(kBlack+i);
      // hMassThrownTopology->SetFillColor(kBlack+i);
      hMassThrownTopology->Rebin(10);

      TString locLegendTitle = hMassThrownTopology->GetTitle();
      locLegendTitle.ReplaceAll("Invariant Mass Topology:","");
      leg->AddEntry(hMassThrownTopology,locLegendTitle,"l");

      if(i==0)
	{
	  hMassThrownTopology->SetTitle("#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus} Invariant Mass; M_{#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}} (GeV)");
	  hMassThrownTopology->Draw();
	}
      else
	hMassThrownTopology->Draw("same");
      hStack->Add(hMassThrownTopology);
    }

  leg->Draw("same");

  TCanvas *ctopomasstack = new TCanvas("ctopomasstack","ctopomasstack",900,600);
  hStack->SetTitle("#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus} Invariant Mass; M_{#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}} (GeV)");
  hStack->Draw();
  leg->Draw("same");

  ctopocent->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopocent.root");
  ctopocent->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopocent.eps");
  ctopomass->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopomass.root");
  ctopomass->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopomass.eps");
  ctopomasstack->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopomasstack.root");
  ctopomasstack->Print("/lustre/hebe/panda/ahamdi/gluex_root_analysis/workdir/dataout/ctopomasstack.eps");

  return;
}
