void payman() {
   const Int_t nbins = 20;
   TH1F *h = new TH1F("h","test",nbins,-3,3);
   h->FillRandom("gaus",1000);
   TCanvas *c1 = new TCanvas("c1","c1",600,800);
   c1->Divide(1,2);
   c1->cd(1);
   h->SetFillColor(50);
   h->DrawCopy("bar2");
   c1->cd(2);
   Double_t stats[5]={0,0,0,0,0};
   h->PutStats(stats); // reset mean value, etc
   for (Int_t i=1;i<=nbins-5;i++)
   h->SetBinContent(i,h->GetBinContent(i+5));
   for (i=nbins-5;i<=nbins;i++) h->SetBinContent(i,0);
   h->Draw("bar2");
}
