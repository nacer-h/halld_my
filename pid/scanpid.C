void scanpid()
{
  TFile *f = new TFile("treepid_011529.root");
  TTree* pid = (TTree*)f->Get("pid");
  TH1F *hdedx = new TH1F("hdedx","hdedx",100,0,2);

  // for(int i =0;i<10000000;++i)
  // {
  // int nentries = pid->GetEntries(Form("dedx*1e5>>hdedx,,1,%d",i));
  pid.Draw("dedx*1e5>>hdedx","","",1,11112465);
  // }
  // TF1 *f2 = new TF1("f2","TMath::DiLog(x)",0,1);
  // f2->SetParameters(10,0.2,0.1);
  // f2->Draw("same");

  // TCanvas *cdedx = new TCanvas("cdedx","cdedx",700,500);
  // cdedx->cd();
  //if(hdedx.GetEntries()>30)

  // hdedx->Draw();

}
