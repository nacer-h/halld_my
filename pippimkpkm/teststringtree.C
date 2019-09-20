void teststringtree()
{
  TFile *f = new TFile("teststring.root","recreate");
  TTree *t = new TTree("ntp", "ntp");

  double height;
  TObjString gender;

  t->Branch("height", &height, "height/D");
  t->Branch("gender", &gender);

  for (int i=0;i<100000;++i)
    {
      if (gRandom->Rndm()>0.5)
	{
	  gender = "male";
	  height = gRandom->Gaus(1.75,0.1);
	}
      else
	{
	  gender = "female";
	  height = gRandom->Gaus(1.65,0.1);
	}

      t->Fill();
    }

  t->Write();
  f->Close();
} 
