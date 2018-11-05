void runpiddedx()
{

  TFile *f = new TFile("hd_root_011529.root");
  TTree* ntp = (TTree*)f->Get("ntp");


  for(int trunc=5; trunc<105; trunc+=5)
  {
    ntp->Process("piddedx.C+",Form("%i",trunc));
  }

}
