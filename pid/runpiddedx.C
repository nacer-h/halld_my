void runpiddedx()
{
  TFile *f = new TFile("treepid_011529.root");
  TTree* pid = (TTree*)f->Get("pid");

  for(int trunc=0; trunc<100; trunc+=5)
  {
    pid->Process("piddedx.C+",Form("%i",trunc));
  }
}
