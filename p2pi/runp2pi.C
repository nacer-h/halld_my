void runp2pi(TString particle)
{
  // TString particle[4] = {"protonflat", "pip","proton","pim"};
  int ntrunca = 0;
  // for(int j = 0; j<4;++j)
  // {
  // if(j == 0) TFile *f = new TFile("/data.local/nacer/halld_my/p2pi/tree_p2pi_proton_17v1h.root");
  // else TFile *f = new TFile(Form("/data.local/nacer/halld_my/p2pi/tree_p2pi_%s_17v1h.root",particle[j].Data()));
  if(particle == "protonflat") TFile *f = new TFile("/data.local/nacer/halld_my/p2pi/tree_p2pi_proton_17v1h.root");
  else TFile *f = new TFile(Form("/data.local/nacer/halld_my/p2pi/tree_p2pi_%s_17v1h.root",particle.Data()));
  TTree* pid = (TTree*)f->Get("pid");

  for(int trunclevel=0; trunclevel<100; trunclevel+=5)
  {
    for(int trunclevel2=0; trunclevel2<100; trunclevel2+=5)
    {
      if(trunclevel + trunclevel2 >99) continue;
      pid->Process("p2pi.C+",Form("%i:%i:%s",trunclevel,trunclevel2,particle.Data()));
      ntrunca++;
    }
  }
  printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++ ntrunca=%i\n", ntrunca);
  // }
}
