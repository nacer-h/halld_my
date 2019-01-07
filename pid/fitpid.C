double getmisid(TF1* fbkg, TF1 *fsig, double &xcut, double eff=0.95)
{
	double xmin, xmax;
	fsig->GetRange(xmin,xmax);
	double integ = 0;
	double step = (xmax - xmin)/4;
	xcut = (xmax + xmin)/2;
	double normsig = fsig->Integral(xmin, xmax);
	double diffinteg = 1;
	// int maxloop = 100;
	while(diffinteg > 1e-6) // && --maxloop>0)
	{
		integ = fsig->Integral(xcut, xmax)/normsig;
		diffinteg = fabs(integ - eff);
		if(integ > eff) xcut += step;
		else xcut -= step;
		step /= 2;
	}

	double misid = fbkg->Integral(xcut,xmax)/fbkg->Integral(xmin,xmax);

	printf("integ = %f | misid = %f | xcut = %f xmin = %f | xmax = %f\n", integ, misid,xcut,xmin,xmax);

	return misid;
}

void fitpid(int trunc = 0)
{
	TFile *fpid = new TFile("fig/pid.root");
	TFile *fpid1 = new TFile("fig/pid1.root","UPDATE");  // file to save in all the plots.

	gStyle->SetOptStat(0);

	TF1 dg("dg","[5]*([9]/([9]+1.0)/[7]*exp(-0.5*((x-[6])/[7])^2)+1.0/([9]+1.0)/[8]*exp(-0.5*((x-[6])/[8])^2))/sqrt(2.0*pi)");

	TString gexp1 = "[0]*((((x-[1])/[2])<=-[3])*(exp(0.5*[3]^2+[3]*((x-[1])/[2]))) + (((x-[1])/[2])>-[3] && ((x-[1])/[2])<=[4])*(exp(-0.5*(x-[1])^2/[2]^2)) + (((x-[1])/[2])>[4])*(exp(0.5*[4]^2-[4]*((x-[1])/[2]))))";
	TString gexp2 = "[5]*((((x-[6])/[7])<=-[8])*(exp(0.5*[8]^2+[8]*((x-[6])/[7]))) + (((x-[6])/[7])>-[8] && ((x-[6])/[7])<=[9])*(exp(-0.5*(x-[6])^2/[7]^2)) + (((x-[6])/[7])>[9])*(exp(0.5*[9]^2-[9]*((x-[6])/[7]))))";

	TF1 fbkg("fbkg",gexp1);
	fbkg.SetLineColor(4); fbkg.SetLineStyle(7);
	TF1 fsig("fsig",gexp1);
	fsig.SetLineColor(8); fsig.SetLineStyle(9);
	TF1 fcomb("fcomb",gexp1+" + "+gexp2);

	TString hnamepthetamisid="hpthetamisid_";
	hnamepthetamisid+=(int)(trunc);
	TH2F *hpthetamisid = new TH2F(hnamepthetamisid, "Momentum vs. #theta vs. mis-id;#theta [deg]; P [GeV/c]; mis-id", 20,20,60,20,0.4,0.8);

	TCanvas *cddedx=new TCanvas("cddedx","cddedx",10,10,1250,950);
	TString cname[20][20];
	int histcnt = 0;
	TPolyMarker *pm;

	for (int p=0;p<20;++p)
	{
		// cddedx->Divide(4,5,0.0001,0.0001);

		for (int tht=0;tht<20;++tht)
		{

			cddedx->cd((histcnt++)%20+1);
			TString hname = Form("hdedxav_%d_%d_%d",p,tht,trunc);

			TH1F *h = (TH1F*) fpid->Get(hname);

			cout <<hname<<" ***** h = "<<h<<endl;

			if (h->Integral(1,100)<10) continue;

			double xmin = h->GetXaxis()->GetXmin();
			double xmax = h->GetXaxis()->GetXmax();

			if (h->GetListOfFunctions() !=0) pm = (TPolyMarker*)h->GetListOfFunctions()->At(0);
			if(pm == 0) continue;
			double npeaks = pm->GetN();
			if(npeaks < 2) continue;
			double peakbkg = pm->GetX()[0];
			double peaksig = pm->GetX()[1];

			printf("trunc = %i | p = %i | theta = %i | npeaks = %d\n", trunc, p, tht, npeaks);

			// double peakbkg = h->GetXaxis()->FindBin(peakbinbkg[0]);
			// double xminbkg = h->GetBinCenter(peakbkg);

			// fcomb.SetParameters(150,-0.4,0.05,1,1,50,0.0,0.05,1,1);
			fcomb.SetParameters(150,-6,0.5,10,1,50,0.0,0.5,1,1);
			if(peakbkg > peaksig) {double tmp=peakbkg;peakbkg=peaksig;peaksig=tmp;}
			if (npeaks>1) fcomb.SetParameters(150,peakbkg,0.05,1,1,50,peaksig,0.05,1,1);
			fcomb.SetParLimits(8,1,5);
			h->Fit("fcomb","","",xmin,xmax);
			h->Draw();

			double *pars = fcomb.GetParameters();
			double *pars5 = &(pars[5]);

			fbkg.SetRange(xmin,xmax);
			fbkg.SetParameters(pars);
			h->GetListOfFunctions()->Add(new TF1(fbkg));

			fsig.SetRange(xmin,xmax);
			fsig.SetParameters(pars5);
			h->GetListOfFunctions()->Add(new TF1(fsig));

			double xcut=0.0;
			hpthetamisid->SetBinContent(tht+1, p+1, getmisid(&fbkg,&fsig, xcut));
			TLine *lxcut = new TLine(xcut, h->GetMinimum(), xcut, h->GetMaximum());
			lxcut->SetLineWidth(1);
			lxcut->SetLineStyle(9);
			lxcut->Draw();
			h->Write();
			cddedx->Update();

			cname[p][tht] = "cddedx_";
			cname[p][tht] += (int)(p);
			cname[p][tht] += "_";
			cname[p][tht]+=tht;
			cname[p][tht]+="_";
			cname[p][tht] += (int)(trunc);
			cddedx->Print(Form("fig/%s.eps",cname[p][tht].Data()),"eps");
			cddedx->Print(Form("fig/%s.root",cname[p][tht].Data()),"root");
		}
		// if(p == 15 & tht == 6)
	}

	TString cnamepthetamisid = "cpthetamisid_";
	cnamepthetamisid +=(int)(trunc);
	TCanvas *cpthetamisid= new TCanvas(cnamepthetamisid,cnamepthetamisid,700,500);
	cpthetamisid->cd();
	gStyle->SetPalette(1);
	// hpthetamisid->GetZaxis()->SetRangeUser(0.00001,1);
	// if(trunc != 95) gPad->SetLogz();
	hpthetamisid->Draw("colz");
	gPad->Update();
	if (hpthetamisid->GetListOfFunctions() !=0) TPaletteAxis *palette = (TPaletteAxis*)hpthetamisid->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.86);
  palette->SetX2NDC(0.90);
	cpthetamisid->Modified();
	hpthetamisid->Write();
	cpthetamisid->Update();
	cpthetamisid->Print(Form("fig/%s.eps",cnamepthetamisid.Data()),"eps");
	cpthetamisid->Print(Form("fig/%s.root",cnamepthetamisid.Data()),"root");

	fpid1->Print();
	fpid1->Close();

}
