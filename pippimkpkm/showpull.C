TH1F *showpull(TH1F *h, TF1 *ff1)
{
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.251);

    pad1->SetBottomMargin(0.0);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.4);
    pad2->SetBorderMode(0);

    pad1->Draw();
    pad2->Draw();

    TH1F *hres = new TH1F("hres", "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    hres->SetStats(0);
    hres->GetXaxis()->SetLabelSize(0.18);
    hres->GetXaxis()->SetTitleSize(0.18);
    hres->GetXaxis()->SetTitleOffset(1.1);
    hres->GetXaxis()->SetTickLength(0.05);

    hres->GetYaxis()->SetLabelSize(0.15);
    hres->GetYaxis()->SetTitleSize(0.15);
    hres->GetYaxis()->SetTitleOffset(0.35);
    hres->GetYaxis()->SetNdivisions(5);

    hres->SetXTitle(h->GetXaxis()->GetTitle());
    hres->SetYTitle("pull [#sigma] ");

    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        if (h->GetBinError(i) != 0)
        {
            hres->SetBinContent(i, (h->GetBinContent(i) - ff1->Eval(h->GetBinCenter(i))) / h->GetBinError(i));
            hres->SetBinError(i, 1.0);
        }
    }
    hres->GetYaxis()->SetRangeUser(-5, 5);
    hres->GetYaxis()->SetNdivisions(505);

    pad2->cd();
    hres->Draw("e");
    pad2->SetGridy();

    // extra
    TLine l;
    l.SetLineColor(4);
    l.DrawLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);

    pad1->cd();

    // pull up the minimum so that a minimum 0 on the y-axis is not cut

    double ymax = 0;
    // this is the maximum for (bin content + bin error)
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        if ((h->GetBinContent(i) + h->GetBinError(i)) > ymax)
            ymax = (h->GetBinContent(i) + h->GetBinError(i));
    double ymin = -h->GetYaxis()->GetLabelSize() * 0.5 * ymax;

    if (gPad->GetFrame()->GetY1() > ymin)
    {
        //h->SetMinimum(-0.05*(h->GetYaxis()->GetXmax()-h->GetXaxis()->GetXmin()));
        h->GetYaxis()->SetRangeUser(ymin, ymax * 1.05);
    }

    return hres;
}