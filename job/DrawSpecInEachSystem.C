{
    TFile* f1=new TFile("AnaFNEH1_rpc.root","read");
    TFile* f2=new TFile("AnaFNEH1_ows.root","read");
    TFile* f3=new TFile("AnaFN_iws.root","read");
    TH1F* h1=f1->Get("h1");
    h1->SetLineColor(kBlue);
    TH1F* h2=f2->Get("h1");
    h2->SetLineColor(kRed);
    TH1F* h3=f3->Get("AD1_PE_RRRRR");
    h3->SetLineColor(kGreen);
    double h1Num=h1->Integral(h1->FindBin(0.7),h1->FindBin(99.));
    double h2Num=h2->Integral(h2->FindBin(0.7),h2->FindBin(99.));
    double h3Num=h3->Integral(h3->FindBin(0.7),h3->FindBin(99.));
    h1->Scale(h2Num/h1Num);
    //h2->Scale(h2Num/h2Num);
    h3->Scale(h2Num/h3Num);
    //h1->Rebin(4);
    //h2->Rebin(4);
    //h3->Rebin(4);
    TCanvas *c=new TCanvas("specInEachSystem","specInEachSystem",800,600);
    h3->Draw();
    h1->Draw("same");
    h2->Draw("same");
    c->SaveAs("specInEachSystem.eps");

}
