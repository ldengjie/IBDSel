{
    char *labels[4] = {"P12E","P13A","P12E+P13A","P14A"};
    double x[4]={1,2,3,4};
    double y1[3][4]={{1.6/100,3.5/100,11.3/100,18.03/100},
                     {9.5/100,31.1/100,11.3/100,38.9/100},
                     {18.3/100,3.3/100,15.6/100,9.3/100}};
    double y2[3][4]={{-37.8/100,-34.8/100,-28.5/100,-27.1/100},
                     {-37.1/100,-44.3/100,-28.5/100,-39.1/100},
                     {-27.96/100,-22.5/100,-35.4/100,-38.1/100}};
    double livetime[3][4]={{190.995*2,198.406+198.406,190.995*2+198.406+198.406,565.436+565.436},
                           {189.646,202.683+202.683,189.646+202.683+202.683,568.03+378.407},
                           {189.786*3,200.726*4,189.786*3+200.726*4,562.451*3+372.685}};
    double firstChiNdf[3][4]={{266.9/173,232./173,210./173,172./173},
                              {128.8/149,168./171,210./173,164.8/173},
                              {26.9/99,68./136,127./153,156./166}};
    double zeroChiNdf[3][4]={{356./174,325./174,383./174,483./174},
                              {165./150,407./172,383./174,673./174},
                              {26.7/100,72./137,134./154,181./167}};
    TCanvas* c[3];
    TGraph* g1[3],g2[3];
    //TText *t1[3][4]; 
    //TText *t2[3][4]; 
    TText *t=new TText(); 
    TString nameStr;
    for( int i=0 ; i<3 ; i++ )
    {

        nameStr=Form("EH%i",i+1);
        c[i]=new TCanvas(nameStr,nameStr,800,300);
        c[i]->SetGridy();
        g1[i]=new TGraph(4,x,y1[i]);
        nameStr=Form("EH%i fast neutron system error",i+1);
        g1[i]->SetTitle(nameStr);
        g1[i]->GetHistogram()->SetMaximum(40./100);
        g1[i]->GetHistogram()->SetMinimum(-50./100);
        g1[i]->SetLineColor(4);
        g1[i]->Draw("AC*");
        g2[i] =new TGraph(4,x,y2[i]);
        g2[i]->SetLineColor(6);
        g2[i]->Draw("C*");
        for( int j=0 ; j<4 ; j++ )
        {
            nameStr=Form("%.2f",firstChiNdf[i][j]);
            t->DrawText(x[j],y1[i][j]+0.02,nameStr);
            nameStr=Form("%.2f",zeroChiNdf[i][j]);
            t->DrawText(x[j],y2[i][j]-0.05,nameStr);
            nameStr=Form("%s(%.0fdays)",labels[j],livetime[i][j]);
            g1[i]->GetHistogram()->GetXaxis()->SetBinLabel(g1[i]->GetHistogram()->FindBin(x[j]),nameStr);
        }
        g1[i]->GetHistogram()->LabelsOption("h","X");
        nameStr=Form("sysErr/EH%i.eps",i+1);
        c[i]->SaveAs(nameStr);
        
    }

/*
    //TText *t = new TText();
    //t->SetTextAlign(32);
    //t->SetTextSize(0.035);
    //t->SetTextFont(72);
    //char *labels[4] = {"P12E","P13A","P12E+P13A","P14A"};
    //for (Int_t i=0;i<4;i++) {
    //t->DrawText(x[i],0,labels[i]);
    //}

    TCanvas* c2=new TCanvas("EH2","EH2",800,300);
    c2->SetGridy();
    TGraph *g21 =new TGraph(4,x,y21);
    g21->SetTitle("EH2 system error");
    g21->GetHistogram()->SetMaximum(40./100);
    g21->GetHistogram()->SetMinimum(-50./100);
    g21->SetLineColor(4);
    for( int i=0 ; i<4 ; i++ )
    {
        g21->GetHistogram()->GetXaxis()->SetBinLabel(g21->GetHistogram()->FindBin(x[i]),labels[i]);
    }
    g21->GetHistogram()->LabelsOption("h","X");
    g21->Draw("AC*");
    TGraph *g22 =new TGraph(4,x,y22);
    g22->SetLineColor(6);
    g22->Draw("C*");


    TCanvas* c3=new TCanvas("EH3","EH3",800,300);
    c3->SetGridy();
    TGraph *g31 =new TGraph(4,x,y31);
    g31->SetTitle("EH3 system error");
    g31->GetHistogram()->SetMaximum(40./100);
    g31->GetHistogram()->SetMinimum(-50./100);
    g31->SetLineColor(4);
    for( int i=0 ; i<4 ; i++ )
    {
        g31->GetHistogram()->GetXaxis()->SetBinLabel(g31->GetHistogram()->FindBin(x[i]),labels[i]);
    }
    g31->GetHistogram()->LabelsOption("h","X");
    g31->Draw("AC*");
    TGraph *g32 =new TGraph(4,x,y32);
    g32->SetLineColor(6);
    g32->Draw("C*");

*/

}
