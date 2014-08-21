#define Ibd_cxx

#include "Ibd.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include  <iostream>
#include  <fstream>
#include  <sstream>
#include  "TGraph.h"
#include  "TGraphQQ.h"
#include  "TPaveLabel.h"
//#include "iomanip.h"
#include  <TFitResult.h>
#include  <TLegend.h>
#include  <TLine.h>
#include  <TGraph.h>
#include <iomanip>
using namespace std;

void Ibd::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();
    //TString DataVerTmp=option(3,4);
    //string DataVerTmp=option.substr(3);
    dataVer=option;
    dataVer=dataVer.substr(3);
    option=option(0,3);
    siteStr=option;
    std::cout<<"dataVer  : "<<dataVer<<endl;
    int daqHistNum=0;
    if( dataVer.find("13")!=string::npos || dataVer.find("14")!=string::npos || dataVer.find("15")!=string::npos || dataVer.find("16")!=string::npos)
    {
        daqHistNum=5;
    } else
    {
        if( dataVer.find("12")!=string::npos || dataVer.find("11")!=string::npos)
        {
            daqHistNum=4;
        }
    }
    std::cout<<"option  : "<<option<<endl;
    //livetime={{85.7705,85.4113,0,0},{98.2744,0,0,0},{112.597,112.596,112.56,0}};
    //livetime0={{85.8167,85.4573,0,0},{98.3095,0,0,0},{112.61,112.609,112.572,0}};
    //double livetime={{85.7705,85.4113},{98.2744},{112.597,112.596,112.56}};
    //double livetime0={{85.8167,85.4573},{98.3095},{112.61,112.609,112.572}};

    for( int i=0 ; i<4 ; i++ )
    {
        UpperNum[i]=0.;
        LowerNum[i]=0.;
    }

    std::cout<<"start to read file "<<endl;
    for( int i=0 ; i<53 ; i++ )
    {
        for( int j=0 ; j<8 ; j++ )
        {
            t_result[i][j]=0.;    
        }
    }
    ifstream infile; 
    string infileName=dataVer+"/"+"result_temp_"+dataVer+".txt";
    infile.open(infileName.c_str(),ios::in);
    int inLineNum=0;
    std::cout<<"get value "<<endl;
    while(inLineNum<53)
    {
        for( int j=0 ; j<8 ; j++ )
        {
            infile>>t_result[inLineNum][j];
        }
        inLineNum++;
        //std::cout<<"inLineNuminLineNum  : "<<inLineNum<<endl;
    }
    infile.close();
    std::cout<<"done "<<endl;

    daqtime[0][0]=t_result[26][0];
    daqtime[0][1]=t_result[26][1];
    daqtime[1][0]=t_result[26][2];
    daqtime[1][1]=t_result[26][3];
    daqtime[2][0]=t_result[26][4];
    daqtime[2][1]=t_result[26][5];
    daqtime[2][2]=t_result[26][6];
    daqtime[2][3]=t_result[26][7];

    daqtime0[0][0]=t_result[2][0];
    daqtime0[0][1]=t_result[2][1];
    daqtime0[1][0]=t_result[2][2];
    daqtime0[1][1]=t_result[2][3];
    daqtime0[2][0]=t_result[2][4];
    daqtime0[2][1]=t_result[2][5];
    daqtime0[2][2]=t_result[2][6];
    daqtime0[2][3]=t_result[2][7];

    epsi_mu[0][0]=t_result[28][0];
    epsi_mu[0][1]=t_result[28][1];
    epsi_mu[1][0]=t_result[28][2];
    epsi_mu[1][1]=t_result[28][3];
    epsi_mu[2][0]=t_result[28][4];
    epsi_mu[2][1]=t_result[28][5];
    epsi_mu[2][2]=t_result[28][6];
    epsi_mu[2][3]=t_result[28][7];

    epsi_mu0[0][0]=t_result[3][0];
    epsi_mu0[0][1]=t_result[3][1];
    epsi_mu0[1][0]=t_result[3][2];
    epsi_mu0[1][1]=t_result[3][3];
    epsi_mu0[2][0]=t_result[3][4];
    epsi_mu0[2][1]=t_result[3][5];
    epsi_mu0[2][2]=t_result[3][6];
    epsi_mu0[2][3]=t_result[3][7];

    epsi_multi[0][0]=t_result[30][0];
    epsi_multi[0][1]=t_result[30][1];
    epsi_multi[1][0]=t_result[30][2];
    epsi_multi[1][1]=t_result[30][3];
    epsi_multi[2][0]=t_result[30][4];
    epsi_multi[2][1]=t_result[30][5];
    epsi_multi[2][2]=t_result[30][6];
    epsi_multi[2][3]=t_result[30][7];

    epsi_multi0[0][0]=t_result[4][0];
    epsi_multi0[0][1]=t_result[4][1];
    epsi_multi0[1][0]=t_result[4][2];
    epsi_multi0[1][1]=t_result[4][3];
    epsi_multi0[2][0]=t_result[4][4];
    epsi_multi0[2][1]=t_result[4][5];
    epsi_multi0[2][2]=t_result[4][6];
    epsi_multi0[2][3]=t_result[4][7];

    livetime[0][0]=daqtime[0][0]*epsi_mu[0][0]*epsi_multi[0][0];
    livetime[0][1]=daqtime[0][1]*epsi_mu[0][1]*epsi_multi[0][1];
    livetime[1][0]=daqtime[1][0]*epsi_mu[1][0]*epsi_multi[1][0];
    livetime[1][1]=daqtime[1][1]*epsi_mu[1][1]*epsi_multi[1][1];
    livetime[2][0]=daqtime[2][0]*epsi_mu[2][0]*epsi_multi[2][0];
    livetime[2][1]=daqtime[2][1]*epsi_mu[2][1]*epsi_multi[2][1];
    livetime[2][2]=daqtime[2][2]*epsi_mu[2][2]*epsi_multi[2][2];
    livetime[2][3]=daqtime[2][3]*epsi_mu[2][3]*epsi_multi[2][3];

    livetime0[0][0]=daqtime0[0][0]*epsi_mu0[0][0]*epsi_multi0[0][0];
    livetime0[0][1]=daqtime0[0][1]*epsi_mu0[0][1]*epsi_multi0[0][1];
    livetime0[1][0]=daqtime0[1][0]*epsi_mu0[1][0]*epsi_multi0[1][0];
    livetime0[1][1]=daqtime0[1][1]*epsi_mu0[1][1]*epsi_multi0[1][1];
    livetime0[2][0]=daqtime0[2][0]*epsi_mu0[2][0]*epsi_multi0[2][0];
    livetime0[2][1]=daqtime0[2][1]*epsi_mu0[2][1]*epsi_multi0[2][1];
    livetime0[2][2]=daqtime0[2][2]*epsi_mu0[2][2]*epsi_multi0[2][2];
    livetime0[2][3]=daqtime0[2][3]*epsi_mu0[2][3]*epsi_multi0[2][3];

    ADNum=daqHistNum-1;
    site=2;
    //if( option.EqualTo("EH1") || option.EqualTo("EH2"))
    if( option=="EH1")
    {
        site=0;
        ADNum=2;
    }
    if( option=="EH2" )
    {
        site=1;
        ADNum=daqHistNum-3;
    }
    for( int i=0 ; i<ADNum ; i++ )
    {
        Num1[i]=0;
        Numo[i]=0;
        Num2[i]=0;
    }
    tNumo=0;
    tNum1=0;
    tNum2=0;
    option+="Ibd_";
    option+=dataVer;
    option+=".root";
    //TString fileName=dataVer+"/"+option;
    TString fileName=dataVer;
    fileName+="/";
    fileName+=option;
    file = new TFile(fileName,"RECREATE");

    // Binning	
    double rnorm = 1.;
    // calculate muon E bins
    const int nBins = 200;
    //const int nBins = 200;
    double xlow = 0.001, xup = 200.;
    BinWt = new  double[nBins+1];//so i can declare it in .h file 
    //double BinWt[nBins+1] = {0.0}; 
    double EBins[nBins+1] = {0.0};
    //EBins[0] = xlow;  //should be this ,but  zhangfh's EBins[0]=0
    double BinWidth = (log10(xup)-log10(xlow))/nBins;

    for(int i=1; i<=nBins; i++) {
        EBins[i] = xlow*pow(10, BinWidth*i);
        BinWt[i] = rnorm/(EBins[i]-EBins[i-1]);
    }
    tFnProEWithrpc=new TH1F("tFnProEWithrpc","Energy of promt signal with using rpc veto",nBins, EBins);
    tFnProEWithoutrpc=new TH1F("tFnProEWithoutrpc","Energy of promt signal without using rpc veto",nBins, EBins);
    tFnProEWithoutrpcUniX=new TH1F("tFnProEWithoutrpcUniX","Energy of promt signal without using rpc veto",198, 0.7,99.7);
    for(int i=0;i<ADNum;i++)
    {
        histname="AD";
        histname+=i+1;
        histname+="t2lastmuonWithrpc";
        t2lastmuonWithrpc[i] = new TH1F(histname,"prompt time to last Muon with using rpc veto",5000,0.,50000.);
        histname="AD";
        histname+=i+1;
        histname+="t2lastmuonWithoutrpc";
        t2lastmuonWithoutrpc[i] = new TH1F(histname,"prompt time to last Muon with using rpc veto",5000,0.,50000.);
        histname="AD";
        histname+=i+1;
        histname+="FnProEWithrpc";
        FnProEWithrpc[i]=new TH1F(histname,"Energy of promt signal with using rpc veto",nBins, EBins);
        histname="AD";
        histname+=i+1;
        histname+="FnProEWithoutrpc";
        FnProEWithoutrpc[i]=new TH1F(histname,"Energy of promt signal without using rpc veto",nBins, EBins);
        histname="AD";
        histname+=i+1;
        histname+="pxy";
        pxy[i] = new TH2F(histname,"pxy",30,-3.,3.,30,-3.,3.);
        histname+="0";
        pxy0[i] = new TH2F(histname,"pxy0",30,-3.,3.,30,-3.,3.);
        histname="AD";
        histname+=i+1;
        histname+="prz";
        prz[i] = new TH2F(histname,"prz",30,0,3.,30,-3.,3.);
        histname+="0";
        prz0[i] = new TH2F(histname,"prz0",30,0,3.,30,-3.,3.);
        histname="AD";
        histname+=i+1;
        histname+="pE";
        pE[i]= new TH1F(histname,"Energy of promt signal",30,0.7 ,12.0);
        histname+="0";
        pE0[i]= new TH1F(histname,"Energy of promt signal 0",30,0.7 ,12.0);
        histname="AD";
        histname+=i+1;
        histname+="dxy";
        dxy[i] = new TH2F(histname,"dxy",30,-3.,3.,30,-3.,3.);
        histname+="0";
        dxy0[i] = new TH2F(histname,"dxy0",30,-3.,3.,30,-3.,3.);
        histname="AD";
        histname+=i+1;
        histname+="drz";
        drz[i] = new TH2F(histname,"drz",30,0,3.,30,-3.,3.);
        histname+="0";
        drz0[i] = new TH2F(histname,"drz0",30,0,3.,30,-3.,3.);
        histname="AD";
        histname+=i+1;
        histname+="dE";
        dE[i]= new TH1F(histname,"Energy of delayed signal",60,6.0,12.0);
        histname+="0";
        dE0[i]= new TH1F(histname,"Energy of delayed signal 0",60,6.0,12.0);
        histname="AD";
        histname+=i+1;
        histname+="intervalT";
        intervalT[i]= new TH1F(histname,"interval time between prompt and delayed signal",400,0.,200.);
        histname+="0";
        intervalT0[i]= new TH1F(histname,"interval time between prompt and delayed signal 0",400,0.,200.);
    }
    std::cout<<"this is begin ,done. "<<endl;
}


void Ibd::SlaveBegin(TTree * /*tree*/)
{

    TString option = GetOption();

}
void Ibd::FillIbd(int itype)
{
    if( itype==1)
    {
        pxy[det-1]->Fill(x[0]/1000,y[0]/1000);
        prz[det-1]->Fill(sqrt(x[0]*x[0]+y[0]*y[0])/1000,z[0]/1000);
        pE[det-1]->Fill(energy[0]);
        dxy[det-1]->Fill(x[1]/1000,y[1]/1000);
        drz[det-1]->Fill(sqrt(x[1]*x[1]+y[1]*y[1])/1000,z[1]/1000);
        dE[det-1]->Fill(energy[1]);
        intervalT[det-1]->Fill(timeInterval);

    }else
    {
        pxy0[det-1]->Fill(x[0]/1000,y[0]/1000);
        prz0[det-1]->Fill(sqrt(x[0]*x[0]+y[0]*y[0])/1000,z[0]/1000);
        pE0[det-1]->Fill(energy[0]);
        dxy0[det-1]->Fill(x[1]/1000,y[1]/1000);
        drz0[det-1]->Fill(sqrt(x[1]*x[1]+y[1]*y[1])/1000,z[1]/1000);
        dE0[det-1]->Fill(energy[1]);
        intervalT0[det-1]->Fill(timeInterval);

    }

}

Bool_t Ibd::Process(Long64_t entry)
{
    GetEntry(entry);

    //std::cout<<"Numo["<<det-1<<"]  : "<<Numo[det-1]<<endl;
    //TODO::4 is not enough,spesically for newAdMuon and newShowerMuon ,maybe need to modify this .
    double t2lastmuonWoutrpc=promptT2Muon[0];
    for( int i=0 ; i<4 ; i++ )
    {
        if( t2lastmuonWoutrpc>promptT2Muon[i] )
        {
            t2lastmuonWoutrpc=promptT2Muon[i];
        }
    }
    double t2lastmuonWrpc=promptT2Muon[4];
    if( t2lastmuonWrpc>t2lastmuonWoutrpc)
    {
        t2lastmuonWrpc=t2lastmuonWoutrpc;
    }
    int bin = FnProEWithrpc[det-1]->FindBin(energy[0]);
    double wt=BinWt[bin];

    if( isIbd==0 )
    {
        Numo[det-1]++;
        //Numo[det-1]++;
        //if( Numo[det-1]-Numo0>100 )
        //{
        //    std::cout<<"Numo0  : "<<Numo0<<endl;
        //    std::cout<<"Numo[det-1]  : "<<Numo[det-1]<<endl;
        //}
        tNumo++;
        FillIbd(0);
        t2lastmuonWithoutrpc[det-1]->Fill(t2lastmuonWoutrpc*1e6);
        FnProEWithoutrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithoutrpc->Fill(energy[0],wt);
        tFnProEWithoutrpcUniX->Fill(energy[0]);
        if( z[0]>0. )
        {
            UpperNum[det-1]++;
        }else
        {
            LowerNum[det-1]++;
        }

    } else if( isIbd==1 )
    {
        Num1[det-1]++;
        tNum1++;
        FillIbd(0);
        FillIbd(1);
        t2lastmuonWithoutrpc[det-1]->Fill(t2lastmuonWoutrpc*1e6);
        t2lastmuonWithrpc[det-1]->Fill(t2lastmuonWrpc*1e6);
        FnProEWithoutrpc[det-1]->Fill(energy[0],wt);
        FnProEWithrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithoutrpc->Fill(energy[0],wt);
        tFnProEWithoutrpcUniX->Fill(energy[0]);
        tFnProEWithrpc->Fill(energy[0],wt);
        if( z[0]>0. )
        {
            UpperNum[det-1]++;
        }else
        {
            LowerNum[det-1]++;
        }
    } else if( isIbd==2 )
    {
        Num2[det-1]++;
        tNum2++;
        FillIbd(1);
        t2lastmuonWithrpc[det-1]->Fill(t2lastmuonWrpc*1e6);
        FnProEWithrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithrpc->Fill(energy[0],wt);
    } else if( isIbd==3 )
    {
        FnProEWithoutrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithoutrpc->Fill(energy[0],wt);
        tFnProEWithoutrpcUniX->Fill(energy[0]);
    } else if( isIbd==4 )
    {
        FnProEWithoutrpc[det-1]->Fill(energy[0],wt);
        FnProEWithrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithoutrpc->Fill(energy[0],wt);
        tFnProEWithoutrpcUniX->Fill(energy[0]);
        tFnProEWithrpc->Fill(energy[0],wt);
    } else if (isIbd==5)
    {
        FnProEWithrpc[det-1]->Fill(energy[0],wt);
        tFnProEWithrpc->Fill(energy[0],wt);
    }else
    {
        //std::cout<<"isIbd is wrong ,please check data. "<<endl;	
    }
    //std::cout<<"this is process ,done. "<<endl;

    return kTRUE;
}

void Ibd::SlaveTerminate()
{

}

void Ibd::Terminate0()
{

}
void Ibd::Terminate()
{
    for( int i=0 ; i<ADNum ; i++ )
    {
        std::cout<<"start to write "<<endl;
        t2lastmuonWithrpc[i]->Write();
        t2lastmuonWithoutrpc[i]->Write();
        FnProEWithrpc[i]->Write();
        FnProEWithoutrpc[i]->Write();
        pxy[i]->Write();
        prz[i]->Write();
        pE[i]->Write();
        dxy[i]->Write();
        drz[i]->Write();
        dE[i]->Write();
        intervalT[i]->Write();
        pxy0[i]->Write();
        prz0[i]->Write();
        pE0[i]->Write();
        dxy0[i]->Write();
        drz0[i]->Write();
        dE0[i]->Write();
        intervalT0[i]->Write();
        std::cout<<"done "<<endl;
        std::cout<<""<<endl;
        std::cout<<"AD"<<i+1<<" Num of IBD                          : ["<<Num2[i]+Num1[i]<<" ("<<((double)(Num2[i]+Num1[i]) -(double)(Numo[i]+Num1[i]) )/(double)(Numo[i]+Num1[i])<<")]  "<<(Num2[i]+Num1[i])/livetime[site][i] <<" /day"<<endl;
        std::cout<<"    Num2: "<<Num2[i]<<" Num1: "<<Num1[i]<<endl;
        std::cout<<"AD"<<i+1<<" Num of IBD (without using rpc veto) : ["<<Numo[i]+Num1[i]<<"]  "<<(Numo[i]+Num1[i])/livetime0[site][i] <<" /day"<<endl;
        std::cout<<"    Numo: "<<Numo[i]<<" Num1: "<<Num1[i]<<endl;
        std::cout<<" "<<endl;
    }

    //fn

    for( int i=1 ; i<=200 ; i++ )
    {
        double nEvt = tFnProEWithrpc->GetBinContent(i);
        tFnProEWithrpc->SetBinError(i, sqrt(nEvt*BinWt[i]));
        double nEvt0 = tFnProEWithoutrpc->GetBinContent(i);
        double nEvt0UniX = tFnProEWithoutrpcUniX->GetBinContent(i);
        tFnProEWithoutrpc->SetBinError(i, sqrt(nEvt0*BinWt[i]));
        tFnProEWithoutrpcUniX->SetBinError(i, sqrt(nEvt0UniX));
    }

    fitHighEdge=50.;
    scaledHighEdge=99.;//here,99.6 should be less than tFnProEWithoutrpcUniX's maxmium 99.7,because >99.7 there are many events .
    if( site==2 )
    {
        reBinNum=3;
    } else
    {
        reBinNum=4;
    }
    TFile* f_ows=new TFile("AnaFNEH1_ows.root","read");
    if( f_ows->IsZombie() )
    {
        std::cout<<"Can't find AnaFNEH1_ows.root ... "<<endl;
    }
    TH1F* h_ows=(TH1F*)f_ows->Get("h1");
    if( !h_ows )
    {
        std::cout<<"Can't find h1 in AnaFNEH1_ows.root ... "<<endl;
    }
    TH1F h1_ows=*h_ows;
    double tnum=tFnProEWithoutrpcUniX->Integral(tFnProEWithoutrpcUniX->FindBin(12.),tFnProEWithoutrpcUniX->FindBin(scaledHighEdge));
    double hnum_ows=h_ows->Integral(h_ows->FindBin(12.),h_ows->FindBin(scaledHighEdge));
    std::cout<<"tnum  : "<<tnum<<endl;
    std::cout<<"hnum_ows  : "<<hnum_ows<<endl;
    h1_ows.Scale(tnum/hnum_ows);

    TFile* f_rpc=new TFile("AnaFNEH1_rpc.root","read");
    if( f_rpc->IsZombie() )
    {
        std::cout<<"Can't find AnaFNEH1_rpc.root ... "<<endl;
    }
    TH1F* h_rpc=(TH1F*)f_rpc->Get("h1");
    if( !h_rpc )
    {
        std::cout<<"Can't find h1 in AnaFNEH1_rpc.root ... "<<endl;
    }
    TH1F h1_rpc=*h_rpc;
    double hnum_rpc=h_rpc->Integral(h_rpc->FindBin(12.),h_rpc->FindBin(scaledHighEdge));
    std::cout<<"hnum_rpc  : "<<hnum_rpc<<endl;
    h1_rpc.Scale(tnum/hnum_rpc);

    TFile* f_iws_MC=new TFile("AnaFNEH1_iws_MC.root","read");
    if( f_iws_MC->IsZombie() )
    {
        std::cout<<"Can't find AnaFNEH1_iws_MC.root ... "<<endl;
    }
    TH1F* h_iws_MC=(TH1F*)f_iws_MC->Get("h1");
    if( !h_iws_MC )
    {
        std::cout<<"Can't find h1 in AnaFNEH1_iws_MC.root ... "<<endl;
    }
    TH1F h1_iws_MC=*h_iws_MC;
    double hnum_iws_MC=h_iws_MC->Integral(h_iws_MC->FindBin(12.),h_iws_MC->FindBin(scaledHighEdge));
    std::cout<<"hnum_iws_MC  : "<<hnum_iws_MC<<endl;
    h1_iws_MC.Scale(tnum/hnum_iws_MC);

    TFile* f_mc=new TFile("AnafilefastMCEH1_noIWS.root","read");
    if( f_mc->IsZombie() )
    {
        std::cout<<"Can't find AnafilefastMCEH1_noIWS.root ... "<<endl;
    }
    TH1F* h_mc=(TH1F*)f_mc->Get("ad1H");
    if( !h_mc )
    {
        std::cout<<"Can't find h1 in AnafilefastMCEH1_noIWS.root ... "<<endl;
    }
    TH1F h1_mc=*h_mc;
    double hnum_mc=h_mc->Integral(h_mc->FindBin(12.),h_mc->FindBin(scaledHighEdge));
    std::cout<<"hnum_mc  : "<<hnum_mc<<endl;
    h1_mc.Scale(tnum/hnum_mc);

    TH1F* h1Q=new TH1F("h1Q","h1Q",174,12.7,99.7);
    TH1F* tFnProEWithoutrpcUniXQ=new TH1F("tFnProEWithoutrpcUniXQ","tFnProEWithoutrpcUniXQ",174,12.7,99.7);
    for( int i=1 ; i<=174 ; i++ )
    {
        h1Q->SetBinContent(i,h1_ows.GetBinContent(i+24));
        tFnProEWithoutrpcUniXQ->SetBinContent(i,tFnProEWithoutrpcUniX->GetBinContent(i+24));
    }
    double res[174], x[174];
    double prob=0.;
    double _chi2=0.;
    int _ndf;
    int _isgood;
    prob=tFnProEWithoutrpcUniXQ->Chi2TestX(h1Q,_chi2,_ndf,_isgood,"UU",res);
    //Graph for Residuals
    for (Int_t i=0; i<174; i++) x[i]= 12.7+i*0.5+0.25;
    TGraph *resgr = new TGraph(174,x,res);
    resgr->GetXaxis()->SetRangeUser(12.7,99.7);
    resgr->GetYaxis()->SetRangeUser(-3.5,3.5);
    resgr->GetYaxis()->SetTitle("Normalized Residuals");
    resgr->SetMarkerStyle(21);
    resgr->SetMarkerColor(2);
    resgr->SetMarkerSize(.9);
    resgr->SetTitle("Normalized Residuals");

    //Quantile-Quantile plot
    TF1 *ff = new TF1("ff","TMath::Gaus(x,0,1)",-10,10);
    TGraphQQ *qqplot = new TGraphQQ(174,res,ff);
    qqplot->SetMarkerStyle(20);
    qqplot->SetMarkerColor(2);
    qqplot->SetMarkerSize(.9);
    qqplot->SetTitle("Q-Q plot of Normalized Residuals");

    //create Canvas
    TCanvas *c1 = new TCanvas("c1","Chistat Plot",10,10,700,600);
    c1->Divide(1,2);
    c1->cd(1);
    string probStr=Form("chi2/ndf : %f/%i   p-value : %f",_chi2,_ndf,prob);
    TPaveLabel* label1=new TPaveLabel(45,3,90,4,probStr.c_str());
    resgr->Draw("APL");
    label1->Draw();
    c1->cd(2);
    qqplot->Draw("AP");
    c1->SaveAs(Form("%s/%s%sfnSpecNorRes.eps",dataVer.c_str(),dataVer.c_str(),siteStr.c_str()));

    file->cd();
    TCanvas* c2=new TCanvas(Form("%sfnSpec",dataVer.c_str()),"c2",600,400);
    //use rpc veto
    //pol0
    cout<<" "<<endl;
    cout<<" "<<endl;
    cout<<">>> log bin ... 12("<<tFnProEWithrpc->FindBin(12.) <<") ~ "<<fitHighEdge <<"("<<tFnProEWithrpc->FindBin(fitHighEdge) <<")"<<endl;
    TF1* f0= new TF1("f0","pol0",12.,fitHighEdge);
    tFnProEWithrpc->Fit(f0,"R+");
    double par0,ipar0;
    f0->GetParameters(&par0);
    ipar0=f0->GetParError(0);
    double NFn0=0.;
    double iNFnsquare0=0.;
    int binlow=tFnProEWithrpc->FindBin(0.7);
    int binup=tFnProEWithrpc->FindBin(12.);
    std::cout<<"binlow  : "<<binlow<<endl;
    std::cout<<"binup  : "<<binup<<endl;
    double bincenter=0.;

    for( int j=binlow ; j<=binup ; j++ )
    {
        bincenter=tFnProEWithrpc->GetBinCenter(j);
        NFn0+=par0/BinWt[j];
        iNFnsquare0+=ipar0*ipar0/(BinWt[j]*BinWt[j]);
    }
    //NFn0=tFnProEWithrpc->Integral(tFnProEWithrpc->FindBin(12),tFnProEWithrpc->FindBin(100),"width")*(12-0.7)/(100-12);
    //iNFnsquare0=sqrt(tFnProEWithrpc->Integral(tFnProEWithrpc->FindBin(12),tFnProEWithrpc->FindBin(100),"width"))*(12-0.7)/(100-12);
    ////NFn00=f00->Integral(0.7,12.0);
    //pol1
    TF1* f1= new TF1("f1","pol1",12.,fitHighEdge);
    tFnProEWithrpc->Fit(f1,"R");
    double par1[2],ipar1[2];
    f1->GetParameters(&par1[0]);
    ipar1[0]=f1->GetParError(0);
    ipar1[1]=f1->GetParError(1);
    double NFn1=0.;
    double iNFnsquare1=0.;
    for( int j=binlow ; j<=binup ; j++ )
    {
        bincenter=tFnProEWithrpc->GetBinCenter(j);
        NFn1+=(par1[0]+par1[1]*bincenter)/BinWt[j];
        iNFnsquare1+=(bincenter*bincenter*ipar1[1]*ipar1[1]+ipar1[0]*ipar1[0])/(BinWt[j]*BinWt[j]);
    }
    //NFn1=(par1[0]*2+12.7*par1[1])*(12-0.7)/2;
    //iNFnsquare1=(ipar1[0]*2+12.7*ipar1[1])*(12-0.7)/2; //this is not appropriate ,from zhangfh.

    //NFn1=f1->Integral(0.7,12.0);

    std::cout<<"NFn0  : "<<NFn0<<endl;
    std::cout<<"NFn1  : "<<NFn1<<endl;
    //double tNum1=tFnProEWithrpcUniX->Integral(tFnProEWithrpcUniX->FindBin(12.),tFnProEWithrpcUniX->FindBin(60));
    //double hNum1_ows=h_ows->Integral(h_ows->FindBin(12.),h_ows->FindBin(60));
    //double fnNum1=h_ows->Integral(h_ows->FindBin(0.7),h_ows->FindBin(12.));
    //fnNum1=fnNum1*tNum1/hNum1_ows;
    //std::cout<<"NFn10-fnNum  : "<<NFn10-fnNum <<endl;
    //std::cout<<"NFn00-fnNum  : "<<NFn00-fnNum<<endl;
    //double fnNum1Err=sqrt(fnNum1+max(abs(NFn1-fnNum1),abs(NFn0-fnNum1))*max(abs(NFn1-fnNum1),abs(NFn0-fnNum1)));
    //std::cout<<"NFn0  : "<<NFn0<<" NFn0-fnNum1  : "<<NFn0-fnNum1<<" (NFn0-fnNum1)/fnNum1 : "<<(NFn0-fnNum1)/fnNum1<<endl;
    //std::cout<<"NFn1  : "<<NFn1<<" NFn1-fnNum1  : "<<NFn1-fnNum1<<" (NFn1-fnNum1)/fnNum1 : "<<(NFn1-fnNum1)/fnNum1<<endl;
    //std::cout<<"fnNum1  : "<<fnNum1<<" stat.err/fnNum1 : "<<sqrt(fnNum1)/fnNum1<<" sys.err/fnNum1 : "<<max(abs(NFn1-fnNum1),abs(NFn0-fnNum1))/fnNum1<<" total.err/fnNum1 : "<<fnNum1Err/fnNum1<<endl;
    //NFn0=f0->Integral(0.7,12.0);

    //without using rpc veto
    cout<<" "<<endl;
    cout<<" "<<endl;
    double tNum=tFnProEWithoutrpcUniX->Integral(tFnProEWithoutrpcUniX->FindBin(12.),tFnProEWithoutrpcUniX->FindBin(scaledHighEdge));
    double hNum_ows=h_ows->Integral(h_ows->FindBin(12.),h_ows->FindBin(scaledHighEdge));
    double fnNum_ows=h_ows->Integral(h_ows->FindBin(0.7),h_ows->FindBin(12.));
    fnNum_ows=fnNum_ows*tNum/hNum_ows;

    double hNum_rpc=h_rpc->Integral(h_rpc->FindBin(12.),h_rpc->FindBin(scaledHighEdge));
    double fnNum_rpc=h_rpc->Integral(h_rpc->FindBin(0.7),h_rpc->FindBin(12.));
    fnNum_rpc=fnNum_rpc*tNum/hNum_rpc;

    double hNum_iws_MC=h_iws_MC->Integral(h_iws_MC->FindBin(12.),h_iws_MC->FindBin(scaledHighEdge));
    double fnNum_iws_MC=h_iws_MC->Integral(h_iws_MC->FindBin(0.7),h_iws_MC->FindBin(12.));
    fnNum_iws_MC=fnNum_iws_MC*tNum/hNum_iws_MC;
    
    double hNum_mc=h_mc->Integral(h_mc->FindBin(12.),h_mc->FindBin(scaledHighEdge));
    double fnNum_mc=h_mc->Integral(h_mc->FindBin(0.7),h_mc->FindBin(12.));
    fnNum_mc=fnNum_mc*tNum/hNum_mc;

    h1_ows.Rebin(reBinNum);
    h1_rpc.Rebin(reBinNum);
    h1_iws_MC.Rebin(reBinNum);
    h1_mc.Rebin(reBinNum);
    tFnProEWithoutrpcUniX->Rebin(reBinNum);
    cout<<">>> uniform bin ... 12("<<tFnProEWithoutrpcUniX->FindBin(12.) <<") ~ "<<fitHighEdge <<"("<<tFnProEWithoutrpcUniX->FindBin(fitHighEdge) <<")"<<endl;
    int binlowUniX=tFnProEWithoutrpcUniX->FindBin(0.7);
    int binupUniX=tFnProEWithoutrpcUniX->FindBin(12.);
    //pol00
    TF1* f00= new TF1("f00","pol0",12.,fitHighEdge);
    TFitResultPtr r0=tFnProEWithoutrpcUniX->Fit(f00,"0R+S");
    //TFitResultPtr r0;
    double chi20=r0->Chi2();
    double ndf0=r0->Ndf();
    double par00,ipar00;
    f00->GetParameters(&par00);
    ipar00=f00->GetParError(0);
    double NFn00=0.;
    double iNFnsquare00=0.;
    for( int j=binlowUniX ; j<=binupUniX ; j++ )
    {
        bincenter=tFnProEWithoutrpcUniX->GetBinCenter(j);
        NFn00+=par00;
        iNFnsquare00+=ipar00*ipar00;
    }

    //NFn00=tFnProEWithoutrpc->Integral(tFnProEWithoutrpc->FindBin(12),tFnProEWithoutrpc->FindBin(100),"width")*(12-0.7)/(100-12);
    //iNFnsquare00=sqrt(tFnProEWithoutrpc->Integral(tFnProEWithoutrpc->FindBin(12),tFnProEWithoutrpc->FindBin(100),"width"))*(12-0.7)/(100-12);
    //pol10
    TF1* f10= new TF1("f10","pol1",12.,fitHighEdge);
    TFitResultPtr r=tFnProEWithoutrpcUniX->Fit(f10,"0RS");
    //TFitResultPtr r;
    double chi2=r->Chi2();
    double ndf=r->Ndf();
    double par10[2],ipar10[2];
    f10->GetParameters(&par10[0]);
    ipar10[0]=f10->GetParError(0);
    ipar10[1]=f10->GetParError(1);
    double NFn10=0.;
    double iNFnsquare10=0.;
    for( int j=binlowUniX ; j<=binupUniX ; j++ )
    {
        bincenter=tFnProEWithoutrpcUniX->GetBinCenter(j);
        NFn10+=(par10[0]+par10[1]*bincenter);
        iNFnsquare10+=(bincenter*bincenter*ipar10[1]*ipar10[1]+ipar10[0]*ipar10[0]);
    }
    //NFn10=(par10[0]*2+12.7*par10[1])*(12-0.7)/2;
    //iNFnsquare10=(ipar10[0]*2+12.7*ipar10[1])*(12-0.7)/2;
    //
    //NFn00=f00->Integral(0.7,12.0);
    //double tNum=tFnProEWithoutrpcUniX->Integral(tFnProEWithoutrpcUniX->FindBin(12.),tFnProEWithoutrpcUniX->FindBin(scaledHighEdge));
    //double hNum_ows=h_ows->Integral(h_ows->FindBin(12.),h_ows->FindBin(scaledHighEdge));
    //double fnNumErr_ows=sqrt(fnNum_ows+max(abs(NFn10-fnNum_ows),abs(NFn00-fnNum_ows))*max(abs(NFn10-fnNum_ows),abs(NFn00-fnNum_ows)));
    //std::cout<<"NFn00  : "<<NFn00<<" NFn00-fnNum_ows  : "<<NFn00-fnNum_ows<<" (NFn00-fnNum_ows)/fnNum_ows : "<<(NFn00-fnNum_ows)/fnNum_ows<<endl;
    //std::cout<<"NFn10  : "<<NFn10<<" NFn10-fnNum_ows  : "<<NFn10-fnNum_ows<<" (NFn10-fnNum_ows)/fnNum_ows : "<<(NFn10-fnNum_ows)/fnNum_ows<<endl;
    //std::cout<<"fnNum_ows  : "<<fnNum_ows<<" stat.err/fnNum_ows : "<<sqrt(fnNum_ows)/fnNum_ows<<" sys.err/fnNum_ows : "<<max(abs(NFn10-fnNum_ows),abs(NFn00-fnNum_ows))/fnNum_ows<<" total.err/fnNum_ows : "<<fnNumErr_ows/fnNum_ows<<endl;

    double fnNum=(fnNum_ows+fnNum_rpc)/2;
    //double fnNum=(fnNum_ows+fnNum_rpc+fnNum_iws_MC)/3;
    multimap<double,string> fnNumMap;
    fnNumMap.insert(make_pair((double)fnNum,"*fnNum"));
    fnNumMap.insert(make_pair((double)fnNum_ows,"fnNum_ows"));
    fnNumMap.insert(make_pair((double)fnNum_rpc,"fnNum_rpc"));
    //fnNumMap.insert(make_pair((double)fnNum_iws_MC,"fnNum_iws_MC"));
    fnNumMap.insert(make_pair((double)NFn00,"NFn00"));
    fnNumMap.insert(make_pair((double)NFn10,"NFn10"));
    int mapSize=fnNumMap.size();
    int mapi=1;
    double fnNumErr=0.;
    for(multimap<double,string>::iterator it=fnNumMap.begin(); it!=fnNumMap.end() ; it++ )
    {
        if( mapi==1||mapi==mapSize )
        {
            fnNumErr=fnNumErr<abs((it->first-fnNum)/fnNum)?abs((it->first-fnNum)/fnNum):fnNumErr;
        }
        nameStr=Form("  %10s  %5.2f  %3.1f%%",it->second.c_str(),it->first,(it->first-fnNum)/fnNum*100); 
        cout<<nameStr<<endl;
    }
    int fnNumBinNum=tFnProEWithoutrpcUniX->FindBin(12.)-tFnProEWithoutrpcUniX->FindBin(0.7)+1;
    //double fnNumLow=(1-fnNumErr)*fnNum/fnNumBinNum;
    //double fnNumHigh=(1+fnNumErr)*fnNum/fnNumBinNum;
    TGraph* gLow=new TGraph();
    TGraph* gHigh=new TGraph();
    //TGraph *grshade = new TGraph();

    for( int i=0 ; i<fnNumBinNum ; i++ )
    {
        double yLow=(h1_ows.GetBinContent(h1_ows.FindBin(0.7)+i)+h1_rpc.GetBinContent(h1_ows.FindBin(0.7)+i))/2*(1-fnNumErr);
        double xLow=h1_ows.GetBinCenter(h1_ows.FindBin(0.7)+i);
        double yHigh=(h1_ows.GetBinContent(h1_ows.FindBin(0.7)+i)+h1_rpc.GetBinContent(h1_ows.FindBin(0.7)+i))/2*(1+fnNumErr);
        double xHigh=h1_ows.GetBinCenter(h1_ows.FindBin(0.7)+i);
        cout<<"xLow yLow xHigh yHigh  : "<<xLow<<" "<<yLow<<" "<<xHigh<<" "<<yHigh<<endl;
        gLow->SetPoint(i,xLow,yLow);
        gHigh->SetPoint(i,xHigh,yHigh);
        //grshade->SetPoint(i,xLow,yHigh);
        //grshade->SetPoint(fnNumBinNum+i,xLow,yLow);
    }
    //grshade->SetFillStyle(3013);
    //grshade->SetFillColor(16);

    gLow->SetLineColor(kGreen+3);
    gHigh->SetLineColor(kGreen+3);
    

    double totallivetime0=0.;
    for( int i=0 ; i<ADNum ; i++ )
    {
        totallivetime0+=livetime0[site][i];
    }
    double totallivetime=0.;
    for( int i=0 ; i<ADNum ; i++ )
    {
        totallivetime+=livetime[site][i];
    }
    cout<<"Rate syst.err  : "<<fnNum*fnNumErr/totallivetime0<<endl;
    cout<<"Rate stat.err  : "<<sqrt(fnNum)/totallivetime0<<endl;
    
    fnNumErr=sqrt(fnNum+(fnNum*fnNumErr)*(fnNum*fnNumErr));
    cout<<"total persent  : "<<Form(" %.2f%% ",fnNumErr/fnNum*100)<<endl;

    std::cout<<"	total0  : "<<fnNum<<" +- "<<fnNumErr\
        <<"   B/S : "<<fnNum/(tNumo+tNum1) \
        <<" rate : "<<fnNum/totallivetime0<<" +- "<< fnNumErr/totallivetime0<<" /day"<<endl;
    //tFnProEWithoutrpcUniX->SetLineColor();
    tFnProEWithoutrpcUniX->SetTitle(Form("%s",siteStr.c_str()));
    tFnProEWithoutrpcUniX->SetStats(kFALSE);
    tFnProEWithoutrpcUniX->GetXaxis()->SetTitle("Prompt signal energy (MeV)");
    tFnProEWithoutrpcUniX->GetYaxis()->SetTitle(Form("Events / %.1f MeV",0.5*reBinNum));
    tFnProEWithoutrpcUniX->SetLineWidth(2);
    tFnProEWithoutrpcUniX->Draw("e");
    h1_ows.SetLineColor(kRed);
    h1_ows.Draw("same");
    h1_rpc.SetLineColor(kBlue);
    h1_rpc.Draw("same");
    h1_iws_MC.SetLineColor(kGreen);
    //h1_iws_MC.Draw("same");
    //TLine lhigh(0.7,fnNumHigh,12.,fnNumHigh);
    //TLine llow(0.7,fnNumLow,12.,fnNumLow);
    //lhigh.SetLineColor(kGreen);
    //llow.SetLineColor(kGreen);
    //lhigh.Draw("same");
    //llow.Draw("same");
    gLow->Draw("same");
    gHigh->Draw("same");
    //grshade->Draw("samef");
    gPad->SetLogy();
        TLegend *legend=new TLegend(.55,.55,.89,.89);
        legend->AddEntry(tFnProEWithoutrpcUniX,"Data","lp");
        legend->AddEntry(&h1_ows,"OWS-tagged","lp");
        legend->AddEntry(&h1_rpc,"RPC-only-tagged","lp");
        //legend->AddEntry(&h1_iws_MC,"IWS-tagged","lp");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->Draw();
    c2->SaveAs(Form("%s/%s%sfnSpecForNuFact2014.eps",dataVer.c_str(),dataVer.c_str(),siteStr.c_str()));
    //c2->Write();
    tFnProEWithrpc->Write();
    tFnProEWithoutrpc->Write();
    tFnProEWithoutrpcUniX->Write();
    //h_ows->Write();
    //delete h_ows;
    //delete c;
    //delete f;
    //file->Close();	

    //write to result_temp.txt
    int ii=0;
    if( site==0 )
    {
        ii=0;
    } else if( site==1 )
    {
        ii=2;
    }else if( site==2 )
    {
        ii=4;
    }else
    {
        std::cout<<"error  : wrong siteNum"<<endl;    
    }
    for( int i=0 ; i<ADNum ; i++ )
    {
        t_result[1][i+ii]=Numo[i]+Num1[i];
        t_result[13][i+ii]=fnNum;
        t_result[14][i+ii]=fnNumErr;
        t_result[15][i+ii]=fnNum/(tNumo+tNum1);
        t_result[16][i+ii]=UpperNum[i]/LowerNum[i];
        t_result[17][i+ii]=fnNum/totallivetime0;
        t_result[18][i+ii]=fnNumErr/totallivetime0;
        /*
           t_result[24][i+ii]=Num2[i]+Num1[i];
           if( Num2[i]+Num1[i]==0 )
           {
           t_result[25][i+ii]=0;
           } else
           {
           t_result[25][i+ii]=((double)(Num2[i]+Num1[i])-(double)(Numo[i]+Num1[i]) )/(double)(Numo[i]+Num1[i]);
           }
           */
        t_result[24][i+ii]=UpperNum[i];
        t_result[25][i+ii]=LowerNum[i];
        t_result[44][i+ii]=(NFn0+NFn1)/2;
        t_result[45][i+ii]=(sqrt((iNFnsquare0+iNFnsquare1)/4+(NFn1-NFn0)*(NFn1-NFn0)));
        t_result[46][i+ii]=((NFn0+NFn1)/2 -(NFn00+NFn10)/2 )/((NFn00+NFn10)/2);
        t_result[47][i+ii]=((NFn0+NFn1)/(tNum2+tNum1))/2;
        t_result[48][i+ii]=0;
        t_result[49][i+ii]=(((NFn0+NFn1)/(tNum2+tNum1))/2 -(NFn00+NFn10)/2/(tNumo+tNum1) )/((NFn00+NFn10)/2/(tNumo+tNum1));
        t_result[50][i+ii]=(NFn0+NFn1)/2/totallivetime;
        t_result[51][i+ii]=(sqrt((iNFnsquare0+iNFnsquare1)/4+(NFn1-NFn0)*(NFn1-NFn0)))/totallivetime;
        t_result[52][i+ii]=( (NFn0+NFn1)/2/totallivetime-(NFn00+NFn10)/2/totallivetime0 )/((NFn00+NFn10)/2/totallivetime0);
    }

    ofstream outfile;
    string outfileName=dataVer+"/"+"result_temp_"+dataVer+".txt";
    outfile.open(outfileName.c_str());
    int outLineNum=0;
    while( outLineNum<53 )
    {
        for( int i=0 ; i<8 ; i++ )
        {
            outfile<< setw(13) <<t_result[outLineNum][i];
            //outfile<< setw(13) <<outStr[outLineNum][i];
            if( i==7 )
            {
                outfile<<endl;
            }
        }
        outLineNum++;
    }
    outfile.close();

    TString outStr[54][9];
    ifstream infile4out; 
    string infileoutName=dataVer+"/"+"result_temp_"+dataVer+".txt";
    infile4out.open(infileoutName.c_str(),ios::in);
    int inLineNum4out=0;
    while(inLineNum4out<53)
    {
        for( int j=1 ; j<9 ; j++ )
        {
            infile4out>>outStr[inLineNum4out][j];
        }
        inLineNum4out++;
    }
    infile4out.close();

    outStr[0][0]="====NoRpc====";//0
    outStr[1][0]="ibdNum";//1
    outStr[2][0]="daqTime";//2
    outStr[3][0]="epsi_mu";//3 
    outStr[4][0]="epsi_multi";//4
    outStr[5][0]="accNum";//5
    outStr[6][0]="+-";//6
    outStr[7][0]="accRate";//7
    outStr[8][0]="+-";//8
    outStr[9][0]="accNum_600";//9
    outStr[10][0]="+-";//10
    outStr[11][0]="accRate_600";//11
    outStr[12][0]="+-";//12
    outStr[13][0]="fnNum";//13
    outStr[14][0]="+-";//14
    outStr[15][0]="fnBs";//15
    outStr[16][0]="U/L";//16
    outStr[17][0]="fnRate";//17
    outStr[18][0]="+-";//18

    outStr[19][0]="LiNum";//18
    outStr[20][0]="+-";//18
    outStr[21][0]="LiRate";//18
    outStr[22][0]="+-";//18

    outStr[23][0]="===WithRpc===";//19
    //outStr[24][0]="ibdNum";//20
    //outStr[25][0]="(-)";//21
    outStr[24][0]="UpperNum";//20
    outStr[25][0]="LowerNum";//21
    outStr[26][0]="daqTime";//22
    outStr[27][0]="(-)";//23
    outStr[28][0]="epsi_mu";//24
    outStr[29][0]="(-)";//25
    outStr[30][0]="epsi_multi";//26
    outStr[31][0]="(-)";//27
    outStr[32][0]="accNum";//28
    outStr[33][0]="+-";//29
    outStr[34][0]="(-)";//30
    outStr[35][0]="accRate";//31
    outStr[36][0]="+-";//32
    outStr[37][0]="(-)";//33
    outStr[38][0]="accNum_600";//34
    outStr[39][0]="+-";//35
    outStr[40][0]="(-)";//36
    outStr[41][0]="accRate_600";//37
    outStr[42][0]="+-";//38
    outStr[43][0]="(-)";//39
    outStr[44][0]="fnNum";//40
    outStr[45][0]="+-";//41
    outStr[46][0]="(-)";//42
    outStr[47][0]="fnBs";//43
    outStr[48][0]="+-";//44
    outStr[49][0]="(-)";//45
    outStr[50][0]="fnRate";//46
    outStr[51][0]="+-";//47
    outStr[52][0]="(-)";//48

    ofstream outfile4out;
    string outfileoutName=dataVer+"/"+"result_"+dataVer+".txt";
    outfile4out.open(outfileoutName.c_str());
    int outLineNum4out=0;
    while( outLineNum4out<53 )
    {
        for( int i=0 ; i<9 ; i++ )
        {
            outfile4out<< setw(13) <<outStr[outLineNum4out][i];
            if( i==8 )
            {
                outfile4out<<endl;
            }
        }
        outLineNum4out++;
    }
    outfile4out.close();
}
