#include <iostream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMath.h>
#include "TStyle.h"
#include "TCanvas.h"

using namespace std;

vector <double>  clusterize(unsigned short *adc, vector <int> &highchans, int nclust, int highThresh, int lowThresh, int event);

void fiberCluster(int runnumber, int board, int highThresh, int lowThresh){

  TString term_cmd = "mkdir -pv "+TString::Format("Plots/run%d",(int)runnumber);
  system(term_cmd.Data());
  
  
  TString term_cmd2 = "mkdir -pv "+TString::Format("Results/run%d",(int)runnumber);
  system(term_cmd2.Data());

  //Initialize results vector
  vector <double> results;
  results.clear();
  results.assign(4,-1);

  //Initialize highest channels vector
  vector <int> highchans;
  results.clear();
  
  int doped=1;
  int dosigr=1;
  int event=0;
  int nclust=0;

  TH1F * h_adc = new TH1F("h_adc", "h_adc", 150, -0.5, 150.5);
  h_adc->GetXaxis()->SetTitle("Arbitrary units");
  
  TH1F * h_cog = new TH1F("h_cog", "h_cog", 33, -0.5, 32.5);
  h_cog->GetXaxis()->SetTitle("cog (mm)");
  
  TH1F * h_cogHigh = new TH1F("h_cogHigh", "h_cogHigh", 33, -0.5, 32.5);
  h_cogHigh->GetXaxis()->SetTitle("cog (mm)");

  TH1F * h_clust = new TH1F("h_clust#", "h_clust#", 11, -0.5, 10.5);
  h_clust->GetXaxis()->SetTitle("clusters number");

  TH1F * h_width = new TH1F("h_width", "h_clust_width", 11, -0.5, 10.5);
  h_width->GetXaxis()->SetTitle("clusters width");

  TH2F * h_2d = new TH2F("h_2d", "h_2d", 33, -0.5, 32.5, 100, 0, 50);
  h_2d->GetXaxis()->SetTitle("cog (mm)");
  h_2d->GetYaxis()->SetTitle("ADC");

  
  TH1F * histo_raw_sigmas = new TH1F("histo_raw_sigmas", "histo_raw_sigmas", 150, 0, 50);


  TString filename = TString::Format("test_herd_20170517_run_%d",(int)runnumber)+TString::Format("_brd%d",(int) board)+".root";
  TString filenamecal = TString::Format("test_herd_20170517_run_calib_%d",(int)runnumber)+TString::Format("_brd%d",(int) board)+".root";
  TString outfile = TString::Format("Results/run%d",(int)runnumber)+TString::Format("/test_herd_20170517_run_%d",(int)runnumber)+TString::Format("_brd%d",(int) board)+"_results.root";

  TString plotdirectory = TString::Format("Plots/run%d",(int)runnumber);


  TFile* f = new TFile(filename,"READ");
  TFile* cal = new TFile(filenamecal,"READ");
  
  TFile* out = new TFile(outfile, "recreate");

  TTree *CLUST = new TTree("FiberClusters","T");
  CLUST->Branch("eventnum", &event);
  CLUST->Branch("nclusters", &nclust);
  CLUST->Branch("FiberClusters", &results);

  cout << endl;
  f->ls();
  cout << endl;

  cout << endl;
  cal->ls();
  cout << endl;


  unsigned short adc[128], padc[128];
  
  TTree* t = (TTree*)f->Get("T");
  t->SetBranchAddress("adc", &adc);
  
  TTree* tc= (TTree*)cal->Get("T");
  tc->SetBranchAddress("adc", &padc);

  int nevents = t->GetEntries();
  int neventscal = tc->GetEntries();
  
  cout << endl;
  cout << "Number of events " << nevents << endl;
  cout << "Number of cal events " << neventscal << endl;


  //Initialize pedestals and sigmas vectors
  vector <double> pedestals, raw_sigmas;
  pedestals.clear();
  pedestals.assign(128,0);
  raw_sigmas.clear();
  raw_sigmas.assign(128,0);  

  //Create graph for pedestals and sigmas
  TGraph *gr_ped=new TGraph();
  TGraph *gr_raw_sigmas=new TGraph();
  gr_ped->SetMarkerStyle(23);
  gr_ped->SetMarkerColor(kBlue);
  gr_ped->SetLineColor(kBlue);
  gr_raw_sigmas->SetMarkerStyle(23);
  gr_raw_sigmas->SetMarkerColor(kBlue);
  gr_raw_sigmas->SetLineColor(kBlue);

  TGraph *gr_event = new TGraph();
  
  //Using events from calibration file to calculate pedestals 
  while(event < neventscal && (doped || dosigr)){
        tc->GetEntry(event);
	
	for(int channel=0; channel < 128; channel++){
	  if (event<1024) {
	    pedestals.at(channel) +=padc[channel];
	  }
	  if ((event==1024) && doped) {
	    // now we finalize the pedestals
	    for (int ch=0; ch<(int)pedestals.size(); ch++) {
	      pedestals.at(ch)/=1024.;
	      gr_ped->SetPoint(gr_ped->GetN(),ch,pedestals.at(ch));
	     }
	     doped=0; 
	    
	  }

	  //Using 2048 events from calibration file to calculate raw sigmas
	  if (event>1024 && event<3072) {
	    raw_sigmas.at(channel) += pow((padc[channel]-pedestals.at(channel)), 2);
	  }
	  if ((event==2048) && dosigr) {
	    for (int ch=0; ch<(int)raw_sigmas.size(); ch++) {
	      raw_sigmas.at(ch)/=2048.;
	      raw_sigmas.at(ch)=sqrt(raw_sigmas.at(ch));
	      //cout << "SIGMA " << raw_sigmas.at(channel) << endl;
	      gr_raw_sigmas->SetPoint(gr_raw_sigmas->GetN(),ch,raw_sigmas.at(ch));
	      histo_raw_sigmas->Fill(raw_sigmas.at(ch));
	    }  
	    dosigr=0;
	  }
	}	
	event++;
  }
  
      // for(int k=0; k < 128; k++){
      // 	cout << raw_sigmas.at(k) << endl;
      // 	cout << pedestals.at(k) << endl;
      //  }

  
  //Event loop
  for(int event_index=0; event_index < nevents; event_index++){
    t->GetEntry(event_index);

    highchans.clear();
    nclust = 0;

    //Pedestal subtraction
    
      // for(int k=0; k < 128; k++){
      // 	adc[k]=adc[k]-pedestals.at(k);
      //  }
    

    for(int k=0; k < 64; k++){
       if(adc[k] > pedestals.at(k)){
       	adc[k]=adc[k]-pedestals.at(k);
       } else{
       	adc[k]=0;
       }
    }

    for(int k=64; k < 128; k++){
      // adc[k]=0;
      if(adc[k] > pedestals.at(k)){
      adc[k]=adc[k]-pedestals.at(k);
      } else{
      adc[k]=0;
      }
    }

      for(int k=0; k < 128; k++){
	gr_event->SetPoint(gr_event->GetN(),k*0.25,adc[k]/500);
      }


      //Fill vector of channels over high Threshold
    for(int j=0; j < 128; j++){
      if(adc[j] > highThresh){
	// Only if no neighboring channel was previously marked as highchan
	if(find(highchans.begin(), highchans.end(), j-1) == highchans.end() || find(highchans.begin(), highchans.end(), j+1) == highchans.end()){	    
	  highchans.push_back (adc[j]);
	  highchans.push_back (j);

	}
      } 
    }



    /*
    for(int i=0; i < highchans.size(); i=i+2){
      cout << "@@@@@@@@" << endl;
      cout << "ADC " << highchans.at(i) << endl;
      cout << "Channel "<< highchans.at(i+1) << endl;
    }
    */
        
    nclust = highchans.size()/2;

    //cout << "Number of Clusters " << nclust << endl;

    
    if(nclust > 0){
      h_clust->Fill(nclust);
    //Cluster calculation
      results.clear();
      results=clusterize(adc, highchans, nclust, highThresh, lowThresh, event_index);
      nclust=results.size()/4;
      
      for(int i=0; i < nclust; i++){			   
	//Fill histograms
	h_adc->Fill(results.at(4*i));
	h_cogHigh->Fill(results.at(4*i+1)*0.25);
	h_cog->Fill(results.at(4*i+2));
	h_width->Fill(results.at(4*i+3));
	h_2d->Fill(results.at(4*i+2),results.at(4*i)); 
      }
    }
    //Fill TBranch with results
    CLUST->Fill();
  }


TCanvas *c = new TCanvas("c", "c", 2000, 2000);
 gStyle->SetOptStat(1111111);
 gStyle->SetOptFit(1111);
 
 c->Divide(3,2);
 
 c->cd(1);
 h_cog->Draw();
 h_cog->Write();

 c->cd(2);
 h_cogHigh->Draw();
 h_cogHigh->Write();
 
 c->cd(3);
 h_width->Draw();
 h_width->Write();

 c->cd(4);
 h_adc->Draw();
 h_adc->Write();

 c->cd(5);
 h_clust->Draw();
 h_clust->Write();

 c->cd(6);
 gPad->SetLogz();
 h_2d->Draw("COLZ");
 h_2d->Write("COLZ");


 TCanvas *c2 = new TCanvas("c2", "c2", 2000, 2000);
 gStyle->SetOptStat(1111111);
 gStyle->SetOptFit(1111);
 
 c2->Divide(2,1);

 c2->cd(1);
 TH1F *frame = gPad->DrawFrame(0,0,128,1000);
 frame->SetTitle("Pedestals; channel; pedestals (adc)");
 frame->GetYaxis()->SetTitleOffset(1.55);
 frame->GetXaxis()->SetNdivisions(4,5,1,kFALSE);
 gr_ped->Draw("l");
 gr_ped->Write("l");
 
 c2->cd(2);
 frame = gPad->DrawFrame(0,-10,128,300);
 frame->SetTitle("Raw sigmas; channel; raw sigmas (adc)");
 frame->GetYaxis()->SetTitleOffset(1.55);
 frame->GetXaxis()->SetNdivisions(4,5,1,kFALSE);
 gr_raw_sigmas->Draw("l");
 gr_raw_sigmas->Write("l");

 c->SaveAs(plotdirectory+"/fiberCluster.png");
 c2->SaveAs(plotdirectory+"/fiberCluster2.png");


 out->Write();
}


vector <double>  clusterize(unsigned short *adc, vector <int> &highchans, int nclust, int highThresh, int lowThresh, int event){

  vector <double> results;
  results.clear();  
    
  double cog=-1;
  double cogN=-1;
  double cogD=-1;
  
  int width=0;
  int L,R=1;

  double maxadc=-1;
  int maxchan=-1;
 
  for(int n=0; n < 2* nclust; n=n+2){
   
    maxadc=highchans.at(n);
    maxchan=highchans.at(n+1);
    
    cogN=0;
    cogD=0;

    
    L=0;
    R=0;
    bool goL = true;
    bool goR = true;
        
    //Find strips higher than low Threshold next to high Threshold ones
    while((goL || goR) && maxchan-L>=0 && maxchan+R <=128){
	
	if(adc[maxchan-L] > lowThresh && goL){
	  cogN+=adc[maxchan-L]*(maxchan-L);  
	  cogD+=adc[maxchan-L];
	  L++;
	  goL=true;
	} else{
	  goL=false;
	}
	
	if(adc[maxchan+R] > lowThresh && goR){
	  cogN+=adc[maxchan+R]*(maxchan+R);  
	  cogD+=adc[maxchan+R];
	  R++;
	  goR=true;
	} else{
	  goR=false;
	}
      }
    
    width=TMath::Abs(R-L)+1;

      cog = (cogN/cogD)*0.25;
      maxadc = cogD;

      results.push_back (maxadc/250);
      results.push_back (maxchan);
      results.push_back (cog);
      results.push_back (width);
  }
  
  return results;
}
