#include <iostream>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>

using namespace std;

void daqEncoding1(string protagonist = 0)
{

  if (protagonist == 0)
  {
    cout << "you must provide a file name!" << endl;
    return;
  }

  string input_path = "/home/gigi/Root/FiberTracker/BT201705/DaqFiles/";
  string output_path = "/home/gigi/Root/FiberTracker/BT201705/RootFiles/";

  string daq_ext = ".daq";
  string root_ext = ".root";

  string input_file_name = input_path + protagonist + daq_ext;
  string output_file_name = output_path + protagonist + root_ext;

  FILE *file = fopen(input_file_name.c_str(), "rb");

  if (file == 0)
  {
    cout << "There was an error in the daq file opening!" << endl;
    return;
  }

  TFile *out_file = new TFile(output_file_name.c_str(), "recreate");

  unsigned int val;
  unsigned char trigger, ev_id;
  unsigned short adc1, adc2, adc[128], padc[128];

  TH1F *h_trigger = new TH1F("h_trigger", "h_trigger", 4, 0, 4);
  TH1F *h_ev_id = new TH1F("h_ev_id", "h_ev_id", 65, 0, 65);
  TH1F *h_adc1 = new TH1F("h_adc1", "h_adc1", 4096, 0, 4096);
  TH1F *h_adc2 = new TH1F("h_adc2", "h_adc2", 4096, 0, 4096);

  TH1F *histo_raw_sigmas = new TH1F("histo_raw_sigmas", "histo_raw_sigmas", 3000, 0., 300.);

  TTree *T = new TTree("T", "T");
  T->Branch("adc", adc, "adc[128]/s");

  int cnt = 0, event = 0;
  vector<double> pedestals, raw_sigmas;
  pedestals.clear();
  pedestals.assign(128, 0);
  raw_sigmas.clear();
  raw_sigmas.assign(128, 0);

  TGraph *gr_ped = new TGraph();
  TGraph *gr_raw_sigmas = new TGraph();
  gr_ped->SetMarkerStyle(23);
  gr_ped->SetMarkerColor(kBlue);
  gr_ped->SetLineColor(kBlue);
  gr_raw_sigmas->SetMarkerStyle(23);
  gr_raw_sigmas->SetMarkerColor(kBlue);
  gr_raw_sigmas->SetLineColor(kBlue);

  int doped = 1, dosigr = 1;
  int prev_ev = 0;
  memset(&adc, 0, sizeof(adc));

  do
  {
    if (feof(file))
      break;

    fread(&val, sizeof(val), 1, file);
    //printf("0x%08x\n",val);
    trigger = val >> 30;
    ev_id = (val >> 24) & 0x3f;
    adc2 = (val >> 12) & 0xfff;
    adc1 = val & 0xfff;
    h_trigger->Fill(trigger);
    h_ev_id->Fill(ev_id);
    h_adc2->Fill(adc2);
    h_adc1->Fill(adc1);

    int channel = (2 * (cnt % 64));
    int event = cnt / 64;

    padc[channel] = adc1;
    padc[channel + 1] = adc2;

    if (event < 1024)
    {
      pedestals.at(channel) += adc1;
      pedestals.at(channel + 1) += adc2;
      //cout<<"event = "<<event<<", channel = "<<channel<<endl;
    }
    if ((event == 1024) && doped)
    {
      // now we finalize the pedestals
      for (int ch = 0; ch < (int)pedestals.size(); ch++)
      {
        pedestals.at(ch) /= 1024.;
        gr_ped->SetPoint(gr_ped->GetN(), ch, pedestals.at(ch));
      }
      doped = 0;
    }

    event -= 1025;
    if ((event >= 0) && (event < 2048))
    {
      raw_sigmas.at(channel) += pow((adc1 - pedestals.at(channel)), 2);
      raw_sigmas.at(channel + 1) += pow((adc2 - pedestals.at(channel + 1)), 2);
    }
    if ((event == 2048) && dosigr)
    {
      for (int ch = 0; ch < (int)raw_sigmas.size(); ch++)
      {
        raw_sigmas.at(ch) /= 2048.;
        raw_sigmas.at(ch) = sqrt(raw_sigmas.at(ch));
        gr_raw_sigmas->SetPoint(gr_raw_sigmas->GetN(), ch, raw_sigmas.at(ch));
        histo_raw_sigmas->Fill(raw_sigmas.at(ch));
      }
      dosigr = 0;
    }

    // cout << hex << "0x" << val << endl ;
    // cout << dec << (int)trigger << " ";
    // cout << dec << (int)ev_id << " ";
    // cout << dec << adc2 << " ";
    // cout << dec << adc1 << endl;

    if (channel == 126)
    {

      for (int in = 0; in < 64; in++)
      {
        adc[in] = padc[-in + 63];
      }

      for (int in = 64; in < 128; in++)
      {
        adc[in] = padc[-in + 191];
      }

      for (int in = 0; in < 128; in++)
      {
        padc[in] = pedestals.at(in);
      }
      T->Fill();
      memset(&adc, 0, sizeof(adc));
      memset(&padc, 0, sizeof(padc));
    }

    cnt++;

  } while (1);

  T->Write(0, TObject::kOverwrite);

  TCanvas *c = new TCanvas("c", "c", 750, 938);
  c->Divide(2, 4);

  c->cd(1);
  h_trigger->Draw();

  c->cd(2);
  h_ev_id->Draw();

  c->cd(3);
  h_adc1->Draw();

  c->cd(4);
  h_adc2->Draw();

  c->cd(5);
  TH1F *frame = gPad->DrawFrame(0, 0, 128, 1000);
  frame->SetTitle("Pedestals; channel; pedestals (adc)");
  frame->GetYaxis()->SetTitleOffset(1.55);
  frame->GetXaxis()->SetNdivisions(4, 5, 1, kFALSE);
  gr_ped->Draw("l");

  c->cd(6);
  frame = gPad->DrawFrame(0, -10, 128, 300);
  frame->SetTitle("Raw sigmas; channel; raw sigmas (adc)");
  frame->GetYaxis()->SetTitleOffset(1.55);
  frame->GetXaxis()->SetNdivisions(4, 5, 1, kFALSE);
  gr_raw_sigmas->Draw("l");

  c->cd(7);
  histo_raw_sigmas->Draw();

  c->Update();

  cout << "We have found " << T->GetEntries() << " events." << endl;
  out_file->Write(0, TObject::kOverwrite);
}
