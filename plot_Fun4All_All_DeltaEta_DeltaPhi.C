#include "TROOT.h"
#include "TClass.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "Riostream.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TMath.h"
#include "math.h"
#include "TColor.h"
#include <vector>
#include <sstream>
#include <algorithm>

#include <cstdlib>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>
#include "TGraph.h"
#include "TGraph2D.h"
#include <algorithm>

int plot_Fun4All_All_DeltaEta_DeltaPhi()
{
	const std::string inFile = "LeptoAna_1000events_All.root";
	const std::string inDirectory = "/gpfs/mnt/gpfs02/phenix/scratch/jlab/Leptoquark/";
	std::string inputFile = inDirectory+inFile;

	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("ntp_leptoquark");

//	const int Nevent = t->GetMaximum("event");
	const int Nevent = 100;
	cout << "Running " << Nevent << " events" << endl;
	const int Nentries1 = t->Draw("isMaxEnergyJet","(isMaxEnergyJet<10)*(calorimeterid<10)","goff");
		Double_t *arr_jet = t->GetV1();
		vector<int> v_jet(Nentries1);
		for(int i = 0; i < Nentries1; i++)
		{
			v_jet[i] = (int)arr_jet[i];
		}
	const int Nentries = t->Draw("towereta:towerphi:towerenergy:event","(isMaxEnergyJet<10)*(calorimeterid<10)","goff");
		Double_t *arr_eta = t->GetV1();
		Double_t *arr_phi = t->GetV2();
		Double_t *arr_e = t->GetV3();
		Double_t *arr_event = t->GetV4();


	if(Nentries1 != Nentries)
	{
		cerr << "ERROR: Draw commands to not return the same dimensions. Check that your logical expressions are the same." << endl;
		return -1;
	}

	vector<double> v_DeltaEta, v_DeltaPhi, v_DeltaTheta, v_Energy;
	vector<double> v_DeltaEta_j2, v_DeltaPhi_j2, v_DeltaTheta_j2, v_Energy_j2;

	for(int i = 0; i < Nevent; i++)
	{
		double Emax = 0;
		int Emax_i = 0;

		double Emax_j2 = 0;
		int Emax_i_j2 = 0;

		for(int j = 0; j < Nentries; j++)
		{
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 1))
			{
				if(t->GetV3()[j] > Emax) 
				{
					Emax = t->GetV3()[j];
					Emax_i = j;
				}
			}
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 2))
			{
				if(t->GetV3()[j] > Emax_j2) 
				{
					Emax_j2 = t->GetV3()[j];
					Emax_i_j2 = j;
				}
			}
		}

		for(int j = 0; j < Nentries; j++)
		{
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 1))
			{
				v_DeltaEta.push_back(t->GetV1()[j] - t->GetV1()[Emax_i]);
				v_DeltaTheta.push_back(2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[j])) - 2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[Emax_i])));
				v_DeltaPhi.push_back(t->GetV2()[j] - t->GetV2()[Emax_i]);
				v_Energy.push_back(t->GetV3()[j]);
			}
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 2))
			{
				v_DeltaEta_j2.push_back(t->GetV1()[j] - t->GetV1()[Emax_i_j2]);
				v_DeltaTheta_j2.push_back(2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[j])) - 2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[Emax_i_j2])));
				v_DeltaPhi_j2.push_back(t->GetV2()[j] - t->GetV2()[Emax_i_j2]);
				v_Energy_j2.push_back(t->GetV3()[j]);
			}
		}
	}


//-----------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	double xmin = -1;
	double xmax = 1;
	double ymin = -0.5;
	double ymax = 0.5;

	std::string title = "isMaxJetEnergy = 1";
	TCanvas *c1 = new TCanvas();
	TH2D *h1 = new TH2D("h1",title.c_str(),40,xmin,xmax,40,ymin,ymax);
	for(int i = 0; (unsigned)i < v_DeltaEta.size(); i++)
	{
		h1->Fill(v_DeltaEta[i],v_DeltaPhi[i], v_Energy[i]/Nevent);
	}
	c1->SetLogz();
	h1->Draw("colz");
//		h1->SetMinimum(1);
//		h1->SetMaximum(1000);
		h1->GetXaxis()->SetTitle("#Delta#eta");
		h1->GetYaxis()->SetTitle("#Delta#phi");
	c1->Update();

//-----------------------------------------------------------------------------------------------------------


	std::string title2 = "isMaxJetEnergy = 2";
	TCanvas *c2 = new TCanvas();
	TH2D *h2 = new TH2D("h2",title2.c_str(),40,xmin,xmax,40,ymin,ymax);
	for(int i = 0; (unsigned)i < v_DeltaEta_j2.size(); i++)
	{
		h2->Fill(v_DeltaEta_j2[i],v_DeltaPhi_j2[i], v_Energy_j2[i]/Nevent);
	}
	c2->SetLogz();
	h2->Draw("colz");
//		h2->SetMinimum(1);
//		h2->SetMaximum(1000);
		h2->GetXaxis()->SetTitle("#Delta#eta");
		h2->GetYaxis()->SetTitle("#Delta#phi");
	c2->Update();


//-----------------------------------------------------------------------------------------------------------


	TCanvas *c3 = new TCanvas();
	TH1D *h3 = h1->ProjectionY();
		h3->SetLineColor(kGreen+3);
	TH1D *h4 = h2->ProjectionY();
		h4->SetLineColor(kBlue);

	TF1 *f3 = new TF1("f3", "gaus", ymin, ymax);
		f3->SetLineColor(kGreen+3);
	h3->Fit("f3","Q R");

	TF1 *f4 = new TF1("f4", "gaus", ymin, ymax);
		f4->SetLineColor(kBlue+2);
	h4->Fit("f4","Q R");

	h3->Draw();
	h4->Draw("SAME");
	TLegend *leg = new TLegend(0.2,0.9,0.7,0.75);
		leg->SetBorderSize(1);
		leg->AddEntry(h3,title.c_str(),"l");
		leg->AddEntry(h4,title2.c_str(),"l");
	leg->Draw("SAME");
	c3->Update();

	cout << endl;
	cout << "*******************************************************" << endl;

	cout << "Ratio of Gaussian Width to Max Energy:" << endl;
	cout << "     for: " << title << " -> " << f3->GetParameter(2) << " / " << h3->GetBinContent(h3->GetMaximumBin()) 
		<< " = " << f3->GetParameter(2) /  h3->GetBinContent(h3->GetMaximumBin()) << endl;
	cout << "     for: " << title2 << " -> " << f4->GetParameter(2) << " / " << h4->GetBinContent(h4->GetMaximumBin()) 
		<< " = " << f4->GetParameter(2) /  h4->GetBinContent(h4->GetMaximumBin()) << endl;
	cout << "     Ratio of ratios: " << (f3->GetParameter(2) /  h3->GetBinContent(h3->GetMaximumBin())) / (f4->GetParameter(2) /  h4->GetBinContent(h4->GetMaximumBin())) << endl;

	cout << "*******************************************************" << endl;
	cout << endl;

	c1->Close();
	c2->Close();
	c3->Close();

//-----------------------------------------------------------------------------------------------------------

	vector<double> v_Ratios;

	for(int i = 0; i < Nevent; i++)
	{
		TH1D *h5 = new TH1D("h5","",20,ymin*2,ymax*2);
		TH1D *h6 = new TH1D("h6","",20,ymin*2,ymax*2);

		for(int j = 0; j < Nentries/100; j++)
		{
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 1))
			{
				h5->Fill(v_DeltaPhi[j], v_Energy[j]);
			}
			if((t->GetV4()[j]-1 == i) && (v_jet[j] == 2))
			{
				h6->Fill(v_DeltaPhi_j2[j], v_Energy_j2[j]);
			}
		}

		TF1 *f5 = new TF1("f5", "gaus", ymin, ymax);
			f3->SetLineColor(kGreen+3);
		h5->Fit("f5","Q");

		TF1 *f6 = new TF1("f6", "gaus", ymin, ymax);
			f6->SetLineColor(kBlue+2);
		h6->Fit("f6","Q");

		TCanvas *c5 = new TCanvas();
		h5->Draw();

		TCanvas *c6 = new TCanvas();
		h6->Draw();

//		cout << endl;
//		cout << "Ratio of Gaussian Width to Max Energy:" << endl;
//		cout << "     for: " << title << " -> " << f5->GetParameter(2) << " / " << h5->GetBinContent(h5->GetMaximumBin()) 
//			<< " = " << f5->GetParameter(2) /  h5->GetBinContent(h5->GetMaximumBin()) << endl;
//		cout << "     for: " << title2 << " -> " << f6->GetParameter(2) << " / " << h6->GetBinContent(h6->GetMaximumBin()) 
//			<< " = " << f6->GetParameter(2) /  h6->GetBinContent(h6->GetMaximumBin()) << endl;
		cout << "     Ratio of ratios for event " << i+1 << " : " << (f5->GetParameter(2) /  h5->GetBinContent(h5->GetMaximumBin())) /
			 (f6->GetParameter(2) /  h6->GetBinContent(h6->GetMaximumBin())) << endl;

		if((f5->GetParameter(2) /  h5->GetBinContent(h5->GetMaximumBin())) /
                         (f6->GetParameter(2) /  h6->GetBinContent(h6->GetMaximumBin())) != NAN)
		v_Ratios.push_back((f5->GetParameter(2) /  h5->GetBinContent(h5->GetMaximumBin())) /
                         (f6->GetParameter(2) /  h6->GetBinContent(h6->GetMaximumBin())));

		delete h5;
		delete h6;
		c5->Close();
		c6->Close();
	}

	TCanvas *c7 = new TCanvas();
	TH1D *h7 = new TH1D("h7","",400,-200,200);
	for(int i = 0; (unsigned)i < v_Ratios.size(); i++)
	{
		h7->Fill(v_Ratios[i]);
	}
	h7->Draw();









	return 0;
}
