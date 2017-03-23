#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"

using namespace std;

int plot_tauOnly_deltaeta_v_deltaphi(std::string inputFile);

int plot_tauOnly_deltaeta_v_deltaphi()
{
	std::string inputFile = "./simdir/tau/g4cemc_eval.root";
	
	plot_tauOnly_deltaeta_v_deltaphi(inputFile);

	return 0;
}


int plot_tauOnly_deltaeta_v_deltaphi(std::string inputFile)
{
	TFile *f = TFile::Open(inputFile.c_str());
	TNtuple *ntp_tower = (TNtuple*)f->Get("ntp_tower");

	//Declaration of leaves types
	Float_t         event;
	Float_t         towerID;
	Float_t         ieta;
	Float_t         iphi;
	Float_t         eta;
	Float_t         phi;
	Float_t         e;
	Float_t         gparticleID;
	Float_t         gflavor;
	Float_t         gnhits;
	Float_t         geta;
	Float_t         gphi;
	Float_t         ge;
	Float_t         gpt;
	Float_t         gvx;
	Float_t         gvy;
	Float_t         gvz;
	Float_t         gembed;
	Float_t         gedep;
	Float_t         efromtruth;

	// Set branch addresses.
	ntp_tower->SetBranchAddress("event",&event);
	ntp_tower->SetBranchAddress("towerID",&towerID);
	ntp_tower->SetBranchAddress("ieta",&ieta);
	ntp_tower->SetBranchAddress("iphi",&iphi);
	ntp_tower->SetBranchAddress("eta",&eta);
	ntp_tower->SetBranchAddress("phi",&phi);
	ntp_tower->SetBranchAddress("e",&e);
	ntp_tower->SetBranchAddress("gparticleID",&gparticleID);
	ntp_tower->SetBranchAddress("gflavor",&gflavor);
	ntp_tower->SetBranchAddress("gnhits",&gnhits);
	ntp_tower->SetBranchAddress("geta",&geta);
	ntp_tower->SetBranchAddress("gphi",&gphi);
	ntp_tower->SetBranchAddress("ge",&ge);
	ntp_tower->SetBranchAddress("gpt",&gpt);
	ntp_tower->SetBranchAddress("gvx",&gvx);
	ntp_tower->SetBranchAddress("gvy",&gvy);
	ntp_tower->SetBranchAddress("gvz",&gvz);
	ntp_tower->SetBranchAddress("gembed",&gembed);
	ntp_tower->SetBranchAddress("gedep",&gedep);
	ntp_tower->SetBranchAddress("efromtruth",&efromtruth);


	TCanvas *c1 = new TCanvas();
	ntp_tower->Draw("eta:phi:e","event==1","colz");
	c1->Update();

	vector<int> vec_event;
	vector<double> vec_eta,vec_phi,vec_e;

	for (int i = 0, N = ntp_tower->GetEntries(); i < N; ++i) 
	{
		ntp_tower->GetEntry(i);
		vec_event.push_back(event);
		vec_eta.push_back(eta);
		vec_phi.push_back(phi);
		vec_e.push_back(e);
	}

	double e_max_prev = 0;
	int i_max = 0;

	int max_event = vec_event[vec_event.size()-1];
	vector<double> delta_eta;
	vector<double> delta_phi;
	double eta_corr = 0;
	double phi_corr = 0;

	for(int j = 0; j < max_event; j++)
	{
		for(int i = 0; i < vec_event.size(); i++)
		{
			if(vec_event[i] == j)
			{
				if(vec_e[i] > e_max_prev)
				{
					i_max = i;
					e_max_prev = vec_e[i];
				}
			}
		}
		for(int i = 0; i < vec_event.size(); i++)
		{
			eta_corr = vec_eta[i_max];
			phi_corr = vec_phi[i_max];

			if(vec_event[i] == j)
			{
				delta_eta.push_back(vec_eta[i] - eta_corr);
				delta_phi.push_back(vec_phi[i] - phi_corr);
					
			}
		}
	}

	TCanvas *c2 = new TCanvas();
	TGraph2D *gr = new TGraph2D(delta_eta.size(), &(delta_eta[0]), &(delta_phi[0]), &(vec_e[0]));
	gr->Draw();
	c2->Update();

	TCanvas *c3 = new TCanvas();
	gr->Draw("colz");


	return 0;
}
