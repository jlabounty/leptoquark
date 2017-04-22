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

using namespace std;

//AVERAGE
double Average(vector<double> v)
{      double sum=0;
       for(int i=0;(unsigned)i<v.size();i++)
               sum+=v[i];
       return sum/v.size();
}
//DEVIATION
double Deviation(vector<double> v, double ave)
{
       double E=0;
       for(int i=0;(unsigned)i<v.size();i++)
               E+=(v[i] - ave)*(v[i] - ave);
       return sqrt(1/static_cast<double>(v.size())*E);
}

int plot_LQGENEP_DiffEta_DiffPhi()
{

	gStyle->SetOptStat(0);

	std::string inputFile = "./outdir/TestOut.14000event.root";
//	std::string inputFile = "./outdir/TestOut.1000event.root";
	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("EICTree");

	int vsize = t->Draw("particles.eta:particles.phi:particles.I:particles.KS","","goff");
	vector<int> v_I, v_KS, v_KF, v_ORIG;
	vector<double> v_phi, v_eta, v_mom;
	vector<double> v_px, v_py, v_pz, v_E, v_m;
	vector<double> v_trans;
	vector<int> v_event;
	double I_old = 1000;
	int event_total = 0;

	for (int i = 0; i < vsize; i++)
	{

		v_eta.push_back(t->GetV1()[i]);
		v_phi.push_back(t->GetV2()[i]);
		v_I.push_back(t->GetV3()[i]);
		v_KS.push_back(t->GetV4()[i]);

		if(I_old > v_I[i]) event_total++;
		I_old = v_I[i];
		v_event.push_back(event_total);
//		if((i % 1000 == 0) && (i > 999)) cout << i << " " << event_total <<  endl;
	}
	t->Draw("particles.id:particles.orig:particles.E:particles.px","","goff");
	for (int i = 0; i < vsize; i++)
	{
		v_KF.push_back(t->GetV1()[i]);
		v_ORIG.push_back(t->GetV2()[i]);
		v_E.push_back(t->GetV3()[i]);
		v_px.push_back(t->GetV4()[i]);
	}

	t->Draw("particles.py:particles.py:particles.m:particles.pt","","goff");
	for (int i = 0; i < vsize; i++)
	{
		v_py.push_back(t->GetV1()[i]);
		v_pz.push_back(t->GetV2()[i]);
		v_m.push_back(t->GetV3()[i]);
		v_trans.push_back(t->GetV4()[i]);
	}
	t->Draw("particles.p:particles.py","","goff");
	for (int i = 0; i < vsize; i++)
	{
		v_mom.push_back(t->GetV1()[i]);
	}

	//--------------------------------------------------------------------------------------------------------------
	
	int k;	
	const int tau_code = 15;	//tau daughters
//	const int tau_code = 1;		//DIS particle daughters
	const int exec_code = 21;
	int mother_tau;
	vector<int> v_daughters_I;
	vector<int> v_daughters_line;
	bool daughter_check = false;
	int Iold = 0, Inew;

	vector<double> v_MotherTau_eta, v_MotherTau_phi;
	double mother_tau_eta_i, mother_tau_phi_i;

	cout << "Size: " << v_event.size() << endl;

	for(int i = 0; (unsigned)i < v_event.size(); i++) 			//loop over all entries in the list
	{
		for(int j = 0; j < event_total; j++)				//loop over all events
		{
			if((v_event[i] == j) && (v_KS[i] != exec_code))		//for each event and not exec_code (21)
			{
				while(k == 0)					//set k = 1 when all daughters found
				{
					if(v_KF[i] == tau_code)
					{
						mother_tau = v_I[i];
						mother_tau_eta_i = v_eta[i];
						mother_tau_phi_i = v_phi[i];
						v_daughters_I.push_back(mother_tau);
						v_daughters_line.push_back(i);
						v_MotherTau_eta.push_back(mother_tau_eta_i);
						v_MotherTau_phi.push_back(mother_tau_phi_i);
					}

					for(int l = 0; (unsigned)l < v_daughters_I.size(); l++)
					{
						if(v_ORIG[i] == v_daughters_I[l]) daughter_check = true;
						if(daughter_check == true)
						{
							v_daughters_I.push_back(v_I[i]);
							v_daughters_line.push_back(i);
							v_MotherTau_eta.push_back(mother_tau_eta_i);
							v_MotherTau_phi.push_back(mother_tau_phi_i);
						}
						daughter_check = false;
					}
					k = 1;
				}
				k = 0;
			}
		}
		Inew = v_I[i];
		if(Inew < Iold) v_daughters_I.clear();
		Iold = Inew;
	}
	
	cout << "Total Daughter Particles: " << v_daughters_line.size() << endl;

	//--------------------------------------------------------------------------------------------------------------

	vector<int> v_d_I, v_d_KS, v_d_KF, v_d_ORIG;
	vector<double> v_d_phi, v_d_eta, v_d_mom;
	vector<double> v_d_DiffPhi, v_d_DiffEta;
	vector<double> v_d_DiffTheta;
	vector<double> v_d_px, v_d_py, v_d_pz, v_d_E, v_d_m;
	vector<double> v_d_trans;
	vector<int> v_d_event;

	for(int i = 0; (unsigned)i < v_I.size(); i++)
	{
		for(int j = 0; (unsigned)j < v_daughters_line.size(); j++)
		{
			if(i == v_daughters_line[j])
			{
				v_d_I.push_back(v_I[i]);		//particle number within each event
				v_d_KS.push_back(v_KS[i]);		//bookkeeping code
				v_d_KF.push_back(v_KF[i]);		//particle data group code
				v_d_ORIG.push_back(v_ORIG[i]);		//particle which spawned this particle
				v_d_event.push_back(v_event[i]);	//event number
				v_d_mom.push_back(v_mom[i]);		//momentum of particle
				v_d_eta.push_back(v_eta[i]);		//eta of particle
				v_d_phi.push_back(v_phi[i]);		//phi of particle
				v_d_DiffEta.push_back(v_eta[i] - v_MotherTau_eta[j]);		//eta of particle - eta of tau
				v_d_DiffTheta.push_back(2*TMath::ATan(TMath::Power(TMath::E(),(-1*v_eta[i]))) - 2*TMath::ATan(TMath::Power(TMath::E(),(-1*v_MotherTau_eta[j]))));		//theta of particle - theeta of tau
				if((v_phi[i]*5 < v_MotherTau_phi[j]) && (v_MotherTau_phi[j] > 4) )
				{
					v_d_DiffPhi.push_back(2*TMath::Pi() + v_phi[i] - v_MotherTau_phi[j]);		//phi of particle - eta of tau
				}
				else if((v_phi[i]/5 > v_MotherTau_phi[j]) && (v_phi[i] > 4))
				{
					v_d_DiffPhi.push_back(v_phi[i] - v_MotherTau_phi[j] - 2*TMath::Pi());		//phi of particle - eta of tau
				}
				else
				{
					v_d_DiffPhi.push_back(v_phi[i] - v_MotherTau_phi[j]);		//phi of particle - eta of tau
				}
				v_d_px.push_back(v_px[i]);		//x momentum
				v_d_py.push_back(v_py[i]);		//y momentum
				v_d_pz.push_back(v_pz[i]);		//z momentum
				v_d_E.push_back(v_E[i]);		//particle energy
				v_d_m.push_back(v_m[i]);		//particle mass
				v_d_trans.push_back(v_trans[i]);	//transverse momentum
			}
		}
	}

	cout << "Daughter vectors created" << endl;

	//--------------------------------------------------------------------------------------------------------------
	
	TFile *f2 = new TFile("DeltaEta_DeltaPhi.root","RECREATE");

	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.4);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
/*
	TCanvas *c1 = new TCanvas();
	TH2F *h  = new TH2F("h","#Delta#eta vs. #Delta#phi",120,-1.5,1.5,120,-1.5,1.5);
	for(int i = 0; (unsigned)i < v_d_eta.size(); i++)
	{
		h->Fill(v_d_DiffEta[i],v_d_DiffPhi[i]);
	}
	h->SetTitle("#Delta#eta vs. #Delta#phi for daughters of #tau from  e-p Leptoquark Events");
	h->GetXaxis()->SetTitle("#Delta#eta");
	h->GetYaxis()->SetTitle("#Delta#phi");
	h->GetZaxis()->SetTitle("Number of Particles");
	h->Draw("COLZ");
	c1->SetRightMargin(0.15);
//	c1->Print("DeltaEta_vs_DeltaPhi.eps");
//	c1->Print("DeltaEta_vs_DeltaPhi.root");
*/	
	//--------------------------------------------------------------------------------------------------------------


	TCanvas *c2 = new TCanvas();

	TTree *t2 = new TTree("tvec","Tree with vectors of tau daughter particles");
		t2->Branch("DIFFeta",&v_d_DiffEta);
		t2->Branch("DIFFphi",&v_d_DiffPhi);
		t2->Branch("DIFFtheta",&v_d_DiffTheta);
		t2->Branch("energy",&v_d_E);
		t2->Branch("I",&v_d_I);
		t2->Branch("KS",&v_d_KS);
		t2->Branch("KF",&v_d_KF);
		t2->Branch("ORIG",&v_d_ORIG);
		t2->Branch("event",&v_d_event);
		t2->Branch("momentum",&v_d_mom);
		t2->Branch("eta",&v_d_eta);
		t2->Branch("phi",&v_d_phi);
		t2->Branch("px",&v_d_px);
		t2->Branch("py",&v_d_py);
		t2->Branch("pz",&v_d_pz);
		t2->Branch("mass",&v_d_m);
		t2->Branch("pTrans",&v_d_trans);
	t2->Fill();	

	TH2F *h2  = new TH2F("h2","#Delta#eta vs. #Delta#phi",40,-1,1,40,-1,1);
		h2->SetTitle("#Delta#eta vs. #Delta#phi for daughters of #tau from  e-p Leptoquark Events");
		h2->GetXaxis()->SetTitle("#Delta#eta");
		h2->GetYaxis()->SetTitle("#Delta#phi");
		h2->GetZaxis()->SetTitle("Particles/Event");
//		h2->SetMaximum(14);
//	t2->Draw("DIFFphi:DIFFeta>>h2","(KF!=1)*(1/14000)","colz same");
	t2->Draw("DIFFphi:DIFFeta>>h2","(KF!=15)*(1/14000)","colz same");
//	t2->Draw("DIFFphi:DIFFeta>>h2","(KF!=15)*(energy/14000)","colz same");
//	t2->Draw("DIFFphi:DIFFeta>>h2","(KF!=1)*(energy/14000)","colz same");
	h2->Draw("colz");

	c2->Update();

	TCanvas *c3 = new TCanvas();
	TH1D *h3 = h2->ProjectionX();
	h3->Draw();

	TCanvas *c4 = new TCanvas();
	TH1D *h4 = h2->ProjectionY();
	h4->Draw();

	TCanvas *c5 = new TCanvas();
	TH2F *h5  = new TH2F("h5","#Delta#theta vs. #Delta#phi",60,-0.5,0.5,60,-0.5,0.5);
		h5->SetTitle("#Delta#theta vs. #Delta#phi for daughters of #tau from  e-p Leptoquark Events");
		h5->GetXaxis()->SetTitle("#Delta#theta");
		h5->GetYaxis()->SetTitle("#Delta#phi");
		h5->GetZaxis()->SetTitle("Average Energy/Event");
	t2->Draw("DIFFphi:DIFFtheta>>h5","(KF!=15)*(energy/14000)","colz same");
	h5->Draw("colz");

	TCanvas *c6 = new TCanvas();
	TH1D *h6 = h5->ProjectionX();
	h6->Draw();


	f2->Write();
	
	//--------------------------------------------------------------------------------------------------------------

	return 0;
}
