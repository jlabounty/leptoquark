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
#include <algorithm>

#include <fstream>
#include <string>
#include <math.h>
#include "TGraph.h"


int Lepto()
{

        int I, KS, KF, ORIG;
        double px, py, pz, E, m;
	double phi, eta, mom;
	double trans;
	vector<double> v_phi, v_eta, v_mom;
	vector<double> v_px, v_py, v_pz, v_E, v_m;
	vector<double> v_trans;
	vector<int> v_event;

	vector<int> v_d_event;
	vector<double> v_d_phi, v_d_eta, v_d_mom;
	vector<double> v_d_px, v_d_py, v_d_pz, v_d_E, v_d_m;
	vector<double> v_d_trans;

	vector<int> v_21_event;
	vector<double> v_21_phi, v_21_eta, v_21_mom;
	vector<double> v_21_px, v_21_py, v_21_pz, v_21_E, v_21_m;
	vector<double> v_21_trans;
	vector<double> v_21_orig;
	vector<int> v_21_I;

	std::ifstream file("TestOut_1000Events.txt");
//	std::ifstream file("TestOut_10_250_LQ13.txt"); //Full number of events!
//	std::ifstream file("TestOut_10_250.txt"); //Full number of events!
//	std::ifstream file("TestOut_20_250.txt"); //Full number of events!
//	std::ifstream file("TestOut_20_325.txt"); //Full number of events!
        std::string str;

	int tau_orig = 0;
	int non_daughter = 0;
	int daughter = 0;
	int mother = 0;
	int high_eta = 0;
	int I_old = 1000, I_new, event_total;
	int neg_eta = 0;
	int problem_eta = 0;

	int line = 0;

        while (std::getline(file, str))
        {
		if(file >> I >> KS >> KF >> ORIG >> px >> py >> pz >> E >> m)
		{
			line++;
			I_new = I;
			if(I_new < I_old) event_total++;
			I_old = I_new;

			if((KS != 21) && (KF == 15) && (ORIG != 8)) cout << "WTF???  Event: " << event_total << " --- Line: " <<  line + 1 << endl;

			if((KS != 21) && (KF == 15) && (ORIG == 8))
			{
				//Calculate the total momentum, transverse momentum, and eta
				mom = sqrt(px**2 + py**2 + pz**2);
				trans = sqrt(px**2 + py**2);
				eta = (0.5)*log((mom+trans)/(mom-trans));
				if(pz < 0)
				{
					eta = -1.0*eta;
					neg_eta++;
				}
				phi = atan(py/px);
				if(px < 0 && py > 0) phi = phi + TMath::Pi();
				if(px < 0 && py < 0) phi = phi - TMath::Pi();

				//Push these values back into vectors
				v_event.push_back(event_total);
				v_mom.push_back(mom);
				v_eta.push_back(eta);
				v_phi.push_back(phi);
				v_px.push_back(px);
				v_py.push_back(py);
				v_pz.push_back(pz);
				v_E.push_back(E);
				v_m.push_back(m);
				v_trans.push_back(trans);

				//Increment counters
				if(eta > 1.1 && eta < 4.0) high_eta++;
				if((eta > 4.0) || (eta < 1.45 && eta > 1.10))
				{
					problem_eta++;
					//cout << eta << endl;
				}
				mother++;
				tau_orig = I;
			}
			else if((KS != 21))
			{
				v_d_event.push_back(event_total);

				mom = sqrt(px**2 + py**2 + pz**2);
				if((KF == 12) || (KF == 14) || (KF == 16) || (KF == 18))
				{
					trans = 0.000000;
					daughter++;
				}
				else
				{
					trans = sqrt(px**2 + py**2);
				}
				eta = (0.5)*log((mom+trans)/(mom-trans));
				if(pz < 0)
				{
					eta = -1.0*eta;
				}
				phi = atan(py/px);
				if(px < 0 && py > 0) phi = phi + TMath::Pi();
				if(px < 0 && py < 0) phi = phi - TMath::Pi();

				v_d_mom.push_back(mom);
				v_d_eta.push_back(eta);
				v_d_phi.push_back(phi);
				v_d_px.push_back(px);
				v_d_py.push_back(py);
				v_d_pz.push_back(pz);
				v_d_E.push_back(E);
				v_d_m.push_back(m);
				v_d_trans.push_back(trans);

			}
			else if(KS != 21)
			{

			}
			else if(KS == 21)
			{

				v_21_event.push_back(event_total);

				mom = sqrt(px**2 + py**2 + pz**2);
				trans = sqrt(px**2 + py**2);
				eta = (0.5)*log((mom+trans)/(mom-trans));
				if(pz < 0)
				{
					eta = -1.0*eta;
				}
				phi = atan(py/px);
				if(px < 0 && py > 0) phi = phi + TMath::Pi();
				if(px < 0 && py < 0) phi = phi - TMath::Pi();

				v_21_mom.push_back(mom);
				v_21_eta.push_back(eta);
				v_21_phi.push_back(phi);
				v_21_px.push_back(px);
				v_21_py.push_back(py);
				v_21_pz.push_back(pz);
				v_21_E.push_back(E);
				v_21_m.push_back(m);
				v_21_trans.push_back(trans);
				v_21_orig.push_back(ORIG);
				v_21_I.push_back(I);
			}

		}
        }

	double ptmiss, pxSUM, pySUM, pxSUM_orig, pySUM_orig;
	double event_count = 0;
	double last_count = 1;
	vector<double> pt_orig, pt_miss;
	for(int i=0; i < v_21_event.size(); i++)
	{
		if(v_21_event[i] > event_count) event_count++;

		if(event_count == last_count)
		{
			if( (v_21_I[i] ==  1) || (v_21_I[i] == 2))
			{
				pxSUM_orig = pxSUM_orig + v_21_pz[i];
				pySUM_orig = 0.000;
			}
			else
			{
				pxSUM = pxSUM + v_21_px[i];
				pySUM = pySUM + v_21_py[i];
			}
		}

		if(event_count > last_count)
		{
			pt_orig.push_back(TMath::Sqrt((pxSUM_orig)**2 + (pySUM_orig)**2));
			pt_miss.push_back(TMath::Sqrt((pxSUM)**2 + (pySUM)**2));
			pxSUM = pySUM = pxSUM_orig = pySUM_orig = 0;
			last_count++;
		}

	}
/*
	for (int i = 0; i < pt_orig.size(); i++)
	{
		cout << pt_orig[i] << "   " << pt_miss[i] << endl;
	}
*/
	cout << endl << "=================================================" << endl << endl;
	cout << "Events: " << event_total << endl;
	cout << "Mother (tau) Particles: " << mother << endl;
	cout << "    Negative Eta (Electron Going): " << neg_eta << endl;
	cout << "    tau in FEMCAL: " << high_eta << endl;
	cout << "    tau in Barrel: " << mother - high_eta << endl;
//	cout << "Daughter Particles: " << daughter << endl;
//	cout << "Non-Daughter Particles: " << non_daughter << endl;
	cout << endl << "Problem eta: " << problem_eta << endl;
	cout << "Problem eta (%): " << ((problem_eta)/(mother))*100 << endl;
	cout << endl << "=================================================" << endl << endl;

	TCanvas *c1 = new TCanvas();
	TH2F *h	 = new TH2F("h","#eta vs. #phi",160,-8,8,80,-4,4);
		for(int i = 0; i < v_eta.size(); i++)
		{
			h->Fill(v_eta[i],v_phi[i]);
		}
		h->SetTitle("#eta vs. #phi of Tau Distribution from e-p Leptoquark Events");
		h->GetXaxis()->SetTitle("#eta");
		h->GetYaxis()->SetTitle("#phi");
	h->Draw("COLZ");

	TCanvas *c3 = new TCanvas();
	TH2F *h2 = new TH2F("h2","p_{x} vs. p_{y}",60,-60,60,60,-60,60);
		for(int i = 0; i < v_px.size(); i++)
		{
			h2->Fill(v_px[i],v_py[i]);
		}
		h2->SetTitle("p_{x} vs. p_{y} of Tau Distribution from e-p Leptoquark Events");
		h2->GetXaxis()->SetTitle("p_{x}");
		h2->GetYaxis()->SetTitle("p_{y}");
	h2->Draw("COLZ");

	TCanvas *c4 = new TCanvas();
	TH2F *h4 = new TH2F("h3","Trans vs. Momentum",100,0,100,100,0,100);
		for(int i = 0; i < v_mom.size(); i++)
		{
			h4->Fill(v_trans[i],v_mom[i]);
		}
		h4->SetTitle("p_{T} vs. Total p of #tau from e-p Leptoquark Events (#tau only)");
		h4->GetXaxis()->SetTitle("p_{T}");
		h4->GetYaxis()->SetTitle("#Sigma p");
		double ding = h4->GetYaxis()->GetXmax();
		double dong = h4->GetXaxis()->GetXmax();
	h4->Draw("COL");

	TH1D * projh4X = h4->ProjectionX();
	TH1D * projh4Y = h4->ProjectionY();

	gStyle->SetPalette(1);

	Float_t rightmax = 1.1*projh4X->GetMaximum();
	Float_t scale = ding/rightmax;
	projh4X->Scale(scale);
	projh4X->SetLineColor(kRed);
	projh4X->SetFillColor(kRed);
	projh4X->SetFillStyle(3001);
	projh4X->Draw("bar SAME");

	Float_t topmax = 1.1*projh4Y->GetMaximum();
	Float_t scale_2 = dong/topmax;
	projh4Y->Scale(scale_2);
	projh4Y->SetFillColor(kBlue);
	projh4Y->SetFillStyle(3001);
	projh4Y->Draw("hbar SAME");

	TCanvas *c5 = new TCanvas();
	h4->Draw("COLZ");

	TCanvas *c55 = new TCanvas();
	projh4X->Draw("bar");

	TCanvas *c6 = new TCanvas();
	TH2F *h6 = new TH2F("h6","p_{T}^{Miss} vs. p_{T}^{original}",500,0,500,500,0,500);
		for(int i = 0; i < pt_orig.size(); i++)
		{
			h6->Fill(pt_orig[i],pt_miss[i]);
		}
		h6->SetTitle("p_{T}^{Miss} vs. p_{T}^{original} Distribution from e-p Leptoquark Events");
		h6->GetXaxis()->SetTitle("p_{T}^{original}");
		h6->GetYaxis()->SetTitle("p_{T}^{Miss}");
	h6->Draw("COLZ");

	TCanvas* c7 = new TCanvas();
	TH1F *h7 = new TH1F("h7", "p_{T}^{Miss}", 80,-0.25,39.75);
		for(int i = 0; i < pt_miss.size(); i++)
		{
			h7->Fill(pt_miss[i]);
		}
		h7->SetTitle("p_{T}^{Miss} Distribution from e-p Leptoquark Events");
		h7->GetYaxis()->SetTitle("Events");
		h7->GetXaxis()->SetTitle("p_{T}^{Miss}");
		h7->SetFillColor(kBlue);
		h7->SetFillStyle(3001);
	h7->Draw();

	TCanvas* c8 = new TCanvas();
	TH1F *h8 = new TH1F("h8", "p_{z-#tau}", 520,-10.5,249.5);
		for(int i = 0; i < v_pz.size(); i++)
		{
			h8->Fill(v_pz[i]);
		}
		h8->SetTitle("p_z Distribution of #tau produced from e-p Leptoquark Events");
		h8->GetYaxis()->SetTitle("Events");
		h8->GetXaxis()->SetTitle("p_z");
		h8->SetFillColor(kBlue);
		h8->SetFillStyle(3001);
	h8->Draw();

        TCanvas *c9 = new TCanvas();
        TH2F *h9 = new TH2F("h9","p_{#tau} vs. #eta_{#tau}",200,-5,5,200,0,250);
                for(int i = 0; i < v_eta.size(); i++)
                {
                        h9->Fill(v_eta[i], TMath::Abs(v_mom[i]));
                }
                h9->SetTitle("p_{#tau} vs. #eta_{#tau} from e-p Leptoquark Events");
                h9->GetYaxis()->SetTitle("p_{#tau}");
                h9->GetXaxis()->SetTitle("#eta");
        h9->Draw("COLZ");


/*
	TCanvas *c5 = new TCanvas();
	TH1F *h5 = new TH1F("h5","Transverse Momentum",200,0,100);
		for(int i = 0; i < v_trans.size(); i++)
		{
			h5->Fill(v_trans[i]);
		}
		h5->SetTitle("p_{T} of Tau Distribution from e-p Leptoquark Events");
		h5->GetXaxis()->SetTitle("p_{T}");
	h5->Draw();
	h5->SetLineColor(kRed);
	h5->SetFillColor(kRed);
	h5->SetFillStyle(3001);
	TH1F *h6 = new TH1F("h6","Total Momentum",200,0,100);
		for(int i = 0; i < v_mom.size(); i++)
		{
			h6->Fill(v_mom[i]);
		}
		h6->GetXaxis()->SetTitle("#Sigma p");
	h6->Draw("SAME");
	h6->SetLineColor(kBlue);
	h6->SetFillStyle(3001);
	h6->SetFillColor(kBlue);
	c5->Update();
*/
	cout << "\a";
        return 0;
}
