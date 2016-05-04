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

int Lepto_test()
{

	int I, KS, KF, ORIG;
	double px, py, pz, E, m;
	double phi, eta, mom;
	double trans;
	vector<int> v_I, v_KS, v_KF, v_ORIG;
	vector<double> v_phi, v_eta, v_mom;
	vector<double> v_px, v_py, v_pz, v_E, v_m;
	vector<double> v_trans;
	vector<int> v_event;

//	std::ifstream file("../TestOut_1000Events.txt");
//	std::ifstream file("TestOut_10_250_LQ13.txt"); //Full number of events!
//	std::ifstream file("../TestOut_10_250.txt"); //Full number of events!
	std::ifstream file("../TestOut_20_250.txt"); //Full number of events!
//	std::ifstream file("TestOut_20_325.txt"); //Full number of events!
	std::string str;

//        int tau_orig = 0;
//        int non_daughter = 0;
//        int daughter = 0;
//        int mother = 0;
//        int high_eta = 0;
        int I_old = 1000, I_new, event_total;
        int neg_eta = 0;
//        int problem_eta = 0;

        int line = 0;


	while (std::getline(file, str))
	{
		if(file >> I >> KS >> KF >> ORIG >> px >> py >> pz >> E >> m)
		{
			line++;
			I_new = I;
			if(I_new < I_old) event_total++;
			I_old = I_new;

			mom = sqrt(pow(px,2) + pow(px,2) + pow(pz,2));
			trans = sqrt(pow(px,2) + pow(py,2));
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
			v_I.push_back(I);			//particle number within each event
			v_KS.push_back(KS);			//bookkeeping code
			v_KF.push_back(KF);			//particle data group code
			v_ORIG.push_back(ORIG);			//particle which spawned this particle
			v_event.push_back(event_total);		//event number
			v_mom.push_back(mom);			//momentum of particle
			v_eta.push_back(eta);			//eta of particle
			v_phi.push_back(phi);			//phi of particle
			v_px.push_back(px);			//x momentum
			v_py.push_back(py);			//y momentum
			v_pz.push_back(pz);			//z momentum
			v_E.push_back(E);			//particle energy
			v_m.push_back(m);			//particle mass
			v_trans.push_back(trans);
		}
	}

	int k;	
	const int tau_code = 15;
	const int exec_code = 21;
	int mother_tau;
	vector<int> v_daughters_I;
	vector<int> v_daughters_line;
	bool daughter_check = false;
	int Iold = 0, Inew;

	cout << "Size: " << v_event.size() << endl;

	for(int i = 0; (unsigned)i < v_event.size(); i++) 				//loop over all entries in the list
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
						v_daughters_I.push_back(mother_tau);
						v_daughters_line.push_back(i);
					}
						for(int l = 0; (unsigned)l < v_daughters_I.size(); l++)
						{
							if(v_ORIG[i] == v_daughters_I[l]) daughter_check = true;
						}
						if(daughter_check == true)
						{
							v_daughters_I.push_back(v_I[i]);
							v_daughters_line.push_back(i);
						}
						daughter_check = false;
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


	vector<int> v_d_I, v_d_KS, v_d_KF, v_d_ORIG;
	vector<double> v_d_phi, v_d_eta, v_d_mom;
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

	//compute delta_eta, delta_phi for each shower
	vector<int> v_event_analysis;
	vector<double> v_delta_eta, v_delta_phi;
	int number_of_events = *std::max_element(v_d_event.begin(), v_d_event.end());
	Iold = Inew = 0;	
	int N = 0;
	double avg_eta, avg_phi, stdev_eta, stdev_phi;
	vector<double> v_oneevent_eta, v_oneevent_phi;

	for(int i = 0; (unsigned)i < v_d_I.size(); i++)
	{
		for (int j = 0; j < number_of_events - 1; j++)
		{
			if(v_d_event[i] == j)
			{
				N++;
				v_oneevent_eta.push_back(v_d_eta[i]);
				v_oneevent_phi.push_back(v_d_phi[i]);
			}
		}
		Inew = v_d_I[i];
		if(Inew < Iold)
		{
			avg_eta = Average(v_oneevent_eta);
			avg_phi = Average(v_oneevent_phi);
			stdev_eta = Deviation(v_oneevent_eta, avg_eta);
			stdev_phi = Deviation(v_oneevent_phi, avg_phi);

			v_delta_eta.push_back(stdev_eta);
			v_delta_phi.push_back(stdev_phi);

			v_oneevent_eta.clear();
			v_oneevent_phi.clear();
		}
		Iold = Inew;
	}

	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.4);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	TCanvas *c1 = new TCanvas();
	TH2F *h  = new TH2F("h","#delta #eta vs. #delta #phi",60,0,3,80,0,4);
	for(int i = 0; (unsigned)i < v_delta_eta.size(); i++)
	{
		h->Fill(v_delta_eta[i],v_delta_phi[i]);
	}
	h->SetTitle("#eta vs. #phi of Tau Distribution from e-p Leptoquark Events");
	h->GetXaxis()->SetTitle("#delta #eta");
	h->GetYaxis()->SetTitle("#delta #phi");
	h->Draw("COLZ");
	c1->SetRightMargin(0.13);
	c1->Print("DeltaEta_vs_DeltaPhi.eps");


	//compute difference in eta, difference in phi for each shower
	vector<double> v_diff_eta, v_diff_phi;
	v_oneevent_eta.clear();
	v_oneevent_phi.clear();
	Iold = Inew = 0;	
	N = 0;
	double tau_eta, tau_phi/*, diff_eta, diff_phi*/;
	for(int i = 0; (unsigned)i < v_d_I.size(); i++)
	{
		Inew = v_d_I[i];
		if(Inew < Iold)
		{
			tau_eta = v_oneevent_eta[0];
			tau_phi = v_oneevent_phi[0];
			
			for(int counter = 1; (unsigned)counter < v_oneevent_eta.size(); counter++)
			{
				v_diff_eta.push_back(v_oneevent_eta[counter] -tau_eta);
				v_diff_phi.push_back(v_oneevent_phi[counter] -tau_phi);
			}

			v_oneevent_eta.clear();
			v_oneevent_phi.clear();
		}
		Iold = Inew;
		for (int j = 0; j < number_of_events - 1; j++)
		{
			if(v_d_event[i] == j)
			{
				N++;
				v_oneevent_eta.push_back(v_d_eta[i]);
				v_oneevent_phi.push_back(v_d_phi[i]);
			}
		}
	}

	TCanvas *c2 = new TCanvas();
	TH2F *h2  = new TH2F("h2","#Delta #eta vs. #Delta #phi",160,-1,1,160,-1,1);
	for(int i = 0; (unsigned)i < v_diff_eta.size(); i++)
	{
		h2->Fill(v_diff_eta[i],v_diff_phi[i]);
	}
	h2->SetTitle("#eta vs. #phi of Tau Distribution from e-p Leptoquark Events");
	h2->GetXaxis()->SetTitle("#delta #eta");
	h2->GetYaxis()->SetTitle("#delta #phi");
	h2->Draw("COLZ");
	c2->SetRightMargin(0.13);
	c2->Print("DiffEta_vs_DiffPhi.eps");


	return 0;
}
