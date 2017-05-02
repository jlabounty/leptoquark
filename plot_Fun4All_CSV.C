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

//AVERAGE
double Average(vector<double> v)
{
	double sum=0;
	for(int i=0;(unsigned)i<v.size();i++)
	sum+=v[i];
	return sum/v.size();
}


//WEIGHTED AVERAGE
double Average(vector<double> v, vector<double> w)
{
	if(v.size() != w.size())
	{
		cerr << "ERROR in Average(v1, v2): Two vectors have different lengths" << endl;
		return nan("");
	}
	double sum=0, sum_w = 0;
	for(int i=0;(unsigned)i<v.size();i++)
	{
		sum+=(v[i]*w[i]);
		sum_w += w[i];
	}
	return sum/sum_w;
}

//DEVIATION
double Deviation(vector<double> v, double ave)
{
	double E=0;
	for(int i=0;(unsigned)i<v.size();i++)
	E+=(v[i] - ave)*(v[i] - ave);
	return sqrt(1/static_cast<double>(v.size())*E);
}

int plot_Fun4All_CSV()
{
        ofstream myfile;
//	const std::string inFile = "LeptoAna_100events_tauOnly.root";
//	const std::string inDirectory = "/direct/phenix+u/jlab/github/forks/macros/macros/g4simulations/";
//	const std::string inFile = "LeptoAna_1000events_tauOnly.root";	const std::string class_string = "tau"; myfile.open("./tauJetSummary_1000.csv");
//	const std::string inFile = "LeptoAna_100events_DISonly.root";	const std::string class_string = "DIS"; myfile.open("./tauJetSummary_100.csv");
//	const std::string inFile = "LeptoAna_1000events_DISonly.root";	const std::string class_string = "DIS"; myfile.open("./DISJetSummary_1000.csv");
//	const std::string inFile = "LeptoAna_1000events_DISonly_30_50.root";	const std::string class_string = "DIS"; myfile.open("./data/DISJetSummary_1000_30_50.csv");
//	const std::string inFile = "LeptoAna_1000events_tauOnly_30_50.root";	const std::string class_string = "tau"; myfile.open("./data/tauJetSummary_1000_30_50.csv");
	const std::string inFile = "LeptoAna_1000events_DISonly_5_275.root";	const std::string class_string = "DIS"; myfile.open("./data/DISJetSummary_1000_5_275_small.csv");
//	const std::string inFile = "LeptoAna_1000events_tauOnly_5_275.root";	const std::string class_string = "tau"; myfile.open("./data/tauJetSummary_1000_5_275_small.csv");
//	const std::string inFile = "LeptoAna_1000events_tauOnly_30_275.root";	const std::string class_string = "tau"; myfile.open("./data/tauJetSummary_1000_30_275_small.csv");

	const std::string inDirectory = "/gpfs/mnt/gpfs02/phenix/scratch/jlab/Leptoquark/";
	std::string inputFile = inDirectory+inFile;

	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("ntp_leptoquark");
//	t->Print();

	const int Nevent = t->GetMaximum("event");
	cout << "Running " << Nevent << " events" << endl;
	const int Nentries = t->Draw("towereta:towerphi:towerenergy:event","isMaxEnergyJet==1","goff");

        

	//output what each of the variables in the csv file will contain"
//	myfile << "# " << "n_Total, n_Above_0p001, n_Above_0p01, n_Above_0p1, n_Above_1, n_Above_1, n_Above_10, eta_avg, eta_std, phi_avg, phi_std,"
//		<< " Delta_eta_avg, Delta_eta_std, Delta_phi_avg, Delta_phi_std, Delta_eta_avg_w, Delta_eta_std_w, Delta_phi_avg_w,"
//		<< " Delta_phi_std_w, towerenergy_sum" << endl;

	//loop over all events
	for(int i = 1; i < Nevent+1; i++)
	{
		vector<double> v_DeltaEta, v_DeltaPhi, v_DeltaTheta; 	//vectors which store the differece in eta/phi/theta from the parent tau/quark
		vector<double> v_Eta, v_Phi, v_TowerEnergy;		//vectors which store absolute eta/phi/theta values
		int n_Above_0p001 = 0, n_Above_0p01 = 0, n_Above_0p1 = 0, n_Above_1 = 0, n_Above_10 = 0, n_Total = 0;	//counters for the number of events above 
															//	a certain energy threshod in an event. 0p001 -> 0.001
		double tower_energy_sum = 0;				//rolling sum of all the energy in a jet.

		double Emax = 0;
		int Emax_i = 0;

		//loop over all entries in t->Draw()
		for(int j = 0; j < Nentries; j++)
		{
			if(t->GetV4()[j] == i)
			{
				//increment counters depending on energy level
				n_Total++;
				tower_energy_sum = tower_energy_sum + t->GetV3()[j];
				if(t->GetV3()[j] > 0.001) n_Above_0p001++;
				if(t->GetV3()[j] > 0.01) n_Above_0p01++;
				if(t->GetV3()[j] > 0.1) n_Above_0p1++;
				if(t->GetV3()[j] > 1.0) n_Above_1++;
				if(t->GetV3()[j] > 10.0) n_Above_10++;
				if(t->GetV3()[j] > Emax) 
				{
					Emax = t->GetV3()[j];
					Emax_i = j;
				}
			}
		}
//		cout << "Maximum energy tower for event " << i << " is " << Emax_i << endl;

		for(int j = 0; j < Nentries; j++)
		{
			if(t->GetV4()[j] == i)
			{
				v_Eta.push_back(t->GetV1()[j]);
				v_Phi.push_back(t->GetV2()[j]);
				v_DeltaEta.push_back(t->GetV1()[j] - t->GetV1()[Emax_i]);
				v_DeltaTheta.push_back(2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[j])) - 2*TMath::ATan(TMath::Power(TMath::E(),-1*t->GetV1()[Emax_i])));
				v_DeltaPhi.push_back(t->GetV2()[j] - t->GetV2()[Emax_i]);
				v_TowerEnergy.push_back(t->GetV3()[j]);
			}
		}

	//compute standard deviations and averages of raw values
	double eta_average = Average(v_Eta);
	double eta_std = Deviation(v_Eta, eta_average);
	double phi_average = Average(v_Phi);
	double phi_std = Deviation(v_Phi, phi_average);

	//avg/std of the difference from tau/quark
	double Delta_eta_average = Average(v_DeltaEta);
	double Delta_eta_std = Deviation(v_DeltaEta, Delta_eta_average);
	double Delta_phi_average = Average(v_DeltaPhi);
	double Delta_phi_std = Deviation(v_DeltaPhi, Delta_phi_average);

	//avg/std of difference from parent quark, but weighted by energy
	double Delta_eta_average_weighted = Average(v_DeltaEta,v_TowerEnergy);
	double Delta_eta_std_weighted = Deviation(v_DeltaEta, Delta_eta_average_weighted);
	double Delta_phi_average_weighted = Average(v_DeltaPhi,v_TowerEnergy);
	double Delta_phi_std_weighted = Deviation(v_DeltaPhi, Delta_phi_average_weighted);

	//output to the csv file
//	myfile	<< n_Total << ", " << n_Above_0p001 << ", " << n_Above_0p01 << ", " << n_Above_0p1 << ", " << n_Above_1 << ", " << n_Above_10 << ", " << eta_average << ", "
//		<< eta_std << ", " << phi_average << ", " << phi_std  << ", " << Delta_eta_average << ", " << Delta_eta_std << ", " << Delta_phi_average << ", " 
//		<< Delta_phi_std << ", " << Delta_eta_average_weighted << ", " << Delta_eta_std_weighted << ", " << Delta_phi_average_weighted << ", " 
//		<< Delta_phi_std_weighted << ", " << tower_energy_sum << ", " << class_string << endl;
        myfile  << n_Total << ", " << n_Above_1 << ", " << n_Above_10 << ", " << Delta_eta_average_weighted << ", " << Delta_eta_std_weighted << ", " 
		<< Delta_phi_average_weighted << ", " << Delta_phi_std_weighted << ", " << tower_energy_sum << ", " << class_string << endl; 

	}


	return 0;
}
