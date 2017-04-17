int DetectorAnalysis()
{

//	std::string inputFile = "./data/LeptoAna.root";
	std::string inputFile = "/direct/phenix+u/jlab/github/sPHENIX/macros/macros/g4simulations/LeptoAna.root";
	TFile *f = TFile::Open(inputFile.c_str());
        TNtuple *t = (TNtuple*)f->Get("ntp_leptoquark");
//	t->Print();i

	const int vsize = t->Draw("towereta:towerphi:towerenergy:event","calorimeterid == 1","goff");
	int nevents = t->GetMaximum("event") + 1;

	double max_energy_tower_i = 0;
	double max_eta = 0, max_phi = 0;
	vector<double> vec_diff_eta(vsize), vec_diff_phi(vsize);

	for(int i = 0; i < nevents; i++)
	{
		for(int j = 0; j < vsize; j++)
		{
			if(i == (int)t->GetV4()[j])
			{
				if(t->GetV3()[j] > max_energy_tower_i)
				{
					max_energy_tower_i = t->GetV3()[j];
					max_eta = t->GetV1()[j];
					max_phi = t->GetV2()[j];
				}
				
			}
		}

		for(int j = 0; j < vsize; j++)
		{
			if(i == (int)t->GetV4()[j])
			{
				vec_diff_eta[j] = (double)( t->GetV1()[j] - max_eta);
				vec_diff_phi[j] = t->GetV2()[j] - max_phi;
			}
		}
		max_energy_tower_i = 0;
	}


	//-----------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas();
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.4);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	TH2F *h  = new TH2F("h","#Delta#eta vs. #Delta#phi",48,-0.6,0.6,48,-0.6,0.6);
	for(int i = 0; (unsigned)i < vec_diff_eta.size(); i++)
	{
		h->Fill(vec_diff_eta[i],vec_diff_phi[i], t->GetV3()[i]);
	}
	h->SetTitle("#Delta#eta vs. #Delta#phi for Towers in Primary Jet (LQ)");
	h->GetXaxis()->SetTitle("#Delta#eta");
	h->GetYaxis()->SetTitle("#Delta#phi");
	h->GetZaxis()->SetTitle("Number of Towers (Energy Weighted)");
	h->Draw("COLZ");
	c1->SetRightMargin(0.15);
//	c1->Print("DeltaEta_vs_DeltaPhi.eps");
//	c1->Print("DeltaEta_vs_DeltaPhi.root");


	//-----------------------------------------------------------------------------------------

	TCanvas *c2 = new TCanvas();
	TH1D *h2 = h->ProjectionX();
	h2->Draw();
	h2->SetTitle("#Delta#eta for Towers in Primary Jet (LQ)");
	c2->Update();

	TCanvas *c3 = new TCanvas();
	TH1D *h3 = h->ProjectionY();
	h3->Draw();
	h3->SetTitle("#Delta#phi for Towers in Primary Jet (LQ)");
	c3->Update();
















	return 0;
}
