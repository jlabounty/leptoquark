int plot_Fun4All_tauOnly_DeltaEta_DeltaPhi()
{
	const std::string inFile = "LeptoAna_10events_tauDaughtersOnly.root";
	const std::string inDirectory = "/direct/phenix+u/jlab/github/forks/macros/macros/g4simulations/";
	std::string inputFile = inDirectory+inFile;

	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("ntp_leptoquark");
//	t->Print();

	const int Nevent = t->GetMaximum("event");
	cout << "Running " << Nevent << " events" << endl;
	const int Nentries = t->Draw("towereta:towerphi:towerenergy:event","isMaxEnergyJet==1","goff");
	vector<double> v_DeltaEta, v_DeltaPhi;

	for(int i = 0; i < Nevent; i++)
	{
		double Emax = 0;
		int Emax_i = 0;
		for(int j = 0; j < Nentries; j++)
		{
			if(t->GetV4()[j] == i)
			{
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
				v_DeltaEta.push_back(t->GetV1()[j] - t->GetV1()[Emax_i]);
				v_DeltaPhi.push_back(t->GetV2()[j] - t->GetV2()[Emax_i]);
			}
		}
	}


//-----------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);

	double xmin = -1;
	double xmax = 1;
	double ymin = -1;
	double ymax = 1;

	std::string title = "";
	TCanvas *c1 = new TCanvas();
	TH2D *h1 = new TH2D("h1",title.c_str(),40,xmin,xmax,40,ymin,ymax);
	for(int i = 0; i < Nentries; i++)
	{
		h1->Fill(v_DeltaEta[i],v_DeltaPhi[i]);
	}
	h1->Draw("colz");
//		h1->SetMinimum(0);
//		h1->SetMaximum(3600);
		h1->GetXaxis()->SetTitle("#Delta#eta");
		h1->GetYaxis()->SetTitle("#Delta#phi");
//		h1->GetYaxis()->SetTitleOffset(1.4);
	c1->Update();

//-----------------------------------------------------------------------------------------------------------


	TCanvas *c2 = new TCanvas();
	TH1D *h2 = h1->ProjectionX();
	h2->Draw();


//-----------------------------------------------------------------------------------------------------------

	return 0;
}
