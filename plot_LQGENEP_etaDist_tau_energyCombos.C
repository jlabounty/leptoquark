int plot_LQGENEP_etaDist_tau_energyCombos()
{

        gStyle->SetOptStat(0);


	int vsize = 0;
        std::string title = "#eta distribution for #tau^{-} produced in LQ Events";
        TCanvas *c1 = new TCanvas();
        double xmin = -2.5;
        double xmax = 7.5;

	//----------------------------------------------------------------------------------------------------------

        std::string inputFile3 = "./outdir/tau_dist_comparison_10000events/e5_p275_LQ.root";
        TFile *f3 = TFile::Open(inputFile3.c_str());
        TTree *t3 = (TTree*)f3->Get("EICTree");

        vsize = t3->Draw("particles.eta:particles.phi:particles.pt","(particles.id==15)*(particles.orig==8)","goff");

        TH1F *h3 = new TH1F("h3",title.c_str(),100,xmin,xmax);
        for(int i = 0; i < vsize; i++)
        {
                h3->Fill(t3->GetV1()[i]);
        }
        h3->Draw("colz");
                h3->GetXaxis()->SetTitle("#eta");
                h3->GetYaxis()->SetTitle("Number of Particles");
                h3->GetYaxis()->SetTitleOffset(1.6);
		h3->SetLineColor(kGreen);
		h3->SetFillColor(kGreen);
		h3->SetFillStyle(3003);
        c1->Update();

	//----------------------------------------------------------------------------------------------------------

        std::string inputFile2 = "./outdir/tau_dist_comparison_10000events/e30_p50_LQ.root";
        TFile *f2 = TFile::Open(inputFile2.c_str());
        TTree *t2 = (TTree*)f2->Get("EICTree");

        vsize = t2->Draw("particles.eta:particles.phi:particles.pt","(particles.id==15)*(particles.orig==8)","goff");

        TH1F *h2 = new TH1F("h2",title.c_str(),100,xmin,xmax);
        for(int i = 0; i < vsize; i++)
        {
                h2->Fill(t2->GetV1()[i]);
        }
        h2->Draw("colz same");
		h2->SetLineColor(kBlue);
		h2->SetFillColor(kBlue);
		h2->SetFillStyle(3003);
        c1->Update();

	//----------------------------------------------------------------------------------------------------------
	
        std::string inputFile = "./outdir/tau_dist_comparison_10000events/e20_p250_LQ.root";
        TFile *f = TFile::Open(inputFile.c_str());
        TTree *t = (TTree*)f->Get("EICTree");

        vsize = t->Draw("particles.eta:particles.phi:particles.pt","(particles.id==15)*(particles.orig==8)","goff");

        TH1F *h1 = new TH1F("h1",title.c_str(),100,xmin,xmax);
        for(int i = 0; i < vsize; i++)
        {
                h1->Fill(t->GetV1()[i]);
        }
        h1->Draw("colz same");
		h1->SetLineColor(kRed);
		h1->SetFillColor(kRed);
		h1->SetFillStyle(3003);
        c1->Update();

	//----------------------------------------------------------------------------------------------------------

	leg = new TLegend(0.55,0.9,0.9,0.75); 
		leg->SetBorderSize(1);
		leg->AddEntry(h1,"20 Gev x 250 Gev","l");
		leg->AddEntry(h2,"30 Gev x 50 Gev","l");
		leg->AddEntry(h3,"5 Gev x 275 Gev","l");
	leg->Draw("same");
        c1->Update();

	return 0;
}
