#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"

using namespace std;

int plot_LQGENEP_tauOnly()
{

	gStyle->SetOptStat(0);

	std::string inputFile = "./outdir/LQGENEP_20_250_q1_q1_1936-5_1000000.root";
//	std::string inputFile = "./outdir/LQGENEP_30_50_q1_q1_1936-5_1000000.root";
//	std::string inputFile = "./outdir/LQGENEP_5_275_q1_q1_1936-5_1000000.root";
//	std::string inputFile = "./outdir/TestOut.100000event.root";
	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("EICTree");

	int vsize = t->Draw("particles.eta:particles.phi:particles.pt:particles.p","(particles.id==15)*(particles.orig==8)","goff");

	double xmin = -2.5;
	double xmax = 7.5;
	double ymin = -0.1;
	double ymax = 2*TMath::Pi()+0.1;

	//-----------------------------------------------------------------------------------------------------------

	std::string title = "";
//	std::string title = "#eta vs. #phi for #tau^{-} produced in LQ Events";
//	TCanvas *c1 = new TCanvas("c1","canvas",1000,750);
	TCanvas *c1 = new TCanvas();
	TH2D *h1 = new TH2D("h1",title.c_str(),50,xmin,xmax,60,ymin,ymax);
	for(int i = 0; i < vsize; i++)
	{
		h1->Fill(t->GetV1()[i],t->GetV2()[i]);
	}
	h1->Draw("colz");
		h1->SetMinimum(0);
		h1->SetMaximum(3600);
		h1->GetXaxis()->SetTitle("#eta");
		h1->GetYaxis()->SetTitle("#phi");
		h1->GetYaxis()->SetTitleOffset(1.4);
	c1->Update();

	//-----------------------------------------------------------------------------------------------------------

	TCanvas *c2 = new TCanvas();
	std::string title2 = "";
//	std::string title2 = "#eta vs. p_{T} for #tau^{-} produced in LQ Events";
	TH2D *h2 = new TH2D("h2",title2.c_str(),50,xmin,xmax,75,0,75);
	for(int i = 0; i < vsize; i++)
	{
		h2->Fill(t->GetV1()[i],t->GetV3()[i]);
	}
	h2->Draw("colz");
		h2->SetMinimum(0);
//		h2->SetMaximum(14500);
		h2->GetXaxis()->SetTitle("#eta");
		h2->GetYaxis()->SetTitle("Transverse Momentum (GeV)");
	c2->Update();

	//-----------------------------------------------------------------------------------------------------------

	TCanvas *c4 = new TCanvas();
//	std::string title4 = "p_{T} vs. #Sigma p for #tau^{-} produced in LQ Events";
	std::string title4 = "";
	TH2D *h4 = new TH2D("h4",title4.c_str(),44,0,66,44,0,66);
	for(int i = 0; i < vsize; i++)
	{
		h4->Fill(t->GetV3()[i],t->GetV4()[i]);
	}
	h4->Draw("colz");
		h4->SetMinimum(0);
//		h4->SetMaximum(22000);
		h4->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
		h4->GetYaxis()->SetTitle("Total Momentum (GeV)");
	c4->Update();

	//-----------------------------------------------------------------------------------------------------------
/*
	TCanvas *c3 = new TCanvas();
		c3->SetGridx();
		c3->SetGridy();
//	gPad->DrawFrame(-11.2,-11.2,11.2,11.2);
	TPad *p = (TPad*)c3->cd();
	p->SetTheta(90.);
	p->SetPhi(0.);
	t->Draw("11:2*TMath::ATan( TMath::E()**(-1*particles.eta))","particles.id==15"," lego2 pol");
	c3->Update();
*/	
	return 0;
}























