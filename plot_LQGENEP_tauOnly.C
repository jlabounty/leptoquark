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

	std::string inputFile = "./outdir/TestOut.1000000event.root";
//	std::string inputFile = "./outdir/TestOut.10000event.root";
	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("EICTree");

	int vsize = t->Draw("particles.eta:particles.phi:particles.pt","(particles.id==15)*(particles.orig==8)","goff");

	double xmin = -2.5;
	double xmax = 7.5;
	double ymin = -0.1;
	double ymax = 2*TMath::Pi()+0.1;

	std::string title = "#eta vs. #phi for #tau^{-} produced in LQ Events";
	TCanvas *c1 = new TCanvas();
	TH2D *h1 = new TH2D("h1",title.c_str(),100,xmin,xmax,100,ymin,ymax);
	for(int i = 0; i < vsize; i++)
	{
		h1->Fill(t->GetV1()[i],t->GetV2()[i]);
	}
	h1->Draw("colz");
		h1->GetXaxis()->SetTitle("#eta");
		h1->GetYaxis()->SetTitle("#phi");
	c1->Update();

	TCanvas *c2 = new TCanvas();
	std::string title2 = "#eta vs. p_{T} for #tau^{-} produced in LQ Events";
	TH2D *h2 = new TH2D("h2",title2.c_str(),100,xmin,xmax,100,0,100);
	for(int i = 0; i < vsize; i++)
	{
		h2->Fill(t->GetV1()[i],t->GetV3()[i]);
	}
	h2->Draw("colz");
		h2->GetXaxis()->SetTitle("#eta");
		h2->GetYaxis()->SetTitle("Transverse Momentum (GeV)");
	c2->Update();
	
	return 0;
}
