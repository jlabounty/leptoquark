int plot_tauOnly_eta_v_phi()
{
	std::string inputFile = "./outdir/TestOut.100000event.root";
	
	plot_tauOnly_eta_v_phi(inputFile);

	return 0;
}


int plot_tauOnly_eta_v_phi(std::string inputFile)
{
	TFile *f = TFile::Open(inputFile.c_str());
	TTree *t = (TTree*)f->Get("EICTree");

	TCanvas *c1 = new TCanvas();
	t->Draw("particles.eta:particles.phi","particles.id==15","colz");

	return 0;
}
