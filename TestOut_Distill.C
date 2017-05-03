int TestOut_Distill()
{

	const std::string inputFile = "/direct/phenix+u/jlab/Leptoquark/Leptoquark_Lorenzo/TestOut.txt";
	std::ifstream fin(inputFile.c_str());
	std::string str;

	std::string foutname = "./outdir/TestOut_Distilled_DIS_10_50.txt";
	std::ofstream fout;
	fout.open(foutname.c_str());

	int I, KS, KF, ORIG;
	const int tau_KF = 1;
	int I_prev = 0;
	int event_count = 0;

	vector<int> tau_orig;			
	while (std::getline(fin, str))
	{
		if(str.length() == 127)
		{
//			cout << str << endl;
			std::istringstream iss(str);

			if(iss >> I >> KS >> KF >> ORIG)
			{
				if(I < I_prev) tau_orig.clear();
				I_prev = I;

				if(I < 3)
				{
					fout << str << endl;
				}
				else if(KS == 21) 
					{ fout << str << endl; }
				else if(TMath::Abs(KF) == tau_KF)
				{
					fout << str << endl;
					tau_orig.push_back(I);
				}
				else
				{
					for(int i = 0; i < tau_orig.size(); i++)
					{
						if(ORIG == tau_orig[i]) 
						{
							fout << str << endl;
							tau_orig.push_back(I);
						}
					}
				}
				
			}
		}
		else
		{
			fout << str << endl;
			if(str.find("Event finished") != std::string::npos) event_count++;
		}
	}

	gSystem->Load("libeicsmear");
	BuildTree(foutname.c_str(),"./outdir/",event_count);

	return 0;
}
