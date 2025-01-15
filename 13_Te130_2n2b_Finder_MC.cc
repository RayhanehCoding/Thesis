// Author : Rayhaneh Dehghani
// Compile the code in Root; using ".L 13_Te130_2n2b_Finder_MC.cc" 
// For the purpose of the Master's degree thesis

// This script is for finding the 2-neutrino 2-beta signal in MC for selected run ranges
// (that describe the detector state in terms of reconstruction and physics)
// Selection cuts will be defined ; energy, Nhits_low, Nhits_high, FV, FitValid, event_index, etc. 




#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "time.h"
#include "TTree.h"
#include "TPad.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>


using namespace std;

int neutrino_finder()
{

	TH1D *Energy2n2b = new TH1D("Te130_2n2b_Energy", "2n2b_Energy", 500, 0.0, 5);
        Energy2n2b -> GetXaxis() -> SetTitle("2n2b_MC_Energy_spectrum");
        Energy2n2b -> GetYaxis() -> SetTitle("event_count");

        TH2D *Energy_posr = new TH2D("2n2bEnergy_posr","2n2bEnergy_posr", 200, 0, 2.7, 500, 0, 6000);
        Energy_posr -> GetXaxis() -> SetTitle("Energy");
        Energy_posr -> GetYaxis() -> SetTitle("posr");
	//-------------------------------------------------------------------------------------------
	TString Directory = "/data/snoplus/production/scint/7.0.9/ntuple/ScintFit_2p2Te130_2n2bRun/";
	TString filename = "ScintFit_2p2Te130_2n2bRun_r";
	TString filename_2 = "_s0_p0.ntuple.root";
	//TString MC_2n2bTe130_path;
	//TString File_name;
	//-------------------------------------------------------------------------------------------
        int NhL=20;
       	int Nhits_low;
	int NhH = 700;
        int Nhits_high ;
	int F = 5500;//[mm]
        int FV_cut = 3300;//[mm]
	
        int nhitsCleaned;
        double EL = 0;//[Mev]
	double energy_low = EL;
	double EH = 2.5;//[Mev]
        double energy_high = EH;
	
	Double_t mcPosr;
        double posr, posx, posy, posz;
        bool fitValid;
	double energy;
        
	int passed_evIndex = 0;
        int passed_fitValid = 0;
        int passed_mcPosr = 0;
        int passed_posr = 0;
	//-------------------------------------------------------------------------------------------
	vector<double> Te2n2b_energy_vec, Te2n2b_nhits_vec;	
		
	for (int Q=300000 ; Q<306498 ; Q=Q+1)
	//for (int Q=300000 ; Q<301067 ; Q=Q+1)
	{
	char numstr[100];
	sprintf(numstr, "%d", Q);
	TString MC_2n2bTe130_path = Directory + filename + numstr + filename_2 ;
		
	ifstream input(MC_2n2bTe130_path);
        if (input.is_open())
        {
	cout << MC_2n2bTe130_path << endl;

	TChain *chain = new TChain ("output");
        chain -> Add(MC_2n2bTe130_path);

	int evIndex;
	chain -> SetBranchAddress("evIndex",&evIndex);//Only for MC files to control re-trigger plotting
	chain -> SetBranchAddress("energy",&energy);
        chain -> SetBranchAddress("posr",&posr);
        chain -> SetBranchAddress("posx",&posx);
        chain -> SetBranchAddress("posy",&posy);
        chain -> SetBranchAddress("posz",&posz);
        chain -> SetBranchAddress("fitValid",&fitValid);
        chain -> SetBranchAddress("nhitsCleaned",&nhitsCleaned);
        chain -> SetBranchAddress("mcPosr",&mcPosr);
       	int entries=chain->GetEntries();
	//-------------------------------------------------------------------------------------------
		for (int i=0; i <= entries; i=i+1)
		{
			chain->GetEntry(i);
			// condition #1
                        if (evIndex <= 0) {
                                passed_evIndex++; // Increment counter if condition is true

                                //condition #2
                                if (fitValid ==1) {
                                        passed_fitValid++; // Increment counter if condition is true

                                        // Condition #3
                                        if (mcPosr <= 6000) {
                                                passed_mcPosr++; // Increment counter if condition is true

                                                // Condition #4
                                                if (posr <= FV_cut) {
                                                        passed_posr++; // Increment counter if condition is true  
			
							Energy2n2b -> Fill(energy);//fill with what ever I want to fill the histogram; i.e energy, posr
						}
					}
				}
			}
		}
		input.close();
               	delete chain;
	}
	else cout <<"file does not exist" << endl;
	}
	cout << "Number of events passing [evIndex] cut = " << passed_evIndex << endl;
        cout << "Number of events passing [evIndex and fitValid] cuts = " << passed_fitValid << endl;
        cout << "Number of events passing [evIndex and fitValid and mcPosr] cuts = " << passed_mcPosr << endl;
        cout << "Number of events passing [evIndex and fitValid and mcPosr and posr] cuts = " << passed_posr << endl;	
        //----------------------------------------------------------------------------------------
	TFile* Te1302n2b_energy_file = new TFile("Te1302n2b_spectrum_2p2_optics.root", "Recreate");//done

	Energy2n2b -> Sumw2();
	//Ploting histograms;
	//-------------------------
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1 -> SetTickx();
	c1 -> SetTicky();
	c1 -> SetGridx();
	c1 -> SetGridy();
	Energy2n2b -> Draw();
	//Po_energy -> SetFillColor(kViolet);
	//Po_energy -> Fit ("Fit_gaus");
	//Energy2n2b -> Scale(1.0/Energy2n2b->Integral());
	Energy2n2b -> Write();
return 0;


}
